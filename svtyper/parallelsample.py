#!/usr/bin/env python3

import sys
import gzip
import argparse
import shutil
import tempfile
from argparse import RawTextHelpFormatter
import os
import multiprocessing
from functools import partial
from subprocess import call
import random
import string


from pysam import AlignmentFile, tabix_compress, tabix_index
from svtyper import singlesample as single
from svtyper import version as svtyper_version


def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
        svtyper\n\
        author: " + svtyper_version.__author__ + "\n\
        version: " + svtyper_version.__version__ + "\n\
        description: Compute genotype of structural variants based on breakpoint depth")

    parser.add_argument('-B', '--bam',
                        type=str, required=True,
                        help='BAM list files')
    parser.add_argument('-i', '--input_vcf',
                        type=str, required=False, default=None,
                        help='VCF input')
    parser.add_argument('-o', '--output_vcf',
                        type=str, required=False, default=None,
                        help='output VCF to write')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='number of threads to use (set at maximum number of available cores)')

    args = parser.parse_args()

    if args.input_vcf is None:
        args.input_vcf = sys.stdin
    if args.output_vcf is None:
        args.output_vcf = sys.stdout

    return args

def fetchId(bamfile):
    """
    Fetch sample id in a bam file
    :param bamfile: the bam file
    :type bamfile: file
    :return: sample name
    :rtype: str
    """
    bamfile_fin = AlignmentFile(bamfile, 'rb')
    name = bamfile_fin.header['RG'][0]['SM']
    bamfile_fin.close()
    return name


def get_bamfiles(bamlist):
    with open(bamlist, "r") as fin:
        bamfiles = [os.path.abspath(f.strip()) for f in fin]
    return bamfiles


def random_string(length=10):
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(length))


def genotype_multiple_samples(bamlist, vcf_in, vcf_out, cores=1):

    bamfiles = get_bamfiles(bamlist)

    if vcf_out == sys.stdout:
        tmp_dir = tempfile.mkdtemp()
    else:
        tmp_dir = os.path.join(os.path.dirname(vcf_out), "tmp" + random_string())
        os.mkdir(tmp_dir)

    if vcf_in == sys.stdin:
        vcf_in = single.dump_piped_vcf_to_file(vcf_in, tmp_dir)
    elif vcf_in.endswith(".gz"):
        vcf_in_gz = vcf_in
        vcf_in = os.path.join(tmp_dir, os.path.basename(vcf_in[:-3]))
        with gzip.open(vcf_in_gz, 'rb') as f_in, open(vcf_in, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    pool = multiprocessing.Pool(processes=cores)
    launch_ind = partial(genotype_single_sample,
                         vcf_in=vcf_in, out_dir=tmp_dir)

    vcf_files = pool.map(launch_ind, bamfiles)

    merge_cmd = "bcftools merge -m id "

    if vcf_out != sys.stdout:
        merge_cmd += "-O z " + " -o " + vcf_out + " "

    merge_cmd += " ".join(vcf_files)

    exit_code = call(merge_cmd, shell=True)

    if exit_code == 0:
        if vcf_out != sys.stdout:
            tabix_index(vcf_out, force=True, preset="vcf")

        shutil.rmtree(tmp_dir)
    else:
        print("Failed: bcftools merge exits with status %d" % exit_code)
        exit(1)
        

def genotype_single_sample(bam, vcf_in, out_dir):
    lib_info_json = bam + ".json"
    sample = fetchId(bam)
    out_vcf = os.path.join(out_dir, sample + ".gt.vcf")
    with open(vcf_in, "r") as inf, open(out_vcf, "w") as outf:
        single.sso_genotype(bam_string=bam,
                            vcf_in=inf,
                            vcf_out=outf,
                            min_aligned=20,
                            split_weight=1,
                            disc_weight=1,
                            num_samp=1000000,
                            lib_info_path=lib_info_json,
                            debug=False,
                            ref_fasta=None,
                            sum_quals=False,
                            max_reads=1000,
                            cores=None,
                            batch_size=None)
    out_gz = out_vcf + ".gz"
    tabix_compress(out_vcf, out_gz, force=True)
    tabix_index(out_gz, force=True, preset="vcf")
    return out_gz


def main():
    args = get_args()
    
    genotype_multiple_samples(args.bam, args.input_vcf, args.output_vcf, cores=args.threads)


def cli():
    try:
        sys.exit(main())
    except IOError as e:
        if e.errno != 32:  # ignore SIGPIPE
            raise


if __name__ == '__main__':
    cli()

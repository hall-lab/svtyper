#!/usr/bin/env python

import logging
import pysam
import json
import argparse, sys, os.path
import math, time
from argparse import RawTextHelpFormatter

import svtyper.version

from svtyper.parsers import Vcf, Variant, Sample
from svtyper.utils import *
from svtyper.statistics import bayes_gt

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
svtyper\n\
author: " + svtyper.version.__author__ + "\n\
version: " + svtyper.version.__version__ + "\n\
description: Compute genotype of structural variants based on breakpoint depth")
    parser.add_argument('-i', '--input_vcf', metavar='FILE', type=argparse.FileType('r'), default=None, help='VCF input (default: stdin)')
    parser.add_argument('-o', '--output_vcf', metavar='FILE',  type=argparse.FileType('w'), default=sys.stdout, help='output VCF to write (default: stdout)')
    parser.add_argument('-B', '--bam', metavar='FILE', type=str, required=True, help='BAM or CRAM file(s), comma-separated if genotyping multiple samples')
    parser.add_argument('-T', '--ref_fasta', metavar='FILE', type=str, required=False, default=None, help='Indexed reference FASTA file (recommended for reading CRAM files)')
    parser.add_argument('-S', '--split_bam', type=str, required=False, help=argparse.SUPPRESS)
    parser.add_argument('-l', '--lib_info', metavar='FILE', dest='lib_info_path', type=str, required=False, default=None, help='create/read JSON file of library information')
    parser.add_argument('-m', '--min_aligned', metavar='INT', type=int, required=False, default=20, help='minimum number of aligned bases to consider read as evidence [20]')
    parser.add_argument('-n', dest='num_samp', metavar='INT', type=int, required=False, default=1000000, help='number of reads to sample from BAM file for building insert size distribution [1000000]')
    parser.add_argument('-q', '--sum_quals', action='store_true', required=False, help='add genotyping quality to existing QUAL (default: overwrite QUAL field)')
    parser.add_argument('--max_reads', metavar='INT', type=int, default=None, required=False, help='maximum number of reads to assess at any variant (reduces processing time in high-depth regions, default: unlimited)')
    parser.add_argument('--split_weight', metavar='FLOAT', type=float, required=False, default=1, help='weight for split reads [1]')
    parser.add_argument('--disc_weight', metavar='FLOAT', type=float, required=False, default=1, help='weight for discordant paired-end reads [1]')
    parser.add_argument('-w', '--write_alignment', metavar='FILE', dest='alignment_outpath', type=str, required=False, default=None, help='write relevant reads to BAM file')
    parser.add_argument('--debug', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--verbose', action='store_true', default=False, help='Report status updates')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.input_vcf == None:
        if not sys.stdin.isatty():
            args.input_vcf = sys.stdin

    # send back the user input
    return args

# methods to grab reads from region of interest in BAM file
def gather_all_reads(sample, chromA, posA, ciA, chromB, posB, ciB, z, max_reads):
    # grab batch of reads from both sides of breakpoint
    read_batch = {}
    read_batch, many = gather_reads(sample, chromA, posA, ciA, z, read_batch, max_reads)
    if many:
        return {}, True

    read_batch, many = gather_reads(sample, chromB, posB, ciB, z, read_batch, max_reads)
    if many:
        return {}, True

    return read_batch, many

def gather_reads(sample,
                 chrom, pos, ci,
                 z,
                 fragment_dict,
                 max_reads):

    # the distance to the left and right of the breakpoint to scan
    # (max of mean + z standard devs over all of a sample's libraries)
    fetch_flank = sample.get_fetch_flank(z)
    chrom_length = sample.bam.lengths[sample.bam.gettid(chrom)]

    many = False
    for i, read in enumerate(sample.bam.fetch(chrom,
                                 max(pos + ci[0] - fetch_flank, 0),
                                 min(pos + ci[1] + fetch_flank, chrom_length))):
        if read.is_unmapped or read.is_duplicate:
            continue

        lib = sample.get_lib(read.get_tag('RG'))
        if lib.name not in sample.active_libs:
            continue

        # read.query_sequence = "*"
        # read.query_qualities = "*"
        if max_reads is not None and i > max_reads:
            many = True
            break

        if read.query_name in fragment_dict:
            fragment_dict[read.query_name].add_read(read)
        else:
            fragment_dict[read.query_name] = SamFragment(read, lib)

    return fragment_dict, many


# ==================================================
# Genotyping function
# ==================================================

def sv_genotype(bam_string,
                vcf_in,
                vcf_out,
                min_aligned,
                split_weight,
                disc_weight,
                num_samp,
                lib_info_path,
                debug,
                alignment_outpath,
                ref_fasta,
                sum_quals,
                max_reads):

    # parse the comma separated inputs
    bam_list = []
    for b in bam_string.split(','):
        if b.endswith('.bam'):
            bam_list.append(pysam.AlignmentFile(b, mode='rb'))
        elif b.endswith('.cram'):
            bam_list.append(pysam.AlignmentFile(b, mode='rc', reference_filename=ref_fasta))
        else:
            sys.stderr.write('Error: %s is not a valid alignment file (*.bam or *.cram)\n' % b)
            exit(1)
            
    min_lib_prevalence = 1e-3 # only consider libraries that constitute at least this fraction of the BAM

    # parse lib_info_path JSON
    lib_info = None
    if lib_info_path is not None and os.path.isfile(lib_info_path):
        lib_info_file = open(lib_info_path, 'r')
        lib_info = json.load(lib_info_file)

    if vcf_in is None:
        sys.stderr.write('Warning: VCF not found.\n')

    # build the sample libraries, either from the lib_info JSON or empirically from the BAMs
    sample_list = list()
    for i in xrange(len(bam_list)):
        if lib_info is None:
            logging.info('Calculating library metrics from %s...' % bam_list[i].filename)
            sample = Sample.from_bam(bam_list[i], num_samp, min_lib_prevalence)
        else:
            logging.info('Reading library metrics from %s...' % lib_info_path)
            sample = Sample.from_lib_info(bam_list[i], lib_info, min_lib_prevalence)

        sample.set_exp_seq_depth(min_aligned)
        sample.set_exp_spanning_depth(min_aligned)
        sample_list.append(sample)
    logging.info('done')

    # diagnostic dump of relevant BAM reads
    if alignment_outpath is not None:
        # create a BAM file of the reads used for genotyping
        out_bam_written_reads = set()
        template_bam = pysam.AlignmentFile(bam_string.split(',')[0], 'rb')
        out_bam = pysam.AlignmentFile(alignment_outpath, 'wb', template_bam)
        template_bam.close()

    # write the JSON for each sample's libraries
    if lib_info_path is not None and not os.path.isfile(lib_info_path):
        logging.info('Writing library metrics to %s...' % lib_info_path)
        lib_info_file = open(lib_info_path, 'w')
        write_sample_json(sample_list, lib_info_file)
        lib_info_file.close()
        logging.info('done')

    # quit early if VCF absent
    if vcf_in is None:
        if alignment_outpath is not None:
            out_bam.close()
        return

    # set variables for genotyping
    z = 3
    split_slop = 3 # amount of slop around breakpoint to count splitters
    in_header = True
    header = []
    breakend_dict = {} # cache to hold unmatched generic breakends for genotyping
    vcf = Vcf()

    # read input VCF
    for line in vcf_in:
        if in_header:
            if line[0] == '#':
                header.append(line)
                if line[1] != '#':
                    vcf_samples = line.rstrip().split('\t')[9:]
                continue
            else:
                in_header = False
                vcf.add_header(header)
                # if detailed:
                vcf.add_custom_svtyper_headers()

                # add the samples in the BAM files to the VCF output
                for sample in sample_list:
                    if sample.name not in vcf.sample_list:
                        vcf.add_sample(sample.name)

                # write the output header
                vcf_out.write(vcf.get_header() + '\n')


        v = line.rstrip().split('\t')
        var = Variant(v, vcf)
        var_length = None # var_length should be None except for deletions
        if not sum_quals:
            var.qual = 0

        # genotype generic breakends
        try:
            svtype = var.get_info('SVTYPE')
        except KeyError:
            sys.stderr.write('Warning: SVTYPE missing at variant %s. Skipping.\n' % (var.var_id))
            vcf_out.write(var.get_var_string() + '\n')
            continue
            
        # print original line if unsupported svtype
        if svtype not in ('BND', 'DEL', 'DUP', 'INV'):
            sys.stderr.write('Warning: Unsupported SVTYPE at variant %s (%s). Skipping.\n' % (var.var_id, svtype))
            vcf_out.write(var.get_var_string() + '\n')
            continue

        if svtype == 'BND':
            if var.info['MATEID'] in breakend_dict:
                var2 = var
                var = breakend_dict[var.info['MATEID']]
                chromA = var.chrom
                chromB = var2.chrom
                posA = var.pos
                posB = var2.pos
                # confidence intervals
                ciA = map(int, var.info['CIPOS'].split(','))
                ciB = map(int, var2.info['CIPOS'].split(','))

                # infer the strands from the alt allele
                if var.alt[-1] == '[' or var.alt[-1] == ']':
                    o1_is_reverse = False
                else: o1_is_reverse = True
                if var2.alt[-1] == '[' or var2.alt[-1] == ']':
                    o2_is_reverse = False
                else: o2_is_reverse = True

                # remove the BND from the breakend_dict
                # to free up memory
                del breakend_dict[var.var_id]
            else:
                breakend_dict[var.var_id] = var
                continue
        else:
            chromA = var.chrom
            chromB = var.chrom
            posA = var.pos
            posB = int(var.get_info('END'))
            # confidence intervals
            ciA = map(int, var.info['CIPOS'].split(','))
            ciB = map(int, var.info['CIEND'].split(','))
            if svtype == 'DEL':
                var_length = posB - posA
                o1_is_reverse, o2_is_reverse =  False, True
            elif svtype == 'DUP':
                o1_is_reverse, o2_is_reverse =  True, False
            elif svtype == 'INV':
                o1_is_reverse, o2_is_reverse = False, False

        # increment the negative strand values (note position in VCF should be the base immediately left of the breakpoint junction)
        if o1_is_reverse: posA += 1
        if o2_is_reverse: posB += 1

        for sample in sample_list:
            # grab reads from both sides of breakpoint
            read_batch, many = gather_all_reads(sample, chromA, posA, ciA, chromB, posB, ciB, z, max_reads)
            if many:
                var.genotype(sample.name).set_format('GT', './.')
                continue

            # initialize counts to zero
            ref_span, alt_span = 0, 0
            ref_seq, alt_seq = 0, 0
            alt_clip = 0

            # ref_ciA = ciA
            # ref_ciB = ciB
            ref_ciA = [0,0]
            ref_ciB = [0,0]

            for query_name in sorted(read_batch.keys()):
                fragment = read_batch[query_name]
                # boolean on whether to write the fragment
                write_fragment = False

                # -------------------------------------
                # Check for split-read evidence
                # -------------------------------------

                # get reference sequences
                for read in fragment.primary_reads:
                    is_ref_seq_A = fragment.is_ref_seq(read, var, chromA, posA, ciA, min_aligned)
                    is_ref_seq_B = fragment.is_ref_seq(read, var, chromB, posB, ciB, min_aligned)
                    if (is_ref_seq_A or is_ref_seq_B):
                        p_reference = prob_mapq(read)
                        ref_seq += p_reference

                        read.set_tag('XV', 'R')
                        write_fragment = True

                # get non-reference split-read support
                for split in fragment.split_reads:

                    split_lr = split.is_split_straddle(chromA, posA, ciA,
                                                       chromB, posB, ciB,
                                                       o1_is_reverse, o2_is_reverse,
                                                       svtype, split_slop)
                    # p_alt = prob_mapq(split.query_left) * prob_mapq(split.query_right)
                    p_alt = (prob_mapq(split.query_left) * split_lr[0] + prob_mapq(split.query_right) * split_lr[1]) / 2.0
                    if split.is_soft_clip:
                        alt_clip += p_alt
                    else:
                        alt_seq += p_alt

                    if p_alt > 0:
                        split.tag_split(p_alt)
                        write_fragment = True

                # -------------------------------------
                # Check for paired-end evidence
                # -------------------------------------

                # tally spanning alternate pairs
                if svtype == 'DEL' and posB - posA < 2 * fragment.lib.sd:
                    alt_straddle = False
                else:
                    alt_straddle = fragment.is_pair_straddle(chromA, posA, ciA,
                                                             chromB, posB, ciB,
                                                             o1_is_reverse, o2_is_reverse,
                                                             min_aligned,
                                                             fragment.lib)

                # check both sides if inversion (perhaps should do this for BND as well?)
                if svtype in ('INV'):
                    alt_straddle_reciprocal = fragment.is_pair_straddle(chromA, posA, ciA,
                                                                        chromB, posB, ciB,
                                                                        not o1_is_reverse,
                                                                        not o2_is_reverse,
                                                                        min_aligned,
                                                                        fragment.lib)
                else:
                    alt_straddle_reciprocal = False

                if alt_straddle or alt_straddle_reciprocal:
                    if svtype == 'DEL':
                        p_conc = fragment.p_concordant(var_length)
                        if p_conc is not None:
                            p_alt = (1 - p_conc) * prob_mapq(fragment.readA) * prob_mapq(fragment.readB)
                            alt_span += p_alt

                            # # since an alt straddler is by definition also a reference straddler,
                            # # we can bail out early here to save some time
                            # p_reference = p_conc * prob_mapq(fragment.readA) * prob_mapq(fragment.readB)
                            # ref_span += p_reference
                            # continue

                            fragment.tag_span(p_alt)
                            write_fragment = True

                    else:
                        p_alt = prob_mapq(fragment.readA) * prob_mapq(fragment.readB)
                        alt_span += p_alt

                        fragment.tag_span(p_alt)
                        write_fragment = True

                # # tally spanning reference pairs
                if svtype == 'DEL' and posB - posA < 2 * fragment.lib.sd:
                    ref_straddle_A = False
                    ref_straddle_B = False
                else:
                    ref_straddle_A = fragment.is_pair_straddle(chromA, posA, ref_ciA,
                                                               chromA, posA, ref_ciA,
                                                               False, True,
                                                               min_aligned,
                                                               fragment.lib)
                    ref_straddle_B = fragment.is_pair_straddle(chromB, posB, ref_ciB,
                                                               chromB, posB, ref_ciB,
                                                               False, True,
                                                               min_aligned,
                                                               fragment.lib)

                if ref_straddle_A or ref_straddle_B:
                    # don't allow the pair to jump the entire variant, except for
                    # length-changing SVs like deletions
                    if not (ref_straddle_A and ref_straddle_B) or svtype == 'DEL':
                        p_conc = fragment.p_concordant(var_length)
                        if p_conc is not None:
                            p_reference = p_conc * prob_mapq(fragment.readA) * prob_mapq(fragment.readB)
                            ref_span += (ref_straddle_A + ref_straddle_B) * p_reference / 2

                            fragment.tag_span(1 - p_conc)
                            write_fragment = True

                # write to BAM if requested
                if alignment_outpath is not None and  write_fragment:
                    for read in fragment.primary_reads + [split.read for split in fragment.split_reads]:
                        out_bam_written_reads = write_alignment(read, out_bam, out_bam_written_reads)

            if debug:
                print '--------------------------'
                print 'ref_span:', ref_span
                print 'alt_span:', alt_span
                print 'ref_seq:', ref_seq
                print 'alt_seq:', alt_seq
                print 'alt_clip:', alt_clip

            # in the absence of evidence for a particular type, ignore the reference
            # support for that type as well
            if (alt_seq + alt_clip) < 0.5 and alt_span >= 1:
                alt_seq = 0
                alt_clip = 0
                ref_seq = 0
            if alt_span < 0.5 and (alt_seq + alt_clip) >= 1:
                alt_span = 0
                ref_span = 0

            if ref_seq + alt_seq + ref_span + alt_span + alt_clip > 0:
                # get bayesian classifier
                if var.info['SVTYPE'] == "DUP": is_dup = True
                else: is_dup = False
                alt_splitters = alt_seq + alt_clip
                QR = int(split_weight * ref_seq) + int(disc_weight * ref_span)
                QA = int(split_weight * alt_splitters) + int(disc_weight * alt_span)
                gt_lplist = bayes_gt(QR, QA, is_dup)
                best, second_best = sorted([ (i, e) for i, e in enumerate(gt_lplist) ], key=lambda(x): x[1], reverse=True)[0:2]
                gt_idx = best[0]

                # print log probabilities of homref, het, homalt
                if debug:
                    print gt_lplist

                # set the overall variant QUAL score and sample specific fields
                var.genotype(sample.name).set_format('GL', ','.join(['%.0f' % x for x in gt_lplist]))
                var.genotype(sample.name).set_format('DP', int(ref_seq + alt_seq + alt_clip + ref_span + alt_span))
                var.genotype(sample.name).set_format('RO', int(ref_seq + ref_span))
                var.genotype(sample.name).set_format('AO', int(alt_seq + alt_clip + alt_span))
                var.genotype(sample.name).set_format('QR', QR)
                var.genotype(sample.name).set_format('QA', QA)
                # if detailed:
                var.genotype(sample.name).set_format('RS', int(ref_seq))
                var.genotype(sample.name).set_format('AS', int(alt_seq))
                var.genotype(sample.name).set_format('ASC', int(alt_clip))
                var.genotype(sample.name).set_format('RP', int(ref_span))
                var.genotype(sample.name).set_format('AP', int(alt_span))
                try:
                    var.genotype(sample.name).set_format('AB', '%.2g' % (QA / float(QR + QA)))
                except ZeroDivisionError:
                    var.genotype(sample.name).set_format('AB', '.')


                # assign genotypes
                gt_sum = 0
                for gt in gt_lplist:
                    try:
                        gt_sum += 10**gt
                    except OverflowError:
                        gt_sum += 0
                if gt_sum > 0:
                    gt_sum_log = math.log(gt_sum, 10)
                    sample_qual = abs(-10 * (gt_lplist[0] - gt_sum_log)) # phred-scaled probability site is non-reference in this sample
                    phred_gq = min(-10 * (second_best[1] - best[1]), 200)
                    var.genotype(sample.name).set_format('GQ', int(phred_gq))
                    var.genotype(sample.name).set_format('SQ', sample_qual)
                    var.qual += sample_qual
                    if gt_idx == 1:
                        var.genotype(sample.name).set_format('GT', '0/1')
                    elif gt_idx == 2:
                        var.genotype(sample.name).set_format('GT', '1/1')
                    elif gt_idx == 0:
                        var.genotype(sample.name).set_format('GT', '0/0')
                else:
                    var.genotype(sample.name).set_format('GQ', '.')
                    var.genotype(sample.name).set_format('SQ', '.')
                    var.genotype(sample.name).set_format('GT', './.')
            else:
                var.genotype(sample.name).set_format('GT', './.')
                var.qual = 0
                var.genotype(sample.name).set_format('GQ', '.')
                var.genotype(sample.name).set_format('SQ', '.')
                var.genotype(sample.name).set_format('GL', '.')
                var.genotype(sample.name).set_format('DP', 0)
                var.genotype(sample.name).set_format('AO', 0)
                var.genotype(sample.name).set_format('RO', 0)
                # if detailed:
                var.genotype(sample.name).set_format('AS', 0)
                var.genotype(sample.name).set_format('ASC', 0)
                var.genotype(sample.name).set_format('RS', 0)
                var.genotype(sample.name).set_format('AP', 0)
                var.genotype(sample.name).set_format('RP', 0)
                var.genotype(sample.name).set_format('QR', 0)
                var.genotype(sample.name).set_format('QA', 0)
                var.genotype(sample.name).set_format('AB', '.')

        # after all samples have been processed, write
        vcf_out.write(var.get_var_string() + '\n')
        if var.info['SVTYPE'] == 'BND':
            var2.qual = var.qual
            var2.active_formats = var.active_formats
            var2.genotype = var.genotype
            vcf_out.write(var2.get_var_string() + '\n')

    # close the files
    vcf_in.close()
    vcf_out.close()
    if alignment_outpath is not None:
        out_bam.close()

    return

def set_up_logging(verbose):
    level = logging.WARNING
    if verbose:
        level = logging.INFO

    logging.basicConfig(format='%(message)s', level=level)

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    set_up_logging(args.verbose)

    if args.split_bam is not None:
        sys.stderr.write('Warning: --split_bam (-S) is deprecated. Ignoring %s.\n' % args.split_bam)

    # call primary function
    sv_genotype(args.bam,
                args.input_vcf,
                args.output_vcf,
                args.min_aligned,
                args.split_weight,
                args.disc_weight,
                args.num_samp,
                args.lib_info_path,
                args.debug,
                args.alignment_outpath,
                args.ref_fasta,
                args.sum_quals,
                args.max_reads)

# --------------------------------------
# command-line/console entrypoint

def cli():
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise

# initialize the script
if __name__ == '__main__':
    cli()

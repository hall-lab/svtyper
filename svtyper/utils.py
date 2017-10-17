from __future__ import print_function

import sys, time, datetime, os, contextlib, tempfile, shutil
from functools import wraps

from svtyper.parsers import SamFragment, Vcf

# ==================================================
# BAM and JSON output
# ==================================================

# write read to BAM file, checking whether read is already written
def write_alignment(read, bam, written_reads, is_alt=None):
    read.query_sequence = None
    read_hash = (read.query_name, read.flag)

    if bam is None or read_hash in written_reads:
        return written_reads
    else:
        bam.write(read)
        written_reads.add(read_hash)
        return written_reads

# dump the sample and library info to a file
def write_sample_json(sample_list, lib_info_file):
    lib_info = {}
    for sample in sample_list:
        s = {}
        s['sample_name'] = sample.name
        s['bam'] = sample.bam.filename
        s['libraryArray'] = []
        s['mapped'] = sample.bam.mapped
        s['unmapped'] = sample.bam.unmapped

        for lib in sample.lib_dict.values():
            l = {}
            l['library_name'] = lib.name
            l['readgroups'] = lib.readgroups
            l['read_length'] = lib.read_length
            l['mean'] = lib.mean
            l['sd'] = lib.sd
            l['prevalence'] = lib.prevalence
            l['histogram'] = lib.hist

            s['libraryArray'].append(l)

        lib_info[sample.name] = s

    # write the json file
    json.dump(lib_info, lib_info_file, indent=4)
    lib_info_file.close()


# ==================================================
# Miscellaneous methods for manipulating SAM alignments
# ==================================================

# get the non-phred-scaled mapq of a read
def prob_mapq(read):
    return 1 - 10 ** (-read.mapping_quality / 10.0)

# ==================================================
# Memoization Helpers
# ==================================================
def memoize_bam_search(func):
    cache = {}

    @wraps(func)
    def wrap(*args):
        # skip the fragment_dict input to gather_reads
        (sample, chrom, lpos, rpos) = args
        inputs = (sample.name, chrom, lpos, rpos)
        if inputs not in cache:
            sys.stderr.write('invoking bam search fn -- {} -- '.format(inputs))
            cache[inputs] = func(*args)
            sys.stderr.write('length: {} \n'.format(len(cache[inputs])))
        else:
            sys.stderr.write('using cached results -- {}\n'.format(inputs))
        return cache[inputs]

    return wrap


@memoize_bam_search
def fetch_reads_from_bam(sample, chrom, left_pos, right_pos):
    return list(sample.bam.fetch(chrom, left_pos, right_pos))


# ==================================================
# logging
# ==================================================

def logit(msg):
    ts = time.strftime("[ %Y-%m-%d %T ]", datetime.datetime.now().timetuple())
    fullmsg = "{} {}".format(ts, msg)
    print(fullmsg, file=sys.stderr)
    sys.stderr.flush()

# ==================================================
# temporary directory handling
# ==================================================

def is_lsf_job():
    rv = True if 'LSB_JOBID' in os.environ else False
    return rv

@contextlib.contextmanager
def cd(newdir, cleanup=lambda: True):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)
        cleanup()

@contextlib.contextmanager
def tempdir():
    root = '/tmp'
    if is_lsf_job():
        root = '/tmp/{}.tmpdir'.format(os.environ['LSB_JOBID'])

    dir_prefix = 'svtyper-{}-'.format(os.getpid())
    dirpath = tempfile.mkdtemp(dir=root, prefix=dir_prefix)

    def cleanup():
        shutil.rmtree(dirpath)

    with cd(dirpath, cleanup):
        yield dirpath


# ==================================================
# vcf helpers
# ==================================================
def vcf_variants(vcf_file):
    with open(vcf_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                yield line

def vcf_headers(vcf_file):
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                yield line
            else:
                break

def vcf_samples(vcf_file):
    samples = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#CHROM'):
                samples = line.rstrip().split('\t')[9:]
                break
    return samples

# ==================================================
# system helpers
# ==================================================
def die(msg):
    sys.exit(msg)

import sys
from functools import wraps

from svtyper.parsers import SamFragment

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
    return 1 - 10 ** (-read.mapq / 10.0)

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

    reads = fetch_reads_from_bam(
        sample,
        chrom,
        max(pos + ci[0] - fetch_flank, 0),
        min(pos + ci[1] + fetch_flank, chrom_length)
    )

    for i, read in enumerate(reads):
        if read.is_unmapped or read.is_duplicate:
            continue

        lib = sample.get_lib(read.get_tag('RG'))
        if lib not in sample.active_libs:
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

@memoize_bam_search
def fetch_reads_from_bam(sample, chrom, left_pos, right_pos):
    return list(sample.bam.fetch(chrom, left_pos, right_pos))

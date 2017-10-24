from __future__ import print_function
import json, sys, os, math, argparse
from itertools import chain

import svtyper.version
from svtyper.parsers import Vcf, Variant, Sample, LiteRead, SamFragment
from svtyper.utils import die, logit, prob_mapq, write_sample_json, tempdir, vcf_headers, vcf_variants, vcf_samples
from svtyper.statistics import bayes_gt

import pysam


def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="\
svtyper\n\
author: " + 'Indraniel Das (idas@wustl.edu)' + "\n\
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
    parser.add_argument('--max_reads', metavar='INT', type=int, default=1000, required=False, help='maximum number of reads to assess at any variant (reduces processing time in high-depth regions, default: 1000)')
    parser.add_argument('--split_weight', metavar='FLOAT', type=float, required=False, default=1, help='weight for split reads [1]')
    parser.add_argument('--disc_weight', metavar='FLOAT', type=float, required=False, default=1, help='weight for discordant paired-end reads [1]')
    parser.add_argument('--debug', action='store_true', help=argparse.SUPPRESS)

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.input_vcf == None:
        if not sys.stdin.isatty():
            args.input_vcf = sys.stdin

    # send back the user input
    return args

def ensure_valid_alignment_file(afile):
    if not (afile.endswith('.bam') or afile.endswith('.cram')):
        die('Error: %s is not a valid alignment file (*.bam or *.cram)\n' % afile)

def open_alignment_file(afile, reference_fasta):
    fd = None
    if afile.endswith('.bam'):
        fd = pysam.AlignmentFile(afile, mode='rb')
    elif afile.endswith('.cram'):
        fd = pysam.AlignmentFile(afile, mode='rc', reference_filename=reference_fasta)
    else:
        die('Error: %s is not a valid alignment file (*.bam or *.cram)\n' % afile)

    return fd

def setup_sample(bam, lib_info_path, reference_fasta, sampling_number, min_aligned):
    fd = open_alignment_file(bam, reference_fasta)

    # only consider libraries that constitute at least this fraction of the BAM
    min_lib_prevalence = 1e-3

    sample = None
    if (lib_info_path is not None) and os.path.exists(lib_info_path):
        # use pre-calculated library metrics
        logit('Reading library metrics from %s...' % lib_info_path)
        with open(lib_info_path, 'r') as f:
            lib_info = json.load(f)
            sample = Sample.from_lib_info(fd, lib_info, min_lib_prevalence)
    else:
        # calculate and include library metrics from bam/cram
        sample = Sample.from_bam(fd, sampling_number, min_lib_prevalence)

    sample.set_exp_seq_depth(min_aligned)
    sample.set_exp_spanning_depth(min_aligned)

    return sample

def dump_library_metrics(lib_info_path, sample):
    sample_list = [sample]
    if (lib_info_path is not None) and (not os.path.exists(lib_info_path)):
        logit('Writing library metrics to %s...' % lib_info_path)
        lib_info_file = open(lib_info_path, 'w')
        write_sample_json(sample_list, lib_info_file)
        lib_info_file.close()
        logit('Finished writing library metrics')

def setup_src_vcf_file(fobj, invcf, rootdir):
    src_vcf = invcf
    if invcf == '<stdin>':
        src_vcf = dump_piped_vcf_to_file(fobj, rootdir)
    return src_vcf

def dump_piped_vcf_to_file(stdin, basedir):
    vcf = os.path.join(basedir, 'input.vcf')
    logit('dumping vcf inputs into a temporary file: {}'.format(vcf))
    line_count = 0
    with open(vcf, 'w') as f:
        for line in stdin:
            print(line, sep='', file=f)
            line_count += 1
    logit('finished temporary vcf dump -- {} lines'.format(line_count))
    return vcf

def init_vcf(vcffile, sample, scratchdir):
    v = Vcf()
    v.filename = vcffile
    hdrs = list(vcf_headers(vcffile))
    v.add_header(hdrs)
    v.add_custom_svtyper_headers()
    vcf_samples_list = vcf_samples(vcffile)
    if sample.name not in vcf_samples_list:
        fname = '<stdin>' if scratchdir in v.filename else v.filename
        msg = ("Note: Did not find sample name : '{}' "
               "in input vcf: '{}' -- adding").format(sample.name, fname)
        logit(msg)
    v.add_sample(sample.name)
    return v

def get_breakpoint_regions(breakpoint, sample, z):
    # the distance to the left and right of the breakpoint to scan
    # (max of mean + z standard devs over all of a sample's libraries)
    fetch_flank = sample.get_fetch_flank(z)

    regions = []
    for side in ('A', 'B'):
        pos = breakpoint[side]['pos']
        ci = breakpoint[side]['ci']
        chrom = breakpoint[side]['chrom']
        chrom_length = sample.bam.lengths[sample.bam.gettid(chrom)]

        left_pos = int(max(pos + ci[0] - fetch_flank, 0))
        right_pos = int(min(pos + ci[1] + fetch_flank, chrom_length))

        regions.append((sample.name, chrom, pos, left_pos, right_pos))

    return regions

def count_reads_in_region(region, sample):
    (sample_name, chrom, pos, left_pos, right_pos) = region
    count = sample.bam.count(chrom, start=left_pos, stop=right_pos, read_callback='all')
    return count

def get_reads_iterator(region, sample):
    (sample_name, chrom, pos, left_pos, right_pos) = region
    iterator = sample.bam.fetch(chrom, left_pos, right_pos)
    return iterator

def retrieve_reads_from_db(breakpoint, sample, z, max_reads):
    (over_threshold, reads) = (False, [])
    (regionA, regionB) = get_breakpoint_regions(breakpoint, sample, z)
    (countA, countB) = ( count_reads_in_region(regionA, sample), count_reads_in_region(regionB, sample) )
    if countA > max_reads or countB > max_reads:
        msg = ("SKIPPING -- Variant '{}' has too many reads\n"
                "\t\t A: {} : {}\n"
                "\t\t B: {} : {}").format(breakpoint['id'], regionA, countA, regionB, countB)
        logit(msg)
        return (overthreshold, reads)

    reads_generator = chain(
        get_reads_iterator(regionA, sample),
        get_reads_iterator(regionB, sample)
    )
    return (over_threshold, reads_generator)

def gather_reads(breakpoint, sample, z, max_reads):
    fragment_dict = {}
    (over_threshold, reads) = retrieve_reads_from_db(breakpoint, sample, z, max_reads)

    for read in reads:
        if read.query_name in fragment_dict:
            fragment_dict[read.query_name].add_read(read)
        else:
            lib = sample.get_lib(read.get_tag('RG'))
            fragment_dict[read.query_name] = SamFragment(read, lib)
    return (fragment_dict, over_threshold)

def make_empty_genotype(variant, sample):
    variant.genotype(sample.name).set_format('GT', './.')
    return variant

def make_detailed_empty_genotype(variant, sample):
    variant.genotype(sample.name).set_format('GT', './.')
    variant.qual = 0
    variant.genotype(sample.name).set_format('GQ', '.')
    variant.genotype(sample.name).set_format('SQ', '.')
    variant.genotype(sample.name).set_format('GL', '.')
    variant.genotype(sample.name).set_format('DP', 0)
    variant.genotype(sample.name).set_format('AO', 0)
    variant.genotype(sample.name).set_format('RO', 0)
    # if detailed:
    variant.genotype(sample.name).set_format('AS', 0)
    variant.genotype(sample.name).set_format('ASC', 0)
    variant.genotype(sample.name).set_format('RS', 0)
    variant.genotype(sample.name).set_format('AP', 0)
    variant.genotype(sample.name).set_format('RP', 0)
    variant.genotype(sample.name).set_format('QR', 0)
    variant.genotype(sample.name).set_format('QA', 0)
    variant.genotype(sample.name).set_format('AB', '.')
    return variant

def check_split_read_evidence(sam_fragment, breakpoint, split_slop, min_aligned):
    (ref_seq, alt_seq, alt_clip) = (0, 0, 0)

    elems = ('chrom', 'pos', 'ci', 'is_reverse')
    (chromA, posA, ciA, o1_is_reverse) = tuple(breakpoint['A'][i] for i in elems)
    (chromB, posB, ciB, o2_is_reverse) = tuple(breakpoint['B'][i] for i in elems)

    # get reference sequences
    for read in sam_fragment.primary_reads:
        is_ref_seq_A = sam_fragment.is_ref_seq(read, None, chromA, posA, ciA, min_aligned)
        is_ref_seq_B = sam_fragment.is_ref_seq(read, None, chromB, posB, ciB, min_aligned)
        if (is_ref_seq_A or is_ref_seq_B):
            p_reference = prob_mapq(read)
            ref_seq += p_reference

    # get non-reference split-read support
    for split in sam_fragment.split_reads:

        svtype = breakpoint['svtype']
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

    return (ref_seq, alt_seq, alt_clip)

def check_paired_end_evidence(fragment, breakpoint, min_aligned):
    (ref_span, alt_span) = (0, 0)
    ref_ciA = [0,0]
    ref_ciB = [0,0]

    elems = ('chrom', 'pos', 'ci', 'is_reverse')
    (chromA, posA, ciA, o1_is_reverse) = tuple(breakpoint['A'][i] for i in elems)
    (chromB, posB, ciB, o2_is_reverse) = tuple(breakpoint['B'][i] for i in elems)
    svtype = breakpoint['svtype']

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
            var_length = breakpoint['var_length']
            p_conc = fragment.p_concordant(var_length)
            if p_conc is not None:
                p_alt = (1 - p_conc) * prob_mapq(fragment.readA) * prob_mapq(fragment.readB)
                alt_span += p_alt

                # # since an alt straddler is by definition also a reference straddler,
                # # we can bail out early here to save some time
                # p_reference = p_conc * prob_mapq(fragment.readA) * prob_mapq(fragment.readB)
                # ref_span += p_reference
                # continue

        else:
            p_alt = prob_mapq(fragment.readA) * prob_mapq(fragment.readB)
            alt_span += p_alt

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
            var_length = breakpoint.get('var_length', None)
            p_conc = fragment.p_concordant(var_length)
            if p_conc is not None:
                p_reference = p_conc * prob_mapq(fragment.readA) * prob_mapq(fragment.readB)
                ref_span += (ref_straddle_A + ref_straddle_B) * p_reference / 2

    return (ref_span, alt_span, ref_ciA, ref_ciB)

def tally_variant_read_fragments(variant, sample, split_slop, min_aligned, breakpoint, sam_fragments, debug):
    # initialize counts to zero
    ref_span, alt_span = 0, 0
    ref_seq, alt_seq = 0, 0
    alt_clip = 0

    ref_ciA = [0,0]
    ref_ciB = [0,0]

    for query_name in sorted(sam_fragments.keys()):
        fragment = sam_fragments[query_name]

        (ref_seq_calc, alt_seq_calc, alt_clip_calc) = \
                check_split_read_evidence(fragment, breakpoint, split_slop, min_aligned)

        ref_seq += ref_seq_calc
        alt_seq += alt_seq_calc
        alt_clip += alt_clip_calc

        (ref_span_calc, alt_span_calc, ref_ciA_calc, ref_ciB_calc) = \
                check_paired_end_evidence(fragment, breakpoint, min_aligned)

        ref_span += ref_span_calc
        alt_span += alt_span_calc
        ref_ciA = [ x + y for x,y in zip(ref_ciA, ref_ciA_calc)]
        ref_ciB = [ x + y for x,y in zip(ref_ciB, ref_ciB_calc)]

    # in the absence of evidence for a particular type, ignore the reference
    # support for that type as well
    if (alt_seq + alt_clip) < 0.5 and alt_span >= 1:
        alt_seq = 0
        alt_clip = 0
        ref_seq = 0
    if alt_span < 0.5 and (alt_seq + alt_clip) >= 1:
        alt_span = 0
        ref_span = 0

    counts = { 'ref_seq' : ref_seq, 'alt_seq' : alt_seq,
               'ref_span' : ref_span, 'alt_span' : alt_span,
               'alt_clip' : alt_clip }

    if debug:
        items = ('ref_span', 'alt_span', 'ref_seq', 'alt_seq', 'alt_clip')
        cmsg = "\n".join(['{}: {}'.format(i, counts[i]) for i in items])
        logit("{} -- read fragment tally counts:\n{}".format(variant.var_id, cmsg))

    return counts

def bayesian_genotype(variant, sample, counts, split_weight, disc_weight, debug):
    is_dup = True if variant.get_svtype() == 'DUP' else False

    elems = ('ref_seq', 'alt_seq', 'alt_clip', 'ref_span', 'alt_span')
    (ref_seq, alt_seq, alt_clip, ref_span, alt_span) = \
        tuple(counts[i] for i in elems)

    # pre-calculations
    alt_splitters = alt_seq + alt_clip
    QR = int(split_weight * ref_seq) + int(disc_weight * ref_span)
    QA = int(split_weight * alt_splitters) + int(disc_weight * alt_span)

    # the actual bayesian calculation and decision
    gt_lplist = bayes_gt(QR, QA, is_dup)
    gt_idx = gt_lplist.index(max(gt_lplist))

    # print log probabilities of homref, het, homalt
    if debug:
        msg = ("{} -- "
               "log probabilities (homref, het, homalt) : "
               "{}").format(variant.var_id, gt_lplist)
        logit(msg)

    variant.genotype(sample.name).set_format('GL', ','.join(['%.0f' % x for x in gt_lplist]))
    variant.genotype(sample.name).set_format('DP', int(ref_seq + alt_seq + alt_clip + ref_span + alt_span))
    variant.genotype(sample.name).set_format('RO', int(ref_seq + ref_span))
    variant.genotype(sample.name).set_format('AO', int(alt_seq + alt_clip + alt_span))
    variant.genotype(sample.name).set_format('QR', QR)
    variant.genotype(sample.name).set_format('QA', QA)
    # if detailed:
    variant.genotype(sample.name).set_format('RS', int(ref_seq))
    variant.genotype(sample.name).set_format('AS', int(alt_seq))
    variant.genotype(sample.name).set_format('ASC', int(alt_clip))
    variant.genotype(sample.name).set_format('RP', int(ref_span))
    variant.genotype(sample.name).set_format('AP', int(alt_span))
    try:
        variant.genotype(sample.name).set_format('AB', '%.2g' % (QA / float(QR + QA)))
    except ZeroDivisionError:
        variant.genotype(sample.name).set_format('AB', '.')

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
        if 1 - (10**gt_lplist[gt_idx] / 10**gt_sum_log) == 0:
            phred_gq = 200
        else:
            phred_gq = abs(-10 * math.log(1 - (10**gt_lplist[gt_idx] / 10**gt_sum_log), 10))
        variant.genotype(sample.name).set_format('GQ', int(phred_gq))
        variant.genotype(sample.name).set_format('SQ', sample_qual)
        variant.qual += sample_qual
        if gt_idx == 1:
            variant.genotype(sample.name).set_format('GT', '0/1')
        elif gt_idx == 2:
            variant.genotype(sample.name).set_format('GT', '1/1')
        elif gt_idx == 0:
            variant.genotype(sample.name).set_format('GT', '0/0')
    else:
        variant.genotype(sample.name).set_format('GQ', '.')
        variant.genotype(sample.name).set_format('SQ', '.')
        variant.genotype(sample.name).set_format('GT', './.')
    
    return variant

def calculate_genotype(variant, sample, z, split_slop, min_aligned, split_weight, disc_weight, breakpoint, max_reads, debug):
    (read_batches, many) = gather_reads(breakpoint, sample, z, max_reads)

    # if there are too many reads around the breakpoint
    if many is True:
        return make_empty_genotype(variant, sample)

    # if there are no reads around the breakpoint
    if bool(read_batches) is False:
        return make_detailed_empty_genotype(variant, sample)

    counts = tally_variant_read_fragments(
        variant,
        sample,
        split_slop,
        min_aligned,
        breakpoint,
        read_batches,
        debug
    )

    total = sum([ counts[k] for k in counts.keys() ])
    if total == 0:
        return make_detailed_empty_genotype(variant, sample)

    variant = bayesian_genotype(variant, sample, counts, split_weight, disc_weight, debug)
    return variant

def genotype_vcf(src_vcf, out_vcf, sample, z, split_slop, min_aligned, sum_quals, split_weight, disc_weight, max_reads, debug):
    # initializations
    bnd_cache = {}
    src_vcf.write_header(out_vcf)
    total_variants = len(list(vcf_variants(src_vcf.filename)))

    for i, vline in enumerate(vcf_variants(src_vcf.filename)):
        v = vline.rstrip().split('\t')
        variant = Variant(v, src_vcf)
        if i % 1000 == 0:
            logit("[ {} | {} ] Processing variant {}".format(i, total_variants, variant.var_id))
        if not sum_quals:
            variant.qual = 0

        if not variant.is_svtype():
            msg = ('Warning: SVTYPE missing '
                   'at variant %s. '
                   'Skipping.\n') % (variant.var_id)
            logit(msg)
            variant.write(out_vcf)
            continue

        if not variant.is_valid_svtype():
            msg = ('Warning: Unsupported SVTYPE '
                   'at variant %s (%s). '
                   'Skipping.\n') % (variant.var_id, variant.get_svtype())
            logit(msg)
            variant.write(out_vcf)
            continue

        breakpoints = src_vcf.get_variant_breakpoints(variant)

        # special BND processing
        if variant.get_svtype() == 'BND':
            if variant.info['MATEID'] in bnd_cache:
                variant2 = variant
                variant = bnd_cache[variant.info['MATEID']]
                del bnd_cache[variant.var_id]
            else:
                bnd_cache[variant.var_id] = variant
                continue

        if breakpoints is None:
            msg = ("Found no breakpoints for variant "
                   "'{}' ({})").format(variant.var_id, variant.get_svtype())
            logit(msg)
            continue

        variant = calculate_genotype(
                variant,
                sample,
                z,
                split_slop,
                min_aligned,
                split_weight,
                disc_weight,
                breakpoints,
                max_reads,
                debug
        )
        variant.write(out_vcf)

        # special BND processing
        if variant.get_svtype() == 'BND':
            variant2.qual = variant.qual
            variant2.active_formats = variant.active_formats
            variant2.genotype = variant.genotype
            variant2.write(out_vcf)

def sso_genotype(bam_string,
                 vcf_in,
                 vcf_out,
                 min_aligned,
                 split_weight,
                 disc_weight,
                 num_samp,
                 lib_info_path,
                 debug,
                 ref_fasta,
                 sum_quals,
                 max_reads):

    # quit early if input VCF is absent
    if vcf_in is None:
        return

    (invcf, outvcf) = (os.path.abspath(vcf_in.name), os.path.abspath(vcf_out.name))
    ensure_valid_alignment_file(bam_string)

    sample = setup_sample(bam_string, lib_info_path, ref_fasta, num_samp, min_aligned)
    dump_library_metrics(lib_info_path, sample)

    # set variables for genotyping
    z = 3
    split_slop = 3 # amount of slop around breakpoint to count splitters

    with tempdir() as scratchdir:
        logit("Temporary scratch directory: {}".format(scratchdir))

        # dump the vcf file into the tmp directory, if we're reading from stdin
        src_vcf_file = setup_src_vcf_file(vcf_in, invcf, scratchdir)

        # create the vcf object
        src_vcf = init_vcf(src_vcf_file, sample, scratchdir)

        # pass through input vcf -- perform actual genotyping
        logit("Genotyping Input VCF")
        genotype_vcf(src_vcf, vcf_out, sample, z, split_slop, min_aligned, sum_quals, split_weight, disc_weight, max_reads, debug)

    sample.close()

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    if args.split_bam is not None:
        sys.stderr.write('Warning: --split_bam (-S) is deprecated. Ignoring %s.\n' % args.split_bam)

    # call primary function
    sso_genotype(args.bam,
                 args.input_vcf,
                 args.output_vcf,
                 args.min_aligned,
                 args.split_weight,
                 args.disc_weight,
                 args.num_samp,
                 args.lib_info_path,
                 args.debug,
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

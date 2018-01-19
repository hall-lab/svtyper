from __future__ import print_function
import json, sys, os, math, argparse, time
import multiprocessing as mp

from cytoolz.itertoolz import partition_all

import svtyper.version
from svtyper.parsers import Vcf, Variant, Sample, SamFragment
from svtyper.utils import die, logit, prob_mapq, write_sample_json, tempdir, vcf_headers, vcf_variants, vcf_samples
from svtyper.statistics import bayes_gt

import pysam


def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="\
svtyper\n\
author: " + 'Indraniel Das (idas@wustl.edu)' + "\n\
version: " + svtyper.version.__version__ + "\n\
description: Compute genotype of structural variants based on breakpoint depth on a SINGLE sample")
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
    parser.add_argument('--cores', type=int, metavar='INT', required=False, default=None, help='number of workers to use for parallel processing')
    parser.add_argument('--batch_size', type=int, metavar='INT', required=False, default=1000, help='number of breakpoints to batch for a parallel processing worker job')

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
    if os.path.basename(invcf) == '<stdin>':
        src_vcf = dump_piped_vcf_to_file(fobj, rootdir)
    return src_vcf

def dump_piped_vcf_to_file(stdin, basedir):
    vcf = os.path.join(basedir, 'input.vcf')
    logit('dumping vcf inputs into a temporary file: {}'.format(vcf))
    line_count = 0
    with open(vcf, 'w') as f:
        for line in stdin:
            print(line, end='', file=f)
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

def collect_breakpoints(vcf):
    breakpoints = []
    for vline in vcf_variants(vcf.filename):
        v = vline.rstrip().split('\t')
        variant = Variant(v, vcf)
        if not variant.has_svtype(): continue
        if not variant.is_valid_svtype(): continue
        brkpts = vcf.get_variant_breakpoints(variant)
        if brkpts is None: continue
        breakpoints.append(brkpts)
    return breakpoints

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

    return tuple(regions)

def count_reads_in_region(region, bam):
    (sample_name, chrom, pos, left_pos, right_pos) = region
    count = bam.count(chrom, start=left_pos, stop=right_pos, read_callback='all')
    return count

def get_reads_iterator(region, bam):
    (sample_name, chrom, pos, left_pos, right_pos) = region
    iterator = bam.fetch(chrom, start=left_pos, stop=right_pos)
    return iterator

def is_over_threshold(bam, variant_id, regions, max_reads):
    over_threshold = False
    (regionA, regionB) = regions
    (countA, countB) = ( count_reads_in_region(regionA, bam), count_reads_in_region(regionB, bam) )
    if countA > max_reads or countB > max_reads:
        over_threshold = True
        msg = ("SKIPPING -- Variant '{}' has a region with too many reads (> {})\n"
                "\t\t A: (sample={} chrom={} center={} leftflank={} rightflank={}) : {}\n"
                "\t\t B: (sample={} chrom={} center={} leftflank={} rightflank={}) : {}").format(
                        variant_id,
                        max_reads,
                        regionA[0], regionA[1], regionA[2], regionA[3], regionA[4],
                        countA,
                        regionB[0], regionB[1], regionB[2], regionB[3], regionB[4],
                        countB,
                )
        logit(msg)
    return over_threshold

def gather_reads(bam, variant_id, regions, library_data, active_libs, max_reads):
    fragment_dict = {}
    over_threshold = is_over_threshold(bam, variant_id, regions, max_reads)

    if over_threshold:
        return (fragment_dict, over_threshold)

    for region in regions:
        for read in get_reads_iterator(region, bam):
            if read.is_unmapped or read.is_duplicate: continue
            lib = library_data[read.get_tag('RG')]
            if lib.name not in active_libs: continue
            if read.query_name in fragment_dict:
                fragment_dict[read.query_name].add_read(read)
            else:
                lib = library_data[read.get_tag('RG')]
                fragment_dict[read.query_name] = SamFragment(read, lib)

    return (fragment_dict, over_threshold)

def blank_genotype_result():
    return {
        'qual' : 0,
        'formats': {
            'GT'  : './.',
            'GQ'  : '.',
            'SQ'  : '.',
            'GL'  : '.',
            'DP'  : 0,
            'AO'  : 0,
            'RO'  : 0,
            'AS'  : 0,
            'ASC' : 0,
            'RS'  : 0,
            'AP'  : 0,
            'RP'  : 0,
            'QR'  : 0,
            'QA'  : 0,
            'AB'  : '.',
        }
    }

def make_empty_genotype_result(variant_id, sample_name):
    gt = blank_genotype_result()
    gt['DP'] = '.'
    return {
        'variant.id' : variant_id,
        'sample.name' : sample_name,
        'genotype' : gt
    }

def make_detailed_empty_genotype_result(variant_id, sample_name):
    return {
        'variant.id' : variant_id,
        'sample.name' : sample_name,
        'genotype' : blank_genotype_result()
    }


def gather_split_read_evidence(sam_fragment, breakpoint, split_slop, min_aligned):
    (ref_seq, alt_seq, alt_clip) = (0, 0, 0)

    elems = ('chrom', 'pos', 'ci', 'is_reverse')
    (chromA, posA, ciA, o1_is_reverse) = [ breakpoint['A'][i] for i in elems ]
    (chromB, posB, ciB, o2_is_reverse) = [ breakpoint['B'][i] for i in elems ]

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

def gather_paired_end_evidence(fragment, breakpoint, min_aligned):
    (ref_span, alt_span) = (0, 0)
    ref_ciA = [0,0]
    ref_ciB = [0,0]

    elems = ('chrom', 'pos', 'ci', 'is_reverse')
    (chromA, posA, ciA, o1_is_reverse) = [ breakpoint['A'][i] for i in elems ]
    (chromB, posB, ciB, o2_is_reverse) = [ breakpoint['B'][i] for i in elems ]
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

def tally_variant_read_fragments(split_slop, min_aligned, breakpoint, sam_fragments, debug):
    # initialize counts to zero
    ref_span, alt_span = 0, 0
    ref_seq, alt_seq = 0, 0
    alt_clip = 0

    ref_ciA = [0,0]
    ref_ciB = [0,0]

    for query_name in sorted(sam_fragments.keys()):
        fragment = sam_fragments[query_name]

        (ref_seq_calc, alt_seq_calc, alt_clip_calc) = \
                gather_split_read_evidence(fragment, breakpoint, split_slop, min_aligned)

        ref_seq += ref_seq_calc
        alt_seq += alt_seq_calc
        alt_clip += alt_clip_calc

        (ref_span_calc, alt_span_calc, ref_ciA_calc, ref_ciB_calc) = \
                gather_paired_end_evidence(fragment, breakpoint, min_aligned)

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
        logit("{} -- read fragment tally counts:\n{}".format(breakpoint['id'], cmsg))

    return counts

def bayesian_genotype(breakpoint, counts, split_weight, disc_weight, debug):
    is_dup = breakpoint['svtype'] == 'DUP'

    elems = ('ref_seq', 'alt_seq', 'alt_clip', 'ref_span', 'alt_span')
    (ref_seq, alt_seq, alt_clip, ref_span, alt_span) = \
        [counts[i] for i in elems]

    # pre-calculations
    alt_splitters = alt_seq + alt_clip
    QR = int(split_weight * ref_seq) + int(disc_weight * ref_span)
    QA = int(split_weight * alt_splitters) + int(disc_weight * alt_span)

    # the actual bayesian calculation and decision
    gt_lplist = bayes_gt(QR, QA, is_dup)
    best, second_best = sorted([ (i, e) for i, e in enumerate(gt_lplist) ], key=lambda(x): x[1], reverse=True)[0:2]
    gt_idx = best[0]

    # print log probabilities of homref, het, homalt
    if debug:
        msg = ("{} -- "
               "log probabilities (homref, het, homalt) : "
               "{}").format(breakpoint['id'], gt_lplist)
        logit(msg)

    result = blank_genotype_result()
    result['formats']['GL'] = ','.join(['%.0f' % x for x in gt_lplist])
    result['formats']['DP'] = int(ref_seq + alt_seq + alt_clip + ref_span + alt_span)
    result['formats']['RO'] = int(ref_seq + ref_span)
    result['formats']['AO'] = int(alt_seq + alt_clip + alt_span)
    result['formats']['QR'] = QR
    result['formats']['QA'] = QA
    # if detailed:
    result['formats']['RS'] = int(ref_seq)
    result['formats']['AS'] = int(alt_seq)
    result['formats']['ASC'] = int(alt_clip)
    result['formats']['RP'] = int(ref_span)
    result['formats']['AP'] = int(alt_span)
    try:
        result['formats']['AB'] = '%.2g' % (QA / float(QR + QA))
    except ZeroDivisionError:
        result['formats']['AB'] = '.'

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
        result['formats']['GQ'] = int(phred_gq)
        result['formats']['SQ'] = sample_qual
        result['qual'] += sample_qual
        if gt_idx == 1:
            result['formats']['GT'] = '0/1'
        elif gt_idx == 2:
            result['formats']['GT'] = '1/1'
        elif gt_idx == 0:
            result['formats']['GT'] = '0/0'
    else:
        result['formats']['GQ'] = '.'
        result['formats']['SQ'] = '.'
        result['formats']['GT'] = './.'
    
    return result

def serial_calculate_genotype(bam, regions, library_data, active_libs, sample_name, split_slop, min_aligned, split_weight, disc_weight, breakpoint, max_reads, debug):
    (read_batches, many) = gather_reads(bam, breakpoint['id'], regions, library_data, active_libs, max_reads)

    # if there are too many reads around the breakpoint
    if many is True:
        return make_empty_genotype_result(breakpoint['id'], sample_name)

    # if there are no reads around the breakpoint
    if bool(read_batches) is False:
        return make_detailed_empty_genotype_result(breakpoint['id'], sample_name)

    counts = tally_variant_read_fragments(
        split_slop,
        min_aligned,
        breakpoint,
        read_batches,
        debug
    )

    total = sum([ counts[k] for k in counts.keys() ])
    if total == 0:
        return make_detailed_empty_genotype_result(breakpoint['id'], sample_name)

    result = bayesian_genotype(breakpoint, counts, split_weight, disc_weight, debug)
    return { 'variant.id' : breakpoint['id'], 'sample.name' : sample_name, 'genotype' : result }

def parallel_calculate_genotype(alignment_file, reference_fasta, library_data, active_libs, sample_name, split_slop, min_aligned, split_weight, disc_weight, max_reads, debug, batch_breakpoints, batch_regions, batch_number):
    logit("Starting batch: {}".format(batch_number))
    bam = open_alignment_file(alignment_file, reference_fasta)

    genotype_results = []
    (skip_count, no_read_count) = (0, 0)
    t0 = time.time()
    for breakpoint, regions in zip(batch_breakpoints, batch_regions):
        (read_batches, many) = gather_reads(bam, breakpoint['id'], regions, library_data, active_libs, max_reads)

        # if there are too many reads around the breakpoint
        if many is True:
            skip_count += 1
            genotype_results.append(make_empty_genotype_result(breakpoint['id'], sample_name))
            continue

        # if there are no reads around the breakpoint
        if bool(read_batches) is False:
            no_read_count += 1
            genotype_results.append(make_detailed_empty_genotype_result(breakpoint['id'], sample_name))
            continue

        counts = tally_variant_read_fragments(
            split_slop,
            min_aligned,
            breakpoint,
            read_batches,
            debug
        )

        total = sum([ counts[k] for k in counts.keys() ])
        if total == 0:
            genotype_results.append(make_detailed_empty_genotype_result(breakpoint['id'], sample_name))
            continue

        result = bayesian_genotype(breakpoint, counts, split_weight, disc_weight, debug)
        genotype_results.append({ 'variant.id' : breakpoint['id'], 'sample.name' : sample_name, 'genotype' : result })

    t1 = time.time()
    logit("Batch {} Processing Elapsed Time: {:.4f} secs".format(batch_number, t1 - t0))
    bam.close()
    return { 'genotypes' : genotype_results, 'skip-count' : skip_count, 'no-read-count' : no_read_count }

def assign_genotype_to_variant(variant, sample, genotype_result):
    variant_id = genotype_result['variant.id']
    sample_name = genotype_result['sample.name']
    outcome = genotype_result['genotype']

    if (variant.var_id != variant_id) or (sample.name != sample_name):
        msg = ("Error: assign_genotype: "
               "Variant/Sample ({}/{}) to genotype result ({}/{}) "
               "mismatch!").format(variant.var_id, sample.name, variant_id, sample_name)
        die(msg)

    if bool(outcome) is False:
        variant.genotype(sample.name).set_format('GT', './.')
    else:
        variant.qual += outcome['qual']
        variant.genotype(sample.name).set_format('GT', outcome['formats']['GT'])
        variant.genotype(sample.name).set_format('GQ', outcome['formats']['GQ'])
        variant.genotype(sample.name).set_format('SQ', outcome['formats']['SQ'])
        variant.genotype(sample.name).set_format('GL', outcome['formats']['GL'])
        variant.genotype(sample.name).set_format('DP', outcome['formats']['DP'])
        variant.genotype(sample.name).set_format('AO', outcome['formats']['AO'])
        variant.genotype(sample.name).set_format('RO', outcome['formats']['RO'])
        # if detailed:
        variant.genotype(sample.name).set_format('AS',  outcome['formats']['AS'])
        variant.genotype(sample.name).set_format('ASC', outcome['formats']['ASC'])
        variant.genotype(sample.name).set_format('RS',  outcome['formats']['RS'])
        variant.genotype(sample.name).set_format('AP',  outcome['formats']['AP'])
        variant.genotype(sample.name).set_format('RP',  outcome['formats']['RP'])
        variant.genotype(sample.name).set_format('QR',  outcome['formats']['QR'])
        variant.genotype(sample.name).set_format('QA',  outcome['formats']['QA'])
        variant.genotype(sample.name).set_format('AB',  outcome['formats']['AB'])
    return variant

def genotype_serial(src_vcf, out_vcf, sample, z, split_slop, min_aligned, sum_quals, split_weight, disc_weight, max_reads, debug):
    # initializations
    bnd_cache = {}
    src_vcf.write_header(out_vcf)
    total_variants = len(list(vcf_variants(src_vcf.filename)))

    # cleanup unused library attributes
    for rg in sample.rg_to_lib:
        sample.rg_to_lib[rg].cleanup()

    for i, vline in enumerate(vcf_variants(src_vcf.filename)):
        v = vline.rstrip().split('\t')
        variant = Variant(v, src_vcf)
        if i % 1000 == 0:
            logit("[ {} | {} ] Processing variant {}".format(i, total_variants, variant.var_id))
        if not sum_quals:
            variant.qual = 0

        if not variant.has_svtype():
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

        result = serial_calculate_genotype(
                sample.bam,
                get_breakpoint_regions(breakpoints, sample, z),
                sample.rg_to_lib,
                sample.active_libs,
                sample.name,
                split_slop,
                min_aligned,
                split_weight,
                disc_weight,
                breakpoints,
                max_reads,
                debug
        )

        variant = assign_genotype_to_variant(variant, sample, result)
        variant.write(out_vcf)

        # special BND processing
        if variant.get_svtype() == 'BND':
            variant2.qual = variant.qual
            variant2.active_formats = variant.active_formats
            variant2.genotype = variant.genotype
            variant2.write(out_vcf)

def apply_genotypes_to_vcf(src_vcf, out_vcf, genotypes, sample, sum_quals):
    # initializations
    bnd_cache = {}
    src_vcf.write_header(out_vcf)
    total_variants = len(list(vcf_variants(src_vcf.filename)))

    for i, vline in enumerate(vcf_variants(src_vcf.filename)):
        v = vline.rstrip().split('\t')
        variant = Variant(v, src_vcf)
        if not sum_quals:
            variant.qual = 0

        if not variant.has_svtype():
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

        # special BND processing
        if variant.get_svtype() == 'BND':
            if variant.info['MATEID'] in bnd_cache:
                variant2 = variant
                variant = bnd_cache[variant.info['MATEID']]
                del bnd_cache[variant.var_id]
            else:
                bnd_cache[variant.var_id] = variant
                continue

        result = genotypes[variant.var_id]

        if result is None:
            msg = ("Found no genotype results for variant "
                   "'{}' ({})").format(variant.var_id, variant.get_svtype())
            logit(msg)
            raise RuntimeError(msg)

        variant = assign_genotype_to_variant(variant, sample, result)
        variant.write(out_vcf)

        # special BND processing
        if variant.get_svtype() == 'BND':
            variant2.qual = variant.qual
            variant2.active_formats = variant.active_formats
            variant2.genotype = variant.genotype
            variant2.write(out_vcf)

def genotype_parallel(src_vcf, out_vcf, sample, z, split_slop, min_aligned, sum_quals, split_weight, disc_weight, max_reads, debug, cores, breakpoint_batch_size, ref_fasta):

    # cleanup unused library attributes
    for rg in sample.rg_to_lib:
        sample.rg_to_lib[rg].cleanup()

    # 1st pass through input vcf -- collect all the relevant breakpoints
    logit("Collecting breakpoints")
    breakpoints = collect_breakpoints(src_vcf)
    logit("Number of breakpoints/SVs to process: {}".format(len(breakpoints)))
    logit("Collecting regions")
    regions = [ get_breakpoint_regions(b, sample, z) for b in breakpoints ]
    logit("Batch breakpoints into groups of {}".format(breakpoint_batch_size))
    breakpoints_batches = list(partition_all(breakpoint_batch_size, breakpoints))
    logit("Batch regions into groups of {}".format(breakpoint_batch_size))
    regions_batches = list(partition_all(breakpoint_batch_size, regions))

    if len(breakpoints_batches) != len(regions_batches):
        raise RuntimeError("Batch error: breakpoint batches ({}) != region batches ({})".format(breakpoints_batches, regions_batches))

    logit("Number of batches to parallel process: {}".format(len(breakpoints_batches)))

    std_args = (
        sample.bam.filename,
        ref_fasta,
        sample.rg_to_lib,
        sample.active_libs,
        sample.name,
        split_slop,
        min_aligned,
        split_weight,
        disc_weight,
        max_reads,
        debug
    )

    pool = mp.Pool(processes=cores)
    results = [pool.apply_async(parallel_calculate_genotype, args=std_args + (b, r, i)) for i, (b, r) in enumerate(zip(breakpoints_batches, regions_batches))]
    results = [p.get() for p in results]
    logit("Finished parallel breakpoint processing")
    logit("Merging genotype results")
    merged_genotypes = { g['variant.id'] : g for batch in results for g in batch['genotypes'] }

    total_variants_skipped = sum([ batch['skip-count'] for batch in results ])
    total_variants_with_no_reads = sum([ batch['no-read-count'] for batch in results ])

    logit("Number of variants skipped (surpassed max-reads threshold): {}".format(total_variants_skipped))
    logit("Number of variants with no reads: {}".format(total_variants_with_no_reads))

    # 2nd pass through input vcf -- apply the calculated genotypes to the variants
    logit("Applying genotype results to vcf")
    apply_genotypes_to_vcf(src_vcf, out_vcf, merged_genotypes, sample, sum_quals)
    logit("All Done!")

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
                 max_reads,
                 cores,
                 batch_size):

    # quit early if input VCF is absent
    if vcf_in is None:
        return

    invcf = os.path.abspath(vcf_in.name)
    full_bam_path = os.path.abspath(bam_string)
    ensure_valid_alignment_file(full_bam_path)

    sample = setup_sample(full_bam_path, lib_info_path, ref_fasta, num_samp, min_aligned)
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

        if cores is None:
            logit("Genotyping Input VCF (Serial Mode)")
            # pass through input vcf -- perform actual genotyping
            genotype_serial(src_vcf, vcf_out, sample, z, split_slop, min_aligned, sum_quals, split_weight, disc_weight, max_reads, debug)
        else:
            logit("Genotyping Input VCF (Parallel Mode)")

            genotype_parallel(src_vcf, vcf_out, sample, z, split_slop, min_aligned, sum_quals, split_weight, disc_weight, max_reads, debug, cores, batch_size, ref_fasta)


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
                 args.max_reads,
                 args.cores,
                 args.batch_size)

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

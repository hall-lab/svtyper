from __future__ import print_function

import time, re, json, sys
from collections import Counter

from svtyper.statistics import mean, stdev, median, upper_mad

# ==================================================
# VCF parsing tools
# ==================================================

class Vcf(object):
    def __init__(self):
        self.file_format = 'VCFv4.2'
        # self.fasta = fasta
        self.reference = ''
        self.sample_list = []
        self.info_list = []
        self.format_list = []
        self.alt_list = []
        self.add_format('GT', 1, 'String', 'Genotype')
        self.header_misc = []
        self.filename = None
        self._bnd_breakpoint_func = None

    def add_custom_svtyper_headers(self):
        self.add_info('SVTYPE', 1, 'String', 'Type of structural variant')
        self.add_format('GQ', 1, 'Integer', 'Genotype quality')
        self.add_format('SQ', 1, 'Float', 'Phred-scaled probability that this site is variant (non-reference in this sample')
        self.add_format('GL', 'G', 'Float', 'Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy')
        self.add_format('DP', 1, 'Integer', 'Read depth')
        self.add_format('RO', 1, 'Integer', 'Reference allele observation count, with partial observations recorded fractionally')
        self.add_format('AO', 'A', 'Integer', 'Alternate allele observations, with partial observations recorded fractionally')
        self.add_format('QR', 1, 'Integer', 'Sum of quality of reference observations')
        self.add_format('QA', 'A', 'Integer', 'Sum of quality of alternate observations')
        self.add_format('RS', 1, 'Integer', 'Reference allele split-read observation count, with partial observations recorded fractionally')
        self.add_format('AS', 'A', 'Integer', 'Alternate allele split-read observation count, with partial observations recorded fractionally')
        self.add_format('ASC', 'A', 'Integer', 'Alternate allele clipped-read observation count, with partial observations recorded fractionally')
        self.add_format('RP', 1, 'Integer', 'Reference allele paired-end observation count, with partial observations recorded fractionally')
        self.add_format('AP', 'A', 'Integer', 'Alternate allele paired-end observation count, with partial observations recorded fractionally')
        self.add_format('AB', 'A', 'Float', 'Allele balance, fraction of observations from alternate allele, QA/(QR+QA)')

    def add_header(self, header):
        for line in header:
            if line.split('=')[0] == '##fileformat':
                self.file_format = line.rstrip().split('=')[1]
            elif line.split('=')[0] == '##reference':
                self.reference = line.rstrip().split('=')[1]
            elif line.split('=')[0] == '##INFO':
                a = line[line.find('<')+1:line.find('>')]
                r = re.compile(r'(?:[^,\"]|\"[^\"]*\")+')
                self.add_info(*[b.split('=')[1] for b in r.findall(a)])
            elif line.split('=')[0] == '##ALT':
                a = line[line.find('<')+1:line.find('>')]
                r = re.compile(r'(?:[^,\"]|\"[^\"]*\")+')
                self.add_alt(*[b.split('=')[1] for b in r.findall(a)])
            elif line.split('=')[0] == '##FORMAT':
                a = line[line.find('<')+1:line.find('>')]
                r = re.compile(r'(?:[^,\"]|\"[^\"]*\")+')
                self.add_format(*[b.split('=')[1] for b in r.findall(a)])
            elif line[0] == '#' and line[1] != '#':
                self.sample_list = line.rstrip().split('\t')[9:]
            elif line.startswith('##fileDate='):
                pass
            else:
                self.header_misc.append(line.rstrip())

    # return the VCF header
    def get_header(self):
        header = '\n'.join(['##fileformat=' + self.file_format,
                            '##fileDate=' + time.strftime('%Y%m%d'),
                            '##reference=' + self.reference] + \
                           [i.hstring for i in self.info_list] + \
                           [a.hstring for a in self.alt_list] + \
                           [f.hstring for f in self.format_list] + \
                           self.header_misc + \
                           ['\t'.join([
                               '#CHROM',
                               'POS',
                               'ID',
                               'REF',
                               'ALT',
                               'QUAL',
                               'FILTER',
                               'INFO',
                               'FORMAT'] + \
                                      self.sample_list
                                  )])
        return header

    def write_header(self, fd=None):
        if fd is None:
            fd = sys.stdout
        print(self.get_header(), file=fd)

    def add_info(self, id, number, type, desc):
        if id not in [i.id for i in self.info_list]:
            inf = self.Info(id, number, type, desc)
            self.info_list.append(inf)

    def add_alt(self, id, desc):
        if id not in [a.id for a in self.alt_list]:
            alt = self.Alt(id, desc)
            self.alt_list.append(alt)

    def add_format(self, id, number, type, desc):
        if id not in [f.id for f in self.format_list]:
            fmt = self.Format(id, number, type, desc)
            self.format_list.append(fmt)

    def add_sample(self, name):
        self.sample_list.append(name)

    # get the VCF column index of a sample
    # NOTE: this is zero-based, like python arrays
    def sample_to_col(self, sample):
        return self.sample_list.index(sample) + 9

    @staticmethod
    def _init_bnd_breakpoint_func():
        bnd_cache = {}

        def _get_bnd_breakpoints(variant):
            if variant.info['MATEID'] in bnd_cache:
                var2 = variant
                var = bnd_cache[variant.info['MATEID']]
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

                breakpoints = {
                    'id' : var.var_id,
                    'svtype' : 'BND',
                    'A' : {'chrom': chromA, 'pos' : posA, 'ci': ciA, 'is_reverse': o1_is_reverse},
                    'B' : {'chrom': chromB, 'pos' : posB, 'ci': ciB, 'is_reverse': o2_is_reverse},
                }

                for k in ('A', 'B'):
                    if breakpoints[k]['is_reverse']:
                        breakpoints[k]['pos'] += 1

                # remove the BND from the bnd_cache
                # to free up memory
                del bnd_cache[var.var_id]

                return breakpoints
            else:
                bnd_cache[variant.var_id] = variant
                return None
         
        return _get_bnd_breakpoints

    @staticmethod
    def _default_get_breakpoints(variant):
        chromA = variant.chrom
        chromB = variant.chrom
        posA = variant.pos
        posB = int(variant.get_info('END'))
        # confidence intervals
        ciA = map(int, variant.info['CIPOS'].split(','))
        ciB = map(int, variant.info['CIEND'].split(','))
        svtype = variant.get_svtype()
        if svtype == 'DEL':
            var_length = posB - posA
            o1_is_reverse, o2_is_reverse =  False, True
        elif svtype == 'DUP':
            o1_is_reverse, o2_is_reverse =  True, False
        elif svtype == 'INV':
            o1_is_reverse, o2_is_reverse = False, False

        if svtype != 'DEL':
            breakpoints = {
                'id' : variant.var_id,
                'svtype' : svtype,
                'A' : {'chrom': chromA, 'pos' : posA, 'ci': ciA, 'is_reverse': o1_is_reverse},
                'B' : {'chrom': chromB, 'pos' : posB, 'ci': ciB, 'is_reverse': o2_is_reverse},
            }
        else:
            breakpoints = {
                'id' : variant.var_id,
                'svtype' : svtype,
                'var_length' : var_length,
                'A' : {'chrom': chromA, 'pos' : posA, 'ci': ciA, 'is_reverse': o1_is_reverse},
                'B' : {'chrom': chromB, 'pos' : posB, 'ci': ciB, 'is_reverse': o2_is_reverse},
            }

        for k in ('A', 'B'):
            if breakpoints[k]['is_reverse']:
                breakpoints[k]['pos'] += 1

        return breakpoints

    def get_variant_breakpoints(self, variant):
        if self._bnd_breakpoint_func is None:
            func = self._init_bnd_breakpoint_func()
            self._bnd_breakpoint_func = func

        breakpoints = None
        if variant.get_svtype() == 'BND':
            func = self._bnd_breakpoint_func
            breakpoints = func(variant)
        else:
            breakpoints = self._default_get_breakpoints(variant)

        return breakpoints

    class Info(object):
        def __init__(self, id, number, type, desc):
            self.id = str(id)
            self.number = str(number)
            self.type = str(type)
            self.desc = str(desc)
            # strip the double quotes around the string if present
            if self.desc.startswith('"') and self.desc.endswith('"'):
                self.desc = self.desc[1:-1]
            self.hstring = '##INFO=<ID=' + self.id + ',Number=' + self.number + ',Type=' + self.type + ',Description=\"' + self.desc + '\">'

    class Alt(object):
        def __init__(self, id, desc):
            self.id = str(id)
            self.desc = str(desc)
            # strip the double quotes around the string if present
            if self.desc.startswith('"') and self.desc.endswith('"'):
                self.desc = self.desc[1:-1]
            self.hstring = '##ALT=<ID=' + self.id + ',Description=\"' + self.desc + '\">'

    class Format(object):
        def __init__(self, id, number, type, desc):
            self.id = str(id)
            self.number = str(number)
            self.type = str(type)
            self.desc = str(desc)
            # strip the double quotes around the string if present
            if self.desc.startswith('"') and self.desc.endswith('"'):
                self.desc = self.desc[1:-1]
            self.hstring = '##FORMAT=<ID=' + self.id + ',Number=' + self.number + ',Type=' + self.type + ',Description=\"' + self.desc + '\">'

class Variant(object):
    def __init__(self, var_list, vcf):
        self.chrom = var_list[0]
        self.pos = int(var_list[1])
        self.var_id = var_list[2]
        self.ref = var_list[3]
        self.alt = var_list[4]
        if var_list[5] == '.':
            self.qual = 0
        else:
            self.qual = float(var_list[5])
        self.filter = var_list[6]
        self.sample_list = vcf.sample_list
        self.info_list = vcf.info_list
        self.info = dict()
        self.format_list = vcf.format_list
        self.active_formats = list()
        self.gts = dict()

        # fill in empty sample genotypes
        if len(var_list) < 8:
            sys.stderr.write('Error: VCF file must have at least 8 columns\n')
            exit(1)
        if len(var_list) < 9:
            var_list.append("GT")

        # make a genotype for each sample at variant
        for s in self.sample_list:
            try:
                s_gt = var_list[vcf.sample_to_col(s)].split(':')[0]
                self.gts[s] = Genotype(self, s, s_gt)
                # import the existing fmt fields
                for j in zip(var_list[8].split(':'), var_list[vcf.sample_to_col(s)].split(':')):
                    self.gts[s].set_format(j[0], j[1])
            except IndexError:
                self.gts[s] = Genotype(self, s, './.')

        self.info = dict()
        i_split = [a.split('=') for a in var_list[7].split(';')] # temp list of split info column
        for i in i_split:
            if len(i) == 1:
                i.append(True)
            self.info[i[0]] = i[1]

    def set_info(self, field, value):
        if field in [i.id for i in self.info_list]:
            self.info[field] = value
        else:
            sys.stderr.write('Error: invalid INFO field, \"' + field + '\"\n')
            exit(1)

    def get_info(self, field):
        return self.info[field]

    def get_info_string(self):
        i_list = list()
        for info_field in self.info_list:
            if info_field.id in self.info.keys():
                if info_field.type == 'Flag':
                    i_list.append(info_field.id)
                else:
                    i_list.append('%s=%s' % (info_field.id, self.info[info_field.id]))
        return ';'.join(i_list)

    def get_format_string(self):
        f_list = list()
        for f in self.format_list:
            if f.id in self.active_formats:
                f_list.append(f.id)
        return ':'.join(f_list)

    def genotype(self, sample_name):
        if sample_name in self.sample_list:
            return self.gts[sample_name]
        else:
            sys.stderr.write('Error: invalid sample name, \"' + sample_name + '\"\n')

    def get_var_string(self):
        s = '\t'.join(map(str,[
            self.chrom,
            self.pos,
            self.var_id,
            self.ref,
            self.alt,
            '%0.2f' % self.qual,
            self.filter,
            self.get_info_string(),
            self.get_format_string(),
            '\t'.join(self.genotype(s).get_gt_string() for s in self.sample_list)
        ]))
        return s

    def write(self, fd=None):
        if fd is None:
            fd = sys.stdout
        print(self.get_var_string(), file=fd)

    def has_svtype(self):
        flag = True
        try:
            svtype = self.get_info('SVTYPE')
        except KeyError:
            flag = False
        return flag

    def get_svtype(self):
        return self.get_info('SVTYPE')

    def is_valid_svtype(self):
        svtype = self.get_svtype()
        flag = svtype in ('BND', 'DEL', 'DUP', 'INV')
        return flag

class Genotype(object):
    def __init__(self, variant, sample_name, gt):
        self.format = dict()
        self.variant = variant
        self.set_format('GT', gt)

    def set_format(self, field, value):
        if field in [i.id for i in self.variant.format_list]:
            self.format[field] = value
            if field not in self.variant.active_formats:
                self.variant.active_formats.append(field)
                # sort it to be in the same order as the format_list in header
                self.variant.active_formats.sort(key=lambda x: [f.id for f in self.variant.format_list].index(x))
        else:
            sys.stderr.write('Error: invalid FORMAT field, \"' + field + '\"\n')
            exit(1)

    def get_format(self, field):
        return self.format[field]

    def get_gt_string(self):
        g_list = list()
        for f in self.variant.active_formats:
            if f in self.format:
                if type(self.format[f]) == float:
                    g_list.append('%0.2f' % self.format[f])
                else:
                    g_list.append(self.format[f])
            else:
                g_list.append('.')
        return ':'.join(map(str,g_list))

# ==================================================
# Library parsing
# ==================================================

# holds a library's insert size and read length information
class Library(object):
    def __init__(self,
                 name,
                 bam,
                 readgroups,
                 read_length,
                 hist,
                 dens,
                 mean,
                 sd,
                 prevalence,
                 num_samp):

        # parse arguments
        self.name = name
        self.bam = bam
        self.num_samp = num_samp
        self.readgroups = readgroups
        self.read_length = read_length
        self.hist = hist
        self.dens = dens
        self.mean = mean
        self.sd = sd
        self.prevalence = prevalence

        # if information is missing, compute it
        if self.read_length is None:
            self.calc_read_length()
        if self.hist is None:
            self.calc_insert_hist()
        if self.dens is None:
            self.calc_insert_density()
        if self.prevalence is None:
            self.calc_lib_prevalence()

        # remove the uneeded attribute (only needed during object initialization)
        del self.bam

    def cleanup(self):
        # the hist attribute is not really needed in the actual genotyping
        if hasattr(self, 'hist'):
            del self.hist

    @classmethod
    def from_lib_info(cls,
                      sample_name,
                      lib_index,
                      bam,
                      lib_info):

        lib = lib_info[sample_name]['libraryArray'][lib_index]

        # convert the histogram keys to integers (from strings in JSON)
        lib_hist = {int(k):int(v) for k,v in lib['histogram'].items()}

        return cls(lib['library_name'],
                   bam,
                   lib['readgroups'],
                   int(lib['read_length']),
                   lib_hist,
                   None,
                   float(lib['mean']),
                   float(lib['sd']),
                   float(lib['prevalence']),
                   0)

    @classmethod
    def from_bam(cls,
                 lib_name,
                 bam,
                 num_samp):

        # get readgroups that comprise the library
        readgroups = []
        for r in bam.header['RG']:
            try:
                in_lib = r['LB'] == lib_name
            except KeyError, e:
                in_lib = lib_name == ''

            if in_lib:
                readgroups.append(r['ID'])

        return cls(lib_name,
                   bam,
                   readgroups,
                   None,
                   None,
                   None,
                   None,
                   None,
                   None,
                   num_samp)

    # calculate the library's prevalence in the BAM file
    def calc_lib_prevalence(self):
        max_count = 100000
        lib_counter = 0
        read_counter = 0

        for read in self.bam.fetch():
            if read_counter == max_count:
                break
            if read.get_tag('RG') in self.readgroups:
                lib_counter += 1
            read_counter += 1

        self.prevalence = float(lib_counter) / read_counter

    # get read length
    def calc_read_length(self):
        max_rl = 0
        counter = 0
        num_samp = 10000
        for read in self.bam.fetch():
            if read.get_tag('RG') not in self.readgroups:
                continue
            if read.infer_query_length() > max_rl:
                max_rl = read.infer_query_length()
            if counter == num_samp:
                break
            counter += 1
        self.read_length = max_rl

    def is_primary(self, read):
        return (not read.is_supplementary and not read.is_secondary)

    # generate empirical histogram of the sample's insert size distribution
    def calc_insert_hist(self):
        counter = 0
        skip = 0
        skip_counter = 0
        mads = 10
        ins_list = []

        # Each entry in valueCounts is a value, and its count is
        # the number of instances of that value observed in the dataset.
        # So valueCount[5] is the number of times 5 has been seen in the data.
        valueCounts = Counter()
        for read in self.bam.fetch():
            if skip_counter < skip:
                skip_counter += 1
                continue
            if (read.is_reverse
                or not read.mate_is_reverse
                or read.is_unmapped
                or read.mate_is_unmapped
                or not self.is_primary(read)
                or read.template_length <= 0
                or read.get_tag('RG') not in self.readgroups):
                continue
            else:
                valueCounts[read.template_length] += 1
                counter += 1
            if counter == self.num_samp:
                break

        if len(valueCounts) == 0:
            sys.stderr.write('Error: failed to build insert size histogram for paired-end reads.\n\
Please ensure BAM file (%s) has inward facing, paired-end reads.\n' % self.bam.filename)
            exit(1)

        # remove outliers
        med = median(valueCounts)
        u_mad = upper_mad(valueCounts, med)
        for x in [x for x in list(valueCounts) if x > med + mads * u_mad]:
            del valueCounts[x]

        self.hist = valueCounts
        self.mean = mean(self.hist)
        self.sd = stdev(self.hist)

    # calculate the density curve for and insert size histogram
    def calc_insert_density(self):
        dens = Counter()
        for i in list(self.hist):
            dens[i] = float(self.hist[i])/self.countRecords(self.hist)
        self.dens = dens

    def countRecords(self, myCounter):
        numRecords = sum(myCounter.values())
        return numRecords

# ==================================================
# Sample parsing
# ==================================================

# holds each sample's BAM and library information
class Sample(object):
    # general constructor
    def __init__(self,
                 name,
                 bam,
                 num_samp,
                 lib_dict,
                 rg_to_lib,
                 min_lib_prevalence,
                 bam_mapped,
                 bam_unmapped):

        self.name = name
        self.bam = bam
        self.lib_dict = lib_dict
        self.rg_to_lib = rg_to_lib
        self.bam_mapped = bam_mapped
        self.bam_unmapped = bam_unmapped

        # get active libraries
        self.active_libs = []
        for lib in lib_dict.values():
            if lib.prevalence >= min_lib_prevalence:
                self.active_libs.append(lib.name)

    # constructor from supplied JSON descriptor
    @classmethod
    def from_lib_info(cls,
                      bam,
                      lib_info,
                      min_lib_prevalence):

        name = bam.header['RG'][0]['SM']
        num_samp = 0
        rg_to_lib = {}
        lib_dict = {}

        try:
            for i in xrange(len(lib_info[name]['libraryArray'])):
                lib = lib_info[name]['libraryArray'][i]
                lib_name = lib['library_name']

                # make library object
                lib_dict[lib_name] = Library.from_lib_info(name,
                                                           i,
                                                           bam,
                                                           lib_info)

                # make a map from readgroup IDs to library objects
                for rg in lib['readgroups']:
                    rg_to_lib[rg] = lib_dict[lib_name]
        except KeyError:
            sys.stderr.write('Error: sample %s not found in JSON library file.\n' % name)
            exit(1)

        return cls(name,
                   bam,
                   num_samp,
                   lib_dict,
                   rg_to_lib,
                   min_lib_prevalence,
                   lib_info[name]['mapped'],
                   lib_info[name]['unmapped'])

    # constructor for empirical distributions
    @classmethod
    def from_bam(cls,
                 bam,
                 num_samp,
                 min_lib_prevalence):

        name = bam.header['RG'][0]['SM']
        rg_to_lib = {}
        lib_dict = {}

        for r in bam.header['RG']:
            try:
                lib_name=r['LB']
            except KeyError, e:
                lib_name=''

            # add the new library
            if lib_name not in lib_dict:
                new_lib = Library.from_bam(lib_name, bam, num_samp)
                lib_dict[lib_name] = new_lib
            rg_to_lib[r['ID']] = lib_dict[lib_name]

        return cls(name,
                   bam,
                   num_samp,
                   lib_dict,
                   rg_to_lib,
                   min_lib_prevalence,
                   bam.mapped,
                   bam.unmapped)

    # get the maximum fetch flank for reading the BAM file
    def get_fetch_flank(self, z):
        return max([lib.mean + (lib.sd * z) for lib in self.lib_dict.values()])

    # return the library object for a specified read group
    def get_lib(self, readgroup):
        return self.rg_to_lib[readgroup]

    # return the expected spanning coverage at any given base
    def set_exp_spanning_depth(self, min_aligned):
        genome_size = float(sum(self.bam.lengths))
        weighted_mean_span = sum([(lib.mean - 2 * lib.read_length + 2 * min_aligned) * lib.prevalence for lib in self.lib_dict.values()])
        exp_spanning_depth = (weighted_mean_span * self.bam_mapped) / genome_size

        self.exp_spanning_depth = exp_spanning_depth

        return

    # return the expected sequence coverage at any given base
    def set_exp_seq_depth(self, min_aligned):
        genome_size = float(sum(self.bam.lengths))
        weighted_mean_read_length = sum([(lib.read_length - 2 * min_aligned) * lib.prevalence for lib in self.lib_dict.values()])
        exp_seq_depth = (weighted_mean_read_length * self.bam_mapped) / genome_size

        self.exp_seq_depth = exp_seq_depth

        return

    def close(self):
        self.bam.close()

# ==================================================
# Class for SAM fragment, containing all alignments
# from a single molecule
# ==================================================

class SamFragment(object):
    def __init__(self, read, lib):
        self.lib = lib
        self.primary_reads = []
        self.split_reads = []
        self.read_set = set()
        self.num_primary = 0
        self.query_name = read.query_name

        self.readA = None
        self.readB = None
        self.ispan = None
        self.ospan = None

        self.add_read(read)

    def is_primary(self, read):
        return (not read.is_supplementary and not read.is_secondary)

    def add_read(self, read):
        # ensure we don't add the same read twice
        read_hash = read.__hash__()
        if read_hash in self.read_set:
            return
        else:
            self.read_set.add(read_hash)

        if self.is_primary(read):
            # add the primary reads
            self.primary_reads.append(read)
            self.num_primary += 1

            # make split candidate and check whether it's valid
            split_candidate = SplitRead(read, self.lib)
            if split_candidate.is_valid():
                self.split_reads.append(split_candidate)

            # complete set of primaries
            if self.num_primary == 2:
                self.readA, self.readB = self.primary_reads

    # tag the read with R (ref), A (alt), or U (unknown) XV tag
    def tag_span(self, p_alt=None):
        if p_alt is None:
            value = 'U'
        elif p_alt > 0:
            value = 'A'
        else:
            value = 'R'
        for read in self.primary_reads:
            if not read.has_tag('XV'):
                read.set_tag('XV', value)

        return

    # get inner span
    def get_ispan(self, min_aligned):
        ispan1 = self.readA.reference_start + min_aligned
        ispan2 = self.readB.reference_end - min_aligned - 1

        return (ispan1, ispan2)

    # get outer span
    def get_ospan(self):
        ospan1 = self.readA.reference_start
        ospan2 = self.readB.reference_end

        return (ospan1, ospan2)

    # returns boolean of whether a single read crosses
    # a genomic reference point with appropriate flanking
    # sequence
    def is_ref_seq(self,
                   read,
                   variant,
                   chrom, pos, ci,
                   min_aligned):

        # check chromosome matching
        # Note: this step is kind of slow
        if read.reference_name != chrom:
            return False

        # ensure there is min_aligned on both sides of position
        if read.get_overlap(max(0, pos - min_aligned), pos + min_aligned) < 2 * min_aligned:
            return False

        return True


    # returns boolean of whether the primary pair of a
    # fragment straddles a genomic point
    def is_pair_straddle(self,
                         chromA, posA, ciA,
                         chromB, posB, ciB,
                         o1_is_reverse, o2_is_reverse,
                         min_aligned,
                         lib):
        if self.num_primary != 2:
            return False

        # check orientation
        if self.readA.is_reverse != o1_is_reverse:
            return False
        if self.readB.is_reverse != o2_is_reverse:
            return False

        # check chromosome matching
        # Note: this step is kind of slow
        if self.readA.reference_name != chromA:
            return False
        if self.readB.reference_name != chromB:
            return False

        # get the inner span
        ispan = self.get_ispan(min_aligned)

        # check orientations and positions
        flank = lib.mean + lib.sd * 3
        if not o1_is_reverse and (ispan[0] > posA + ciA[1] or ispan[0] < posA + ciA[0] - flank):
            return False
        if o1_is_reverse and (ispan[0] < posA + ciA[0] or ispan[0] > posA + ciA[1] + flank):
            return False
        if not o2_is_reverse and (ispan[1] > posB + ciB[1] or ispan[1] < posB + ciB[0] - flank):
            return False
        if o2_is_reverse and (ispan[1] < posB + ciB[0] or ispan[1] > posB + ciB[1] + flank):
            return False

        return True

    # calculate the probability that a read pair is concordant at a breakpoint,
    # given the putative variant size and insert distribution of the library.
    def p_concordant(self, var_length=None):
        # a priori probability that a read-pair is concordant
        disc_prior = 0.05
        conc_prior = 1 - disc_prior

        ospan = self.get_ospan()

        # outer span length
        ospan_length = abs(ospan[1] - ospan[0])

        # if no variant length (such as in the case of a non-deletion variant)
        # default to mean plus 3 stdev
        z = 3
        if var_length is None:
            var_length = self.lib.mean + self.lib.sd * z

        try:
            p = float(self.lib.dens[ospan_length]) * conc_prior / (conc_prior * self.lib.dens[ospan_length] + disc_prior * (self.lib.dens[ospan_length - var_length]))
        except ZeroDivisionError:
            p = None

        return p > 0.5

# ==================================================
# Class for a split-read, containing all alignments
# from a single chimeric read
# ==================================================

# each SplitRead object has a left and a right SplitPiece
# (reads with more than 2 split alignments are discarded)
class SplitRead(object):
    def __init__(self, read, lib):
        self.query_name = read.query_name
        self.read = read
        self.lib = lib
        self.sa = None
        self.q1 = None
        self.q2 = None
        self.is_soft_clip = False

    # the piece of each split alignemnt
    class SplitPiece(object):
        def __init__(self,
                     chrom,
                     reference_start,
                     is_reverse,
                     cigar,
                     mapq):
            self.chrom = chrom
            self.reference_start = reference_start
            self.reference_end = None
            self.is_reverse = is_reverse
            self.cigar = cigar
            self.mapping_quality = mapq
            self.left_query = None
            self.right_query = None

            # get query positions
            self.query_pos = self.get_query_pos_from_cigar(self.cigar, self.is_reverse)

        # get the positions of the query that are aligned
        @staticmethod
        def get_query_pos_from_cigar(cigar, is_reverse):
            query_start = 0
            query_end = 0
            query_length = 0

            # flip if negative strand
            if is_reverse:
                cigar = cigar[::-1]

            # iterate through cigartuple
            for i in xrange(len(cigar)):
                k, n = cigar[i]
                if k in (4,5): # H, S
                    if i == 0:
                        query_start += n
                        query_end += n
                        query_length += n
                    else:
                        query_length += n
                elif k in (0,1,7,8): # M, I, =, X
                    query_end += n
                    query_length +=n

            d = QueryPos(query_start, query_end, query_length);
            return d

        def set_reference_end(self, reference_end):
            self.reference_end = reference_end

    def is_clip_op(self, op):
        '''
        whether the CIGAR OP code represents a clipping event
        '''
        return op == 4 or op == 5

    # check if passes QC, and populate with necessary info
    def is_valid(self,
                 min_non_overlap = 20,
                 min_indel = 50,
                 max_unmapped_bases = 50):
        # check for SA tag
        if not self.read.has_tag('SA'):
            # Include soft-clipped reads that didn't generate split-read alignments
            if self.is_clip_op(self.read.cigar[0][0]) or self.is_clip_op(self.read.cigar[-1][0]):
                clip_length = max(self.read.cigar[0][1] * self.is_clip_op(self.read.cigar[0][0]), self.read.cigar[-1][1] * self.is_clip_op(self.read.cigar[-1][0]))
                # Only count if the longest clipping event is greater than the cutoff and we have mapped a reasonable number of bases
                if clip_length > 0 and (self.read.query_length - self.read.query_alignment_length) <= max_unmapped_bases:
                    a = self.SplitPiece(self.read.reference_name,
                                        self.read.reference_start,
                                        self.read.is_reverse,
                                        self.read.cigar,
                                        self.read.mapping_quality)
                    a.set_reference_end(self.read.reference_end)
                    b = self.SplitPiece(None,
                                        1,
                                        self.read.is_reverse,
                                        self.read.cigar,
                                        0)
                    b.set_reference_end(1)
                    self.set_order_by_clip(a, b)
                    self.is_soft_clip = True
                    return True
                else:
                    return False

            return False

        # parse SA tag
        sa_list = self.read.get_tag('SA').rstrip(';').split(';')
        if len(sa_list) > 1:
            return False
        else:
            self.sa = sa_list[0].split(',')
            mate_chrom = self.sa[0]
            mate_pos = int(self.sa[1]) - 1 # SA tag is one-based, while SAM is zero-based
            mate_is_reverse = self.sa[2] == '-'
            mate_cigar = self.cigarstring_to_tuple(self.sa[3])
            mate_mapq = int(self.sa[4])

        # make SplitPiece objects
        a = self.SplitPiece(self.read.reference_name,
                            self.read.reference_start,
                            self.read.is_reverse,
                            self.read.cigar,
                            self.read.mapping_quality)
        a.set_reference_end(self.read.reference_end)

        b = self.SplitPiece(mate_chrom,
                            mate_pos,
                            mate_is_reverse,
                            mate_cigar,
                            mate_mapq)
        b.set_reference_end(self.get_reference_end_from_cigar(b.reference_start, b.cigar))

        # set query_left and query_right splitter by alignment position on the reference
        # this is used for non-overlap and off-diagonal filtering
        # (query_left and query_right are random when on different chromosomes
        if self.read.reference_name == mate_chrom:
            if self.read.pos > mate_pos:
                self.query_left = b
                self.query_right = a
            else:
                self.query_left = a
                self.query_right = b
        else:
            self.set_order_by_clip(a, b)

        # check non-overlap
        if self.non_overlap() < min_non_overlap:
            return False

        # check off-diagonal distance and desert
        # only relevant when split pieces are on the same chromosome and strand
        if (self.query_left.chrom == self.query_right.chrom
            and self.query_left.is_reverse == self.query_right.is_reverse):
            # use end diagonal on left and start diagonal on right since
            # the start and end diags might be different if there is an
            # indel in the alignments
            if self.query_left.is_reverse:
                left_diag = self.get_start_diagonal(self.query_left)
                right_diag = self.get_end_diagonal(self.query_right)
                ins_size = right_diag - left_diag
            else:
                left_diag = self.get_end_diagonal(self.query_left)
                right_diag = self.get_start_diagonal(self.query_right)
                ins_size = left_diag - right_diag
            if abs(ins_size) < min_indel:
                return False

            # check for desert gap of indels
            desert = self.query_right.query_pos.query_start - self.query_left.query_pos.query_end - 1
            if desert > 0 and desert - max(0, ins_size) > max_unmapped_bases:
                return False

        # passed all checks. valid split-read
        return True

    # reference position where the alignment would have started
    # if the entire query sequence would have aligned
    @staticmethod
    def get_start_diagonal(split_piece):
        sclip = split_piece.query_pos.query_start
        if split_piece.is_reverse:
            sclip = split_piece.query_pos.query_length - split_piece.query_pos.query_end
        return split_piece.reference_start - sclip

    # reference position where the alignment would have ended
    # if the entire query sequence would have aligned
    @staticmethod
    def get_end_diagonal(split_piece):
        query_aligned = split_piece.query_pos.query_end
        if split_piece.is_reverse:
            query_aligned = split_piece.query_pos.query_length - split_piece.query_pos.query_start
        return split_piece.reference_end - query_aligned


    # adapted from Matt Shirley (http://coderscrowd.com/app/public/codes/view/171)
    @staticmethod
    def cigarstring_to_tuple(cigarstring):
        cigar_dict = {'M':0, 'I':1,'D':2,'N':3, 'S':4, 'H':5, 'P':6, '=':7, 'X':8}
        pattern = re.compile('([MIDNSHPX=])')
        values = pattern.split(cigarstring)[:-1] ## turn cigar into tuple of values
        paired = (values[n:n+2] for n in xrange(0, len(values), 2)) ## pair values by twos
        return [(cigar_dict[pair[1]], int(pair[0])) for pair in paired]

    @staticmethod
    def get_reference_end_from_cigar(reference_start, cigar):
        '''
        This returns the coordinate just past the last aligned base.
        This matches the behavior of pysam's reference_end method
        '''
        reference_end = reference_start
        
        # iterate through cigartuple
        for i in xrange(len(cigar)):
            k, n = cigar[i]
            if k in (0,2,3,7,8): # M, D, N, =, X
                reference_end += n
        return reference_end

    def non_overlap(self):
        # get overlap of aligned query positions
        overlap = self.get_query_overlap(self.query_left.query_pos.query_start,
                                         self.query_left.query_pos.query_end,
                                         self.query_right.query_pos.query_start,
                                         self.query_right.query_pos.query_end)

        # get minimum non-overlap
        left_non_overlap = 1 + self.query_left.query_pos.query_end - self.query_left.query_pos.query_start - overlap
        right_non_overlap = 1 + self.query_right.query_pos.query_end - self.query_right.query_pos.query_start - overlap
        non_overlap = min(left_non_overlap, right_non_overlap)

        return non_overlap

    def get_query_overlap(self, s1, e1, s2, e2):
        o = 1 + min(e1, e2) - max(s1, s2)
        return max(0, o)

    @staticmethod
    def check_split_support(split, chrom, pos, is_reverse, split_slop):
        if split.chrom != chrom:
            return False

        if is_reverse:
            coord = split.reference_start
        else:
            coord = split.reference_end

        if (coord > pos + split_slop
                or coord < pos - split_slop):
            return False
        return True

    def is_split_straddle(self,
                          chromA, posA, ciA,
                          chromB, posB, ciB,
                          o1_is_reverse, o2_is_reverse,
                          svtype, split_slop):

        # arrange the SV breakends from left to right
        if (chromA != chromB
            or (chromA == chromB and posA > posB)):
            chrom_left = chromB
            pos_left = posB
            ci_left = ciB
            is_reverse_left = o2_is_reverse
            chrom_right = chromA
            pos_right = posA
            ci_right = ciA
            is_reverse_right = o1_is_reverse
        else:
            chrom_left = chromA
            pos_left = posA
            ci_left = ciA
            is_reverse_left = o1_is_reverse
            chrom_right = chromB
            pos_right = posB
            ci_right = ciB
            is_reverse_right = o2_is_reverse


        # check split chromosomes against variant
        left_split = False
        right_split = False

        if (not self.is_soft_clip) or svtype == 'DEL' or svtype == 'INS':
            left_split = self.check_split_support(self.query_left,
                    chrom_left,
                    pos_left,
                    is_reverse_left,
                    split_slop)
            right_split = self.check_split_support(self.query_right,
                    chrom_right,
                    pos_right,
                    is_reverse_right,
                    split_slop)
        elif svtype == 'DUP':
            left_split = self.check_split_support(self.query_left,
                    chrom_right,
                    pos_right,
                    is_reverse_right,
                    split_slop)
            right_split = self.check_split_support(self.query_right,
                    chrom_left,
                    pos_left,
                    is_reverse_left,
                    split_slop)
        elif svtype == 'INV':
                # check all possible sides
            left_split_left = self.check_split_support(self.query_left,
                    chrom_left,
                    pos_left,
                    is_reverse_left,
                    split_slop)
            left_split_right = self.check_split_support(self.query_left,
                    chrom_right,
                    pos_right,
                    is_reverse_right,
                    split_slop)
            left_split = left_split_left or left_split_right
            right_split_left = self.check_split_support(self.query_right,
                    chrom_left,
                    pos_left,
                    is_reverse_left,
                    split_slop)
            right_split_right = self.check_split_support(self.query_right,
                    chrom_right,
                    pos_right,
                    is_reverse_right,
                    split_slop)
            right_split = right_split_left or right_split_right

        return (left_split, right_split)

    # tag the read with R (ref), A (alt), or U (unknown) XV tag
    def tag_split(self, p_alt=None):
        if p_alt is None:
            value = 'U'
        elif p_alt > 0:
            value = 'A'
        else:
            value = 'R'

        self.read.set_tag('XV', value)

        return

    def set_order_by_clip(self, a, b):
        '''
        Determine which SplitPiece is the leftmost based
        on the side of the longest clipping operation
        '''
        if self.is_left_clip(a.cigar):
            self.query_left = b
            self.query_right = a
        else:
            self.query_left = a
            self.query_right = b

    def is_left_clip(self, cigar):
        '''
        whether the left side of the read (w/ respect to reference) is clipped.
        Clipping side is determined as the side with the longest clip.
        Adjacent clipping operations are not considered
        '''
        left_tuple = cigar[0]
        right_tuple = cigar[-1]
        left_clipped = self.is_clip_op(left_tuple[0])
        right_clipped = self.is_clip_op(right_tuple[0])

        return (left_clipped and not right_clipped) or (left_clipped and right_clipped and left_tuple[1] > right_tuple[1])


# structure to hold query position information
class QueryPos (object):
    """
    struct to store the start and end positions of query CIGAR operations
    """
    def __init__(self, query_start, query_end, query_length):
        self.query_start = int(query_start)
        self.query_end = int(query_end)
        self.query_length  = int(query_length)


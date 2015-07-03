#!/usr/bin/env python

import argparse, sys, copy, gzip
import math, time, re
import numpy
from scipy import stats
from collections import Counter
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.2 $"
__date__ = "$Date: 2014-04-28 14:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
sv_classifier.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: classify structural variants")
    parser.add_argument('-i', '--input', metavar='VCF', dest='vcf_in', type=argparse.FileType('r'), default=None, help='VCF input [stdin]')
    parser.add_argument('-g', '--gender', metavar='FILE', dest='gender', type=argparse.FileType('r'), required=True, default=None, help='tab delimited file of sample genders (male=1, female=2)\nex: SAMPLE_A\t2')
    parser.add_argument('-a', '--annotation', metavar='BED', dest='ae_path', type=str, default=None, help='BED file of annotated elements')
    parser.add_argument('-f', '--fraction', metavar='FLOAT', dest='f_overlap', type=float, default=0.9, help='fraction of reciprocal overlap to apply annotation to variant [0.9]')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.vcf_in == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.vcf_in = sys.stdin
    # send back the user input
    return args

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

    # return the VCF header
    def get_header(self):
        header = '\n'.join(['##fileformat=' + self.file_format,
                            '##fileDate=' + time.strftime('%Y%m%d'),
                            '##reference=' + self.reference] + \
                           [i.hstring for i in self.info_list] + \
                           [a.hstring for a in self.alt_list] + \
                           [f.hstring for f in self.format_list] + \
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
            sys.stderr.write('\nError: VCF file must have at least 8 columns\n')
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
            sys.stderr.write('\nError: invalid INFO field, \"' + field + '\"\n')
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
            sys.stderr.write('\nError: invalid sample name, \"' + sample_name + '\"\n')

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
            sys.stderr.write('\nError: invalid FORMAT field, \"' + field + '\"\n')
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

# test whether variant has read depth support
def has_depth_support(var, gender):
    slope_threshold = 0.1
    rsquared_threshold = 0.1
    
    if 'CN' in var.active_formats:
        # allele balance list
        ab_list = []
        for s in var.sample_list:
            ab_str = var.genotype(s).get_format('AB')
            if ab_str == '.':
                ab_list.append(-1)
                continue

            ab_list.append(float(ab_str))

        # populate read-depth list, accounting for sample gender
        rd_list = []
        for s in var.sample_list:
            if (var.chrom == 'X' or var.chrom == 'Y') and gender[s] == 1:
                rd_list.append(float(var.genotype(s).get_format('CN')) * 2)
            else:
                rd_list.append(float(var.genotype(s).get_format('CN')))

        rd = numpy.array([ab_list, rd_list])

        # remove missing genotypes
        rd = rd[:, rd[0]!=-1]

        # ensure non-uniformity in genotype and read depth
        if len(numpy.unique(rd[0,:])) > 1 and len(numpy.unique(rd[1,:])) > 1:
            # calculate regression
            (slope, intercept, r_value, p_value, std_err) = stats.linregress(rd)
            # print slope, intercept, r_value, var.info['SVTYPE'], var.var_id

            # # write the scatterplot to a file
            # f = open('data/%s_%s_%sbp.txt' % (var.info['SVTYPE'], var.var_id, var.info['SVLEN']), 'w')
            # numpy.savetxt(f, numpy.transpose(rd), delimiter='\t')
            # f.close()
            
            if r_value ** 2 < rsquared_threshold:
                return False

            if var.info['SVTYPE'] == 'DEL':
                slope = -slope

            if slope < slope_threshold:
                return False

            return True
    return False

def to_bnd(var):
    var1 = copy.deepcopy(var)
    var2 = copy.deepcopy(var)

    # update svtype
    var1.info['SVTYPE'] = 'BND'
    var2.info['SVTYPE'] = 'BND'

    # update variant id
    var1.info['EVENT'] = var.var_id
    var2.info['EVENT'] = var.var_id
    var1.var_id = var.var_id + "_1"
    var2.var_id = var.var_id + "_2"
    var1.info['MATEID'] = var2.var_id
    var2.info['MATEID'] = var1.var_id
    
    # update position
    var2.pos = var.info['END']

    # update CIPOS and CIEND
    var2.info['CIPOS'] = var.info['CIEND']
    var2.info['CIEND'] = var.info['CIPOS']
    var2.info['CIPOS95'] = var.info['CIEND95']
    var2.info['CIEND95'] = var.info['CIPOS95']

    # delete svlen and END
    del var1.info['SVLEN']
    del var2.info['SVLEN']
    del var1.info['END']
    del var2.info['END']

    # add SECONDARY to var2
    var2.info['SECONDARY'] = True

    if var.info['SVTYPE'] == 'DEL':
        var1.alt = 'N[%s:%s[' % (var.chrom, var.info['END'])
        var2.alt = ']%s:%s]N' % (var.chrom, var.pos)

    elif var.info['SVTYPE'] == 'DUP':
        var1.alt = ']%s:%s]N' % (var.chrom, var.info['END'])
        var2.alt = 'N[%s:%s[' % (var.chrom, var.pos)
    return var1, var2

def reciprocal_overlap(a, b):
    # catch divide by zero error
    if a[1] == a[0] or b[1] == b[0]:
        return 0
    overlap = float(min(a[1], b[1]) - max(a[0], b[0]))
    return min(overlap / (a[1] - a[0]), overlap / (b[1] - b[0]))

def annotation_intersect(var, ae_dict, threshold):
    best_overlap = 0
    best_feature = ''
    slop = 100
    
    if var.chrom in ae_dict:
        var_start = var.pos - 1
        var_end = int(var.info['END'])
        i = 0
        while 1:
            # bail if end of dict
            if i >= len(ae_dict[var.chrom]):
                break
            feature = ae_dict[var.chrom][i]
            if feature[0] - slop < var_end:
                if feature[1] + slop > var_start:
                    overlap = reciprocal_overlap([var_start, var_end], feature)
                    if overlap > best_overlap:
                        best_overlap = overlap
                        best_feature = feature[2]
            else:
                break
            i += 1
        if best_overlap >= threshold:
            return best_feature

    return None

# primary function
def sv_classify(vcf_in, gender_file, ae_dict, f_overlap):
    vcf_out = sys.stdout
    vcf = Vcf()
    header = []
    in_header = True

    gender = {}
    # read sample genders
    for line in gender_file:
        v = line.rstrip().split('\t')
        gender[v[0]] = int(v[1])

    for line in vcf_in:
        if in_header:
            if line[0] == '#':
                header.append(line)
                continue
            else:
                in_header = False
                vcf.add_header(header)
                # write the output header
                vcf_out.write(vcf.get_header() + '\n')

        # parse variant line
        v = line.rstrip().split('\t')
        var = Variant(v, vcf)

        # check intersection with mobile elements
        if ae_dict is not None and var.info['SVTYPE'] in ['DEL']:
            ae = annotation_intersect(var, ae_dict, f_overlap)
            if ae is not None:
                if ae.startswith('SINE') or ae.startswith('LINE'):
                    ae = 'ME:' + ae
                var.alt = '<DEL:%s>' % ae
                var.info['SVTYPE'] = 'DEL:%s' % ae
                vcf_out.write(var.get_var_string() + '\n')
                continue

        # annotate based on read depth
        if var.info['SVTYPE'] in ['DEL', 'DUP']:
            if has_depth_support(var, gender):
                # write variant
                vcf_out.write(var.get_var_string() + '\n')
            else:
                for m_var in to_bnd(var):
                    vcf_out.write(m_var.get_var_string() + '\n')

        # variant is not DEL or DUP, just write out original
        else:
            vcf_out.write(var.get_var_string() + '\n')
    vcf_out.close()
    return


# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # load the annotated elements
    ae_dict = None
    if args.ae_path is not None:
        if args.ae_path.endswith('.gz'):
            ae_bedfile = gzip.open(args.ae_path, 'rb')
        else:
            ae_bedfile = open(args.ae_path, 'r')
        ae_dict = {}

        for line in ae_bedfile:
            v = line.rstrip().split('\t')
            if len(v) < 4:
                continue
            v[1] = int(v[1])
            v[2] = int(v[2])
            if v[0] in ae_dict:
                ae_dict[v[0]].append(v[1:])
            else:
                ae_dict[v[0]] = [v[1:]]

    # call primary function
    sv_classify(args.vcf_in, args.gender, ae_dict, args.f_overlap)

    # close the files
    args.vcf_in.close()
    args.gender.close()
    if args.ae_path is not None:
        ae_bedfile.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 

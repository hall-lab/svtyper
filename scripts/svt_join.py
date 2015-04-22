#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2015-04-13 14:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
svt_join.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Join genotyped VCFs from multiple samples")
    # parser.add_argument('-a', '--argA', metavar='argA', type=str, required=True, help='description of argument')
    # parser.add_argument('-b', '--argB', metavar='argB', required=False, help='description of argument B')
    # parser.add_argument('-c', '--flagC', required=False, action='store_true', help='sets flagC to true')
    parser.add_argument('-m', '--master', type=argparse.FileType('r'), default=None, help='VCF file to set first 8 columns of variant info [first file in vcf_list]')
    parser.add_argument('-q', '--sum_quals', required=False, action='store_true', help='Sum QUAL scores of input VCFs as output QUAL score')
    parser.add_argument('vcf_list', metavar='vcf', nargs='*', type=argparse.FileType('r'), default=None, help='VCF file(s) to join')

    # parse the arguments
    args = parser.parse_args()

    if len(args.vcf_list) < 1:
        parser.print_help()
        exit(1)

    # send back the user input
    return args

# primary function
def svt_join(master, sum_quals, vcf_list):
    header = []
    # draw header and variant info from first
    # VCF in the list
    # master = vcf_list[0]

    if master is None:
        master = open(vcf_list[0].name)

    sample_list = []

    # print header
    while 1:
        master_line = master.readline()
        if not master_line:
            break
        if master_line[:2] != "##":
            break
        print (master_line.rstrip())

    # get sample names
    # master_v = master_line.rstrip().split('\t')
    # for sample in master_v[9:]:
    #     sample_list.append(sample)
    for vcf in vcf_list:
        while 1:
            line = vcf.readline()
            if not line:
                break
            if line[:2] == "##":
                continue
            if line[0] == "#":
                line_v = line.rstrip().split('\t')
                for sample in line_v[9:]:
                    sample_list.append(sample)
                break
    print '\t'.join(master_line.rstrip().split('\t')[:9] + sample_list)
    
    # iterate through VCF body
    while 1:
        master_line = master.readline()
        if not master_line:
            break
        master_v = master_line.rstrip().split('\t')
        master_chrom = master_v[0]
        master_pos = master_v[1]

        out_v = master_v[:9] # output array of fields
        qual = float(out_v[5])

        # sys.stdout.write( '\t'.join(master_v) )
        for vcf in vcf_list:
            line = vcf.readline()
            if not line:
                sys.stderr.write('\nError: VCF files differ in length\n')
                exit(1)
            line_v = line.rstrip().split('\t')
            line_chrom = line_v[0]
            line_pos = line_v[1]
            
            # ensure that each VCF position agrees with the master
            if (master_chrom != line_chrom or
                master_pos != line_pos):
                sys.stderr.write('\nError: variant in %s (%s:%s) conflicts with master (%s:%s)\n' %
                                 (vcf.name, line_chrom, line_pos, master_chrom, master_pos))
                exit(1)

            qual += float(line_v[5])
            out_v = out_v + line_v[9:]
            # sys.stdout.write( '\t' + '\t'.join(line.rstrip().split('\t')[9:]) )
        if sum_quals:
            out_v[5] = qual
        sys.stdout.write( '\t'.join(map(str, out_v)) + '\n')
        # sys.stdout.write('\n')
    
    # close files
    master.close()
    for vcf in vcf_list:
        vcf.close()
    
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    svt_join(args.master, args.sum_quals, args.vcf_list)

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise

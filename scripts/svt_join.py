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
    parser.add_argument('vcf_list', metavar='vcf', nargs='*', type=argparse.FileType('r'), default=None, help='VCF file(s) to join')

    # parse the arguments
    args = parser.parse_args()

    if len(args.vcf_list) < 1:
        parser.print_help()
        exit(1)

    # send back the user input
    return args

# primary function
def svt_join(vcf_list):
    header = []
    # draw header and variant info from first
    # VCF in the list
    master = vcf_list[0]
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
    master_v = master_line.rstrip().split('\t')
    for sample in master_v[9:]:
        sample_list.append(sample)
    for vcf in vcf_list[1:]:
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
        out_v = [] # output array of fields
        master_line = master.readline()
        if not master_line:
            break
        out_v = master_line.rstrip().split('\t')
        qual = float(out_v[5])

        # sys.stdout.write( '\t'.join(master_v) )
        for vcf in vcf_list[1:]:
            line = vcf.readline()
            if not line:
                sys.stderr.write('\nError, VCF files differ in length\n')
                exit(1)
            line_v = line.rstrip().split('\t')
            qual += float(line_v[5])
            out_v = out_v + line_v[9:]
            # sys.stdout.write( '\t' + '\t'.join(line.rstrip().split('\t')[9:]) )
        out_v[5] = qual
        sys.stdout.write( '\t'.join(map(str, out_v)) + '\n')
        # sys.stdout.write('\n')
        
    for vcf in vcf_list:
        vcf.close()
    
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    svt_join(args.vcf_list)

    # close the input file
    # args.input.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise

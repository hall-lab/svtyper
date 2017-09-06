#!/usr/bin/env bash

set -e

# # make truth set and sort the BAMs
# BAM=/gscmnt/gc2802/halllab/sv_aggregate/MISC/realigned_BAMs/NA12878/NA12878.bam
# BAM_BASE=`basename $BAM .bam`
# python -m svtyper.core \
#     -i data/example.vcf \
#     -B $BAM \
#     -l data/$BAM_BASE.bam.json \
#     -w data/$BAM_BASE.target_loci.bam \
#     > data/example.gt.vcf
#
# sambamba sort data/NA12878.target_loci.bam

# run test
python -m svtyper.core \
    -i data/example.vcf \
    -B data/NA12878.target_loci.sorted.bam \
    -l data/NA12878.bam.json \
    > out.vcf

diff -I '^##fileDate=' data/example.gt.vcf out.vcf

#!/usr/bin/env bash

# # make truth set and sort the BAMs
# ../svtyper \
#     -i truth.vcf \
#     -B NA12893.bam \
#     -S NA12893.splitters.bam \
#     --dump NA12893.target_loci \
#     > truth.gt.vcf
# sambamba sort NA12893.target_loci.bam
# sambamba sort NA12893.target_loci.splitters.bam

# # find duplicates in BAM file
# sambamba view NA12893.target_loci.sorted.bam | sort -k1,1 | zapdups -v
# sambamba view NA12893.target_loci.splitters.sorted.bam | sort -k1,1 | zapdups -v

# run test
../svtyper \
    -i truth.vcf \
    -B NA12893.target_loci.sorted.bam \
    -S NA12893.target_loci.splitters.sorted.bam \
    -l NA12893.target_loci.lib_info.json \
    > test.gt.vcf

diff -I '^#' truth.gt.vcf test.gt.vcf
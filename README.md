SVTyper
=======
[![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/hall-lab/svtyper/master/LICENSE)
[![Build Status](https://travis-ci.org/hall-lab/svtyper.svg?branch=master)](https://travis-ci.org/hall-lab/svtyper)

Bayesian genotyper for structural variants

## Example

```sh
svtyper \
    -i sv.vcf \
    -B sample.bam \
    -l sample.bam.json \
    > sv.gt.vcf
```

## Overview

SVTyper performs breakend genotyping of structural variants using paired-end and split-read alignments from BWA-MEM.

## Installation

### Requirements
- Python 2.7 (or newer)
- Pysam 0.8.1 (or newer)

Clone the repository

```sh
git clone git@github.com:hall-lab/svtyper.git
```

Test the installation

```sh
cd svtyper/test

../svtyper \
    -i example.vcf \
    -B NA12878.target_loci.sorted.bam \
    -l NA12878.bam.json
    > test.vcf
```

## Citation

C Chiang, R M Layer, G G Faust, M R Lindberg, D B Rose, E P Garrison, G T Marth, A R Quinlan, and I M Hall. SpeedSeq: ultra-fast personal genome analysis and interpretation. Nat Meth (2015). doi:10.1038/nmeth.3505.

http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3505.html

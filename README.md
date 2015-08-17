svtyper
=======

Bayesian genotyper for structural variants

### Example workflow

#### Data
```
wget http://colbychiang.com/hall/cshl_sv_2014/data/NA12878.20.bam
wget http://colbychiang.com/hall/cshl_sv_2014/data/NA12878.20.bam.bai
wget http://colbychiang.com/hall/cshl_sv_2014/data/NA12878.20.splitters.bam
wget http://colbychiang.com/hall/cshl_sv_2014/data/NA12878.20.splitters.bam.bai
wget http://colbychiang.com/hall/cshl_sv_2014/data/NA12878.20.vcf.gz
```

#### Genotype with SVTyper
```
zcat NA12878.20.vcf.gz \
    | ./svtyper \
        -B NA12878.20.bam \
        -S NA12878.20.splitters.bam \
        > NA12878.20.gt.vcf
```
#### Warning
2015-06-17  
At the moment, SVTyper only works properly when BAMs were aligned using the -M flag to mark shorter reads as secondary alignments. This causes errors like [#7](https://github.com/cc2qe/svtyper/issues/7). We are working on a fix for this.

# dedup
A small script for NGS data deduplication, optimized for cfDNA sequencing data

# Dependency
This script depends on `pysam` module, please install it by `pip` or `easy_install` first

# Usage
```shell
./dedup.py -1 in.sorted.bam -o out.bam > dedup.log
```
Be sure that the in.sorted.bam is sorted

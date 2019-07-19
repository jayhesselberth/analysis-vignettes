#! /usr/bin/env bash

#BSUB -J rpfs_by_codon
#BSUB -e %J.err
#BSUB -o %J.out

python3 rpfs_by_codon.py \
    | gzip -c > rpfs.tsv.gz

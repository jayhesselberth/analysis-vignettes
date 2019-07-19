#! /usr/bin/env python3

#BSUB -J rpfs_by_codon
#BSUB -e %J.err
#BSUB -o %J.out

'''
Analyze ribosome profiling data.

Goal is to visualize the "ribosome shadow".

Strategy: XXX

Data is from Guydosh et al (2014) Mol Cell

Output is a TSV that we can plot in R (or matplotlib, if you want
to wait).
'''

from pandas import DataFrame
from pybedtools import BedTool
from pyfaidx import Fasta

from collections import defaultdict, Counter

import gzip
import sys
import pdb

def main():

    # fasta must be compressed with `bgzip` (from samtools dist), not gzip
    genome = Fasta("data/sacCer1.fa.bgz")

    # gene positions and ribosome-protected fragments ("rpfs")
    genes = BedTool("data/sgdGenes.sacCer1.bed.gz").sort()
    rpfs = BedTool("data/SRR1042864.pos.bg.gz").sort()

    # generate codons for each gene
    codons = load_gene_codons(genes)

    # map rpf signal onto codons
    codon_signal = map_signals(codons, rpfs)

    codon_table = make_codon_table()

    print(">> summing signals for codons ...", file=sys.stderr)
    codon_sums = defaultdict(dict)
    for ivl in codon_signal:

        # signal is column 7 
        signal = ivl.fields[6]
        if signal == '.': continue

        codon_seq = genome[ivl.chrom][ivl.start:ivl.end].seq
        codon_aa = codon_table[codon_seq]

        # dict[gene][codon_pos] = (codon_aa, signal)
        codon_sums[ivl.name][ivl.score] = (codon_aa, signal)

        pdb.set_trace()

    print(">> calcualting signal at codon offsets ...", file=sys.stderr)
    codon_offset_signal = defaultdict(Counter)
    for gene, codon_pos in codon_sums.items():
        for offset in range(-20, 20):
            codon_offset = codon_pos + offset
            try:
                signal_offset = codon_sums[gene][codon_offset]
            except KeyError:
                continue
            codon, signal_cur = codon_sums[gene][codon_pos]
            codon_offset_signal[codon][offset] += signal_offset

    for codon, offset in codon_offset_signal.items():
        signal = codon_offset_signal[codon][offset]
        print([codon, offset, signal], sep = "\t")

def make_codon_table():
    tbl = dict()

    with gzip.open("data/codon-table.txt.gz", 'rb') as ct:
        for line in ct:
            # decode the bytes
            line = line.decode("utf-8")

            if line.startswith("#"): continue

            fields = line.strip().split('\t')
            aa = fields[0]; codon = fields[1];

            tbl[codon] = aa

    return tbl

def map_signals(codons, rpfs):
    print(">> mapping rpf signal to codons ... ", file=sys.stderr, end='')

    res = codons.map(rpfs, c=4, o='sum')

    print("done", file=sys.stderr)

    return res

def load_gene_codons(genes):
    '''create BedTool of codon positions via intermediate DataFrame'''
    print(">> loading codons ... ", file=sys.stderr, end='')

    codons_df = DataFrame(
        codon_positions(genes),
        columns=['chrom', 'start', 'end',
        'name', 'score', 'strand'])

    res = BedTool.from_dataframe(codons_df).sort()

    print("done", file=sys.stderr)

    return res

def codon_positions(bt):
    ''' generator to make codons for each gene. `score` holds
    the codon number in each gene (starting at 1).'''
    for ivl in bt:

        # skip negative strand genes and genes with more than one exon
        if ivl.strand == "-" or int(ivl.fields[9]) > 1:
            continue

        for pos, codon_start in enumerate(range(ivl.start, ivl.end, 3)):
            yield {
                'chrom': ivl.chrom,
                'start': codon_start,
                'end': codon_start + 3,
                'name': ivl.name,
                'score': pos
            }

if __name__ == "__main__":
    main()


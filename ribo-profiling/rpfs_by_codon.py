#! /usr/bin/env python3

'''
Analyze ribosome profiling data. Goal is to visualize the "ribosome shadow".

Data is from Guydosh et al (2014) Mol Cell. Sample is WT yeast treated with
3-aminotriazole for 2 hours to inhibit histidine synthesis, depleting histidinyl-tRNAs
and stalling ribosomes at histidine codons.

Output are TSV that we can plot in R (or matplotlib, if you want
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
    print(">> mapping rpf signal to codons ... ", file=sys.stderr)
    codon_signal = codons.map(rpfs, c=4, o='sum')

    codon_table = make_codon_table()

    print(">> summing signals for codons ...", file=sys.stderr)
    codon_sums = defaultdict(dict)
    for ivl in codon_signal:

        # signal is column 7 
        signal = ivl.fields[6]
        if signal == '.': continue

        codon_seq = genome[ivl.chrom][ivl.start:ivl.end].seq

        # dict[gene][codon_pos] = (codon_aa, signal)
        codon_sums[ivl.name][int(ivl.score)] = (codon_seq, int(signal))

    print(">> calcualting signal at codon offsets ...", file=sys.stderr)

    codon_offset_signal = defaultdict(Counter)
    for gene, codon_data in codon_sums.items():
        for codon_pos, data in codon_sums[gene].items():

            codon_seq, signal = data

            for offset in range(-10, 11):
                codon_offset = codon_pos + offset
                try:
                    _, signal_offset = codon_sums[gene][codon_offset]
                except KeyError:
                    continue

                codon_offset_signal[codon_seq][offset] += signal_offset

    for codon_seq in codon_offset_signal:
        codon_aa = codon_table[codon_seq]
        for offset, signal in codon_offset_signal[codon_seq].items():
            print('\t'.join(map(str, [codon_seq, codon_aa, offset, signal])))

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


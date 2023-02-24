#!/usr/bin/env python

import sys
from Bio import SeqIO

def split(fastafile, prefix, nseqs):
    ofc = 1
    osq = 0
    outfile = "{}-{:03d}.fa".format(prefix, ofc)
    with open(fastafile, "r") as f:
        out = open(outfile, "w")
        sys.stdout.write(outfile + "\n")
        try:
            rec = SeqIO.parse(f, "fasta")
            for seq in rec:
                SeqIO.write(seq, out, "fasta")
                osq += 1
                if osq == nseqs:
                    out.close()
                    ofc += 1
                    outfile = "{}-{:03d}.fa".format(prefix, ofc)
                    out = open(outfile, "w")
                    sys.stdout.write(outfile + "\n")
                    osq = 0
        finally:
            out.close()

def main(args):
    if len(args) < 3:
        sys.stdout.write("""Usage: splitFasta.py input.fa prefix nseqs

Split a fasta file `input.fa' into multiple fasta files, each one containing `nseqs' sequences.
The output files have the form PREFIX-NNN.fa, where PREFIX is the second argument, and NNN is
a progressive number, padded with 0 to the left. 

The program writes the names of all output files to standard output, so they can be collected
in a file with:

  splitFasta.py input.fa prefix nseqs > outfiles.txt

""")
        sys.exit(1)
    split(args[0], args[1], int(args[2]))

if __name__ == "__main__":
    main(sys.argv[1:])

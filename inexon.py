#!/usr/bin/env python3

####################################################################################################
#                                                                                                  #
# inexon.py                                                                                        #
#                                                                                                  #
# Author: Katharina Hoff                                                                           #
# Program: None, was customized for Monica Sheffer's project                                       #
# E-Mail: katharina.hoff@uni-greifswald.de                                                         #
# Date: 20.5.2021                                                                                  #
#                                                                                                  #
####################################################################################################

__author__ = "Katharina J. Hoff"
__copyright__ = "Copyright 2021. All rights reserved."
__credits__ = "Matthis Ebel and Mario Stanke"
__version__ = "1.0.1"
__email__ = "katharina.hoff@uni-greifswald.de"
__status__ = "development"

import argparse
import re

''' parse command line options '''

parser = argparse.ArgumentParser(description="Determines for a position in a genome and a transcript id \
whether the position is in a coding exon (CDS) of that transcript. If yes, the number of the \
exon in the transcript is returned. The numbering starts from 1 and goes upstream to downstream.",
                                 epilog="""
SYNOPSIS

inexon --txid NM_004958.3 --chr chr1 --pos 11257000 --gff alltx.gtf

    --gff         Specifies the name of a file in GTF format that contains the transcripts'
                  coordinates.
    --txid        a transcript id
    --chr         name of the genomic sequence, e.g. contig, scaffold, chromosome id
    --pos         1-based sequence coordinate

  output:
    For each position and transcript combination the output is either -1 if the position is not
    inside a coding exon of that transcript or it is the number of the CDS exon that contains
    that position.    
""",
                                 formatter_class=argparse.RawTextHelpFormatter)

mandatoryArgs = parser.add_argument_group("mandatory arguments")

mandatoryArgs.add_argument("--gff", "--gtf",
                           dest="gff",
                           metavar="FILE",
                           help="Specifies the name of a file in GTF format that contains the transcripts' coordinates.",
                           required=True)


posfileArgs = parser.add_argument_group("position file arguments")

posfileArgs.add_argument("--posfile",
                         dest="posfile",
                         metavar="FILE",
                         help="""Name of an optional file that contains (multiple) input positions in tab separated format, e.g.
chr1     11248000
chr1     3244100""")

args = parser.parse_args()


''' read position file / create query list '''

positions = []

if (args.posfile):
    try:
        with open(args.posfile, "r") as fh:
            for line in fh:
                line = line.strip("\n")
                fields = line.split("\t")
                positions.append(tuple([fields[0], fields[1]]))
    except IOError:
        print("Could not open file with positions: " + args.posfile + "!")


''' read the GTF file '''

gtf = {}
try:
    with open(args.gff, "r") as fh:
        for line in fh:
            # skip comments and track lines
            if re.search("^#|track ", line, re.I):
                continue
            fields = line.split("\t")
            # skip non-exon features
            if (not re.search("CDS", fields[2], re.I)):
                continue
            seqname = str(fields[0])
            start = int(fields[3])
            end = int(fields[4])
            transcriptID = str(
                re.search('transcript_id "?([^;"]+)"?', fields[8]).group(1))
            if (seqname not in gtf):
                gtf[seqname] = []
            gtf[seqname].append(tuple([start, end, transcriptID]))
except IOError:
    print("Could not open GTF file: " + args.gff + "!")

# sort the exons
for seqname in gtf.keys():
        gtf[seqname].sort(key=lambda elem: elem[0])


''' check positions '''
for position in positions:
    seqname = position[0]
    pos = position[1]
    if (seqname not in gtf):
        print("Sequence", seqname, "is not there.")
        continue
    for exon in gtf[seqname]:
        if (int(pos) >= int(exon[0])) and (int(pos) <= int(exon[1])):
            print(seqname, "\t", str(pos), "\t", exon[2])
            break

#!/usr/bin/env python3

"""
Slices fasta files
"""

# Imports --------------------------------------------------------------------------------------------------
import argparse
from PyIsland import Islands

# Args -----------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="")

parser.add_argument(
    "-n", "--nucl", help="bucket/path/to/fasta", type=str, required=True
)

parser.add_argument("-f", "--gff", help="bucket/path/to/gff", type=str, required=True)

parser.add_argument(
    "-g", "--gbk", help="bucket/path/to/output/", type=str, required=True
)

parser.add_argument("-c", "--contig", help="contig accession", type=str, required=True)

parser.add_argument(
    "-s", "--start", help="start coordinate, 1-based", type=int, required=True
)

parser.add_argument(
    "-e", "--stop", help="stop coordinate, 1-based", type=int, required=True
)


args = parser.parse_args()

#  ------------------------------------------------------------------------------------------------------------------------


#  ------------------------------------------------------------------------------------------------------------------------


# make a gbk file
if __name__ == "__main__":
    Islands.make_gbk(args.nucl, args.gff, args.gbk, args.contig, args.start, args.stop)

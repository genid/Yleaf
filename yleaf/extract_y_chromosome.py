#!/usr/bin/env python
"""
Script for extracting Y chromosome data from a full reference genome.

Developed for Yleaf by Alaina Hardie, @trianglegrrl

License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.
"""

import argparse
import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
LOG = logging.getLogger(__name__)


def extract_y_chromosome(input_file: Path, output_file: Path):
    """
    Extract Y chromosome data from a full reference genome.

    Args:
        input_file: Path to the full reference genome file
        output_file: Path where the Y chromosome data will be saved
    """
    LOG.info(f"Extracting Y chromosome from {input_file}")

    with open(output_file, "w") as fo:
        with open(input_file) as fi:
            record = False
            y_chrom_found = False

            for line in fi:
                if line.startswith(">chrY") or line.startswith(">Y"):
                    record = True
                    y_chrom_found = True
                    LOG.info("Found Y chromosome")
                    fo.write(line)
                elif record:
                    if line.startswith(">"):
                        break
                    fo.write(line)

    if not y_chrom_found:
        LOG.error("No Y chromosome found in the reference file!")
        return False

    LOG.info(f"Y chromosome successfully extracted to {output_file}")
    return True


def main():
    parser = argparse.ArgumentParser(description="Extract Y chromosome from a reference genome")
    parser.add_argument("-i", "--input", required=True, type=str,
                      help="Input reference genome file (.fa, .fasta, or .fna)")
    parser.add_argument("-o", "--output", required=True, type=str,
                      help="Output file for Y chromosome (.fa, .fasta, or .fna)")
    args = parser.parse_args()

    input_file = Path(args.input)
    output_file = Path(args.output)

    if not input_file.exists():
        LOG.error(f"Input file {input_file} does not exist")
        sys.exit(1)

    valid_extensions = ['.fa', '.fasta', '.fna']
    if input_file.suffix.lower() not in valid_extensions:
        LOG.warning(f"Input file extension {input_file.suffix} is not a standard FASTA extension (.fa, .fasta, or .fna)")

    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Extract Y chromosome
    if not extract_y_chromosome(input_file, output_file):
        sys.exit(1)


if __name__ == "__main__":
    main()
#!/usr/bin/python
"""
Code for downloading the reference genome and extracting the specific y-chromomse data

Developed at Erasmus Medical Center Department of Genetic Identification

License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

Autor: Diego Montiel Gonzalez
Extensively modified by: Bram van Wersch
"""

import os
from pathlib import Path
import urllib.request
import gzip
import shutil
import urllib.request
import logging

from src import yleaf_constants

LOG: logging = logging.getLogger("yleaf_logger")


def main(
    choice: str
):
    if choice == yleaf_constants.HG19:
        reference_choice = [yleaf_constants.HG19]
    elif choice == yleaf_constants.HG38:
        reference_choice = [yleaf_constants.HG38]
    else:
        reference_choice = [yleaf_constants.HG19, yleaf_constants.HG38]

    # running downloader
    for dir_name in reference_choice:
        install_genome_files(dir_name)


def install_genome_files(
    reference_choice: str
):
    dir_path = yleaf_constants.DATA_FOLDER / reference_choice
    LOG.info(f"Starting with preparing {reference_choice}...")

    if not dir_path.exists():
        try:
            os.mkdir(dir_path)
        except IOError:
            LOG.error("Failed to create directory for output files")
            raise SystemExit("Failed to create directory for output files")

    ref_file = Path(dir_path / yleaf_constants.FULL_REF_FILE)
    ref_gz_file = Path(str(ref_file) + ".gz")
    try:
        if not ref_file.exists() and not ref_gz_file.exists():
            LOG.info(f"Downloading the {reference_choice} genome...")
            urllib.request.urlretrieve(f"http://hgdownload.cse.ucsc.edu/goldenPath/{reference_choice}"
                                       f"/bigZips/{reference_choice}.fa.gz", ref_gz_file)
        if not ref_file.exists():
            LOG.info("Unpacking the downloaded archive...")
            with gzip.open(ref_gz_file, 'rb') as f_in:
                with open(ref_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(ref_gz_file)
        ychrom_file = dir_path / yleaf_constants.Y_REF_FILE
        if not ychrom_file.exists():
            LOG.info("Writing Ychromosomal data")
            get_ychrom_data(ref_file, ychrom_file)
    # try and cleanup when user aborts, this attempts to not leave have downloaded files
    except KeyboardInterrupt:
        try:
            os.remove(ref_gz_file)
            os.remove(ref_file)
        # skip on IOerrors and such
        finally:
            raise


def get_ychrom_data(
    full_data_path: Path,
    yhcrom_file: Path
):
    with open(yhcrom_file, "w") as fo:
        with open(full_data_path) as fi:
            record = False
            for line in fi:
                if line == ">chrY\n":
                    record = True
                    fo.write(line)
                elif record:
                    if line.startswith(">"):
                        break
                    fo.write(line)

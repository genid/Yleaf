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

from yleaf import yleaf_constants

LOG: logging = logging.getLogger("yleaf_logger")


def _make_progress_hook(label: str):
    last_pct = [-1]
    def hook(block_num, block_size, total_size):
        if total_size <= 0:
            return
        pct = min(100, int(block_num * block_size * 100 / total_size))
        if pct >= last_pct[0] + 5:
            last_pct[0] = pct
            mb = block_num * block_size / (1024 * 1024)
            total_mb = total_size / (1024 * 1024)
            LOG.info(f"Downloading {label}: {pct}% ({mb:.0f} / {total_mb:.0f} MB)")
    return hook


def main(
    choice: str
):
    if choice == yleaf_constants.HG19:
        reference_choice = [yleaf_constants.HG19]
    elif choice == yleaf_constants.HG38:
        reference_choice = [yleaf_constants.HG38]
    elif choice == yleaf_constants.T2T:
        reference_choice = [yleaf_constants.T2T]
    else:
        reference_choice = [yleaf_constants.HG19, yleaf_constants.HG38]

    # running downloader
    for dir_name in reference_choice:
        install_genome_files(dir_name)


def install_genome_files(
    reference_choice: str
):
    LOG.info(f"Starting with preparing {reference_choice}...")

    if reference_choice == yleaf_constants.HG19:
        ref_file = yleaf_constants.HG19_FULL_GENOME
    elif reference_choice == yleaf_constants.T2T:
        ref_file = yleaf_constants.T2T_FULL_GENOME
    else:
        ref_file = yleaf_constants.HG38_FULL_GENOME

    ref_gz_file = Path(str(ref_file) + ".gz")
    try:
        if (not ref_file.exists() or os.path.getsize(ref_file) < 100) and not ref_gz_file.exists():

            LOG.info(f"Downloading {reference_choice} reference genome (~3 GB). This may take several minutes...")
            # T2T (CHM13v2.0 + HG002 Y) is hosted on UCSC as "hs1"
            ucsc_name = "hs1" if reference_choice == yleaf_constants.T2T else reference_choice
            urllib.request.urlretrieve(
                f"http://hgdownload.cse.ucsc.edu/goldenPath/{ucsc_name}/bigZips/{ucsc_name}.fa.gz",
                ref_gz_file,
                reporthook=_make_progress_hook(reference_choice),
            )
        if not ref_file.exists() or os.path.getsize(ref_file) < 100:
            LOG.info(f"Decompressing {reference_choice} reference genome, please wait...")
            with gzip.open(ref_gz_file, 'rb') as f_in:
                with open(ref_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(ref_gz_file)

        if reference_choice == yleaf_constants.HG19:
            ychrom_file = yleaf_constants.HG19_Y_CHROMOSOME
        elif reference_choice == yleaf_constants.T2T:
            ychrom_file = yleaf_constants.T2T_Y_CHROMOSOME
        else:
            ychrom_file = yleaf_constants.HG38_Y_CHROMOSOME

        if not ychrom_file.exists() or os.path.getsize(ychrom_file) < 100:
            LOG.debug("Writing Ychromosomal data")
            get_ychrom_data(ref_file, ychrom_file)
    # try and cleanup when user aborts, this attempts to not leave half downloaded files
    except KeyboardInterrupt:
        try:
            os.remove(ref_gz_file)
            os.remove(ref_file)
        # skip on IOerrors and such
        finally:
            raise


_FULL_GENOME_PATHS = {
    yleaf_constants.HG19: yleaf_constants.HG19_FULL_GENOME,
    yleaf_constants.HG38: yleaf_constants.HG38_FULL_GENOME,
    yleaf_constants.T2T:  yleaf_constants.T2T_FULL_GENOME,
}

_YCHROM_PATHS = {
    yleaf_constants.HG19: yleaf_constants.HG19_Y_CHROMOSOME,
    yleaf_constants.HG38: yleaf_constants.HG38_Y_CHROMOSOME,
    yleaf_constants.T2T:  yleaf_constants.T2T_Y_CHROMOSOME,
}


def download_chry_fasta(reference_choice: str) -> Path:
    """Return the chrY FASTA for the given build.

    Resolution order (identical for all builds):
    1. chrY FASTA already present → return immediately.
    2. Full genome already downloaded → extract chrY from it.
    3. Neither present → download full genome (same path as BAM mode) then extract.
    """
    dest = _YCHROM_PATHS[reference_choice]

    if dest.exists() and dest.stat().st_size > 100:
        return dest

    dest.parent.mkdir(parents=True, exist_ok=True)

    full_genome = _FULL_GENOME_PATHS[reference_choice]
    if not (full_genome.exists() and full_genome.stat().st_size > 100):
        LOG.info(f"Full genome for {reference_choice} not found, downloading now...")
        install_genome_files(reference_choice)

    LOG.info(f"Extracting chrY FASTA from {reference_choice} full genome reference...")
    get_ychrom_data(full_genome, dest)
    LOG.info(f"chrY FASTA ready: {dest}")
    return dest


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

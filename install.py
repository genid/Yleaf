#!/usr/bin/python

import subprocess
import os
from pathlib import Path
import urllib.request
import gzip
import shutil
import urllib.request
import sys


_HG19 = "hg19"
_HG38 = "hg38"
_BOTH = "both"
_EXPECTED_BWA_EXTENSIONS = ("0123", "amb", "ann", "bwt.2bit.64", "pac")


def main():
    # allow to give an argv argument
    args = sys.argv[1:]
    options = [_HG19, _HG38, _BOTH]

    # choose what to install
    choice = None
    for option in options:
        if option in args:
            choice = option
    if choice is None:
        while True:
            choice = input(f'Please choose a build Genome; {"/".join(options)}: ')
            choice = choice.lower()
            if choice not in options:
                print("Invalid option selected.")
            else:
                break

    if choice == _HG19:
        dir_names = [_HG19]
    elif choice == _HG38:
        dir_names = [_HG38]
    else:
        dir_names = [_HG19, _HG38]

    # running downloader
    for dir_name in dir_names:
        install_genome_files(dir_name)


def install_genome_files(dir_name):
    dir_path = Path(dir_name)
    print(f"Starting with preparing {dir_name}...")

    if not dir_path.exists():
        try:
            os.mkdir(dir_path)
        except IOError:
            raise IOError("Failed to create directory for output files")
    if not Path(dir_path / f"{dir_name}.fa").exists() and not Path(f"{dir_name}.fa.gz").exists():
        print(f"Downloading the {dir_name} genome...")
        urllib.request.urlretrieve(f"http://hgdownload.cse.ucsc.edu/goldenPath/{dir_name}/bigZips/{dir_name}.fa.gz",
                                   f"{dir_name}.fa.gz")
    full_fasta_file = dir_path / f"{dir_name}.fa"
    if not full_fasta_file.exists():
        print("Unpacking the downloaded archive...")
        with gzip.open(f"{dir_name}.fa.gz", 'rb') as f_in:
            with open(dir_path / f"{dir_name}.fa", 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(f"{dir_name}.fa.gz")
    ychrom_file = dir_path / "chrY.fa"
    if not ychrom_file.exists():
        print("Writing Ychromosomal data")
        get_ychrom_data(full_fasta_file, ychrom_file)


def get_ychrom_data(full_data_path: Path, yhcrom_file: Path):
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


if __name__ == '__main__':
    main()

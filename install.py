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
    if not Path(dir_path / f"{dir_name}.fa").exists():
        print("Unpacking the downloaded archive...")
        with gzip.open(f"{dir_name}.fa.gz", 'rb') as f_in:
            with open(dir_path / f"{dir_name}.fa", 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(f"{dir_name}.fa.gz")

    # downloading and preparing snp data. Is comparatively quick
    snp_file = Path(dir_path / "snp_data.txt.gz")
    unzipped_snp_file = Path(dir_path / "snp_data.txt")
    filtered_snp_file = Path(dir_path / "snp_data_filtered.csv")

    if not snp_file.exists() and not unzipped_snp_file.exists() and not filtered_snp_file.exists():
        print("Downloading snp data...")
        if dir_name == _HG38:
            urllib.request.urlretrieve("https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/"
                                       "chr_rpts/chr_Y.txt.gz", filename=str(snp_file))
        else:
            urllib.request.urlretrieve("https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/"
                                       "chr_rpts/chr_Y.txt.gz", filename=str(snp_file))
    if not unzipped_snp_file.exists() and not filtered_snp_file.exists():
        with gzip.open(snp_file, 'rb') as f_in:
            with open(unzipped_snp_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(snp_file)
    if not filtered_snp_file.exists():
        filter_snp_data(unzipped_snp_file, filtered_snp_file)
    print(f"Finished downloading all for {dir_name}")


def filter_snp_data(unzipped_snp_file, filtered_snp_file):
    print("Filtering out usefull SNPs...")
    with open(unzipped_snp_file) as fi:
        # ignore the beatifull consistent great header
        skip_header_amnt = 7
        [fi.readline() for _ in range(skip_header_amnt)]
        with open(filtered_snp_file, "w") as fo:
            fo.write("rs_nr,location,minor_allele,frequency\n")
            for index_2, line in enumerate(fi):
                values = line.split('\t')
                try:
                    rs_nr = values[0]
                    location = values[11]
                    freq_info = values[-1] if ":" in values[-1].strip() else None
                except IndexError:
                    print(f"Failed to read line {index_2 + skip_header_amnt}.")
                    continue
                if freq_info is None:
                    continue
                try:
                    minor_allele, _, freq = freq_info.strip().split(":")
                except IndexError:
                    # missing information
                    continue
                # minor allele is not a single nucleotide
                if len(minor_allele) > 1:
                    continue
                fo.write(f"{rs_nr},{location},{minor_allele},{freq}\n")
    os.remove(unzipped_snp_file)


if __name__ == '__main__':
    main()

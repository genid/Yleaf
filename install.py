#!/usr/bin/python

import subprocess
import os
from pathlib import Path
import urllib.request
import gzip
import shutil


_HG19 = "hg19"
_HG38 = "hg38"
_BOTH = "both"
_EXPECTED_BWA_EXTENSIONS = ("0123", "amb", "ann", "bwt.2bit.64", "pac")

options = [_HG19, _HG38, _BOTH]
while True:
    choice = input(f'Please choose a build Genome; {"/".join(options)}: ')
    choice = choice.lower()
    if choice not in options:
        print("Invalid option selected.")
    break

if choice == _HG19:
    dir_names = [_HG19]
elif choice == _HG38:
    dir_names = [_HG38]
else:
    dir_names = [_HG19, _HG38]

for index in range(len(dir_names)):
    files_already_present = False
    dir_name = dir_names[index]
    dir_path = Path(dir_name)
    print(f"Starting with preparing {dir_name}")

    if not dir_path.exists():
        try:
            os.mkdir(dir_path)
        except IOError:
            raise IOError("Failed to create directory for output files")
    if not Path(dir_path / f"{dir_name}.fa").exists() and not Path(f"{dir_name}.fa.gz").exists():
        print(f"Downloading the {dir_name} genome")
        urllib.request.urlretrieve(f"http://hgdownload.cse.ucsc.edu/goldenPath/{dir_name}/bigZips/{dir_name}.fa.gz",
                                   f"{dir_name}.fa.gz")
    else:
        files_already_present = True
        print("Already found existing archive. Using existing files.")
    if not Path(dir_path / f"{dir_name}.fa").exists():
        print("Unpacking the downloaded archive")
        with gzip.open(f"{dir_name}.fa.gz", 'rb') as f_in:
            with open(dir_path / f"{dir_name}.fa", 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(f"{dir_name}.fa.gz")
    else:
        files_already_present = True
    print(f"Finished downloading all for {dir_name}")

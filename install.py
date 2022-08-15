#!/usr/bin/python

import subprocess
import os
from pathlib import Path
import urllib.request
import zipfile

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
    if not Path(dir_path / f"{dir_name}.fa").exists() and not Path(dir_path / f"{dir_name}.fa.gz").exists():
        print(f"Downloading the {dir_name} genome")
        urllib.request.urlretrieve(f"http://hgdownload.cse.ucsc.edu/goldenPath/{dir_name}/bigZips/{dir_name}.fa.gz",
                                   f"{dir_name}.fa.gz")
    else:
        files_already_present = True
        print("Already found existing archive. Using existing dowloaded file...")
    if not Path(dir_path / f"{dir_name}.fa").exists():
        print("Unpacking the downloaded archive")
        with zipfile.ZipFile(f"{dir_name}.fa.gz", "r") as zip_ref:
            zip_ref.extractall(dir_path)
    else:
        files_already_present = True
    have_to_index = False
    for extension in _EXPECTED_BWA_EXTENSIONS:
        found = False
        for file in dir_path.iterdir():
            if file.name.endswith(extension):
                found = True
                break
        if not found:
            have_to_index = True
            break
    if have_to_index:
        print("Starting with bwa-mem2 indexing for faster fastq runs")
        cmd = "bwa-mem2 index -a bwtsw hg19/hg19.fa"
        process = subprocess.Popen(cmd, shell=True)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            print(f"Failed to create bwa-mem2 index, with message {stderr.decode('utf-8')}")
            raise SystemExit("Failed command execution")
    else:
        files_already_present = True
    if files_already_present:
        print("All required files have been downloaded and prepared. Certain files that are already present have not"
              f" been remade. If you experience issues running Yleaf in fastq mode then consider deleting the"
              f" {dir_name} folder entirely and running install again.")
    else:
        print(f"Finished downloading all for {dir_name}")
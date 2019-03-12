#!/usr/bin/python

import subprocess
import argparse
import os
import sys
import errno

if not os.path.isdir("tmp"):
    os.mkdir("tmp")

def force_symlink(file1, file2):
    try:
        os.symlink(file1, file2)
    except OSError as e:
        if e.errno == errno.EEXIST:
            os.remove(file2)
            os.symlink(file1, file2)

homedir = os.getenv("HOME")
homebin = os.path.join(homedir, "bin")
parser = argparse.ArgumentParser(description="Install yleaf")
parser.add_argument("--prefix", default=homebin, help="Folder to put executables in")
args = parser.parse_args()

# get hg19 data
if os.path.isdir("index_hg19"):
    print("Reference genome hg19 does exist, \
I assume this is because you have installed \
this software before, before, and I will skip this step. \
In case of doubt, delete it and run install.py again.")
else:
    if os.path.isfile("hg19.tar.gz"):
        print("Extracting from existing archive...")
    else:
        print("Downloading hg19 reference genome...")
        os.system("wget -O hg19.tar.gz \
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz")

    os.system("gunzip hg19.tar.gz")
    os.system("tar -xf hg19.tar")
    os.system("mkdir index_hg19")
    os.system("mv *.fa index_hg19")
    os.system("cat index_hg19/*.fa >> hg19.fa")
    os.system("rm index_hg19/*.fa")
    os.system("mv hg19.fa index_hg19")
    
    print "Building index from hg19 genome, this may take some time..."
    cmd = "bwa index -a bwtsw index_hg19/hg19.fa" 
    subprocess.call(cmd, shell=True)

# make sure prefix dir is there
if not os.path.exists(args.prefix):
    try:
        os.mkdir(args.prefix)
    except:
        sys.stderr.write("Failed to create directory {0}! \n".format(args.prefix))

# get current location of scripts
app_folder = os.path.dirname(os.path.realpath(__file__))
print("yleaf executables are in {0}".format(app_folder))

# link to executables
print("Linking to executables...")
yleaf_py = os.path.join(app_folder, "Yleaf.py")
yleaf_py_link = os.path.join(args.prefix, "Yleaf.py")
#clean_tree_pr = os.path.join(app_folder, "clean_tree_printout.sh")
#clean_tree_pr_link = os.path.join(args.prefix, "clean_tree_printout.sh")
force_symlink(yleaf_py, yleaf_py_link)
#force_symlink(clean_tree_pr, clean_tree_pr_link)

# make sure prefix dir is in PATH
paths = os.getenv("PATH").split(":")
if not args.prefix in paths:
    bashprofile = os.path.join(homedir + ".profile")
    with open(bashprofile, "a") as profilefh:
        profilefh.write("""\nPATH="{0}:$PATH" """.format(args.prefix))

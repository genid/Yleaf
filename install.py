#!/usr/bin/python

import subprocess
import argparse
import os
import sys
import errno

flag = True
while flag:
    print("Please choose a build Genome reference to download and process with BWA-MEM")
    choice = input("hg19/hg38: ")            
    if str(choice) == "hg19":                
        # get hg19 data
        if os.path.isdir("hg19"):
            print("Reference genome hg19 does exist, \
        I assume this is because you have installed \
        this software before, this step will be skip. \
        In case of doubt, delete it and run install.py again.")
        else:
            if os.path.isfile("hg19.tar.gz"):
                print("Extracting from existing archive...")
                os.system("gunzip hg19.tar.gz")
                os.system("tar -xf hg19.tar")
                os.system("mkdir hg19")
                os.system("mv *.fa hg19")
                os.system("cat hg19/*.fa >> hg19.fa")
                os.system("rm hg19/*.fa")
                os.system("mv hg19.fa hg19")                
                print("Building index from hg19 genome, this may take some time...")
                cmd = "bwa index -a bwtsw hg19/hg19.fa" 
                subprocess.call(cmd, shell=True)
                flag = False  
            else:
                print("Downloading hg19 reference genome...")
                os.system("wget -O hg19.tar.gz \
        http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz")

                os.system("gunzip hg19.tar.gz")
                os.system("tar -xf hg19.tar")
                os.system("mkdir hg19")
                os.system("mv *.fa hg19")
                os.system("cat hg19/*.fa >> hg19.fa")
                os.system("rm hg19/*.fa")
                os.system("mv hg19.fa hg19")                
                print("Building index from hg19 genome, this may take some time...")
                cmd = "bwa index -a bwtsw hg19/hg19.fa" 
                subprocess.call(cmd, shell=True)
                flag = False  

    elif str(choice) == "hg38":                
        # get hg38 data
        if os.path.isdir("hg38"):
            print("Reference genome hg38 does exist, \
        I assume this is because you have installed \
        this software before, before, and I will skip this step. \
        In case of doubt, delete it and run install.py again.")
        else:
            if os.path.isfile("hg38.fa.gz"):
                print("Extracting from existing archive...")
                os.system("gunzip hg38.fa.gz")    
                os.system("mkdir hg38")    
                os.system("mv hg38.fa hg38")                
                print("Building index from hg38 genome, this may take some time...")
                cmd = "bwa index -a bwtsw hg38/hg38.fa" 
                subprocess.call(cmd, shell=True)        
                flag = False                
            else:
                print("Downloading hg38 reference genome...")
                os.system("wget -O hg38.fa.gz \
        http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz")

                os.system("gunzip hg38.fa.gz")    
                os.system("mkdir hg38")    
                os.system("mv hg38.fa hg38")
                
                print("Building index from hg38 genome, this may take some time...")
                cmd = "bwa index -a bwtsw hg38/hg38.fa" 
                subprocess.call(cmd, shell=True)        
                flag = False                
    else:
        print("Please type hg19 or hg38")                               

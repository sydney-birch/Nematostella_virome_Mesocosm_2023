#! /usr/bin/env python3

#import modules 
import argparse
import subprocess
import os
import gzip
import time
import shutil

#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("-a", help="path to fastqc output")

args = parser.parse_args()

#This script iterates through the directory its placed in and runs fastqc - give it the path to where you want the output

def run_fastqc (output_path):

    for item in os.scandir():
        #print("iterating on item: ", item) #for debug
        if item.is_file():
            if "NH_" in item.name: 
                current_sample = item.name
                print("current sample: ", current_sample)
                result = subprocess.run("fastqc {0} -o {1}".format(current_sample, output_path), shell=True)
                print("moving to next sample")
                
    print("Fastqc is complete!")
    
#call funciton
fastqc_result = run_fastqc(args.a)    


## Run in dir with fastq files - move script and slurm to fastqc dir after running

#pwd: /scratch/sbirch1/Nematostella_transcriptomics/raw_reads_2-19-24/21081FL-06-01-01_S1_L001_R1_001.fastq.gz+_ds.c38327b14c2c4b9e82c49d5a8328b70e/Fastq-021324-XP-fcB-22HGNTLT3-L001-I8I8

#./1_fastqc.py -a ../../../1_fastqc/before_trim

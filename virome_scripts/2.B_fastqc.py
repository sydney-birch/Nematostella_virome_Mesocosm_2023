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
parser.add_argument("-b", help="path to fastq files")

args = parser.parse_args()

#This script iterates through the directory given and runs fastqc on the paired reads - give it the path to where you want the output

def run_fastqc (input_path, output_path):
    print("input path: ", input_path)
    os.chdir(input_path)
    
    for item in os.scandir():
        #print("iterating on item: ", item) #for debug
        if item.is_file():
            if "P.fq.gz" in item.name: 
                current_sample = item.name
                print("current sample: ", current_sample)
                result = subprocess.run("fastqc {0} -o {1}".format(current_sample, output_path), shell=True)
                print("moving to next sample")
                
    print("Fastqc is complete!")
    
#call funciton
fastqc_result = run_fastqc(args.b, args.a)    


## Run in fastqc dir

#./2.B_fastqc.py -a after_trim -b ../2_trimmomatic

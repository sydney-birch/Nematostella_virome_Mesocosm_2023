#! /usr/bin/env python3

#import modules 
import argparse
import subprocess
import os


#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("-a", help="name of dir with fastq files")
#parser.add_argument("-b", help="name of new file {0}-read_counts.txt")

args = parser.parse_args()

## This script goes through a given directory and counts the number of reads in a fastq file 

 
def count_reads(input_dir):
    os.chdir(input_dir)
    #with open("{0}-read_counts.txt".format(file_name), "w") as out_handle:
    for item in os.scandir():
        if ".fq.gz" in item.name:
            sample_name = item.name
            #print("On Sample: ", sample_name)
            result = subprocess.run("echo -n \"{0} number of reads: \"".format(sample_name), shell=True)
            result = subprocess.run("zcat {0} | grep @ | wc -l".format(sample_name), shell=True)        

#call funciton
result = count_reads(args.a)

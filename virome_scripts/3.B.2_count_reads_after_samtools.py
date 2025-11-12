#! /usr/bin/env python3

#import modules 
import argparse
import subprocess
import os


#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
#parser.add_argument("-a", help="name of dir with fastq files")
#parser.add_argument("-b", help="name of new file {0}-read_counts.txt")

args = parser.parse_args()

## This script goes through a given directory after samtools and goes into each sample dir and counts the number of reads in a fastq file 
## Run in dir with all the sample dirs

#Step 1: prep files - make a list of your samples which each have their own dir (remove any unwanted dirs)
 
def sample_list_prep ():
    sample_list = []

    for item in os.scandir():
        #print("iterating on item: ", item) #for debug
        if item.is_dir():
            sample_name = item.name
            print("sample Name: ", sample_name)
        
            sample_list.append(sample_name)
            
    #sample_list.remove(rem_dir)
    #sample_list.remove(rem_b)
    return sample_list

#call funciton
sample_list_result = sample_list_prep()
print("structure of sample_list_result: ", sample_list_result) #for debug


print("moving to next function")
 
def count_reads():
    for dir_name in sample_list_result:
        os.chdir(dir_name)
    
        for item in os.scandir():
            if ".fastq.gz" or ".fq.gz" in item.name:
                sample_name = item.name
                #print("On Sample: ", sample_name)
                result = subprocess.run("echo -n \"{0} number of reads: \"".format(sample_name), shell=True)
                result = subprocess.run("zcat {0} | grep @ | wc -l".format(sample_name), shell=True) 
        os.chdir("..")   

#call funciton
result_1 = count_reads()

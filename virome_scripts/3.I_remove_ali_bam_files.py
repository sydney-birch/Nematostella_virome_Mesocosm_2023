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
#parser.add_argument("-r", help="name of dirs not to include")
#parser.add_argument("-a", help="name of dirs not to include")
#parser.add_argument("-a", help="name of next dir for mapping")
#parser.add_argument("-b", help="number in sam alignement file to use in samtools(if using second alignment sam file - use 2)")

args = parser.parse_args()

## This script is for the samtools processing of MAPPED READS - it takes the sam files from the first HISAT2 and turns it 
## into a bam file then creates a fastq file of only sequences that did map 

##This script should be run in the 2nd ali dir - it will use the sam file from the 1st alignment to get the mapped reads
	## sam file structure: Nematostella_genome_ali_1_{1}.sam

#Step 1: prep files - make a list of your samples which each have their own dir 
 
def sample_list_prep ():
    sample_list = []

    for item in os.scandir():
        #print("iterating on item: ", item) #for debug
        if item.is_dir():
            sample_name = item.name
            print("sample Name: ", sample_name)
        
            sample_list.append(sample_name)
            
    #print("group db: ", group_db)#for debug

    return sample_list

#call funciton
sample_list_result = sample_list_prep()
print("structure of sample_list_result: ", sample_list_result) #for debug

print("moving to next function")


#Step 2: Run sam tools  

def run_rm_ali_file():
    for dir_name in sample_list_result:
        os.chdir(dir_name)
        print ("cwd: ", os.getcwd())
        for item in os.scandir():
            print("current item: ", item)
            if item.name == "align_reads.bam": 
                print("will remove current item: ", item.name)
                #mv to RNA_filt dir for alignment
                result = subprocess.run("rm align_reads.bam", shell=True)
        os.chdir("..")
        print("cwd: ", os.getcwd())

    print("All Samtools processing is complete")

#call run_samtools function
result_1 = run_rm_ali_file()

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
#parser.add_argument("--rem", "-r", nargs="?", help="name of dirs not to include")
#parser.add_argument("-b", help="name of dirs not to include")
#parser.add_argument("-a", help="name of next dir for mapping")
parser.add_argument("-b", help="number in sam alignement file to use in samtools(if using second alignment sam file - use 2)")

args = parser.parse_args()

## This script is for the samtools processing - it takes the sam files from HISAT2 and turns it 
## into a bam file then creates a fastq file of only sequences that did not map 

##This script should be run in the dir that was created by the hisat2 py script
	## For example if hisat2 was run in the 1st ali dir to generate the 2nd aligement - it created a 2nd_ali dir
	## Run this script in that 2nd_ali dir which has all sub dirs and within those dirs contains the sam file 
	## in this example it would start with Nematostella_genome_ali_2

#Step 1: prep files - make a list of your samples which each have their own dir (remove any unwanted dirs)
 
def sample_list_prep ():
    sample_list = []
 #   os.mkdir("fastqc_results")
    for item in os.scandir():
        #print("iterating on item: ", item) #for debug
        if item.is_dir():
            sample_name = item.name
            print("sample Name: ", sample_name)
        
            sample_list.append(sample_name)
            
    #print("group db: ", group_db)#for debug
    #sample_list.remove(rem_dir)
    #sample_list.remove(rem_b)
    return sample_list

#call funciton
sample_list_result = sample_list_prep()
print("structure of sample_list_result: ", sample_list_result) #for debug

print("moving to next function")


#Step 2: Run sam tools  

def run_samtools(ali_num):
    for dir_name in sample_list_result:
        os.chdir(dir_name)
        print ("cwd: ", os.getcwd())
        
        #run first part
        result = subprocess.run("samtools view Nematostella_genome_ali_{0}_*.sam -f 0x4 -h -b -o unalign_reads.bam".format(ali_num), shell=True)
        #Run second part 
        result = subprocess.run("samtools collate -u -O unalign_reads.bam | \
        samtools fastq -1 unaligned_{0}_{1}_ali_paired1.fq.gz -2 unaligned_{0}_{1}_ali_paired2.fq.gz -0 /dev/null -s /dev/null -n".format(dir_name, ali_num), shell=True)
        print("samtools complete on {0}".format(dir_name))
        os.chdir("..")
        print("cwd: ", os.getcwd())

    print("All Samtools processing is complete")

#call run_samtools function
result_1 = run_samtools(args.b)

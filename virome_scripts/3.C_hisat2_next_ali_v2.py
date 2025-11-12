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
parser.add_argument("-c", help="name of next dir for mapping - if making 3rd ali - the next dir will be 4th_ali")
parser.add_argument("-b", help="number to name new sam alignement file (if making 3rd alignment sam file - use 3)")

args = parser.parse_args()

#Run this in the 1st ali dir (3_HISAT2/1st_ali_11-9-23 there will be sub dirs in each) - 
#the output will be a sam file from the 2nd ali - move this into the 2nd ali dir (3_HISAT2/2nd_ali_12-4-23/)

## This script is for the next hisat2 mapping (not the initial - any mapping after) - it uses the two fastq files
## generated from samtools (all reads that didn't map) and maps them again to the genome. It outputs a samfile 
## and the next main dir for mapping and moves the sam file into the new mapping dir with the subdir of the sample 


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


#Step 2: Run Hisat2 - loop through the list you just created - go into each dir and run the 
        #hisat2 code and move back to the main dir and create 

def run_hisat2(ali_num,next_dir_name):
    os.mkdir("../{0}".format(next_dir_name))

    for dir_name in sample_list_result:
        os.chdir(dir_name)
        print ("cwd: ", os.getcwd())
        #Run hisat2
        #result = subprocess.run("echo \"starting {0}\" ".format(dir_name), shell=True)
        result = subprocess.run("hisat2 -q -p 12 -x ../../1st_ali_1-23-25/Nematostella_genome -1 *_paired1.fq.gz -2 *_paired2.fq.gz -S Nematostella_genome_ali_{0}_{1}.sam".format(ali_num,dir_name), shell=True)
        print("hisat2 complete on {0}".format(dir_name))
        print("cwd after hisat2 b4 mv: ", os.getcwd())
        #make new dir for next round of mapping and move sam alignement file into that dir 
        os.mkdir("../../{0}/{1}".format(next_dir_name,dir_name))
        result = subprocess.run("mv Nematostella_genome_ali_{0}_{1}.sam ../../{2}/{1}".format(ali_num,dir_name,next_dir_name), shell=True)
        print("file move complete")
        os.chdir("..")
    print("All hisat2 is complete")

#call run_fastqc function
result_1 = run_hisat2(args.b,args.c)


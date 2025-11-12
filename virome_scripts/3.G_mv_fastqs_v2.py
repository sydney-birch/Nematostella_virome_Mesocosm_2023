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
parser.add_argument("-c", help="name of main dir to move fastq files into (this dir should have a fastq_files dir within it)")
parser.add_argument("-b", help="number to name new sam alignement file (if making 3rd alignment sam file - use 3)")

args = parser.parse_args()


#Step 1: prep files - make a list of your samples which each have their own dir (remove any unwanted dirs)
 
def sample_list_prep ():
    sample_list = []

    for item in os.scandir():
        #print("iterating on item: ", item) #for debug
        if item.is_dir():
            sample_name = item.name
            print("sample Name: ", sample_name)
        
            sample_list.append(sample_name)

    return sample_list

#call funciton
sample_list_result = sample_list_prep()
print("structure of sample_list_result: ", sample_list_result) #for debug


print("moving to next function")


#Step 2: move fastq files into RNA filt dir into the unmapped_reads dir

def mv_fq_files(ali_num,next_dir_name):
    #os.mkdir("../{0}".format(next_dir_name))
    for dir_name in sample_list_result:
        os.chdir(dir_name)
        print("cwd: ", os.getcwd())           
        #os.mkdir("../../{0}/{1}".format(next_dir_name,dir_name))
        result = subprocess.run("mv unaligned_{1}_{0}_ali_paired1.fq.gz ../../../{2}/unmapped_fastq_files/".format(ali_num,dir_name,next_dir_name), shell=True)
        result = subprocess.run("mv unaligned_{1}_{0}_ali_paired2.fq.gz ../../../{2}/unmapped_fastq_files/".format(ali_num,dir_name,next_dir_name), shell=True)
        print("file move complete")
        os.chdir("..")
        
#call run_fastqc function
result = mv_fq_files(args.b,args.c)   


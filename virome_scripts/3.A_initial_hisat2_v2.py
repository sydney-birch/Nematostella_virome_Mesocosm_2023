#! /usr/bin/env python3

#import modules 
import argparse
import subprocess
import os
import gzip


#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("-a", help="name of dir with fastqs after trimmommatic - make sure names are in the following format NH_T0-NH_B4-filtered_2P.fq")
parser.add_argument("-b", help="number to name new sam alignement file (if making 3rd alignment sam file - use 3)")
parser.add_argument("-c", help="name of next dir for mapping")
parser.add_argument("-d", help="path back from fastq dir to 1st ali dir")


args = parser.parse_args()

##Run this in the 1st ali dir (3_HISAT2/1st_ali_11-9-23 there will not be sub dirs - this scripts makes them) 
## It outputs a samfile and the next main dir for mapping and moves the sam file into the new mapping dir with the subdir of the sample 

#NH_T0-SC_B5_filtered_1P.fq.gz
#0     1   2   3        4

#Step 1: prep files - make a list of your samples which each have their own dir (remove any unwanted dirs)
 
def sample_list_prep (input_path, back_path):
    sample_list = []
    print("input path: ", input_path)
    os.chdir(input_path)
    
    for item in os.scandir():
        #print("iterating on item: ", item) #for debug
        if item.is_file():
            if "filtered_1P.fq" in item.name: 
                temp_list= []
                sp_item = item.name.split("_")
                temp_list.append(sp_item[0])
                temp_list.append(sp_item[1])
                temp_list.append(sp_item[2])
                new_name = "_".join(temp_list) # ex: NH_T0-SC_B5
                print("sample Name: ", new_name)
                sample_list.append(new_name)
    
    print("cwd: ", os.getcwd())
    os.chdir(back_path)
    print("cwd after chdir: ", os.getcwd()) 
               
    return sample_list

#call funciton
sample_list_result = sample_list_prep(args.a, args.d)
print("structure of sample_list_result: ", sample_list_result) #for debug


print("moving to next function")


#Step 2: Run Hisat2 - loop through the list you just created - make sample dirs and run the 
        #hisat2 code and move back to the main dir 

def run_initial_hisat2(ali_num,trim_dir,next_dir):
    os.mkdir("../{0}".format(next_dir))
    for dir_name in sample_list_result:
        os.mkdir("{0}".format(dir_name))
        os.chdir(dir_name)
        print ("cwd: ", os.getcwd())
        
        #Run hisat2
        result = subprocess.run("hisat2 -q -p 12 -x ../Nematostella_genome -1 ../{0}/{1}_filtered_1P.fq.gz -2 ../{0}/{1}_filtered_2P.fq.gz -S Nematostella_genome_ali_{2}_{1}.sam".format(trim_dir,dir_name,ali_num), shell=True)
        print("hisat2 complete on {0}".format(dir_name))
        print("cwd after hisat2 b4 mv: ", os.getcwd())
        
        #make new dir for next round of mapping and move sam alignement file into that dir
        os.mkdir("../../{0}/{1}".format(next_dir,dir_name))
        result = subprocess.run("mv Nematostella_genome_ali_{0}_{1}.sam ../../{2}/{1}".format(ali_num,dir_name,next_dir), shell=True)
        print("file move complete")
        os.chdir("..")
    print("All hisat2 is complete")

#call run_fastqc function
result_1 = run_initial_hisat2(args.b, args.a, args.c)
 
#Run in 1st_ali dir
#./3.A_initial_hisat2_v2.py -a ../../2_trimmomatic -b 1 -c 2nd_ali -d ../3_HISAT2/1st_ali_2-20-24/

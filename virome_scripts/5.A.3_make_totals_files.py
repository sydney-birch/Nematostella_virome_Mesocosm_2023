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
parser.add_argument("-a", help="path to input dir that has all fastq files (in .fq or .fastq), should have paired1 and paired2 files")
parser.add_argument("-b", help="aligned or filtered (aligned = aligned_NH_T0-ME_B5_1_ali_paired1.fastq 0-6 items) or (filtered = filtered_unaligned_NH_T0-NS_B1_4_ali_paired1.fq 0-7 items)")
parser.add_argument("-c", help="path from fastq files back to start dir")
parser.add_argument("-d", help="name ID for final total file ({0}-total_R1.fastq)")

args = parser.parse_args()

# Plan of attack
	#First make 2 list of files - iterate through input dir with fastq files - use it to make two files where the R1s are in one list and R2 in the other (same locations)

#aligned_NH_T0-ME_B5_1_ali_paired1.fastq
#   0     1   2    3 4  5   6
   
#filtered_unaligned_NH_T0-NS_B1_4_ali_paired1.fq
#   0         1      2  3    4  5  6    7  

	#Second - go through and concatenate 

R1_files = []
R2_files = []

#Step 1 make lists of R1 and R2 files to be concatenated 
def file_prep(input_dir, arg_flag):
    if arg_flag == "aligned":
        os.chdir(input_dir)

        for item in os.scandir():
            #print("iterating on item: ", item) #for debug
            if "paired1" in item.name:
                R1_name = item.name
                #print("R1 Name: ", R1_name)
                R1_files.append(R1_name)
            
                temp_list= []
                sp_name = item.name.split("_")
                #print("split name: ", sp_name) 
                #temp_list.append("filtered")
                temp_list.append(sp_name[0])
                temp_list.append(sp_name[1])
                temp_list.append(sp_name[2])
                temp_list.append(sp_name[3])
                temp_list.append(sp_name[4])
                temp_list.append(sp_name[5])
                temp_list.append("paired2.fq.gz")
                R2_name = "_".join(temp_list)
                #print("New R2 Name: ", R2_name)
                R2_files.append(R2_name)
                
    elif arg_flag == "filtered":
        os.chdir(input_dir)

        for item in os.scandir():
            #print("iterating on item: ", item) #for debug
            if "paired1" in item.name:
                R1_name = item.name
                #print("R1 Name: ", R1_name)
                R1_files.append(R1_name)
            
                temp_list= []
                sp_name = item.name.split("_")
                #print("split name: ", sp_name) 
                #temp_list.append("filtered")
                temp_list.append(sp_name[0])
                temp_list.append(sp_name[1])
                temp_list.append(sp_name[2])
                temp_list.append(sp_name[3])
                temp_list.append(sp_name[4])
                temp_list.append(sp_name[5])
                temp_list.append(sp_name[6])
                temp_list.append("paired2.fq.gz")
                R2_name = "_".join(temp_list)
                #print("New R2 Name: ", R2_name)
                R2_files.append(R2_name)

#call funciton
result = file_prep(args.a,args.b)
print("structure of R1 list: ", R1_files) #for debug
print("structure of R2 list: ", R2_files) #for debug


#Step 2) Concatenate the winning reads into Total R1 and Total R2 files (zipped) 
def make_totals_files(assembly_dir_path, name_id, R1_list, R2_list):
    #os.chdir("{0}".format(assembly_dir_path))
    #print("cwd after cd: ", os.getcwd()) #for debug
    
    #open two zipped files to write to - total_R1 and total_R2
    with gzip.open("{0}-total_R1.fastq".format(name_id), "wb") as out_handel_1:
        with gzip.open("{0}-total_R2.fastq".format(name_id), "wb") as out_handel_2:

            for idx in range(len(R1_list)): 
                match_name_1 = R1_list[idx]
                print("idx: ", idx)
                print("match name 1: ", match_name_1)
                match_name_2 = R2_list[idx]
                print("match name 2: ", match_name_2)

                #write total r1 file - read in each line and write out each line
                with gzip.open("{0}".format(match_name_1), "rb") as in_handel_1:
                    print("writing for file:", match_name_1)
                    for line in in_handel_1:
                        out_handel_1.write(line)
                        
                #write total r2 file - read in each line and write out each line
                with gzip.open("{0}".format(match_name_2), "rb") as in_handel_2:
                    print("writing for file", match_name_2)
                    for line in in_handel_2:
                        out_handel_2.write(line)
                    
    result = subprocess.run("mv {0}-total_R* {1}".format(name_id,assembly_dir_path), shell=True)
    print("total R1 and R2 files have been made and moved to assembly dir!!")

#Call function 
make_totals_files(args.c, args.d, R1_files, R2_files)


#how to execute:
#./5.A_make_totals_files.py -a ../4_RNA_filt/mapped_fastq_files -b aligned -c ../../5_assemblies/

#./5.A_make_totals_files.py -a ../4_RNA_filt/filtered_unmapped_fastq_files -b filtered -c ../../5_assemblies/


#! /usr/bin/env python3

#import modules 
import argparse
import subprocess
import os
import gzip


#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
#parser.add_argument("-a", help="name of output dir")
#parser.add_argument("-b", help="name of file to index for mapping")
parser.add_argument("-c", help="name of dir that has input fastq files")
#parser.add_argument("-d", help="path back from fastq dir to 1st ali dir")


args = parser.parse_args()


#This script loops through a given directory of fastq files and runs spades rnaviral - will create 
#93 assemblies - run on orion for 696 hours = 29 days


#filtered_unaligned_SC_T96-NC.M_B2_4_ali_paired2.fq.gz 
#    0        1      2     3     4  5 6     7
 
#Step 1: Run spades by looping through each file in fastq dir (input dir)  
 
def run_spades(input_dir):
    os.chdir(input_dir)
    print("current dir: ", os.getcwd())
    for item in os.scandir():
        if "paired1.fq.gz" in item.name:
            print("file name currently on: ", item.name)
 
            #get new name - rename without all the extra info stuff:
            temp_list= []
            old_name_list=[]
            sp_name = item.name.split("_")
            print("split name: ", sp_name) 
            temp_list.append(sp_name[2])
            temp_list.append(sp_name[3])
            temp_list.append(sp_name[4])

            old_name_list.append(sp_name[0])
            old_name_list.append(sp_name[1])
            old_name_list.append(sp_name[2])
            old_name_list.append(sp_name[3])
            old_name_list.append(sp_name[4])
            old_name_list.append(sp_name[5])
            old_name_list.append(sp_name[6])
            
            new_name = "_".join(temp_list)
            print("New Name: ", new_name)
            
            old_name = "_".join(old_name_list)
            print("Old Name: ", old_name)
        
            #Run spades - hard code output dir (5_assemblies) 
            result = subprocess.run("spades.py --rnaviral -1 {0}_paired1.fq.gz -2 {0}_paired2.fq.gz -o ../../5_assemblies/{1}_spades_output".format(old_name, new_name), shell=True)
            print("{0} assembly completed!!".format(new_name))
        
#call function
spades_result = run_spades(args.c)

#./5.A_run_spades.py -c ../4_RNA_filt/rRNA_filtered_unmapped_fastq_files

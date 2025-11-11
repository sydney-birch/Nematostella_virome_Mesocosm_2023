#! /usr/bin/env python3

#import modules 
import argparse
import subprocess
import os
import gzip


#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("-a", help="file that has the names of the admera names and my names to make dictionary")
parser.add_argument("-b", help="directory to navigate to with fastq files")

args = parser.parse_args()

# renames the fastq files from new sequencing from admera 

#make a dictionary of the term from admera to the name that I want {key:Value} {admera name:My name} 
#iterate through each file in the dir 
#split the name on _ --> search in dictionary for {0} --> join with dictionary value (my name), {spots: 1,2,3,4} with _

#  21081FL-06-02-11_S27_L002_R1_001.fastq.gz


#Step 1: make dictionary from input file of file names to change 

names_db = {}

with open(args.a, "r") as in_handle:
    for line in in_handle:
        line = line.rstrip()
        sp_line = line.split("\t")
        #print("my name: ", sp_line[0])
        #print("name to change: ", sp_line[1])
        
        names_db.setdefault(sp_line[1],sp_line[0])
        
    print("dictionary: ", names_db)

#Step 2: iterate through each file in the directory and change the names    

def name_change (dir_name, name_db):
    os.chdir(dir_name)
    for item in os.scandir():
        print("iterating on item: ", item) #for debug
        if "Undetermined" in item.name:
            continue

        if item.is_file():
            sample_name = item.name
            print("sample Name: ", sample_name)
            
            sp_name = sample_name.split("_")
            print ("split name: ", sp_name)
            
            print("value of key: ", name_db[sp_name[0]])
            temp_list = []
            temp_list.append(name_db[sp_name[0]])
            temp_list.append(sp_name[1])
            temp_list.append(sp_name[2])
            temp_list.append(sp_name[3])
            temp_list.append(sp_name[4])
            print("temp list: ", temp_list)
            
            new_name = "_".join(temp_list)
            print("New Name: ", new_name)
            
            print ("previous name: ", sample_name)
            result = subprocess.run("mv {0} {1}".format(sample_name, new_name), shell=True)
            print("name successfully changed: ", new_name)
            
#call run_fastqc function
result = name_change(args.b, names_db) 

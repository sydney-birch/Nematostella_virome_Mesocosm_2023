#! /usr/bin/env python3

#import modules 
import argparse
import subprocess
import os
import gzip


#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("-a", help="name of output dir")
parser.add_argument("-b", help="name of file to index for mapping")
parser.add_argument("-c", help="path to dir that has input fastq files")
parser.add_argument("-d", help="flag for T0 or T96")

args = parser.parse_args()

#Step 1: make new dir for output (filtered data) and run index

def prep_work(output_dir,index_file,time_point_flag):
    #os.mkdir("{0}".format(output_dir)) #comment out after first run in dir
    if time_point_flag == "T0":
        print("T0 index is being made")
        #result = subprocess.run("hisat2-build {0} T0_index".format(index_file), shell=True) #already made
    elif time_point_flag =="T96":    
        print("T96 index is being made")
        #result = subprocess.run("hisat2-build {0} T96_index".format(index_file), shell=True) #already made
 
#call function
prep_work_result = prep_work(args.a, args.b, args.d)

 
#Step 2: Run hisat2 by looping through each file  
 
def run_hisat2(input_dir,output_dir,time_point_flag):
    os.chdir(input_dir)
    if time_point_flag == "T0":
        for item in os.scandir():
            if "Undetermined" in item.name:
                print("this was an undetermined - skipped: ", item.name)
                continue
                
            elif "R1" in item.name:
                if "T0" in item.name: 
                    print("T0 file currently on: ", item.name)
 
                    #get new name - rename without all the extra info stuff:
                    temp_list= []
                    temp_2_list= []
                    sp_name = item.name.split("_")
                    print("split name: ", sp_name) 
                
                    temp_list.append(sp_name[0])
                    temp_list.append(sp_name[1])
                    temp_list.append(sp_name[2])

                    new_name = "_".join(temp_list)
                    print("New Name: ", new_name)
                    #make full name 
                    temp_2_list.append(sp_name[0])
                    temp_2_list.append(sp_name[1])
                    temp_2_list.append(sp_name[2])
                    temp_2_list.append(sp_name[3])
                    temp_2_list.append(sp_name[4])
                    
                    full_name = "_".join(temp_2_list)
                    print("New Name: ", full_name)

                    #Run hisat2 on T0
                    result = subprocess.run("hisat2 -q -p 12 -x ../6_BLAST_unmapped/T0_mapping/T0_index -1 {0}_R1_001.fastq.gz -2 {0}_R2_001.fastq.gz -S ../6_BLAST_unmapped/T0_mapping/{1}/T0_RAW_mapped_{2}.sam".format(full_name, output_dir,new_name), shell=True)
    
    elif time_point_flag == "T96":
        for item in os.scandir():
            if "Undetermined" in item.name:
                print("this was an undetermined - skipped: ", item.name)
                continue
                
            elif "R1" in item.name:
                if "T96" in item.name: 
                    print("T96 file currently on: ", item.name)
 
                    #get new name - rename without all the extra info stuff:
                    temp_list= []
                    temp_2_list= []
                    sp_name = item.name.split("_")
                    print("split name: ", sp_name) 
                
                    temp_list.append(sp_name[0])
                    temp_list.append(sp_name[1])
                    temp_list.append(sp_name[2])

                    new_name = "_".join(temp_list)
                    print("New Name: ", new_name)
                    
                    #make full name 
                    temp_2_list.append(sp_name[0])
                    temp_2_list.append(sp_name[1])
                    temp_2_list.append(sp_name[2])
                    temp_2_list.append(sp_name[3])
                    temp_2_list.append(sp_name[4])
                    
                    full_name = "_".join(temp_2_list)
                    print("New Name: ", full_name)
        
                    #Run hisat2 on T96
                    result = subprocess.run("hisat2 -q -p 12 -x ../6_BLAST_unmapped/T96_mapping/T96_index -1 {0}_R1_001.fastq.gz -2 {0}_R2_001.fastq.gz -S ../6_BLAST_unmapped/T96_mapping/{1}/T96_RAW_mapped_{0}.sam".format(full_name, output_dir,new_name), shell=True)
          
#call function
hisat2_result = run_hisat2(args.c, args.a, args.d)


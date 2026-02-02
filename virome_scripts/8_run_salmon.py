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
parser.add_argument("-a", help="path to fastq files")
parser.add_argument("-b", help="path from fastq to salmon dir)")
parser.add_argument("-c", help="filtered = 7 spots, aligned = 6 spots")

args = parser.parse_args()


#Step 1: prep files - make a list of your samples which each have their own dir (remove any unwanted dirs)
 
# filtered_unaligned_NH_T14-ME_B1_4_ali_paired1.fq
#     0         1     2    3    4 5  6    7
     
#aligned_NH_T0-ME_B4_1_ali_paired2.fastq
#   0     1   2    3 4  5    6  
 
def sample_list_prep(fastq_dir,back_dir,flag):
    sample_list = []
    os.chdir(fastq_dir)
    if flag == "filtered":
        for item in os.scandir():
            #print("iterating on item: ", item) #for debug
            if "paired1" in item.name:
                print("file name currently on: ", item.name)
 
                #get new name - rename without all the extra info stuff:
                temp_list= []
                sp_name = item.name.split("_")
                print("split name: ", sp_name) 
                #temp_list.append("filtered")
                temp_list.append(sp_name[0])
                temp_list.append(sp_name[1])
                temp_list.append(sp_name[2])
                temp_list.append(sp_name[3])
                temp_list.append(sp_name[4])
                temp_list.append(sp_name[5])
                temp_list.append(sp_name[6])
                new_name = "_".join(temp_list)
                print("New Name: ", new_name)
                sample_list.append(new_name)
        os.chdir(back_dir)
        
    elif flag == "aligned":
        for item in os.scandir():
            #print("iterating on item: ", item) #for debug
            if "paired1" in item.name:
                print("file name currently on: ", item.name)
 
                #get new name - rename without all the extra info stuff:
                temp_list= []
                sp_name = item.name.split("_")
                print("split name: ", sp_name) 
                #temp_list.append("filtered")
                temp_list.append(sp_name[0])
                temp_list.append(sp_name[1])
                temp_list.append(sp_name[2])
                temp_list.append(sp_name[3])
                temp_list.append(sp_name[4])
                temp_list.append(sp_name[5])
                new_name = "_".join(temp_list)
                print("New Name: ", new_name)
                sample_list.append(new_name)
        os.chdir(back_dir)
    
    return sample_list

#call funciton
sample_list_result = sample_list_prep(args.a,args.b,args.c)
print("structure of sample_list_result: ", sample_list_result) #for debug



#Step 2: Run salmon

def run_salmon(fastq_dir):
    for sample in sample_list_result:

        result = subprocess.run("salmon quant -p 12 --seqBias --gcBias -i mapped_index -l A -1 {0}/{1}_paired1.fq.gz -2 {0}/{1}_paired2.fq.gz -o {1}".format(fastq_dir,sample), shell=True)
        print("finished: {0}".format(sample))

    print("Salmon is complete")

#call run_samtools function
result_1 = run_salmon(args.a)

# filtered unmapped
#./8_run_salmon.py -a ../4_RNA_filt/filtered_unmapped_fastq_files/ -b ../../8_salmon/ -c filtered

# mapped
#./8_run_salmon.py -a ../4_RNA_filt/mapped_fastq_files / -b ../../8_salmon/ -c aligned

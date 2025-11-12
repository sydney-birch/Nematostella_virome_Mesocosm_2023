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
parser.add_argument("-a", help="input dir with sam files")
#parser.add_argument("-b", help="number in sam alignement file to use in samtools(if using second alignment sam file - use 2)")

args = parser.parse_args()



#Step 1: prep files - make a list of your samples which each have their own dir 
 
def sample_list_prep (input_dir):
    sample_list = []
    os.chdir(input_dir)

    for item in os.scandir():
        #print("iterating on item: ", item) #for debug
        
        sample_name = item.name
        print("sample Name: ", sample_name)
        
        sample_list.append(sample_name)

    return sample_list

#call funciton
sample_list_result = sample_list_prep(args.a)
print("structure of sample_list_result: ", sample_list_result) #for debug



print("moving to next function")


#Step 2: Run sam tools  

def run_samtools():
    for sample in sample_list_result:
        #os.chdir(dir_name)
        #print ("cwd: ", os.getcwd())
        print("on ssample: ", sample)
        
        temp_list= []
        sp_name = sample.split("_")
        print("split name: ", sp_name) 

        temp_list.append(sp_name[0])
        temp_list.append(sp_name[1])
        temp_list.append(sp_name[4])
        temp_list.append(sp_name[5])
        temp_list.append(sp_name[6])

        new_name = "_".join(temp_list)
        print("New Name: ", new_name)
        
        #run first part - makes bam file of all MAPPED reads
        result = subprocess.run("samtools view {0} -f 0x2 -h -b -o align_reads-{1}.bam".format(sample, new_name), shell=True)
        #Run second part 
        result = subprocess.run("samtools collate -u -O align_reads-{0}.bam | \
        samtools fastq -1 {0}_paired1.fq.gz -2 {0}_paired2.fq.gz -0 /dev/null -s /dev/null -n".format(new_name), shell=True)
        print("samtools complete on {0}".format(sample))


    print("All Samtools processing is complete")

#call run_samtools function
result_1 = run_samtools()

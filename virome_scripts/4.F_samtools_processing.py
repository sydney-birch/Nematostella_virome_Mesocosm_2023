#! /usr/bin/env python3

#import modules 
import argparse
import subprocess
import os

#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
#parser.add_argument("--rem", "-r", nargs="?", help="name of dirs not to include")
#parser.add_argument("-b", help="name of dirs not to include")
parser.add_argument("-a", help="name of dir with sam files")
#parser.add_argument("-b", help="name of output dir)")

args = parser.parse_args()

## This script is for the samtools processing - it takes the sam files from HISAT2 and turns it 
## into a bam file then creates a fastq file of only sequences that did not map 

##This script should be run in the dir that was created by the hisat2 py script
	## For example if hisat2 was run in the 1st ali dir to generate the 2nd aligement - it created a 2nd_ali dir
	## Run this script in that 2nd_ali dir which has all sub dirs and within those dirs contains the sam file 
	## in this example it would start with Nematostella_genome_ali_2

#Step 1: prep files - make a list of your samples which each have their own dir (remove any unwanted dirs)
 
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


#filtered_unaligned_NH_T0-NH_B2_4_ali.sam
#0            1      2    3   4  5  6
#Step 2: Run sam tools  

def run_samtools(sample_list):

    for sample in sample_list:
        print("on sample: {0}".format(sample))
        #make new name
        temp_list= []
        sp_name = sample.split("_")
        print("split name: ", sp_name) 
        #temp_list.append("filtered")
        temp_list.append(sp_name[0])
        temp_list.append(sp_name[1])
        temp_list.append(sp_name[2])
        temp_list.append(sp_name[3])
        temp_list.append(sp_name[4])
        #temp_list.append(sp_name[5])
        temp_list.append("ali")
        new_name = "_".join(temp_list)
        print("New Name: ", new_name)
        
        #run first part
        result = subprocess.run("samtools view {0} -f 0x4 -h -b -o unalign_reads-{1}.bam".format(sample,new_name), shell=True)
        #Run second part 
        result = subprocess.run("samtools collate -u -O unalign_reads-{0}.bam | \
        samtools fastq -1 {0}_paired1.fq.gz -2 {0}_paired2.fq.gz -0 /dev/null -s /dev/null -n".format(new_name), shell=True)
        print("samtools complete on {0}".format(sample))
        

    print("All Samtools processing is complete")

#call run_samtools function
result_1 = run_samtools(sample_list_result)

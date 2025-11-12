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
parser.add_argument("-c", help="name of dir that has input fastq files")
#parser.add_argument("-d", help="path back from fastq dir to 1st ali dir")


args = parser.parse_args()


#First have script make a new dir to output filtered data (filtered_unmapped_fastq_files)

#next have script run the index on the rnadb (output it into the new dir?)

# give script dir to fastq files - loop through each file and run hisat2 

#the output name:      (Current name: unaligned_NH_T0-NH_B5_4_ali_paired1.fastq) 
	#split name on _ make a list and append 
		#first position unaligned, Second position - filtered and then the rest 

#next script sam tools
	
#unaligned_NH_T0-NH_B5_4_ali_paired1.fastq
#     0     1    2   3 4  5      6



#Step 1: make new dir for output (filtered data) and run index

def prep_work(output_dir,index_file):
    os.mkdir("{0}".format(output_dir))
    result = subprocess.run("hisat2-build {0} rRNA_db_index".format(index_file), shell=True)
 
#call function
prep_work_result = prep_work(args.a, args.b)
 
 
 
#Step 2: Run hisat2 by looping through each file  
 
def run_hisat2(input_dir,output_dir):
    os.chdir(input_dir)
    for item in os.scandir():
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
            #temp_list.append(sp_name[6])
            new_name = "_".join(temp_list)
            print("New Name: ", new_name)
        
        #Run hisat2 
        result = subprocess.run("hisat2 -q -p 12 -x ../rRNA_db_index -1 {0}_paired1.fq.gz -2 {0}_paired2.fq.gz -S ../{1}/filtered_{0}.sam".format(new_name, output_dir), shell=True)
        
#call function
hisat2_result = run_hisat2(args.c, args.a)


#./4.C_filtered_hisat2.py -a filtered_unmapped_fastq_files -b rRNAdb_ed1.fa -c unmapped_fastq_files

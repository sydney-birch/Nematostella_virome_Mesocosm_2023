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

def prep_work(output_dir,index_file,time_point_flag):
    os.mkdir("{0}".format(output_dir))
    if time_point_flag == "T0":
        print("T0 index being made")
        #result = subprocess.run("hisat2-build {0} T0_index".format(index_file), shell=True)
    elif time_point_flag =="T96":   
        print("T96 index being made") 
        #result = subprocess.run("hisat2-build {0} T96_index".format(index_file), shell=True)
 
#call function
prep_work_result = prep_work(args.a, args.b, args.d)
 
 
 
#Step 2: Run hisat2 by looping through each file  
 
def run_hisat2(input_dir,output_dir,time_point_flag):
    os.chdir(input_dir)
    if time_point_flag == "T0":
        for item in os.scandir():
            if "paired1" in item.name:
                if "T0" in item.name: 
                    print("T0 file currently on: ", item.name)
 
                    #get new name - rename without all the extra info stuff:
                    temp_list= []
                    sp_name = item.name.split("_")
                    print("split name: ", sp_name) 
                
                    temp_list.append(sp_name[0])
                    temp_list.append(sp_name[1])
                    temp_list.append(sp_name[2])
                    temp_list.append(sp_name[3])
                    temp_list.append(sp_name[4])
                    temp_list.append(sp_name[5])
                
                    new_name = "_".join(temp_list)
                    print("New Name: ", new_name)
        
                    #Run hisat2 on T0
                    result = subprocess.run("hisat2 -q -p 12 -x ../../6_BLAST_unmapped/T0_mapping/T0_index -1 {0}_ali_paired1.fq.gz -2 {0}_ali_paired2.fq.gz -S ../../6_BLAST_unmapped/T0_mapping/{1}/T0_mapped_{0}.sam".format(new_name, output_dir), shell=True)
    
    elif time_point_flag == "T96":
        for item in os.scandir():
            if "paired1" in item.name:
                if "T96" in item.name: 
                    print("T96 file currently on: ", item.name)
 
                    #get new name - rename without all the extra info stuff:
                    temp_list= []
                    sp_name = item.name.split("_")
                    print("split name: ", sp_name) 
                
                    temp_list.append(sp_name[0])
                    temp_list.append(sp_name[1])
                    temp_list.append(sp_name[2])
                    temp_list.append(sp_name[3])
                    temp_list.append(sp_name[4])
                    temp_list.append(sp_name[5])
                
                    new_name = "_".join(temp_list)
                    print("New Name: ", new_name)
        
                    #Run hisat2 on T96
                    result = subprocess.run("hisat2 -q -p 12 -x ../../6_BLAST_unmapped/T96_mapping/T96_index -1 {0}_ali_paired1.fq.gz -2 {0}_ali_paired2.fq.gz -S ../../6_BLAST_unmapped/T96_mapping/{1}/T96_mapped_{0}.sam".format(new_name, output_dir), shell=True)
          
#call function
hisat2_result = run_hisat2(args.c, args.a, args.d)



#./5.A.5_hisat2_total_assemblies.py -a T0_total_mapping_fastqs -b total-SC-T0_assembly.fa_FIX.fa -c ../../4_RNA_filt/rRNA_filtered_unmapped_fastq_files/ -d T0
#./5.A.5_hisat2_total_assemblies.py -a T96_total_mapping_fastqs -b total-SC-T96_assembly.fa_FIX.fa -c ../../4_RNA_filt/rRNA_filtered_unmapped_fastq_files/ -d T96

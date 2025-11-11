#! /usr/bin/env python3

#import modules 
import argparse
import subprocess
import os
import gzip


#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("-a", help="path to raw reads from trimm dir")
parser.add_argument("-b", help="path from raw reads dir back to trimmomatic output dir")
#parser.add_argument("-a", help="path to fastqc output")

args = parser.parse_args()

#iterate through dir with fastq files - look for R1 - make a dictionary R1 and R2??
#iterate through collection (dictionary) - run trimmommatic that way 

#Run in trimmomatic dir - output files into trimommatic dir


#Step 1: prep files - make a list of your samples which each have their own dir (remove any unwanted dirs)
 
# NH_T14-MA_B1_S25_L002_R1_001.fastq.gz
# 0    1     2  3   4    5  6 
 
def sample_db_prep (input_path, back_path):
    sample_db = {}
    
    #print("input path: ", input_path)
    os.chdir(input_path)
    
    for item in os.scandir():
        #print("iterating on item: ", item) #for debug
        if item.is_file():
            if "SC_" in item.name: 
                if "_R1_" in item.name:
                    #make r1 file object
                    r1_file = item.name
                    
                    #make temp list - add R2 as value to db 
                    temp_list= []
                    sp_item = item.name.split("_")
                    temp_list.append(sp_item[0])
                    temp_list.append(sp_item[1])
                    temp_list.append(sp_item[2])
                    temp_list.append(sp_item[3])
                    temp_list.append(sp_item[4])
                    temp_list.append("R2")
                    temp_list.append(sp_item[6])
                    #print("temp list: ", temp_list)
                
                    r2_file = "_".join(temp_list)
                
                    print("R2 Name: ", r2_file)
                    print("R1 Name: ", r1_file)
                    #Add R1 File as key to db and R2 as value
                    sample_db.setdefault(r1_file, r2_file)
    
    #print("completed db: ", sample_db)
    print("cwd: ", os.getcwd())
    os.chdir(back_path)
    print("cwd after chdir: ", os.getcwd())
    
            
    return sample_db

#call funciton
sample_db_result = sample_db_prep(args.a, args.b)
print("structure of sample_db_result: ", sample_db_result) #for debug


print("Moving to next function to run trimmomatic")


#Step 2: Run Trimmomatic - loop through the list you just created - and run trimmomatic on 
#both the R1 and R2 files - output the info to the trimmomatic dir - pathway given (args.a)

def run_trimmomatic(names_db,input_path):
    for sample_key in names_db:
        print("key, should be R1: ", sample_key)
        print("value, should be R2: ", names_db[sample_key])
        #get new name - rename without all the extra info stuff:
        temp_list= []
        sp_name = sample_key.split("_")
        #print("split name: ", sp_name)
        temp_list.append(sp_name[0])
        temp_list.append(sp_name[1])
        temp_list.append(sp_name[2])
        new_name = "_".join(temp_list)
        print("New Name: ", new_name)
        
        #print("path to fastqs: ", input_path)
        
        #Run trimmomatic
        result = subprocess.run("java -jar /apps/pkg/trimmomatic/0.39/trimmomatic-0.39.jar PE -phred33 -threads 4 {0}/{1} {0}/{2} -baseout {3}_filtered.fq.gz ILLUMINACLIP:/apps/pkg/trimmomatic/0.39/adapters/TruSeq3-PE-2.fa:1:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36".format(input_path, sample_key,names_db[sample_key], new_name), shell=True)
        
    print("All trimmomatic is complete")

#call run_fastqc function
result_1 = run_trimmomatic(sample_db_result, arg

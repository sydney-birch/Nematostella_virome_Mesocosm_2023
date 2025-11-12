#! /usr/bin/env python3

#import modules 
import argparse
import subprocess
import os
import gzip


#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="name of dir that has input fastq files")

args = parser.parse_args()


#filtered_unaligned_SC_T96-NC.M_B2_4_ali_paired2.fq.gz 
#    0        1      2     3     4  5 6     7
 

#make a dictionary - iterate through input dir split the name take the [2]SC and [3]T0-FL.F join - make key
#add that and [4]B5 as the value 

#then, iterate through dict - get length of list (number of reps)
#if statement - if length = 3 then run spades command - index list with three
# if length == 4 then run spades command - index list with three

#Step 1: setup structure of samples - fill in names dictionary 
names_db = {}
 
def make_names_db(input_dir):
    os.chdir(input_dir)
    print("current dir:", os.getcwd())
    for item in os.scandir():
        if "paired1.fq.gz" in item.name:
            #print("file name currently on: ", item.name)
 
            #make temp lists to get names for the keys and values to add:
            key_list= []
            value_list=[]
            
            sp_name = item.name.split("_")
            #print("split name: ", sp_name) 
            
            #populate key list
            key_list.append(sp_name[2])
            key_list.append(sp_name[3])

            #populate value list
            value_list.append(sp_name[2])
            value_list.append(sp_name[3])
            value_list.append(sp_name[4])

            #get key and value names
            key_name = "_".join(key_list)
            #print("Key Name: ", key_name)
            
            value_name = "_".join(value_list)
            #print("Value Name: ", value_name)
            
            #add to dictionary
            names_db.setdefault(key_name,[]) #set key name and make list as value
            names_db[key_name].append(value_name) #add value name to value list 
        
    print("Final names dictionary: ", names_db)
    return names_db
    
#call function
names_db_result = make_names_db(args.i)


#Step 2: Run spades 
def run_spades():
    print("current dir:", os.getcwd())
    
    #iterate through dictionary and run spades
    for key in names_db:
        #get length of list
        num_samples = len(names_db[key])
        print("number of samples: ", num_samples)
        #run spades
        if num_samples == 3:
            print("in run 3")
            result = subprocess.run("spades.py --rnaviral --pe1-1 filtered_unaligned_{0}_4_ali_paired1.fq.gz --pe2-1 filtered_unaligned_{1}_4_ali_paired1.fq.gz --pe3-1 filtered_unaligned_{2}_4_ali_paired1.fq.gz --pe1-2 filtered_unaligned_{0}_4_ali_paired2.fq.gz --pe2-2 filtered_unaligned_{1}_4_ali_paired2.fq.gz --pe3-2 filtered_unaligned_{2}_4_ali_paired2.fq.gz -o ../../5_assemblies/{3}_output".format(names_db[key][0],names_db[key][1],names_db[key][2],key), shell=True)
            print("{0} assembly completed!!".format(key)) 
              
        elif num_samples == 4:
            print("in run 4")
            result = subprocess.run("spades.py --rnaviral --pe1-1 filtered_unaligned_{0}_4_ali_paired1.fq.gz --pe2-1 filtered_unaligned_{1}_4_ali_paired1.fq.gz --pe3-1 filtered_unaligned_{2}_4_ali_paired1.fq.gz --pe4-1 filtered_unaligned_{3}_4_ali_paired1.fq.gz --pe1-2 filtered_unaligned_{0}_4_ali_paired2.fq.gz --pe2-2 filtered_unaligned_{1}_4_ali_paired2.fq.gz --pe3-2 filtered_unaligned_{2}_4_ali_paired2.fq.gz --pe4-2 filtered_unaligned_{3}_4_ali_paired2.fq.gz -o ../../5_assemblies/{4}_output".format(names_db[key][0],names_db[key][1],names_db[key][2],names_db[key][3],key), shell=True)
            print("{0} assembly completed!!".format(key))
    print("All assemblies completed")

#call function
result = run_spades()

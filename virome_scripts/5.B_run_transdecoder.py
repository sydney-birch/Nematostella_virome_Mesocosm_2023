#! /usr/bin/env python3

#import modules 
import argparse
import subprocess
import os
import gzip

#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()

parser.add_argument("-c", help="name of output dir for final nuc fastas")
parser.add_argument("-b", help="name of output dir for final prot fastas")

args = parser.parse_args()

#Run this in the 5_assembly dir after all assemblies have been created. Each assembly has
#its own dir - this will go into each dir - rename the fasta file, run it thru transdecoder, 
#rename the prot fasta, and move both nuc and prot fastas to respective dirs outside of the 
#spades sample dir_name

 #SC_T0-ME.M_B2_spades_output
 #0     1     2


#Step 1: prep files - make a list of your samples which each have their own dir 
def sample_list_prep ():
    sample_list = []

    for item in os.scandir():
        #print("iterating on item: ", item) #for debug
        if item.is_dir():
    
            sample_name = item.name
            print("sample Name: ", sample_name)
            
            sample_list.append(sample_name)

    return sample_list

#call funciton
sample_list_result = sample_list_prep()
print("structure of sample_list_result: ", sample_list_result) #for debug


print("moving to next function")


#Step 2: Run transdecoder - loop through the list you just created - go into each dir 
		#change name of nuc fasta to sample name, run transdecoder, change name of prot fasta
		#and move nuc fasta to new nuc fasta dir and prot fasta to prot fasta dir 
def run_transdecoder(nuc_dir,prot_dir):
    os.mkdir("{0}".format(nuc_dir))
    os.mkdir("{0}".format(prot_dir))
    for dir_name in sample_list_result:
        os.chdir(dir_name)
        print ("cwd: ", os.getcwd())
        
        #Get sample part of name to rename fasta file
        sp_name = dir_name.split("_")
        print ("split name: ", sp_name)
            
        temp_list = []
        temp_list.append(sp_name[0])
        temp_list.append(sp_name[1])
        temp_list.append(sp_name[2])
        print("temp list: ", temp_list)
            
        new_name = "_".join(temp_list)
        print("New Name: ", new_name)
        
        #Rename scaffolds.fasta to actual name
        result = subprocess.run("mv scaffolds.fasta {0}-nuc.fa".format(new_name), shell=True)    
        
        #Run transdecoder
        result = subprocess.run("TransDecoder.LongOrfs -t {0}-nuc.fa".format(new_name), shell=True)
        print("transdecoder complete on {0}-nuc.fa".format(dir_name))
        
        #move nuc fasta to fastas dir
        result = subprocess.run("cp {0}-nuc.fa ../{1}".format(new_name,nuc_dir), shell=True) 
        
        #rename and move prot fasta to prot_fastas dir
        result = subprocess.run("mv {0}-nuc.fa.transdecoder_dir/longest_orfs.pep {0}-nuc.fa.transdecoder_dir/{0}-prot.fa".format(new_name), shell=True)
        result = subprocess.run("cp {0}-nuc.fa.transdecoder_dir/{0}-prot.fa ../{1}".format(new_name,prot_dir), shell=True)
    
        print("cwd after transdecoder: ", os.getcwd())

        os.chdir("..")
    print("All translation and renaming is complete")

#call run_transdecoder function
result_1 = run_transdecoder(args.c, args.b)

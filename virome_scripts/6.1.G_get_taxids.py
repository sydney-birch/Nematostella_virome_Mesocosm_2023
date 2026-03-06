#! /usr/bin/env python3

#import modules 
import argparse
import subprocess
import os


#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input, the ACCID file")
parser.add_argument("-d", help="the Refseq database")
parser.add_argument("-b", help="base name to use in temp and final name")

args = parser.parse_args()

#first edit down the basename - split the current name on a "-", add to a temp list and join with "-"
temp_list=[]

split_name= args.b.split("-")
temp_list.append(split_name[0])
temp_list.append(split_name[1])

#join the list to make the shorten base name
base_name="-".join(temp_list)
print("base name: ", base_name)
 
#Open the accid input file, iterate through each line and search the accid in the refseq database - append search result to new file
with open(args.i, "r") as in_handle:
    for accid in in_handle:
        accid = accid.strip()
        refseq = args.d
        #print("current accid: ", accid)
        result = subprocess.run("grep {0} {1} >> {2}-taxids_accids".format(accid, refseq, base_name), shell=True)

#open the new file created of taxid and accids and open a new output file that will have just taxids
with open("{0}-taxids_accids".format(base_name), "r") as in_handle_2:
    with open("{0}-taxids".format(base_name), "w") as out_handle: 
        #iterate through each line, split on tab, and then write the taxid to new output file 
        for line in in_handle_2:
            sp_line = line.split("\t")
            out_handle.write(sp_line[0]+"\n")
        

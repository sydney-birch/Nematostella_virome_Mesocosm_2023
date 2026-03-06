#! /usr/bin/env python3

#import modules 
import argparse
import subprocess
import os


#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input, the file with taxids one per line")

args = parser.parse_args()

#first edit down the basename - split the current name on a "-", add to a temp list and join with "-"
temp_list=[]

print("args.i: ", args.i)
split_name= args.i.split("/")
sp_2_name = split_name[1].split("-")
temp_list.append(sp_2_name[0])
temp_list.append(sp_2_name[1])

#join the list to make the shorten base name
base_name="-".join(temp_list)
print("base name: ", base_name)
 
#Run taxonkit line
#code_line = "cat {0} | ./taxonkit reformat --data-dir TAXONKIT_DB -I 1 -F -P -f \"{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}\" >> {1}_linage.txt".format(args.i,base_name)
#print("code line: ", code_line) 
#result = subprocess.run("cat {0} | ./taxonkit reformat --data-dir TAXONKIT_DB -I 1 -F -P -f \"{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}\" >> {1}_linage.txt".format(args.i,base_name), shell=True)
result = subprocess.run("cat {0} | ./taxonkit reformat --data-dir TAXONKIT_DB -I 1 -F -P -f \"{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}\\t{{t}}\" >> {1}_linage.txt".format(args.i,base_name), shell=True)

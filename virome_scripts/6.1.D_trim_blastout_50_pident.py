#! /usr/bin/env python3

#import modules 
import argparse
import subprocess
import os


#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input the blastout file")

args = parser.parse_args()

with open(args.i, "r") as in_handle:
    with open("{0}_50_pi".format(args.i), "w") as out_handle: 
        for line in in_handle:
            sp_line = line.split("\t")  
            pident = float(sp_line[2])
            if pident >= 50: 
                out_handle.write(line)

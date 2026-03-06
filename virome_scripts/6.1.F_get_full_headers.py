#! /usr/bin/env python3

#import modules 
import argparse
import subprocess
import os


#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("-a", help="input file with headers(output from get_accessions)")
parser.add_argument("-f", help="input file of the transcriptome with full header names (contains spaces) ")
parser.add_argument("-o", help="name of output file")


args = parser.parse_args()

#First iterate through file with headers only (from get accessions script) - populate the keys with accession ids 
headers_dict={}

def make_dict(headers_file):

    with open("{0}".format(headers_file), "r") as in_handel_1:
        for header in in_handel_1:
            header = header.rstrip()
            #print("current header: {0}".format(header))
            #add header to dictionary
            headers_dict.setdefault(header)
        #print ("header dictionary: {0}".format(headers_dict))
            
#call funciton
result = make_dict(args.a)   


#next open up the transcriptome with full header names - iterate through each line in file and iterate through dict - if match then add as value in dict
def get_headers(fasta_file,output_file):
    #populate dict with new headers
    with open("{0}".format(fasta_file), "r") as in_handel_1:
        for line in in_handel_1:
            for header_key in headers_dict:
                if line.startswith(">{0}".format(header_key)):
                    headers_dict[header_key] = line.rstrip()
        
        #print ("populated header dictionary: {0}".format(headers_dict))
   
    #write new headers to output file    
    with open("{0}".format(output_file), "w") as out_handel_1:
        for header_key in headers_dict:
            correct_header = headers_dict[header_key].replace(">"," " )
            correct_header = correct_header.strip()
            #print("correct header: {0}".format(correct_header))
            out_handel_1.write("{0}\n".format(correct_header))
            
            
#call funciton
result = get_headers(args.f,args.o) 

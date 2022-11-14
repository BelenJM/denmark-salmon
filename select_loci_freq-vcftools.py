#!/usr/bin/env python

# Overall aim: extract the information from .frq.count
# from vcftools, to create the input for signasel (Hubert et al.)

# Opens the txt file, and copies the file into a common txt file
# Opens the txt file and modifies the input file
# from __future__ import print_function

import sys
import os
import gzip

#from Bio.SeqIO.QualityIO import FastqGeneralIterator

input_file1 = sys.argv[1] # input one is the loci to search
input_file2 = sys.argv[2] # input two should be the name of the river to use for files to search from
input_file3 = sys.argv[3] # input three is the output name

#allele = {}

directory = r'C:\Users\bmen\Desktop\salmon\allele_freq_ALL_vcftools_bayescan23092020'
# "a" makes it appending
# any text that gets written in the file
#common_file = open(output_file1, "a")

#name_ID = input_file1.split('_')[1:4]

def str_join(*args):
    return ''.join(map(str, args))
    
# with open(input_file2, "rU") as second_file:
with open(input_file1, "r") as loci, open(input_file3, "a") as out: 
    index = 0
    for l in loci:
        print(l)
        line = l.strip("\n").split("\t") # removes the \n and then splits the line into the elements separated by \t
        index += 1
        print(line[0])
        #new_minor = []
        #minor = []
        for filename in os.listdir(directory):
            filename_full = str_join("5e5.", input_file2)
            if filename.startswith(filename_full):
                file = os.path.join(filename)
                print("The name is:", file)
                
                # printing Ne and generation
                if file == "5e5.SK1913.frq":
                    year = "1913"
                    Ne = float(1729.9)
                    gen = float(1)
                if file == "5e5.SK1930early.frq":
                    year = "1930early"
                    Ne = float(76.9)
                    gen = float(6)
                if file == "5e5.SK1930late.frq":
                    year = "1930late"
                    Ne = float(286)
                    gen = float(7)
                if file == "5e5.SK1945.frq":
                    year = "1945"
                    Ne = float(40.4)
                    gen = float(8)
                if file == "5e5.SK1955.frq":
                    year = "1955"
                    Ne = float(132.2)
                    gen = float(11)
                if file == "5e5.SK1990.frq":
                    year = "1990"
                    Ne = float(18.9)
                    gen = float(20)
                if file == "5e5.SK1990_minusAdmix.frq":
                    year = "1990_minusAdmix"
                    Ne = float(18.9)
                    gen = float(20)
                if file == "5e5.SK1999.frq":
                    year = "1999"
                    Ne = float(68.6)
                    gen = float(22)
                if file == "5e5.SK1999_minusAdmix.frq":
                    year = "1999_minusAdmix"
                    Ne = float(68.6)
                    gen = float(22)
                if file == "5e5.SK2008.frq":
                    year = "2008"
                    Ne = float(129.3)
                    gen = float(25)
                if file == "5e5.SK2015.frq":
                    year = "2015"
                    Ne = float(273.4)
                    gen = float(27)
                
                with open(file, "r") as AC_file: # opening the allele counts file as text using "rt"
                    next(AC_file)

                    for AC in AC_file:
                        AC_line = AC.strip("\n").split("\t") # removes the \n and then splits the line into the elements separated by \t
                        #print(type(maf_line))
                        #f=[i.strip() for i in maf_line]
                        #print(maf_line)[0]
                        #for i in maf_line:
                        #print(maf_line[0])
                        if AC_line[0] == line[0]:
                            if AC_line[1] == line[1]:
                                #print(maf_line[1])
                                no_ind = int(AC_line[3])/2
                                #print(no_ind)
                                allele1 = AC_line[4][0]
                                #print(maf_line[4][2:-1])
                                allele1_count = AC_line[4][2:]
                                allele2 = AC_line[5][0]
                                allele2_count = AC_line[5][2:]
                                # for pcadapt:
                                out.write("%i\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (index, year, Ne, gen, file[6:10], line[0], line[1], no_ind, allele1, allele1_count, allele2, allele2_count))
                                #minor = new_minor
                                break
                            else:
                                continue
                        else:
                            continue

out.close()
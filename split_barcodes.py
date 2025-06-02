#!/usr/bin/env python

import subprocess, sys, gzip, os, re
import mysql.connector

from glob import glob
#from time import sleep
import re
import argparse
from argparse import RawTextHelpFormatter
from datetime import datetime

# Option for setting the barcode length - can default to the length in Sierra but we want to be able to set it too.
# Option (for Jon's NEB UMI data) to take remaining sequence from I1 file and add it the read ID of read 1 files. 
# Check the TrAEL format as that's how the deduplication works
# option to pass in a csv samplesheet if we don't want to get the sample info from Sierra.

# single barcoded
# python3 ../../split_barcodes.py 20250404_AV240405_AV_B_PRG6210_PE75_04042025

# dual barcoded
# python3 ../../split_barcodes.py 20250127_AV240405_AV_B_MT6176_ET6177_SE75_27012025 lane1_NoIndex_L001_R1.fastq.gz lane1_NoIndex_L001_I1.fastq.gz lane1_NoIndex_L001_I2.fastq.gz

fhs = {}           # storing the filehandles for all output files - dictionary of filehandles where key is sample barcode
fhsR2 = {}
#paired_end = False
double_coded = False
prepath = "/bi/group/bioinf/Laura_B/python_barcode_splitting/*/"
path_from_run_folder = ""

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description = '''For demultiplexing fastq files''')
parser.add_argument('run_folder', type=str, default="", help='run folder name')
parser.add_argument('--sample_sheet', type=str, default="", help='[Optional] sample sheet will be pulled from Sierra by default')
parser.add_argument('--lane_number', type=str, default="1", help='Lane number on flow cell i.e. L001 will be 1. Default: 1')

#parser.add_argument('--verbose', default=False, action='store_true', help='verbose processing')
parser.add_argument('--i1_umi', default=False, action='store_true', help='If UMI is present after the barcode in the I1 file. Set barcode length if this is specified')
parser.add_argument('--barcode_length', type=int, default=0, help='If barcode length differs from actual length of sequences in the index file(s)')

args=parser.parse_args()

run_folder = args.run_folder


def main():

    print(datetime.now())
    
    file_location = f"{prepath}{run_folder}{path_from_run_folder}"
   
    lane_number = args.lane_number

   #paired_end = False

	# open logfile
    sample = "temp_name"

    log_filename = f"{sample}_indexing.log"
    fhs["log"] = open(log_filename, mode = "w")

    try:
        expected_barcodes_list = get_expected_barcodes(run_folder, lane_number)
        expected_barcodes = expected_barcodes_list[0]
        double_coded = expected_barcodes_list[1]
        print(f"barcodes are {expected_barcodes}. \nIs this a double coded library? {double_coded}")

       # split_fastqs(R1, R2, I1, I2, expected_barcodes, double_coded)
        split_fastqs(file_location, expected_barcodes, double_coded)

    except Exception as err:
        print(f"\n !! Couldn't get expected barcodes !!  Is the run folder correct? {run_folder}")
        print(err)

    close_filehandles()
    print(datetime.now())

#----------------------------------------------
#  do the splitting
#----------------------------------------------
def split_fastqs(file_location, expected_barcodes, double_coded):
  
    #filename = R1
    #filename2 = R2

    R1 = get_R1(file_location)
    R2 = get_R2(file_location)
    I1 = get_I1(file_location)
    I2 = get_I2(file_location)

    #if (double_coded and I2 is not None):
    # dual barcoded in Sierra but couldn't find an I2 file

    if (R2 is not None):
        paired_end = True
    else:
        paired_end = False

    for key in expected_barcodes:

        new_filenameR1 = f"laneX_{key}_{expected_barcodes[key]}_R1.fastq.gz"
		#new_filename = f"{sample}_{sample_level_barcode}_{suffix}"
        open_filehandles(new_filenameR1, key)

        if paired_end:
            new_filenameR2 = f"laneXX_{key}_{expected_barcodes[key]}_R2.fastq.gz"
            open_filehandles(new_filenameR2, key)

    #with gzip.open(filename) as cf, gzip.open(filename2) as cf2:
    with gzip.open(R1) as cf, gzip.open(I1) as i1, gzip.open(I2) as i2:
    #with gzip.open(filename) as cf, gzip.open(filename2) as cf2, gzip.open(I1) as i1, gzip.open(I2) as i2:
		# unassigned_count = 0 # count the number of reads that don't match supplied barcodes
		# unpaired_count = 0 # count the number of R2 barcodes that don't match R1
        line_count = 0
        barcode = ""

        #while True:
        while line_count <= 4000: 
            readID_R1  = cf.readline().decode().strip()
            seq_R1     = cf.readline().decode().strip()
            line3_R1   = cf.readline().decode().strip()
            qual_R1    = cf.readline().decode().strip()

            shortID_R1 = readID_R1.split(" ")[0]
			
            if not qual_R1:
                break
			
            if paired_end:
                readID_R2  = cf2.readline().decode().strip()
                seq2_R2    = cf2.readline().decode().strip()
                line3_R2   = cf2.readline().decode().strip()
                qual_R2    = cf2.readline().decode().strip()
                shortID_R2 = readID2.split(" ")[0]

            readID_I1  = i1.readline().decode().strip()
            seq_I1     = i1.readline().decode().strip()
            line3_I1   = i1.readline().decode().strip()
            qual_I1    = i1.readline().decode().strip()
            shortID_I1 = readID_I1.split(" ")[0]

            if double_coded:
                readID_I2  = i2.readline().decode().strip()
                seq_I2     = i2.readline().decode().strip()
                line3_I2   = i2.readline().decode().strip()
                qual_I2    = i2.readline().decode().strip()
                shortID_I2 = readID_I2.split(" ")[0]

                barcode = f"{seq_I1}_{seq_I2}"
            
            else:
                barcode = seq_I1

            if barcode in expected_barcodes.keys():
                #print(f"Found it!! {barcode} has the name {expected_barcodes[barcode]}")
                readID_R1 = f"{readID_R1} {barcode}"

                # this one's quicker - is it because we're not writing out so many times?
                fhs[barcode].write (("\n".join([readID_R1, seq_R1, line3_R1, qual_R1]) + "\n").encode())

            #else:
                #print(f"Couldn't find this {barcode}")

            line_count += 1

        #     # I don't think that we should need to check this
            if shortID_R1 != shortID_I1:
                err_msg = f"\n!! IDs do not match for read {line_count}, exiting... !!\n"
                print(err_msg)
                fhs["log"].write(err_msg)
                exit()

def open_filehandles(fname, sample_level_barcode):
	#print (f"Opening filehandle for {sample_level_barcode} and {fname}")
	fhs[sample_level_barcode] = gzip.open (fname,mode='wb',compresslevel=3)


def close_filehandles():
	for name in fhs.keys():
		fhs[name].close() 

#---------------------
# quick barcode check
#---------------------
def get_barcodes_I1(run_folder, n_bars_to_check):

    try:
        os.chdir(f"/primary/{run_folder}/Unaligned/Project_External/Sample_lane1/")

        bc_check_cmd = f"zcat lane1_NoIndex_L001_I1.fastq.gz | head -n 400000 | awk 'NR % 4 == 2' | sort | uniq -c | sort -k 1 -n -r | head -n {n_bars_to_check} > found_barcodes_L001_I1.txt"

        try:
            subprocess.run(bc_check_cmd, shell=True, executable="/bin/bash")
            
        except Exception as err:
            print(f"\n !! Couldn't run barcode checking command !!")
            print(err)
            
    except Exception as err:
            print(f"\n !! Couldn't run barcode checking command !!")
            print(err)


#---------------------
# quick barcode check
#---------------------
def get_barcodes_I2(run_folder, n_bars_to_check):

    try:
        os.chdir(f"/primary/{run_folder}/Unaligned/Project_External/Sample_lane1/")
        bc_check_cmd = f"zcat lane1_NoIndex_L001_I2.fastq.gz | head -n 100000 | awk 'NR % 4 == 2' | sort | uniq -c | sort -k 1 -n -r | head -n {n_bars_to_check} > found_barcodes_L001_I2.txt"

        try:
            subprocess.run(bc_check_cmd, shell=True, executable="/bin/bash")
            
        except Exception as err:
            print(f"\n !! Couldn't run barcode checking command !!")
            print(err)
            
    except Exception as err:
        print(f"\n !! Couldn't run barcode checking command !!")
        print(err)


#---------------------------------------------------
# remove any unwanted characters from sample names
#---------------------------------------------------
def clean_sample_name(sample_name):    
    return re.sub(r'[^a-zA-Z0-9.\-_]+_?', '_', sample_name)


#---------------------------------------------------------
# Get the expected barcodes from sample sheet in Sierra
#---------------------------------------------------------
def get_expected_barcodes(run_folder, lane_number):

    try:
        #os.chdir(f"/primary/{run_folder}/Unaligned/Project_External/Sample_lane1/")
        double_coded = False
        barcode_dict = {}
        count = 0

        cnx = mysql.connector.connect(user='sierrauser', password='', host='bilin2.babraham.ac.uk', database='sierra')
        cursor = cnx.cursor()

        # lane_number will mostly be 1 for aviti and miseq - but we are starting to get some 2s
        #query = (f"select barcode.5_prime_barcode,barcode.3_prime_barcode,barcode.name from run,flowcell,lane,barcode WHERE run.run_folder_name = '{run_folder}' and run.flowcell_id=flowcell.id AND run.flowcell_id=lane.flowcell_id AND lane.lane_number=1 AND lane.sample_id = barcode.sample_id")
        #query = (f"select barcode.5_prime_barcode,barcode.3_prime_barcode,barcode.name from run,flowcell,lane,barcode WHERE run.run_folder_name = '{run_folder}' and run.flowcell_id=flowcell.id AND run.flowcell_id=lane.flowcell_id AND lane.lane_number='{lane_number}' AND lane.sample_id = barcode.sample_id")

        query = (
            f"select barcode.5_prime_barcode,barcode.3_prime_barcode,barcode.name from run,flowcell,lane,barcode " 
            f"WHERE run.run_folder_name = '{run_folder}' and run.flowcell_id=flowcell.id AND run.flowcell_id=lane.flowcell_id " 
            f"AND lane.lane_number='{lane_number}' AND lane.sample_id = barcode.sample_id"
        )

        cursor.execute(query)

        for (row) in cursor:

            bc1 = row[0].strip()
            bc2 = row[1].strip()
            # sample name is always the 3rd field, the 2nd is empty if it is single barcoded
            sample_name = row[2].strip()
            sample_name = clean_sample_name(sample_name)

            if count == 0:                
                if  bc2 == "":
                    print("barcode 2 is empty, assuming single coded")
                    double_coded = False
                else:
                    print("This is a double coded library.")
                    double_coded = True   

            if(double_coded):
                barcode_seq = f"{bc1}_{bc2}"
            else:
                barcode_seq = bc1
                
            barcode_dict[barcode_seq] = sample_name
            count += 1

        cnx.close()
        print(f"\nFound {count} expected barcodes in Sierra\n")
        return([barcode_dict, double_coded])
        
    except Exception as err:
        print(f"\n !! Couldn't get expected barcodes !!")
        print(err)

#------------------------------
# Locate the Read 1 fastq file
#------------------------------
def get_R1(noIndex_location):
    
    R1_location = f"{noIndex_location}/*NoIndex*R1.fastq.gz"
    R1 = glob(R1_location)

    if len(R1) == 1:
        return(R1[0])     
    else:
        print(f"!! Something went wrong with locating the R1 file here: {R1_location}, exiting... !!\n")
        exit()


#------------------------------
# Locate the Read 2 fastq file
#------------------------------
def get_R2(noIndex_location):
    
    R2_location = f"{noIndex_location}/*NoIndex*R2.fastq.gz"
    R2 = glob(R2_location)

    if len(R2) == 1:
        #paired_end = True
        return(R2[0])   
    elif len(R2) == 0:
        #print(f"Couldn't find an R2 file here: {R2_location}. \nAssuming a single barcoded library.")
        print(f"Couldn't find an R2 file, assuming a single barcoded library.")
        return(None)    
    else:
        print(f"!! Something went wrong with locating the R2 file here: {R2_location}, exiting... !!\n")
        exit()


#-------------------------------
# Locate the first barcode file
#-------------------------------
def get_I1(noIndex_location):
    
    I1_location = f"{noIndex_location}/*NoIndex*I1.fastq.gz"
    I1 = glob(I1_location)

    if len(I1) == 1:
        return(I1[0])     
    else:
        print(f"!! Something went wrong with locating the I1 file here: {I1_location}, exiting... !!\n")
        exit()


#--------------------------------
# Locate the second barcode file
#--------------------------------
def get_I2(noIndex_location):
    
    I2_location = f"{noIndex_location}/*NoIndex*I2.fastq.gz"
    I2 = glob(I2_location)

    if len(I2) == 1:
        dual_barcoded = True
        return(I2[0])   
    elif len(I2) == 0:
        print(f"Couldn't find an I2 file here: {I2_location}. \nAssuming a single barcoded library.")
        return(None)    
    else:
        print(f"!! Something went wrong with locating the I2 file here: {I2_location}, exiting... !!\n")
        exit()



if __name__ == "__main__":
    main()

#!/usr/bin/env python3

"""
Written by Devon Gregory

Last edited 1-31-23
"""

import sys
import os
import datetime
import argparse

parser = argparse.ArgumentParser(
        description='Pull GISAID samples from a fasta file that were collected in the past 90 days'
)

parser.add_argument(
    'in_file',
    type=argparse.FileType('r'),
    help='fasta file'
)

args = parser.parse_args()

cutoff = datetime.date.today() - datetime.timedelta(days=90)

passed = False
passed_dict = {}
cur_id = ""
first_date = datetime.date.today()
last_date = cutoff
for line in args.in_file:
    if line.startswith('>'):
        try:
            date = datetime.date.fromisoformat(line.split('|')[2].strip())
        except IndexError:
            continue
        except ValueError:
            try:
                year, month, day = line.split('|')[2].strip().split("-")
                if int(day) == 0:
                    day = "01"
                if int(month) == 0:
                    month = "01"
                if int(day) == 0:
                    year = "01"
                date = datetime.date.fromisoformat(f"{year}-{month}-{day}")
            except IndexError:
                continue
        if date >= cutoff:
            cur_id = line.strip()
            passed_dict[cur_id] = ""
            passed = True
            if date > last_date:
                last_date = date
            if date < first_date:
                first_date = date
            
        else:
            passed = False
            
                
            
    elif passed:
        passed_dict[cur_id] += line.strip()
    
if passed_dict:
    #with open(f"{args.in_file.name}_{first_date}--{last_date}.fasta", "w") as out_fh:    
    with open("last90.fasta", "w") as out_fh:
        for entry in passed_dict:
            out_fh.write(entry)
            out_fh.write("\n")
            out_fh.write(passed_dict[entry])
            out_fh.write("\n")
else:
    print("No entries fround for the last 90 days")


#!/usr/bin/env python3

"""
Written by Devon Gregory

Last edited 2-20-23
"""

import sys
import os
import argparse
import calendar
import time


SARS2_CDS =  {"ORF1a"    : [266, 13468],
               "ORF1b"    : [13468, 21555],
               "S"        : [21563, 25384],
               "ORF3a"    : [25393, 26220],
               "E"        : [26245, 26472],
               "M"        : [26523, 27191],
               "ORF6"     : [27202, 27387],
               "ORF7a"    : [27394, 27759],
               "ORF7b"    : [27756, 27887],
               "ORF8"     : [27894, 28259],
               "N"        : [28274, 29533],
               "ORF10"    : [29558, 29674]}


def arg_parser():
    """
    Parses the command line arguments using pythin argparse.
    returns arguments object
    """

    parser = argparse.ArgumentParser(description='process MSA fasta with first entry as a reference for vanriance information')

    parser.add_argument("in_tsv", type=argparse.FileType('r'), help="name of tsv with NextClade output")

    parser.add_argument('-m', nargs='*', type=str, dest='mut_file', action='append', help='Name of text file(s) with list of mutations of interest')
    
    parser.add_argument('-n', type=int, choices=[0, 1], dest='n', default=1, help='look for last new mutation file (-n 0 to disable)')

    parser.add_argument('-o', type=str, dest='out_file', help='Name of tabular output file')

    parser.add_argument('-p', action='store_true', dest='p', help='activate mutation search')
    
    parser.add_argument('-d', type=int, dest='days', default=7, help='number of consecutive days a mutation must be above the threshold to be tagged as of potential interest')

    parser.add_argument('-t', type=float, dest='freq_threshold', default=.2, help='frequency threshold a mutation must be above to be tagges as of potential interest')

    args = parser.parse_args()

    return args


def decode_gisaid_id(gisaid_id):
    """
    Pulls location and date information from a gisaid name.
    Parameters:
    gisaid_id - string of a gisaid virus name - str
    Functionality:
        splits appart the full gisaid name to pull the location and date
        Returns str with the location and date of the sample
    """

    try:
        loc = gisaid_id.split('/')[1]
        if loc == "env":
            loc = gisaid_id.split('/')[2]
        loc = loc.replace(" ", "_").replace("-", "_")
        if loc == "USA":
            loc += f"-{gisaid_id.split('/')[2].split('-')[0]}"

        date = gisaid_id.split('|')[2].strip()
        ident = f"{loc}|{date}"
    except IndexError:
        ident = gisaid_id

    return ident


def parse_input(tsv_file):
    """
    parses a nextclade.tsv output file to get the pango and mutation information for each sample.
    Parameters:
    tsv_file - nextclade.tsv outpt - file object
    Functionality:
        goes through each line of the file and pulls the sample information
        and collapses duplicate sequences from the same location and date
        returns dict of samples with the count of their occurance
    """
    seq_dict = {}
    first_line = tsv_file.readline()
    tsv_file.seek(0)
    if not first_line.startswith("index\tseqName\tclade\t"):
        print(f"{tsv_file.name} not recognized as nextclade output tsv")
        exit(1)
    else:
        total_reads = -1
        for line in tsv_file:
            split_line = line.split("\t")
            samp_info = ";;".join((
                                    split_line[3], # pango
                                    split_line[24], # subs
                                    split_line[25], # dels
                                    split_line[26], # ins
                                    split_line[31], # priv Rev
                                    split_line[32], # priv labeled
                                    split_line[33], # priv unlabeled
                                    split_line[29], # aa dels
                                    split_line[30], # aa ins
                                    split_line[27], # frameshift
                                    split_line[28], # aasubs
                                    ))
            total_reads += 1
            ident = decode_gisaid_id(split_line[1])
            try:
                seq_dict[ident]
            except KeyError:
                seq_dict[ident] = { samp_info : 1 }
            else:
                try:
                    seq_dict[ident][samp_info] += 1
                except KeyError:
                    seq_dict[ident][samp_info] = 1


    return seq_dict, total_reads


def get_muts(list_file):
    """
    parses a text file with a mutation of interest on each line
    Parameters:
    list_file - name of file with a list of mutation of interest - str
    Functionality:
        goes through each line of the file and pulls a mutation if one is present
        returns list of the mutations
    """
    muts = []
    with open(list_file, "r") as mut_list:
        for line in mut_list:
            if line.strip():
                muts.append(line.strip())
    return muts


def mut_list_aa_match(nt_mut_list, aa_subs):

    new_muts = []
    for mut in nt_mut_list.split(","):
        if mut:
            position = int(mut.split("|")[0].strip("ATCG-"))
            AAchange = ""
            for cds in SARS2_CDS:
                if position >= SARS2_CDS[cds][0] and position <= SARS2_CDS[cds][1]:
                    AApos = str(int(((position - SARS2_CDS[cds][0])//3)+1))
                    for AAmut in aa_subs.split(","):
                        if cds in AAmut:
                            if ":" in AAmut and AApos == AAmut.split(":")[1][1:-1]:
                                AAchange = AAmut
                    break
            if AAchange:
                new_muts.append(f"{mut}|{AAchange}")
            else:
                new_muts.append(mut)

    return ",".join(new_muts)


def del_list_aa_match(nt_del_list, aa_del_list):
    new_list = []
    for mut in nt_del_list.split(","):
        if mut:
            try:
                start_pos, end_pos = mut.split("-")
            except ValueError:
                start_pos, end_pos = mut, mut
            start_pos, end_pos = int(start_pos), int(end_pos)
            start_cds = ""
            end_cds = ""
            for cds in SARS2_CDS:
                if start_pos >= SARS2_CDS[cds][0] and start_pos <= SARS2_CDS[cds][1]:
                    start_cds = cds
                if end_pos >= SARS2_CDS[cds][0] and end_pos <= SARS2_CDS[cds][1]:
                    end_cds = cds
            if start_cds and end_cds and start_cds == end_cds:
                aa_start_pos = int(((start_pos - SARS2_CDS[end_cds][0])//3)+1)
                aa_end_pos = int(((end_pos - SARS2_CDS[end_cds][0])//3)+1)
                AAchange = []
                for aamut in aa_del_list.split(","):
                    if end_cds in aamut:
                        position = ""
                        for c in aamut.split(":")[-1]:
                            if c.isnumeric():
                                position += c
                        position = int(position)
                        if position >= aa_start_pos and position <= aa_end_pos:
                            AAchange.append(aamut)
                if AAchange:
                    new_list.append(f"{mut}|{';'.join(AAchange)}")
                else:
                    new_list.append(mut)
            elif start_cds or end_cds:
                new_list.append(f"{mut}|orf(s) disrupted:{start_cds}{end_cds}")

    return ",".join(new_list)

def aa_match(entry_dict):
    aa_matched_dict = {}
    for dayloc in entry_dict:
        if dayloc == "seqName":
            aa_matched_dict[dayloc] = entry_dict[dayloc]
        else:
            for entry in entry_dict[dayloc]:
                split_entry = entry.split(";;")
                split_entry[5] = mut_list_aa_match(split_entry[5], split_entry[-1])
                split_entry[6] = mut_list_aa_match(split_entry[6], split_entry[-1])
                split_entry[2] = del_list_aa_match(split_entry[2], split_entry[7])
                try:
                    aa_matched_dict[dayloc][";;".join(split_entry)] = entry_dict[dayloc][entry]
                except KeyError:
                    aa_matched_dict[dayloc] = {";;".join(split_entry) : entry_dict[dayloc][entry]}

    return aa_matched_dict


def get_date(filename):
    if filename.endswith("-new-mutations.txt"):
        return filename.split("-new-mutations")[0].replace("_","")
    elif filename.startswith("new_mutations_"):
        return filename.split("_")[-2].replace("_","")



def main():

    args = arg_parser()

    timestamp = time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime(time.time()))

    (collapsed_entries, total_entries) = parse_input(args.in_tsv)

    collapsed_entries = aa_match(collapsed_entries)

    unique_day_loc_seq = 0

    muts_of_interest = []

    new_muts_of_interest = []


    if args.mut_file:
        for file_list in args.mut_file:
            for file in file_list:
                muts_of_interest += get_muts(file)

    if args.n == 1:
        parse_path = os.getcwd()
        if os.path.isdir("/work/results/"):
            parse_path = "/work/results/"
        old_new_mut_files = []
        for file in os.listdir(parse_path):
            if file.endswith(".txt") and ( file.endswith("-new-mutations.txt") or file.startswith("new_mutations_")):
                old_new_mut_files.append(os.path.join(parse_path, file))
        if old_new_mut_files:
            old_new_mut_files = sorted(old_new_mut_files, key=lambda x: get_date(x))
            muts_of_interest += get_muts(old_new_mut_files[-1])
        

    if args.p:

        """
        parses the dict of collapsed nextclade samples to get information of the occurance of individual mutations
        """
        mut_loc_time_dict = {}
        samp_counts = {}
        for loc_time in collapsed_entries:
            try: # ignores samples that aren't from gisaid
                loc, date = loc_time.split("|")
            except ValueError:
                continue
            date = int("".join(date.split("-")))

            for muts in collapsed_entries[loc_time]:

                try: # get the count of samples at each location
                    samp_counts[loc] += collapsed_entries[loc_time][muts]
                except KeyError:
                    samp_counts[loc] = collapsed_entries[loc_time][muts]
                try: # get the count of samples on each date at each location
                    samp_counts["|".join((loc, str(date)))] += collapsed_entries[loc_time][muts]
                except KeyError:
                    samp_counts["|".join((loc, str(date)))] = collapsed_entries[loc_time][muts]
                split_mut_types = muts.split(";;")
                # only looking at indels, and private non-reversion subs
                for mut in ",".join((split_mut_types[2], split_mut_types[3], split_mut_types[5], split_mut_types[6],)).split(","):
                    if mut:
                        try:
                            mut_loc_time_dict[mut]
                        except KeyError:
                            mut_loc_time_dict[mut] = {}

                        try:
                            mut_loc_time_dict[mut][loc]
                        except KeyError:
                            mut_loc_time_dict[mut][loc] = {}

                        try:
                            mut_loc_time_dict[mut][loc][date] += collapsed_entries[loc_time][muts]
                        except KeyError:
                            mut_loc_time_dict[mut][loc][date] = collapsed_entries[loc_time][muts]

        for mut in mut_loc_time_dict: # test if mutations meet the requirements to be considered a new potential interesting mutation
            total_mut_count = 0
            for loc in mut_loc_time_dict[mut]:

                mut_loc_count = 0
                loc_cal = calendar.Calendar()
                consecutive_days = 0
                years = list(mut_loc_time_dict[mut][loc].keys())

                for year in range((min(years)//10000), (max(years)//10000)+1):
                    for month in range(1, 13):
                        for day in loc_cal.itermonthdays(year, month):
                            if not day == 0:
                                try:
                                    date = (year * 10000) + ( month * 100) + day
                                    mut_loc_count += mut_loc_time_dict[mut][loc][date]
                                    total_mut_count += mut_loc_time_dict[mut][loc][date]
                                    if (mut_loc_time_dict[mut][loc][date] / samp_counts["|".join((loc, str(date)))] ) >= args.freq_threshold:
                                        consecutive_days += 1

                                except KeyError:
                                    if consecutive_days > args.days and not mut in new_muts_of_interest:
                                        new_muts_of_interest.append(mut)
                                    consecutive_days = 0

                if consecutive_days > args.days and not mut in new_muts_of_interest:
                    new_muts_of_interest.append(mut)

    if new_muts_of_interest: # write new mutation of interest to file
        with open(f"new_mutations_{timestamp}.txt", "w") as new_mut_out:
            for mut in new_muts_of_interest:
                if not mut in muts_of_interest:
                    muts_of_interest.append(mut)
            for mut in muts_of_interest:
                new_mut_out.write(f"{mut}\n")


    if args.out_file:
        out_f = args.out_file
    else:
        out_f = f"{args.in_tsv.name}.col.tsv"
    mut_background_dict = {}
    with open(out_f, "w") as output:
        output.write("seqName\tcount\t")
        output.write('\t'.join(list(collapsed_entries['seqName'].keys())[0].split(';;')[:-1]))
        if muts_of_interest:
            output.write("\tmutations of interest")
        output.write("\n")
        for day_loc in sorted(collapsed_entries.keys()):
            if not day_loc == "seqName":
                unique_day_loc_seq += len(collapsed_entries[day_loc])
                for seq in collapsed_entries[day_loc]:
                    output.write(f"{day_loc}\t{collapsed_entries[day_loc][seq]}")
                    seq_entries = seq.split(';;')
                    for i in range(0, len(seq_entries)-1):

                        output.write(f"\t")
                        if i < 4 or i > 6:
                            output.write(f"{seq_entries[i]}")
                        elif seq_entries[i].startswith("private"):
                            output.write(f"{seq_entries[i]}")
                        else:
                            if seq_entries[i]:
                                for mut in seq_entries[i].split(","):
                                    output.write(mut)
                    if muts_of_interest:
                        output.write("\t")
                        for mut in muts_of_interest:
                            if mut in seq:
                                output.write(f"{mut}, ")
                                try:
                                    mut_background_dict[mut]
                                except KeyError:
                                    mut_background_dict[mut] = {seq_entries[0] : 1}
                                else:
                                    try:
                                        mut_background_dict[mut][seq_entries[0]] += 1
                                    except KeyError:
                                        mut_background_dict[mut][seq_entries[0]] = 1

                    output.write("\n")



    print(f"{total_entries} collapsed to {len(collapsed_entries)-1} day locations and {unique_day_loc_seq} sequence entries")
    if mut_background_dict:
        with open(f"Mutation_Clade_Backgrounds_{timestamp}.txt", "w") as mut_clades:
            for mut in mut_background_dict:
                mut_clades.write(f"{mut}:")
                for clade in mut_background_dict[mut]:
                    mut_clades.write(f" {clade} : {mut_background_dict[mut][clade]},")
                mut_clades.write("\n")

if __name__ == "__main__":
    main()


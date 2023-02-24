#!/bin/bash

while getopts "d:t:n:m:" opt; do
    case $opt in
	d)
	    DAYS="--parser_days $OPTARG"
	    ;;
	t)
	    FREQ="--parser_freq $OPTARG"
	    ;;
	n)
	    NSEQS="--num_seqs $OPTARG"
	    ;;
	m)
	    MUTS="--parser_mut_list $OPTARG"
	    ;;
    esac
done
shift $((OPTIND - 1))

FASTA=$1
GDIR="/work/data/"

if [[ -z $FASTA ]];
then
  echo "Usage: tracker.sh [options] fastafile"
  echo
  echo "Options:"
  echo "  -n N | Split input fasta file in blocks of N sequences each (default: 100000)"
  echo "  -d D | Days to track mutations of interest (default: 14)."
  echo "  -t T | Minimum frequency to track mutations of interest (default: 0.2)."
  echo "  -m M | List of mutations of interest."
  exit 1
fi

if [[ -d $GDIR ]];
then
  GDIR="--nc_data $GDIR"
else
  echo "Error: genome data directory not found."
  exit 1
fi

RUNID=$(date --rfc-3339=date)

if [[ -z $MUTS ]];
then
  if [[ ! -z "$(ls /work/results/*-new-mutations.txt)" ]];
  then
    previous_mutations=$(ls -t /work/results/*-new-mutations.txt | head -1)
    MUTS="--parser_mut_list $previous_mutations"
  fi
fi

nextflow run -resume /work/bin/tracker.nf --name $RUNID --fastafile `pwd`/$FASTA $DAYS $FREQ $NSEQS $MUTS $GDATA



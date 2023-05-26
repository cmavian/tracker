#!/usr/bin/env nextflow

/* Parameters */

params.name = ""
params.fastafile = ""
params.num_seqs = 100000
params.nc_data = "data/sars-cov-2"
params.parser_mut_list = ""
params.parser_days = 14
params.parser_freq = 0.2

/* Processes */

process SplitFasta {
	executor "local"

	output:
	path 'fragment-*.fa' into nextclade_ch

	script:
	"""
#!/bin/bash
#python3.9 /work/bin/NinetyDay.py ${params.fastafile} 
python3.9 /work/bin/splitFasta.py ${params.fastafile} fragment ${params.num_seqs} #last90.fasta
"""
}

Channel
	.from nextclade_ch
	.flatten()
	.set { nextclade_one }

process Nextclade {
	memory "5G"
	time "5:00:00"
	cpus 8

	input:
	path frag from nextclade_one

	output:
	path 'ncoutput/nextclade.tsv' into merge_nextclade_ch

	script:
	"""

nextclade dataset get --name 'sars-cov-2' --output-dir 'data/sars-cov-2'

nextclade run -j 8 -D ${params.nc_data} --output-all=ncoutput $frag
"""
}

Channel 
	.from merge_nextclade_ch
	.collectFile(keepHeader: true, skip: 1, name: 'nextclade.tsv')
	.set { merge_nextclade }

process ParseNextclade {
	executor "local"

	input:
	val filename from merge_nextclade

	output:
	path 'nextclade.tsv'
	path 'parsed.tsv'
	path 'Mutation_Clade_Backgrounds_*.txt', optional: true
	path 'new_mutations_*.txt', optional: true

	publishDir 'results', pattern: '*.tsv', mode: 'copy', saveAs: { fn -> "${params.name}-${fn}" }
	publishDir 'results', pattern: 'Mutation_Clade*.txt', mode: 'copy', saveAs: { fn -> "${params.name}-mutation-clades.txt" }
	publishDir 'results', pattern: 'new_mutations_*.txt', mode: 'copy', saveAs: { fn -> "${params.name}-new-mutations.txt" }

	script:
	"""
#!/bin/bash

cp $filename .
python3.9 /work/bin/NC_parse.py -p -m ${params.parser_mut_list} -d ${params.parser_days} -t ${params.parser_freq} -o parsed.tsv $filename
"""
}


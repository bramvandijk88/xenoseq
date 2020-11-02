# Xenoseq
Pipeline to compare multiple short-read libraries, 
and assemble contigs that only occur in one library.

## Overview
Xenoseq is pipeline of (meta)genomic tools to extract
contigs that are unique to one dataset. To use it, you
have to submit a short-read library as a query file (in
which it will look for unique contigs), and one or more 
subject files (which will be used as a reference for 
which reads are unique). Both these files need to be 
supplied in fasta format. If you have uncorrected fastq
files, you have to first quality trim/merge/convert them
yourself (using e.g. prinseq-lite), or use the xenoseq_prep
example script to use FLASh and prinseq for this step.

## Install: conda setup (required) 
To use xenoseq, a conda environment is provided. 
(see docs.conda.io/projects/conda/en/latest/) 
To install this conda environment, simply use: 

```bash
conda env create -f conda_env.yml
conda activate xenoseq
```

If you prefer using other virtual environments, this 
yml-file lists all the required packages and their versions.

## Install: add to path (optional)

If you want to use the pipeline scripts from anywhere on 
your machine, you have to add the xenoseq directory (where
you cloned the repository) to your global PATH variable: 

```bash
export PATH='$PATH:<XENOSEQ_DIR>'
```

e.g. if you cloned / downloaded xenoseq into your home dir:

```bash
export PATH='$PATH:~/xenoseq'
```

If you don't want to do this each time you login to a new
terminal, add the export PATH code to you ~/.bashrc 
(or ~/.profile) file:

```bash
echo export PATH='$PATH:~/<XENOSEQ_DIR>' >> ~./bashrc
```

## Running 
To test whether everything is working, try using the example
programL

```
> ./xenoseq_example
```

This example will use mock reads found in samples/reads and 
will search for (artficially added) unique contigs:

Your unique contigs will be in the file xenoseq_contigs.fasta
in the provided output directory (e.g. example_out)". This 
example also benchmarks the quality of this in silico dataset,
which despite the relatively low coverage still shows how 100%
of the unique contigs are retrieved. (however, not over their
entire length)

If you want to run only the main pipeline which takes raw reads
and assembles unique contigs, this is an example command:

```bash
xenoseq \
	-s subject_library_1.fasta \
	-s subject_library_2.fasta \
	-s subject_library_n.fasta \
	-q query_library.fasta \
	-o example_out/unique_contigs.fasta 

	(-c 4) 		# Nr. of threads
	(-megahit) 	# Megahit instead of spades
	(-blastn) 	# Blastn instead of BWA mem
```

Or see the help page:

```
xenoseq -h
```

## Provided scripts / binaries:
xenoseq_bin/flash  	 	| FLASh software for merging paired-end reads
xenoseq_bin/prinseq-lite 	| Quality trim, remove duplicates, remove adapters
xenoseq_bin/seek_uniq_seq.py 	| Python-code to use BLAST to detect reads (this is slow, and no longer the default)
xenoseq_prep 	 		| BASH-script for prepping paired-end reads with FLASh and prinseq
xenoseq		 		| BASH-script for main pipeline.
xenoseq_example  		| BASH-script for example of xenoseq (output in this dir/example_out) 

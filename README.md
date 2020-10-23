# Xenoseq
Pipeline to compare multiple short-read libraries, and assemble contigs that only occur in one library. 

## Overview
Xenoseq is pipeline of (meta)genomic tools to extract contigs that are unique to one dataset. To use it, you have to submit a short-read library as a query file (in which it will look for unique contigs), one or more subject files (which will be used as a reference for which reads are unique) in fasta format. If you have uncorrected fastq files, you either have to quality trim / merge them yourself, or use the xenoseq_prep command to use FLASh and PRINSEQ to merge and quality check reads. 

## Install 
To use xenoseq, a conda environment is provided. (see https://docs.conda.io/projects/conda/en/latest/) To install your conda environment, simply use: 

> conda env create -f conda_env.yml
> conda activate xenoseq

If you prefer using other virtual environments, this yml-file will provide you with all the required packages and their version.

## Running 
To test whether everything is working, try using the example program './xenoseq_example', which will use mock reads found in samples/reads and will search for
(artficially added) unique contigs:

> ./xenoseq_example

Your unique contigs will be in the file "xenoseq_contigs.fasta". This example also benchmarks the quality of this (low coverage) in silico dataset to show how 100% of the unique contigs are retrieved, albeit not over their entire length.

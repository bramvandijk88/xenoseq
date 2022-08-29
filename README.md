# Xenoseq

```
	 __    __                                                             
	/  |  /  |                                                            
	AT |  GC |  ______   _______    ______    _______   ______    ______  
	GG  \/AG/  /      \ /       \  /      \  /       | /      \  /      \ 
 	TA  TT<  /GATACT  |TAAATGG  |/CCGTAA  |/AATAAAS/ /ATTTCT  |/ATGTTA  |
	  TGAC  \ AA    GA |TG |  GG |AG |  AG |TT      \ TT    CG |GA |  AT |
	 AA /AT  |GATCCCGT/ TA |  AA |CG \__AG | GGAACT  |TACGGGTA/ GT \__AT |
	AA |  GC |CC       |GA |  GT |GC    GT/ /     GG/ TA       |AG    AT |
	GG/   TG/  GTAGGCC/ CA/   TT/  TAAATG/  ATGCGCG/   ATGCAAT/  AGGGTTT |
	                                                                  AA |
	                                                                  AA |
 	                                                                  AA/ 
```
Xenoseq is a simple bioinformatic pipeline to find sequences that appear to be newly
introduced into a community. The input are (sets of) query and reference samples,
which the pipeline will use to detect:

* unique_contigs.fasta; sequences in query not present in the reference
* xenotypic_contigs.fasta; the subset of unique contigs that can be linked to *another* reference
* xenotypic_coverage\*.txt; text files describing the coverage of the xenotypic sequences in other samples

## Overview
Xenoseq wraps read trimming (fastp), assembly (megahit), read mapping (BWA),
read filtering (samtools), and local alignment (blast), to detect putative evidence
of horizontal transfer between communities. This tool was used in a recent publication
(<link>) to detect the movement of MGEs and nanobacteria in compost communities.

## Install: conda setup (required)
To use xenoseq, a conda environment file is provided to install the above mentioned
dependencies (for information on conda, see docs.conda.io/projects/conda/en/latest/)

To install the xenoseq environment, simply use:

```
conda env create -f environment.yml
conda activate xenoseq
```

## Install: add to path (optional)

To run the pipeline, you either have to provide the full path to
the xenoseq binary (e.g. /home/user/XENOSEQ_DIR/xenoseq, or
add the xenoseq directory (where you cloned the repository)
to your global PATH variable:

```
export PATH='$PATH:<XENOSEQ_DIR>'
```

e.g. if you cloned / downloaded xenoseq into your home dir:

```
export PATH='$PATH:~/XENOSEQ_DIR'
```

If you don't want to do this each time you login to a new
terminal, add the export PATH code to you ~/.bashrc
(or ~/.profile) file:

```
echo export PATH='$PATH:~/<XENOSEQ_DIR>' >> ~/.bashrc
```

## Metadata file

Usage of the pipeline requires a file listing all the "queries" and "subjects" in a
text file. The pipeline will look for unique reads in the query files by comparing
them to the corresponding subjects. These reads will be assembled into contigs. A
subset of these contigs will be "xenotypic", i.e. having a foreign origin, by
aligning the sequences to the remaining subjects, providing extra evidence for
horizontal transfer of viral sequences or genes.

The metadata should look like this:

```
# SAMPLES (Query, Reference)
Horizontal1	Ancestral1
Horizontal2	Ancestral2
Horizontal3	Ancestral3
Horizontal4	Ancestral4
Vertical1	Ancestral1
Vertical2	Ancestral2
Vertical3	Ancestral3
Vertical4	Ancestral4
```

To run the pipeline on the example data, use:
```
> ./xenoseq -m example_metadata.tsv -o Xenoseq_example -t
```

This example will use mock reads found in samples/reads and will search for xenotypic contigs in simulated data. Your unique/xenotypic contigs will then be
stored in Xenoseq_example/<QUERY>. By default, xenoseq will assume your reads will
be stored in /samples/reads, and will assume the files correspond to the names
in the metadata file with the suffix \"\_R1.fq\" and \"\_R2.fq\" (these options
can be changed with -p and -r). In other words; if you have your own reads you need to
modify the metadata (and optionally, modify the read path/prefix). 

The usage of the trace option (-t) will ensure all samples will be mapped back 
to the unique contigs and generate coverage statistics. 

## All options
```bash

Usage:
         xenoseq -m <meta_data_tsv> -o <output_dir> -c <num_cores> -l -t
Mandatory:
        -m/--metadata           File containing the metadata (tsv file with query-reference sets)
Optional options:
        -p/--path_to_reads <STRING>     Path to reads for samples in metadata
        -r/--read_suffix <STRING>       Read suffix for paired files in metadata (e.g. _R*.fq for using _R1.fq and _R2.fq)
        -l/--link                       After detecting unique contigs, attempt to link them to other reference samples.
        -t/--trace                      After detecting xenotypic contigs, trace them across all samples.
        -c/--cores <INT>                Number of threads to use in parallisable parts of the pipeline
        -o/--output <STRING>            Output directory to put all the data
        -L/--alignment_length           Minimal alignment length to link unique sequences to other reference samples.
        -P/--alignment_pid              Minimal percent identity to link unique sequences to other reference samples.
```

If you want to modify any of the options (e.g. filtering thresholds, quality trimming),
you can modify the relevant subroutines of xenoseq given in xenoseq_bin/functions.sh.
Generally, changing these options is not required.
```
xenoseq -h
```

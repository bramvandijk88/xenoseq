#!/bin/bash

function trimheaders() {
	cat $1 | awk -v max=${2} "{if(\$0~\">\"){printf(\"%s\n\",substr(\$1,1,max))}else{print \$0}}";
}
function makeidsuniq() {
	perl -pe 's/$/_$seen{$_}/ if ++$seen{$_}>1 and /^>/; ' $1 ;
}
function concatfasta() {
	grep -v "^>" $1 | awk 'BEGIN { ORS=""; print ">Concatenated_Sequence\n" } { print } END { print "\n"} '
}

function parse_metadata() {
	query_samples=()
	reference_samples=()
	output_dirs=()
	# ancestral=$(cat $1 | cut -f2)
	while read p; do
		if [ ! ${p::1} == "#" ]; then
			query=$(echo $p | cut -d ' ' -f1)
			query_samples+=($query)
			ref=$(echo $p | cut -d ' ' -f2)
			reference_samples+=($ref)
			out=$(echo $p | cut -d ' ' -f3)
			output_dirs+=($out)
		fi
	done < ${1}
}

function trim() {
	r1=$1
	r2=$2
	out=$3
	if [ -f $out ]; then
		command=""
		echo -e "[xenoseq_trim     $(date +%d-%m_%H:%M:%S)] ${ORANGE}NOTE: no job created, using pre-existing $out.${NC}"
	else
		command="${SCRIPT_PATH}/xenoseq_bin/xenoseq_prep -r1 $r1 -r2 $r2 -c $cores -o $out;"
		echo -e "[xenoseq_cmd      $(date +%d-%m_%H:%M:%S)] ${GREY}$command${NC}"
		#eval $command
		#success $? "Trimming (fastp)"
	fi
}

function samtools_coverage() {
	bam_name=$1
	if [ -f ${1}_coverage.txt ]; then
		command=""
		echo -e "[xenoseq_cover    $(date +%d-%m_%H:%M:%S)] ${ORANGE}NOTE: using pre-existing ${1}_coverage.txt${NC}"
	else
		command="samtools coverage ${1}.sorted.bam"
		echo -e "[xenoseq_cmd      $(date +%d-%m_%H:%M:%S)] ${GREY}$command${NC}"
		$command > ${1}_coverage.txt
		success $? "Samtools coverage"
	fi
}

function link_contig (){
	if [ -f $3 ]; then 
		echo -e "[xenoseq_link     $(date +%d-%m_%H:%M:%S)] ${ORANGE}NOTE: no job created, using pre-existing ${3}${NC}"
	else
		blastn -query $1 -db $2 -outfmt '6 qseqid sseqid pident length slen' \
			| awk -v blength=$blength -v bpid=$bpid '{if($4>blength && $3 > bpid)print $1}' | sort | uniq | sort -nr | awk -v r="$ref" '{print $0"\t"r}' >> $3
		success $? "Blast linking"
	fi
}

# assembly with megahit
function assemble_mh() {
	reads=$1
	minlen=$2
	out=$3

	if [ -f $out/final.contigs.fa ]; then
		command=""
		echo -e "[xenoseq_mega     $(date +%d-%m_%H:%M:%S)] ${ORANGE}NOTE: no job created, using pre-existing $out/final.contigs.fa${NC}"
	else
		rm -rf $out
		command="megahit -r $reads -t $cores --min-contig-len $minlen --presets meta-large -o $out;"
		echo -e "[xenoseq_cmd      $(date +%d-%m_%H:%M:%S)] ${GREY}$command${NC}"
		# $command 2> ${output}/${reference_samples[$i]}/logs/megahit.log
		# success $? "Assembly (megahit)"
	fi
}

# Make burrows-wheeler index
function bwa_index()  {
	if [ -f $1.bwt ]; then
		command=""
		echo -e "[xenoseq_indx     $(date +%d-%m_%H:%M:%S)] ${ORANGE}NOTE: no job created, using pre-existing $1.bwt${NC}"
	else
		command="bwa index $1;"
		echo -e "[xenoseq_cmd      $(date +%d-%m_%H:%M:%S)] ${GREY}$command${NC}"
		#$command 2> ${output}/${reference_samples[$i]}/logs/bwa_index.log
		#success $? "Bwa index"
	fi
}

# Make blastdb
function blastdb() {
	if [ -f $1.nhr ]; then
		echo -e "[xenoseq_blasdb   $(date +%d-%m_%H:%M:%S)] ${ORANGE}NOTE: no job created, using pre-existing blastdb $1.nhr${NC}"
	else
		command="makeblastdb -in $1 -dbtype nucl;"
		echo -e "[xenoseq_cmd      $(date +%d-%m_%H:%M:%S)] ${GREY}$command${NC}"
		# $command > ${output}/${reference_samples[$i]}/logs/blastdb.log 2>&1
		# success $? "Make blastdb"
	fi
}

# Map reads using BWA MEM
function bwa_map()  {
	if [ -f $3.sorted.bam ]; then
		echo -e "[xenoseq_map      $(date +%d-%m_%H:%M:%S)] ${ORANGE}NOTE: using pre-existing file $3.sorted.bam${NC}"
	else
		command="bwa mem -t $cores $1 $2"
		echo -e "[xenoseq_cmd      $(date +%d-%m_%H:%M:%S)] ${GREY}$command${NC}"
		$command > $3.bam 2> ${output}/${reference_samples[$i]}/logs/bwa_map.log
		success $? "Bwa mapping"
		command="samtools sort -@ $cores $3.bam -o $3.sorted.bam"
		echo -e "[xenoseq_cmd      $(date +%d-%m_%H:%M:%S)] ${GREY}$command${NC}"
		$command 2> ${output}/${reference_samples[$i]}/logs/samtools.log
		success $? "Samtools sort"
	fi
}

function unmapped_to_fasta() {
	if [ -f $output/${query_samples[$i]}/reads/unique_reads.fasta ]; then
		echo -e "[xenoseq_filter   $(date +%d-%m_%H:%M:%S)] ${ORANGE}NOTE: using pre-existing file $output/${query_samples[$i]}/reads/unique_reads.fasta${NC}"
	else
		echo -e "[xenoseq_filter   $(date +%d-%m_%H:%M:%S)] Samtools: extracting non-hits "
		samtools index $1
		command="samtools view -f 4 $1"
		$command > $output/${query_samples[$i]}/read_mapping/Unmapped_to_${2}.bam 2> ${output}/logs/samtools_filter.log
		success $? "samtools filter"
		if [ -s $output/${query_samples[$i]}/read_mapping/Unmapped_to_${2}.bam ]; then
			echo -e "[xenoseq_filter   $(date +%d-%m_%H:%M:%S)] Samtools: exporting fasta file for non-hits... "
			command="samtools fasta $output/${query_samples[$i]}/read_mapping/Unmapped_to_${2}.bam"
			$command > $output/${query_samples[$i]}/reads/unique_reads.fasta 2> ${output}/logs/samtools_fasta.log
			success $? "samtools fasta"
		fi
	fi
}

# Function to check if previous command was succesful
function success() {
	if [ ! $1 -eq 0 ]; then
  		echo -e "${RED}Error in xenoseq ($2). See relevant log files in ${output}/${query_samples[$i]}/logs (exit 1).${NC}" >&2
  		exit 1;
	fi
}

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



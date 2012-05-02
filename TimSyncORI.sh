#!/bin/bash

# expects 3 args: fastq of read1, fastq of read 2, and prefix for output files


# Only uses the first part (up to the first space) of a sequence name and assumes this is identical for read 1 and read 2 from the same fragment
# Will not carry through to the output the second part of an Illumina header >> decided to do it this way as it is unclear whether this will always be present + it affects the names of reads if it is not parsed out of the header.
# Expects non-compressed files >> compressed files could be made an option
# outputs 4 files
declare fastq1=$1
declare fastq2=$2
declare prefix=$3

# call with file name
function countFastqRecords(){
declare lines=`wc -l $1 | grep -o "[ 0-9 ]*" | tr -d " "`
echo $((${lines}/4))
}

function countFastqRecordsToScreen(){
    echo $1 `countFastqRecords ${1}`
}

echo
date
echo
echo "Input files number of fq records"
countFastqRecordsToScreen $fastq1 
countFastqRecordsToScreen $fastq2 

# Inputs are file name and read number
# Create a tab separated fastq format with the second column as the read: 1 or 2
# FS is space to only get first part of name
function tabularFq(){
cat $1 | awk -v read=$2 'BEGIN{FS=" "; OFS="\t";};{printf $1"\t"; if((NR+3)%4==0){printf read"\t"}; if((NR%4)==0){print ""}};END{}'}

# declare temporary files
declare tempFastq1=${fastq1}_temp_R1.txt
declare tempFastq2=${fastq2}_temp_R2.txt

# tabularise and write to temporary files for read 1 and 2
tabularFq $fastq1 1 > $tempFastq1
tabularFq $fastq2 2 > $tempFastq2

# Set up the output file names here 
declare pairedPrefix=${prefix}_paired_
declare singlePrefix=${prefix}_single_

# classify into 4 files: paired read1, paired read 2, single read 1 and single read2
cat ${tempFastq1} ${tempFastq2} | sort -k1,2 | awk -v pairedPrefix=${pairedPrefix} -v singlePrefix=${singlePrefix} 'BEGIN{OFS="\t"; FS="\t"; prevName=0; prevLine=0; prevRead=0; currentName=0; currentLine=0; currentRead=0};{currentName=$1; currentLine=$0; currentRead=$2; if(prevName!=0){if(prevName==currentName){file=pairedPrefix"r"prevRead".txt"; print prevLine > file; file=pairedPrefix"r"currentRead".txt";print currentLine > file; currentName=0}else{file=singlePrefix"r"prevRead".txt"; print prevLine > file}}; prevName=currentName; prevLine=currentLine; prevRead=currentRead};END{}'

# transform back from tabular to fastq
# need to drop the read info in the second column
function tabToFastq(){
cat $1 | awk 'BEGIN{OFS="\t"; FS="\t";};{print $1; print $3; print $4; print $5};END{}'

}

tabToFastq ${pairedPrefix}r1.txt >  ${pairedPrefix}r1.fq
tabToFastq ${pairedPrefix}r2.txt >  ${pairedPrefix}r2.fq
tabToFastq ${singlePrefix}r1.txt >  ${singlePrefix}r1.fq
tabToFastq ${singlePrefix}r2.txt >  ${singlePrefix}r2.fq

echo
echo "Output files number of fq records"
countFastqRecordsToScreen ${pairedPrefix}r1.fq
countFastqRecordsToScreen ${pairedPrefix}r2.fq
countFastqRecordsToScreen ${singlePrefix}r1.fq
countFastqRecordsToScreen ${singlePrefix}r2.fq
echo
echo "Paired output files should have same number of reads and sum of single and paired in output files should match corresponding input file"
echo
date

# cleanup
# Remove the tabular files of inputs
rm -f $tempFastq1 $tempFastq2

# Remove the tabular files of the classified reads
rm -f ${pairedPrefix}r1.txt ${pairedPrefix}r2.txt ${singlePrefix}r1.txt ${singlePrefix}r2.txt
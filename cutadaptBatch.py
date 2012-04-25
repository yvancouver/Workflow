# April 2012 Yvan Strahm
## One need to pass the R1 and R2 adapters. R1 is  usually the Truseq Adapter index
## R1 adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC first A added (ligation) and stop just before index
## R2 adapter is the TruSeq universal adapter in reverse complement (sure about that?)
## R2 adapter AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
## the script should traverse a directory found the read and call adapt with the correct option
## cutadapt -m 20 -a adapter -o fastq.clipped.gz original.fastq.gz
#
import re
import sys
import os

for arg in sys.argv:
    print arg , " is a ", type(arg)

r1_adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
r2_adapter = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"

root = sys.argv[1]
m = sys.argv[2]


def cutadaptMe():
    file = open("/Users/yvans/Home/Dropbox/travail/BRCA12/WORKFLOWS/WorkFlowGeneral.sh")
    while 1:
        line = file.readline()
        #print dirname
        if not line:
            break
        if (re.search("export READS1", line)):
            d.write(line.strip()+R1+"\n")
        elif (re.search("export READS2", line)):
            d.write(line.strip()+R2+"\n")
        elif (re.search("export DIR", line)):
            d.write(line.strip()+folder_write_workflow+"\n")
        elif (re.search("export WORKINGDIR", line)):
            d.write(line.strip()+folder_write_workflow+"/Analysis/\n")
        elif (re.search("export RG", line)):
            d.write(line.strip()+"\"@RG\\tID:subset\\tPL:ILLUMINA\\tSM:subset\"\n")
        else:
            d.write(line)

i = 1  

for dirname, dirnames, filenames in os.walk(root, topdown=True):
    for dir in dirnames:
        #print "DIRNAME\t",dirname
        #print "DIR\t",dir
        #print "FILENAMES\t",filenames
        for file in filenames:
            if (re.search("fastq.gz", file)):
                #print dirname,dir,file
                fileToBeClipped=str(os.path.join(dirname,dir,file))
                print fileToBeClipped
                cutadaptMe(fileToBeClipped,adapter,m)
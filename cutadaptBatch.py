# April 2012 Yvan Strahm
## One need to pass the R1 and R2 adapters. R1 is  usually the Truseq Adapter index
## R1 adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC first A added (ligation) and stop just before index
## R2 adapter is the TruSeq universal adapter in reverse complement (sure about that?)
## R2 adapter AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
## the script should traverse a directory found the read and call adapt with the correct option
## cutadapt -m 20 -a adapter -o fastq.clipped.gz original.fastq.gz
## args should be ['cutadapt', '-m', '20', '-a', 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC', '-o', '/Users/yvans/Home/Analysis/Cardio_analysis_01.17.2012/Project_Diag-Cardiomyopathy-2011-12-09/Sample_KA-001/KA-001_ATCACG_L001_R1_001.clipped_m20.fastq.gz', '/Users/yvans/Home/Analysis/Cardio_analysis_01.17.2012/Project_Diag-Cardiomyopathy-2011-12-09/Sample_KA-001/KA-001_ATCACG_L001_R1_001.fastq.gz']

import re
import sys
import os
#import subprocess
from subprocess import Popen, PIPE, STDOUT
import shlex

for arg in sys.argv:
    print arg , " is a ", type(arg)

r1_adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
r2_adapter = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
root = sys.argv[1]
m = sys.argv[2]


def cutadaptMe(file,adapter,m):
    #from subprocess import Popen, PIPE, STDOUT
    print file
    print adapter
    print m
    cmd = "cutadapt -m "+ m + " -a "+ adapter +" -o "+ file[:-8]+"clipped_m"+m+".fastq.gz" + "  "+file
    print "CMD\t", cmd
    args = shlex.split(cmd)
    print "ARGS\t",args
    print "\n\n"
    job = Popen(args, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    output = job.communicate()
    results=file[:-8]+".report"
    result_handle=open(results,"w+")
    result_handle.write(output)

for dirname, dirnames, filenames in os.walk(root, topdown=True):
    for dir in dirnames:
        #print "DIRNAME\t",dirname
        #print "DIR\t",dir
        #print "FILENAMES\t",filenames
        for file in filenames:
            if (re.search("R1_001.fastq.gz", file)):
                #print dirname,dir,file
                adapter = r1_adapter
                fileToBeClipped=str(os.path.join(dirname,file))
                print "FILETOBECLIPPED1 ",fileToBeClipped
                cutadaptMe(fileToBeClipped,adapter,m)
            elif (re.search("R2_001.fastq.gz", file)):
                #print dirname,dir,file
                adapter  = r2_adapter
                fileToBeClipped=str(os.path.join(dirname,file))
                print "FILETOBECLIPPED2 ",fileToBeClipped
                cutadaptMe(fileToBeClipped,adapter,m)
#!/usr/bin/python
# April 2012 Yvan Strahm
## One need to pass the R1 and R2 adapters. R1 is  usually the Truseq Adapter index
## R1 adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC first A added (ligation) and stop just before index
## R2 adapter is the TruSeq universal adapter in reverse complement (sure about that?)
## R2 adapter AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
## the script should traverse a directory found the read and call adapt with the correct option
## cutadapt -m 32 -a adapter -o fastq.clipped.gz original.fastq.gz
## -m 32 will remove every ry smaller than 32 nt, bwa seed is 32 by default
## one should give an gunziped file and retrieve a uncompress one, because TimSync is (for now) taking ony uncompress file
## args should be ['cutadapt', '-m', '32', '-a', 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC', '-o', '/Users/yvans/Home/Analysis/Cardio_analysis_01.17.2012/Project_Diag-Cardiomyopathy-2011-12-09/Sample_KA-001/KA-001_ATCACG_L001_R1_001.clipped_m20.fastq.gz', '/Users/yvans/Home/Analysis/Cardio_analysis_01.17.2012/Project_Diag-Cardiomyopathy-2011-12-09/Sample_KA-001/KA-001_ATCACG_L001_R1_001.fastq.gz']
## have a look here http://argparse.googlecode.com/svn/trunk/doc/argparse-vs-optparse.html
##                  http://docs.python.org/library/argparse.html#module-argparse
#

# Import system package

import re
import sys
import os


# Import parsing package
import argparse
import textwrap

# Import process control package 
from subprocess import Popen, PIPE, STDOUT
import shlex


parser = argparse.ArgumentParser(prog='cutadaptBatch',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
                                                            cutadaptBatch is use to clip a gunziped fastq file with, for now, some Illumina adapters
                                                            This was made in order to clean some Halo produced data.
                                                            Author: Yvan Strahm (yvan.strahm@gmail.com)
                                                            License: GPL 3.0 (http://www.gnu.org/licenses/gpl-3.0.txt)
                                                            '''
                                                            ),
                                 epilog="If any questions contact me good, luck!"
                                 )
parser.add_argument('-m', help='match region length, I use generally 32 [Default 32]', required=False, nargs='?', default=32, type=int, action="store")
parser.add_argument('-d', help='fastq.gz containing dir [Required]', required=True, nargs='?')
parser.add_argument('-r1', help='Read1 adapter [Default = AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC', default="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC")
parser.add_argument('-r2', help='Read2 adapter [Default = AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT',required=False, default="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT")
values = parser.parse_args()

def testArgs(values):
    if values.m == None:
        values.m = 32
    if values.d == None:
        sys.exit("please enter a directory containing the reads")
    print dir(values)
    print "m type Option\t", type(values.m)
    print "m  Option\t", values.m
    print "dir Option\t", values.d
    print values.r1
testArgs(values)

def cutadaptMe(f,adapter,m):
    cmd = "cutadapt -m "+ m + " -a "+ adapter +" -o "+ f[:-8]+"clipped_m"+m+".fastq.gz" + "  "+f
    args = shlex.split(cmd)
    job = Popen(args, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    output = job.communicate()
    results=f[:-8]+".report"
    result_handle=open(results,"w+")
    result_handle.write(output[0])

for dirname, dirnames, filenames in os.walk(values.d, topdown=True):
    for f in filenames:
        if (re.search("R1_001.fastq.gz", f)):
            adapter = r1_adapter
            fileToBeClipped=str(os.path.join(dirname,f))
            print "FILETOBECLIPPED1 ",fileToBeClipped
            cutadaptMe(fileToBeClipped,adapter,values.m)
        elif (re.search("R2_001.fastq.gz", f)):
            adapter  = r2_adapter
            fileToBeClipped=str(os.path.join(dirname,f))
            print "FILETOBECLIPPED2 ",fileToBeClipped
            cutadaptMe(fileToBeClipped,adapter,values.m)
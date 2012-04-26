# April 2012 Yvan Strahm
## One need to pass the R1 and R2 adapters. R1 is  usually the Truseq Adapter index
## R1 adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC first A added (ligation) and stop just before index
## R2 adapter is the TruSeq universal adapter in reverse complement (sure about that?)
## R2 adapter AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
## the script should traverse a directory found the read and call adapt with the correct option
## cutadapt -m 20 -a adapter -o fastq.clipped.gz original.fastq.gz
## args should be ['cutadapt', '-m', '20', '-a', 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC', '-o', '/Users/yvans/Home/Analysis/Cardio_analysis_01.17.2012/Project_Diag-Cardiomyopathy-2011-12-09/Sample_KA-001/KA-001_ATCACG_L001_R1_001.clipped_m20.fastq.gz', '/Users/yvans/Home/Analysis/Cardio_analysis_01.17.2012/Project_Diag-Cardiomyopathy-2011-12-09/Sample_KA-001/KA-001_ATCACG_L001_R1_001.fastq.gz']
## have a look here http://argparse.googlecode.com/svn/trunk/doc/argparse-vs-optparse.html
##                  http://docs.python.org/library/argparse.html#module-argparse
#

import re
import sys
import os
import argparse
import textwrap

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
parser.add_argument('-m', help='match region length, I use genreally 20', required=True, nargs='?', default=20, type=int)
parser.add_argument('-d', help='fastq.gz containing dir', required=True, nargs='?')
parser.add_argument('-r1', help='Read1 adapter',required=False)
parser.add_argument('-r2', help='Read2 adapter',required=False)
values = parser.parse_args()
print values.m
print values.d

def testArgs(values):
    pass
testArgs(values)

sys.exit("on s'arretes ici pour l'instant")

r1_adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
r2_adapter = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
root = sys.argv[1]
m = sys.argv[2]


def cutadaptMe(file,adapter,m):
    cmd = "cutadapt -m "+ m + " -a "+ adapter +" -o "+ file[:-8]+"clipped_m"+m+".fastq.gz" + "  "+file
    args = shlex.split(cmd)
    job = Popen(args, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    output = job.communicate()
    results=file[:-8]+".report"
    result_handle=open(results,"w+")
    result_handle.write(output[0])

for dirname, dirnames, filenames in os.walk(root, topdown=True):
    for file in filenames:
        if (re.search("R1_001.fastq.gz", file)):
            adapter = r1_adapter
            fileToBeClipped=str(os.path.join(dirname,file))
            print "FILETOBECLIPPED1 ",fileToBeClipped
            cutadaptMe(fileToBeClipped,adapter,m)
        elif (re.search("R2_001.fastq.gz", file)):
            adapter  = r2_adapter
            fileToBeClipped=str(os.path.join(dirname,file))
            print "FILETOBECLIPPED2 ",fileToBeClipped
            cutadaptMe(fileToBeClipped,adapter,m)
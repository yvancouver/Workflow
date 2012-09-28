#!/usr/bin/python
# Aout 2012 Yvan Strahm yvan.strahm@gmail.com

import sys
import pybedtools

import argparse
import textwrap

import CollectCov

parser = argparse.ArgumentParser(prog='bedtools_coverage',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
                                                            bedtools_coverage is use to return a bed file containing region with a desired minimal coverage
                                                            This was made in order to get a better idea of the coverage from Halo produced data.
                                                            Author: Yvan Strahm (yvan.strahm@gmail.com)
                                                            License: GPL 3.0 (http://www.gnu.org/licenses/gpl-3.0.txt)
                                                            '''
                                                            ),
                                 epilog="If any questions contact me! Good luck!"
                                 )
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
parser.add_argument('-t', help = 'Runs the testdoc module', required=False, action="store_true")
parser.add_argument('-a', help = 'bam file', required=False)
parser.add_argument('-b', help = 'bed file', required=False)
parser.add_argument('-c', help = 'coverage file, if it is provided no coverageBed will be called', required=False)
parser.add_argument('-n', help = 'minimal coverage required', required = True, type = int)

values = parser.parse_args()

def test_values(values):
# test if the a b pair is passed properly 
    if (values.a and values.b):
        try:
            with open(values.a) as f1:
                try:
                    with open(values.b) as f2:pass
                except IOError as e2:
                    print 'Please check that the bed file exists'
        except IOError as e1:
            print 'Please check that the bam file exists'
    elif values.c:
        pass
        #print values.c
        #print type(values.c)
    else:
        sys.exit("Please, check the path of the bam and/or the bed files")
        
def runCoverageBed(a_file,b_file):

    a = pybedtools.BedTool(a_file)
    if a._isbam:
        print "Is a ",a.fn," a bam file ? ", a._isbam
    else:
        message = "Check your bam file"+str(a.fn)
# for now use to bed file to developp the Collectcov script 
# but a should be a bam fiel produced by the pipeline
        #sys.exit(message)
    b = pybedtools.BedTool(b_file)
    #print "Is b ",b.fn," b bam file ? ", b._isbam
    #print"command coverageFile = a.coverage(b,d=True)"
    coverageFile = a.coverage(b,d=True)
    return coverageFile

test_values(values)

if (values.a and values.b):
    coverageFile = runCoverageBed(values.a,values.b)
    CollectCov.CollectCov(coverageFile,values.n)
else:
    CollectCov.CollectCov(values.c,values.n)

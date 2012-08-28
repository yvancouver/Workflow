#!/usr/bin/python
# Aout 2012 Yvan Strahm yvan.strahm@gmail.com

import sys
import pybedtools

import argparse
import textwrap

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
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument('-t', help='Runs the testdoc module', required=False, action="store_true")



#a = pybedtools.BedTool('/Volumes/ToveExFat/120615/040_KA005/020_refineAlignment/030_BQRecalGATK/all.realigned.markDup.baseQreCali.bam')
#a = pybedtools.BedTool(sys.argv[1])
#a = pybedtools.BedTool('/Volumes/ToveExFat/120615/040_KA005/020_refineAlignment/030_BQRecalGATK/all.realigned.markDup.baseQreCali1000b.bam')
#print type(a)
#print "Is a ",a.fn," a bam file ? ", a._isbam


#b = pybedtools.BedTool(sys.argv[2])
#b = pybedtools.BedTool('/Users/yvans/Home/Dropbox/travail/BED_GFF_INTERVALS/Valided_And_Correct/CARDIO/Galaxy_GFFtoBED_EGS179.r150.readableregion_b37_sorted_GeneName.bed')
#print type(b)
#print "Is b ",b.fn," a bam file ? ",b._isbam
#print"command c = a.coverage(b,d=True)"

#c = a.coverage(b,d=True)
c = pybedtools.BedTool('/Users/yvans/Home/Dropbox/travail/meetings/data/zero_coverage_agilent_halo_bait.bed')

'''''
1    1    50    NR_047544_exon_0_0_chr1_156052369_f    0    +    1    0
1    1    50    NR_047544_exon_0_0_chr1_156052369_f    0    +    2    0
1    1    50    NR_047544_exon_0_0_chr1_156052369_f    0    +    29    0
1    1    50    NR_047544_exon_0_0_chr1_156052369_f    0    +    30    0
7    1000    2000    NM_001077653_exon_7_0_chr7_35293105_r    0    -    446    0
7    1000    2000    NM_001077653_exon_7_0_chr7_35293105_r    0    -    447    0
7    1000    2000    NM_001077653_exon_7_0_chr7_35293105_r    0    -    459    0
7    1000    2000    NM_001077653_exon_7_0_chr7_35293105_r    0    -    460    0
7    1000    2000    NM_001166220_exon_5_0_chr7_35293105_r    0    -    446    0
7    1000    2000    NM_001166220_exon_5_0_chr7_35293105_r    0    -    447    0
7    1000    2000    NM_001166220_exon_5_0_chr7_35293105_r    0    -    459    0
7    1000    2000    NM_001166220_exon_5_0_chr7_35293105_r    0    -    460    0
'''
def CollectCov(coverageResult,cov):
    chrom = ""
    start = 0
    stop = 0
    feature = "None"
    colors = "\t255,0,0"
     
## Get the chrom, start and stop
## From the example above we can see that chrom, start, stop, feature,score, and strand stay the same only posinfeat change
## One should return an line with summarizing the "strech" of the feature satisfying the wanted coverage.
## In the first example it should return
## 1    1    30    NR_047544_exon_0_0_chr1_156052369_f    0    +
## 7    1446 1460    NM_001077653_exon_7_0_chr7_35293105_r    0    -
## 7    1446 1460    NM_001166220_exon_5_0_chr7_35293105_r    0    -
    for entry in  coverageResult[0:len(coverageResult)]:
        if int(entry[7]) <= int(cov):

            if feature != entry.name:
## new feature so have to get his name and the start position (column 6)
                if feature == "None":
                    if entry.start == 1:
                        start = entry.start+int(entry[6])-2
                    else:    
                        start = entry.start+int(entry[6])-1
                    feature = entry.name
                if feature != "None" and feature != entry.name:
                    print chrom+"\t"+str(start)+"\t"+str(stop)+"\t"+feature+"\t"+entry.score+"\t"+entry.strand+"\t"+str(start)+"\t"+str(stop)+colors
                    if entry.start == 1:
                        start = entry.start+int(entry[6])-2
                    else:    
                        start = entry.start+int(entry[6])-1
                    feature = entry.name
            else:
                chrom = entry[0]
                if entry.start == 1:
                    stop = entry.start + int(entry[6])-1
                else:
                    stop = entry.start + int(entry[6])
                
    last_entry = coverageResult[len(coverageResult)-1]
    print last_entry.chrom+"\t"+str(start)+"\t"+str(stop)+"\t"+last_entry.name+"\t"+last_entry.score+"\t"+last_entry.strand+"\t"+str(start)+"\t"+str(stop)+colors

CollectCov(c,0)

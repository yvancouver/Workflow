import sys
import pybedtools

a = pybedtools.BedTool('/Volumes/ToveExFat/120615/040_KA005/020_refineAlignment/030_BQRecalGATK/all.realigned.markDup.baseQreCali.bam')
#a = pybedtools.BedTool('/Volumes/ToveExFat/120615/040_KA005/020_refineAlignment/030_BQRecalGATK/all.realigned.markDup.baseQreCali1000b.bam')
#print type(a)
#print "Is a ",a.fn," a bam file ? ", a._isbam

b = pybedtools.BedTool(sys.argv[1])
#print type(b)
#print "Is b ",b.fn," a bam file ? ",b._isbam
#print"command c = a.coverage(b,d=True)"

c = a.coverage(b,d=True)
#c = pybedtools.BedTool('/Volumes/ToveExFat/120615/040_KA005/CoverageYvan/040_KA005_uscsRefSeq33genes_b37_zeroTest.bed')

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
    
## Get the chrom, start and stop
## From the example above we can see that chrom, start, stop, feature,score, and strand stau the same only posinfeat change
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
                    print chrom+"\t"+str(start)+"\t"+str(stop)+"\t"+feature+"\t"+entry.score+"\t"+entry.strand
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
    print last_entry.chrom+"\t"+str(start)+"\t"+str(stop)+"\t"+last_entry.name+"\t"+last_entry.score+"\t"+last_entry.strand
CollectCov(c,0)



'''print entry.chrom
                print entry.start
                print entry.stop
                print entry.name
                print entry.fields
                print entry.score
                print entry.strand
                print entry.fields[6]
                print entry.fields[7]
                
    chrom = entry[0]
    start = entry[1]
    stop = entry[2]
    feature = entry[3]
    score = entry[4]
    strand = entry[5]
    posinfeat = entry[6]
    coverage = entry[7]
    print chrom,"\t",start,"\t",stop,"\t",feature,"\t",score,"\t",strand,"\t",posinfeat,"\t",coverage
'''     
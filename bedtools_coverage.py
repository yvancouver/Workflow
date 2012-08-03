import sys
import pybedtools

#a = pybedtools.BedTool('/Volumes/ToveExFat/120615/040_KA005/020_refineAlignment/030_BQRecalGATK/all.realigned.markDup.baseQreCali.bam')
a = pybedtools.BedTool('/Volumes/ToveExFat/120615/040_KA005/020_refineAlignment/030_BQRecalGATK/all.realigned.markDup.baseQreCali1000b.bam')
#print type(a)
#print "Is a ",a.fn," a bam file ? ", a._isbam

b = pybedtools.BedTool(sys.argv[1])
#print type(b)
#print "Is b ",b.fn," a bam file ? ",b._isbam
#print"command c = a.coverage(b,d=True)"

#c = a.coverage(b,d=True)
c = pybedtools.BedTool('/Volumes/ToveExFat/120615/040_KA005/CoverageYvan/040_KA005_uscsRefSeq33genes_b37_zeroTest.bed')
'''for i in c[0:len(c)]:
    if int(i[7]) == 0:
        print str(i).strip()
'''
#chrom = ""
#start = ""
#stop = ""
feature = "None"
#score = ""
#strand = ""
posinfeat = ""
coverage = ""

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
def CollectCov(coverageResult,cov,feature):
    i = 0
    chrom = ""
    start = 0
    stop = 0
    score = ""
    strand = ""
    lengthFeature = 0
    #print "chrom\tstart\tstop\tfeature\t\t\t\t\tscore\tstrand"
    for entry in  coverageResult[0:len(coverageResult)]:
        i += 1
        #print "at line ",str(i),"entry is ",entry
        if int(entry[7]) <= int(cov):
## Get the chrom, start and stop
## From the example above we can see that chrom, start, stop, feature,score, and strand stau the same only posinfeat change
## One should return an line with summarizing the "strech" of the feature satisfying the wanted coverage.
## In the first example it should return
## 1    1    30    NR_047544_exon_0_0_chr1_156052369_f    0    +
## 7    1446 1460    NM_001077653_exon_7_0_chr7_35293105_r    0    -
## 7    1446 1460    NM_001166220_exon_5_0_chr7_35293105_r    0    -

            if feature != entry.name:
                print "entry_start",entry.start
                print "entry_6", entry[6]
            
                start = entry.start+int(entry[6])-1
                current_feature = entry.name
                feature = entry.name
                if feature != "None":
                    print chrom,"\t",str(start-1),"\t",str(stop),"\t",current_feature,"\t",entry.score,"\t",entry.strand
                #print "current :",current_feature
                #print chrom,"\t",str(start-1),"\t",str(stop),"\t",current_feature,"\t",entry.score,"\t",entry.strand
            else:
                #print chrom,"\t",str(start-1),"\t",str(stop),"\t",current_feature,"\t",entry.score,"\t",entry.strand
                lengthFeature = int(entry[6])
                chrom = entry[0]
                stop = entry.start+int(entry[6])
                score = entry.score
                strand = entry.strand
                current_feature = entry.name
    
    last_entry = coverageResult[len(coverageResult)-1]
    print last_entry.chrom,"\t",str(last_entry.start+lengthFeature),"\t",str(stop),"\t",last_entry.name,"\t",last_entry.score,"\t",last_entry.strand
CollectCov(c,0,feature)



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
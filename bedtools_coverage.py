import sys
import pybedtools

a = pybedtools.BedTool('/Volumes/ToveExFat/120615/040_KA005/020_refineAlignment/030_BQRecalGATK/all.realigned.markDup.baseQreCali.bam')
print type(a)
print "Is a ",a.fn," a bam file ? ", a._isbam

b = pybedtools.BedTool(sys.argv[1])
print type(b)
print "Is b ",b.fn," a bam file ? ",b._isbam
print"command c = a.coverage(b,d=True)"

c = a.coverage(b,d=True)

for i in c[0:len(c)]:
    print i

chrom = ""
start = ""
stop = ""
feature = ""
score = ""
stand = ""
posinfeat = ""
coverage = ""

###################################################################################################################
# chrom    start        stop        feature                                score    stand    posinfeat  coverage  #
# 1        156052368    156052454   NR_047544_exon_0_0_chr1_156052369_f    0        +        1          0         #
#                                    .                                                                            #
#                                    .                                                                            #
#                                    .                                                                            #
# 1        156052368    156052454   NR_047544_exon_0_0_chr1_156052369_f   0         +        86         0         #
# 1        156052738    156052975   NR_047544_exon_1_0_chr1_156052739_f   0         +        129        0         #
#                                    .                                                                            #
#                                    .                                                                            #
#                                    .                                                                            #
# 1        156052738    156052975   NR_047544_exon_1_0_chr1_156052739_f    0        +        216        0         #
# 1        156084460    156085065   NM_005572_exon_0_0_chr1_156084461_f    0        +        29         0         #
#                                    .                                                                            #
#                                    .                                                                            #
#                                    .                                                                            #
# 1        156084460    156085065   NM_170708_exon_0_0_chr1_156084461_f    0        +        54         0         #
# 1        156084503    156085065   NR_047544_exon_3_0_chr1_156084504_f    0        +        1          0         #
# 1        156084503    156085065   NR_047544_exon_3_0_chr1_156084504_f    0        +        2          0         #
#                                    .                                                                            #
#                                    .                                                                            #
#                                    .                                                                            #
# 7        35293104    35293711     NM_001077653_exon_7_0_chr7_35293105_r    0      -        446        0         #
# 7        35293104    35293711     NM_001077653_exon_7_0_chr7_35293105_r    0      -        447        0         #
# 7        35293104    35293711     NM_001077653_exon_7_0_chr7_35293105_r    0      -        448        0         #
#                                    .                                                                            #
#                                    .                                                                            #
#                                    .                                                                            #
# 7        35293104    35293711     NM_001166220_exon_5_0_chr7_35293105_r    0      -        459        0         #
# 7        35293104    35293711     NM_001166220_exon_5_0_chr7_35293105_r    0      -        460        0         #
###################################################################################################################
def CollectCov(coverageResult,cov):
    for entry in  coverageResult[1:len(coverageResult)]:
        if int(entry[7]) <= int(cov):
## Get the chrom, start and stop
            if chrom != entry[0] and start != entry[1] and stop !=entry[3] and feature != entry[4]:
                chrom = entry[0]
                start = entry[1]
                stop = entry[2]
                feature = entry[4]
                print entry

#CollectCov(c,0)
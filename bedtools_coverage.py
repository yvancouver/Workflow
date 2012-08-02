import pybedtools

import sys

a = pybedtools.BedTool('/Volumes/ToveExFat/120615/040_KA005/020_refineAlignment/030_BQRecalGATK/all.realigned.markDup.baseQreCali.bam')

b = pybedtools.BedTool(sys.argv[1])

b.coverage(a,bam='abam')

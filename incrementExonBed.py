#!/usr/bin/python
# Aout 2012 Yvan Strahm yvan.strahm@gmail.com

import sys
import pybedtools

a = pybedtools.BedTool('/Users/yvans/Home/Dropbox/travail/BED_GFF_INTERVALS/Valided_And_Correct/CARDIO/Galaxy_GFFtoBED_EGS179.r150.readableregion_b37_sorted_GeneName.bed')
feature = "None"
for entry in  a[0:len(a)]:
	if feature != entry.name:
		i = 1
		print entry.chrom+"\t"+str(entry.start)+"\t"+str(entry.stop)+"\t"+entry.name+"_"+str(i)+"\t"+entry.score+"\t"+entry.strand
		feature = entry.name
	else:
		i += 1
		print entry.chrom+"\t"+str(entry.start)+"\t"+str(entry.stop)+"\t"+entry.name+"_"+str(i)+"\t"+entry.score+"\t"+entry.strand
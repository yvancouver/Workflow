#!/usr/bin/python
# Aout 2012 Yvan Strahm yvan.strahm@gmail.com
#'chrom', 'count', 'deparse_attrs', 'end', 'fields', 'file_type', 'length', 'name', 'o_amt', 'o_end', 'o_start', 'score', 'start', 'stop', 'strand'
import sys
import re
import pybedtools

seed = "None"
list_of_exons = []

def massage_the_list(list):
    print "new list"
    for i in list_of_exons:
            print str(i).rstrip()
    print len(list)
    print type(list)
    
# First get the length of the list in order to get the number of exon
    print len(list)
# determine if the gene is on the reverse or the forward stand
# search reverse
    if (re.search("_r$", entry.name)):
        #print "FOUNDrevere"
        print entry.name
    if (re.search("_f$", entry.name)):
        #print "FOUNDforward"
        print entry.name

    
bedfile = pybedtools.BedTool(sys.argv[1])
for entry in  bedfile[0:len(bedfile)]:
    #if (re.search(pattern+"_cds_[0-9_a-z]+_(r|f)", entry.name)):
    pattern = "(NM_"+re.escape(seed)+")"+"(_cds_)([0-9_a-z]+)(_r|f)"
    match = re.search(pattern, entry.name)
    if match:
#        print "all","NM_000364_cds_0_0_chr1_201328338_r"
        #print "0",match.group(0)
        list_of_exons.append(entry)
        
        massage_the_list(list_of_exons)
        list_of_exons = []
#        print "1",match.group(1)
#        print "2",match.group(2)
#        print "3",match.group(3)
#        print "4",match.group(4)
#        print match.group(5)
#        sys.exit()
    else:
        #print entry.name
        list_of_exons.append(entry)
        pre_pattern = entry.name.split("_")
        seed = pre_pattern[1]
        #sys.exit(pattern)
'''

entry = bedfile[0]
print "chr",entry.chrom,"start",entry.start,"stop",entry.stop,"name",entry.name,"score",entry.score,"strand",entry.strand,"o_start",entry.o_start,"o_amt",entry.o_amt
print dir(bedfile[0])
'''
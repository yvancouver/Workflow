## Get the chrom, start and stop
## From the example above we can see that chrom, start, stop, feature,score, and strand stay the same only posinfeat change
## One should return an line with summarizing the "strech" of the feature satisfying the wanted coverage.
## coverageBed -a test2.bed -b test1_1.bed -d
## should return
## chr    start    stop   exon in b  score strand start    stop     color code
## 7    10050    10060    Pos1_ex1    0    +    10050    10060    255,0,0
## 7    10080    10100    Pos1_ex1    0    +    10080    10100    255,0,0
## 7    10200    10250    Pos1_ex2    0    +    10200    10250    255,0,0
## 7    10300    10310    Pos1_ex2    0    +    10300    10310    255,0,0
## 7    10400    10500    Pos1_ex3    0    +    10400    10500    255,0,0

def CollectCov(coverageResult,cov):
    import pybedtools
    coverage_object = pybedtools.BedTool(coverageResult)
    chrom = ""
    start = 0.1
    stop = 0.1
    feature = "None"
    colors = "\t255,0,0"
    i = 1
    notSee = True
    for entry in  coverage_object[0:len(coverage_object)]:
        if int(entry[7]) <= cov:
            #print "F:t",feature == entry.name
            #print "So:t",int(entry[6]) == stop+1
            if feature != entry.name and stop != int(entry[6]):
                #print dir(entry)
                if feature == "None" and notSee:
# Need to remember the start of the region
# which is entry.start + start
                    notSee = False
                    feature = entry.name
                    chr_start = entry.start
                    chr_stop = entry.stop
                    strand = entry.strand
                    score = entry.score
                    start = int(entry[6])
                    stop = int(entry[6])
                    chrom = entry.chrom
                    print "1a\tline:",i,"\t",chrom,"\t",chr_start,"\t",chr_stop,"\t",feature,"\t",score,"\t",strand,"\t",stop                    
                    #continue
                elif feature != entry.name:
                    print "1b:\tline:",i,"\t",chrom,"\t",chr_start+start,"\t",chr_start+stop,"\t",feature,"\t",score,"\t",strand,"\t",chr_start+start,"\t",chr_start+stop,"\t",colors
                    print "1c:\tline:",i,"\t",entry
                    feature = entry.name
                    stop = int(entry[6])
                    print "1d:\tline:",i,"\t",chrom,"\t",chr_start+start,"\t",chr_start+stop,"\t",feature,"\t",score,"\t",strand,"\t",chr_start+start,"\t",chr_start+stop,"\t",colors
                    
                #print feature , entry.start, start
                i += 1
                #print"\t\tjumped line 51"
                continue
# if the feature is the same AND the base is next one just remember it
# but start is becoming stop
            if feature == entry.name and int(entry[6]) == stop+1 :
                stop = int(entry[6])
                print "2:\tline:",i,"\t", stop
                i += 1
                #print"\t\tjumped line 59"
                continue
# if the feature change or the base is not the next one
# print the previous and record the new values
            if feature == entry.name or int(entry[6]) != stop+1 :
                #print "3a_ori:\tline:",i-1,"\t",chrom,"\t",chr_start,"\t",chr_stop,"\t",feature,"\t",score,"\t",strand,"\t",stop
                print "3a:\tline:",i,"\t",chrom,"\t",chr_start+start,"\t",chr_start+stop,"\t",feature,"\t",score,"\t",strand,"\t",chr_start+start,"\t",chr_start+stop,"\t",colors
                #print "3b:\tline:",i,"\t",entry
                start = int(entry[6])
                stop = int(entry[6])
                #print "3c:\tline:",i,"\t",chrom,"\t",chr_start+start,"\t",chr_start+stop,"\t",feature,"\t",score,"\t",strand,"\t",chr_start+start,"\t",chr_start+stop,"\t",colors
                #print "jumped line 69"
                continue
            elif feature != entry.name or int(entry[6]) == stop+1 :
                #print "4a_ori:\tline:",i-1,"\t",chrom,"\t",chr_start,"\t",chr_stop,"\t",feature,"\t",score,"\t",strand,"\t",stop
                print "4a:\tline:",i,"\t",chrom,"\t",chr_start+start,"\t",chr_start+stop,"\t",feature,"\t",score,"\t",strand,"\t",chr_start+start,"\t",chr_start+stop,"\t",colors
                #print "4b:\tline:",i,"\t",entry
                feature = entry.name
                start = int(entry[6])
                stop = int(entry[6])
                print "4c:\tline:",i,"\t",chrom,"\t",chr_start+start,"\t",chr_start+stop,"\t",feature,"\t",score,"\t",strand,"\t",chr_start+start,"\t",chr_start+stop,"\t",colors

            i += 1
    last_entry = coverage_object[len(coverage_object)-1]
    print "last entry_ori:\t\t",last_entry.chrom,"\t",last_entry.start,"\t",last_entry.stop,"\t",last_entry.name,"\t",last_entry.score,"\t",last_entry.strand,"\t",start,"\t",stop,colors
            

## Get the chrom, start and stop
## From the example above we can see that chrom, start, stop, feature,score, and strand stay the same only posinfeat change
## One should return an line with summarizing the "strech" of the feature satisfying the wanted coverage.
## coverageBed -a test2.bed -b test1_1.bed -d
## should return check
## chr    start    stop   exon in b  score strand start    stop     color code
## 7    10050    10060    Pos1_ex1    0    +    10050    10060    255,0,0
## 7    10080    10100    Pos1_ex1    0    +    10080    10100    255,0,0
## 7    10200    10250    Pos1_ex2    0    +    10200    10250    255,0,0
## 7    10300    10310    Pos1_ex2    0    +    10300    10310    255,0,0
## 7    10400    10500    Pos1_ex3    0    +    10400    10500    255,0,0
#
# Need to record when the feature is the same BUT the base position has increased by more then 1
# or
# the feature has changed
#
# In the first case start and stop need to be reset (see condition 3)
#
# the second case requires the feature, start and stop to be record (condition 1)

def CollectCov(coverageResult,cov):
    import pybedtools
    coverage_object = pybedtools.BedTool(coverageResult)
    chrom = ""
    start = 0.1
    stop = 0.1
    feature = "None"
    colors = "255,0,0"
    i = 1
    notSee = True
    for entry in  coverage_object[0:len(coverage_object)]:
        if int(entry[7]) <= cov:
            if feature != entry.name and stop != int(entry[6]):
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
                    #continue
                elif feature != entry.name:
                    print chrom+"\t"+str(chr_start+start-1)+"\t"+str(chr_start+stop)+"\t"+feature+"\t"+score+"\t"+strand+"\t"+str(chr_start+start-1)+"\t"+str(chr_start+stop)+"\t"+colors
                    chrom = entry.chrom
                    feature = entry.name
                    chr_start = entry.start
                    start = int(entry[6])
                    stop = int(entry[6])
                    
                i += 1
                continue
# if the feature is the same AND the base is next one just remember it
# but start is becoming stop
            if feature == entry.name and int(entry[6]) == stop+1 :
                stop = int(entry[6])
                i += 1
                continue
# if the feature change or the base is not the next one
# print the previous and record the new values
            if feature == entry.name or int(entry[6]) != stop+1 :
                print chrom+"\t"+str(chr_start+start-1)+"\t"+str(chr_start+stop)+"\t"+feature+"\t"+score+"\t"+strand+"\t"+str(chr_start+start-1)+"\t"+str(chr_start+stop)+"\t"+colors
                start = int(entry[6])
                stop = int(entry[6])
                continue
            elif feature != entry.name or int(entry[6]) == stop+1 :
                print chrom,"HOHOHProblems\t",chr_start+start,"\t",chr_start+stop,"\t",feature,"\t",score,"\t",strand,"\t",chr_start+start-1,"\t",chr_start+stop,"\t",colors
                chrom = entry.chrom
                feature = entry.name
                start = int(entry[6])
                stop = int(entry[6])
                print "4c:\tline:",i,"\t",chrom,"\t",chr_start+start,"\t",chr_start+stop,"\t",feature,"\t",score,"\t",strand,"\t",chr_start+start,"\t",chr_start+stop,"\t",colors

            i += 1
    last_entry = coverage_object[len(coverage_object)-1]
    print last_entry.chrom+"\t"+str(last_entry.start+start-1)+"\t"+str(last_entry.start+stop)+"\t"+last_entry.name+"\t"+last_entry.score+"\t"+last_entry.strand+"\t"+str(last_entry.start+start-1)+"\t"+str(last_entry.start+stop)+"\t"+colors        
    
    #print "\t\t\t7    10051    10056    Pos1_ex1    0    +    10051    10056    255,0,0\n\t\t\t7    10081    10086    Pos1_ex1    0    +    10081    10086    255,0,0\n\t\t\t7    10206    10211    Pos1_ex2    0    +    10206    10211    255,0,0\n\t\t\t8    10406    10409    Pos1_ex3    0    +    10406    10409    255,0,0"
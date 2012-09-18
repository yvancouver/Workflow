## Get the chrom, start and stop
## From the example above we can see that chrom, start, stop, feature,score, and strand stay the same only posinfeat change
## One should return an line with summarizing the "strech" of the feature satisfying the wanted coverage.
## 7    10050    10059    Pos1_ex1    0    +    50    0
## 7    10081    10100    Pos1_ex1    0    +    80    0
## 7    10200    10250    Pos1_ex2    0    +    51    0
## 7    10300    10310    Pos1_ex2    0    +    151    0
## 7    10400    10500    Pos1_ex3    0    +    51    0


def CollectCov(coverageResult,cov):
    import pybedtools
    coverage_object = pybedtools.BedTool(coverageResult)
    #chrom = ""
    start = 0
    stop = 0
    feature = "None"
    colors = "\t255,0,0"

    for entry in  coverage_object[0:len(coverage_object)]:
        #if int(entry[7]) <= int(cov):
        print cov, feature
        if int(entry[7]) <= cov:
        
            if feature != entry.name:
## new feature so have to get his name and the start position (column 6)
                if feature == "None":
                    if entry.start == 1:
                        print "start entry",entry
                        start = entry.start+int(entry[6])-2
                    else:    
                        start = entry.start+int(entry[6])-1
                    feature = entry.name
                if feature != "None" and feature != entry.name:
                    #print chrom+"\t"+str(start)+"\t"+str(stop)+"\t"+feature+"\t"+entry.score+"\t"+entry.strand+"\t"+str(start)+"\t"+str(stop)+colors
                    if entry.start == 1:
                        start = entry.start+int(entry[6])-2
                    else:    
                        start = entry.start+int(entry[6])-1
                    feature = entry.name
            else:
                #chrom = entry[0]
                if entry.start == 1:
                    stop = entry.start + int(entry[6])-1
                else:
                    stop = entry.start + int(entry[6])
                
    last_entry = coverage_object[len(coverage_object)-1]
    print "last entry:",last_entry.chrom+"\t"+str(start)+"\t"+str(stop)+"\t"+last_entry.name+"\t"+last_entry.score+"\t"+last_entry.strand+"\t"+str(start)+"\t"+str(stop)+colors
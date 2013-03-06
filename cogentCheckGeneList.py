'''
Chromosome Name,Gene Start (bp),Gene End (bp),Associated Gene Name,RefSeq mRNA [e.g. NM_001195597],Strand,Band,Ensembl Gene ID,Ensembl Transcript ID,Transcript Start (bp),Transcript End (bp)
17,42976510,42982758,CCDC103,NM_213607,1,q21.31,ENSG00000167131,ENST00000417826,42977142,42982758
3,180320646,180588793,CCDC39,NM_181426,-1,q26.33,ENSG00000145075,ENST00000442201,180332376,180397288
17,78010435,78074412,CCDC40,NM_017950,1,q25.3,ENSG00000141519,ENST00000397545,78010435,78074412
16,84178865,84212373,DNAAF1,NM_178452,1,q23.3,ENSG00000154099,ENST00000378553,84178922,84211524
14,50091892,50101948,DNAAF2,NM_018139,-1,q21.3,ENSG00000165506,ENST00000298292,50091892,50101948
14,50091892,50101948,DNAAF2,NM_001083908,-1,q21.3,ENSG00000165506,ENST00000406043,50091892,50101948
19,55670031,55678090,DNAAF3,NM_178837,-1,q13.42,ENSG00000167646,ENST00000391719,55670031,55678018
7,21582833,21941457,DNAH11,NM_003777,1,p15.3,ENSG00000105877,ENST00000328843,21582833,21941457
5,13690440,13944652,DNAH5,NM_001369,-1,p15.2,ENSG00000039139,ENST00000265104,13690440,13944652
9,34457412,34520982,DNAI1,NM_012144,1,p13.3,ENSG00000122735,ENST00000242317,34458833,34520982
17,72270386,72311023,DNAI2,NM_001172810,1,q25.1,ENSG00000171595,ENST00000582036,72270386,72311023
17,72270386,72311023,DNAI2,NM_023036,1,q25.1,ENSG00000171595,ENST00000311014,72270429,72311023
14,74111578,74165604,DNAL1,NM_031427,1,q24.3,ENSG00000119661,ENST00000553645,74111702,74163197
7,766338,829190,HEATR2,NM_017802,1,p22.3,ENSG00000164818,ENST00000297440,766338,826112
16,70841281,71264625,HYDIN,NM_001270974,-1,q22.2,ENSG00000157423,ENST00000393567,70841281,71264592
8,133584320,133687838,LRRC6,NM_012472,-1,q24.22,ENSG00000129295,ENST00000250173,133584449,133687813
7,37888199,37940003,NME8,NM_016616,1,p14.1,ENSG00000086288,ENST00000199447,37888199,37940003
6,116937642,116954148,RSPH4A,NM_001010892,1,q22.1,ENSG00000111834,ENST00000229554,116937650,116954141
6,116937642,116954148,RSPH4A,NM_001161664,1,q22.1,ENSG00000111834,ENST00000368581,116937642,116954148
6,43612783,43640337,RSPH9,NM_001193341,1,p21.1,ENSG00000172426,ENST00000372165,43612783,43640336
6,43612783,43640337,RSPH9,NM_152732,1,p21.1,ENSG00000172426,ENST00000372163,43612783,43640336
19,48799714,48825130,CCDC114,NM_144577,-1,q13.33,ENSG00000105479,ENST00000315396,48799714,48823332
2,26624784,26679579,CCDC164,NM_145038,1,p23.3,ENSG00000157856,ENST00000288710,26624784,26679579
'''
import sys
import csv
import os

from cogent.db.ensembl import HostAccount, Species, Genome
import sqlalchemy as sql
#from sqlalchemy.engine import result
def getExons(gene):
    print gene
    #for ex in gene.Transcripts.Exons:
    #    print ex

def CheckRefSeq(ENSG,ENST,RefSeq,start,end):
    # gene fetch
    gene = human.getGeneByStableId(StableId=ENSG)

    # get the db tables we need
    external_db = human.CoreDb.getTable("external_db")
    object_xref = human.CoreDb.getTable("object_xref")
    xref = human.CoreDb.getTable("xref")

    # get the external db ID for refseq mrna
    refseq_mrna_id = sql.select([external_db.c.external_db_id],
                                external_db.c.db_name.like('RefSeq_mRNA')).execute().fetchone()
    
    # query for the entire lot
    print
    print "Querying for mRNA REFSEQ entries for all transcripts"
    # get all the transcript IDs
    tr_ids = [t._table_rows["transcript"][0] for t in gene.Transcripts]
    query = sql.select([object_xref, xref],
                       sql.and_(xref.c.xref_id==object_xref.c.xref_id,
                       object_xref.c.ensembl_id.in_(tr_ids),
                       xref.c.external_db_id == refseq_mrna_id[0]))
    
    results = query.execute().fetchall()
    print len(results)
    for res in results:
        #print dir(res)
        #print "row\t",type(res._row[8])
        #print result.keys
        #print result.values
        #print "\titems\t", res.items
        if str(RefSeq) != res._row[8]:
            #print RefSeq, " = ", res._row[8]
            print "NEXT"
            next
        else:
            print RefSeq, " = ", res._row[8]
            print "SAME\n Need to get exons"
            getExons(res)
            exit
#Define the Ensembl connections details
Release = 70

if 'ENSEMBL_ACCOUNT' in os.environ:
    host, username, password = os.environ['ENSEMBL_ACCOUNT'].split()
    account = HostAccount(host, username, password)
else:
    account = None

human = Genome(Species='human',Release=Release, account=account)
print human.Release
print human.GeneralRelease

# Get the csv lines
## Watch out for norwegian characters and empty space
with open('/Users/yvans/Home/workspace/Workflow/TestFiles/Pre_Primair_ciliaer_dyskinesi_mart_export.csv', 'rb') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in csvreader:
        if row[7] != "Ensembl Gene ID":
            print "Querying"
            print(row[7],row[8],row[4],row[1],row[2])
            CheckRefSeq(row[7],row[8],row[4],row[1],row[2])

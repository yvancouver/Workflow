import os
import sqlalchemy as sql
from cogent.db.ensembl import HostAccount, Genome

#account = HostAccount(*os.environ['ENSEMBL_ACCOUNT'].split())
if 'ENSEMBL_ACCOUNT' in os.environ:
    host, username, password = os.environ['ENSEMBL_ACCOUNT'].split()
    account = HostAccount(host, username, password)
else:
    account = None

human = Genome('human', Release=69, account=account)

# BRCA1
gene = human.getGeneByStableId(StableId="ENSG00000167131")

# get the db tables we need
external_db = human.CoreDb.getTable("external_db")
object_xref = human.CoreDb.getTable("object_xref")
xref = human.CoreDb.getTable("xref")

# get the external db ID for refseq mrna
refseq_mrna_id = sql.select([external_db.c.external_db_id],
                            external_db.c.db_name.like('RefSeq_mRNA')).execute().fetchone()


# query for a specific transcript ID
print "Querying for mRNA REFSEQ entries for one transcript"
query = sql.select([object_xref, xref],
        sql.and_(xref.c.xref_id==object_xref.c.xref_id,
        object_xref.c.ensembl_id == 1345831,
        xref.c.external_db_id == refseq_mrna_id[0]))

results = query.execute().fetchall()
print len(results)
for result in results:
    print result

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
for result in results:
    print result
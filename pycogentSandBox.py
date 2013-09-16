## Specifying a Host and Account

import os
Release = 71
from cogent.db.ensembl import HostAccount
if 'ENSEMBL_ACCOUNT' in os.environ:
     host, username, password = os.environ['ENSEMBL_ACCOUNT'].split()
     account = HostAccount(host, username, password)
else:
	account = None
print account

## What Species Are Available?

#from cogent.db.ensembl import Species
#print Species

## Interrogating a Genome

from cogent.db.ensembl import HostAccount, Genome
human = Genome(Species='human', Release=Release, account=account)
print human

'''
A Note on Coordinate Systems

The positions employed on Ensembls web-site, and in their MySQL database differ 
from those used internally by cogent.db.ensembl.In all cases where you are querying 
cogent.db.ensembl objects directly inputting nucleotide positions you can indicate 
you are using Ensembl coordinates by setting ensembl_coord=True.
If you are explicitly passing in a cogent.db.ensembl region, that argument has no effect.
'''

## Selecting Gene
#Via StableID
brca1 = human.getGeneByStableId(StableId='ENSG00000012048')
print brca1.Description

#Or gene symbol
genes = human.getGenesMatching(Symbol='brca2')
for gene in genes:
     if gene.Symbol.lower() == 'brca2':
         break

brca2 = gene # so we keep track of this reference for later on
print "Symbol\t", brca2.Symbol
print "Descr.\t", brca2.Description
print "gene\t", brca2
print "loc.\t", brca2.Location
print "length\t" ,len(brca2)

'''
Each location is directly tied to the parent genome and the coordinate above also 
shows the coordinates' type (chromosome in this case), name (13), start, end and strand. 
The start and end positions are python indices and will differ from the Ensembl indices 
in that start will be the Ensembl index - 1. This is because python counts from 0, not 1.
In querying for regions using a specific set of coordinates, it is possible to put 
in the Ensembl coordinates (demonstrated below).
'''

##Transcripts
'''
transcripts = brca2.Transcripts
for transcript in transcripts:
	print transcript
	print transcript.StableId
	print len(transcript.Exons)

transcript = brca2.CanonicalTranscript
for exon in transcript.Exons:
     print exon, exon.Location
'''


##Variations
#print human.getDistinct('effect')

effects =  human.getDistinct('effect')
for effect in effects:
	print effect
'''
nsyn_variants = human.getVariation(Effect='non_synonymous_codon',like=True)
for nsyn_variant in nsyn_variants:
	print "nsyn_variant\t",nsyn_variant
	break	#break will break out of the smallest enclosing loop
    
#print len(nsyn_variants)
'''
#brca2_snps = human.getFeatures(feature_types='synonymous_variant', region=brca2)
brca2_snps = human.getFeatures(region=brca2)
print brca2_snps
#for snp in brca2_snps:
#	print snp

brca2_snps = human.getFeatures(feature_types='variation',region=brca2)
for snp in brca2_snps:
     if 'non_synonymous_codon' in snp.Effect:
         break
print snp
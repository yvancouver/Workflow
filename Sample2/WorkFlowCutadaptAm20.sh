#!/bin/bash
#Where should I call the script???? Start in the main Folder above the passed_filter
##
# executables

# Declare files and binaries

export LOG="$PWD/$0_`date '+%F_%T'`.log"
echo `date` > timeUsed.txt
echo $LOG

# BWA
# Version: 0.5.9-r16
export BWA=/Users/yvans/Home/bin/bwa-0.5.9/bwa

#SAMTOOLS
#Version: 0.1.18 (r982:295)
export SAMTOOLS=/Users/yvans/Home/bin/samtools-0.1.18/samtools

##
# files and directories, This need to be edited, could this be passed by the config.cfg file?
##
export DB=/Users/yvans/Home/bin/GATK_resource_bundle_from_Ying_17_01_2012/1.2/b37/bwa_v5.10/human_g1k_v37_decoy.fasta
export GATK=/Users/yvans/Home/bin/GenomeAnalysisTK-1.3-24-gc8b1c92/
export PICARD=/Users/yvans/Home/bin/picard-tools-1.61/picard-tools-1.61/
export DBSNP=/Users/yvans/Home/bin/GATK_resource_bundle_from_Ying_17_01_2012/1.2/b37/dbsnp_132.b37.vcf
export HAPMAP=/Users/yvans/Home/bin/GATK_resource_bundle_from_Ying_17_01_2012/1.2/b37/hapmap_3.3.b37.sites.vcf
export OMNI=/Users/yvans//Home/bin/GATK_resource_bundle_from_Ying_17_01_2012/1.2/b37/1000G_omni2.5.b37.sites.vcf

#The reads which have passed QC are located in the passed_filter
# DIR is the top of the Reads containing folder
# WORKINGDIR is where this particular project will be done
# should the working directory name have datestamp?

export DIR=/Users/yvans/Home/Analysis/BRCA_analysis_01.12.2012/Sample_Diag-HaloBRCA1A-test-2
export WORKINGDIR=/Users/yvans/Home/Analysis/BRCA_analysis_01.12.2012/Sample_Diag-HaloBRCA1A-test-2/CutadaptAm20_Analysis/

# Are the reads from the default location or not? Should I passe dit by absolute path??
export READS1=/Users/yvans/Home/Analysis/BRCA_analysis_01.12.2012/Sample_Diag-HaloBRCA1A-test-2/passed_filter/CutATest2_paired_r1.fq
export READS2=/Users/yvans/Home/Analysis/BRCA_analysis_01.12.2012/Sample_Diag-HaloBRCA1A-test-2/passed_filter/CutATest2_paired_r2.fq

#
## RG line for the GATK consistency like that "@RG\tID:\tPL:ILLUMINA\tSM:"
## ID = Read group identier. Each @RG line must have a unique ID. The value of ID is used in the RG tags of alignment records. Must be unique among all read groups in header section. Read group IDs may be modied when merging SAM files in order to handle collisions.
## PL = Platform/technology used to produce the reads. Valid values: CAPILLARY, LS454, ILLUMINA, SOLID, HELICOS, IONTORRENT and PACBIO.
## SM = Sample. Use pool name where a pool is being sequenced.
#
export RG="@RG\tID:sample2CutAdaptAm20\tPL:ILLUMINA\tSM:sample2"

echo -e "`date`
Files and binaries for this analysis
BWA:\t\t$BWA
SAMTOOLS:\t$SAMTOOLS 
DB:\t\t$DB 
GATK:\t\t$GATK 
PICARD:\t\t$PICARD
DBSNP:\t\t$DBSNP
OMNI:\t\t$OMNI
HAPMAP:\t\t$HAPMAP
DIR:\t\t$DIR 
WORKINGDIR:\t$WORKINGDIR 
READS1:\t\t$READS1
READS2:\t\t$READS2 
RG:\t\t$RG

"> $LOG

#
## Create working directory
#


if [ ! -d $WORKINGDIR ] ; then
	mkdir -p $WORKINGDIR/010_alignment ;
	cd  $WORKINGDIR/010_alignment ;
else
	cd  $WORKINGDIR/010_alignment ;
fi


# Align Reads 1 against ref. genome
 echo -e "at `date`
 \taligning first bwa">>$LOG
 $BWA aln -t 3 $DB $READS1 > aln1.sai 2> aln1.log
 if [ ! -s aln1.sai ] ; then
 	echo "first index build failed check path and files" >>$LOG
 	exit
 fi
 
 # Align Reads 2
 echo -e "at `date`
 \taligning second bwa">>$LOG
 $BWA aln -t 3 $DB $READS2 > aln2.sai 2>aln2.log
 if [ ! -s aln2.sai ] ; then
 	echo "second index build failed check path and files" >>$LOG
 	exit
 fi
 echo -e "at `date`
 	\tfinished bwa aln">>$LOG
 	
 #
 ## Create sam and bam files
 ## Need to add read group with the -r option otherwise my version of GATK will complain
 
 #test if all the files are present
 echo -e "at `date`
 \t starting sampe job" >> $LOG
 if [[ ! -e aln1.sai || ! -e aln2.sai ]] ; then
 	echo " .sai files do not exist" >> $LOG ;
 	exit
 elif [ ! -e $DB ] ; then
 	echo " DB $DB does not exist" >> $LOG ;
 	exit
 else 
 	$BWA sampe -r $RG $DB aln1.sai aln2.sai $READS1 $READS2 >aln.sam 2> sampe.log;
 	echo -e "at `date`
 	\tfinished sampe job" >> $LOG
 fi


#
## For now I am using samtools but should I use picard to be more carefull?????
## But what picard commands? It is also recommended by the GATK to use Picard 
## SortSam.jar SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT INPUT=file.sam OUTPUT=file.bam??

## Creating a bam file
#
if [[ ! -e aln.sam || ! -s aln.sam ]] ;then 
	echo "aln.sam do not exist or is empty have a look in sampe.log" >> $LOG ;
	exit
elif [ ! -e $DB.fai ] ; then
	echo "$DB.fai could not be found" >> $LOG ;
	exit
else
	echo -e "at `date`
	\t starting bam building and indexing" >> $LOG ;
	java -Xmx4g -jar $PICARD/SortSam.jar \
	I=aln.sam \
	O=aln.posiSrt.bam \
	SO=coordinate \
	VALIDATION_STRINGENCY=SILENT \
	2>SortSam.txt

	java -Xmx4g -jar $PICARD/BuildBamIndex.jar \
	INPUT=aln.posiSrt.bam \
	2>BuildBamIndex.txt
fi

#
## cleaning steps
#
rm -rf aln1.sai aln2.sai aln.sam


#
## Calculate mapped reads
#
if [[ ! -e aln.posiSrt.bam || ! -s aln.posiSrt.bam ]] ; then
	echo "the aln.posiSrt.bam file is missing or empty" >> $LOG ;
	exit
else
	echo -e "at `date`
	\tproducing flagstats" >> $LOG ;
	$SAMTOOLS flagstat aln.posiSrt.bam > flagstats.txt
fi

#
## Check is mapping worked
#
echo -e "at `date`
	\tproducing mapped.txt" >> $LOG ;

$SAMTOOLS view -X aln.posiSrt.bam | awk 'BEGIN{unmapped=0; unique=0; ambiguous=0; FS="\t"};\
($1 !~ /^@/) \
{if($2 !~ /.*u.*/){if($5==0){ambiguous=ambiguous+1}else{unique=unique+1}}else{unmapped=unmapped+1}};\
END{total=ambiguous + unique + unmapped; \
print "Unmapped: " unmapped; \
print "Mapping uniquely: " unique; \
print "Mapping ambiguously: "  ambiguous; \
print "Total: "  total;}'> mapped.txt

if [[ ! -s mapped.txt || ! -e mapped.txt ]] ; then 
	echo "Something wrong happened with mapped.txt, Have a look " >> $LOG;
	exit
fi
	
#Create a pdf of insert size

if [ ! -e aln.posiSrt.bam ] ; then
	echo "No more aln.posiSrt.bam check your files" >> $LOG ;
	exit 
else
	echo -e "at `date`
	\tproducing CollectInsertSizeMetrics.log" >> $LOG ;
	java -jar $PICARD/CollectInsertSizeMetrics.jar \
	INPUT=aln.posiSrt.bam \
	OUTPUT=insertSizeMetrics.txt \
	HISTOGRAM_FILE=insertSizeHistogram.pdf \
	2>CollectInsertSizeMetrics.txt
fi

if [ ! -e insertSizeMetrics.txt ] ; then
	echo "Check CollectInsertSizeMetrics.txt, something went wrong with CollectInsertSizeMetrics" >> $LOG;
	exit
fi

## Realign the reads 
# Create targets
# If everything went ok then the errRealignerTargetCreator should be empty (zero byte)

echo -e "at `date`
	\tproducing aln.posiSrt.intervals" >> $LOG ;
java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R $DB \
-o aln.posiSrt.intervals \
-I aln.posiSrt.bam \
--known $DBSNP \
2>errRealignerTargetCreator > realignerTargetCreatorInfo.txt

#
## Test if the error file is not empty
## -s == NOT empty
## ! -s == empty
#

if [ -s errRealignerTargetCreator ] ; then
	echo "RealignerTargetCreator produced some errors please have a look at the errRealignerTargetCreator and realignerTargetCreatorInfo.txt" >>$LOG;
	exit
fi

#
## realign the reads
## need to use the maxreads here because of the number of reads
## If everything went ok then the errIndelRealigner should be empty (zero byte)
#
echo -e "at `date`
	\tproducing aln.posiSrt.realigned.bam" >> $LOG ;
java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar  \
-T IndelRealigner \
-maxReads 1000000 \
-I ../010_alignment/aln.posiSrt.bam \
-R $DB \
-targetIntervals aln.posiSrt.intervals \
-o aln.posiSrt.realigned.bam \
-compress 0 \
2>errIndelRealigner > indelRealignerInfo.txt

# Test if the error file is empty
if [ -s errIndelRealigner ] ; then
	echo "IndelRealigner produced some errors please have a look at the errIndelRealigner and indelRealignerInfo.txt" >>$LOG;
	exit
fi

## realigning
# Normally errFixMate contains INFO on the realignent and if something wet wrong the word exeption should appear
# Should I stop the analysis if the test failed?
echo -e "at `date`
	\tproducing aln.posiSrt.fixedMate.bam" >> $LOG ;
java -Xmx4g -jar $PICARD/FixMateInformation.jar \
INPUT=aln.posiSrt.realigned.bam \
OUTPUT=aln.posiSrt.fixedMate.bam \
SO=coordinate \
VALIDATION_STRINGENCY=SILENT \
CREATE_INDEX=true \
2>errFixMate


# I don't know how to catch an error produced by PICARD
# Try
if grep -E "ERROR|Exception" errFixMate > /dev/null ; then 
	echo "FixMateInformation reported an error have a look at ErrFixMate"; 
	exit
fi;

##### Base quality score recalibration

## count covariates before recalibration
# errCountCovariatesPre should be empty or containing ERROR
echo -e "at `date`
	\tproducing table.recal_data.pre.csv" >> $LOG ;
java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar \
-l INFO \
-R $DB \
-knownSites $DBSNP \
-I aln.posiSrt.fixedMate.bam \
-T CountCovariates \
-cov ReadGroupCovariate \
-cov QualityScoreCovariate \
-cov DinucCovariate \
-cov CycleCovariate \
-recalFile table.recal_data.pre.csv \
-nt 3 \
2>errCountCovariatesPre > countCovariatesPreInfo.txt

if [ -s errCountCovariatesPre ] ; then 
	echo "CountCovariates PRE produced some errors please have a look at the errCountCovariatesPre and countCovariatesPreInfo.txt" >>$LOG;
	exit
fi

## recalibration
echo -e "at `date`
	\tproducing aln.posiSrt.baseQreCali.bam" >> $LOG ;
java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar \
-l INFO \
-R $DB \
-I aln.posiSrt.fixedMate.bam \
-T TableRecalibration \
--out aln.posiSrt.baseQreCali.bam \
-recalFile table.recal_data.pre.csv \
-OQ \
-pQ 5 \
2>errTableRecalibration > tableRecalibrationInfo.txt

if [ -s errTableRecalibration ] ; then 
	echo "TableRecalibration produced some errors please have a look at the errTableRecalibration and tableRecalibrationInfo.txt" >>$LOG;
	exit
fi

echo -e "\tHow will samtools name this index"
$SAMTOOLS index aln.posiSrt.baseQreCali.bam


## count covariates after recalibration
echo -e "at `date`
	\tproducing aln.posiSrt.baseQreCali.bam" >> $LOG ;
java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar \
-l INFO \
-R $DB \
-knownSites $DBSNP \
-I aln.posiSrt.baseQreCali.bam \
-T CountCovariates \
-cov ReadGroupCovariate \
-cov QualityScoreCovariate \
-cov DinucCovariate \
-cov CycleCovariate \
-recalFile table.recal_data.post.csv \
-nt 3 \
2>errCountCovariatesPost > countCovariatesPostInfo.txt

if [ -s errCountCovariatesPost ] ; then 
	echo "CountCovariates Post produced some errors please have a look at the errCountCovariatesPost and countCovariatesPostInfo.txt" >>$LOG;
	exit
fi

## check the improvement
mkdir pre_recal
mkdir post_recal

echo -e "at `date`
	\tproducing plotPreInfo.txt" >> $LOG ;
java -Xmx4g -jar $GATK/AnalyzeCovariates.jar \
-recalFile table.recal_data.pre.csv \
-outputDir pre_recal/   \
-ignoreQ 5 \
2>errPlotpre > plotPreInfo.txt

if [ -s errPlotpre ] ; then 
	echo "AnalyzeCovariates produced some errors please have a look at the errPlotpre and plotPreInfo.txt" >>$LOG;
	exit
fi

echo -e "at `date`
	\tproducing  plotPostInfo.txt" >> $LOG ;
java -Xmx4g -jar $GATK/AnalyzeCovariates.jar \
-recalFile table.recal_data.post.csv \
-outputDir post_recal/ \
-ignoreQ 5 \
2>errPlotpost > plotPostInfo.txt

if [ -s errPlotpost ] ; then 
	echo "AnalyzeCovariates produced some errors please have a look at the errPlotpost and plotPostInfo.txt" >>$LOG;
	exit
fi

#This need to be fixed??? Which bed file should I used?? See With Jungbai If he is willing to share.
################################ checking coverage ################################

# Getting a simple coverage stat on coverage
echo -e "at `date`
	\tproducing  coverage data" >> $LOG ;
coverageBed -d -abam aln.posiSrt.baseQreCali.bam -b ~/Home/Dropbox/travail/BRCA12/BRCA.hg19.end.ampregion.bed >coverageEachBase.bed
coverageBed -hist -abam aln.posiSrt.baseQreCali.bam -b ~/Home/Dropbox/travail/BRCA12/BRCA.hg19.end.ampregion.bed >coverageEachBase.hist.bed

# Using picard
echo -e "at `date`
	\tproducing hsMetrics.txt" >> $LOG ;
java -jar $PICARD/CalculateHsMetrics.jar \
INPUT=aln.posiSrt.baseQreCali.bam \
OUTPUT=hsMetrics.txt \
BAIT_INTERVALS=/Users/yvans/Home/Dropbox/travail/BRCA12/BRCA.hg19.end.ampregion.intervals  \
TARGET_INTERVALS=/Users/yvans/Home/Dropbox/travail/BRCA12/BRCA.hg19.end.ampregion.intervals  \
2> errCalculateHsMetrics > CalculateHsMetrics.txt

if grep -E "ERROR|Exception" errCalculateHsMetrics > /dev/null ; then 
	echo "FixMateInformation reported an error have a look at errCalculateHsMetrics and CalculateHsMetrics.txt"; 
	exit
fi;
# Again you will want to get this into a readable format

paste <(cat hsMetrics.txt | tail -n 4 | head -n 1 | tr "\t" "\n") <(cat hsMetrics.txt | tail -n 3 | head -n 1 | tr "\t" "\n")>hsMetricsReformated.txt

mkdir $WORKINGDIR/020_realignment
cd $WORKINGDIR/020_realignment

#Part2
################################ Variant Calling ################################

## snp
echo -e "at `date`
	\tproducing snps.raw.vcf" >> $LOG ;
java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar \
-R $DB \
-T UnifiedGenotyper \
-I ../010_alignment/aln.posiSrt.baseQreCali.bam \
--dbsnp $DBSNP \
-o snps.raw.vcf \
-stand_call_conf 30 \
-stand_emit_conf 4 \
-baq CALCULATE_AS_NECESSARY \
-L /Users/yvans/Home/Dropbox/travail/BRCA12/BRCA.hg19.end.ampregion.intervals \
2>errSnpCalling > snpCallingInfo.txt

if [ -s errSnpCalling ] ; then 
	echo "UnifiedGenotyper produced some errors please have a look at the errSnpCalling and snpCallingInfo.txt" >>$LOG;
	exit
fi

## indel
echo -e "at `date`
	\tproducing indel.raw.vcf " >> $LOG ;
java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar \
-R $DB \
-T UnifiedGenotyper \
-glm INDEL \
-I ../010_alignment/aln.posiSrt.baseQreCali.bam \
--dbsnp ~/Home/bin/GATK_resource_bundle_from_Ying_17_01_2012/1.2/b37/dbsnp_132.b37.vcf \
-o indel.raw.vcf \
-L /Users/yvans/Home/Dropbox/travail/BRCA12/BRCA.hg19.end.ampregion.intervals \
2>errIndelCalling > indelCallingInfo.txt

if [ -s errIndelCalling ] ; then 
	echo "UnifiedGenotyper INDEL produced some errors please have a look at the errIndelCalling and indelCallingInfo.txt" >>$LOG;
	exit
fi

# Check the result files
# - How many variants are called:
echo -e "at `date`
	\tproducing grep output " >> $LOG ;
grep -v "#" snps.raw.vcf |wc -l > snps.raw.vcf.out
if [ ! -s snps.raw.vcf.out ] ; then
	echo "WARNING snps.raw.vcf.out is empty " >>$LOG;
fi
grep -v "#" indel.raw.vcf |wc -l > indel.raw.vcf.out

if [ ! -s indel.raw.vcf.out ] ; then
	echo "WARNING indel.raw.vcf.out is empty " >>$LOG;
fi
# -How many variants are passed basic filter
grep -v "#" snps.raw.vcf |grep PASS|wc -l > snps.raw.vcf.pass.out
if [ ! -s snps.raw.vcf.pass.out ] ; then 
	echo "WARNING snps.raw.vcf.pass.out is empty " >>$LOG;
fi
grep -v "#" indel.raw.vcf |grep PASS|wc -l > indel.raw.vcf.pass.out
if [ ! -s indel.raw.vcf.pass.out ] ; then
	echo "WARNING indel.raw.vcf.pass.out is empty " >>$LOG;
fi
# -How many variants are reported in dbSNP
grep –vP "^#" *.vcf | grep –cP "\trs" > variantsDBSNPs.out
if [ ! -s variantsDBSNPs.out ] ; then
	echo "WARNING variantsDBSNPs.out is empty " >>$LOG;
fi

################################ Variant hard filtration ################################

## snp
echo -e "at `date`
	\tproducing snp.filter.vcf " >> $LOG ;
java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar \
-R $DB \
-T VariantFiltration \
-o snp.filter.vcf \
--variant snps.raw.vcf \
--filterExpression "QD < 2.0" \
--filterExpression "MQ < 40.0" \
--filterExpression "FS > 60.0" \
--filterExpression "HaplotypeScore > 13.0" \
--filterExpression "MQRankSum < -12.5" \
--filterExpression "ReadPosRankSum < -8.0" \
--filterName "QDFilter" \
--filterName "MQFilter" \
--filterName "FSFilter" \
--filterName "HaplotypeScoreFilter" \
--filterName "MQRankSumFilter" \
--filterName "ReadPosRankSumFilter" \
2>errSnpFilter > snpFilterInfo.txt
if [ -s errSnpFilter ] ; then
	echo "VariantFiltration produced some errors please have a look at the errSnpFilter and snpFilterInfo.txt" >>$LOG;
	exit
fi

## indel
echo -e "at `date`
	\tproducing indel.filter.vcf " >> $LOG ;
java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R $DB \
-o indel.filter.vcf \
--variant indel.raw.vcf \
--filterExpression "QD < 2.0" \
--filterExpression "ReadPosRankSum < -20.0" \
--filterExpression "InbreedingCoeff < -0.8" \
--filterExpression "FS > 200.0" \
--filterName "QDFilter" \
--filterName "ReadPosRankSumFilter" \
--filterName "InbreedingCoeffFilter" \
--filterName "FSFilter" \
2>errIndelFilter > indelFilterInfo.txt

if [ -s errIndelFilter ] ; then
	echo "VariantFiltration produced some errors please have a look at the errIndelFilter and indelFilterInfo.txt" >>$LOG;
	exit
fi

############################ variant quality score recalibration ################################
## to little data to run this file
## Variant Recalibrator
## for halo brca too few targets to do it
# 
# echo -e "at `date`
# 	\tproducing variantRecalibratorInfo.txt  " >> $LOG ;
# java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar \
# -T VariantRecalibrator \
# -R $DB \
# -input snps.raw.vcf \
# -resource:hapmap,known=false,training=true,truth=true,prior=15.0 /Users/yvans/Home/bin/GATK_resource_bundle_from_Ying_17_01_2012/1.2/b37/hapmap_3.3.b37.sites.vcf \
# -resource:omni,known=false,training=true,truth=false,prior=12.0 /Users/yvans//Home/bin/GATK_resource_bundle_from_Ying_17_01_2012/1.2/b37/1000G_omni2.5.b37.sites.vcf \
# -resource:dbsnp,known=true,training=false,truth=false,prior=8.0 /Users/yvans/Home/bin/GATK_resource_bundle_from_Ying_17_01_2012/1.2/b37/dbsnp_132.b37.vcf \
# -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
# -recalFile output.recal \
# -tranchesFile output.tranches \
# -rscriptFile output.plots.R \
# --maxGaussians 4 \
# --target_titv 3.2 \
# 2>errVariantRecalibrator > variantRecalibratorInfo.txt 
# 
# if [ -s errVariantRecalibrator ] ; then
# 	echo "VariantRecalibrator produced some errors please have a look at the errVariantRecalibrator and variantRecalibratorInfo.txt" >>$LOG;
# 	exit
# fi
# 
# ## Apply Recalibrator
# echo -e "at `date`
# 	\tproducing snps.filter.vcf " >> $LOG ;
# java -Xmx4g -jar  $GATK/GenomeAnalysisTK.jar \
# -T ApplyRecalibration \
# -R $DB \
# -input snps.raw.vcf \
# --ts_filter_level 99.0 \
# -tranchesFile output.tranches \
# -recalFile output.recal \
# -o snps.filter.vcf \
# > applyRecalibratorInfo.txt

#if [ -s errVariantRecalibrator ] then ;
#	echo "VariantRecalibrator produced some errors please have a look at the errVariantRecalibrator and variantRecalibratorInfo.txt" >>$LOG;
#	exit
#fi
############################ variant evaluation by GATK ################################
echo -e "at `date`
	\tproducing output.eval.gatkreport" >> $LOG ;
java -Xmx2g -jar GenomeAnalysisTK.jar \
-R $DB \
-T VariantEval 
-o output.eval.gatkreport \
--eval:set1 snp.filter.vcf \
--dbsnp $DBSNP

############################ variant annotation by annovar ################################
echo -e "at `date`
	\tstarting annovar" >> $LOG ;
perl ./convert2annovar.pl snp.filter.vcf -format vcf4 -includeinfo > snps.filter.avinput
./summarize_annovar.pl snps.filter.avinput -buildver hg19 -verdbsnp 132 /Users/yvans/Home/bin/annovar_2011Sep11/humandb -outfile sumSNP

echo -e "at `date`
	\tfinished annovar" >> $LOG ;

echo `date` >> timeUsed.txt
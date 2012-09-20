##
# Should the different DB and soft should be declared in this script
##

# Cutadapt
# Version ?
CUT=/volumm
mCUT=32
adapter1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
adapter2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT


# BWA
# Version: 0.5.9-r26-dev
export BWA=/Volumes/data.odin/software/mac/bwa/bwa-0.5.10

#SAMTOOLS
#Version: 0.1.18 (r982:295)
export SAMTOOLS=/Volumes/data.odin/software/mac/samtools/samtools-0.1.18/samtools

export DB=/Users/yvans/Home/bin/GATK_resource_bundle_from_Ying_17_01_2012/1.2/b37/bwa_v5.10/human_g1k_v37_decoy.fasta
export GATK=/Users/yvans/Home/bin/GenomeAnalysisTK-1.4-37-g0b29d54/
export PICARD=/Volumes/data.odin/software/mac/picard/picard-tools-1.62/picard-tools-1.62/
export DBSNP=/Volumes/data.odin/common/GATK_resource_bundle/1.2/b37/dbsnp_132.b37.vcf
export HAPMAP=/Volumes/data.odin/common/GATK_resource_bundle/1.2/b37/hapmap_3.3.b37.sites.vcf
export OMNI=/Volumes/data.odin/common/GATK_resource_bundle/1.2/b37/1000G_omni2.5.b37.sites.vcf
export Mills_1000G_gold=/Volumes/data.odin/common/GATK_resource_bundle/1.2/b37/Mills_and_1000G_gold_standard.indels.b37.sites.vcf
export ANNOVAR=/Volumes/data.odin/software/mac/annovar/annovar_2012May25/

DB=/Volumes/data.odin/common/GATK_resource_bundle/1.2/b37/bwa_v5.10/human_g1k_v37_decoy.fasta
# What version should we use??
GATK=/Users/yvans/Home/bin/GenomeAnalysisTK-1.4-37-g0b29d54/
#
PICARD=/Volumes/data.odin/software/mac/picard/picard-tools-1.62/picard-tools-1.6
DBSNP=/Volumes/data.odin/common/GATK_resource_bundle/1.2/b37/dbsnp_132.b37.vcf"
HAPMAP=/Users/yvans/Home/bin/GATK_resource_bundle_from_Ying_17_01_2012/1.2/b37/hapmap_3.3.b37.sites.vcf
OMNI=/Users/yvans//Home/bin/GATK_resource_bundle_from_Ying_17_01_2012/1.2/b37/1000G_omni2.5.b37.sites.vcf
Mills_Gold_1000=/Volumes/data.odin/common/GATK_resource_bundle/1.2/b37/

# hsMetrics
BAIT=
TARGET=

RG="@RG\tID:1\tPL:ILLUMINA\tSM:"$PREFIX"_paired_a_m32"
# LOG
export LOG="$PWD/$0_`date '+%F_%T'`.log"
export TimeUsed=timeUsed.txt
echo `date` > $TimeUsed
echo $LOG

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
ANNOVAR:\t\t$ANNOVAR
Mills_1000g:\t\t$Mils_1000G_gold
DIR:\t\t$DIR 
WORKINGDIR:\t$WORKINGDIR 
READS1:\t\t$READS1
READS2:\t\t$READS2 
RG:\t\t$RG

"> $LOG

READS1=
READS2=


date

# test if the two files are different. Other test should follow, like are they really paired end file 
# XXXXX_R1_XXX.fastq.gz and XXXX_R2_XXX.fastq.gz

if [ $READS1 != $READS2 ] ; then  echo NOT  ; fi

echo $READS1
echo $READS2

echo

echo "cutadapt r1"
time $CUT -m $mCUT -a $adapter1 \
$READS1 > ${READS1%.fastq.gz}_Cutadapt_a_m32.fastq \
2> ${READS1%.fastq.gz}_cutadapt_R1_a_m32.log

echo
echo "cutadapt r2"
time $CUT -m $mCUT  -a $adapter2 \
$READS2 > ${READS2%.fastq.gz}_Cutadapt_a_m32.fastq \
2> ${READS2%.fastq.gz}_cutadapt_R2_a_m32.log

#dont forget to change the line number

export PREFIX=${READS1%_L00X_R1_001.fastq.gz}

echo "$PREFIX"

echo
echo "syncing"
/Users/yvans//Home/workspace/Workflow/TimSync.sh \
${READS1%.fastq.gz}_Cutadapt_a_m32.fastq \
${READS2%.fastq.gz}_Cutadapt_a_m32.fastq \
$PREFIX > ${READS2%.fastq.gz}_TimSync.log

gzip ${READS1%.fastq.gz}_Cutadapt_a_m32.fastq
gzip ${READS2%.fastq.gz}_Cutadapt_a_m32.fastq

##
# don't forget to add some lines in order to clean after the sync, delete or compress files
##

echo
echo "mapping r1"
time $BWA aln -t 3 $DB \
"$PREFIX"_paired_r1.fq >"$PREFIX"_paired_r1.sai \
2> "$PREFIX"_paired_r1.sai.log

echo
echo "mapping r2"
time $BWA aln -t 3 $DB \
"$PREFIX"_paired_r2.fq > "$PREFIX"_paired_r2.sai \
2> "$PREFIX"_paired_r2.sai.log

echo
echo "pairing"
time $BWA sampe -r $RG $DB \
"$PREFIX"_paired_r1.sai "$PREFIX"_paired_r2.sai "$PREFIX"_paired_r1.fq "$PREFIX"_paired_r2.fq > "$PREFIX"_paired.sam 2>sampe.log

gzip "$PREFIX"_paired_r1.fq
gzip "$PREFIX"_paired_r2.fq

echo
echo "Picard sort"
time java -Xmx8g -jar $PICARD/SortSam.jar \
I="$PREFIX"_paired.sam \
O="$PREFIX"_paired.bam \
SO=coordinate 2>errSortSam

echo "Check for the bam file validity"

java -Xmx8g -jar $PICARD/ValidateSamFile.jar \
INPUT="$PREFIX"_paired.bam \
MODE=SUMMARY \
MAX_OUTPUT=null
>first_bam_validity.txt

echo
echo "Picard Fixmate"
time java -Xmx8g -jar $PICARD/FixMateInformation.jar \
I="$PREFIX"_paired.bam \
O="$PREFIX"_pairedFixed.bam \
2>errFixMate

rm -f "$PREFIX"_paired.sam
rm -f "$PREFIX"_paired.bam
rm -f "$PREFIX"_paired_r1.fq
rm -f "$PREFIX"_paired_r2.fq

echo
echo "Picard indexing"
time java -Xmx8g -jar $PICARD/BuildBamIndex.jar \
INPUT="$PREFIX"_pairedFixed.bam \
VALIDATION_STRINGENCY=SILENT \
2>errBuildBamIndex

echo
echo "flagstat"
time $SAMTOOLS flagstat "$PREFIX"_pairedFixed.bam > flagstats.txt

echo
echo "mapped.txt"
time /Users/yvans/Home/bin/samtools-0.1.18/samtools view -X "$PREFIX"_pairedFixed.bam | awk 'BEGIN{unmapped=0; unique=0; ambiguous=0; FS="\t"};\
	($1 !~ /^@/) \
	{if($2 !~ /.*u.*/){if($5==0){ambiguous=ambiguous+1}else{unique=unique+1}}else{unmapped=unmapped+1}};\
	END{total=ambiguous + unique + unmapped; \
	print "Unmapped: " unmapped; \
	print "Mapping uniquely: " unique; \
	print "Mapping ambiguously: "  ambiguous; \
	print "Total: "  total;}'> mapped.txt

echo
echo "InsertMetrics"
time java -Xmx8g -jar $PICARD/CollectInsertSizeMetrics.jar \
INPUT="$PREFIX"_pairedFixed.bam \
OUTPUT=insertSizeMetrics.txt \
HISTOGRAM_FILE=insertSizeHistogram.pdf \
2>CollectInsertSizeMetrics.txt

echo
echo "mark duplicates"
time java -Xmx8g -jar $PICARD/MarkDuplicates.jar \
INPUT="$PREFIX"_pairedFixed.bam \
OUTPUT="$PREFIX"_pairedFix.dup.bam \
METRICS_FILE=duplicates_metrics.txt \
REMOVE_DUPLICATES=false 2>errMarkDuplicates

## Remove dup bam
rm -f "$PREFIX"_pairedFix.dup.bam

echo
echo "Realign"
time java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R $DB \
-o "$PREFIX"_paired.intervals \
-I "$PREFIX"_pairedFixed.bam \
--known $DBSNP \
> realignerTargetCreatorInfo.txt 2>errRealignerTargetCreator 

echo
echo "Indel realigner"
time java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar \
-T IndelRealigner \
-maxReads 1000000 \
-I "$PREFIX"_pairedFixed.bam \
-R $DB \
-targetIntervals "$PREFIX"_paired.intervals \
-o "$PREFIX"_paired.realigned.bam \
-compress 0  > indelRealignerInfo.txt 2>errIndelRealigner


rm -f "$PREFIX"_pairedFixed.bam

echo
echo "count variant before"
time java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar \
-l INFO \
-R $DB \
-knownSites $DNSNP \
-I "$PREFIX"_paired.realigned.bam \
-T CountCovariates -cov ReadGroupCovariate \
-cov QualityScoreCovariate \
-cov DinucCovariate \
-cov CycleCovariate \
-recalFile table.recal_data.pre.csv \
-nt 3  > countCovariatesPreInfo.txt 2>errCountCovariatesPre

echo
echo " recalibration"
time java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar \
-l INFO \
-R $DB \
-I "$PREFIX"_paired.realigned.bam \
-T TableRecalibration \
--out "$PREFIX"_paired.baseQreCali.bam \
-recalFile table.recal_data.pre.csv -OQ -pQ 5 \
> tableRecalibrationInfo.txt 2>errTableRecalibration

rm -f "$PREFIX"_paired.realigned.bam

echo
echo "count variant after"
time java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar \
-l INFO \
-R $DB \
-knownSites $DBSNP \
-I "$PREFIX"_paired.baseQreCali.bam \
-T CountCovariates \
-cov ReadGroupCovariate \
-cov QualityScoreCovariate \
-cov DinucCovariate \
-cov CycleCovariate \
-recalFile table.recal_data.post.csv \
-nt 3  > countCovariatesPostInfo.txt 2>errCountCovariatesPost

echo "Check for the final bam file validity"

java -Xmx8g -jar $PICARD/ValidateSamFile.jar \
INPUT="$PREFIX"_paired.baseQreCali.bam \
MODE=SUMMARY \
MAX_OUTPUT=null
>final_bam_validity.txt

echo
echo "create pre and post folder"
mkdir pre_recal
mkdir post_recal

echo
echo "evaluation of recalibration"
echo
echo "pre"
time java -Xmx8g -jar $GATK/AnalyzeCovariates.jar \
-recalFile table.recal_data.pre.csv \
-outputDir pre_recal/ \
-ignoreQ 5 2>errPlotpre > plotPreInfo.txt

echo
echo "post"
time java -Xmx8g -jar $GATK/AnalyzeCovariates.jar \
-recalFile table.recal_data.post.csv \
-outputDir post_recal/ \
-ignoreQ 5 2>errPlotpost > plotPostInfo.txt

echo
echo "hsMetrics"
time java -Xmx8g -jar $PICARD/CalculateHsMetrics.jar \
INPUT="$PREFIX"_paired.baseQreCali.bam \
OUTPUT=hsMetrics.txt \
TARGET_INTERVALS=/Users/yvans/Home/Analysis/Cardio_analysis_01.17.2012/SingleSample/010_alignment/CardioAllRefSeqExons_ucsc_corr20.interval_list  \
BAIT_INTERVALS=/Users/yvans/Home/Analysis/Cardio_analysis_01.17.2012/SingleSample/010_alignment/CardioEGS179.r150.readableregion_corr_sort.interval_list \
2> errCalculateHsMetrics > CalculateHsMetrics.txt

echo
echo "reformat hsmetrics"
/Users/yvans/Home/workspace/Workflow/reFormathsMetrics.sh hsMetrics.txt >hsMetrics.csv

echo
echo "unified genotyper SNPs"
time java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar \
-R $DB \
-T UnifiedGenotyper \
-glm SNP \
-A AlleleBalance \
-I "$PREFIX"_paired.baseQreCali.bam \
--dbsnp /Users/yvans/Home/bin/GATK_resource_bundle_from_Ying_17_01_2012/1.2/b37/dbsnp_132.b37.vcf \
-o snps.raw.vcf \
-stand_call_conf 50 \
-stand_emit_conf 10 \
-dcov 1000 \
-L /Users/yvans/Home/Analysis/Cardio_analysis_01.17.2012/SingleSample/010_alignment/CardioAllRefSeqExons_ucsc_corr20.interval_list \
> snpCallingInfo.txt  2>errSnpCalling

echo
echo "unified genotyper INDELs"
time java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar \
-R $DB \
-T UnifiedGenotyper \
-glm INDEL \
-A AlleleBalance \
-I "$PREFIX"_paired.baseQreCali.bam \
--dbsnp /Users/yvans/Home/bin/GATK_resource_bundle_from_Ying_17_01_2012/1.2/b37/dbsnp_132.b37.vcf \
-o indel.raw.vcf \
-stand_call_conf 50 \
-stand_emit_conf 10 -dcov 1000 \
-L /Users/yvans/Home/Analysis/Cardio_analysis_01.17.2012/SingleSample/010_alignment/CardioAllRefSeqExons_ucsc_corr20.interval_list \
> indelCallingInfo.txt  2>errIndelCalling

echo
echo "hard filtration SNPs"
time java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar \
-R $DB \
-T VariantFiltration \
-o snps.filter.vcf \
--variant snps.raw.vcf \
--filterExpression "QD < 2.0" \
--filterExpression "MQ < 34.0" \
--filterExpression "FS > 60.0" \
--filterExpression "MQRankSum < -12.5" \
--filterName "QDFilter" \
--filterName "MQFilter" \
--filterName "FSFilter" \
--filterName "MQRankSumFilter" \
> snpFilterInfo.txt 2>errSnpFilter

echo
echo "hard filtration INDELs"
time java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R $DB \
-o indel.filter.vcf \
--variant indel.raw.vcf \
--filterExpression "QD < 2.0"  \
--filterExpression "FS > 200.0" \
--filterName "QDFilter" \
--filterName "FSFilter"  \
> indelFilterInfo.txt 2>errIndelFilter

echo
echo "GATK eval"
time java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar \
-R $DB \
-T VariantEval -o output.eval.gatkreport \
--eval:set1 snps.filter.vcf \
--dbsnp /Users/yvans/Home/bin/GATK_resource_bundle_from_Ying_17_01_2012/1.2/b37/dbsnp_132.b37.vcf \
2>errVarinantEval > variantEval.txt

echo
echo "annovar"
time perl /Volumes/data.odin/software/mac/annovar/annovar_2012May25/convert2annovar.pl snps.filter.vcf -format vcf4 -includeinfo > snps.filter.avinput
time perl /Volumes/data.odin/software/mac/annovar/annovar_2012May25/summarize_annovar.pl snps.filter.avinput -buildver hg19 -verdbsnp 132 /Volumes/data.odin/software/mac/annovar/annovar_2012May25/humandb/ -outfile sumSNP
time perl /Volumes/data.odin/software/mac/annovar/annovar_2012May25/ indel.filter.vcf -format vcf4 -includeinfo > indel.filter.avinput
time perl /Volumes/data.odin/software/mac/annovar/annovar_2012May25/annotate_variation.pl indel.filter.avinput -buildver hg19 /Volumes/data.odin/software/mac/annovar/annovar_2012May25/humandb/ -outfile sumINDEL

date
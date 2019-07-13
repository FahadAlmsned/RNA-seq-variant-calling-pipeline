#############################################################################################
# !/bin/bash
#
# AUTHOR: Fahad Almsned
# email: almsned.fahad@gmail.com
#
# RNA-Seq variant calling pieline accoring to GATK (4.1.2.0) best practices.
# https://www.broadinstitute.org/gatk/guide/article?id=3891
#
# USAGE: bash VC_pipeline.sh <Input_Reads1.fq> <Input_Reads2.fq> <output_basename>
#
# Assumes STAR aligner and samtools are under path
#
# '--sjdbOverhang' is 'ReadLength-1'
############################################################################################

start=`date +%s`

gunzip *

pwd=$(pwd)

fwd=$pwd/$1
rev=$pwd/$2
bn=$3

#Create an output directory
opdir=$pwd/$bn"_processed"
mkdir $opdir

#Path to reference genome
fasta=/mnt/e/VariantCalling/ref_genome/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta

#1st pass
genomeDir=/mnt/e/VariantCalling/thezone/SJ/hg38
mkdir -p $genomeDir
cd $genomeDir

echo -e "["$(date)"]\t1st_Pass_Genome_Generate .."
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $fasta

runDir=/mnt/e/VariantCalling/thezone/SJ/1pass
mkdir -p $runDir
cd $runDir

echo -e "["$(date)"]\t1st_pass_Aligning .."
STAR --genomeDir $genomeDir --readFilesIn $fwd $rev --runThreadN 8

#2 pass
genomeDir=/mnt/e/VariantCalling/thezone/SJ/hg38_2pass
mkdir -p $genomeDir
cd $genomeDir

echo -e "["$(date)"]\t2nd_Pass_Genome_Generate .."
STAR --runMode genomeGenerate --genomeDir $genomeDir --runThreadN 8 --sjdbOverhang 99 --sjdbFileChrStartEnd /mnt/e/VariantCalling/thezone/SJ/1pass/SJ.out.tab --genomeFastaFiles $fasta

runDir=/mnt/e/VariantCalling/thezone/SJ/2pass
mkdir $runDir
cd $runDir

echo -e "["$(date)"]\2nd_pass_Aligning .."
STAR --genomeDir $genomeDir --readFilesIn $fwd $rev --runThreadN 8 --outFileNamePrefix $bn"_"

#----------------------------------------------
mv $bn"_Aligned.out.sam" $pwd
cd $pwd
rm -r /mnt/e/VariantCalling/thezone/SJ/
#-----------------------------------

#Picard+samtools
Picard=/mnt/e/VariantCalling/tools/picard.jar

echo -e "["$(date)"]\tPicarding_AddOrReplaceReadGroups .."
java -jar $Picard AddOrReplaceReadGroups \
I=$bn"_Aligned.out.sam" \
O=$bn".bam" \
SO=coordinate \
RGID=4 \
RGLB=lib1 \
RGPL=illumina \
RGPU=unit1 \
RGSM=20

echo -e "["$(date)"]\tPicarding_MarkDuplicates .."
java -jar $Picard MarkDuplicates \
I=$bn".bam" \
O=$bn"_dedupped.bam" \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT \
M=output.metrics

cp $fasta $pwd

echo -e "["$(date)"]\tPicarding_CreateSequenceDictionary .."
java -jar $Picard CreateSequenceDictionary \
R=resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
O=resources_broad_hg38_v0_Homo_sapiens_assembly38.dict


echo -e "["$(date)"]\tsamtooling_creating the fasta index file .."
samtools faidx resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta

#------------------------------------
#GATK
GATK=/mnt/e/VariantCalling/tools/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar
dbSNP146=/mnt/e/VariantCalling/ref_genome/dbsnp_146.hg38.vcf
millsIndels=/mnt/e/VariantCalling/ref_genome/Mills_and_1000G_gold_standard.indels.hg38.vcf

echo -e "["$(date)"]\tGATK_SplitNCigarReads .."
java -jar $GATK SplitNCigarReads \
-R resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
-I $bn"_dedupped.bam" \
-O $bn"_split.bam"

#java -jar $GATK IndexFeatureFile \
#-F $dbSNP146 $millsIndels $v0_1000G_phase1


#GATK Base recalibration (highly recommended, but does not work without known SNP data)
#Perform BQSR

echo -e "["$(date)"]\tGATK_BQSR .."
java -jar $GATK BaseRecalibrator \
-R resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
-O recal_data.table \
-I $bn"_split.bam" \
--known-sites $dbSNP146 \
--known-sites $millsIndels

echo -e "["$(date)"]\tGATK_ApplyBQSR .."
java -jar $GATK ApplyBQSR \
-R resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
-I $bn"_split.bam" \
--bqsr-recal-file recal_data.table \
-O $bn"_recal.bam"

echo -e "["$(date)"]\tGATK_HaplotypeCaller .."
java -jar $GATK HaplotypeCaller \
-R resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
-I $bn"_recal.bam" \
-stand-call-conf 20.0 \
-O $bn".vcf"

mv $bn".vcf" $opdir

#Filter variants
#echo -e "["$(date)"]\tFiltering Variants.."
#java -jar $GATK VariantFiltration \
#-R Homo_sapiens.GRCh38.dna.fasta \
#-V $bn".vcf" \
#-window 35 -cluster 3 \
#--filter-expression "QD < 2.0" --filter-name QD2 \
#--filter-expression "FS > 30.0" --filter-name FS30 \
#-O $opdir/$bn"_filtered.vcf"



rm output.metrics
rm recal_data.table
rm resources_broad_hg38_v0_Homo_sapiens_assembly38*
rm $bn"_Aligned.out.sam"
rm $bn".bam"
rm $bn"_dedupped.bam"
rm $bn"_dedupped.bai"
rm $bn"_split.bam"
rm $bn"_split.bai"
rm $bn"_recal.bam"
rm $bn"_recal.bai"
rm $bn".vcf.idx" 




echo -e "["$(date)"]\tDONE!"

end=`date +%s`

runtime=$((end-start))

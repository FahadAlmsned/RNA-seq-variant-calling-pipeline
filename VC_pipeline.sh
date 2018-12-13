#!/bin/bash
#
# AUTHOR: Fahad Almsned
#email: almsned.fahad@gmail.com
#
# RNA-Seq variant calling pieline accoring to GATK Best practices.
# https://www.broadinstitute.org/gatk/guide/article?id=3891
#
# Make sure to guzip all files using guzip *.fq.gz
#
# Call with following arguments
# bash VC_pipeline.sh <Input_Reads1.fq> <Input_Reads2.fq> <output_basename>
#
# Assumes STAR aligner is under path
#
#

pwd=$(pwd)

fwd=$pwd/$1
rev=$pwd/$2
bn=$3

#Create an output directory
opdir=$pwd/$bn"_processed"
mkdir $opdir

#Path to reference genome
fasta=/mnt/c/Users/ale6aly/Documents/VariantCalling/ref_genome/Homo_sapiens.GRCh38.dna.fasta

#1st pass
genomeDir=/mnt/c/Users/ale6aly/Documents/VariantCalling/thezone/SJ/hg38
mkdir -p $genomeDir
cd $genomeDir

echo -e "["$(date)"]\t1st_Pass_Genome_Generate .."
STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir $genomeDir \
--genomeFastaFiles $fasta

runDir=/mnt/c/Users/ale6aly/Documents/VariantCalling/thezone/SJ/1pass
mkdir -p $runDir
cd $runDir

echo -e "["$(date)"]\t1st_pass_Aligning .."
STAR --genomeDir $genomeDir \
--readFilesIn $fwd $rev \
--runThreadN 8

#2 pass
genomeDir=/mnt/c/Users/ale6aly/Documents/VariantCalling/thezone/SJ/hg38_2pass
mkdir -p $genomeDir
cd $genomeDir

echo -e "["$(date)"]\t2nd_Pass_Genome_Generate .."
STAR --runMode genomeGenerate \
--genomeDir $genomeDir \
--runThreadN 8 \
--sjdbOverhang 101 \
--sjdbFileChrStartEnd /mnt/c/Users/ale6aly/Documents/VariantCalling/thezone/SJ/1pass/SJ.out.tab \
--genomeFastaFiles $fasta

runDir=/mnt/c/Users/ale6aly/Documents/VariantCalling/thezone/SJ/2pass
mkdir $runDir
cd $runDir

echo -e "["$(date)"]\2nd_pass_Aligning .."
STAR --genomeDir $genomeDir \
--readFilesIn $fwd $rev \
--runThreadN 8 \
--outFileNamePrefix $bn"_"

#----------------------------------------------
mv $bn"_Aligned.out.sam" $pwd
cd $pwd
rm -r /mnt/c/Users/ale6aly/Documents/VariantCalling/thezone/SJ
#-----------------------------------

#Picard+samtools
Picard=/mnt/c/Users/ale6aly/Documents/VariantCalling/tools/picard.jar
samtools=/mnt/c/Users/ale6aly/Documents/VariantCalling/tools/samtools/bin/samtools

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

rm $bn"_Aligned.out.sam"

echo -e "["$(date)"]\tPicarding_MarkDuplicates .."
java -jar $Picard MarkDuplicates \
I=$bn".bam" \
O=$bn"_dedupped.bam" \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT \
M=output.metrics

rm $bn".bam"

cp $fasta $pwd

echo -e "["$(date)"]\tPicarding_CreateSequenceDictionary .."
java -jar $Picard CreateSequenceDictionary \
R=Homo_sapiens.GRCh38.dna.fasta \
O=Homo_sapiens.GRCh38.dna.dict

echo -e "["$(date)"]\tsamtooling_creating the fasta index file .."
$samtools faidx Homo_sapiens.GRCh38.dna.fasta

#------------------------------------
#GATK
GATK=/mnt/c/Users/ale6aly/Documents/VariantCalling/tools/gatk-4.0.10.1/gatk-package-4.0.10.1-local.jar
millsIndels=/mnt/c/Users/ale6aly/Documents/VariantCalling/ref_genome/Mills_and_1000G_gold_standard.indels.hg38.vcf
dbSNP146=/mnt/c/Users/ale6aly/Documents/VariantCalling/ref_genome/dbsnp_146.hg38.vcf

echo -e "["$(date)"]\tGATK_SplitNCigarReads .."
java -jar $GATK SplitNCigarReads \
-R Homo_sapiens.GRCh38.dna.fasta \
-I $bn"_dedupped.bam" \
-O $bn"_split.bam"

rm $bn"_dedupped.bam"
rm $bn"_dedupped.bai"


#Perform BQSR
echo -e "["$(date)"]\tPerforming BQSR.."
java -jar $GATK -T BaseRecalibrator \
-I $bn"_split.bam" \
-R $fasta \
-knownSites $millsIndels \
-knownSites $dbSNP138 \
-o $bn"_recal.table"


Print recalibrated reads
echo -e "["$(date)"]\tPrinting recalibrated reads.."
java -d64 -jar $gatk -T PrintReads -R $ref -I $opdir/$bn"_processed.bam" -nct 50 -BQSR $opdir/$bn"_recal.table" -o $opdir/$bn"_recal.bam" 2>$opdir/$bn.BQSR2.log

rm $opdir/$bn"_processed.bam"
rm $opdir/$bn"_processed.bai"


echo -e "["$(date)"]\tGATK_HaplotypeCaller .."
java -jar $GATK HaplotypeCaller \
-R Homo_sapiens.GRCh38.dna.fasta \
-I $bn"_split.bam" \
-stand-call-conf 20.0 \
-O $bn".vcf"

rm $bn"_split.bam"
rm $bn"_split.bai"


#Filter variants
echo -e "["$(date)"]\tFiltering Variants.."
java -jar $GATK VariantFiltration \
-R Homo_sapiens.GRCh38.dna.fasta \
-V $bn".vcf" \
-window 35 -cluster 3 \
--filter-expression "QD < 2.0" --filter-name QD2 \
--filter-expression "FS > 30.0" --filter-name FS30 \
-O $opdir/$bn"_filtered.vcf" 

echo -e "["$(date)"]\tDONE!"


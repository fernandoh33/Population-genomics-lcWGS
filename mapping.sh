#!/bin/bash
#SBATCH --account=your_account
#SBATCH --time=0-24:00
#SBATCH --ntasks=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=60

mkdir bam.files trimmed.reads out.fastqc out.fastqc.trimmed

export JAVA_TOOL_OPTIONS="-Xms256m -Xmx2g"

RAW=fastq.files
TRIMDIR=trimmed.samples
GENOME=reference/your_reference.fasta
BAMDIR=bam.files
S1=_R1.fastq.gz
S2=_R2.fastq.gz

module load StdEnv/2023
module load gcc/12.3
module load fastp/0.23.4
module load samtools/1.20
module load fastqc/0.12.1
module load bamtools/2.5.2
module load bwa/0.7.18
module load qualimap/2.3

#Index the Genome
bwa index $GENOME

#Create Sequence List
#ls $RAW | grep $S1 | sed 's/_R1.fastq.gz//' | sort -u > list.samples.txt

#Remove adapters, duplicates, poly Gs, and low quality bases
for i in $(cat list.samples.txt);
do fastp -i $RAW/$i$S1 -I $RAW/$i$S2 -o $TRIMDIR/$i$S1 -O $TRIMDIR/$i$S2 --cut_right --dedup -h $TRIMDIR/$i'.html' -g -w 60;
done

#Run FastQC on trimmed reads
for i in $(cat list.samples.txt);
do fastqc $TRIMDIR/$i$S1 -o out.fastqc.trimmed;
fastqc $TRIMDIR/$i$S2 -o out.fastqc.trimmed;
done

#Run Alignment, index, and bam statistics
for i in $(cat list.samples.txt);
do bwa mem -t 60 $GENOME $TRIMDIR/$i$S1 $TRIMDIR/$i$S2|samtools view -bh|samtools sort -T tmp -@ 60 -o $BAMDIR/$i'.bam';
samtools index -@ 60 $BAMDIR/$i'.bam';
bamtools stats -in $BAMDIR/$i'.bam' > $BAMDIR/$i'_bamstats.txt';
qualimap bamqc -bam $BAMDIR/${i}'.bam' -outdir $BAMDIR/${i}_qualimap -outfile ${i}_qualimap_report.txt -nt 60;
done

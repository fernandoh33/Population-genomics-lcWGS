#!/bin/bash
#SBATCH --account=your_account
#SBATCH --time=0-24:00
#SBATCH --ntasks=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=24

#ngsRelate website: https://github.com/ANGSD/NgsRelate
#installing ngsRelate
git clone --recursive https://github.com/SAMtools/htslib
git clone https://github.com/ANGSD/ngsRelate
cd htslib/;make -j2;cd ../ngsRelate;make HTSSRC=../htslib/;
#First, generate a file with allele frequencies (geno.lowld.mafs.gz) and a file with genotype likelihoods (geno.lowld.beagle.gz)
#for ngsRelate, snps are assumed to be independent, so we should only include sites in low ld (-rg)
#unlinked positions is a file with one site per row as: Chr1:1234
#adjust -nInd, -minInd, and REF variable
REF=path_to_your_reference/reference.fasta

angsd -bam bam.list -fai $REF.fai -nInd 20 -doMajorMinor 1 -doPost 1 -doMaf 1 -doGlf 2 -out geno.lowld -gl 2 -minMapQ 30 -minQ 20 -minMaf 0.05 -SNP_pval 1e-6 -minInd 10 -rf unlinked.positions
#Then we extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)
zcat geno.lowld.mafs.gz | cut -f5 |sed 1d > allele.freq
#run NgsRelate
#adjust -n (number of individuals) and full_path/
full_path/ngsRelate  -G geno.lowld.beagle.gz -n 20 -f allele.freq  -O out.relatedness -z samples.list

#ngsRelate website: https://github.com/ANGSD/NgsRelate
#installing ngsRelate
git clone --recursive https://github.com/SAMtools/htslib
git clone https://github.com/ANGSD/ngsRelate
cd htslib/;make -j2;cd ../ngsRelate;make HTSSRC=../htslib/;

#First, we need to generate a file with allele frequencies (geno.file.mafs.gz) and a file with genotype likelihoods (geno.file.beagle.gz)
#adjust -nInd and -minInd

angsd -bam bam.list -fai $REF.fai -nInd 20 -doMajorMinor 1 -doPost 1 -doMaf 1 -doGlf 2 -out geno.file -gl 2 -minMapQ 30 -minQ 20 -minMaf 0.05 -SNP_pval 1e-6 -minInd 10

### Then we extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)

zcat geno.file.mafs.gz | cut -f5 |sed 1d > allele.freq

### run NgsRelate
#adjust -n (number of individuals)

ngsrelate  -G geno.file.beagle.gz -n 20 -f allele.freq  -O out.relatedness

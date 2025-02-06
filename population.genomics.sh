#!/bin/bash
#SBATCH --account=your_account
#SBATCH --time=0-24:00
#SBATCH --ntasks=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=24

#Population genomics analyses using low coverage WGS
#First, we need to install three programs: angsd for most of the analyses, ngsLD and prune_graph for ld pruning
#installing angsd
git clone https://github.com/ANGSD/angsd.git 
cd angsd
make
#installing ngsLD
#ngsLD website: https://github.com/EFox16/ngsLD-Tutorial
git clone https://github.com/fgvieira/ngsLD.git
cd ngsLD
make
make test
make clean
#installing prune_graph
git clone https://github.com/fgvieira/prune_graph.git
cd prune_graph
cargo build --release
cargo test
#Once the three programs are installed and working, we can start the analyses 
REF=path_to_your_reference/reference.fasta
#generate a file with genotype likelihoods of high quality sites (variant + invariant) for estimating diversity statistics
#if no ancestral reference is available, the sfs is folded, and the reference is used as ancestral, -ref need also be provided for adjusting base qualities (-baq 1)
#adjust -nInd -minInd
angsd -bam bam.list -nInd 20 -minInd 10 -doSaf 1 -baq 1 -anc $REF -ref $REF -GL 2 -P 16 -minMapQ 30 -minQ 20 -out all.sites
#finding the global estimate
realSFS out.saf.idx -P 16 -fold 1 > out.sfs
#calculate theta per site
realSFS saf2theta out.saf.idx -sfs out.sfs -outname out
#calculate theta per window, the example is window size 100kb with steps of 20kb but can be changed of course
thetaStat do_stat out.thetas.idx -win 100000 -step 20000  -outnames theta.w100.s20
#generate a file with genotype likelihoods of high quality snps with maf > 0.05 as invariant sites are useless for ld estimation
#adjust -nInd and -minInd 
angsd -bam bam.list -fai $REF.fai -nInd 20 -doMajorMinor 1 -doPost 1 -doMaf 1 -doGlf 2 -out geno.file -gl 2 -minMapQ 30 -minQ 20 -minMaf 0.05 -SNP_pval 1e-6 -minInd 10
#make a list of positions from the .mafs file
zcat geno.file.mafs.gz | cut -f 1,2 | tail -n +2 > geno.file_positions.txt
#get positions in low LD
#adjust full_path, --nInd, and --n_sites
full_path/ngsLD --n_threads 30 --verbose 1 --n_ind 20 --n_sites 10000 --geno geno.file.beagle.gz --probs --pos geno.file_positions.txt --max_kb_dist 50 --min_maf 0.05 --extend_out| sort -k 1,1Vr -k 2,2V > testLD.ld
full_path/prune_graph --header --in testLD.ld --weight-field "r2" --weight-filter "dist <= 50000 && r2 >= 0.2" --out unlinked.positions
#generate the beagle file with snps filtered for low LD
angsd -b bam.list -rf unlinked.positions -out geno.lowld -doMajorMinor 1 -doPost 1 -doMaf 1 -doGlf 2 -gl 2 -fai $REF.fai
#run pca using pcangsd
pcangsd -b geno.lowld.beagle.gz -t 16 -o out.pcangsd
#run ngsadmix, ngsadmix is within angsd, so no need to install
#the output files with admixture proportions are out.admixture.K*.qopt
for i in `seq 1 10`;do NGSadmix -likes geno.lowld.beagle.gz -K $i -outfiles out.admixure.K$i -P 10;done
#in R
samples.info=read.table("sample.info.txt",header=F) #sample.info.txt is a file with the name of the samples in a single column, with no header
out.pca=as.matrix(read.table("~/OneDrive/Downloads/out.pcangsd.cov"))
eigenvalues=eigen(out.pca)
plot(eigenvalues$vectors)
out.pca=cbind(samples.info, eigenvalues$vectors[,1:10)
colnames(out.pca)=c("sample","PC1",PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

#!/bin/bash
#SBATCH --account=your-account
#SBATCH --time=0-12:00
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --cpus-per-task=8

module load StdEnv/2023
module load angsd/0.940

#get headers from beagle and mafs files, headers are the same for all chromosomes so any chromosome file is useful
zcat geno.for.ld/Scaffold_10.vimineum.for.ld.maf0.1.minInd0.5.beagle.gz | head -1 > header.beagle
zcat geno.for.ld/Scaffold_10.vimineum.for.ld.maf0.1.minInd0.5.mafs.gz | head -1 > header.mafs

# stiltgrass.chrs is a list of chromosomes, one per line, the format of ngsLD is chr_position and the beagle format is chr:position so we need to convert
for chr in $(cat stiltgrass.chrs);
do sed 's/:/_/g' ld.estimation/$chr.unlinked.positions > tmp.$chr.unlinked.positions;
zcat geno.for.ld/$chr.vimineum.for.ld.maf0.1.minInd0.5.beagle.gz | awk 'NR==FNR{a[$1]; next} $1 in a' tmp.$chr.unlinked.positions - > tmp.$chr.filtered.txt;
cat header.beagle tmp.$chr.filtered.txt |bgzip -c > ld.estimation/$chr.vimineum.lowLD.maf.0.1.beagle.gz;
done;

# making a unique beagle file for all chrs is now doable as we have a few thousands of positions instead of millons
cat header.beagle tmp*.filtered.txt | bgzip -c > ld.estimation/all.chrs.vimineum.lowLD.maf.0.1.beagle.gz
rm tmp*

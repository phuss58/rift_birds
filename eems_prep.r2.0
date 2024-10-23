library(vcfR)

options(scipen=999)

filepath <- "/common/cooperlab/phuss58/albertine/06_plink/"
savepath <- "/common/cooperlab/phuss58/albertine/06_eems/"

source(paste0(savepath,"bed2diffs.r"))

## Batis

vcf <- read.vcfR(paste0(filepath,"batis_nothin.recode.vcf"))

genotypes <- vcfR2genlight(vcf) 

genotype_matrix <-as.matrix(genotypes)

diff_matrix <- bed2diffs_v2(genotypes = genotype_matrix)

write.table(diff_matrix, paste0(savepath,"batis_nothin.diffs"), col.names = F, row.names = F, quote = F, sep='\t')

## Phylloscopus

vcf <- read.vcfR("phylloscopus_nothin.recode.vcf")

genotypes <- vcfR2genlight(vcf) 

genotype_matrix <-as.matrix(genotypes)

diff_matrix <- bed2diffs_v2(genotypes = genotype_matrix)

write.table(diff_matrix, "phyllscopus_nothin.diffs", col.names = F, row.names = F, quote = F, sep='\t')

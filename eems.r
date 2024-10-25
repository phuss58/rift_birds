library(vcfR)

options(scipen=999)

filepath <- "/common/cooperlab/phuss58/albertine/05_filtered_vcf/"
savepath <- "/common/cooperlab/phuss58/albertine/06_eems/"

source(paste0(savepath,"bed2diffs.r"))

## Batis

vcf <- read.vcfR(paste0(filepath,"batis_nothin.vcf"))

genotypes <- vcfR2genlight(vcf) 

genotype_matrix <-as.matrix(genotypes)

diff_matrix <- bed2diffs_v2(genotypes = genotype_matrix)

write.table(diff_matrix, paste0(savepath,"batis_nothin.diffs"), col.names = F, row.names = F, quote = F, sep='\t')

## Phylloscopus

vcf <- read.vcfR("phylloscopus_nothin.vcf")

genotypes <- vcfR2genlight(vcf) 

genotype_matrix <-as.matrix(genotypes)

diff_matrix <- bed2diffs_v2(genotypes = genotype_matrix)

write.table(diff_matrix, paste0(savepath, "phylloscopus_nothin.diffs"), col.names = F, row.names = F, quote = F, sep='\t')


## Chamaetylas

vcf <- read.vcfR(paste0(filepath,"chamaetylas_nothin.vcf"))

genotypes <- vcfR2genlight(vcf) 

genotype_matrix <-as.matrix(genotypes)

diff_matrix <- bed2diffs_v2(genotypes = genotype_matrix)

write.table(diff_matrix, paste0(savepath,"chamaetylas_nothin.diffs"), col.names = F, row.names = F, quote = F, sep='\t')


## Cossypha
vcf <- read.vcfR(paste0(filepath,"cossypha_nothin.vcf"))

genotypes <- vcfR2genlight(vcf) 

genotype_matrix <-as.matrix(genotypes)

diff_matrix <- bed2diffs_v2(genotypes = genotype_matrix)

write.table(diff_matrix, paste0(savepath,"cossypha_nothin.diffs"), col.names = F, row.names = F, quote = F, sep='\t')

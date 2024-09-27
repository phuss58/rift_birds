# load bcftools
module load bcftools

# ex

vcftools --vcf ${workdir}/04_vcf/CM001988.1.vcf --max-missing 1.0 --min-alleles 2 --max-alleles 2 --max-maf 0.49 --mac 2 --recode --recode-INFO-all --out ${workdir}/test.file

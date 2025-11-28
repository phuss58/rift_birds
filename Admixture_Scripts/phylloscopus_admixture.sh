## example run with phylloscopus

# combine all chromosome vcfs into one vcf
# pick one file, use headers

workdir=/common/cooperlab/phuss58/albertine/
admixture=/common/cooperlab/phuss58/albertine/06_plink/dist/admixture_linux-1.3.0/admixture

module load plink/1.90
module load vcftools/0.1

cd ${workdir}/05_filtered_vcf

grep "#" phylloscopus_structure_CM001988.1_structure_nowindow.recode.vcf > phylloscopus_nothin.vcf
for i in $( ls phylloscopus_structure_* ); # for i in all files
    do grep -v "#" $i >> phylloscopus_nothin.vcf; # get the lines, combine in file, skip headers
    done

# did it work?
wc -l phylloscopus_nothin.vcf

# make chromosome map for the vcfs
grep -v "#" phylloscopus_nothin.vcf | # get the file headers
    cut -f 1 | # cut the headers
    uniq | # get unique values
    awk '{print $0"\t"$0}' > phylloscopus_chrom_map.txt # copy and place in new file

# run vcftools for the combined vcfs
# 5kbp thinning
# vcftools --vcf 5kbpthin.vcf  --plink --chrom-map chrom_map.txt --out 5kbpthin 
# no thinning
vcftools --vcf phylloscopus_nothin.vcf  --plink --chrom-map phylloscopus_chrom_map.txt --out ${workdir}/06_plink/phylloscopus_nothin 


# convert  with plink
# 5kbp thinning
# plink --file ${workdir}/05_structure/5kbpthin --recode12 --allow-extra-chr --out ${workdir}/05_structure/5kbpthin_plink
# no thinning
plink --file ${workdir}/06_plink/phylloscopus_nothin --recode12 --allow-extra-chr --out ${workdir}/06_plink/phylloscopus_nothin_plink

# run pca on each dataset
# 5kbp thinning
# plink --file ${workdir}/05_structure/5kbpthin_plink --pca --allow-extra-chr --out ${workdir}/05_structure/5kbpthin_plink_pca
# no thinning
plink --file ${workdir}/06_plink/phylloscopus_nothin_plink --pca --allow-extra-chr --out ${workdir}/06_plink/phylloscopus_nothin_plink_pca


# run admixture on each dataset
# 5kbp thinning
# for K in 3; do admixture --cv ${workdir}/05_structure/5kbpthin_plink.ped $K; done
# no thinning
for K in 2 3 4 5; do $admixture --cv ${workdir}/06_plink/phylloscopus_nothin_plink.ped $K | tee ${workdir}/log${K}.out; done
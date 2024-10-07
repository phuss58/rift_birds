# running admixture

for K in 3; do admixture --cv ${workdir}/05_structure/5kbpthin_plink.ped $K; done
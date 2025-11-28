admixture=/common/cooperlab/phuss58/albertine/06_plink/dist/admixture_linux-1.3.0/admixture



# running admixture

for K in 1 2 3 4 5 6; do $admixture --cv ${workdir}/05_structure/5kbpthin_plink.ped $K; done

for K in 1 2 3 4 5 5 6; do $admixture --cv ${workdir}/06_plink/nothin_plink.ped $K | tee log${K}.out; done
# Gather gene fpkm from cufflinks 
cd ${cufflinks_dir}

gene_list=(...)

for gene in ${gene_list[*]}
do 
grep -w ${gene} *WT*/genes.fpkm_tracking > gene_2017_subset/WT_${gene}_fpkm.txt
grep -w ${gene} *KO*/genes.fpkm_tracking > gene_2017_subset/KO_${gene}_fpkm.txt
done
#verify using: ls gene_2017_subset | wc -l

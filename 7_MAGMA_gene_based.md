```bash

## MAGMA gene-based analysis on 692 LOF variants and INT 7 Phylum, 63 Genus, 120 species and 181 pathways. Sample shown for MAGMA analysis from one Bacterial species 

./magma --annotate --snp-loc W_SNPeff_LOF_05.bim --gene-loc NCBI37.3.gene.loc --out W_SNPeff_LOF_05

./magma --bfile W_SNPeff_LOF_05 --gene-annot W_SNPeff_LOF_05.genes.annot --pheno file=meta_pheno_magma.txt use=m1  --covar file=COV_MAGMA.txt include=Group-Batch  --out W_MAGMA_LOF05_m1


```bash

## MAGMA gene-based analysis on high impact variants 

./magma --annotate --snp-loc W_clean_05.bim --gene-loc NCBI37.3.gene.loc --out W_ALL_05

#!/bin/bash

# Extract phenotype names (excluding FID and IID)
phenotypes=$(head -n1 Pheno_All_forGWAS_noFID.txt | cut -f3-)

# Convert phenotype names into an array
IFS=$'\t' read -r -a pheno_array <<< "$phenotypes"

# Loop over each phenotype
for phenotype in "${pheno_array[@]}"; do
  echo "Analyzing phenotype: $phenotype"

  ./magma \
    --bfile W_SNPeff_LOF_05 \
    --gene-annot W_SNPeff_LOF_05.genes.annot \
    --pheno file=Pheno_All_forGWAS_noFID.txt use="$phenotype" \
    --covar file=META_COV_noFID.txt include=Age-PC10 \
    --out W_LOF05_"$phenotype"
done





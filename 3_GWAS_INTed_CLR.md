```bash

## GWAS on 7 Phylum and 63 Genus from 16S INTed data using Plink2

plink2 --bfile W_16S_rmrelat --pheno 16S_INT_final.txt --covar COV_16S_FINAL.txt --covar-name Sex Age Batch Pheno PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 --covar-variance-standardize --glm hide-covar cols=+a1freq,+machr2 --out INT

sort -k 14 -g INT.Genus41.glm.linear > Genus41_INT.linear

awk '$14 != "NA"' FS=' ' Genus41_INT.linear > Genus41_INTf.linear

awk '{if ($7 >=0.01) print $0}' Genus41_INTf.linear > Genus41_INT01.linear

awk '{if ($7 >=0.05) print $0}' Genus41_INTf.linear > Genus41_INT05.linear

## GWAS on 120 Species INT using Plink2

plink2 --bfile W_meta_rmrelat --pheno mOTU_INT_final.txt --covar COV_META_FINAL.txt --covar-name Sex Age Batch Pheno PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 --covar-variance-standardize --glm hide-covar cols=+a1freq,+machr2 --out INT

sort -k 14 -g INT.m1.glm.linear > m1_INT.linear

awk '$14 != "NA"' FS=' ' m1_INT.linear > m1_INTf.linear

awk '{if ($7 >=0.05) print $0}' m1_INTf.linear > m1_INT05.linear

## GWAS on 181 Pathways INT using Plink2

plink2 --bfile W_meta_rmrelat --pheno path_INT_final.txt --covar COV_META_FINAL.txt --covar-name Sex Age Batch Pheno PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 --covar-variance-standardize --glm hide-covar cols=+a1freq,+machr2 --out INT

sort -k 14 -g INT.K1.glm.linear > K1_INT.linear

awk '$14 != "NA"' FS=' ' K1_INT.linear > K1_INTf.linear

awk '{if ($7 >=0.05) print $0}' K1_INTf.linear > K1_INT05.linear


## GWAS on 7 Phulum and 63 CLR transformed using Plink2*

plink2 --bfile W_16S_rmrelat --pheno pheno_16S.txt --covar covar_16S_FINAL.txt --covar-name Sex Age Batch Pheno PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 --covar-variance-standardize --glm hide-covar cols=+a1freq,+machr2 --out CLR

sort -k 14 -g CLR.Phy7.glm.linear > Phy7_CLR.linear

awk '$14 != "NA"' FS=' ' Phy7_CLR.linear > Phy7_CLRf.linear

awk '{if ($7 >=0.05) print $0}' Phy7_CLRf.linear > Phy7_CLR05.linear

## GWAS on 120 Species CLR transformed using Plink2

plink2 --bfile W_meta_rmrelat --pheno mOTU_CLR_FINAL_JC.txt --covar META_COV.txt --covar-name Sex Age Batch Pheno PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 --covar-variance-standardize --glm hide-covar cols=+a1freq,+machr2 --out CLR

sort -k 14 -g CLR.m1.glm.linear > m1_CLR.linear

awk '$14 != "NA"' FS=' ' m1_CLR.linear > m1_CLRf.linear

awk '{if ($7 >=0.05) print $0}' m1_CLRf.linear > m1_CLR05.linear


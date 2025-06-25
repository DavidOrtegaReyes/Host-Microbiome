# Interplay between host genetics and gut microbiome composition in the Japanese population

This repository contains the codes used for a study on the relationship between host genetics and gut microbiota composition in the Japanese population. The study aimed to investigate how host genetic variation influences the gut microbiota relative abundance composition

A. Contents of the folder (list of scripts) in order of analysis: 
1_GWAS_QC.md
2_Microbiome_relative_abundance_filtering_transformation.md
2.5_Shapiro_Wilk.md
3_GWAS_Host_Microbiome.md
3.1_GWAS_Binary.md
3.2_GWAS_Quantitative.md
4_variant_annotation.md
5_PheWAS_INTed.md
6_MAGMA_gene_based.md
7_replication_study_Ishida_method.md   

B. Description of the scripts:
1. How quality control was conducted for GWAS 
2. Filtering, normalization and transformation of relative abundance from shotgun metagenomic sequencing mOTUs
2.1. Shapiro Wilk test for identification of binary phenotypes
3. GWAS between host genotype and relative abundance composition using PLINK2
3.1 Script for binary GWAS
3.2 Script for quantitative GWAS 
4. Annotation of functional varints on vcf files using SNPeff
5. PheWAS between loss of function variants and all filtered taxa
6. Gene-based analysis between loss of function variants and all filtered taxa using MAGMA
7. Ishida et al. 16S rRNA method replication using QIIME2

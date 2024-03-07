```bash

## Annotation of genotype: SNPeff

#Variant annotation: used Qced Plink file
#Convert to VCF:

plink --bfile W_clean_Auto --recode vcf --make-bed --out W_clean

bgzip W_clean.vcf

tabix -p vcf  W_clean.vcf.gz

# Select and filter variants:

/home/david/jdk-16.0.2/bin/java -Djava.io.tmpdir=/home/david/WGS_IGT/tmp -jar ~/gatk-4.2.5.0/gatk-package-4.2.5.0-local.jar SelectVariants --max-indel-size 10 -V W_clean.vcf.gz -O W_annot_selected.vcf

cat W_annot_selected.vcf | /home/david/jdk-16.0.2/bin/java -Djava.io.tmpdir=/home/david/WGS_IGT/tmp -jar SnpSift.jar filter " (( exists INDEL ) & ( QUAL >= 20 )) |  ( QUAL >= 30 ) " > W_annot.filtered.vcf

# Annotate variants:

java -Djava.io.tmpdir=/Volumes/HD-PCFSU3-A/tmp -jar snpEff.jar -v -lof GRCh37.75 W_annot_selected.vcf > W_annotated_SNPeff_lof5.vcf

# Extract high and moderate effect variants: 

java -Djava.io.tmpdir=/Volumes/HD-PCFSU3-A/tmp -jar SnpSift.jar filter "(ANN[*].IMPACT has 'HIGH') | (ANN[*].IMPACT has 'MODERATE')" W_annotated_SNPeff_lof5.vcf  > W_annotated_SNPeff_HIGH_MODERATE5.vcf

# Extract missense and nonsense variants:** 

java -Djava.io.tmpdir=/Volumes/HD-PCFSU3-A/tmp -jar SnpSift.jar filter "(EFF[*].FUNCLASS has 'MISSENSE') | (EFF[*].FUNCLASS has 'NONSENSE')" W_annotated_SNPeff_lof5.vcf  > W_annotated_SNPeff_FUNCLASS5.vcf

# Extract LOF variants:

java -Djava.io.tmpdir=/Volumes/HD-PCFSU3-A/tmp -jar SnpSift.jar filter "(exists LOF[*].PERC)" W_annotated_SNPeff_lof5.vcf > W_annotated_SNPeff_LOF_5.vcf

# Change LOF-extracted VCFs to Plink format with different MAFs for further analyses:

plink --vcf W_annotated_SNPeff_LOF_5.vcf.gz --recode --make-bed --out W_SNPeff_LOF

plink --bfile W_SNPeff_LOF --maf 0.01 --make-bed --out W_SNPeff_LOF_01

plink --recodeA --bfile W_SNPeff_LOF_01 --extract SNPs_LOF_01.txt --out W_SNPeff_LOF_01_SNPs

plink --recodeA --bfile W_SNPeff_LOF_05 --out W_SNPeff_LOF_05_SNPs
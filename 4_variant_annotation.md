```bash

## Annotation of genotype: SNPeff

#Variant annotation: used Qced Plink file
#Convert to VCF:

plink --bfile W_clean_Auto --recode vcf --make-bed --out W_clean

bgzip W_clean.vcf

tabix -p vcf  W_clean.vcf.gz

# Select and filter variants:

/home/david/jdk-16.0.2/bin/java -Djava.io.tmpdir=tmp -jar ~/gatk-4.2.5.0/gatk-package-4.2.5.0-local.jar SelectVariants --max-indel-size 10 -V W_clean.vcf.gz -O W_annot_selected.vcf

cat W_annot_selected.vcf | java -Djava.io.tmpdir=tmp -jar SnpSift.jar filter " (( exists INDEL ) & ( QUAL >= 20 )) |  ( QUAL >= 30 ) " > W_annot.filtered.vcf

# Annotate variants:

java -Djava.io.tmpdir=tmp -jar snpEff.jar -v -lof GRCh37.75 W_annot_selected.vcf > W_annotated_SNPeff_lof5.vcf

# Extract high and moderate effect variants: 

java -Djava.io.tmpdir=/tmp -jar SnpSift.jar filter "(ANN[*].IMPACT has 'HIGH')" W_annotated_SNPeff_lof5.vcf  > W_annotated_SNPeff_HIGH.vcf


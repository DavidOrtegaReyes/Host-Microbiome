```bash

## Pre QC

# Filter out low quality reads:

for i in /home/david/WGS_IGT/W35PASS/*.vcf.gz ; do bcftools view -f PASS $i -Oz -o $i.out; done

for i in /home/david/WGS_IGT/W36PASS/*.vcf.gz ; do bcftools view -f PASS $i -Oz -o $i.out; done

for file in *.hard-filtered.vcf.gz.out
do
mv $file ${file%-}.vcf.gz
done

for i in *hard.vcf.gz ; do tabix -p vcf $i; done

# Concatenate individual VCFs:

bcftools merge *hard.vcf.gz -Oz -o W36_merged.vcf.gz

Tools/bcftools merge *hard.vcf.gz -Oz -o W35_merged.vcf.gz

tabix -p vcf  W35_merged.vcf.gz

tabix -p vcf  W36_merged.vcf.gz

# Merge W35 and W36:

bcftools merge W35_merged_FINAL.vcf.gz W36_merged_FINAL.vcf.gz --missing-to-ref -Oz -o W_merged_FINAL.vcf.gz

#Remove duplicated individuals:

vcftools --remove dupsamples.txt --gzvcf W_merged_FINAL.vcf.gz --recode --stdout | gzip -c > W_dedup.vcf.gz

tabix -p vcf W_dedup.vcf.gz

# Lift down hg38 to hg19:
jdk-16.0.2/bin/java -Djava.io.tmpdir=/home/david/WGS_IGT/tmp -jar picard.jar LiftoverVcf I=W_dedup.vcf.gz O=W_merged_hg19.vcf.gz CHAIN=hg38ToHg19.over.chain REJECT=rejected_variants.vcf R=ucsc.hg19.fasta TMP_DIR=/home/david/WGS_IGT/tmp MAX_RECORDS_IN_RAM=10000

# Rename chromosomes:
bcftools annotate --rename-chrs /home/david/DKDanalysis/ChrX/chr_names.txt W_merged_hg19.vcf.gz  -Oz -o W_hg19_chr.vcf.gz

tabix -p vcf W_hg19_chr.vcf.gz

# Include SNPs for dbSNP:
bcftools annotate -c ID -a /home/gluster/CommonData/PublicData/ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b149_GRCh37p13/VCF/All_20161121.vcf.gz W_hg19_chr.vcf.gz > W_merged_FINAL.vcf.gz

tabix -p vcf W_merged_FINAL.vcf.gz

# Recode to PLINK:
plink --vcf W_merged_FINAL.vcf.gz --recode --allow-extra-chr 0 --make-bed --out W_merged_FINAL

plink --bfile W_merged_FINAL --split-x b37 --make-bed --out W_split_FINAL

plink --bfile W_split_FINAL --impute-sex .3 .3 --make-bed --out W_split

# Change to chr:pos
awk '{if($2 == "."){print $1,$1":"$4,$3,$4,$5,$6} else {print $0}}' W_split_is.bim > W_split_is_chrpos.bim

mv W_split_is_chrpos.bim W_split_is.bim


## Per-sample QC

# Summary of data:

plink --bfile W_split --make-bed --allow-extra-chr

# Gender check:

plink --bfile W_split_is --set-hh-missing --make-bed --out W_merged_cleanerx

plink --bfile W_merged_cleanerx --check-sex

grep "PROBLEM" plink.sexcheck | awk '{print$1,$2}'> sex_discrepancy.txt

plink --bfile W_merged_cleanerx --remove sex_discrepancy.txt --make-bed --out W_merged_sexcheck

# Identify samples (individuals) of poor quality:
#(i) call rate (the proportion of SNPs at which a genotype has been called)
#(ii) heterozygosity (i.e. the proportion of SNPs at which the called genotype is a heterozygote).

plink --bfile W_merged_cleanerx --missing --out plink.missing

plink --bfile W_merged_cleanerx --exclude inversion.txt --range --indep-pairwise 50 5 0.2 --out W_indepSNP

plink --bfile W_merged_cleanerx --extract W_indepSNP.prune.in --het --out R_check

# Remove individuals that do not meet sample call rate QC criteria (<95%):

R --vanilla --slave --args plink.missing.imiss R_check.het 0.95 0 0.4 < sampleqc.remove.R > remove.sampleqc.txt

head remove.sampleqc.txt 

plink --bfile W_merged_cleanerx --remove remove.sampleqc.txt --make-bed --out W_merged_rmcr

# Remove individuals that do not meet sample heterozygosity QC criteria (4 SD from the mean):
#The following code generates a list of individuals who deviate more than 4 standard deviations from the heterozygosity rate mean:

Rscript --no-save heterozygosity_outliers_list.R

sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt

plink --bfile W_merged_cleanerx --remove het_fail_ind.txt --make-bed --out W_rmhet

# Evaluate population structure
# Identify individuals of divergent ancestry by combining the genotypes of the study population with genotypes of a reference dataset consisting of individuals from known ethnicities (1000 genomes phase 3). Principal component analysis (PCA) is then used to detect population structure.

#Match study genotypes and reference data
#Filter reference and study data for non A-T or G-C SNPs

awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" || $5$6 == "AT" || $5$6 == "TA")  {print $2}' W_rmhet.bim  > W_nf.ac_gt_snps

awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" || $5$6 == "AT" || $5$6 == "TA")  {print $2}'  hapmap3.bim  > all_hapmap.ac_gt_snps

plink --bfile  hapmap3 --exclude all_hapmap.ac_gt_snps --allow-extra-chr --make-bed --out all_hapmap.no_ac_gt_snps

plink --bfile  W_rmhet --exclude W_nf.ac_gt_snps --make-bed --out W_nf.no_ac_gt_snps

# Prune study data
# Prune for variants in LD with an r2> 0.2 in a 50kb window, using plink --indep-pairwise to compute the LD-variants; exclude range is used to remove genomic ranges of known high-LD structure (highld.txt)

plink --bfile  W_nf.no_ac_gt_snps --exclude range  highld.txt --indep-pairwise 500 5 0.05 --out W_nf.no_ac_gt_snps

plink --bfile  W_nf.no_ac_gt_snps --extract W_nf.no_ac_gt_snps.prune.in --make-bed --out W_nf.pruned

# Filter reference data for the same SNP set as in study

plink --bfile  all_hapmap.no_ac_gt_snps --extract W_nf.no_ac_gt_snps.prune.in --allow-extra-chr --make-bed --out all_hapmap.pruned

# Check and correct chromosome mismatch

awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} ($2 in a && a[$2] != $1) {print a[$2],$2}' W_nf.pruned.bim all_hapmap.pruned.bim | sed -n '/^[XY]/!p' > all_hapmap.toUpdateChr

plink --bfile all_hapmap.pruned --update-chr all_hapmap.toUpdateChr 1 2 --make-bed --out all_hapmap.updateChr

# Position mismatch

awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} ($2 in a && a[$2] != $4)  {print a[$2],$2}' W_nf.pruned.bim all_hapmap.pruned.bim > all_hapmap.toUpdatePos

# Possible allele flips

Check if non-matching allele codes are a case of allele flips.

awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' W_nf.pruned.bim all_hapmap.pruned.bim > all_hapmap.toFlip

# Upate positions and flip alleles

plink --bfile all_hapmap.updateChr --update-map all_hapmap.toUpdatePos 1 2 --flip all_hapmap.toFlip --make-bed --out all_hapmap.flipped

# Remove mismatches

awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' W_nf.pruned.bim all_hapmap.flipped.bim > all_hapmap.mismatch

plink --bfile all_hapmap.flipped --exclude all_hapmap.mismatch --make-bed --out all_hapmap.clean

# Merge study genotypes and reference data*

plink --bfile W_nf.pruned  --bmerge all_hapmap.clean.bed all_hapmap.clean.bim all_hapmap.clean.fam  --make-bed --out W_nf.merge.all_hapmap

#PCA analysis

plink --bfile W_nf.merge.all_hapmap --cluster --pca 10 --maf 0.05 --out MDS_Wmb_hapmap

# Convert population codes into superpopulation codes

$awk '{print$2,$2,$7}' relationships_w_pops_121708.txt > race_1kG.txt

sed 's/CHB/ASN/g' race_1kG.txt>race_1kG2.txt

sed 's/ASW/AFR/g' race_1kG2.txt>race_1kG3.txt

sed 's/CEU/EUR/g' race_1kG3.txt>race_1kG4.txt

sed 's/GIH/AMR/g' race_1kG4.txt>race_1kG5.txt

sed 's/CHD/ASN/g' race_1kG5.txt>race_1kG6.txt

sed 's/MKK/AFR/g' race_1kG6.txt>race_1kG7.txt

sed 's/LWK/AFR/g' race_1kG7.txt>race_1kG8.txt

sed 's/TSI/EUR/g' race_1kG8.txt>race_1kG9.txt

sed 's/MEX/AMR/g' race_1kG9.txt>race_1kG10.txt

sed 's/YRI/AFR/g' race_1kG10.txt>race_1kG11.txt

sed 's/FIN/EUR/g' race_1kG11.txt>race_1kG12.txt

sed 's/JPT/JPT/g' race_1kG12.txt>race_1kG13.txt

sed 's/PUR/AMR/g' race_1kG13.txt>race_1kG14.txt

# Create a racefile of study data

awk '{print$1,$2,"OWN"}' W_nf.pruned.fam > racefile_own.txt

# Concatenate racefiles

cat race_1kG14.txt racefile_own.txt | sed -e '1i\FID IID race' > racefile.txt

# Exclude ethnic outliers
# Select individuals in data below cut-off thresholds. 

awk '{ if ($3 >0.012 && $4 >-0.005) print $1,$2 }' MDS_Wmb_hapmap.eigenvec > MDS_Wmb_hapmap_popclean.txt

plink --bfile  W_rmhet --remove-fam  MDS_Wmb_hapmap_popclean.txt --make-bed --out W_popclean

# Assess genetic relatedness

plink --bfile W_clean_allsamples --maf 0.05 --geno 0.05 --hwe 0.000001 --indep-pairwise 50 5 0.2

# Calculate relatedness between each pair of samples:

plink --bfile W_clean_allsamples --extract plink.prune.in --genome --min 0.25

plink --bfile W_clean_allsamples --missing

# Generate and remove a list containing the sample from each pair with the lowest call rate

R --vanilla --slave --args plink.imiss plink.genome < removeRelated.R > remove.related.txt

plink --bfile W_meta_clean_Auto --remove remove.related.txt --make-bed --out W_rmrelat_allsamples


## Per-marker QC

plink --bfile W_16S_sampleclean --chr X --make-bed --out W_16S_sampleclean_ChrX

plink --bfile W_16S_sampleclean --autosome --make-bed --out W_16S_sampleclean_Auto

plink --bfile W_meta_sampleclean_Auto --geno 0.05 --hwe 0.000001 --make-bed --out W_meta_clean_Auto


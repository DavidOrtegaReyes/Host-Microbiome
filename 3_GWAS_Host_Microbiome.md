```bash


awk '{OFS="\t"; $2 = $1":"$4":"$5":"$6; print}' W_meta_clean_Auto.bim > W_meta_clean_Auto_F.bim

for chr in {1..22}
do
    plink2 \
    --bfile W_meta_clean_Auto \
    --chr $chr \
    --maf 0.005 \
    --make-bed \
    --out W_Auto_meta_chr$chr
done

#Run linear regression

chmod +x GWAS_MB_RUN.sh

./GWAS_MB_RUN.sh


#Run logistic regression

chmod +x GWAS_MB_RUN_Binary.sh

./GWAS_MB_RUN_Binary.sh


#Quantitative GWAS sumstats extract lead variants

for file in *_sorted.txt; do
    # Extract the base filename without the '_sorted.txt' suffix
    base="${file%_sorted.txt}"

    # Remove '#' from the header line and save to a temporary file
    sed '1s/^#//' "$file" > "${base}_modified.txt"

    # Run the pick.sh script with the modified file
    bash pick.sh -i "${base}_modified.txt" -I 3 -c 1 -p 2 -P 14 -o "${base}_lead.txt"

    # Run awk on the output file to create the .tab file
    awk -v OFS='\t' '{$1=$1; print}' "${base}_lead.txt" > "${base}_lead.tab"

    # Optional: Remove the temporary modified file
    rm "${base}_modified.txt"
done


#Extract only common variants
for file in *_sumstats.txt; do
    awk '$7 >= 0.05' "$file" > "${file%_sumstats.txt}_sorted_MAF05.txt"
done

#extract by pval
for file in *_MAF05.txt; do
  echo "Processing $file..."
  awk 'NR>1 && $15 <= 1E-05' "$file" > "${file%_MAF05.txt}_sug.pval.out"
  echo "Filtered file saved as ${file%_MAF05.txt}_sug.pval.out"
done


Extract_lead
for file in *_sorted.txt; do
    # Extract the base filename without the '_sorted.txt' suffix
    base="${file%_sorted_sug.pval.out.txt}"

    # Run the pick.sh script with the modified file
    bash pick.sh -i "$file" -I 3 -c 1 -p 2 -P 15 -o "${base}_maf05_lead.txt"

    # Run awk on the output file to create the .tab file
    awk -v OFS='\t' '{$1=$1; print}' "${base}_maf05_lead.txt" > "${base}_maf05_lead.tab"
done
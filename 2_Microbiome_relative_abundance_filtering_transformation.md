```bash

## Filtering on Phylum, Genus, Species and Bacterial Pathways
# Filtered relative abundance. Prevalence 10% and abundance 0.001%, at phylum, genus, species and bacterial pathways
# Phylum from 16S rRNA sequencing: 15 phyla —> after filtering: 7 phyla
# Genus from 16S rRNA sequencing: 142 genera —> after filtering: 63 genera
# Species from shotgun metagenome sequencing: 811 species —> after filtering: 120 species
# Bacterial pathways from shotgun metagenome mOTUs: 198 bacterial pathways —> after filtering: 181 bacterial pathways

## Two types of transformation and normalizetion of relative abundance data
# Centered log ratio transfromation (CLR) on filtered reltive abundances

Rscript CLR_transf.R 

options = commandArgs(trailingOnly = TRUE)
taxonomy_file = options[1]

if(!require(compositions)){
  install.packages("compositions")
  library(compositions)
}

tax = read.table(taxonomy_file,header=F,sep="\t")



clr = clr(tax,index = "clr")

all_div = data.frame(clr = clr)


write.table(all_div, file = "pheno_CLR.txt",sep="\t",row.names = F,quote = F)


# Inverse normal transformation (INT) on filtered reltive abundances 
#prior natural log transformation

Rscript Inverse_norm_trans.R

options = commandArgs(trailingOnly = TRUE)
taxonomy_file = options[1]
output_file = options[2]


tax = read.table(taxonomy_file,header=T,row.names=1,sep="\t")




corrected_data = apply(tax,2,function(x){
  numPhenos = length(which(!is.na(x))) 
  
  
  quantilePheno = (rank(x, na.last="keep", ties.method="random")-0.5)/numPhenos
  phenoIRNT = qnorm(quantilePheno)
  return(phenoIRNT)
})
corrected_data = as.data.frame(t(corrected_data))
corrected_data2 = data.frame(rownames(corrected_data),corrected_data)
colnames(corrected_data2)[1] = "-"


write.table(corrected_data2, file = output_file,sep="\t",row.names = F,quote = F)
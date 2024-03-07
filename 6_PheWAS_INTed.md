```bash

## PheWAS analysis on 692 LOF variants and INT 7 Phylum, 63 Genus, 120 species and 181 pathways. Sample shown for PheWAS analysis at Species level

library(PheWAS)

genotypes=read.table("W_SNPeff_LOF_05_SNPs.raw",header=TRUE)[,c(-2:-6)]
names(genotypes)[1]="IID"

phenotypes=read.csv("pheno_sp_phewas.csv", header=TRUE)


covariates=read.table("COV_META_PHEWAS.txt",header=TRUE)

results=phewas(phenotypes=phenotypes,genotypes=genotypes,covariates=covariates,cores=5)


results_d=addPhecodeInfo(results)

#List the significant results

sig_results <- results_d[results_d$bonferroni&!is.na(results_d$p),]
sig_table<- write.table(results, file = "sigLOF01sp.txt", sep = "\t",
                        row.names = TRUE, col.names = NA)

top<- results[order(results$p)[1:200],]
top_table<- write.table(top, file = "topsig_LOF01_sp.txt", sep = "\t",
                        row.names = TRUE, col.names = NA)


jpeg("PheWAS_LOF01_sp.jpeg")
phenotypeManhattan(results, annotate.phenotype.description=F, suggestive.line = 0.05, significant.line = 0.00001716,
                   OR.size = F, OR.direction = F,
                   y.axis.interval = 5, use.color=F,x.group.labels=F,x.phenotype.labels=T, y.axis.label = expression(-log[10](italic(p))))
dev.off()
```bash
##Two-sample MR between taxa relative abundance and traits from BBJ

##Exposure data:

escher_exp_dat <- read_exposure_data("UA_sug.txt", sep= "	")
head(escher_exp_dat)
escher_exp_dat <- clump_data(escher_exp_dat, pop="EAS")


##Outcome data:

ubiq_dat <- read_outcome_data(
  snps = escher_exp_dat$SNP,
  filename = "Dl.linear",
  sep = "	",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval"
)
head(ubiq_dat)

##Harmonize data:
dat <- harmonise_data(
  exposure_dat = escher_exp_dat, 
  outcome_dat = ubiq_dat
)

##MR:
res <- mr(dat)
res
pdf(file = "/Users/davidorey/Desktop/MR_UAvsDl.pdf",   
    width = 4, 
    height = 4) 


p1 <- mr_scatter_plot(res, dat)
p1[[1]]

dev.off()


## Forest plot:
pdf(file = "/Users/davidorey/Desktop/UAvsDl_forestplot.pdf",   
    width = 4, 
    height = 4) 

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]
dev.off()
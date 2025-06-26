```bash
#causal mediation analysis for 1000SNPs (rs671+999control SNPs) *2groups (stratified by genotypes)  
1000SNP list：random_rs671freq_SNPs
==============================================================================================================
R
library(dplyr)
library(mediation)
library(ggplot2)
library(ggsci)
library(reshape2)

file_list <- paste("random_rs671freq_SNPs")
SNP_list <- read.table(file=file_list,header=F)
SNP_list <- unlist(SNP_list)

finres <- c()
SNPres <- c()

for (A in c(1:length(SNP_list))){
SNP <- SNP_list[A]
file_name1 <- paste0(",SNP,"_major.txt")
file_name2 <- paste0(",SNP,"_minorcombo.txt")

d1 <- read.table(file=file_name1,header=T)
model.m1 <- glm(blom_IL10~INT_Xylose+age+M0F1,data=d1,family=gaussian(link="identity"))
model.y1 <- glm(HOMAIR~INT_Xylose+age+M0F1+blom_IL10,data=d1,family=gaussian(link="identity"))
out.l1 <- mediate(model.m1,model.y1,treat="INT_Xylose",mediator="blom_IL10",boot=T)
res1 <- c(summary(out.l1)$tau.coef,summary(out.l1)$tau.p,summary(out.l1)$d0,summary(out.l1)$d0.p)


d2 <- read.table(file=file_name2,header=T)
model.m2 <- glm(blom_IL10~INT_Xylose+age+M0F1,data=d2,family=gaussian(link="identity"))
model.y2 <- glm(HOMAIR~INT_Xylose+age+M0F1+blom_IL10,data=d2,family=gaussian(link="identity"))
out.l2 <- mediate(model.m2,model.y2,treat="INT_Xylose",mediator="blom_IL10",boot=T)
res2 <- c(summary(out.l2)$tau.coef,summary(out.l2)$tau.p,summary(out.l2)$d0,summary(out.l2)$d0.p)


result <- c(res1[3],res1[4],res2[3],res2[4],res1[3]/res2[3])
result <- matrix(data=result,nrow=1,ncol=5)
SNPres <- cbind(SNP,result)
finres <- rbind(finres,SNPres)
}

write.table(finres,file='rs671_similarfrqSNPs_ACMEratio.mod.txt', row.names=F, col.names=F, sep=" ",quote=FALSE)
plot <- read.table('rs671_similarfrqSNPs_ACMEratio.mod.txt',header=F)

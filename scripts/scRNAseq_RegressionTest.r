library(tidyr)
library(broom)
library(data.table)
library(dplyr)
library(RNOmni)

source("./SingleCell_RNAseq/eQTL/qq_plot.r")

mydata <- fread('counts.input.tsv')

TEST = 'meanHF' #sumNonsyn

mydata %>%
  group_by(variable) %>%
  do(tidy(lm(value ~ rankNorm(meanHF) + experiment + donor, .))) %>%
  mutate(Beta = as.character(round(estimate, 10)), "Pvalue" = round(p.value, 10), SE = round(std.error, 10)) %>%
  select(term, Beta, SE, "Pvalue") %>%
  as.data.frame() -> myresult
write.table(file = paste("myresult_", TEST, "tsv", sep=""), myresult, sep="\t", quote=FALSE, row.names=FALSE)

myresult <- subset(myresult, myresult$term == 'rankNorm(meanHF)')
myresult$geneIDlink <- gsub('_.*','',myresult$variable)
myresult_pvalue <- subset(myresult, myresult_ann$Pvalue < 0.001)
myresult_pvalue <- myresult_pvalue[order(myresult_pvalue$Pvalue), ]

write.table(file = paste("myresult_", TEST, "p0001.tsv", sep=""), myresult_pvalue, sep="\t", quote=FALSE, row.names=FALSE)

png(paste("myresult_", TEST, "png", sep=""))
qq.plot(myresult$Pvalue)
dev.off()




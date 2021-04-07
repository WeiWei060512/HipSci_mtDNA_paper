library(openxlsx)
library(tidyr)
library(data.table)
library(dplyr)
library(textshape)
library(MAST)

### DE analysis using MAST ###

mastTest <- function(mycounts, myvar, mycells, donor_name, var_pos){
    message(donor_name, '\t', var_pos)
    mydonor_cells <- subset(mycells, mycells$donor_short_id == donor_name)
    mydonor_counts <- as.matrix(mycounts %>% dplyr::select(rownames(mydonor_cells)))
    mydonor_var <- subset(myvar, (myvar$donor_short_id == donor_name & myvar$POS.x == var_pos))
    mydonor_cells$genotype <- ifelse(rownames(mydonor_cells) %in% mydonor_var$cell_name, "1","0")
    sca <- FromMatrix((mydonor_counts), mydonor_cells, mygene_list)
    cond<-factor(colData(sca)$genotype)
    cond<-relevel(cond,"0")
    colData(sca)$condition<-cond
    zlmCond <- zlm(~condition + cell_stage_short + pct_counts_endogenous, sca)
    summaryCond <- summary(zlmCond, doLRT='condition1') 
    summaryDt <- summaryCond$datatable  
    fcHurdle <- merge(summaryDt[contrast=='condition1' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast=='condition1' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
    fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    setorder(fcHurdle, fdr)
    write.table(file = paste0('MASTTestperDonorAllStages_', donor_name, 'pos', var_pos, '.output.tsv'), fcHurdle, sep='\t', row.names=TRUE, quote=FALSE)
}


myvar <- fread('./DEanalysis/outputMAST_Cov/single_experiments/scRNA_snps.txt')
mygene_list <- read.delim('gene_list.txt')
mygene_list$primerid <- mygene_list$geneIDs
mycounts <- readRDS('scRNA_counts.rds')
mycells <- fread('scRNA_meta.txt')
mycells <- subset(mycells, mycells$cell_name %in% colnames(mycounts))
mycells <- column_to_rownames(mycells, 'cell_name')
mycells$wellKey <- rownames(mycells)
mydonor_pos <- read.xlsx('./DEanalysis/outputMAST_Cov/single_experiments/donor_varpos_list.xlsx') # a list contains donor id and the variant tested

for(i in 1:nrow(myinput)) {
    mastTest(mycounts, myvar, mycells, myinput[i,1], myinput[i,2])
}


##### summary output ######
files <- list.files(path="./DEanalysis/outputMAST_Cov/single_experiments/", pattern="*.output.tsv", full.names=TRUE, recursive=FALSE)
lapply(files, function(x) {
    mydata <- read.delim(x)
    mydata <- subset(mydata, mydata$fdr < 0.1)
    mydata <- mydata[!grepl("RPL|RPS|MT-|MRP|MRPS", mydata$primerid),]
    if (nrow(mydata) > 0){
        mydata['files'] = x
        write.table(file=paste0(x, '.fdrbelow10per.tsv'), mydata, sep='\t', row.names=TRUE, quote=FALSE)
    }
})


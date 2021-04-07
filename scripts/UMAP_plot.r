library(M3C)

##### UMAP plots ######
output_dir <- './'
donor <- 'eipl_1'

filePath <- paste0(output_dir, donor, 'cluster.rds')
mydonor <- readRDS(filePath)
donor <- gsub('.*/','',filePath)
donor <- gsub('cluster.rds','',donor)

pdf(paste0(output_dir, donor, '_sumHF.cluster.pdf')) # width=8, height=5
umap(mydonor$data,labels=(mydonor$sumHF),controlscale=TRUE,scale=1,dotsize=3,printheight = 15, printwidth = 22)
dev.off()

pdf(paste0(output_dir, donor, '_stage.cluster.pdf')) 
umap(mydonor$data,labels=(mydonor$cell_stage_short),controlscale=TRUE,scale=3,dotsize=3)
dev.off()

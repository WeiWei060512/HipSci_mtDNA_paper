
### this script was used to process variance component analysis

library(data.table)
library(statmod); 
library(variancePartition)

##### define highly variable genes ##### 

highlyVar_genes <- function(mycounts, geneList) {
    means <- rowMeans(mycounts)
    vars <- apply(mycounts,1,var)
    cv2 <- vars/means^2
    minMeanForFit <- unname( quantile( means[ which( cv2 > .3 ) ], .95 ) )
    useForFit <- means >= minMeanForFit # & spikeins
    fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ),cv2[useForFit] )
    a0 <- unname( fit$coefficients["a0"] )
    a1 <- unname( fit$coefficients["a1tilde"])
    fit$coefficients
    ##Rank genes by the significance of deviation from the fit
    afit <- a1/means+a0
    varFitRatio <- vars/(afit*means^2)
    varorder <- order(varFitRatio,decreasing=T)
    myCounts <- mycounts[varorder,]
    geneList <- geneList[varorder,]
    rownames(myCounts) <- geneList
    return (myCounts)
}

##### variance component analysis ##### 

varComp <- function(mycounts, mymeta, mycounts_geneIDs, myoutput) {
    form <- ~ nonsyn.mt + (1|donor_id) + (1|experiment) + (1|cell_stage)
    varPart <- fitExtractVarPartModel(mycounts[1:4000, ], form, mymeta)
    rownames(varPart) <- mycounts_geneIDs[1:4000]
    pdf(paste0('boxplot_', myoutput, '.pdf'))
    vp <- sortCols( varPart )
    plotVarPart( vp )
    dev.off()
    write.table(file = paste0('varPart_', myoutput, '.tsv'), varPart, row.names=TRUE, sep='\t', quote=FALSE)
}

varComp_nonsyn <- function(mycounts_sub, mymeta_sub, mycounts_geneIDs, myoutput) {
    form <- ~ nonsyn.mt
    varPart <- fitExtractVarPartModel(mycounts_sub[1:4000, ], form, mymeta_sub)
    rownames(varPart) <- mycounts_geneIDs[1:4000]
    pdf(paste0('boxplot_', myoutput, '.pdf'))
    vp <- sortCols( varPart )
    plotVarPart( vp )
    dev.off()
    write.table(file = paste0('varPart_', myoutput, '.tsv'), varPart, row.names=TRUE, sep='\t', quote=FALSE)
}


mycounts <- readRDS('scRNA_counts.rds')
mymeta <- fread('scRNA_meta.txt')
geneList <- fread('gene_list.txt')

mydonor_list <- c('bima_1', 'ciwj_2', 'eesb_1', 'eipl_1', 'fafq_1', 'fiaj_3', 'garx_2', 'hajc_1', 'hayt_1', 'hecn_3', 'hoik_1', 'iisa_3', 'iiyk_4', 'je    jf_2', 'joxm_1', 'juuy_2', 'kajh_3', 'keui_1', 'kolf_2', 'kuco_1', 'letw_1', 'liqa_1', 'lise_3', 'meue_4', 'miaj_6', 'naah_2', 'naju_    1', 'oilg_3', 'pahc_4', 'pipw_5', 'puie_5', 'qaqx_1', 'qayj_3', 'qehq_3', 'quls_2', 'rayr_1', 'rutc_2', 'seru_1', 'sohd_3', 'sojd_3',     'toco_5', 'tolg_6', 'uilk_3', 'vass_1', 'vazt_1', 'vuna_3', 'wuye_2', 'xojn_3', 'yelp_3', 'yoch_6', 'zerv_8')

cell_stage = 'defendo' # 'iPSC','mesendo'
for (mydonor in mydonor_list){
    print (mydonor)
    mymeta_sub <- subset(mymeta, mymeta$donor_id == mydonor)
    mycounts_sub <- mycounts %>% dplyr::select(mymeta_sub$cell_name)
    mycounts_sub <- highlyVar_genes(mycounts_sub, geneList)
    mycounts_geneIDs <- row.names(mycounts_sub)
    varComp_nonsyn(mycounts_sub, mymeta_sub, mycounts_geneIDs, paste0(mydonor,'_', cell_stage))
}   

##### plots #####
varComp_out_sub <- subset(varComp_out, varComp_out$donor_id %in% c('joxm_1'))
ggplot(varComp_out_sub, aes(x=cell_stage, y=nonsyn.mt, fill=cell_stage)) + geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill='white')+
    theme_classic() + ylim(0,1) + scale_fill_manual(values=c('#b85042','#12a4d9','#ffc13b'))

#RaceID4 
#create folders
dir.create("data")
dir.create("data/counts")
dir.create("data/timepoints")
dir.create("plots")
dir.create("plots/heatmaps")
dir.create("plots/others")
dir.create("plots/tsne")

library(tidyverse)
library(viridis)
library(RaceID)
library(Seurat)
library(Matrix)

date = Sys.Date()

non_micr <- unlist(c(read_csv("data/non-microglia-cell-ids.csv")[2], 
                     read_csv("data/more-non-microglia-cell-ids.csv")[2], 
                     read_csv("data/more-more-non-microglia-cell-ids.csv")[2], 
                     read_csv("data/more-more-non-microglia-cell-ids.csv")[2]))

#add hypothalamic microglia
hyth <- list()
hyth[[1]] <- read.csv('data/counts/A761.all_samples.gencode_genomic.corrected_merged.csv', stringsAsFactors = F, sep = '\t', row.names = "GENEID")

hyth[[2]] <- read.csv('data/counts/P215_A326.all_samples.gencode_genomic.corrected_merged.csv', stringsAsFactors = F, sep = '\t', row.names = "GENEID")
hyth <- lapply(hyth, function(x) {
    rownames(x) = gsub("_.*", "", rownames(x))
    x = x[,!colnames(x) %in% non_micr]
})

#generate data frame with timepoints
timepoints <- data.frame(diet=c("1d", "3d"), pattern=c("(9wks_CT|9wks_HFD1d)", "(14wks_HFD3d|14wks_NCD)"), stringsAsFactors = F)

save(timepoints, file="data/timepoints.Robj")

for (i in 1:nrow(timepoints)) {
#RaceID
sc <- SCseq(hyth[[i]][,grepl(timepoints[[2]][i], colnames(hyth[[i]]))])

# filtering of expression data
sc <- filterdata(sc, 
                 mintotal=500,
)

sc <- CCcorrect(sc, 
                dimR = T, 
                nComp = 20,
                CGenes = c('Jun',
                           'Fos',
                           'Zfp36',
                           'Atf3',
                           'Hspa1a|Hspa1b',
                           'Dusp1',
                           'Egr1',
                           'Malat1'))

    sc <- compdist(sc,metric="pearson")
    sc <- clustexp(sc) 
    sc <- findoutliers(sc)
    sc <- comptsne(sc)
    sc <- compfr(sc,knn=10)

#Save sc file
save(sc, file = paste0('data/timepoints/sc', '-with-mt-genes-', timepoints[[1]][i],'.Robj'))

}

#analysis without mt genes
hyth <- lapply(hyth, function(x) {
    x = x[!grepl("mt-", rownames(x)),]
})

for (i in 1:nrow(timepoints)) {
    #RaceID
    sc <- SCseq(hyth[[i]][,grepl(timepoints[[2]][i], colnames(hyth[[i]]))])
    
    # filtering of expression data
    sc <- filterdata(sc, 
                     mintotal=500,
    )
    
    sc <- CCcorrect(sc, 
                    dimR = T, 
                    nComp = 20,
                    CGenes = c('Jun',
                               'Fos',
                               'Zfp36',
                               'Atf3',
                               'Hspa1a|Hspa1b',
                               'Dusp1',
                               'Egr1',
                               'Malat1'))
    
    sc <- compdist(sc,metric="pearson")
    sc <- clustexp(sc) 
    sc <- findoutliers(sc)
    sc <- comptsne(sc)
    sc <- compfr(sc,knn=10)
    
    #Save sc file
    save(sc, file = paste0('data/timepoints/sc', '-without-mt-genes-', timepoints[[1]][i],'.Robj'))
    
}


plotexpmap(sc,"Mrc1",logsc=F,fr=F)
plotexpmap(sc,"Lyve1",logsc=F,fr=F)
plotexpmap(sc,"Cd163",logsc=F,fr=F)
plotexpmap(sc,"Tmem119",logsc=F,fr=F)
plotexpmap(sc,"Cx3cr1",logsc=F,fr=F)
plotexpmap(sc,"Hexb",logsc=F,fr=F)
plotexpmap(sc,"Ptprc",logsc=F,fr=F)
plotexpmap(sc,"Cd3e",logsc=F,fr=F)
plotexpmap(sc,"Itgam",logsc=F,fr=F)
plotexpmap(sc,"Cd8a",logsc=F,fr=F)
plotexpmap(sc,"Cd4",logsc=F,fr=F)
plotexpmap(sc,"H2-Aa",logsc=F,fr=F)
plotexpmap(sc,"Zbtb46",logsc=F,fr=F)
plotexpmap(sc,"Mog",logsc=F,fr=F)
plotexpmap(sc,"Mbp",logsc=F,fr=F)
plotexpmap(sc,"Gfap",logsc=F,fr=F)
plotexpmap(sc,"Wfdc17",logsc=F,fr=F)
plotexpmap(sc,"Cd79a",logsc=F,fr=F)
plotexpmap(sc,"Cst3",logsc=F,fr=F)
plotexpmap(sc,"Nkg7",logsc=F,fr=F)
plotexpmap(sc,"Rho",logsc=F,fr=F)
plotexpmap(sc,"S100a11",logsc=F,fr=F)
plotexpmap(sc,"Fn1",logsc=F,fr=F)
plotexpmap(sc,"Gfap",logsc=F,fr=F)
plotexpmap(sc,"Cd209a",logsc=F,fr=F)
plotexpmap(sc,"Rho",logsc=F,fr=F)

plotexpmap(sc,"Mrc1",logsc=F,fr=T)
plotexpmap(sc,"Lyve1",logsc=F,fr=T)
plotexpmap(sc,"Cd163",logsc=F,fr=T)
plotexpmap(sc,"Tmem119",logsc=F,fr=T)
plotexpmap(sc,"Cx3cr1",logsc=F,fr=T)
plotexpmap(sc,"Ptprc",logsc=F,fr=T)
plotexpmap(sc,"Cd3e",logsc=F,fr=T)
plotexpmap(sc,"Itgam",logsc=F,fr=T)
plotexpmap(sc,"Cd8a",logsc=F,fr=T)
plotexpmap(sc,"H2-Aa",logsc=F,fr=T)
plotexpmap(sc,"Zbtb46",logsc=F,fr=T)
plotexpmap(sc,"Ly6c2",logsc=F,fr=T)
plotexpmap(sc,"Cd177",logsc=F,fr=T)
plotexpmap(sc,"Igkc",logsc=F,fr=T)
plotexpmap(sc,"Wfdc17",logsc=F,fr=T)
plotexpmap(sc,"Plp1",logsc=F,fr=T)
plotexpmap(sc,"Mog",logsc=F,fr=T)
#plot marker genes


#identify myeloid cells
#write.csv(data.frame(names(sc@cpart[sc@cpart %in% c(7)])), 'data/non-microglia-cell-ids.csv')
write.csv(data.frame(names(sc@cpart[sc@cpart %in% c(6,10)])), 'data/more-more-more-non-microglia-cell-ids.csv')

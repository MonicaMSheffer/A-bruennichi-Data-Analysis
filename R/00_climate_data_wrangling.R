#This script takes GPS data for each spider we collected, turns it into a spatial data file so it can work with the bioclim files (downloaded from worldclim), extracts the data for each bioclimatic variable for each site, adds it into a table with the spider data, then runs a PCA on the environmental data per collecting site, to output a final file of the adult spider phenotypes and the environmental principle components that correspond to their collecting site.

#load libraries
library(sf);library(stars);library(raster);library(mapview);library(tidyverse);library(raster);library(tmap);library(spData);library(factoextra);library(tidyverse);library(dplyr);library(openxlsx)


#Need to take non-spatial spatial data from a GPS point in a .csv file and turn it into a spatial file
spider_csv = read.csv("data/rad_phenotypes_completed_sortedLat.csv")
#View(spider_csv)
spider_noNA = dplyr::filter(spider_csv, !is.na(GPS_N) & !is.na(GPS_E))
#View(spider_noNA)
spider_sf = st_as_sf(spider_noNA,
                     coords = c("GPS_E", "GPS_N"), # has to be long, lat order for this projection
                     crs = 4326) # this is for WGS84
# save our new spatial file
# write_sf(spider_sf, 'output/spider_occ.gpkg')

#import information on bioclim variables
bioclim_table=read.csv("data/bioclim_info.csv")

#import spatial file for spider sites
spider <- read_sf('output/spider_occ.gpkg')

#extract bioclim data for each individual
bioclim_list=list()
for (i in 1:length(bioclim_table$BioVar)) {
  filepath=paste0("data/bioclim/wc2.1_30s_bio_",i,".tif")
  world_bioclim = raster(filepath)
  bioclim_list[[i]] = raster::extract(world_bioclim, spider)
}
names(bioclim_list) = bioclim_table$varDescrip
spider=cbind(spider,do.call(cbind,bioclim_list))

#we are focusing our analysis for this manuscript on the European range of the species; remove sites in Japan and Russia
spider_final<-spider[which(spider$Country!="Japan"&spider$Country!="Russia"),]

st_write(spider_final,dsn="output/spider_allBioClim.csv",append=FALSE) #write to new file

#to visualize data, look at any of the bioclim variables interactively, one at a time
options(viewer=NULL)
mapview(spider_final["annMeanTemp"]) 

#make plots of each bioclim variable, save as PDF for paper/supplement
source("R/myplot.R")
pdf("output/bioclim_plots_all.pdf")
for (i in 1:length(bioclim_table$varDescrip)){
  obj=bioclim_table$varDescrip[i]
  print(obj)
  myplot(obj,spider_final)
}
dev.off()

###Use climate data for a PCA to summarize variation into components for use in the rest of the data analysis
dat=read.csv("output/spider_allBioClim.csv")

#dat file has a separate row for each individual, but we just need one value per site, so filter by distinct population numbers
dat_unique = dat %>% distinct(population_number,.keep_all = TRUE)

#we only need bioclimatic data for sites in Europe, so filter out Japan and Russia
dat_filtered=dat_unique[which(dat_unique$Country!="Russia" & dat_unique$Country!="Japan"),]

#group Estonia, Latvia and Lithuania into "Baltic" for cleaner plotting
dat_filtered$groups=dat_filtered$Country
dat_filtered$groups[which(dat_filtered$Country=="Estonia"|dat_filtered$Country=="Lativa"|dat_filtered$Country=="Lithuania")]="Baltic"


#calculate PCA
res.pca_noJnoR=prcomp(dat_filtered[,19:37],scale=T)

#vizualize PCA 
fviz_eig(res.pca_noJnoR)
fviz_pca_biplot(res.pca_noJnoR,axes=c(1,2),col.ind=dat_filtered$groups,repel=T,addEllipses = TRUE,mean.point=F,geom="point",shapes=c(13,14,15,16,17,18,19))
fviz_pca_biplot(res.pca_noJnoR,axes=c(1,3),col.ind=dat_filtered$groups,repel=T,addEllipses = TRUE,mean.point=F,geom="point",shapes=c(13,14,15,16,17,18,19))
fviz_pca_biplot(res.pca_noJnoR,axes=c(1,4),col.ind=dat_filtered$groups,repel=T,addEllipses = TRUE,mean.point=F,geom="point",shapes=c(13,14,15,16,17,18,19))
fviz_pca_biplot(res.pca_noJnoR,axes=c(1,5),col.ind=dat_filtered$groups,repel=T,addEllipses = TRUE,mean.point=F,geom="point",shapes=c(13,14,15,16,17,18,19))



fviz_pca_biplot(res.pca_noJnoR,axes=c(1,3),col.ind=dat_filtered$groups,repel=T,addEllipses = TRUE,mean.point=F,geom="point")
fviz_pca_biplot(res.pca_noJnoR,axes=c(2,3),col.ind=dat_filtered$groups,repel=T,addEllipses = TRUE,mean.point=F,geom="point")
fviz_pca_biplot(res.pca_noJnoR,axes=c(3,4),col.ind=dat_filtered$groups,repel=T,addEllipses = TRUE,mean.point=F,geom="point")
fviz_pca_biplot(res.pca_noJnoR,axes=c(4,5),col.ind=dat_filtered$groups,repel=T,addEllipses = TRUE,mean.point=F,geom="point")

#look at loadings of specific PCs
sort(res.pca_noJnoR[["rotation"]][,"PC1"])
sort(res.pca_noJnoR[["rotation"]][,"PC2"])
sort(res.pca_noJnoR[["rotation"]][,"PC3"])
sort(res.pca_noJnoR[["rotation"]][,"PC4"])
sort(res.pca_noJnoR[["rotation"]][,"PC5"])

#make table of loadings
pc1_loadings=t(as.data.frame(res.pca_noJnoR[["rotation"]][,"PC1"]))
pc2_loadings=t(as.data.frame(res.pca_noJnoR[["rotation"]][,"PC2"]))
pc3_loadings=t(as.data.frame(res.pca_noJnoR[["rotation"]][,"PC3"]))
pc4_loadings=t(as.data.frame(res.pca_noJnoR[["rotation"]][,"PC4"]))
pc5_loadings=t(as.data.frame(res.pca_noJnoR[["rotation"]][,"PC5"]))
pc_loadings=as.data.frame(rbind(pc1_loadings,pc2_loadings,pc3_loadings,pc4_loadings,pc5_loadings))
pc_loadings$PC=c("PC1","PC2","PC3","PC4","PC5")
write.xlsx(pc_loadings,'output/pc_loadings.xlsx')

#add PC data for the first 5 PCs into the adult spider data table
results_pca=as.data.frame(res.pca_noJnoR[["x"]])
popname_results=cbind(dat_filtered$population_number,results_pca)
popname_results$`dat_filtered$population_number`
pca_index=match(dat$population_number,popname_results$`dat_filtered$population_number`)
#View(popname_results[pca_index,])

dat_final=cbind(dat[,1:18],popname_results[pca_index,])
write.csv(dat_final,"output/phenotype_environmentPCA_table.csv")

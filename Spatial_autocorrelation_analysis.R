library(data.table)
library(ggcorrplot)
library(colorspace)
library(Hmisc)
library(rgdal)
library(spdep)
library(spatialreg)
library(car)
library(cAIC4)
require(lme4)
require(randomForest)
library(sp)

#1. Spatial Regression for Viral Diversity and Environmental Factors:
meta <- read.csv("meta_scale.csv")
latlon <- meta[,c(3,4)]
coordinates(latlon) <- ~lon + lat
nb <- dnearneigh(latlon, d1=0, d2=26)
w <- nb2listw(nb, style = "W")

efalpha <- read.csv("efalpha.csv")
efalpha <- efalpha[efalpha$effect_size >= 0.4,]
meta <- meta[,c("shannon_species",efalpha$column[c(2:15,17:22)])]

spearman_res_ori <- NULL
for (i in c(2:21)){
col <- colnames(meta)[i]
print(col)
slm <- lagsarlm(shannon_species ~ get0(col),data=meta,listw=w)
summary(slm)
spearman_res <- cor.test(resid(slm), meta[!is.na(meta[,i]),i], method="spearman")
print(spearman_res)
spearman_original <- cor.test(meta$shannon_species, meta[,i], method="spearman")
print(spearman_original)
spearman_res_ori <- rbind(spearman_res_ori,c(colnames(meta)[i],spearman_res$estimate,spearman_res$p.value,spearman_original$estimate,spearman_original$p.value))
}

#2.	Spatial Regression for Viral and Microbial Community Diversity:
microbialdiversity <- read.csv("nonpareil_diversity.csv")
viraldiversity <- read.csv("meta_scale.csv")
viraldiversity <- viraldiversity[,c("Sample.id","Biome","shannon_species","lat","lon")]
viraldiversity <- viraldiversity[pmatch(viraldiversity$Sample.id,microbialdiversity$sample),]
vh <- data.frame(cbind(viraldiversity$lat,viraldiversity$lon,viraldiversity$shannon_species, microbialdiversity$diversity))
vh <- na.omit(vh)
latlon <- vh[,c(2,1)]

colnames(latlon) <- c("lon","lat")
coordinates(latlon) <- ~lon + lat
nb <- dnearneigh(latlon, d1=0, d2=26)
w <- nb2listw(nb, style = "W")

slm <- lagsarlm(X3 ~ X4,data=vh,listw=w)
spearman_res <- cor.test(resid(slm), vh$X4, method="spearman")
spearman_original <- cor.test(vh$X3, vh$X4,method="spearman")
print(spearman_res)
print(spearman_original)

#3. Random Forest Model Residuals:
set.seed(11)
meta <- read.csv("meta_rfmodel.csv",row.names = 1)
shannon.rf <- randomForest(shannon.red ~ lat+ref_6+bio4+ai_et0+contrast+ref_4+shannon+bio5+bio17+silt2_1+sand2_1, 
                           data = meta,na.action = na.omit)
shannon.rf.pred <- predict(shannon.rf, meta)

latlon <- meta[!is.na(shannon.rf.pred),c(4,3)]
coordinates(latlon) <- ~lon + lat
w <- knn2nb(knearneigh(latlon, k = 5))
w <- nb2listw(w, style = "W")

residues <- meta$shannon.red[!is.na(shannon.rf.pred)]- shannon.rf.pred[!is.na(shannon.rf.pred)]
moran_test <- moran.test(residues, listw = w)
print(moran_test)

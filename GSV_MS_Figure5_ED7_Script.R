#Model establishment
# prediction models
library(data.table)
library(car)
library(cAIC4)

## split datasets
## select main factors 
## compute multiple covariate groups
metadata <- read.csv("envdata_v3.csv",row.names = 1,as.is = TRUE)
metadata.1 <- read.csv("metadata_20210107.csv",row.names = 1,as.is = TRUE)
quanti.1 = metadata[,6:85]
quanti.2 = metadata.1[,6:5]
quanti <- cbind(quanti.1,quanti.2)

rm.names <- c("New_consensus_full_class_9", "New_consensus_full_class_7", "New_nppmean" , 
  "New_Gppmean" ,"New_slope","New_Mean_of_annual_predicted_ECe_1980_2018" ,
  "New_Mean_of_annual_predicted_ESP_1980_2018")

quanti.rm <- na.omit(quanti[ ,!colnames(quanti) %in% rm.names])
tree <- hclustvar(X.quanti = quanti.rm[,1:73], X.quali = quanti.rm[,74:75])
stab <- stability(tree, B = 200) 

## cluster number k = 5
tree.cut5 <- cutreevar(tree, k = 5)

## cluster number k = 11
tree.cut15 <- cutreevar(tree, k = 15)
melt(sort(tree.cut15$cluster))->a
write.csv(a, 'cov.txt',sep = "\t")
env.rep <- quanti.rm[ , match(1:11,tree.cut11$cluster)]

## transfer metadata into sp::SpatialPointsDataFrame
metadata <- read.csv("../Metadata/metadata_new.csv",row.names = 1,as.is = TRUE)
shannon <- fread(file = "shannon.csv")
abund <- fread(file = "vir_abundance.txt")
shannon.new <- shannon[pmatch(metadata$Sample.id , shannon$V1), ]
abund.new <- abund[pmatch(metadata$Sample.id , abund$V1), ]

xy <- metadata[ , c(4,3)]
env.meta <- cbind(metadata[ , c(2:15)],shannon.new$x, abund.new$x)
colnames(env.meta) <- c("project", "lat","long","northness","continent",
                        "mat","map","altitude","uvi","npp","soil_pH",
                        "moisture","soc","evapotion","shannon","richness")

env.meta <- env.meta[env.meta$richness < 15000 , ]

shannon.red <- tapply(env.meta$shannon, 
                      INDEX = paste(env.meta$lat,env.meta$long), 
                      mean)

richness.red <- tapply(env.meta$richness, 
                       INDEX = paste(env.meta$lat, env.meta$long), 
                       mean)

env.meta.red <- cbind(env.meta[!duplicated(xy),c(2:14)],
                      shannon= shannon.red,
                      richness = richness.red)

env.meta.red[,5:13] <- scale(env.meta.red[,5:13])
env.meta.red$richness <- log10(env.meta.red$richness+1)

## least squre model
require(lme4)
shannon.lmer <- lmer(shannon~altitude+(1|altitude)+(1|northness)+(1|npp)+(1|uvi)+(1|mat)+
                       (1|map)+(1|evapotion)+(1|soil_pH)+(1|soc)+(1|evapotion)+
                       northness+npp+uvi+(northness*npp)+(npp*northness)+(northness*uvi)+(npp*uvi)+(uvi*npp)+(uvi*northness)+
                       mat+map+evapotion+(mat*map)+(mat*evapotion)+(map*mat)+(map*evapotion)+(evapotion*mat)+(evapotion*map)+
                       soil_pH+soc+moisture+(soil_pH*soc)+(soil_pH*moisture)+(soc*soil_pH)+(soc*moisture)+(moisture*soil_pH)++(moisture*soc),  
                     data = env.meta.red,control = lmerControl(calc.derivs = FALSE))


shannon.lmer <- lmer(shannon~altitude+(1|altitude)+(1|northness)+(1|npp)+(1|uvi)+(1|mat)+
                       (1|map)+(1|evapotion)+(1|soil_pH)+(1|soc)+(1|evapotion)+
                       northness+npp+uvi+(northness*npp)+(npp*northness)+(northness*uvi)+(npp*uvi)+(uvi*npp)+(uvi*northness)+
                       mat+map+evapotion+(mat*map)+(mat*evapotion)+(map*mat)+(map*evapotion)+(evapotion*mat)+(evapotion*map)+
                       soil_pH+soc+moisture+(soil_pH*soc)+(soil_pH*moisture)+(soc*soil_pH)+(soc*moisture)+(moisture*soil_pH)+dp+(moisture*soc),  
                     data = env.meta.red,control = lmerControl(calc.derivs = FALSE))


vif(shannon.lmer)#5.0，3.0，6.0

shannon.lmer <- lmer(shannon~altitude+(1|altitude)+(1|northness)+(1|npp)+(1|uvi)+(1|mat)+
                       (1|map)+(1|evapotion)+(1|soil_pH)+(1|soc)+(1|evapotion)+
                       northness+npp+uvi+(northness*npp)+(npp*northness)+(northness*uvi)+(npp*uvi)+(uvi*npp)+(uvi*northness)+
                       mat+map+(mat*map)+(map*mat)+
                       soil_pH+soc+moisture+(soil_pH*soc)+(soil_pH*moisture)+(soc*soil_pH)+(soc*moisture)+(moisture*soil_pH)++(moisture*soc),  
                     data = env.meta.red,control = lmerControl(calc.derivs = FALSE))
vif(shannon.lmer)
shannon.lmer <- lmer(shannon~altitude+(1|altitude)+(1|npp)+(1|uvi)+(1|mat)+
                       +(1|soil_pH)+(1|soc)+
                       northness+npp+uvi+(northness*npp)+(npp*northness)+(northness*uvi)+(npp*uvi)+(uvi*npp)+(uvi*northness)+
                       mat+
                       soil_pH+soc+moisture+(soil_pH*soc)+(soil_pH*moisture)+(soc*soil_pH)+(soc*moisture)+(moisture*soil_pH)++(moisture*soc), 
                     data = env.meta.red,control = lmerControl(calc.derivs = FALSE))
vif(shannon.lmer)
cAIC(shannon.lmer)#去除后变小，删

shannon.lmer <- lmer(shannon~altitude+(1|altitude)+(1|npp)+(1|uvi)+(1|mat)+
                       (1|soil_pH)+(1|moisture)+northness+npp+uvi+
                       (npp*uvi)+mat+soil_pH+moisture,
                     data = env.meta.red, 
                     control = lmerControl(calc.derivs = FALSE))
cAIC(shannon.lmer)

# cross-validation
env.meta.red <- na.omit(env.meta.red)

shannon.pred.10f <- NULL
for(i in 1:437){
  shannon.lmer <- lmer(shannon~0+altitude+(1|altitude)+(1|npp)+(1|uvi)+(1|mat)+
                         (1|soil_pH)+(1|moisture)+npp+uvi+northness+
                         (npp*uvi)+mat+soil_pH+moisture,
                       data = env.meta.red[-i, ], 
                       control = lmerControl(calc.derivs = FALSE))
  
  shannon.pred <- predict(shannon.lmer, 
                          env.meta.red[i,],
                          re.form = NULL,
                          allow.new.levels = TRUE,
                          type = "response")
  shannon.pred.10f <- c(shannon.pred.10f, shannon.pred)
  print(i)
}

plot(env.meta.red[,14],
     shannon.pred.10f,
     xlim=c(1,8),
     ylim=c(1,8))
abline(0,1)

chisq.test(env.meta.red[,14],
           richness.pred.10f[-c(1:10)])

## richness

richness.lmer <- lmer(richness~altitude+(1|altitude)+(1|northness)+(1|npp)+(1|uvi)+(1|mat)+
                        (1|map)+(1|soil_pH)+(1|soc)+
                        northness+npp+uvi+(northness*npp)+(npp*northness)+(northness*uvi)+(npp*uvi)+(uvi*npp)+(uvi*northness)+
                        mat+map+(mat*map)+
                        soil_pH+soc+moisture+(soil_pH*soc)+(soil_pH*moisture)+(soc*soil_pH)+(soc*moisture)+(moisture*soil_pH)++(moisture*soc),  
                      data = env.meta.red,control = lmerControl(calc.derivs = FALSE))

vif(richness.lmer)

richness.lmer <- lmer(richness~altitude+
                        (1|altitude)+
                        (1|npp)+
                        (1|uvi)+
                        (1|mat)+
                        (1|soil_pH)+
                        (1|soc)+
                        northness+
                        npp+
                        uvi+
                        (northness*npp)+
                        (northness*uvi)+
                        (npp*uvi)+
                        mat+
                        #map+
                        soil_pH+
                        soc+
                        moisture+
                        (soil_pH*soc)+
                        (soil_pH*moisture)+
                        (soc*moisture),
                      data = env.meta.red,
                      control = lmerControl(calc.derivs = FALSE))

vif(richness.lmer)
cAIC(richness.lmer)

richness.lmer <- lmer(richness~(1|altitude)+
                        (1|npp)+
                        (1|uvi)+
                        (1|mat)+
                        (1|soil_pH)+
                        (1|soc)+
                        northness+
                        npp+
                        uvi+
                        (northness*npp)+
                        (npp*uvi)+
                        soil_pH+
                        soc+
                        moisture+
                        (soil_pH*soc)+
                        (soil_pH*moisture) +
                        (soc*moisture),
                      data = env.meta.red, 
                      control = lmerControl(calc.derivs = FALSE))
cAIC(richness.lmer)

richness.pred.10f <- NULL
for(i in 1:437){
  richness.lmer <- lmer(richness~(1|altitude)+(1|npp)+(1|uvi)+(1|mat)+(1|soil_pH)+
                          northness+npp+uvi+(northness*npp)+(npp*uvi)+
                          soil_pH+soc+moisture+(soil_pH*soc)+(soil_pH*moisture)+(soc*moisture),
                        data = env.meta.red, 
                        control = lmerControl(calc.derivs = FALSE))
  
  richness.pred <- predict(richness.lmer, 
                           env.meta.red[i,],
                           re.form = NULL,
                           allow.new.levels = TRUE,
                           type = "response")
  richness.pred.10f <- c(richness.pred.10f, richness.pred)
  print(i)
}

plot(env.meta.red[,15],
     richness.pred.10f,
     cex = .2,
     xlim=c(1,4),
     ylim=c(1,4)
)
abline(0,1)

rmse(shannon.lmer, env.meta.red)
rsquare(shannon.lmer, env.meta.red)

rmse(richness.lmer, env.meta.red)
rsquare(richness.lmer, env.meta.red)

## random forest
require(randomForest)

shannon.rf <- randomForest(shannon ~ lat+long+mat+map+altitude+uvi+npp+soil_pH+moisture+soc+evapotion, 
                           data = env.meta.red,
                           na.action = na.omit)
shannon.rf.pred <- predict(shannon.rf, env.meta.red)

shannon.rf.pred.10f <- NULL

for(i in 1:437){
  shannon.rf <- randomForest(shannon ~ lat+long+mat+map+altitude+uvi+npp+soil_pH+moisture+soc+evapotion, 
                             data = env.meta.red[-i,],
                             na.action = na.omit)
  shannon.rf.pred <- predict(shannon.rf, env.meta.red[i,])
  shannon.rf.pred.10f <- c(shannon.rf.pred.10f, shannon.rf.pred)
}

plot(env.meta.red[,14],
     shannon.rf.pred.10f,
     xlim=c(0,8),
     ylim=c(0,8))
abline(0,1)

rmse(model = shannon.rf,
     data = na.omit(env.meta.red))

rsquare(model = shannon.rf,
        data = na.omit(env.meta.red))

richness.rf <- randomForest(richness~lat+long+mat+map+altitude+uvi+npp+soil_pH+moisture+soc+evapotion, 
                            data = env.meta.red,
                            na.action=na.omit)

richness.rf.pred <- predict(richness.rf, env.meta.red)

richness.rf.pred.10f <- NULL
set.seed(333)
seq <- sample(1:518)
for(i in 0:9){
  richness.rf <- randomForest(richness ~ lat+long+mat+map+altitude+uvi+npp+soil_pH+moisture+soc+evapotion, 
                              data = env.meta.red[-seq[c(1:50+i*50)],],
                              na.action = na.omit)
  richness.rf.pred <- predict(richness.rf, env.meta.red[seq[c(1:50+i*50)],])
  richness.rf.pred.10f <- c(richness.rf.pred.10f, richness.rf.pred)
}

plot(env.meta.red[names(richness.rf.pred.10f),15],
     richness.rf.pred.10f,
     xlim=c(0,5),
     ylim=c(0,5))
abline(0,1)


#Calculation
require(lme4)
require(modelr)
library(data.table)

env.meta <- read.csv("$meta.csv",row.names = 1,as.is = TRUE)

index.lmer <- lmer(output_index ~ formula
                   data = env.meta, 
                   control = lmerControl(calc.derivs = FALSE))

factor.a <- fread("$factor.a.csv") #Input all the environmental layers
factor.a <- data.frame(factor.a)
factor.b <- fread("$factor.b.csv") 
factor.b <- data.frame(factor.b)

a <- matrix(nrow=18000,ncol=36000) 

for (i in 1:36000){
  env.meta <- data.frame(factor.a=factor.a[,i],factor.b=factor.b[,i])#Merge column i of every factor into a whole array in data.frame format.
                         
                         colnames(env.meta)<-c("factor.a","factor.b")
                         
                         abundance.pred <- predict(index.lmer, #Model generated in the above
                                                   env.meta, 
                                                   re.form = NULL,
                                                   allow.new.levels = TRUE,
                                                   type = "response")
                         a[,i] <- abundance.pred
}

a[is.na(a)==TRUE] <- 0
fwrite(a,"cal_index.csv",col.names=F,sep=",")

#map Visualization
library(sf)
library(raster)
library(dplyr)
library(spData)
library(spDataLarge)
library(tmap)

nz <- raster("index.tif")
tm_shape(nz)+
  tm_raster(style = "fisher")

#Validation-bootstrapping procedure
library(boot)
library(data.table)
require(lme4)
require(modelr)

#generate env.meta.all
factor.a <- fread("$factor.a.csv") #Input all the environmental layers
factor.a <- as.matrix(factor.a)
factor.a <- as.vector(factor.a)
factor.b <- fread("$factor.b.csv") 
factor.b <- as.matrix(factor.b)
factor.b <- as.vector(factor.b)

env.meta.all <- cbind(factor.a,factor.b)
colnames(env.meta.all) <- c("factor.a","factor.b")

#input env.meta
env.meta <- read.csv("$meta.csv",row.names = 1,as.is = TRUE)
#Environmental information of all the sample sites, and every row represents a factor and the index needed to calculate. Be sure that rownames should be similar with names in the input model.

#Stratified bootstrapping
##function generation
cal <- function(env.meta,ind){
  d <- env.meta[ind,]
  index.lmer <- lmer(index~factor.a+factor.b,
                     data = d, control = lmerControl(calc.derivs = FALSE))#Input the model generated above.
  index.pred <- predict(index.lmer, 
                        env.meta.all, # env.meta.all is global environmental data in format as "meta_filter.csv"
                        re.form = NULL,
                        allow.new.levels = TRUE,
                        type = "response")
  result <- t(index.pred)
  write.table(result,"result.csv",append = T,col.names=F,sep=",")
}
##tansfer into vector
env.meta$northness<-as.factor(env.meta$factor.x)#The column used as stratification category.
##bootstrap
test.boot <- boot(env.meta,cal,R = 100,strata = env.meta.new$factor.x,ncpus=28)
#R——bootstrap iteration times  strata——The column used as stratification category.  ncpus——threads number

data <- read.csv("result.csv")         
sde <- apply(data,2,sd)
mean <- apply(data,2,mean)
final <- sde/mean

matrix(final,nrow=18000,ncol=36000)
fwrite(a,"bootstrap_result.csv",col.names=F,sep=",")
#map generation step as described above.

#Validation-PCA
library(tidyverse)
library(ggfortify)
library(ggrepel)
library(sp)

env.meta <- read.csv("$meta.csv",row.names = 1,as.is = TRUE)
#Environmental information of all the sample sites, and every row represents a factor and the index needed to calculate. Be sure that rownames should be similar with names in the input model.


#Calculate PCA result of sampled data
pcaOutput = prcomp(meta, center = TRUE,scale. = TRUE,retx=TRUE)
summary(pcaOutput)

#Chosse axises
pcaOutput$sdev^2
#Choose the axises that contributed more than 80% 

b <- rep(0,648000000)
#Transfer all the data into PCA
pcameta <- scale(env.meta.all, pcaOutput$center, pcaOutput$scale) %*% pcaOutput$rotation #env.all.meta generated by the steps as described in 'Stratified bootstrapping procedure' part.
pcaScores<- pcaOutput$x
for (i in 2:6){
  j <- 1
  while (j <= i-1){
    pc_chull = chull(pcaScores[,c(j,i)])#chull calculate
    pc_Coordinates = pcaScores[pc_chull,c(j,i)]
    a <- point.in.polygon(pcameta[,j], pcameta[,i], pc_Coordinates[,1], pc_Coordinates[,2], mode.checked=FALSE)#classify if the dot is in the chull
    a[(a>1)] <- 1
    b <- b+a
    j <- j+1
  }
}

matrix(b,nrow=18000,ncol=36000)
fwrite(a,"pca_result.csv",col.names=F,sep=",")
#Map generation step as described above.

#Extended data Figure8a
library(ggplot2)
library(ggcorrplot)
library(corrplot)
library(tidyverse)

env.var <- readxl::read_xlsx("Supplementary_Table4.xlsx", col_types = "numeric")
env.var <- as.data.frame(env.var[,5:76])
tree <- hclustvar(X.quanti = env.var)
tree.cut14 <- cutreevar(tree, k = 14)
melt(sort(tree.cut14$cluster))


tree.phy <- as.phylo(tree)

g.tree <- ggtree(tree.phy, ladderize = TRUE)

label.df <- g.tree$data[1:72, c(4,6,7)]

sel <- rep("italic", 72)
names(sel) <- g.tree$data$label[1:72]
sel[label.df$label %in% tolower(row.names(drv))] <- "bold"

label.df <- label.df%>%
        mutate(cluster = sort(tree.cut14$cluster)[g.tree$data$label[1:72]]) %>%
        mutate(select = sel)

label.df$label[8] <- "soil_pH"

ggtree(tree.phy, ladderize = TRUE) +
        geom_point(data = label.df[label.df$select == "bold",], 
                   color = "darkred",
                   aes(x = 4.63, y, )) +
        geom_label(data = label.df[label.df$select == "bold",], 
                   color = "darkred",
                   size = 3,
                   nudge_x = .1,fontface = "bold",
                   hjust = 0,fill = "white",
                   aes(x = 4.6, y, label = label)) +
        geom_text(
                aes(x = 4.6,
                        y,
                        label = label,
                        color = factor(cluster)
                ),
                fontface = "italic",size = 2.5,
                data = label.df[label.df$select != "bold",],
                inherit.aes = FALSE,
                nudge_x = .1,
                hjust = 0
        ) +
        xlim(0, 6) +
        guides(color = "none", size = "none") +
        geom_vline(xintercept = 4.05,linetype = 2,color = "darkred")

ggsave("variance_select.pdf",height = 9,width = 6)

#Extended data Figure8d&e
# load files

drv <- read.csv("DRV.csv",header = TRUE,row.names = 1)
drv[is.na(drv)] <- 0

drv[,3:16] <- -(drv[,3:16] + t(drv[,3:16]))

row.names(drv) <- tolower(row.names(drv))
colnames(drv) <- tolower(colnames(drv))

ggcorrplot(drv[,16:1],
           #method = "circle",
           pch.cex = .6,
           type = "full",
           colors = c(pal_npg(alpha = .8)(5)[2],"white",pal_npg(alpha = .8)(5)[5]),
           show.legend = FALSE,
           show.diag = FALSE)

ggsave("DRV_sq.pdf", height = 5,width = 5)

richness <- read.csv("richness.csv",header = TRUE,row.names = 1)
richness[is.na(richness)] <- 0
richness[,3:16] <- - richness[,3:16] - t(richness[,3:16])

row.names(richness) <- tolower(row.names(richness))
colnames(richness) <- tolower(colnames(richness))

ggcorrplot(richness[,16:1],
           #method = "circle",
           type = "full", pch.cex = .6,
           colors = c(pal_npg(alpha = .8)(5)[2],"white",pal_npg(alpha = .8)(5)[5]),
           show.legend = FALSE,
           show.diag = FALSE)

ggsave("richness_sq.pdf", height = 5,width = 5)




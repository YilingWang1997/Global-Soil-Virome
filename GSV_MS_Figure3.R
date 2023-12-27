library(reshape2)
library(Hmisc)
library(data.table)
library(ggcorrplot)
library(igraph)
library(shapefiles)
library(ggbeeswarm)
library(tidyverse)
library(ggraph)
library(ggsci)
library(colorspace)
library(ggpubr)
library(geosphere)
library(metacom)
library(data.table)

#fig3a meta-community
pal <- c("#BE564A",
         "#F9B483",
         "#B4C2CF",
         "#78BAD2",
         "#576289",
         "#6C958F",
         "#6B9725",
         "#BED992",
         "#FEBE44",
         "#FDEACD",
         "#7C7C7C",
         "#C9C9C9"
)

metadata <- read.csv("../Metadata/metadata_20210107.csv",)
rownames(metadata) <- metadata$Sample.id
vc.table <- fread(file = "otutable_family.csv",stringsAsFactors = FALSE)
vc.01mat <- vc.table[,-1]
vc.01mat[vc.table[,-1] != 0 ] <- 1 
vc.dom <- as.matrix(vc.01mat[rowSums(vc.01mat) > 10, ])
vc.ord <- OrderMatrix(vc.dom)

Turnover(vc.ord,method='r1',sims=100, scores=1, binary=TRUE)
#name          stat
#1    turnover  3.653519e+10
#2           z -6.685896e+02
#3           p  0.000000e+00
#4     simMean  7.831441e+10
#5 simVariance  6.248858e+07
#6 method = r1            NA
BoundaryClump(vc.ord)
#   name       stat
#1 index   11.71937
#2     p    0.00000
#3    df 2190.00000
Coherence(vc.ord,method = "r1", sims=10,scores = 1)
#name         stat
#1      embAbs 2506530.0000
#2           z    -547.1011
#3           p       0.0000
#4     simMean 3866724.3000
#5 simVariance    2486.1846
#6 method = r1           NA

  
sample.ord <- colnames(vc.ord)
continent <- metadata[sample.ord ,7]
land <- metadata[sample.ord ,6]

row.names(vc.ord) <- 1:dim(vc.dom)[1]
colnames(vc.ord) <- 1:dim(vc.dom)[2]
vc.heat.df <- gather(data.frame(vc.ord)) %>%
  mutate(vc = rep(1:dim(vc.dom)[1],  dim(vc.dom)[2])) %>%
  mutate(sample = rep(1:dim(vc.dom)[2], each = dim(vc.dom)[1])) %>%
  mutate(land = rep(land, each = dim(vc.dom)[1])) %>%
  mutate(cont = rep(continent, each = dim(vc.dom)[1])) %>%
  filter(value != 0)

ggplot(vc.heat.df, aes(y = vc, x = sample)) +
  geom_point(shape = 15,
             alpha= .5,
             size = .01,
             color = pal[1]) +
  theme_light() + xlab("1,859 samples") + ylab("17,700 family-level vOTUs") +
  theme(aspect.ratio = .618,
        axis.title = element_text(face = "bold",size = 10))

ggsave(filename = "metacom.pdf", width = 100,height = 70,units = "mm")

#fig 3b beta & environment factor

library(vegan)

sample <- read.csv("1799biomcont.csv")
bc <- read.csv("bc_dist_1214.csv",row.names = 1)
meta <- read.csv("meta_scale_new_all_no0_1845_1208.csv",row.names = 1)
meta <- meta[meta$Sample.id%in%sample$Sample.id,]


infos <- NULL
for (i in c(3:86,90:91)){
dist<- dist(meta[,i])
result <- mantel(bc,dist,na.rm=TRUE)
statistic <- result$statistic
signif <- result$signif
info <- data.frame(colnames(meta)[i],statistic,signif)
infos <- rbind(infos,info)
}

infoall <- infos
infoall <- infoall[infoall$statistic>=0.1665,]

efbc <- read.csv("bc_effectsize_esp.tsv",sep="\t")
efbc2 <- read.csv("bc_effectsize.tsv",sep="\t")
efbc <- rbind(efbc,efbc2)
efbc <- efbc[efbc$effect_size>=0.25,]

infoall <- infoall[infoall$colnames.meta..i. %in% efbc$column,]

for (i in infoall$colnames.meta..i.){
  dist<- dist(meta[,i])
  result <- mantel(bc,dist,na.rm=TRUE)
  print(i)
  print(result)
}

mantel <- read.csv("mantelresult.csv")
mantel <- mantel[mantel$statistic>=0.1665,]
mantel$start <- 0.16
ggplot(data = mantel,aes(x=reorder(env,statistic),xend = reorder(env,statistic) , y=start,yend = statistic)) + 
  #geom_bar(stat = 'identity',fill="#9cc6dd",width=0.7)+
  geom_segment(aes(color="#9cc6dd"),size=8,show.legend = FALSE)+
  coord_flip()+xlab("")+ylab("Mantel r")+
  scale_color_manual(values = c("#9cc6dd"))+
  geom_errorbar(aes(ymin=statistic,ymax=statistic+quantiles,width=0.2))+
  #geom_hline(yintercept = 0.1,linetype="dotted")+
  ylim(0.16,0.19)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.text.x=element_text(angle=60,size=10,hjust=1,face = "bold"),axis.text.y=element_text(size=9,face = "bold"))

ggsave("mantel.pdf",width=1.8,height=5)

#fig 3c alpha & environment factor
library(data.table)
library(ggcorrplot)
library(colorspace)
library(Hmisc)

meta <- read.csv("meta_scale_new_all_no0_1799_1214.csv")

efalpha1 <- read.csv("espece.csv.deciles_shannonspecies.tsv",sep = "\t")
efalpha2 <- read.csv("meta_scale_new_all_no0_1845_deciles_orderbybc_shannonspecies_nameedit.tsv",sep = "\t")
efalpha <- rbind(efalpha1,efalpha2)
efalpha <- efalpha[efalpha$effect_size >= 0.4,]
meta <- meta[,c("shannon_species",efalpha$column[c(2:15,17:22)])]


cor.shannon <- NULL
for(i in c(2:21)){
  corv <- rcorr(meta$shannon_species, meta[,i], type="spearman")
  res <- c(corv$r[1,2],corv$P[1,2])
  cor.shannon <- rbind(cor.shannon,round(res, 3))
}

env.cor <- cor(na.omit(meta[,c(2:21)], method = "spearman"))
env.cor.p <- cor_pmat(na.omit(meta[,c(2:21)], method = "spearman"))

cor.df <- data.frame(value = env.cor[lower.tri(env.cor)],
                     x = c(1,1:2,1:3,1:4,1:5,1:6,1:7,1:8,1:9,1:10,1:11,1:12,1:13,1:14,1:15,1:16,1:17,1:18,1:19),
                     y = c(19,
                           rep(18,2),
                           rep(17,3),
                           rep(16,4),
                           rep(15,5),
                           rep(14,6),
                           rep(13,7),
                           rep(12,8),
                           rep(11,9),
                           rep(10,10),
                           rep(9,11),
                           rep(8,12),
                           rep(7,13),
                           rep(6,14),
                           rep(5,15),
                           rep(4,16),
                           rep(3,17),
                           rep(2,18),
                           rep(1,19)
                           ))

cornet <- data.frame(x = c(1:20, 1:20),
                     y = c(20:1, 20:1),
                     xend = rep(18,40),
                     yend = rep(18,40),
                     value = c(cor.shannon[,1],cor.shannon[,1]),
                     group = rep(c("Shannon's H’"),40)
)

label.df <- data.frame(x = c(1:20,  18),
                       y = c(20:1, 18),
                       env = c("ESP","Reflectance_6","Sand (0-5cm)","Sand (5-15cm)","Longitude","MTWM","Silt (0-5cm)",
                               "Latitude","TS","TAR","MTDQ","Reflectance_4","RDQ","RWQ",
                               "PDQ","Silt (5-15cm)","EVI Heterogeneity","OCD (5-15cm)","Reflectance_2",
                               "MTCQ","Shannon's H`"))

label.df <- data.frame(x = c(1:20,  18),
                       y = c(20:1, 18),
                       env = c("Exchangeable sodium percentage","Reflectance_6","Sand (0-5cm)","Sand (5-15cm)","Longitude","Max temperature of warmest month","Silt (0-5cm)",
                               "Latitude","Temperature seasonality","Temperature annual range","Mean temperature of driest quarter","Reflectance_4","Radiation of driest quarter","Radiation of warmest quarter",
                               "Precipitation of driest quarter","Silt (5-15cm)","EVI Heterogeneity","Organic carbon density (5-15cm)","Reflectance_2",
                               "Mean temperature of coldest quarter","Shannon's H`"))

label.df$color <- c(1,3,1,1,2,4,1,2,4,4,4,3,4,4,4,1,1,1,3,4,1)

ggplot(na.omit(cor.df), aes(x = x, y = y, alpha = abs(value))) +
  geom_point(aes(color = value),
             size = 7,
             shape = 16) +
  theme_void() +
  scale_color_continuous_divergingx(palette = "RdBu",
                                    name = "Spearman's rho ") +
  scale_size_continuous(range = c(0.1, 3),
                        limits = c(0.052, .5),
                        breaks = c(.1, .2, .3,.4),
                        name = "Mantel r") +
  theme(aspect.ratio = 11/15,
        legend.position = "right") +
  geom_curve(data = cornet, aes(
    x = x,
    y = y,
    xend = xend,
    yend = yend,
    size = abs(value)
  ),inherit.aes = FALSE,
  curvature = 0.2,
  color = "grey",
  alpha = .6,
  angle = 20,
  lineend = "round") +
  xlim(-1,20) +
  geom_label(data = label.df, 
             aes(
               x = x,
               y = y,
               label = env,
               fill = color), 
             nudge_x = -.4,
             inherit.aes = FALSE, 
             fontface = "bold", 
             hjust = 0,
             size = 3,fill = grey(.9,alpha = .5)) +
  guides(alpha = FALSE)


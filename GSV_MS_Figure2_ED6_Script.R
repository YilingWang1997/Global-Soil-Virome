# package loading
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
library(vegan)
library(dplyr)
library(ggplot2)

# color pallate
pal <- c("#C75E48","#FDEACD","#FEB27A","#445481","#5D9089","#Feb928","#747474","#6BB9D2")



# Figure 2
## alpha-diversity (figure 2a)

alphabetical <- c("Agricultural Land","Artificial Surfaces","Bare Land","Forest","Grassland","Shrubland","Tundra","Wetland")
orderalpha <- c("Agricultural Land","Artificial Surfaces","Bare Land","Wetland","Grassland","Tundra","Forest","Shrubland")
metadata <-  read.csv("meta_scale_new_all_no0_1799_1213.csv")
meta <- read.csv("Supplementary_Table1_1213.csv")
metadata <- merge(metadata,meta,by="Sample.id")
pal_alpha<- pal[pmatch(orderalpha,alphabetical)]

alpha.div_sta <-
  data.frame( sample = metadata$Sample.id,
              cont = metadata$Continent,
              land = metadata$Biome,
              shannon_spe = metadata$shannon_species,
              shannon_ge = metadata$shannon_genus,
              shannon_fa = metadata$shannon_family,
              seqdepth = metadata$Size.bp.
  )

library(agricolae)
aov.shannon<-aov(shannon_spe ~ land,data = alpha.div_sta)
summary(aov.shannon)
result<-LSD.test(aov.shannon,"land",p.adj = "bonferroni")
result$groups
aov.shannon<-aov(shannon_norm ~ land,data = alpha.div_sta)
summary(aov.shannon)
result<-LSD.test(aov.shannon,"land",p.adj = "bonferroni")
result$groups

library(agricolae)
aov.shannon<-aov(shannon_ge ~ land,data = alpha.div_sta)
summary(aov.shannon)
result<-LSD.test(aov.shannon,"land",p.adj = "bonferroni")
result$groups

library(agricolae)
aov.shannon<-aov(shannon_fa ~ land,data = alpha.div_sta)
summary(aov.shannon)
result<-LSD.test(aov.shannon,"land",p.adj = "bonferroni")
result$groups

alpha.div_spe <-
  data.frame(land = metadata$Biome,
             shannon_spe = metadata$shannon_species,
             level = "Species Level"
  )

colnames(alpha.div_spe) <- c("land","shannon","level")

alpha.div_ge <-
  data.frame(land = metadata$Biome,
             shannon_ge = metadata$shannon_genus,
             level = "Genus Level"
  )

colnames(alpha.div_ge) <- c("land","shannon","level")

alpha.div_fa <-
  data.frame(land = metadata$Biome,
             shannon_fa = metadata$shannon_family,
             level = "Family Level"
  )

colnames(alpha.div_fa) <- c("land","shannon","level")

alpha.div <- rbind(alpha.div_spe,alpha.div_ge,alpha.div_fa)

alpha.div$land <-
  factor(
    alpha.div$land,
    levels = c(
      "Agricultural Land",
      "Artificial Surfaces",
      "Bare Land",
      "Wetland",
      "Grassland",
      "Tundra",
      "Forest",
      "Shrubland"
    ),
    labels = c(
      "Agricultural Land",
      "Artificial Surfaces",
      "Bare Land",
      "Wetland",
      "Grassland",
      "Tundra",
      "Forest",
      "Shrubland"
    )
  )

alpha.div$level <- 
  factor(
    alpha.div$level,
    levels = c(
      "Species Level",
      "Genus Level",
      "Family Level"
    ),
    labels = c(
      "Species Level",
      "Genus Level",
      "Family Level"
    ))


ggplot(na.omit(alpha.div),
       aes(x = land, y = shannon, color = land)) +
  geom_boxplot(alpha = 1, outlier.colour = NA,
               draw_quantiles = .5
  )+geom_beeswarm(size = .5, alpha = .2)+
  facet_grid(~level,
             switch = "y",
             scales = "free",
             space = "free_x") +
  theme_bw() + ylab("Shannon") + xlab("") +
  #scale_y_log10()+
  scale_color_manual(values = pal_alpha)+
  # scale_y_continuous(breaks = seq(0,5,1),
  #                   limits = c(1,4),
  #labels = expression(10^0,10^1,10^2,10^3,10^4,10^5)
  #                  )+
  theme(legend.position = "none", 
        
        panel.grid = element_blank(),
        #aspect.ratio = 8,
        axis.text.x = element_text(angle = -65,
                                   vjust = 0.5, 
                                   hjust = 0,
                                   size = 9),
        strip.text.x = element_text(size=9),
        strip.text.y = element_text(vjust = 0),
        strip.placement = "outside",
        strip.background.y = element_rect(fill = NA, color = NA))

### figure 2b & Extended figure 6

diversity<- read.csv("nonpareil_diversity.csv",row.names = 1)
diversity$index <- gsub("_R1","",diversity$index)
diversity$index <- gsub("_paired","",diversity$index)
diversity$index <- gsub("-1a","",diversity$index)
diversity$index <- gsub("-2a","",diversity$index)
diversity$index <- gsub("_paired_1","",diversity$index)
diversity$index <- gsub("_1","",diversity$index)
diversity$index <- gsub(".clean","",diversity$index)

viraldiversity <- read.csv("meta_scale_new_all_no0_1799_1214.csv")
viraldiversity <- viraldiversity[,c("Sample.id","Biome","shannon_species")]
viraldiversity <- viraldiversity[pmatch(diversity$index,viraldiversity$Sample.id),]

data <- cbind(diversity,viraldiversity)
data <- data[data$result.6.>=22 & data$result.6. <=25,]
plot(data$shannon_species,data$result.6.)

rare <- read.csv("shannon_refraction_2022.csv",row.names = 1)
rownames(rare) <- gsub("_paired__merge.csv","",rownames(rare))
rownames(rare) <- gsub("_merge.csv","",rownames(rare))
rownames(rare) <- gsub("rhiz","rice",rownames(rare))
rownames(rare) <- gsub("rhio","rice",rownames(rare))
rare$Sample.id <- rownames(rare)

data <- merge(data,rare,by="Sample.id")

data$Biome <-
  factor(
    data$Biome,
    levels = c(
      "Agricultural Land",
      "Artificial Surfaces",
      "Bare Land",
      "Wetland",
      "Grassland",
      "Tundra",
      "Forest",
      "Shrubland"
    ),
    labels = c(
      "Agricultural Land",
      "Artificial Surfaces",
      "Bare Land",
      "Wetland",
      "Grassland",
      "Tundra",
      "Forest",
      "Shrubland"
    )
  )

alphabetical <- c("Agricultural Land","Artificial Surfaces","Bare Land","Forest","Grassland","Shrubland","Tundra","Wetland")
pal <- c("#C75E48","#FDEACD","#FEB27A","#445481","#5D9089","#Feb928","#747474","#6BB9D2")
orderalpha <- c("Agricultural Land","Artificial Surfaces","Bare Land","Wetland","Grassland","Tundra","Forest","Shrubland")
pal_alpha<- pal[pmatch(orderalpha,alphabetical)]

ggplot(data,aes(x=result.6.,y=shannon_species,color=Biome))+geom_point(size=1.5,shape=16,alpha=0.5)+scale_color_manual(values=pal_alpha)+
  theme_bw()+xlab("Metagenome Nd \n(Microbial community)")+ylab("Shannon\n(Viral community, Species level)")+
  theme(legend.position = "none",
        panel.grid = element_blank())

ggplot(data,aes(x=result.6.,y=shannon_species,color=Biome))+geom_point()+scale_color_npg()+geom_smooth(aes(group=Biome),method="lm")
ggplot(data,aes(x=result.6.,y=x,color=Biome))+geom_point()+scale_color_npg()+geom_smooth(aes(group=Biome),method="lm")

lmresult <- lm(data$shannon_species~data$x)

cor.test(data$shannon_species,data$x,method="spearman")
for (i in 1:8){
target <- orderalpha[i]
data_sel <- data[data$Biome==target,]
print(cor.test(data_sel$shannon_species,data_sel$x,method="spearman"))}

summary(lmresult)
ggplot(data[data$x!=0,],aes(x=shannon_species,y=x,))+geom_point(color="darkgrey",size=2.5,shape=16,alpha=0.2)+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank())


ggplot(data[data$x!=0,],aes(x=x,y=shannon_species,color=Biome))+geom_point()+scale_color_manual(values=pal_alpha)+geom_smooth(aes(group=Biome),method="lm")+
  theme_bw()+xlab("Shannon\n(Subsampled reads, Species level)")+ylab("Shannon\n(Total reads, Species level)")+
  theme(panel.grid = element_blank())


target <- orderalpha[1]
data_sel <- data[data$Biome==target,]
plot1 <- ggplot(data_sel[data_sel$x!=0,],aes(x=x,y=shannon_species))+geom_point(color=pal_alpha[1])+geom_smooth(method="lm",color=pal_alpha[1])+
         theme_bw()+xlab("Shannon\n(Subsampled reads, Species level)")+ylab("Shannon\n(Total reads, Species level)")+
         theme(panel.grid = element_blank())
target <- orderalpha[2]
data_sel <- data[data$Biome==target,]
plot2 <- ggplot(data_sel[data_sel$x!=0,],aes(x=x,y=shannon_species))+geom_point(color=pal_alpha[2])+geom_smooth(method="lm",color=pal_alpha[2])+
  theme_bw()+xlab("Shannon\n(Subsampled reads, Species level)")+ylab("Shannon\n(Total reads, Species level)")+
  theme(panel.grid = element_blank())
target <- orderalpha[3]
data_sel <- data[data$Biome==target,]
plot3 <- ggplot(data_sel[data_sel$x!=0,],aes(x=x,y=shannon_species))+geom_point(color=pal_alpha[3])+geom_smooth(method="lm",color=pal_alpha[3])+
  theme_bw()+xlab("Shannon\n(Subsampled reads, Species level)")+ylab("Shannon\n(Total reads, Species level)")+
  theme(panel.grid = element_blank())
target <- orderalpha[4]
data_sel <- data[data$Biome==target,]
plot4 <- ggplot(data_sel[data_sel$x!=0,],aes(x=x,y=shannon_species))+geom_point(color=pal_alpha[4])+geom_smooth(method="lm",color=pal_alpha[4])+
  theme_bw()+xlab("Shannon\n(Subsampled reads, Species level)")+ylab("Shannon\n(Total reads, Species level)")+
  theme(panel.grid = element_blank())
target <- orderalpha[5]
data_sel <- data[data$Biome==target,]
plot5 <- ggplot(data_sel[data_sel$x!=0,],aes(x=x,y=shannon_species))+geom_point(color=pal_alpha[5])+geom_smooth(method="lm",color=pal_alpha[5])+
  theme_bw()+xlab("Shannon\n(Subsampled reads, Species level)")+ylab("Shannon\n(Total reads, Species level)")+
  theme(panel.grid = element_blank())
target <- orderalpha[6]
data_sel <- data[data$Biome==target,]
plot6 <- ggplot(data_sel[data_sel$x!=0,],aes(x=x,y=shannon_species))+geom_point(color=pal_alpha[6])+geom_smooth(method="lm",color=pal_alpha[6])+
  theme_bw()+xlab("Shannon\n(Subsampled reads, Species level)")+ylab("Shannon\n(Total reads, Species level)")+
  theme(panel.grid = element_blank())
target <- orderalpha[7]
data_sel <- data[data$Biome==target,]
plot7 <- ggplot(data_sel[data_sel$x!=0,],aes(x=x,y=shannon_species))+geom_point(color=pal_alpha[7])+geom_smooth(method="lm",color=pal_alpha[7])+
  theme_bw()+xlab("Shannon\n(Subsampled reads, Species level)")+ylab("Shannon\n(Total reads, Species level)")+
  theme(panel.grid = element_blank())
target <- orderalpha[8]
data_sel <- data[data$Biome==target,]
plot8 <- ggplot(data_sel[data_sel$x!=0,],aes(x=x,y=shannon_species))+geom_point(color=pal_alpha[8])+geom_smooth(method="lm",color=pal_alpha[8])+
  theme_bw()+xlab("Shannon\n(Subsampled reads, Species level)")+ylab("Shannon\n(Total reads, Species level)")+
  theme(panel.grid = element_blank())

library(cowplot)
#plot_grid(plot1,plot2,plot3,plot4,plot5,plot6,plot7,plot8,nrow=2)
plot_grid(plot1,plot2,plot3,plot4,nrow=1)
plot_grid(plot5,plot6,plot7,plot8,nrow=1)

### NMDS (figure 2c) 

meta <- read.csv("Supplementary_Table1_1213.csv")
bc_dist <- read.csv("bc_dist_1214.csv",row.names = 1)
dist.1 <- na.exclude(bc_dist)
dist.2 <- na.exclude(t(dist.1))
nmds <- metaMDS(dist.2, k = 2,)
nmdsforg <- nmds$points
nmds.df <- nmdsforg %>%
  mutate(land = eco$Biome) %>%
  mutate(cont = eco$Continent) %>%
  arrange(land) 

eco <- read.csv("1799biomcont.csv",row.names = 2)

eco <- eco[rownames(nmdsforg),]

nmds.df.remain <- nmds.df[nmds.df$MDS1<=-0.0003&nmds.df$MDS1>=-0.0018&nmds.df$MDS2<=0.0011&nmds.df$MDS2>=-0.0001,]
write.csv(nmds.df.remain,"nmds.df.remain.csv")

nmds.df <- read.csv("nmds.df.remain.csv")
pal<- c("#C75E48","#FDEACD","#FEB27A","#445481","#5D9089","#Feb928","#747474","#6BB9D2")

eco.seq <- factor(nmds.df$land, 
                  levels = c("Agricultural Land","Forest","Wetland","Grassland",       
                             "Tundra" ,"Artificial Surfaces","Bare Land","Shrubland"))
nmds.df <- nmds.df[order(eco.seq),]

ggplot(na.omit(nmds.df), aes(
  x = MDS1,
  y = MDS2,
  color = land
)) +
  geom_point(size = .8,
             shape = 16,
             alpha = .5) +
  xlim(-0.0018,-0.0003) + ylim(-0.0001,0.0011) +
  theme_bw() + xlab("NMDS 1") + ylab("NMDS 2")+
  scale_color_manual(values = pal) +
  guides(colour = guide_legend(title = "Biome",
                               override.aes = list(alpha = 1))) +
  #stat_ellipse(level = .6) +
  #facet_wrap( ~ sep, scales = "free", ncol = 2) +
  theme(aspect.ratio = 1,
        legend.position = "right",#c(.3,.7),
        axis.title = element_text(size = 8),
        legend.key.size = unit(x = 2, "mm"),
        legend.background = element_blank(), 
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7,face = "bold"),
        panel.grid = element_blank())
ggsave("nmdsland.pdf")

cont.seq <- factor(nmds.df$cont, 
                   levels = c("North America","Asia","Australia","Europe","South America","Africa"))

nmds.df <- nmds.df[order(cont.seq),]

ggplot(na.omit(nmds.df), aes(
  x = MDS1,
  y = MDS2,
  color = cont
)) +
  geom_point(size = 1.5,
             shape = 16,
             alpha = .5) +
  xlim(-0.0018,-0.0003) + ylim(-0.0001,0.0011) +
  theme_bw() + xlab("NMDS 1") + ylab("NMDS 2")+
  scale_color_manual(values = pal)+
  guides(colour = guide_legend(title = "Ecosystem",nrow = 4,
                               override.aes = list(alpha = 1))) +
  #stat_ellipse(level = .6) +
  #facet_wrap( ~ sep, scales = "free", ncol = 2) +
  theme(aspect.ratio = 1,
        legend.position = "bottom",#c(.3,.7),
        axis.title = element_text(size = 8),
        legend.key.size = unit(x = 2, "mm"),
        legend.background = element_blank(), 
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7,face = "bold"),
        panel.grid = element_blank())

ggsave("nmdscont.pdf")


library(vegan)
otu <- read.csv("OTU_GSV_average_2022_1859.csv",row.names = 1)
otu <- t(otu)
eco <- read.csv("1799biomcont.csv")
rownames(eco) <- eco$Sample.id

remain  <- read.csv("nmds.df.remain.csv")
eco <- eco[remain$X,]
otu <- otu[rownames(eco),]

Biome <- eco$Biome
adonis(otu~Biome, otu, permutations = 999, distance = 'bray')

Continent <- eco$Continent
adonis(otu~Continent, otu, permutations = 999, distance = 'bray')




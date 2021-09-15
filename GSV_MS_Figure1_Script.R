# package loading
library(rgdal)
library(ggplot2)
library(tidyr)
library(ggalluvial)
library(stringr)
library(reshape2)
library(ggforce)
library(UpSetR)
library(rPython)

#Fig1a
sample <- read.csv("ST1_Sampleinfo_GSV.csv")
sample2 <- unite(sample, "Biome_Continent", Biome, Continent, sep = ",", remove = FALSE)
summary <- table(sample2$Biome_Continent)
summary <- data.frame(summary)
library(stringr)
summary2 <- str_split_fixed(summary$Var1, pattern = ",",2)
summary2 <- data.frame(summary2)
summary3 <- cbind(summary2,summary$Freq)
colnames(summary3) <- c("Biome","Continent","Value")

gap <- data.frame(Biome=c("g1","g2","g3","g4","g5","g6"),Continent=c("gap1","gap2","gap3","gap4","gap5","gap6"),Value=rep(30,6))
summary4 <- rbind(summary3,gap)
  
summary4$Biome <- factor(summary4$Biome,levels = c("Wetland","g6","Tundra","g5","Shrubland","g4","Grassland","g3","Forest","g2","Bare Land","g1","Agricultural Land"))

summary4$Continent <- factor(summary4$Continent,levels = c("gap6","South America","gap5","North America","gap4","Europe","gap3","Australia","gap2","Asia","gap1","Africa"))

cbPalette <- c("#747474","#FFFFFF00","#FADE9A","#FFFFFF00","#5D9089","#FFFFFF00","#445481","#FFFFFF00","#6BB9D2","#FFFFFF00","#FEB27A","#FFFFFF00","#C75E48")

ggplot(data = summary4,
       aes(y = Value,axis1 = Continent, axis2 = Biome)) +
  #scale_x_discrete(limits = c("continent", "ecosystem"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = Biome),width = 0,knot.pos = 1/6,alpha=1) +
  #geom_stratum() + 
  geom_text(stat = "stratum", label.strata = TRUE) +
  theme_minimal()+
  scale_fill_manual(values = cbPalette)+
  #scale_x_continuous(breaks = 1:2)+
  coord_flip()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

#Fig1b
# read shapefile
wmap <- readOGR(dsn="ne_110m_land", layer="ne_110m_land")

# convert to dataframe
wmap_df <- fortify(wmap)

# create a blank ggplot theme
theme_opts <-list(theme(panel.grid.minor = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.background = element_blank(),
                        plot.background = element_rect(fill="white"),
                        panel.border = element_blank(),
                        axis.line = element_blank(),
                        axis.text.x = element_blank(),
                        axis.text.y = element_blank(),
                        axis.ticks = element_blank(),
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank(),
                        plot.title = element_text(size=22,hjust = .5)))

# plot map
sample <- read.csv("ST1_Sampleinfo_GSV.csv")
sample.public <- sample[(sample$Type=="Public"),]
sample.public <- data.frame(lat=sample.public$Lat,lon=sample.public$Lon,size=sample.public$Size,biome=sample.public$Biome)
sample.inhouse <- sample[(sample$Type=="In-house"),]
sample.inhouse <- data.frame(lat=sample.inhouse$Lat,lon=sample.inhouse$Lon,size=sample.inhouse$Size,biome=sample.inhouse$Biome)

ggplot(wmap_df, aes(long,lat, group=group)) + 
  geom_polygon(fill = "gray90") + 
  coord_equal() + 
  theme_opts+
  #guides(fill=guide_legend(title="Abundance"))+
  geom_point(data=sample.public,
             aes(lon, lat,color=biome,size=size),
             alpha=.8,
             inherit.aes = FALSE)+
  geom_point(data=sample.inhouse,
             aes(lon, lat,size=size,color=biome),
             shape=17,
             alpha=.8,
             inherit.aes = FALSE)+
  scale_size_continuous(range = c(1,4))+
  scale_color_manual(values=c("#C17360","#EAB88D","#7CABBE","#525E80","#66837F","#DAC798","#7E7E7E"))+
  theme(legend.position = "none")

#Fig1d
expressionInput <- c(GSV = 37841, PIGEON= 34424,GOV2=17912,GPD=18009,IMGVR_Soil = 27539,GVD=512,refseq=944,`PIGEON&GSV` = 9239,`GOV2&GSV` = 504,`GPD&GSV` = 270,`IMGVR_Soil&GSV` = 22541, `GVD&GSV` = 185,`refseq&GSV`=60)
upset(fromExpression(expressionInput),nsets = 7)

#Fig1e
gtdbtaxa <- read.csv("gtdb_taxanomy_rep_31910.csv")#GTDB taxanomy information
result <- read.csv("GSV_hostlinkage_summary_taxa_1873.csv")#Host information summary
gtdbtaxa$hostNCBIName <- gsub("_genomic.fna","",gtdbtaxa$hostNCBIName)
result_taxa <- merge(result,gtdbtaxa,by.x="file",by.y="hostNCBIName")
phyluminfo <- data.frame(table(result_taxa$hostPhylum))
phyluminfo <- phyluminfo[(phyluminfo$Freq>=1),]
write.csv(phyluminfo,"phyluminfo_1873.csv")
classinfo <- data.frame(table(result_taxa$hostClass))
classinfo <- classinfo [(classinfo $Freq>=1),]

phylumclass<-data.frame(table(result_taxa[,9:10]))
phylumclass<-phylumclass[(phylumclass$Freq>=1),]

phyluminfo <- phyluminfo[order(phyluminfo$Freq,decreasing=T), ]

phylumclass <- merge(phylumclass,phyluminfo,by.x="hostPhylum",by.y="Var1")

phylumclass <- phylumclass[order(-phylumclass$Freq.y,-phylumclass$Freq.x),]

phylumclass$hostClass <- as.character(phylumclass$hostClass)
phylumclass$hostClass[(phylumclass$Freq.x<=20)] <- "Others"
library(dplyr)
phylumclass$a <- paste(phylumclass$hostPhylum,phylumclass$hostClass)
phylumclass_new <- phylumclass %>% 
  dplyr::group_by(a) %>% 
  dplyr::summarise(Freq.x=sum(Freq.x))
library(tidyr)
phylumclass_new <- separate(phylumclass_new,a,c("phylum","class"),sep=" ",remove=T)
phylumclass_new <- merge(phylumclass_new,phyluminfo,by.x="phylum",by.y="Var1")

result_taxa$hostSuperkingdom <- as.character(result_taxa$hostSuperkingdom)
for (i in 1:dim(phylumclass_new)[1]){
  name <- phylumclass_new$phylum[i]
  pos <-  grep(name,result_taxa$hostPhylum)
  phylumclass_new$kingdom[i] <- result_taxa$hostSuperkingdom[pos[1]]
}

phylumclass_new <- phylumclass_new[order(phylumclass_new$kingdom,-phylumclass_new$Freq,-phylumclass_new$Freq.x),]

write.csv(phylumclass_new,"phylumclass_1873.csv")

phylumclass <- read.csv("phylumclass_1873_me.csv")

phylumclass$ymax <- phylumclass$Freq.x

for (i in 1:dim(phylumclass)[1]){
  phylumclass$ymax[i] <- sum(phylumclass$Freq.x[1:i])
}

phylumclass$ymin<- c("0",  phylumclass$ymax[1:(dim(phylumclass)[1]-1)])
phylumclass$ymin<-round(as.numeric(phylumclass$ymin),3)

compute_angle = function(perc){
  angle = -1
  if(perc < 0.25) # 1st q [90,0]
    angle = 90 - (perc/0.25) * 90
  else if(perc < 0.5) # 2nd q [0, -90]
    angle = (perc-0.25) / 0.25 * -90
  else if(perc < 0.75) # 3rd q [90, 0]
    angle = 90 - ((perc-0.5) / 0.25 * 90)
  else if(perc < 1.00) # last q [0, -90]
    angle = ((perc -0.75)/0.25) * -90
  # Or even more compact, but less readable
  if(perc < 0.5) # 1st half [90, -90]
    angle = (180 - (perc/0.5) * 180) - 90
  else # 2nd half [90, -90]
    angle = (90 - ((perc - 0.5)/0.5) * 180)
  return(angle)
}

phylumclass_pop = phylumclass %>%
  mutate(running=cumsum(Freq.x), pos=running - Freq.x/2) %>% group_by(1:n()) %>% # to compute row by row
  mutate(angle=compute_angle((running - Freq.x/2) / 35463))

phylumclass_pop$order <- rownames(phylumclass_pop)

phylumclass_pop$order <- as.numeric(phylumclass_pop$order)
order <- tapply(phylumclass_pop$order,phylumclass_pop$phylum,mean)
ymax <- tapply(phylumclass_pop$ymax,phylumclass_pop$phylum,max)
ymin <- tapply(phylumclass_pop$ymin,phylumclass_pop$phylum,min)
Freq.x <- tapply(phylumclass_pop$Freq,phylumclass_pop$phylum,mean)
phylum_pop <- cbind(ymax,ymin,Freq.x,order)
phylum_pop <- data.frame(phylum_pop)
phylum_pop$phylum <- rownames(phylum_pop)
phylum_pop$y <- (phylum_pop$ymax+phylum_pop$ymin)/2

phylum_pop <- phylum_pop[order(phylum_pop$order),]

phylum_pop = phylum_pop %>%
  mutate(running=cumsum(Freq.x), pos=running - Freq.x/2) %>% group_by(1:n()) %>% # to compute row by row
  mutate(angle=compute_angle((running - Freq.x/2) / 35463))

ggplot(phylumclass_pop) + 
  geom_rect(aes(fill=class, ymax=ymax, ymin=ymin, xmax=10, xmin=6)) + 
  geom_rect(aes(fill=phylum, ymax=ymax, ymin=ymin, xmax=6, xmin=2.5)) + 
  geom_rect(aes(fill=kingdom, ymax=ymax, ymin=ymin, xmax=10.4, xmin=10.3)) +
  geom_text(aes(x=8, y=(ymax+ymin)/2, label=class,angle=angle),size=2.5)+
  geom_text(data=phylum_pop,aes(x=4.3, y=y, label=phylum,angle=angle),size=2.5)+
  xlim(c(0,11)) + 
  theme(aspect.ratio=1) + 
  coord_polar(theta="y")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.ticks = element_blank()) +   
  theme(panel.border = element_blank()) +
  theme(axis.text = element_blank(),legend.position = "none")

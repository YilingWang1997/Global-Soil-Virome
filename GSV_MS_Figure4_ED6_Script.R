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
library(geosphere)

# color pallate
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

# data file loading
metadata <- read.csv("update_data/metadata_2021_1873.csv",as.is = TRUE) # metadata
shannon <- read.csv("update_data/shannon_2021_1873.csv",row.names = 1)
shannon <- shannon[metadata$Sample.id,] # arrange
simpson <- read.csv("update_data/simpson_2021_1873.csv",row.names = 1)
simpson <- simpson[metadata$Sample.id,] # arrange
richness <- read.csv("update_data/richness_2021_1873.csv",row.names = 1)
richness <- richness[metadata$Sample.id,] # arrange
shannon.vc <- read.csv("update_data/shannon_vc_2021.csv",row.names = 1) # VC level diversity
shannon.vc <- shannon.vc[metadata$Sample.id,] # arrange 
pcoa <- read.csv("update_data/pcoa_93.csv",row.names = 1) # viruses present in more than 5% 
pcoa.vc <- read.csv("update_data/pcoa_2021_vc.csv",row.names = 1) # VC level
pcoa <- pcoa[metadata$Sample.id,] # arrange 
pcoa.vc <- pcoa.vc[metadata$Sample.id,] # arrange 
nmds <- read.csv("update_data/nmds_93.csv",row.names = 1) # viruses present in more than 5% 
nmds <- nmds[metadata$Sample.id,]
ecodist <- read.csv("update_data/ecodist_2021_1873.csv")[,-1] # Jaccard distance among ecosystems
areadist <- read.csv("update_data/areadist_2021_1873.csv")[,-1] # Jaccard distance among continents
alldist <- read.csv("update_data/all_dist_2021_1873.csv",row.names = 1) # Jaccard distance among all samples 
vc.table <- fread(file = "GSV_vc_averageabundance_10kb_2021.csv",stringsAsFactors = FALSE)

# VC meta-community

vc.01mat <- vc.table[,-1]
vc.01mat[vc.table[,-1] != 0 ] <- 1 
vc.dom <- as.matrix(vc.01mat[rowSums(vc.01mat) > 18, ])

vc.ord <- OrderMatrix(vc.dom)
sample.ord <- colnames(vc.ord)
continent <- metadata[sample.ord ,6]
land <- metadata[sample.ord ,5]

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
        theme_light() + xlab("1873 samples") + ylab("1778 dominant VC") +
        theme(aspect.ratio = .618,
              axis.title = element_text(face = "bold",size = 8))

ggsave("abs_pre.pdf", height = 60, width = 70,units = "mm")

class.df <- data.frame(land = land,
                       cont = continent,
                       x = 1:1873)

class.df$land <- factor(class.df$land, 
                        levels = c("Agricultural land","Forest","Grassland","Wetland",
                                   "Tundra","Shrubland","Bare Land")[7:1])

ggplot(class.df[!is.na(class.df$land), ], 
       aes(x = x, y = land, fill = land)) +
        geom_tile()+
        theme_minimal()+
        scale_fill_manual(values = pal)+
        theme(aspect.ratio = .42, 
              axis.title = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_text(face = "bold",size = 8,color = "black"),
              panel.grid = element_blank(),
              legend.position = "none",
              legend.key.height = unit(1,"mm"))

ggsave(filename = "land_dis.pdf", width = 70,height = 50,units = "mm")

class.df$cont <- factor(class.df$cont,
                        levels = c("Africa", "South America", 
                                   "Europe", "Asia", "North America"))


ggplot(class.df[!is.na(class.df$cont), ], 
       aes(x = x, y = cont, fill = cont)) +
        geom_tile()+
        theme_minimal()+
        scale_fill_manual(values = pal)+
        theme(aspect.ratio = .3, 
              axis.title = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_text(face = "bold",size = 8, color = "black"),
              panel.grid = element_blank(),
              legend.position = "none",
              legend.key.height = unit(1,"mm"))

ggsave(filename = "cont_dis.pdf", width = 70,height = 50,units = "mm")


# global network (figure 3a)

all.mat <- 1 - alldist[metadata$Sample.id, metadata$Sample.id]
all.mat[is.na(all.mat)] <- 0
all.mat[all.mat < .25] <- 0 

g.all <-
        graph_from_adjacency_matrix(
                as.matrix(all.mat),
                mode = "undirected",
                weighted = TRUE,
                diag = FALSE
        )

ly <- data_frame(x = metadata$lon,
                 y = metadata$lat)

library(rgdal)
wmap <- readOGR(dsn="ne_110m_land", layer="ne_110m_land")
wmap_df <- fortify(wmap)

netmap <- ggraph(g.all, layout = ly) +
        geom_polygon(data = wmap_df, 
                     aes(long, lat, 
                         group=group), 
                     fill = grey(.75))+
        geom_node_point(size = .05, 
                        alpha = .3, 
                        color = pal[7])+
        geom_edge_diagonal(#aes(alpha = weight), 
                       color = pal[1],
                       alpha = .5,
                       width = .05,
                       lineend = "round") +
        theme_graph() + 
        theme(legend.position = "none") +
        coord_equal(ylim = c(-55,85))

# matnetplot (Figure 3b)
cor.shannon <- NULL
for(i in c(3,4,7:15)+1){
        corv <- rcorr(shannon, metadata[,i], type="spearman")
        res <- c(corv$r[1,2],corv$P[1,2])
        cor.shannon <- rbind(cor.shannon,round(res, 3))
}

cor.abund <- NULL
for(i in c(3,4,7:15)+1){
        corv <- rcorr(richness, metadata[,i], type="spearman")
        res <- c(corv$r[1,2],corv$P[1,2])
        cor.abund <- rbind(cor.abund, round(res, 3))
}

div.cor.df <- rbind(cor.abund[,1],cor.shannon[,1])
div.p.df <- rbind(cor.abund[,2],cor.shannon[,2])

env.cor <- cor(na.omit(metadata[,c(3,4,7:15)+1]), method = "spearman")
env.cor.p <- cor_pmat(na.omit(metadata[,c(3,4,7:15)+1]), method = "spearman")

cor.df <- data.frame(value = env.cor[lower.tri(env.cor)],
                     x = c(1,1:2,1:3,1:4,1:5,1:6,1:7,1:8,1:9,1:10),
                     y = c(10,
                           rep(9,2),
                           rep(8,3),
                           rep(7,4),
                           rep(6,5),
                           rep(5,6),
                           rep(4,7),
                           rep(3,8),
                           rep(2,9),
                           rep(1,10)))

cornet <- data.frame(x = c(1:11, 1:11),
                     y = c(11:1, 11:1),
                     xend = c(rep(8,11), rep(10,11)),
                     yend = c(rep(10,11), rep(8,11)),
                     value = c(div.cor.df[1,],div.cor.df[1,]),
                     group = rep(c("Richness", "Shannon's H’"),11)
)

direct <- rep(3,22)
direct[cornet$value < 0] <- 4
cornet$direct <- direct

label.df <- data.frame(x = c(1:11, 10, 8),
                       y = c(11:1, 8,10),
                       env = c("Latitude","Longitude","MAP","MAT","Altitude",
                               "UVI","NPP","Soil pH","Moisture",
                               "SOC", "Evaporation","Richness", "Shannon's H`"))

cornet$value[c(div.p.df[1,],div.p.df[2,]) > .05] <- NA

ggplot(na.omit(cor.df), aes(x = x, y = y, alpha = abs(value))) +
        geom_point(aes(color = value),
                   size = 7,
                   shape = 16) +
        theme_void() +
        scale_color_continuous_divergingx(palette = "Fall",
                                          name = "Spearman's rho ") +
        scale_size_continuous(range = c(0.1, 3),
                              limits = c(0.001, .3),
                              breaks = c(.05, .1, .2),
                              name = "Mantel r") +
        theme(aspect.ratio = 12/15,
              legend.position = "right") +
        geom_curve(data = cornet, aes(
                x = x,
                y = y,
                xend = xend,
                yend = yend,
                size = abs(value)
        ),inherit.aes = FALSE,
        curvature = 0.2,
        color = pal[11],
        alpha = .6,
        angle = 20,
        lineend = "round") +
        xlim(-1,14) +
        geom_label(data = label.df, 
                   aes(
                           x = x,
                           y = y,
                           label = env), 
                   nudge_x = -.4,
                   inherit.aes = FALSE, 
                   fontface = "bold", 
                   hjust = 0,
                   size = 3,fill = grey(.9,alpha = .5)) +
        guides(alpha = FALSE)

ggsave(filename = "matnet.pdf",height = 4,width = 5)

##Extended Data Fig. 5

meta.data.env <- metadata[,c(3,4,7:15)+1]
colnames(meta.data.env) <- c("Latitude","Longitude","MAT","MAP","Altitude",
                             "UVI","NPP", "pH","Moisture",
                             "SOC", "Evaporation")
meta.melt <- melt(meta.data.env)

richness.df <- data.frame(meta.melt, r = richness, s = shannon)
colnames(richness.df) <- c("variable","value","r","s")
richness.df$variable <- factor(richness.df$variable, 
                               labels =c("Latitude","Longitude",
                                         "MAT(degree*C)","MAP (mm)","Altitude (m)",
                                         "UVI (mW%.%m^-2)",
                                         "NPP(kg%.%m^-2%.%year^-1)", 
                                         "Soil~~pH",
                                         "Moisture (mm)","SOC (kg%.%m^-2)", 
                                         "Evaporation (mm)") 
)

ggplot(richness.df, aes(x= value, y = s)) +
        geom_point(alpha = .3, size = 1) +
        facet_wrap(~variable, switch = "x", 
                   strip.position = "bottom",
                   labeller = "label_parsed",
                   scales = "free",ncol = 4) +
        theme_bw()+ xlab("") +ylab("Shannon's H") +
        geom_smooth(method = "lm", 
                    formula = y ~ poly(x, 4), 
                    se = FALSE,
                    color = pal[1])+
        theme(panel.grid = element_blank(),
              aspect.ratio = 1,
              strip.background = element_blank(),
              strip.placement = "outside")

ggplot(richness.df, aes(x= value, y = log10(r))) +
        geom_point(alpha = .3, size = 1) +
        facet_wrap(~variable, switch = "x", 
                   strip.position = "bottom",
                   labeller = "label_parsed",
                   scales = "free",ncol = 4) +
        theme_bw()+ xlab("") +ylab("Richness (log10)") +
        geom_smooth(method = "lm", 
                    formula = y ~ poly(x, 4), 
                    se = FALSE,
                    color = pal[1])+
        theme(panel.grid = element_blank(),
              aspect.ratio = 1,
              strip.background = element_blank(),
              strip.placement = "outside")

# distance decay (Fig. 3c)

dist <- fread(input = "bc_dist_rm.csv")

dist.new <- data.frame(dist[,-1])
row.names(dist.new) <- dist$V1
colnames(dist.new) <- dist$V1
pos <- pmatch(metadata$Sample.id, colnames(dist.new))
dist.new1 <- dist.new[pos[!is.na(pos)] , pos[!is.na(pos)]]

coord.df <- metadata[, c(5,4)]
row.names(coord.df) <- metadata$Sample.id

geo.dist <-  distm(coord.df[row.names(dist.new1) , ])

bc.dist <- as.dist(dist.new1, diag = FALSE,upper = TRUE)
geo.dist1 <- as.dist(geo.dist,diag = FALSE,upper = TRUE)

dist.decay.df <- data.frame(bc = c(bc.dist), geo = c(geo.dist1)/1000)
g.dist <- ggplot(dist.decay.df, aes(x = geo, y = bc)) +
        geom_point(size = .1, color = pal[5]) +
        theme_bw() +
        theme(aspect.ratio = .5,panel.grid = element_blank()) +
        scale_x_continuous(breaks = c(0,5e+3,1e+4,1.5e+4,2e+4),
                           limits = c(0, 2.0e+4),
                           labels = expression(0, 
                                               0.5%*%10^4, 
                                               1.0%*%10^4, 
                                               1.5%*%10^4, 
                                               2.0%*%10^4))+
        xlab("Geographic distance (km)") + 
        ylab("Bray–Curtis dissimilarity") 

ggsave(g.dist ,filename = "g_dist.tiff",width = 4,height = 3,compression = "lzw")




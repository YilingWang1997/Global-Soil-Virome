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
vctable <- fread("GSV_vc_averageabundance_10kb_2021.csv")

# Figure 2
## alpha-diversity (figure 2a)
alpha.div <-
        data.frame(land = metadata$LandcoverClass,
                   continent = metadata$Continent,
                   richness =  richness,
                   shannon = shannon,
                   simpson = simpson)

alpha.div$land <-
        factor(
                alpha.div$land,
                levels = c(
                        "Agricultural land",
                        "Wetland",
                        "Grassland",
                        "Forest",
                        "Shrubland",
                        "Bare Land",
                        "Tundra"
                ),
                labels = c(
                        "Agriculture",
                        "Wetland",
                        "Grassland",
                        "Forest",
                        "Shrubland",
                        "Bare land",
                        "Tundra"
                )
                
        )

alpha.div$continent <- factor(
        alpha.div$continent,
        labels = c("Africa",
                   "Asia",
                   "Australia",
                   "Europe",
                   "North\nAmerica",
                   "South\nAmerica"),
        levels = c("Africa",
                   "Asia",
                   "Australia",
                   "Europe",
                   "North America",
                   "South America")
)

## Richness 
ggplot(na.omit(alpha.div),
       aes(x = land, y = log10(richness), color = land)) +
        geom_boxplot(alpha = 1, outlier.colour = NA,
                     draw_quantiles = .5
        ) +
        geom_beeswarm(size = .5, alpha = .2)+
        facet_grid(~continent, 
                   switch = "y",
                   scales = "free",
                   space = "free_x") +
        theme_bw() + ylab("Richness") + xlab("")+
        #scale_y_log10()+
        scale_color_manual(values = pal)+
        scale_y_continuous(breaks = seq(0,5,1),
                           limits = c(1,4),
                           labels = expression(10^0,10^1,10^2,10^3,10^4,10^5))+
        theme(legend.position = "none", 
              panel.grid = element_blank(),
              aspect.ratio = 10,
              axis.text.x = element_text(angle = -65,
                                         vjust = 0.5, 
                                         hjust = 0,
                                         size = 7),
              strip.text.x = element_text(size=7),
              strip.text.y = element_text(vjust = 0),
              strip.placement = "outside",
              strip.background.y = element_rect(fill = NA, color = NA))

## Shannon
ggplot(na.omit(alpha.div),
       aes(x = land, y = shannon, color = land)) +
        geom_boxplot(alpha = 1, outlier.colour = NA,
                     draw_quantiles = .5
        ) +
        geom_beeswarm(size = .5, alpha = .2)+
        facet_grid(~continent, 
                   switch = "y",
                   scales = "free",
                   space = "free_x") +
        theme_bw() + ylab("Shannon") + xlab("")+
        #scale_y_log10()+
        scale_color_manual(values = pal)+
        theme(legend.position = "none", 
              panel.grid = element_blank(),
              aspect.ratio = 10,
              axis.text.x = element_text(angle = -65,
                                         vjust = 0.5, 
                                         hjust = 0,
                                         size = 7),
              strip.text.x = element_text(size=7),
              strip.text.y = element_text(vjust = 0),
              strip.placement = "outside",
              strip.background.y = element_rect(fill = NA, color = NA))


## Simpson
ggplot(na.omit(alpha.div),
       aes(x = land, y = simpson, color = land)) +
        geom_boxplot(alpha = 1, outlier.colour = NA,
                     draw_quantiles = .5
        ) +
        geom_beeswarm(size = .5, alpha = .2)+
        facet_grid(~continent, 
                   switch = "y",
                   scales = "free",
                   space = "free_x") +
        theme_bw() + ylab("simpson") + xlab("")+
        scale_y_continuous(limits = c(.95,1)) +
        scale_color_manual(values = pal) +
        theme(legend.position = "none", 
              panel.grid = element_blank(),
              aspect.ratio = 10,
              axis.text.x = element_text(angle = -65,
                                         vjust = 0.5, 
                                         hjust = 0,
                                         size = 7),
              strip.text.x = element_text(size=7),
              strip.text.y = element_text(vjust = 0),
              strip.placement = "outside",
              strip.background.y = element_rect(fill = NA, color = NA))

### NMDS (figure 2b-c) 

nmds.df <- nmds %>%
        mutate(land = metadata$LandcoverClass) %>%
        mutate(cont = metadata$Continent) %>%
        arrange(land) 

eco.seq <- factor(nmds.df$land, 
                  levels = c("Agricultural land","Forest","Grassland",       
                             "Shrubland","Tundra" ,"Wetland","Bare Land"))

nmds.df <- nmds.df[order(eco.seq),]

ggplot(na.omit(nmds.df), aes(
        x = MDS1,
        y = MDS2,
        color = land
)) +
        geom_jitter(size = 2,
                    shape = 16,
                    alpha = .5) +
        xlim(-.25,.25) + ylim(-.25,.25) +
        theme_bw() + xlab("NMDS 1") + ylab("NMDS 2")+
        scale_color_manual(values = pal) +
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

cont.seq <- factor(nmds.df$cont, 
                   levels = c("North America","Asia","Australia","Europe","South America","Africa"))

nmds.df <- nmds.df[order(cont.seq),]

ggplot(na.omit(nmds.df), aes(
        x = MDS1,
        y = MDS2,
        color = cont
)) +
        geom_point(size = 2,
                   shape = 16,
                   alpha = .5) +
        xlim(-.25,.25) + ylim(-.25,.25) +
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


# similarity among ecosystems (figure 2d)
g.eco <-
        graph_from_adjacency_matrix(
                as.matrix(1 - ecodist),
                mode = "undirected",
                weighted = TRUE,
                diag = FALSE
        )
V(g.eco)$name <- c("Agriculture","Bare\nland","Forest","Grassland",       
                   "Shrubland","Tundra" ,"Wetland")

g.eco <- g.eco %>%
        delete_edges(c(1:19)[E(g.eco)$weight < .01])

ly.eco <- layout.circle(g.eco)

ggraph(g.eco, layout = ly.eco[c(5,7,3,4,1,6,2),]) +
        geom_edge_diagonal(
                aes(edge_width = weight),
                alpha = .5,
                color = pal_npg(alpha = .5)(4)[4],
                lineend = "round"
        ) +
        geom_node_point(size = 8, 
                        color = pal[1:7]) +
        geom_node_text(
                aes(label = name),
                nudge_y = c(.2,-.3,.3,.3,-.3,-.3,.3),#[c(5,7,3,4,1,6,2)],
                nudge_x = c(-.5,0,0,-.5,0,0,0),#[c(5,7,3,4,1,6,2)],
                size = 3, 
                family = "", 
                lineheight = .8
        ) +
        scale_edge_width_continuous(range = c(.5, 1.5), 
                                    name = "Jaccard\nsimilarity") +
        theme_graph(base_family = 'Helvetica') + 
        xlim(-2, 2) + ylim(-1.5, 1.5) +
        theme(aspect.ratio = .75,text = element_text(family = "Helvetica"))

ggsave(filename = "figures/subnet_eco.pdf", height = 4, width = 5)

# similarity among continents (figure 2e)
g.area <-
        graph_from_adjacency_matrix(
                as.matrix(1 - areadist),
                mode = "undirected",
                weighted = TRUE,
                diag = FALSE
        )
V(g.area)$name <- c("Africa","Asia","Australia","Europe","North\nAmerica","South\nAmerica")

ly.area <- layout.circle(g.area)

ggraph(g.area, layout = ly.area[c(3,2,1,4,5,6),]) +
        geom_edge_diagonal2(aes(edge_width = weight), 
                            alpha = .5, 
                            color = pal_npg(alpha = .6)(6)[6],
                            lineend = "round") +
        geom_node_point(size = 8, color = pal[1:6]) +
        geom_node_text(
                aes(label = name),
                size = 3, 
                nudge_y = c(.3,.3,0,0,-.3,-.3),
                nudge_x = c(0,0,.5,-.5,0,0),
                family = "Helvetica", 
                lineheight = .8
        ) +
        scale_edge_width_continuous(range = c(.3, 1.5), 
                                    name = "Jaccard\nsimilarity") +
        theme_graph(base_family = 'Helvetica') + xlim(-2, 2) + ylim(-1.2, 1.2) +
        theme(aspect.ratio = .6,legend.position = "none")



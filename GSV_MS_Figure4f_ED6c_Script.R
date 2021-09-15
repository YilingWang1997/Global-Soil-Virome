#co-occurrence network
vir.present <- fread("vir_present.txt")
vir.table <- fread(file = "OTU_GSV_average_2115.csv")
vir.table.clean <- t(vir.table[vir.present$x > 105, -1])
vir.01 <- vir.table.clean
vir.01[vir.table.clean != 0] <- 1


# ecosystem combine
vir.table.clean

require(igraph)
require(Hmisc)

vir.cor <- cor(vir.table.clean, method = "spearman")
vir.cor.p <- Hmisc::rcorr(vir.table.clean, type =  "spearman")

# network enhancement
require(neten)
vir.cor.en <- neten::Network_Enhancement(vir.cor)

# p-adjust
vir.cor.p.adj <- matrix(p.adjust(vir.cor.p$P, method = "fdr"),nrow = 1885)

# adj matrix 

vir.mat <- matrix(0, nrow = 1885,ncol = 1885)
vir.mat[vir.cor.en > 1 & vir.cor.p.adj < 0.05] <- 1

# generate network
g.vir.conet <-
  graph_from_adjacency_matrix(vir.mat, 
                              mode = "undirected", 
                              weighted = TRUE)

V(g.vir.conet)$name <- c(vir.table[vir.present$x > 105, 1])

write.graph(g.vir.conet, 
            file = "vir_conet.graphml",
            format = "graphml")

g.vir.conet <- read.graph(file = "vir_conet_pos.graphml",
                          format = "graphml")

ly.vir <- data.frame(V(g.vir.conet)$x, V(g.vir.conet)$y)

require(ggraph)

g.eco1 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(eco1))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[3])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "figures/g_eco1.pdf",
       height = 4,
       width = 4)


g.eco2 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(eco2))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[3])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "figures/g_eco2.pdf",
       height = 4,
       width = 4)

g.eco3 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(eco3))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[3])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "figures/g_eco3.pdf",
       height = 4,
       width = 4)

g.eco4 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(eco4))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[3])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "figures/g_eco4.pdf",
       height = 4,
       width = 4)

g.eco5 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(eco5))) +
  theme_graph(background = "black",base_family = 'Helvetica') +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[3])) +
  theme(aspect.ratio = 1, 
        title = element_text(color = "white"),
        legend.position = "none")

ggsave(filename = "figures/g_eco5.pdf",
       height = 4,
       width = 4)

g.eco6 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(eco6))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[3])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "figures/g_eco6.pdf",
       height = 4,
       width = 4)

g.eco7 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(eco7))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[3])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "figures/g_eco7.pdf",
       height = 4,
       width = 4)

g.area1 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(area1))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[5])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "figures/g_area1.pdf",
       height = 4,
       width = 4)

g.area2 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(area2))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[5])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "figures/g_area2.pdf",
       height = 4,
       width = 4)

g.area3 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(area3))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[5])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "figures/g_area3.pdf",
       height = 4,
       width = 4)

g.area4 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(area4))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[5])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "figures/g_area4.pdf",
       height = 4,
       width = 4)

g.area5 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(color = pal_npg(alpha = .6)(6)[5]) +
  theme_graph(background = "black") +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "figures/g_area5.pdf",
       height = 4,
       width = 4)

g.area6 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(area6))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[5])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "figures/g_area6.pdf",
       height = 4,
       width = 4)

### degree distribution
deg.dis <- table(degree(g.vir.conet))

ggplot(melt(deg.dis), aes(x = Var1, y = value)) +
  geom_bar(stat = "identity", fill = pal_npg()(6)[3]) +
  scale_x_continuous(n.breaks = 8) +
  theme_bw() + xlab("Degree") + ylab("Frequency") +
  theme(
    panel.grid = element_blank(),
    aspect.ratio = .5,
    axis.title = element_text(size = 9)
  )

ggsave(filename = "figures/degree_dis.pdf", 
       width = 3, 
       height = 2)

vir.com.eco <- apply(vir.01, MARGIN = 2, tapply, metadata$LandcoverClass, sum)
vir.com.area <- apply(vir.01, MARGIN = 2, tapply, metadata$Continent, sum)

vir.com.eco[vir.com.eco != 0] <- 1
vir.com.area[vir.com.area != 0] <- 1

V(g.vir.conet)$eco1 <- vir.com.eco[1, ]
V(g.vir.conet)$eco2 <- vir.com.eco[2, ]
V(g.vir.conet)$eco3 <- vir.com.eco[3, ]
V(g.vir.conet)$eco4 <- vir.com.eco[4, ]
V(g.vir.conet)$eco5 <- vir.com.eco[5, ]
V(g.vir.conet)$eco6 <- vir.com.eco[6, ]
V(g.vir.conet)$eco7 <- vir.com.eco[7, ]

V(g.vir.conet)$area1 <- vir.com.area[1, ]
V(g.vir.conet)$area2 <- vir.com.area[2, ]
V(g.vir.conet)$area3 <- vir.com.area[3, ]
V(g.vir.conet)$area4 <- vir.com.area[4, ]
V(g.vir.conet)$area5 <- vir.com.area[5, ]
V(g.vir.conet)$area6 <- vir.com.area[6, ]

V(g.vir.conet)$size <- degree(g.vir.conet)

# simulated network

ba.net <- ba.game(n = 1885, m = 36)
er.net <- erdos.renyi.game(n = 1885,p.or.m = 67445,type = "gnm")

deg.df <- rbind(data.frame(melt(table(degree(ba.net))), net = "BA model"),
                data.frame(melt(table(degree(er.net))), net = "ER model"),
                data.frame(melt(table(degree(g.vir.conet))), net = "GSV network"))

ggplot(deg.df, aes(x = Var1, y = value/67445, color = net)) +
  geom_point(alpha = .5)+
  facet_grid(~net, scales = "free") +
  scale_x_log10()+
  scale_y_log10(breaks = c(.00001, .0001,.001,.01), 
                labels = expression(10^-5, 10^-4,10^-3,10^-2)) +
  scale_color_manual(values = pal[c(1,2,6)])+
  theme_bw()+ 
  xlab(expression(italic(k))) + 
  ylab(expression(italic(p[k]))) +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        legend.position = "none")
ggsave("figures/deg_pro.pdf",width = 6,height = 3)
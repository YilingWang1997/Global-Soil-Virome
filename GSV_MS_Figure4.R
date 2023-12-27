#fig4a,b, network
require(igraph)
require(Hmisc)
require(ggraph)
require(ggsci)
library(reshape2)

otu <- read.csv("OTU_GSV_average_2022_1859.csv",row.names = 1)
otu <- t(otu)
eco <- read.csv("1799biomcont.csv")
rownames(eco) <- eco$Sample.id
otu <- otu[rownames(eco),]
vir.com.eco <- apply(otu, MARGIN = 2, tapply,eco$Biome, sum)
rownames(vir.com.eco.t)<-gsub("=",".",rownames(vir.com.eco.t))
vir.com.eco <- data.frame(vir.com.eco)
vireco1020 <- read.csv("vir.com.eco_1020.csv",row.names=1)
vir.com.eco2 <- vir.com.eco[,colnames(vireco1020)]

write.csv(vir.com.eco2,"vir.com.eco_1230.csv")

vir.com.eco <- read.csv("vir.com.eco_1230.csv",row.names = 1)
vir.com.eco[vir.com.eco !=0] <-1

nodesinfo <- read.csv("nodes_co.graphml.csv")
nodesinfo <- nodesinfo[,c("Id","v_name")]

g.vir.conet <- read.graph(file = "co_pos.graphml",
                          format = "graphml")

colnames(vir.com.eco) <- nodesinfo$Id[pmatch(nodesinfo$v_name,colnames(vir.com.eco))]
vir.com.eco <- vir.com.eco[,get.vertex.attribute(g.vir.conet)$id]

V(g.vir.conet)$eco1 <- vir.com.eco[1, ]
V(g.vir.conet)$eco2 <- vir.com.eco[2, ]
V(g.vir.conet)$eco3 <- vir.com.eco[3, ]
V(g.vir.conet)$eco4 <- vir.com.eco[4, ]
V(g.vir.conet)$eco5 <- vir.com.eco[5, ]
V(g.vir.conet)$eco6 <- vir.com.eco[6, ]
V(g.vir.conet)$eco7 <- vir.com.eco[7, ]
V(g.vir.conet)$eco8 <- vir.com.eco[8, ]

ly.vir <- data.frame(V(g.vir.conet)$x, V(g.vir.conet)$y)

eco11 <- vir.com.eco[1, ]
g.eco1 <- ggraph(g.vir.conet, layout = ly.vir) +
  geom_node_point(aes(color = factor(eco11))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[3])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "subnetwork/g_eco1.pdf",
       height = 4,
       width = 4)

eco22 <- vir.com.eco[2, ]
g.eco2 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(eco22))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[3])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "subnetwork/g_eco2.pdf",
       height = 4,
       width = 4)

eco33 <- vir.com.eco[3, ]
g.eco3 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(eco33))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[3])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "subnetwork/g_eco3.pdf",
       height = 4,
       width = 4)

eco44 <- vir.com.eco[4, ]
g.eco4 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(eco44))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[3])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "subnetwork/g_eco4.pdf",
       height = 4,
       width = 4)

eco55 <- vir.com.eco[5, ]
g.eco5 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(eco55))) +
  theme_graph(background = "black",base_family = 'Helvetica') +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[3])) +
  theme(aspect.ratio = 1, 
        title = element_text(color = "white"),
        legend.position = "none")

ggsave(filename = "subnetwork/g_eco5.pdf",
       height = 4,
       width = 4)

eco66 <- vir.com.eco[6, ]
g.eco6 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(eco66))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[3])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "subnetwork/g_eco6.pdf",
       height = 4,
       width = 4)

eco77 <- vir.com.eco[7, ]
g.eco7 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(eco77))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[3])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "subnetwork/g_eco7.pdf",
       height = 4,
       width = 4)

eco88 <- vir.com.eco[8, ]
g.eco88 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(eco88))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[3])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "subnetwork/g_eco8.pdf",
       height = 4,
       width = 4)


vir.com.area <- read.csv("vir.com.cont.csv",row.names = 1)
vir.com.area[vir.com.area !=0] <-1

colnames(vir.com.area) <- nodesinfo$Id[pmatch(nodesinfo$v_name,colnames(vir.com.area))]
vir.com.area <- vir.com.area[,get.vertex.attribute(g.vir.conet)$id]

V(g.vir.conet)$area1 <- vir.com.area[1, ]
V(g.vir.conet)$area2 <- vir.com.area[2, ]
V(g.vir.conet)$area3 <- vir.com.area[3, ]
V(g.vir.conet)$area4 <- vir.com.area[4, ]
V(g.vir.conet)$area5 <- vir.com.area[5, ]
V(g.vir.conet)$area6 <- vir.com.area[6, ]

area11 <- vir.com.area[1, ]
g.area1 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(area11))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[5])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "subnetwork/g_area1.pdf",
       height = 4,
       width = 4)

area22 <- vir.com.area[2, ]
g.area2 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(area22))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[5])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "subnetwork/g_area2.pdf",
       height = 4,
       width = 4)

area33 <- vir.com.area[3, ]
g.area3 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(area33))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[5])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "subnetwork/g_area3.pdf",
       height = 4,
       width = 4)

area44 <- vir.com.area[4, ]
g.area4 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(area44))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[5])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "subnetwork/g_area4.pdf",
       height = 4,
       width = 4)

area55 <- vir.com.area[5, ]
g.area4 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(area55))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[5])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "subnetwork/g_area5.pdf",
       height = 4,
       width = 4)

area66 <- vir.com.area[6, ]
g.area6 <- ggraph(g.vir.conet, layout =  ly.vir) +
  geom_node_point(aes(color = factor(area66))) +
  theme_graph(background = "black") +
  scale_color_manual(values = c(grey(level = .8, alpha = .1), 
                                pal_npg(alpha = .6)(6)[5])) +
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(filename = "subnetwork/g_area6.pdf",
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




# simulated network
pal <- c("#C75E48","#F9B483","#6C958F")
ba.net <- ba.game(n = 5117,m = 25)
er.net <- erdos.renyi.game(n = 5117,p.or.m = 125256,type = "gnm")

deg.df <- rbind(data.frame(melt(table(degree(ba.net))), net = "BA model"),
                data.frame(melt(table(degree(er.net))), net = "ER model"),
                data.frame(melt(table(degree(g.vir.conet))), net = "GSV network"))

ggplot(deg.df, aes(x = Var1, y = value/125256, color = net)) +
  geom_point(alpha = .5)+
  facet_grid(~net, scales = "free") +
  scale_x_log10()+
  scale_y_log10(breaks = c(.00001, .0001,.001,.01), 
                labels = expression(10^-5, 10^-4,10^-3,10^-2)) +
  scale_color_manual(values = pal)+
  theme_bw()+ 
  xlab(expression(italic(k))) + 
  ylab(expression(italic(p[k]))) +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        legend.position = "none")
ggsave("deg_pro.pdf",width = 6,height = 3)

#fig3c&d
virus <- read.csv("virus_forbipartite.csv",row.names = 1)
host <- read.csv("c_host_forbipartite_realabun.csv",row.names = 1)

sample <- read.csv("../SupplementaryTable/Supplementary_Table1.csv")
min(sample$Size)
sample <- sample[,c("Sample.id","Size")]

sample <- sample[pmatch(colnames(host),sample$Sample.id),]
host_norm <- t(t(host)/sample$Size)*mean(sample$Size)

virus <- virus[,colnames(virus) %in% colnames(host_norm)]
virus <- virus[,pmatch(colnames(host_norm),colnames(virus))]

data <- rbind(host_norm,virus)

links <- read.csv("../../MAG_virus/viruslinkages_MAG.csv")
cdb <- read.csv("Cdb.csv")
cdb$genome <- gsub(".fa","",cdb$genome)
links <- merge(links,cdb,by.x="bin",by.y="genome")

uniquelink <- unique(links[,c("secondary_cluster","virus")])
uniquelink <- uniquelink[grep("GSV_",uniquelink$virus),]

rp <- NULL
for (i in 1:dim(uniquelink)[1]){
  print(i)
    h <- data[uniquelink$secondary_cluster[i],]
    v <- data[uniquelink$virus[i],]
    hv <- cbind(t(h),t(v))
    hv <- hv[hv[,1]!=0 & hv[,2] !=0 ,]
    if (length(dim(hv))>0){
      if (dim(hv)[1] >= 3){
    r <- lm(log10(hv[,2])~log10(hv[,1]))
    res_su <- data.frame(summary(r)$coefficients)
    rp <- rbind(rp,data.frame(uniquelink[i,],res_su$Estimate[2],res_su$Pr...t..[2],dim(hv)[1]))}
    else{}}
    else{}
}

#rp2 <- cbind(uniquelink,rp)
#write.csv(rp,"rp_new.csv")
########两张图在这里#########3
library(ggplot2)
rp <- read.csv("rp_new.csv")
rp2_f <- rp[rp$res_su.Pr...t...2.<0.05,]
rp2_f <- rp2_f[rp2_f$dim.hv..1.>=18,]
ggplot(rp2_f, aes(x=res_su.Estimate.2.)) + 
  geom_histogram(aes(y=..density..),binwidth=0.1,alpha=0.2,color="#3760A6",fill="#3760A6")+
  geom_density(alpha=.4,fill="#3760A6",color="#3760A6")+
  geom_vline(xintercept=0,linetype="dashed",color="grey")+
  xlab("Pearson correlation coefficient (r) \nof virus vs. host abundances")+
  ylab("Density")+
  theme_classic(base_size = 17,base_line_size = .3)


hist(rp2_f$res_su.Estimate.2.,
     nclass = 100 ,col = "#3760A6",
     border = "white",xlab = "Pearson correlation coefficient (r) \nof virus vs. host abundances",ylab="Frequence",
     main = "",cex.axis=1.5,cex.lab=1.5)

rhv <- NULL
for (i in 1:dim(uniquelink)[1]){
  print(i)
  h <- data[uniquelink$secondary_cluster[i],]
  v <- data[uniquelink$virus[i],]
  hv <- cbind(t(h),t(v))
  hv <- hv[hv[,1]!=0 & hv[,2] !=0 ,]
  if (length(dim(hv))>0){
  if (dim(hv)[1] >= 3){
  result <- lm(log10(hv[,1])~log10(hv[,2]/hv[,1]))
  res_su <- data.frame(summary(result)$coefficients)
  rhv <- rbind(rhv,data.frame(uniquelink[i,],res_su$Estimate[2],res_su$Pr...t..[2],dim(hv)[1]))}
  else{}}
  else{}
}
write.csv(rhv,"rhv_new.csv")

rhv <- read.csv("rhv_new.csv")
rhv_f <- rhv[rhv$res_su.Pr...t...2.<0.05,]
rhv_f <- rhv[rhv$dim.hv..1. >= 18,]

rhv_f <- rhv_f[rhv_f$res_su.Pr...t...2.<0.05,]
hist(rhv_f$res_su.Estimate.2.)
ggplot(rhv_f[rhv_f$res_su.Estimate.2.<=0.5 & rhv_f$res_su.Estimate.2.>=-2,], aes(x=res_su.Estimate.2.)) + 
  geom_histogram(aes(y=..density..),binwidth=0.1,alpha=0.2,color="#E19B26",fill="#E19B26")+
  geom_density(alpha=.4,fill="#E19B26",color="#E19B26")+
  geom_vline(xintercept=0,linetype="dashed",color="#E19B26")+
  xlab("Pearson correlation coefficient (r) \nof VHRs vs. host abundances")+
  ylab("Density")+
  theme_classic(base_size = 17,base_line_size = .3)

hist(rhv_f$res_su.Estimate.2.[rhv_f$res_su.Estimate.2.<=5])
hist(rhv_f$res_su.Estimate.2.[rhv_f$res_su.Estimate.2.<=0.5 & rhv_f$res_su.Estimate.2.>=-2],
     nclass = 100 ,col = "darkgrey",
     border = "white",xlab = "Pearson correlation coefficient (r) \nof VHRs vs. host abundances",ylab="Frequence",main = "")


maginfo <- read.csv("../../MAG_virus/mag95info.csv")
maginfo <- maginfo[maginfo$cluster %in% rownames(host_norm),]

points <- NULL
for (i in 1:dim(rhv_f)[1]){
  print(i)
  h <- data[rhv_f$secondary_cluster[i],]
  v <- data[rhv_f$virus[i],]
  hv <- data.frame(cbind(t(h),t(v)))
  hv <- hv[hv[,1]!=0 & hv[,2] !=0 ,]
  hv$hlog <- log10(hv[,1])
  hv$hvlog <- log10(hv[,2]/hv[,1])
  hv$group <- rep(i,dim(hv)[1])
  hv$host <- rep(rhv_f$secondary_cluster[i],dim(hv)[1])
  hv$virus <- rep(rhv_f$virus[i],dim(hv)[1])
  colnames(hv) <- c("h","v","hlog","hvlog","group","host","virus")
  points <- rbind(points,hv)}


points <- NULL
for (i in 1:dim(uniquelink)[1]){
  print(i)
  h <- data[uniquelink$secondary_cluster[i],]
  v <- data[uniquelink$virus[i],]
  hv <- data.frame(cbind(t(h),t(v)))
  hv <- hv[hv[,1]!=0 & hv[,2] !=0 ,]
  hv$hlog <- log10(hv[,1])
  hv$hvlog <- log10(hv[,2]/hv[,1])
  hv$group <- rep(i,dim(hv)[1])
  hv$host <- rep(uniquelink$secondary_cluster[i],dim(hv)[1])
  hv$virus <- rep(uniquelink$virus[i],dim(hv)[1])
  colnames(hv) <- c("h","v","hlog","hvlog","group","host","virus")
  points <- rbind(points,hv)}

write.csv(points,"vhr-h.csv")
ggplot(points, aes(hlog, hvlog, color = group)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

test <- lm(formula = hvlog~hlog,points)

gf <- data.frame(table(points$group))
points_filter <- points[points$group %in% gf$Var1[gf$Freq>=18],]
ggplot(points_filter, aes(log10(h), log10(v))) +
  geom_point() +
  geom_smooth(aes(group=group),method = "lm", se = FALSE)

ggplot(points_filter, aes(h, v)) +
  geom_point() +
  geom_smooth(aes(group=group),method = "lm", se = FALSE)

ggplot(test, aes(log10(h), log10(v))) +
  geom_point() +
  geom_smooth(aes(group=group),method = "lm", se = FALSE)



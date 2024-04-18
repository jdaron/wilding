library(reshape)
library(ggplot2)
library(ggthemes)
library(tidyverse)
library("factoextra")
library("FactoMineR")
library(ggrepel)
library(plyr)
library(reshape2)
library(scatterpie)
library(sp)
library(RColorBrewer)
library(grid)
library(gtable)
library("gridExtra")
library(measurements)

setwd("~/bioInf/wilding/github/wilding/paper/input/figure1/")

##########################
### Figure 1:
##########################

### Figure 1A: map Africa and Gabon
###===============

# 1. Create map plot object
worldmap <- borders("world",
                    colour = "gray20",
                    fill = "white"
)

colour = "#fefefe"
fill = "gray50"

th =  theme(
  panel.grid.major = element_blank(),
  panel.background = element_rect(fill = "aliceblue"),
  axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  axis.line=element_blank(),
#  legend.position="none",
  panel.grid.minor=element_blank()
)

p = ggplot() + worldmap +
  coord_quickmap() +
  coord_sf(xlim = c(-20, 60), ylim = c(-25, 40), expand = FALSE) +
  th
p

# 2. Load meta data and process
ag = read.table("ag1000g.samplesLocation.txt", header = T, dec= ".")
head(ag)

ag.d = dcast(data=ag, country+site_name+lat+lon~sp, value.var = "nbsamples")
ag.d[is.na(ag.d)] = 0
ag.d$nbsamples=ag.d$col+ag.d$gam+ag.d$mixte
ag.d$col=(ag.d$col/(ag.d$col+ag.d$gam+ag.d$mixte))
ag.d$gam=(ag.d$gam/(ag.d$col+ag.d$gam+ag.d$mixte))
ag.d$mixte=(ag.d$mixte/(ag.d$col+ag.d$gam+ag.d$mixte))

wil = read.table("~/bioInf/wilding/github/wilding/paper/figuresScripts/wilding.samplesLocation.txt", header = T, dec= ".")
wil.d = dcast(data=wil, country+site_name+lat+lon~sp, value.var = "nbsamples")
wil.d$nbsamples=wil.d$col
wil.d$col=1
wil.d$gam=0
wil.d$mixte=0

# 3. Get radius
maxr = 2
ag.d$radius = (ag.d$nbsamples/max(ag.d$nbsamples))*maxr
wil.d$radius = (wil.d$nbsamples/max(ag.d$nbsamples))*maxr
ag.d$radius[which(ag.d$nbsamples<10)]=NA

pf = p + geom_scatterpie(aes(x=lon, y=lat, group=site_name, r=radius), data=ag.d,
                    cols=c("col", "gam", "mixte"), color=NA) +
  geom_scatterpie(aes(x=lon, y=lat, group=site_name, r=radius), data=wil.d,
                      cols=c("col", "gam", "mixte"), color=NA) +
  coord_sf(xlim = c(-20, 60), ylim = c(-30, 30), expand = FALSE)
pf  

ggsave("worldmap.pdf", units = "mm", width = 120, height = 120, useDingbats=FALSE) # export PDF 7x9

alt = getData("alt", country="Gabon", path = tempdir())
plot(alt)

pdf("gabon.pdf", width = conv_unit(100, "mm", "inch"), height = conv_unit(100, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
plot(alt)
points(wil.d$lon[1], wil.d$lat[1], pch = 19, cex=0.5)
points(wil.d$lon[3], wil.d$lat[3], pch = 19, cex=0.5)
dev.off()

### Figure 1B: plot La Lope map ecosystems
###===============
#Load libraries
library (sf)
library(tidyverse)
library(readr)
library(ggspatial)

## layers
pnl_hab <- st_read("~/bioInf/wilding/github/wilding/paper/input/figure1/SEGC_burn_zones/Habitat_Types_SEGC.shp")  %>%
  st_transform(4326) %>%   ## transfort to WGS84
  mutate(Savanna = factor(SAVANNA))

pnl_roads <- st_read("~/bioInf/wilding/github/wilding/paper/input/figure1/SEGC_burn_zones/Roads_Lope.shp")  %>%
  st_transform(4326) 

pnl_rail <- st_read("~/bioInf/wilding/github/wilding/paper/input/figure1/SEGC_burn_zones/Railroad_Lope.shp")  %>%
  st_transform(4326) 

pnl_routes <- st_read("~/bioInf/wilding/github/wilding/paper/input/figure1/SEGC_burn_zones/Routes_SEGC.shp")  %>%
  st_transform(4326) 

pnl_village <- st_read("~/bioInf/wilding/github/wilding/paper/input/figure1/SEGC_burn_zones/Villages_Lope.shp")  %>%
  st_transform(4326) 

# In pnl_hab, Savana code, 0 = forest, 1 = savanna, 2 = water
col_map_elts = c("0"="seagreen", "1"="khaki1", "2"="dodgerblue1", # for landcover 
                 "Sylvatic"="black", "Domestic"= "black",
                 "village"="gray90", "main_road"="gray50", 
                 "rails"="black", "sec_roads"="rosybrown2", "shape"=23)

base_map = ggplot() + 
  geom_sf(data= pnl_hab, aes(fill = Savanna), show.legend = FALSE) +
  geom_sf(data=pnl_village, aes(fill="village")) +
  geom_sf(data= pnl_routes[pnl_routes$TYPE=="Road",], aes(col="sec_roads")) + 
  geom_sf(data=pnl_rail, aes(linetype="twodash", col="rails"), size=1.5)+
  geom_sf(data=pnl_roads, aes(col="main_road"), linetype="solid", size=1.6)+
  #  geom_point(data = host_pref_map, aes(x=DECIMAL_LONGITUDE, y=DECIMAL_LATITUDE, shape="21", fill = HABITAT), color="white", size=3.5) + 
  coord_sf(xlim = c(11.55, 11.675), ylim = c(-0.24, -0.05), expand = FALSE)  + 
  annotation_scale(location="br") +
  annotation_north_arrow(location="br", 
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         pad_y = unit(0.6, "cm"), pad_x = unit(0.5, "cm"),
                         style = north_arrow_fancy_orienteering(text_size = 8)) + 
  scale_color_manual(name="", values = col_map_elts) + 
  scale_fill_manual(name="", values = col_map_elts) +
  scale_linetype_manual(name="", values="twodash") +
  scale_shape_manual(name="", values = 21, breaks = "21") +
  theme_void() + 
  theme(legend.position = "none")

pdf("lalope_ecosystem.pdf", width = conv_unit(100, "mm", "inch"), height = conv_unit(1
                                                                                                                                             00, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
base_map
dev.off()

### Figure 1C: PCA
###===============
# info
infowi <- read.table("wilding.samples.meta.txt", h=T)
infowi = infowi[,c("id", "population")]
infoag <- read.table("ag1000g.samples.meta.txt", h=T, sep="\t")
colnames(infoag)[1] = "id"
infoag = infoag[,c("id", "population")]

info = rbind(infoag, infowi)

# colors
colag <- read.table("ag1000g.samples.colors.txt", h=T, sep="\t")
colag$shape = "ag"
colwi <- read.table("wilding.samples.colors.txt", h=T, sep="\t")
colwi$shape = "wi"

cols = rbind(colag, colwi)
cols$order = 1:dim(cols)[1]

group.colors = as.character(cols$colors[cols$order])
names(group.colors) = factor(cols$population, levels=cols$population[cols$order])
group.shape = c(rep(21, dim(colag)[1]), rep(24, dim(colwi)[1]))
names(group.shape) = factor(cols$population, levels=cols$population[cols$order])

### PCA AG1000G and wilding
prefix = "ag1000g_phase2_col.wilding.merged.biallelic.3.pca"
inF = paste(prefix, ".eigenval.txt", sep="")
eigenval <- data.frame(read.table(inF, header=FALSE, skip=0, sep=" "))
sum(eigenval)
eigenval = eigenval$V1*100
eigenval = round(eigenval, digits = 2)

inF = paste(prefix, ".eigenvec.txt", sep="")
eigenvec <- data.frame(read.table(inF, header=FALSE, skip=0, sep="\t"))
eigenvec = eigenvec[,1:7]
colnames(eigenvec) = c("id", "EV1", "EV2", "EV3", "EV4", "EV5", "EV6")

merge = c()
merge = merge(eigenvec, info, by="id")
head(merge)

mat = data.frame()
mat <- data.frame(
  id = merge$id,
  pop = merge$population,
  #  location = merge$location,
  #  site = merge$site,
  EV1 = merge$EV1,    # the first eigenvector
  EV2 = merge$EV2,    # the second eigenvector
  EV3 = merge$EV3,    # the second eigenvector
  EV4 = merge$EV4,    # the second eigenvector
  EV5 = merge$EV5,    # the second eigenvector
  EV6 = merge$EV6,    # the second eigenvector
  stringsAsFactors = FALSE
)

mat$pop = factor(mat$pop, levels=cols$population[cols$order])

xlab = paste("EV1 (", round(eigenval[1], digits = 2), "%)", sep = "")
ylab = paste("EV2 (", round(eigenval[2], digits = 2), "%)", sep = "")

size = 2
alpha = 1
stroke = 0

p2 = ggplot() +
  geom_point(data=mat, aes(x=EV1, y=EV2, color=pop, fill=pop, shape=pop), size=size, alpha=alpha, stroke = stroke) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  scale_shape_manual(values = group.shape) +
  xlab(xlab) +
  ylab(ylab) +
  theme_classic() +
  theme(legend.position="none")

p2

xlab = paste("EV1 (", round(eigenval[1], digits = 2), "%)", sep = "")
ylab = paste("EV3 (", round(eigenval[3], digits = 2), "%)", sep = "")

p3 = ggplot() +
  geom_point(data=mat, aes(x=EV1, y=EV3, color=pop, fill=pop, shape=pop), size=size, alpha=alpha, stroke = stroke) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  scale_shape_manual(values = group.shape) +
  xlab(xlab) +
  ylab(ylab) +
  theme_classic() +
  theme(legend.position="none")
p3

xlab = paste("EV2 (", round(eigenval[2], digits = 2), "%)", sep = "")
ylab = paste("EV3 (", round(eigenval[3], digits = 2), "%)", sep = "")

p4 = ggplot() +
  geom_point(data=mat, aes(x=EV2, y=EV3, color=pop, fill=pop, shape=pop), size=size, alpha=alpha, stroke = stroke) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  scale_shape_manual(values = group.shape) +
  xlab(xlab) +
  ylab(ylab) +
  theme_classic() +
  theme(legend.position="none")
p4

pdf("pca.ev1.ev2.ev3.pdf", width = conv_unit(240, "mm", "inch"), height = conv_unit(80
                                                                                                                                            , "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(p2, p3, p4, ncol = 3, nrow = 1)
dev.off()

xlab = paste("EV1 (", round(eigenval[1], digits = 2), "%)", sep = "")
ylab = paste("EV4 (", round(eigenval[4], digits = 2), "%)", sep = "")

p5 = ggplot() +
  geom_point(data=mat, aes(x=EV1, y=EV4, color=pop, fill=pop, shape=pop), size=size, alpha=alpha, stroke = stroke) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  scale_shape_manual(values = group.shape) +
  xlab(xlab) +
  ylab(ylab) +
  theme_classic() +
  theme(legend.position="none")
p5

# 3D plot
library("plot3D")
head(mat)
df_cols = data.frame(pop = names(group.colors), cols=unname(group.colors))
mat = merge(mat, df_cols, by="pop")
mat = mat[order(mat$pop),]
mat$pop = factor(mat$pop)
mat$shape = 19
mat$shape[which(mat$pop%in%c("LBVwil", "LPdom", "LPfor"))] = 17

head(mat)
tail(mat)

# width = conv_unit(1100, "mm", "inch"), height = conv_unit(620, "mm", "inch"),
pdf(file="pca.3D.pdf", useDingbats=FALSE)
scatter3D(mat$EV1, mat$EV2, mat$EV3,
          pch = mat$shape,
          colvar = as.numeric(mat$pop), 
          col = unique(as.character(mat$cols)),
          ticktype = "detailed",   bty ="b2", 
          type = "h",
          cex=1.5,
          theta = -50, phi = 10)
dev.off()

### Figure 1E: Admixture
###=====================
rm(order) # remove vector order cause it creates issue bellow
####### K=3
# 1.1 Load data
colsFrame = c("#B11217", "#FBA081", "#009e60")
# colsFrame = c("#B11217", "#FBA081", "#3a75c4")

coord = read.table("popIndex.txt")
colnames(coord) = c("popId", "popName")

inF = "ag1000gCol_wilding.admixture.K3.txt"
ad = read.table(inF)
colnames(ad)[1] = c("popId")
colnames(ad)[2:dim(ad)[2]] = paste(rep("k", dim(ad)[2]-1), 1:(dim(ad)[2]-1), sep="")
head(ad)

# 2.2 Plot bargraph admixture
# 2.2.1 Reorder individu within each population to make barplot pretty
ad$ind = 1:dim(ad)[1]
reorder = c() # reorder barplot
for(i in unique(ad$popId)){
  tmp = ad[which(ad$popId==i),]
  if(dim(tmp)[1]>1){
    reorder = c(reorder, tmp$ind[do.call(order, rev(as.list(tmp[,2:(dim(tmp)[2]-1)])))] )
  } else{
    reorder = c(reorder, tmp$ind)
  }
}
length(reorder)
length(ad$popId)
ad = ad[match(reorder, ad$ind),]
ad$newOrder = 1:dim(ad)[1]
head(ad)

# 2.2.2 Add spacer betwen each pop
spacerSize = 5
tmp = ad
for(i in unique(tmp$popId)){
  print(i)
  tmp$newOrder[-which(tmp$popId<i)] = tmp$newOrder[-which(tmp$popId<i)] + spacerSize
}
ad = tmp

# 2.3 Make dataframe for population label
labs = matrix(data=NA, nrow = length(unique(ad$popId)), ncol = 3)
for(i in unique(ad$popId)){
  md = median(ad$newOrder[which(ad$popId==i)])
  labs[i, ] = c(md, i, -0.1)
}
colnames(labs) = c("xcoord", "popId", "ycoord")
labs = data.frame(labs)
labs = merge(labs, coord, by="popId")

# 2.4 Plot bargraph
ad.m = melt(ad, id.vars =c("newOrder", "popId"), measure.vars = colnames(tmp[,2:(dim(tmp)[2]-2)]))
colnames(ad.m) = c("newOrder", "popId", "k", "p")

cols = colsFrame[1:length(grep("k[0-9]", colnames(ad)))] # select colors
names(cols) = colnames(ad)[grep("k[0-9]", colnames(ad))]
ad.m$p[which(ad.m$p==0)] = NA
head(ad.m)

p1 = ggplot() +
  geom_bar(ad.m, mapping = aes(x=newOrder, y=p, color=k, fill=k), stat = "identity", position = "stack") +
  geom_text(data = labs, aes(x=xcoord, y=-0.05, label=popName), angle = 45, size=3, hjust=1 ) +
  ylim(c(-0.3,1.1)) +
  scale_color_manual(values = cols) + 
  scale_fill_manual(values = cols) +
  coord_fixed(100, expand = TRUE) +
  theme(
    legend.position="none",
    axis.line=element_blank(),axis.text.x=element_blank(),
    axis.text.y=element_blank(),axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),plot.background=element_blank())
p1

####### K=6
# 1.1 Load data
colsFrame = c("#B11217", "#FBA081", "#009e60", "#fcd116", "#E42F28", "#FC684A")
#colsFrame = c("#B11217", "#FBA081", "#3a75c4", "#fcd116", "#E42F28", "#FC684A")

coord = read.table("popIndex.txt")
colnames(coord) = c("popId", "popName")

inF = "ag1000gCol_wilding.admixture.K6.txt"
ad = read.table(inF)
colnames(ad)[1] = c("popId")
colnames(ad)[2:dim(ad)[2]] = paste(rep("k", dim(ad)[2]-1), 1:(dim(ad)[2]-1), sep="")
head(ad)

# 2.2 Plot bargraph admixture
# 2.2.1 Reorder individu within each population to make barplot pretty
ad$ind = 1:dim(ad)[1]
reorder = c() # reorder barplot
for(i in unique(ad$popId)){
  tmp = ad[which(ad$popId==i),]
  if(dim(tmp)[1]>1){
    reorder = c(reorder, tmp$ind[do.call(order, rev(as.list(tmp[,2:(dim(tmp)[2]-1)])))] )
  } else{
    reorder = c(reorder, tmp$ind)
  }
}
length(reorder)
length(ad$popId)
ad = ad[match(reorder, ad$ind),]
ad$newOrder = 1:dim(ad)[1]
head(ad)

# 2.2.2 Add spacer betwen each pop
spacerSize = 5
tmp = ad
for(i in unique(tmp$popId)){
  print(i)
  tmp$newOrder[-which(tmp$popId<i)] = tmp$newOrder[-which(tmp$popId<i)] + spacerSize
}
ad = tmp

# 2.3 Make dataframe for population label
labs = matrix(data=NA, nrow = length(unique(ad$popId)), ncol = 3)
for(i in unique(ad$popId)){
  md = median(ad$newOrder[which(ad$popId==i)])
  labs[i, ] = c(md, i, -0.1)
}
colnames(labs) = c("xcoord", "popId", "ycoord")
labs = data.frame(labs)
labs = merge(labs, coord, by="popId")

# 2.4 Plot bargraph
ad.m = melt(ad, id.vars =c("newOrder", "popId"), measure.vars = colnames(tmp[,2:(dim(tmp)[2]-2)]))
colnames(ad.m) = c("newOrder", "popId", "k", "p")

cols = colsFrame[1:length(grep("k[0-9]", colnames(ad)))] # select colors
names(cols) = colnames(ad)[grep("k[0-9]", colnames(ad))]
head(ad.m)
ad.m$p[which(ad.m$p==0)] = NA

p2 = ggplot() +
  geom_bar(ad.m, mapping = aes(x=newOrder, y=p, color=k, fill=k), stat = "identity", position = "stack") +
  geom_text(data = labs, aes(x=xcoord, y=-0.05, label=popName), angle = 45, size=3, hjust=1 ) +
  ylim(c(-0.3,1.1)) +
  scale_color_manual(values = cols) + 
  scale_fill_manual(values = cols) +
  coord_fixed(100, expand = TRUE) +
  theme(
    legend.position="none",
    axis.line=element_blank(),axis.text.x=element_blank(),
    axis.text.y=element_blank(),axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),plot.background=element_blank())
p2

pdf("wilding.ag1000g.col.admix_barplot.pdf", width = conv_unit(200, "mm", "inch"), height = conv_unit(100, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(p1, p2, ncol = 1, nrow = 2)
dev.off()

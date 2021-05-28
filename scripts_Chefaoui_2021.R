


##%######################################################%##
#                                                          #
####          SCRIPTS FOR ANALYSIS AND FIGURES          ####
####    OF PAPER: Chefaoui (2021)Seasonal variations    ####
####       of waterbird ecological networks under       ####
####                different saltworks                 ####
####         management. Ecological Informatics         ####
#                                                          #
##%######################################################%##


# The code is presented as individual chunks:


##%######################################################%##
#                                                          #
####             "betalinkr_multi" analysis             ####
#                                                          #
##%######################################################%##

install.packages("bipartite") 


# we load data of bird abundance for each season as matrixes:

load("/home/rosa/Documents/Aves/Aves.R/kk_summer.RData")
load("/home/rosa/Documents/Aves/Aves.R/kk_Autumn.RData")
load("/home/rosa/Documents/Aves/Aves.R/kk_Winter.RData")
load("/home/rosa/Documents/Aves/Aves.R/kk_Spring.RData")

# vector of types of saltpans
load("KK_sites.rda") 


# 1. betalinkr _ more than two webs _ taking ALL SITES pooled:
betalinkr_multi(webs2array(abun2_summer, abun2_Autumn,abun2_Winter, abun2_Spring),partitioning="commondenom",  index="jaccard", binary=F)

#####################
## 2. betalinkr using types of saltworks separately:

# we prepare data:
tipos_Autumn_abun<- aggregate(x=abun2_Autumn, list(site_type$Type), FUN=sum)
tipos_Spring_abun<- aggregate(x=abun2_Spring, list(site_type$Type), FUN=sum)
tipos_Summer_abun<- aggregate(x=abun2_summer, list(site_type$Type), FUN=sum)
tipos_Winter_abun<- aggregate(x=abun2_Winter, list(site_type$Type), FUN=sum)

rownames(tipos_Autumn_abun)<- tipos_Autumn_abun[,1]
tipos_Autumn_abun<- tipos_Autumn_abun[, -1]
tipos_Autumn_abun<-as.matrix(tipos_Autumn_abun)

rownames(tipos_Spring_abun)<- tipos_Spring_abun[,1]
tipos_Spring_abun<- tipos_Spring_abun[, -1]
tipos_Spring_abun<-as.matrix(tipos_Spring_abun)

rownames(tipos_Summer_abun)<- tipos_Summer_abun[,1]
tipos_Summer_abun<- tipos_Summer_abun[, -1]
tipos_Summer_abun<-as.matrix(tipos_Summer_abun)

rownames(tipos_Winter_abun)<- tipos_Winter_abun[,1]
tipos_Winter_abun<- tipos_Winter_abun[, -1]
tipos_Winter_abun<-as.matrix(tipos_Winter_abun)

# we calculate betalinkr_multi:
betalinkr_multi(webs2array(tipos_Summer_abun, tipos_Autumn_abun, tipos_Winter_abun, tipos_Spring_abun),partitioning="commondenom",  index="jaccard", binary=F )


##%######################################################%##
#                                                          #
####           PDI: Paired Differences Index            ####
#                                                          #
##%######################################################%##
install.packages("bipartite") 


# we load data of bird abundance for each season as matrixes:

load("/home/rosa/Documents/Aves/Aves.R/kk_summer.RData")
load("/home/rosa/Documents/Aves/Aves.R/kk_Autumn.RData")
load("/home/rosa/Documents/Aves/Aves.R/kk_Winter.RData")
load("/home/rosa/Documents/Aves/Aves.R/kk_Spring.RData")

# vector of types of saltpans
load("KK_sites.rda") 


# we prepare data:
tipos_Autumn_abun<- aggregate(x=abun2_Autumn, list(site_type$Type), FUN=sum)
tipos_Spring_abun<- aggregate(x=abun2_Spring, list(site_type$Type), FUN=sum)
tipos_Summer_abun<- aggregate(x=abun2_summer, list(site_type$Type), FUN=sum)
tipos_Winter_abun<- aggregate(x=abun2_Winter, list(site_type$Type), FUN=sum)

rownames(tipos_Autumn_abun)<- tipos_Autumn_abun[,1]
tipos_Autumn_abun<- tipos_Autumn_abun[, -1]
tipos_Autumn_abun<-as.matrix(tipos_Autumn_abun)

rownames(tipos_Spring_abun)<- tipos_Spring_abun[,1]
tipos_Spring_abun<- tipos_Spring_abun[, -1]
tipos_Spring_abun<-as.matrix(tipos_Spring_abun)

rownames(tipos_Summer_abun)<- tipos_Summer_abun[,1]
tipos_Summer_abun<- tipos_Summer_abun[, -1]
tipos_Summer_abun<-as.matrix(tipos_Summer_abun)

rownames(tipos_Winter_abun)<- tipos_Winter_abun[,1]
tipos_Winter_abun<- tipos_Winter_abun[, -1]
tipos_Winter_abun<-as.matrix(tipos_Winter_abun)


## Firstly, we remove singletons (species with only one observation) from abundance data: 

# Summer season: we remove Netta_rufin, Galli_chlor, Ardea_ciner, Calid_ferru
tipos_Summer_abun_sin<- tipos_Summer_abun[,-c(15,16,18,19)]

# Winter season: we remove C.ale N.pha T.neb A.acu V.van
tipos_Winter_abun_sin<- tipos_Winter_abun[,-c(4,13,18,19,20)]

# Autumn season: we remove A.cly
tipos_Autumn_abun_sin<- tipos_Autumn_abun[,-c(16)]

# Spring season: we remove N.pha P.squ P.leu A.cly
tipos_Spring_abun_sin<- tipos_Spring_abun[,-c(15,16,19,20)]


## then, we calculate PDI: 

PDI_Summer<-PDI(tipos_Summer_abun_sin) 

PDI_Autumn<-PDI(tipos_Autumn_abun_sin) 

PDI_Winter<-PDI(tipos_Winter_abun_sin)

PDI_Spring<-PDI(tipos_Spring_abun_sin) 



##%######################################################%##
#                                                          #
####    BIPARTITE NETWORK VISUALIZATION with igraph     ####
#                                                          #
##%######################################################%##


# part of the code for igraph network has been adapted from Katya Ognyanova's code:
# https://kateto.net/network-visualization


install.packages("igraph") 

# As an example, we load data of Summer season as a matrix:

load("/home/rosa/Documents/Aves/Aves.R/kk_summer.RData") # Abundance data per site
load("KK_sites.rda") # vector of types of saltpans

## converting the raw data to a WEIGHTED igraph network object:
net_abun_summer <- graph_from_incidence_matrix(abun2_summer, weighted= T)

# we search for vertices:
V(net_abun_summer) # the first 52 vertices are sites, the rest are species

## colours for vertices according to abandoned and active sites (for sites), 
# and IUCN categories (for species):

# Active saltpan=1 : "cyan3"
# Abandoned saltpan=2 : "goldenrod"
# LC=3 : "green3"
# NT=4 : "olivedrab3"
# VU=5 : "yellow"
# EN=6 : "lightsalmon1"
# CR=7 : "red3"

# we create vector of colours for Summer season:
colors_summer<- append(as.factor(site_type$Type),c(3,3,3,5,5,5,5,
                                                   3,4,4,3,5,7,5,
                                                   5,3,3,3,5,5))

colors_summer<- as.factor(colors_summer)


###### EDGES COLORS #########
# Let's color the edges of the graph based on their source node color.
# We'll get the starting node for each edge with "ends()".
edge_st_summer <- ends(net_abun_summer, es=E(net_abun_summer), names=F)[,1]

## vector of active and inactive saltpans needed for edge.color:

active<-site_type[site_type$Type=="ACT",]
edge.col_summer<-ifelse (edge_st_summer %in% active$ID ,"Act", "Inact")
edge.col_summer<-as.factor(edge.col_summer)


plot(net_abun_summer,
     layout=layout.bipartite,
     edge.arrow.size=1,
     edge.arrow.width=2,
     vertex.frame.color= (c("cyan4", "darkgoldenrod", "green3"
                            ,"olivedrab3"
                            ,"yellow"
                            ,"red3")[colors_summer]),
     vertex.shape = (c("rectangle", "circle")[V(net_abun_summer)$type+1]),
     vertex.size=5,
     vertex.size2=23,
     vertex.color=(c("cyan3", "goldenrod", "green3"
                     ,"olivedrab3"
                     ,"yellow"
                     ,"red3")[colors_summer]),
     vertex.label=V(net_abun_summer)$name,
     vertex.label.family="Helvetica",
     vertex.label.color="black",
     vertex.label.cex=0.7,
     edge.lty=1,
     edge.color= (c("#00CDCDB4", "#DAA520B4") [edge.col_summer]),
     edge.width= E(net_abun_summer)$weight/7,
     margin=c(-0.1,0,-0.1,0), 
     asp=0.2) 
title("Summer",cex.main=1.4,col.main="bisque4", adj=0.1, line=0.1)


##%######################################################%##
#                                                          #
####          NETWORK VISUALIZATION: VISWEB      
#                                                          #
##%######################################################%##

install.packages("bipartite") 
# As an example, we load data of Summer season as a matrix:

load("/home/rosa/Documents/Aves/Aves.R/kk_summer.RData") # Abundance data per site

## Visweb for Summer season:

visweb(abun2_summer, circles=TRUE,type="nested", boxes=T,  labsize=2.5, circle.max=3,
       textsize=3, 
       text="interaction",
       circle.col="#00000064",	
       circle.min=0.3 ,	#minimal size of circles, use to rescale circles appropriately, default is 0.2
       outerbox.border="blue")


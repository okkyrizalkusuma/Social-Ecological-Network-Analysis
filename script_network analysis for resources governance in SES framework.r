source("dia.r")
source("community.r")

## Read model specification
data <- model.dia("Directed_Connectivity of RG.dia")
data[order(data$From), ] # CEK DUPLIKASI

## Examine unweighted adjacency matrix
A <- adjacency.matrix(data, labels = TRUE)
A

## igraph
library(igraph)
# weight=TRUE
inet <- graph.adjacency(t(A), mode = "directed") 

class(inet) # igraph
names(inet[2])
summary(inet)
gorder(inet)
gsize(inet)


## Add attributes to the network, vertices, or edges
V(inet)$name

V(inet)$Unit<-c("RA", "RA", "RA", "RA", "RA", "RA", "RA", "RA",
                "RG", "RG", "RG", "RG", "RG",
                "RS", "RS", "RS", "RS", "RS", "RS",
                "RU", "RU", "RU", "RU", "RU"
                )
				
table(V(inet)$Unit)

#V(net)$degree <- degree(net)

## layout
l <- layout_in_circle(inet)
l <- layout_with_lgl(inet)
 

# Generate colors based on Unit type:
colrs <- adjustcolor(c("skyblue", "orange", "yellow", "pink"), 1)
V(inet)$color <- colrs[as.factor(V(inet)$Unit)]

atribut <- as.data.frame(cbind(grup = V(inet)$Unit, nama = V(inet)$name))
atribut[order(atribut$grup), ]

##################################################################
##################################################################
dev.new(width=50,height=50,unit="in")

par(mar = c(3,0,3,0))
plot(inet, 
     edge.arrow.size=.1, 
     vertex.color=NULL, 
     vertex.size=15,
     vertex.label = NULL,
     vertex.frame.color="gray", 
     vertex.label.color="black",
     vertex.label.cex=0.8, 
     vertex.label.dist=0.7, 
     edge.curved=NULL,
     edge.color = "blue",
     #layout = MCoords, 
     edge.arrow.size = 0.3)

plot(inet,
     edge.arrow.size=.1,
     edge.color ='blue',
     main = 'Connectivity of Resource Governance (RG)',
     vertex.size=18,
     vertex.label.size=10,
     vertex.lebel.color="black",
     vertex.lebel.dist=NULL,
     vertex.lebel.cex=0.9,
     #layout=in_circle(inet))
     vertex.color=c("yellow", "orange")[1+(V(inet)$Unit=="RG")])


#legend("topright", c("RA","RG","RS","RU"), pch=21, title = "SES Components", col=colrs, pt.bg=colrs, pt.cex=2, cex=1, bty="n", ncol=1)

#title("Resources Governance Connectivity", 
	  #line = 1, sub = 'SEZ Tanjung Lesung, Banten Province - Indonesia')


#SAVE
dev.print(tiff,"Model Dasar.tiff",res=600, compression="lzw+p",height=7,width=7,units="in")

## Node degrees
bulat <-layout_with_lgl(inet)
deg <- degree(inet, mode="all")
sort(deg)

par(mar = c(3,0,3,0))	#satu-satu
plot(inet, layout = bulat,
     edge.arrow.size=.2,
     vertex.size=deg*1.2,
     vertex.label.cex=0.7,
     vertex.label.dist=0.3,  
     vertex.label.color="black")

dev.print(tiff,"Degree RG.tiff",res=600, compression="lzw+p",height=7,width=7,units="in")


#Hubs and authorities 
hs <- hub_score(inet, weights=NA)$vector
as <- authority_score(inet, weights=NA)$vector

par(mfrow = c(1,2), 
    mar = c(0,2,0,2),
	oma = c(1,1,0,1))

plot(inet, 
edge.arrow.size=.1, 
vertex.size=hs*20, 
vertex.label.dist=0.5,
vertex.label.color="black",
vertex.label.cex=1,
layout=l)

title(expression(bolditalic("Hubs")), 
sub = expression(bolditalic("The Size of Node Indicate The Magnitude of Hubs Value")),
cex.main = 1, cex.sub = 0.7,
line = -3)


plot(inet, 
edge.arrow.size=.1, 
vertex.size=as*20,
vertex.label.dist=0.5,
vertex.label.color="black", 
vertex.label.cex=1,
layout=l)

title(expression(bolditalic("Authorities")), 
sub = expression(bolditalic("The Size of Node Indicate The Magnitude of Authorities Value")),
cex.main = 1.3, cex.sub = 0.7,
line = -3)

#SAVE
dev.print(tiff,"Hubs_Auth.tiff",res=600,compression="lzw+p",height=8,width=8,units="in")


#Community detection and Dendogram
cluster_walktrap(jack)
jack <- graph.data.frame(data, directed = F)
cluster_fast_greedy(jack)
cfg <- cluster_fast_greedy(jack)
plot(cfg,
     jack,
     vertex.size = 13,
     vertex.label.cex = 0.8)

cfg <- cluster_fast_greedy(jack)

#SAVE
dev.print(tiff,"Community Detection.tiff",res=600,compression="lzw+p",height=5,width=10,units="in")

dendPlot(cfg, mode="hclust")

#SAVE
dev.print(tiff,"Dendogram.tiff",res=300,compression="lzw+p",height=5,width=5,units="in")





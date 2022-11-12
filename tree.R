library(tidyverse)
library("ape")
library("Biostrings")
library(adegenet)
library(pegas)
library(pheatmap)
library(phangorn)
library(stats)
library(ade4)
library(ggplot2)
library(factoextra)




rm(list = ls())
t <- read.csv("blib_sats.csv")



t <- t%>%arrange(sampleGroup)
groupt <- t%>%group_by(sampleGroup)
samplet <- slice_sample(groupt, n =5)
t <- t%>%select(name,sampleGroup,sample_number,A0,A1,A2,A3,A4,A5,A6,A7,A8,B0,B1,B2,B3,B4,B5,B6,B7,B8,)





t <- t%>%
  filter(str_length(A0) == 486)%>%
  filter(str_length(A1) == 486)%>%
  filter(str_length(A2) == 486)%>%
  filter(str_length(A3) == 486)%>%
  filter(str_length(A4) == 486)%>%
  filter(str_length(A5) == 486)%>%
  filter(str_length(A6) == 486)%>%
  filter(str_length(A7) == 486)%>%
  filter(str_length(A8) == 486)%>%
  filter(str_length(B0) == 486)%>%
  filter(str_length(B1) == 486)%>%
  filter(str_length(B2) == 486)%>%
  filter(str_length(B3) == 486)%>%
  filter(str_length(B4) == 486)%>%
  filter(str_length(B5) == 486)%>%
  filter(str_length(B6) == 486)%>%
  filter(str_length(B7) == 486)%>%
  filter(str_length(B8) == 486)






tdistinct <- t%>%distinct(A0,A1,A2,A3,A4,A5,A6,A7,A8,B0,B1,B2,B3,B4,B5,B6,B7,B8, .keep_all = TRUE)
#tdistinct <- t%>%distinct(name, .keep_all = TRUE)

tlast <- tdistinct%>%filter(sampleGroup >= max(sampleGroup)-7)
tlast <- tlast%>%unite(A0,A1,A2,A3,A4,A5,A6,A7,A8,B0,B1,B2,B3,B4,B5,B6,B7,B8, sep = "")


#distinctLoci <- tdistinct%>%select(A0:B8)





fullDNA <- DNAStringSet(c(tlast$A0))
names(fullDNA) <- paste0("sg",tlast$sampleGroup,"sn",tlast$sample_number)
fullDNA

fullBin <- as.DNAbin(fullDNA)












D <- dist.dna(as.matrix(fancybin), model = "TN93", as.matrix = TRUE)

#dLastLoci <- dist.gene(lastLoci, method = "pairwise")
#dDistinctLoci <- dist.gene(distinctLoci, method = "pairwise")

D.pca <- prcomp(D, scale = TRUE)
fviz_eig(D.pca)


fviz_pca_ind(D.pca,
             geom = "point",
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


temp <- as.data.frame(as.matrix(D))
table.paint(temp, cleg = 0, clabel.row=0.1, clabel.col=0.1)
heatmap(as.matrix(D))

fullBin2 <- as.phyDat(fullBin)
class(fullBin2)

tre.ini <- nj(dist.dna(fullBin,model="raw"))

tre.ini


parsimony(tre.ini, fullBin2)


tre.pars <- optim.parsimony(tre.ini, fullBin2)

plot(tre.pars, type = "phylogram", show.tip=TRUE, edge.width=0.5, cex= 0.6,use.edge.length = TRUE)
plot(root(tre.pars,7),type = "phylogram", show.tip=TRUE, edge.width=0.5, cex= 0.6,use.edge.length = TRUE)
plot(root(tre.pars,7),type = "cladogram", show.tip=TRUE, edge.width=0.5, cex= 0.6)
plot(root(tre.pars,7),type = "fan", show.tip=TRUE, edge.width=0.5, cex= 0.6)
plot(root(tre.pars,7),type = "unrooted", show.tip=TRUE, edge.width=0.5, cex= 0.6)
plot(root(tre.pars,7),type = "radial", show.tip=TRUE, edge.width=0.5, cex= 0.6)

pml(tre.ini, dna2, k=4)




tre <- bionj(D)

tre <- ladderize(tre)
tre2 <- bionj(dLastLoci)
tre2 <- ladderize(tre2)



plot(tre, show.tip.label  = TRUE, show.node.label = TRUE, use.edge.length = TRUE, edge.width = 0.5, node.width = 0.5, cex = 0.5, node.depth = 0.1, type = "phylogram")

plot(tre2, show.tip.label  = TRUE, show.node.label = TRUE, use.edge.length = TRUE, edge.width = 0.5, node.width = 0.5, cex = 0.5, node.depth = 0.1, type = "phylogram")


tre2 <- root(genetre, outgroup = 41)
tre2 <- ladderize(tre2)

tre <- root(tre, 41)
tre <- ladderize(tre)
  
plot(tre,type = "phylogram", show.tip=TRUE,use.edge.length = TRUE, edge.width=0.5, cex= 0.6)
plot(tre,type = "phylogram", show.tip=TRUE,use.edge.length = FALSE, edge.width=0.5, cex= 0.6)

plot(tre,type = "cladogram", show.tip=TRUE, use.edge.length = FALSE, edge.width=0.5, cex= 0.6)

plot(tre,type = "fan", show.tip=TRUE, use.edge.length = FALSE, edge.width=0.5, cex= 0.6)
plot(tre,type = "unrooted", show.tip=TRUE, use.edge.length = FALSE)
plot(tre,type = "radial", show.tip=TRUE, use.edge.length = FALSE, edge.width=0.5, cex= 0.6)



plot.phylo(tre2, show.tip.label  = TRUE, show.node.label = TRUE, use.edge.length = TRUE, edge.width = 1, node.width = 1, cex = 0.5)



x <- as.vector(D)
y <- as.vector(as.dist(cophenetic(tre2)))
plot(x, y, xlab="original pairwise distances", ylab="pairwise distances on the tree", main="Is NJ appropriate?", pch=20, col=transp("black",.1), cex=3)
abline(lm(y~x), col="red")

cor(x,y)^2

tre3 <- as.phylo(hclust(D,method="average"))
y <- as.vector(as.dist(cophenetic(tre3)))
plot(x, y, xlab="original pairwise distances", ylab="pairwise distances on the tree", main="Is UPGMA appropriate?", pch=20, col=transp("black",.1), cex=3)
abline(lm(y~x), col="red")

cor(x,y)^2

myBoots <- boot.phylo(tre2, as.matrix(lastA0), function(e) root(bionj(dist.dna(e, model = "TN93")),1))

myBoots

plot(tre2, show.tip=FALSE, edge.width=0.2)
title("NJ tree + bootstrap values")
nodelabels(myBoots)

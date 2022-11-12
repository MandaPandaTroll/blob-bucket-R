#Genome analysis November 2022

rm(list = ls())

{
  library(tidyverse)
  library(ape) 
  library(ggplot2)
  library(adegenet)
  library(Biostrings)
  library(insect)
  library(MSA2dist)
  library(pheatmap)
  library(PopGenome)
  library(DECIPHER)
  library(factoextra)
  library(pegas)
  library(phangorn)
  library(stats)
  library(ade4)



  rm(list=ls())
  gc()
}


{
  popData <- read_csv("Pop_data.csv")
  blib_raw <- read_csv("blib_sats.csv")%>%select(time,name, sampleGroup, sample_number, testA, testB) 
  
  blib_raw$name <- as.numeric(gsub("blib",'',blib_raw$name))
  
  blib_raw <- blib_raw%>%arrange(sampleGroup)
  
  blib_unique <- distinct(blib_raw, name, .keep_all = TRUE)
  
  blib_unique$name <- paste0(
    blib_unique$sampleGroup, '_',blib_unique$sample_number,"_",blib_unique$name)
  
  sum(duplicated(blib_raw$name))
  sum(duplicated(blib_unique$name))
  rm(blib_raw)
  
  blib_sample <- blib_unique%>%group_by(sampleGroup)%>%slice_sample(n = 5)
}

sample_groups_vector <- blib_unique%>%select(sampleGroup)
sample_groups_vector <- c(sample_groups_vector$sampleGroup)


blib_split <- blib_unique%>%select(testA,testB)%>%split( sample_groups_vector)
blib_split_list_A <- list()
blib_split_list_B <- list()

for(i in 0:length(blib_split)){
  u <- as.character(i-1)
  blib_split_list_A[i] <- blib_split[[u]][1]
  blib_split_list_B[i] <- blib_split[[u]][1]
}
AS <- ( lapply(blib_split_list_A, DNAStringSet))

thebin <- lapply(AS, as.DNAbin)

nucdivs_vector <- sapply(thebin,nuc.div)
nucdivs_df <- data.frame(time = seq(to = max(blib_unique$time), from = 50, by = 50),diversity = nucdivs_vector )
#nucdivs_df <- nucdivs_df[1:nrow(nucdivs_df)-1,]

nucdivs_df <- nucdivs_df%>%mutate(diversity_normalised = (diversity-min(diversity))/(max(diversity)-min(diversity)) )

nucdivs_df <- nucdivs_df%>%mutate(blibN = popData%>%filter(t %% 50 == 0)%>%select(blibN))
nucdivs_df <- nucdivs_df%>%mutate(blibN_normalised = (blibN-min(blibN))/(max(blibN)-min(blibN)))
nucdivs_g <- nucdivs_df%>%select(time,diversity_normalised,blibN_normalised)%>%gather(key = "variable", value = "value", -time)

ggplot(nucdivs_g,aes(x = time, y = value, colour = variable))+geom_point(, size = 0.1)+geom_smooth(span = 0.1)


library("ggpubr")





cor(nucdivs_df$blibN,nucdivs_df$diversity)
cor.test(nucdivs_df$diversity,nucdivs_df$blibN)


ggplot(nucdivs_g,aes(x = time, y = value, colour = variable))+geom_point(, size = 0.1)+geom_smooth(span = 0.1)

ggplot(nucdivs_df,aes(blibN,diversity))+geom_smooth(span = 0.1)

firstGroup <- blib_unique%>%filter(sampleGroup == 0)
firstGroupA_string <- DNAStringSet(firstGroup$testA)
firstGroupB_string <- DNAStringSet(firstGroup$testB)

sampleA <- DNAStringSet(blib_sample$testA)
names(sampleA) <- paste0("A",blib_sample$sampleGroup,blib_sample$sample_number)
BrowseSeqs(sampleA, highlight = 0)
sampleD <- dist.dna(as.DNAbin(sampleA, as.matrix = TRUE), model = "TN93")

sample_phylo <- bionj(sampleD)
sample_bin <- as.phyDat(as.DNAbin(sampleA))
sample.pars <- optim.parsimony(sample_phylo, sample_bin)

rootedsamplephylo <- compute.brlen(sample.pars)
rootedsamplephylo <- rtt(rootedsamplephylo,blib_sample$sampleGroup)

plot(ladderize(rootedsamplephylo), type = "phylogram", show.tip=TRUE, edge.width=0.1, cex= 0.1,use.edge.length = TRUE)

names(firstGroupA_string) <- paste0("A",firstGroup$name)
names(firstGroupB_string) <- paste0("B",firstGroup$name)
blib_last <- blib_unique%>%filter(sampleGroup >= max(blib_unique$sampleGroup -3))

blib_last$name <- paste0(blib_last$sampleGroup, '_', blib_last$sample_number)


A_string <- DNAStringSet(blib_last$testA)
B_string <- DNAStringSet(blib_last$testB)
names(A_string) <- paste0("A",blib_last$name)
names(B_string) <- paste0("B",blib_last$name)
A_string <- unique(A_string)
B_string <- unique(B_string)




blib_unique$sampleGroup <- as.factor(blib_unique$sampleGroup)

testt <- blib_unique %>%
  group_by(sampleGroup, testA)%>%group_walk(testA,DNAStringSetList)
  

A_string_full <- DNAStringSet(blib_unique$testA)
names(A_string_full) <- paste0("A",blib_unique$name)


A_bin_full <- as.DNAbin(unique(A_string_full))

segsites <- seg.sites(A_bin_full)

A_bin <- (as.DNAbin(A_string))
B_bin <- as.DNAbin(B_string)
AB_bin <- join(A_bin,B_bin)
nuc.div(A_bin)
#write.FASTA(A_bin, "A_bin.fasta")
#write.FASTA(B_bin, "B_bin.fasta")



A_D <- dist.dna(A_bin, model = "TN93", as.matrix = TRUE)
B_D <- dist.dna(B_bin, model = "TN93", as.matrix = TRUE)
AB_D <- dist.dna(AB_bin, model = "TN93", as.matrix = TRUE)




A_D.pca <- prcomp(A_D, scale = TRUE)
B_D.pca <- prcomp(B_D, scale = TRUE)
AB_D.pca <- prcomp(AB_D, scale = TRUE)

fviz_eig(A_D.pca)
fviz_eig(B_D.pca)
fviz_eig(AB_D.pca)

fviz_pca_ind(A_D.pca,
             geom = "point",
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_ind(B_D.pca,
             geom = "point",
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_ind(AB_D.pca,
             geom = "point",
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

temp_A <- as.data.frame(as.matrix(A_D))
temp_B <- as.data.frame(as.matrix(B_D))
temp_AB <- as.data.frame(as.matrix(AB_D))

table.paint(temp_A, cleg = 0, clabel.row=0.1, clabel.col=0.1)
table.paint(temp_B, cleg = 0, clabel.row=0.1, clabel.col=0.1)
table.paint(temp_AB, cleg = 0, clabel.row=0.1, clabel.col=0.1)

heatmap(as.matrix(A_D))
heatmap(as.matrix(B_D))
heatmap(as.matrix(AB_D))


A_bin2 <- as.phyDat(A_bin)
B_bin2 <- as.phyDat(B_bin)
AB_bin2 <- as.phyDat(AB_bin)
class(A_bin2)


treA.ini <- bionj(dist.dna(A_bin,model="TN93"))
treB.ini <- bionj(dist.dna(B_bin,model="TN93"))
treAB.ini <- bionj(dist.dna(AB_bin,model="TN93"))



tre.ini

parsimony(treA.ini, A_bin2)
parsimony(treB.ini, B_bin2)
parsimony(treAB.ini, AB_bin2)

treA.pars <- optim.parsimony(treA.ini, A_bin2)
treB.pars <- optim.parsimony(treB.ini, B_bin2)
treAB.pars <- optim.parsimony(treAB.ini, AB_bin2)



plot(treA.pars, type = "phylogram", show.tip=TRUE, edge.width=0.5, cex= 0.3,use.edge.length = TRUE)

chronA <- chronos(treA.ini, lambda = 1, model = "correlated", quiet = FALSE,
                   calibration = makeChronosCalib(treA.ini),
                   control = chronos.control())
plot(chronA, type = "phylogram", show.tip=TRUE, edge.width=0.5, cex= 0.3,use.edge.length = TRUE)

plot(ladderize(root(chronA, out = c("A95_298","A95_91","A95_47","A95_323","A95_13","A95_31","A95_149"))), type = "phylogram", show.tip=TRUE, edge.width=0.5, cex= 0.3,use.edge.length = FALSE)


rootre <- root(tre.pars, "Ag44s24")
rootre <- ladderize(rootre)
plot(rootre, type = "radial", show.tip=TRUE, edge.width=0.5, cex= 0.6,use.edge.length = FALSE)


chroney <- chronos(tre.ini, lambda = 1, model = "correlated", quiet = FALSE,
                   calibration = makeChronosCalib(tre.ini),
                   control = chronos.control())
plot(chroney)
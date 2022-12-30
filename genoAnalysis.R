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
  names(popData)[1] <- "time"
  
  blib_raw <- bind_rows(lapply(list.files(pattern = "blib_sats*"), read.csv))
  blib_raw <- distinct(blib_raw, name, .keep_all = TRUE)
  blibPheno <- read_csv("Blib_genetics.csv")
  blibPheno <- distinct(blibPheno, name, .keep_all = TRUE)
  #blibPheno <- blibPheno%>%select(-sample_number)
  
 
 
  blib_raw <- blib_raw%>%select(time,name, sampleGroup, sample_number, testA, testB)
  
  blib_unique <- distinct(blib_raw, name, .keep_all = TRUE)
  
  
  
  blib_unique <- blib_unique%>%arrange(sampleGroup)
  
  
  blib_unique <- blib_unique%>%filter(str_length(testA) == 4374 )%>%filter(str_length(testB) == 4374 )
  blib_unique$name <- paste0(
    blib_unique$sampleGroup, '_',blib_unique$sample_number,"_",blib_unique$name)
  
  sum(duplicated(blib_raw$name))
  sum(duplicated(blib_unique$name))
  
  
  blib_sample <- blib_unique%>%group_by(sampleGroup)%>%filter(n() >= 10)%>% slice_sample(n = 10)
  blib_sample$name <- paste0(
    blib_sample$sampleGroup, '_',blib_sample$sample_number)
  
  rm(blib_raw)
}

plot(blibPheno$generation,blibPheno$Heterozygosity)
{
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
  
blib_split_list_A <- blib_split_list_A[lapply(blib_split_list_A,length)>0]
blib_split_list_B <- blib_split_list_B[lapply(blib_split_list_B,length)>0]
  AS <- ( lapply(blib_split_list_A, DNAStringSet))
  
  thebin <- lapply(AS, as.DNAbin)

  nucdivs_vector <- sapply(thebin,nuc.div)
  nucdivs_df <- data.frame(time = seq(from = 50, to = 50*(length(nucdivs_vector)), by = 50),diversity = nucdivs_vector )
 # nucdivs_df <- nucdivs_df[1:nrow(nucdivs_df)-1,]
  
  nucdivs_df <- merge(popData[,1:2], nucdivs_df,by = "time")
  nucdivs_df <- nucdivs_df%>%mutate(diversity_normalised = (diversity-min(diversity))/(max(diversity)-min(diversity)) )
  
  #nucdivs_df <- nucdivs_df%>%mutate(blibN = popData%>%select(blibN))
  #nucdivs_df$blibN <- nucdivs_df$blibN$blibN
  nucdivs_df <- nucdivs_df%>%mutate(blibN_normalised = (blibN-min(blibN))/(max(blibN)-min(blibN)))
  nucdivs_df <- tibble(nucdivs_df)
  nucdivs_g <- nucdivs_df%>%select(time,diversity_normalised,blibN_normalised)%>%gather(key = "variable", value = "value", -time)
  
  }

blibPheno <- left_join(blibPheno,popData,by = "time")
ggplot(blibPheno, aes(time,generation))+geom_point(alpha = 0.3)+geom_smooth()

plot(blibPheno$time,blibPheno$generation)

gtimemodel <- lm(generation~time, data = blibPheno)
abline(gtimemodel, col = "red")


plot(nucdivs_df$blibN,nucdivs_df$diversity);abline(lm(diversity~blibN, data = nucdivs_df), col = "red")
#plot(lm(diversity~blibN, data = nucdivs_df))
plot(nucdivs_df$time,nucdivs_df$diversity)

ggplot(nucdivs_g,aes(x = time/(60), y = value, colour = variable))+geom_line()+geom_point( size = 0.3)

ggplot(nucdivs_df,aes(x = time/60))+geom_smooth(aes(y = diversity))+geom_point(aes(y = diversity))

#library("ggpubr")





cor(nucdivs_df$blibN,nucdivs_df$diversity)
cor.test(nucdivs_df$diversity,as.numeric(nucdivs_df$blibN))




ggplot(nucdivs_df,aes(blibN,diversity))+geom_smooth()+geom_point()

#firstGroup <- blib_unique%>%filter(sampleGroup == 0)
#firstGroupA_string <- DNAStringSet(firstGroup$testA)
#firstGroupB_string <- DNAStringSet(firstGroup$testB)


sampleA <- DNAStringSet(blib_sample$testA)

names(sampleA) <- paste0("A","_",blib_sample$name)
sampleA <- unique(sampleA)
BrowseSeqs(sampleA, highlight = 0)
sampleD <- dist.dna(as.DNAbin(sampleA, as.matrix = TRUE), model = "TN93")


sampleD.pca <- prcomp(sampleD, scale = TRUE)
fviz_eig(sampleD.pca)


fviz_pca_ind(sampleD.pca,
             geom = "point",
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

temp <- as.data.frame(as.matrix(sampleD))

table.paint(temp, cleg = 0, clabel.row=0.1, clabel.col=0.1)


sample_phylo <- bionj(sampleD)
sample_bin <- as.phyDat(as.DNAbin(sampleA))
sample.pars <- optim.parsimony(sample_phylo, sample_bin)

rootedsamplephylo <- compute.brlen(sample.pars)
rootedsamplephylo <- rtt(rootedsamplephylo,blib_sample$sampleGroup)

plot(ladderize(rootedsamplephylo), type = "phylogram", show.tip=TRUE, edge.width=0.1, cex= 0.1,use.edge.length = TRUE)

names(firstGroupA_string) <- paste0("A",firstGroup$name)
names(firstGroupB_string) <- paste0("B",firstGroup$name)





blib_last <- blib_unique%>%filter(sampleGroup >= max(blib_unique$sampleGroup))

blib_last$name <- paste0(blib_last$sampleGroup, '_', blib_last$sample_number)


A_string <- DNAStringSet(blib_last$testA)
B_string <- DNAStringSet(blib_last$testB)
names(A_string) <- paste0("A",blib_last$name)
names(B_string) <- paste0("B",blib_last$name)
A_string <- unique(A_string)
B_string <- unique(B_string)


BrowseSeqs(A_string, highlight = 0)




A_bin <- (as.DNAbin(A_string))
A_phyd <- as.phyDat(A_bin)
B_bin <- as.DNAbin(B_string)
AB_bin <- join(A_bin,B_bin)
AB_phyd <- as.phyDat(AB_bin)
nuc.div(AB_bin)


Dist_JC=dist.ml(A_phyd,  model = "JC69")
NJ =NJ(Dist_JC) 
plot(ladderize(NJ), type = "fan", show.tip=TRUE, edge.width=0.5, cex= 0.3,use.edge.length = TRUE)



#Calculate likelihood of NJ tree
ML_START=pml(NJ, A_phyd) 

#Create maximum-likelihood tree from NJ tree

ML_JC=optim.pml(ML_START,model="JC",  optNni = T) #The simple algorithm 



#root, ladderize and plot maximum-likelihood tree

ML_JC <- ladderize(ML_JC$tree)
ML_JC <- root(ML_JC,out = "A67_28")
plot(ML_JC,use.edge.length=T, type ="phylogram", show.tip=TRUE, edge.width=0.5, cex= 0.2) 

#The complex algorithm 
ML_GTR=optim.pml(ML_START,model="GTR",  optNni = T)  

ML_GTR <- root(ML_GTR$tree,out = "A67_160")

plot(ML_GTR,use.edge.length=T, type ="phylogram", show.tip=TRUE, edge.width=0.5, cex= 0.2) 
plot(ML_GTR,use.edge.length=F, type ="fan", show.tip=TRUE, edge.width=0.5, cex= 0.2) 



#Compare ML trees

AIC(ML_JC, ML_GTR)


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



treA.ini

parsimony(treA.ini, A_bin2)
parsimony(treB.ini, B_bin2)
parsimony(treAB.ini, AB_bin2)

treA.pars <- optim.parsimony(treA.ini, A_bin2)
treB.pars <- optim.parsimony(treB.ini, B_bin2)
treAB.pars <- optim.parsimony(treAB.ini, AB_bin2)



plot(ladderize(treA.pars), type = "phylogram", show.tip=TRUE, edge.width=0.5, cex= 0.3,use.edge.length = TRUE)

chronB <- chronos(treB.ini, lambda = 1, model = "correlated", quiet = FALSE,
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



#Spicy output code
{
  blib_raw$name <- paste0("G_",blib_raw$sampleGroup,"_S_",blib_raw$sample_number)
  
  
  library(Biostrings)
  library(ape)
  
  
  Haplo_A.dnastring <- DNAStringSet(blib_raw$testA)  
  Haplo_B.dnastring <- DNAStringSet(blib_raw$testB) 
  names(Haplo_A.dnastring) <- paste0("Ch_A_", blib_raw$name)
  names(Haplo_B.dnastring) <- paste0("Ch_B_", blib_raw$name)
  write.csv(blib_raw,"blib_raw_verylong.csv")
  rm(blib_raw)
  Haplo_A.dnabin <- as.DNAbin(Haplo_A.dnastring)
  Haplo_B.dnabin <- as.DNAbin(Haplo_B.dnastring)
  
  
  
  
  Haplos <-insect::join(Haplo_A.dnabin,Haplo_B.dnabin)
  
  write.FASTA(Haplo_A.dnabin, "blib_raw_dnabin_chA_verylong.fas")
  write.FASTA(Haplo_B.dnabin, "blib_raw_dnabin_chB_verylong.fas")
}

install.packages("tidyverse")
install.packages("Biostrings")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

BiocManager::install("ape")

BiocManager::install("DECIPHER")


if(!require(scales)){
  install.packages("scales", dependencies=TRUE)
  library(scales)
}



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
library(DECIPHER)

rm(list=ls())
gc()

blib_raw <- read_csv("blib_sats.csv")%>%select(time,name, sampleGroup, sample_number, testA, testB) 

blib_raw$name <- as.numeric(gsub("blib",'',blib_raw$name))

blib_raw <- blib_raw%>%arrange(sampleGroup)
  
  
  testA_df <- blib_raw%>%select(name,sampleGroup,sample_number,testA)
  testB_df <- blib_raw%>%select(name,sampleGroup,sample_number,testB)
  
  A_DNAstring <- DNAStringSet(c(testA_df$testA))
  B_DNAstring <- DNAStringSet(c(testB_df$testB))
  
  
  names(A_DNAstring) <- paste0("A", '_', testA_df$sampleGroup, '_', testA_df$sample_number)
  names(B_DNAstring) <- paste0("B", '_', testB_df$sampleGroup, '_', testB_df$sample_number)
  
  #names(A_DNAstring) <- paste0("A", "g",testA_df$sampleGroup,"s",testA_df$sample_number)
  #names(B_DNAstring) <- paste0("B", "g",testB_df$sampleGroup,"s",testB_df$sample_number)
  
  A_DNAstring <- unique(A_DNAstring)
  B_DNAstring <- unique(B_DNAstring)
  

 
 
 

  
  
pairwiseAlignment(A_DNAstring[1], A_DNAstring[5754 ],
                              patternQuality=PhredQuality(22L),
                              subjectQuality=PhredQuality(22L),
                              type="global",
                              substitutionMatrix=NULL, fuzzyMatrix=NULL,
                              gapOpening=10, gapExtension=4,
                              scoreOnly=FALSE)
  

alignall <- DNAMultipleAlignment(A_DNAstring)
alignall
detail(alignall)
consenso <- consensusMatrix(alignall, baseOnly = TRUE)

BrowseSeqs(A_DNAstring, highlight = 0)
tranA <- translate(A_DNAstring)
BrowseSeqs(tranA, highlight = 0)


  
  
  blib_filtered <- blib_raw%>%
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
  

  
  



rm(list=ls())
blib_raw <- read_csv("blib_sats.csv")

group_blib <- blib_raw%>%filter(sampleGroup %% 5 == 0)
group_blib <- group_blib%>%group_by(sampleGroup)
group_blib <- group_blib%>%filter(n() > 10)
group_blib <- slice_sample(group_blib, n = 10)

testA_df_sample <- group_blib%>%select(time,name,sampleGroup,sample_number,testA)
testB_df_sample <- group_blib%>%select(time,name,sampleGroup,sample_number,testB)



Asample_DNAstring <- DNAStringSet(c(testA_df_sample$testA))
Bsample_DNAstring <- DNAStringSet(c(testB_df_sample$testB))



names(Asample_DNAstring) <- paste0("A", '_', testA_df_sample$sampleGroup, '_', testA_df_sample$sample_number)
names(Bsample_DNAstring) <- paste0("B", '_', testB_df_sample$sampleGroup, '_', testB_df_sample$sample_number)

Asample_DNAstring <- unique(Asample_DNAstring)
Bsample_DNAstring <- unique(Bsample_DNAstring)



#ABsample <- DNAStringSet(c( A = Asample_DNAstring, B = Bsample_DNAstring))
#ABsample <- unique(ABsample)


BrowseSeqs(Asample_DNAstring, highlight = 0)
tranA <- translate(Asample_DNAstring)
BrowseSeqs(tranA, highlight = 0)

tranAB <- translate(ABsample)
BrowseSeqs(tranAB, highlight = 0)


fullBin = as.DNAbin(Asample_DNAstring)

D <- dist.dna(fullBin, model = "raw", as.matrix = TRUE)
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

tre.ini <- nj(dist.dna(fullBin,model="TN93"))

tre.ini
df_names <- data.frame("names" = Asample_DNAstring@ranges@NAMES)
match_df(testA_df_sample,df_names)
parsimony(tre.ini, fullBin2)

tre.pars <- optim.parsimony(tre.ini, fullBin2)

plotTreeTime(tre.pars,)

plot(tre.pars, type = "phylogram", show.tip=TRUE, edge.width=0.5, cex= 0.6,use.edge.length = TRUE)

rootre <- root(tre.pars, "Ag44s24")
rootre <- ladderize(rootre)
plot(rootre, type = "radial", show.tip=TRUE, edge.width=0.5, cex= 0.6,use.edge.length = FALSE)



rm(list=ls())
blib_raw <- read_csv("blib_sats.csv")


blib_last <- blib_raw%>%arrange(sampleGroup,sample_number)%>%filter(sampleGroup >= max(sampleGroup-3))



lastA <- blib_last%>%select(name,sampleGroup,sample_number,testA)
lastB <- blib_last%>%select(name,sampleGroup,sample_number,testB)


Alast_DNAstring <- DNAStringSet(c(lastA$testA))
Blast_DNAstring <- DNAStringSet(c(lastB$testB))

names(Alast_DNAstring) <- paste0("A", '_', lastA$sampleGroup, '_', lastA$sample_number)
names(Blast_DNAstring) <- paste0("B", '_', lastB$sampleGroup, '_', lastB$sample_number)

Alast_DNAstring <- unique(Alast_DNAstring)
Blast_DNAstring <- unique(Blast_DNAstring)
ABLast <- DNAStringSet(c(Alast_DNAstring,Blast_DNAstring))


# First use a fixed substitution matrix
mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
globalAlign <-
  pairwiseAlignment(Alast_DNAstring, Blast_DNAstring, substitutionMatrix = mat,
                    gapOpening = 5, gapExtension = 2)
localAlign <-
  pairwiseAlignment(Alast_DNAstring, Blast_DNAstring, type = "local", substitutionMatrix = mat,
                    gapOpening = 5, gapExtension = 2)
overlapAlign <-
  pairwiseAlignment(Alast_DNAstring, Blast_DNAstring, type = "overlap", substitutionMatrix = mat,
                    gapOpening = 5, gapExtension = 2)



consensusString(Alast_DNAstring, ambiguityMap=IUPAC_CODE_MAP,
                threshold=0.25, shift=0L, width=NULL)

tranA <- translate(Alast_DNAstring)
BrowseSeqs(tranA, highlight = 0)

AbinAA <- as.AAbin(tranA)

fullBin <- as.DNAbin(Alast_DNAstring)
findMutations(fullBin)

D <- dist.dna(fullBin, model = "TN93", as.matrix = TRUE)
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

tre.ini <- nj(dist.dna(fullBin,model="TN93"))

tre.ini


parsimony(tre.ini, fullBin2)

tre.pars <- optim.parsimony(tre.ini, fullBin2)

plot(tre.pars, type = "phylogram", show.tip=TRUE, edge.width=0.5, cex= 0.5,use.edge.length = TRUE, show.node.label = TRUE)

rootre <- root(tre.pars, "A_327_1")
rootre <- ladderize(rootre)
plot(rootre, type = "fan", show.tip=TRUE, edge.width=0.5, cex= 0.6,use.edge.length = TRUE)


Alast_AA <- as.matrix(translate(Alast_DNAstring))

Daa <- dist.gene(Alast_AA)

Daa.pca <- prcomp(Daa, scale = TRUE)
fviz_eig(Daa.pca)


fviz_pca_ind(Daa.pca,
             geom = "point",
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


temp <- as.data.frame(as.matrix(Daa))
table.paint(temp, cleg = 0, clabel.row=0.1, clabel.col=0.1)
heatmap(as.matrix(Daa))



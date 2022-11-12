if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


 

BiocManager::install("MSA2dist")
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
rm(list=ls())
gc()
}


{
blib_raw <- read_csv("blib_sats.csv")%>%select(time,name, sampleGroup, sample_number, testA, testB) 

blib_raw$name <- as.numeric(gsub("blib",'',blib_raw$name))

blib_raw <- blib_raw%>%arrange(sampleGroup)

blib_unique <- distinct(blib_raw, name, .keep_all = TRUE)

blib_unique$name <- paste0(
  blib_unique$sampleGroup, '_',blib_unique$sample_number,"_",blib_unique$name)

sum(duplicated(blib_raw$name))
sum(duplicated(blib_unique$name))
rm(blib_raw)
}

firstGroup <- blib_unique%>%filter(sampleGroup == 0)
firstGroupA_string <- DNAStringSet(firstGroup$testA)
firstGroupB_string <- DNAStringSet(firstGroup$testB)



names(firstGroupA_string) <- paste0("A",firstGroup$name)
names(firstGroupB_string) <- paste0("B",firstGroup$name)
blib_last <- blib_unique%>%filter(sampleGroup >= max(blib_unique$sampleGroup))

blib_last$name <- paste0(blib_last$sampleGroup, '_', blib_last$sample_number)


A_string <- DNAStringSet(blib_last$testA)
B_string <- DNAStringSet(blib_last$testB)
names(A_string) <- paste0(blib_last$name)
names(B_string) <- paste0(blib_last$name)
A_string <- unique(A_string)
B_string <- unique(B_string)



A_bin <- (as.DNAbin(A_string))
B_bin <- as.DNAbin(B_string)

write.FASTA(A_bin, "A_bin.fasta")
write.FASTA(B_bin, "B_bin.fasta")


{OG_blib_genes_map <- FindGenes(firstGroupAB_string, minGeneLength = 30, processors = NULL, showPlot = TRUE)
  OG_blib_genes_sequence <- (ExtractGenes(OG_blib_genes_map,firstGroupAB_string))
  BrowseSeqs(unique(OG_blib_genes_sequence))
  aa <- ExtractGenes(OG_blib_genes_map, firstGroupAB_string, type="AAStringSet")
  BrowseSeqs((aa))
  aa
  }









labels(A_bin)
A_genind <- DNAbin2genind(A_bin)
B_genind <- DNAbin2genind(B_bin)


pool(A_genind,B_genind)

samplegroupsizeVector <- group_size(blib_unique)

numGroupsVector <- (seq(0,max(blib_unique$sampleGroup), 1))
sampleGroupLabels <- rep(c(numGroupsVector), times = samplegroupsizeVector)

sampleGroup_Ind_df <- tibble(sampleGroup = sampleGroupLabels, individual = (blib_unique$name), chrom = '_')

sampleGroup_Ind_df <- rbind(sampleGroup_Ind_df,sampleGroup_Ind_df)
sampleGroup_Ind_df$chrom[1:(nrow(sampleGroup_Ind_df))/2] <- 'A'
sampleGroup_Ind_df$chrom[((nrow(sampleGroup_Ind_df)/2)+1):(nrow(sampleGroup_Ind_df))] <- 'B'



sampleGroup_Ind_df$sampleGroup <- as.factor(sampleGroup_Ind_df$sampleGroup)
sampleGroup_Ind_df$individual <- as.factor(sampleGroup_Ind_df$individual)




strata(A_genind) <- sampleGroup_Ind_df
strata(A_genind)
setPop(A_genind) <- ~sampleGroup
nPop(A_genind)
indNames(A_genind)
A_genpop <- genind2genpop(A_genind)

distmat_A_samplegroups <- dist.genpop(A_genpop)
heatmap(as.matrix(distmat_A_samplegroups))

library(pegas)
hw.test(A_genind)

Hs(A_genind)
Hs(B_genind)



rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggcorrplot)
library("PopGenome")
library(ape)
library("gridExtra")
popSizes <- read_csv("Pop_data.csv")



blib_harmon_vector <- 1/popSizes$blibN

sumharmon <- sum(blib_harmon_vector)

temp <- sumharmon/nrow(popSizes)

ne <- 1/temp

plot(popSizes$t, popSizes$blibN)

popSizes_scaled <- popSizes
popSizes_scaled$blibN <- (popSizes_scaled$blibN-min(popSizes_scaled$blibN))/(max(popSizes_scaled$blibN)-min(popSizes_scaled$blibN))


popSizes_scaled$blobN <- (popSizes_scaled$blobN-min(popSizes_scaled$blobN))/(max(popSizes_scaled$blobN)-min(popSizes_scaled$blobN))

popSizes_scaled$blybN <- (popSizes_scaled$blybN-min(popSizes_scaled$blybN))/(max(popSizes_scaled$blybN)-min(popSizes_scaled$blybN))


popSizes_scaled$meanHeterozygosity_blib <- (popSizes_scaled$meanHeterozygosity_blib-min(popSizes_scaled$meanHeterozygosity_blib))/(max(popSizes_scaled$meanHeterozygosity_blib)-min(popSizes_scaled$meanHeterozygosity_blib))


popSizes_scaled$meanHeterozygosity_blob <- (popSizes_scaled$meanHeterozygosity_blob-min(popSizes_scaled$meanHeterozygosity_blob))/(max(popSizes_scaled$meanHeterozygosity_blob)-min(popSizes_scaled$meanHeterozygosity_blob))


popSizes_scaled$meanHeterozygosity_blyb <- (popSizes_scaled$meanHeterozygosity_blyb-min(popSizes_scaled$meanHeterozygosity_blyb))/(max(popSizes_scaled$meanHeterozygosity_blyb)-min(popSizes_scaled$meanHeterozygosity_blyb))



popG <- popSizes %>%
  select(t, blibN, blobN, blybN, meanHeterozygosity_blib,meanHeterozygosity_blob,meanHeterozygosity_blyb) %>%
  gather(key = "species", value = "value", -t)

popG_scaled <- popSizes_scaled %>%
  select(t, blibN, blobN, blybN, meanHeterozygosity_blib,meanHeterozygosity_blob,meanHeterozygosity_blyb) %>%
  gather(key = "species", value = "value", -t)

ggplot(popG, aes(x = t, y = value)) + 
  geom_point(aes(color = species), size = 0.6)+
  scale_color_manual(values = c("green", "magenta","#FFCC66","forestgreen","purple","orange" ))+theme_grey()


blibHetG <- popSizes_scaled %>% select(t,blibN,meanHeterozygosity_blib)%>% gather(key = "key", value = "value", -t)

ggplot(blibHetG, aes(x = t, y = value)) +
  geom_point(aes(color = key), size = 0.6, alpha = 0.5)

plot()
justPops <- popSizes %>%
  select( blibN, blobN, blybN) 




corr <- cor(justPops)
ggcorrplot(corr)


blobPheno <- read_csv("Blob_genetics.csv")

blobPheno_byGenerations <- tibble(blobPheno)

blobPheno_byGenerations <- blobPheno_byGenerations[order(blobPheno_byGenerations$generation),]








blob_colors <- blobPheno_byGenerations%>%select(generation, redAllele1, redAllele2, greenAllele1, greenAllele2, blueAllele1, blueAllele2)



blob_colors <- blob_colors%>%mutate(bluePheno = (blueAllele1+blueAllele2)/2, greenPheno = (greenAllele1+greenAllele2)/2,redPheno = (redAllele1+redAllele2)/2)

blob_colorPhenos <- blob_colors%>%select(generation,bluePheno,greenPheno,redPheno)

blob_colorPhenos_g <- blob_colorPhenos%>%gather(key = "trait", value = "val", -generation)

ggplot(blob_colorPhenos_g, aes(x = generation, y = val))+geom_count(aes(color = trait))+scale_color_manual(values = c("blue","green","red" ))

blob_colors_g <- blob_colors%>%gather(key = "locus", value = "val", -generation)

ggplot(blob_colors_g, aes(x = generation, y = val))+geom_count(aes(color = locus, alpha = 0.01))+scale_color_manual(values = c("blue", "darkblue","green","forestgreen","red","red4" ))+theme_bw()


{
blob_scaled_traits <- blobPheno_byGenerations

blob_scaled_traits$moveAllele1 <- (blob_scaled_traits$moveAllele1-min(blob_scaled_traits$moveAllele1))/(max(blob_scaled_traits$moveAllele1)-min(blob_scaled_traits$moveAllele1))

blob_scaled_traits$moveAllele2 <- (blob_scaled_traits$moveAllele2-min(blob_scaled_traits$moveAllele2))/(max(blob_scaled_traits$moveAllele2)-min(blob_scaled_traits$moveAllele2))

blob_scaled_traits$redAllele1 <- (blob_scaled_traits$redAllele1-min(blob_scaled_traits$redAllele1))/(max(blob_scaled_traits$redAllele1)-min(blob_scaled_traits$redAllele1))

blob_scaled_traits$redAllele2 <- (blob_scaled_traits$redAllele2-min(blob_scaled_traits$redAllele2))/(max(blob_scaled_traits$redAllele2)-min(blob_scaled_traits$redAllele2))


#blob_scaled_traits$greenAllele1 <- (blob_scaled_traits$greenAllele1-min(blob_scaled_traits$greenAllele1))/(max(blob_scaled_traits$greenAllele1)-min(blob_scaled_traits$greenAllele1))


#blob_scaled_traits$greenAllele2 <- (blob_scaled_traits$greenAllele2-min(blob_scaled_traits$greenAllele2))/(max(blob_scaled_traits$greenAllele2)-min(blob_scaled_traits$greenAllele2))


blob_scaled_traits$blueAllele1 <- (blob_scaled_traits$blueAllele1-min(blob_scaled_traits$blueAllele1))/(max(blob_scaled_traits$blueAllele1)-min(blob_scaled_traits$blueAllele1))


blob_scaled_traits$blueAllele2 <- (blob_scaled_traits$blueAllele2-min(blob_scaled_traits$blueAllele2))/(max(blob_scaled_traits$blueAllele2)-min(blob_scaled_traits$blueAllele2))


blob_scaled_traits$lifeLengthAllele1 <- (blob_scaled_traits$lifeLengthAllele1-min(blob_scaled_traits$lifeLengthAllele1))/(max(blob_scaled_traits$lifeLengthAllele1)-min(blob_scaled_traits$lifeLengthAllele1))

blob_scaled_traits$lifeLengthAllele2 <- (blob_scaled_traits$lifeLengthAllele2-min(blob_scaled_traits$lifeLengthAllele2))/(max(blob_scaled_traits$lifeLengthAllele2)-min(blob_scaled_traits$lifeLengthAllele2))


blob_scaled_traits$lookDistAllele1 <- (blob_scaled_traits$lookDistAllele1-min(blob_scaled_traits$lookDistAllele1))/(max(blob_scaled_traits$lookDistAllele1)-min(blob_scaled_traits$lookDistAllele1))

blob_scaled_traits$lookDistAllele2 <- (blob_scaled_traits$lookDistAllele2-min(blob_scaled_traits$lookDistAllele2))/(max(blob_scaled_traits$lookDistAllele2)-min(blob_scaled_traits$lookDistAllele2))


blob_scaled_traits$turnTorqueAllele1 <- (blob_scaled_traits$turnTorqueAllele1-min(blob_scaled_traits$turnTorqueAllele1))/(max(blob_scaled_traits$turnTorqueAllele1)-min(blob_scaled_traits$turnTorqueAllele1))

blob_scaled_traits$turnTorqueAllele2 <- (blob_scaled_traits$turnTorqueAllele2-min(blob_scaled_traits$turnTorqueAllele2))/(max(blob_scaled_traits$turnTorqueAllele2)-min(blob_scaled_traits$turnTorqueAllele2))


blob_scaled_traits$e2repAllele1 <- (blob_scaled_traits$e2repAllele1-min(blob_scaled_traits$e2repAllele1))/(max(blob_scaled_traits$e2repAllele1)-min(blob_scaled_traits$e2repAllele1))

blob_scaled_traits$e2repAllele2 <- (blob_scaled_traits$e2repAllele2-min(blob_scaled_traits$e2repAllele2))/(max(blob_scaled_traits$e2repAllele2)-min(blob_scaled_traits$e2repAllele2))

blob_scaled_traits$sizeAllele1 <- (blob_scaled_traits$sizeAllele1-min(blob_scaled_traits$sizeAllele1))/(max(blob_scaled_traits$sizeAllele1)-min(blob_scaled_traits$sizeAllele1))

blob_scaled_traits$sizeAllele2 <- (blob_scaled_traits$sizeAllele2-min(blob_scaled_traits$sizeAllele2))/(max(blob_scaled_traits$sizeAllele2)-min(blob_scaled_traits$sizeAllele2))

blob_scaled_traits <- blob_scaled_traits%>%mutate(population_size = popSizes_scaled$blobN)
}

blob_scaled_traits_g <- blob_scaled_traits%>%select(generation, moveAllele1,moveAllele2,redAllele1,redAllele2,blueAllele1,blueAllele2,lifeLengthAllele1,lifeLengthAllele2,lookDistAllele1,lookDistAllele2,turnTorqueAllele1,turnTorqueAllele2,e2repAllele1,e2repAllele2,sizeAllele1,sizeAllele2)%>%gather(key = "locus", value = "val", -generation)

ggplot(blob_scaled_traits_g, aes(x = generation, y = val))+geom_smooth(aes(color = locus, alpha = 0.01))+theme_bw()




blybPheno <- read_csv("Blyb_genetics.csv")

blybPheno_byGenerations <- tibble(blybPheno)

blybPheno_byGenerations <- blybPheno_byGenerations[order(blybPheno_byGenerations$generation),]








blyb_colors <- blybPheno_byGenerations%>%select(generation, redAllele1, redAllele2, greenAllele1, greenAllele2, blueAllele1, blueAllele2)



blyb_colors <- blyb_colors%>%mutate(bluePheno = (blueAllele1+blueAllele2)/2, greenPheno = (greenAllele1+greenAllele2)/2,redPheno = (redAllele1+redAllele2)/2)

blyb_colorPhenos <- blyb_colors%>%select(generation,bluePheno,greenPheno,redPheno)

blyb_colorPhenos_g <- blyb_colorPhenos%>%gather(key = "trait", value = "val", -generation)

ggplot(blyb_colorPhenos_g, aes(x = generation, y = val))+geom_count(aes(color = trait))+scale_color_manual(values = c("blue","green","red" ))

blyb_colors_g <- blyb_colors%>%gather(key = "locus", value = "val", -generation)

ggplot(blyb_colors_g, aes(x = generation, y = val))+geom_count(aes(color = locus, alpha = 0.01))+scale_color_manual(values = c("blue", "darkblue","green","forestgreen","red","red4" ))+theme_bw()


{
  blyb_scaled_traits <- blybPheno
  
  blyb_scaled_traits$moveAllele1 <- blyb_scaled_traits$moveAllele1-blyb_scaled_traits$moveAllele1[1]
  
  blyb_scaled_traits$moveAllele2 <- blyb_scaled_traits$moveAllele2-blyb_scaled_traits$moveAllele2[1]
  
  blyb_scaled_traits$redAllele1 <- blyb_scaled_traits$redAllele1-blyb_scaled_traits$redAllele1[1]
  
  blyb_scaled_traits$redAllele2 <- blyb_scaled_traits$redAllele2-blyb_scaled_traits$redAllele2[1]
  
  
  blyb_scaled_traits$greenAllele1 <- blyb_scaled_traits$greenAllele1-blyb_scaled_traits$greenAllele1[1]
  
  
  blyb_scaled_traits$greenAllele2 <- blyb_scaled_traits$greenAllele2-blyb_scaled_traits$greenAllele2[1]
  
  
  blyb_scaled_traits$blueAllele1 <- blyb_scaled_traits$blueAllele1-blyb_scaled_traits$blueAllele1[1]
  
  
  blyb_scaled_traits$blueAllele2 <- blyb_scaled_traits$blueAllele2-blyb_scaled_traits$blueAllele2[1]
  
  
  blyb_scaled_traits$lifeLengthAllele1 <- blyb_scaled_traits$lifeLengthAllele1-blyb_scaled_traits$lifeLengthAllele1[1]
  
  blyb_scaled_traits$lifeLengthAllele2 <- blyb_scaled_traits$lifeLengthAllele2-blyb_scaled_traits$lifeLengthAllele2[1]
  
  
  blyb_scaled_traits$lookDistAllele1 <- blyb_scaled_traits$lookDistAllele1-blyb_scaled_traits$lookDistAllele1[1]
  
  blyb_scaled_traits$lookDistAllele2 <- blyb_scaled_traits$lookDistAllele2-blyb_scaled_traits$lookDistAllele2[1]
  
  
  blyb_scaled_traits$turnTorqueAllele1 <- blyb_scaled_traits$turnTorqueAllele1-blyb_scaled_traits$turnTorqueAllele1[1]
  
  blyb_scaled_traits$turnTorqueAllele2 <- blyb_scaled_traits$turnTorqueAllele2-blyb_scaled_traits$turnTorqueAllele2[1]
  
  
  blyb_scaled_traits$e2repAllele1 <- blyb_scaled_traits$e2repAllele1-blyb_scaled_traits$e2repAllele1[1]
  
  blyb_scaled_traits$e2repAllele2 <- blyb_scaled_traits$e2repAllele2-blyb_scaled_traits$e2repAllele2[1]
  
  blyb_scaled_traits$sizeAllele1 <- blyb_scaled_traits$sizeAllele1-blyb_scaled_traits$sizeAllele1[1]
  
  blyb_scaled_traits$sizeAllele2 <- 
    blyb_scaled_traits$sizeAllele2-blyb_scaled_traits$sizeAllele2[1]

  #OWO

  blyb_scaled_traits$moveAllele1 <- ((blyb_scaled_traits$moveAllele1-min(blyb_scaled_traits$moveAllele1))/(max(blyb_scaled_traits$moveAllele1)-min(blyb_scaled_traits$moveAllele1))) 
  
  blyb_scaled_traits$moveAllele2 <- ((blyb_scaled_traits$moveAllele2-min(blyb_scaled_traits$moveAllele2))/(max(blyb_scaled_traits$moveAllele2)-min(blyb_scaled_traits$moveAllele2))) 
  
  blyb_scaled_traits$redAllele1 <- ((blyb_scaled_traits$redAllele1-min(blyb_scaled_traits$redAllele1))/(max(blyb_scaled_traits$redAllele1)-min(blyb_scaled_traits$redAllele1))) 
  
  blyb_scaled_traits$redAllele2 <- ((blyb_scaled_traits$redAllele2-min(blyb_scaled_traits$redAllele2))/(max(blyb_scaled_traits$redAllele2)-min(blyb_scaled_traits$redAllele2))) 
  
  
blyb_scaled_traits$greenAllele1 <- ((blyb_scaled_traits$greenAllele1-min(blyb_scaled_traits$greenAllele1))/(max(blyb_scaled_traits$greenAllele1)-min(blyb_scaled_traits$greenAllele1))) 
  
  
blyb_scaled_traits$greenAllele2 <- ((blyb_scaled_traits$greenAllele2-min(blyb_scaled_traits$greenAllele2))/(max(blyb_scaled_traits$greenAllele2)-min(blyb_scaled_traits$greenAllele2))) 
  
  
  blyb_scaled_traits$blueAllele1 <- ((blyb_scaled_traits$blueAllele1-min(blyb_scaled_traits$blueAllele1))/(max(blyb_scaled_traits$blueAllele1)-min(blyb_scaled_traits$blueAllele1))) 
  
  
  blyb_scaled_traits$blueAllele2 <- ((blyb_scaled_traits$blueAllele2-min(blyb_scaled_traits$blueAllele2))/(max(blyb_scaled_traits$blueAllele2)-min(blyb_scaled_traits$blueAllele2))) 
  
  
  blyb_scaled_traits$lifeLengthAllele1 <- ((blyb_scaled_traits$lifeLengthAllele1-min(blyb_scaled_traits$lifeLengthAllele1))/(max(blyb_scaled_traits$lifeLengthAllele1)-min(blyb_scaled_traits$lifeLengthAllele1))) 
  
  blyb_scaled_traits$lifeLengthAllele2 <- ((blyb_scaled_traits$lifeLengthAllele2-min(blyb_scaled_traits$lifeLengthAllele2))/(max(blyb_scaled_traits$lifeLengthAllele2)-min(blyb_scaled_traits$lifeLengthAllele2))) 
  
  
  blyb_scaled_traits$lookDistAllele1 <- ((blyb_scaled_traits$lookDistAllele1-min(blyb_scaled_traits$lookDistAllele1))/(max(blyb_scaled_traits$lookDistAllele1)-min(blyb_scaled_traits$lookDistAllele1))) 
  
  blyb_scaled_traits$lookDistAllele2 <- ((blyb_scaled_traits$lookDistAllele2-min(blyb_scaled_traits$lookDistAllele2))/(max(blyb_scaled_traits$lookDistAllele2)-min(blyb_scaled_traits$lookDistAllele2))) 
  
  
  blyb_scaled_traits$turnTorqueAllele1 <- ((blyb_scaled_traits$turnTorqueAllele1-min(blyb_scaled_traits$turnTorqueAllele1))/(max(blyb_scaled_traits$turnTorqueAllele1)-min(blyb_scaled_traits$turnTorqueAllele1))) 
  
  blyb_scaled_traits$turnTorqueAllele2 <- ((blyb_scaled_traits$turnTorqueAllele2-min(blyb_scaled_traits$turnTorqueAllele2))/(max(blyb_scaled_traits$turnTorqueAllele2)-min(blyb_scaled_traits$turnTorqueAllele2))) 
  
  
  blyb_scaled_traits$e2repAllele1 <- ((blyb_scaled_traits$e2repAllele1-min(blyb_scaled_traits$e2repAllele1))/(max(blyb_scaled_traits$e2repAllele1)-min(blyb_scaled_traits$e2repAllele1))) 
  
  blyb_scaled_traits$e2repAllele2 <- ((blyb_scaled_traits$e2repAllele2-min(blyb_scaled_traits$e2repAllele2))/(max(blyb_scaled_traits$e2repAllele2)-min(blyb_scaled_traits$e2repAllele2))) 
  
  blyb_scaled_traits$sizeAllele1 <- ((blyb_scaled_traits$sizeAllele1-min(blyb_scaled_traits$sizeAllele1))/(max(blyb_scaled_traits$sizeAllele1)-min(blyb_scaled_traits$sizeAllele1))) 
  
  blyb_scaled_traits$sizeAllele2 <- 
    ((blyb_scaled_traits$sizeAllele2-min(blyb_scaled_traits$sizeAllele2))/(max(blyb_scaled_traits$sizeAllele2)-min(blyb_scaled_traits$sizeAllele2))) 
  
  
  blyb_scaled_traits[is.na(blyb_scaled_traits)] <- 1
  
  blyb_scaled_traits[,10:27] <- blyb_scaled_traits[,10:27] -1
  
}

blyb_scaled_traits_g <- blyb_scaled_traits%>%select(time, moveAllele1,moveAllele2,redAllele1,redAllele2,greenAllele1,greenAllele2,blueAllele1,blueAllele2,lifeLengthAllele1,lifeLengthAllele2,lookDistAllele1,lookDistAllele2,turnTorqueAllele1,turnTorqueAllele2,e2repAllele1,e2repAllele2,sizeAllele1,sizeAllele2)%>%gather(key = "locus", value = "val", -time)

traitplot <- ggplot(blyb_scaled_traits_g, aes(x = time, y = val))+geom_smooth(aes(color = locus), se = FALSE, span = 0.2)+theme_bw()

traitplot
popplot <- ggplot(popSizes_scaled, aes(x = t, y = blybN))+geom_point()



grid.arrange(traitplot, popplot, nrow = 2)

uniq <- unique(blyb_scaled_traits_g)


ggplot(uniq, aes(x = time, y = val))+geom_jitter(aes(color = locus))

rm(list = ls())
blibPheno <- read_csv("Blib_genetics.csv")

blibPheno_byGenerations <- tibble(blibPheno)

blibPheno_byGenerations <- blibPheno_byGenerations[order(blibPheno_byGenerations$generation),]








blib_colors <- blibPheno_byGenerations%>%select(generation, redAllele1, redAllele2, greenAllele1, greenAllele2, blueAllele1, blueAllele2)



blib_colors <- blib_colors%>%mutate(bluePheno = (blueAllele1+blueAllele2)/2, greenPheno = (greenAllele1+greenAllele2)/2,redPheno = (redAllele1+redAllele2)/2)

blib_colorPhenos <- blib_colors%>%select(generation,bluePheno,greenPheno,redPheno)

blib_colorPhenos_g <- blib_colorPhenos%>%gather(key = "trait", value = "val", -generation)

ggplot(blib_colorPhenos_g, aes(x = generation, y = val))+geom_count(aes(color = trait))+scale_color_manual(values = c("blue","green","red" ))

blib_colors_g <- blib_colors%>%gather(key = "locus", value = "val", -generation)

ggplot(blib_colors_g, aes(x = generation, y = val))+geom_count(aes(color = locus, alpha = 0.01))+scale_color_manual(values = c("blue", "darkblue","cyan","green","forestgreen","darkgreen","red","red4","darkred" ))+theme_bw()


{
  blib_scaled_traits <- blibPheno
  
  blib_scaled_traits$moveAllele1 <- blib_scaled_traits$moveAllele1-blib_scaled_traits$moveAllele1[1]
  
  blib_scaled_traits$moveAllele2 <- blib_scaled_traits$moveAllele2-blib_scaled_traits$moveAllele2[1]
  
  blib_scaled_traits$redAllele1 <- blib_scaled_traits$redAllele1-blib_scaled_traits$redAllele1[1]
  
  blib_scaled_traits$redAllele2 <- blib_scaled_traits$redAllele2-blib_scaled_traits$redAllele2[1]
  
  
  blib_scaled_traits$greenAllele1 <- blib_scaled_traits$greenAllele1-blib_scaled_traits$greenAllele1[1]
  
  
  blib_scaled_traits$greenAllele2 <- blib_scaled_traits$greenAllele2-blib_scaled_traits$greenAllele2[1]
  
  
  blib_scaled_traits$blueAllele1 <- blib_scaled_traits$blueAllele1-blib_scaled_traits$blueAllele1[1]
  
  
  blib_scaled_traits$blueAllele2 <- blib_scaled_traits$blueAllele2-blib_scaled_traits$blueAllele2[1]
  
  
  blib_scaled_traits$lifeLengthAllele1 <- blib_scaled_traits$lifeLengthAllele1-blib_scaled_traits$lifeLengthAllele1[1]
  
  blib_scaled_traits$lifeLengthAllele2 <- blib_scaled_traits$lifeLengthAllele2-blib_scaled_traits$lifeLengthAllele2[1]
  
  
  
  blib_scaled_traits$turnTorqueAllele1 <- blib_scaled_traits$turnTorqueAllele1-blib_scaled_traits$turnTorqueAllele1[1]
  
  blib_scaled_traits$turnTorqueAllele2 <- blib_scaled_traits$turnTorqueAllele2-blib_scaled_traits$turnTorqueAllele2[1]
  
  
  blib_scaled_traits$e2repA <- blib_scaled_traits$e2repA-blib_scaled_traits$e2repA[1]
  
  blib_scaled_traits$e2repB <- blib_scaled_traits$e2repB-blib_scaled_traits$e2repB[1]
  
  
  
  #OWO
  
  blib_scaled_traits$moveAllele1 <- ((blib_scaled_traits$moveAllele1-min(blib_scaled_traits$moveAllele1))/(max(blib_scaled_traits$moveAllele1)-min(blib_scaled_traits$moveAllele1))) 
  
  blib_scaled_traits$moveAllele2 <- ((blib_scaled_traits$moveAllele2-min(blib_scaled_traits$moveAllele2))/(max(blib_scaled_traits$moveAllele2)-min(blib_scaled_traits$moveAllele2))) 
  
  blib_scaled_traits$redAllele1 <- ((blib_scaled_traits$redAllele1-min(blib_scaled_traits$redAllele1))/(max(blib_scaled_traits$redAllele1)-min(blib_scaled_traits$redAllele1))) 
  
  blib_scaled_traits$redAllele2 <- ((blib_scaled_traits$redAllele2-min(blib_scaled_traits$redAllele2))/(max(blib_scaled_traits$redAllele2)-min(blib_scaled_traits$redAllele2))) 
  
  
  blib_scaled_traits$greenAllele1 <- ((blib_scaled_traits$greenAllele1-min(blib_scaled_traits$greenAllele1))/(max(blib_scaled_traits$greenAllele1)-min(blib_scaled_traits$greenAllele1))) 
  
  
  blib_scaled_traits$greenAllele2 <- ((blib_scaled_traits$greenAllele2-min(blib_scaled_traits$greenAllele2))/(max(blib_scaled_traits$greenAllele2)-min(blib_scaled_traits$greenAllele2))) 
  
  
  blib_scaled_traits$blueAllele1 <- ((blib_scaled_traits$blueAllele1-min(blib_scaled_traits$blueAllele1))/(max(blib_scaled_traits$blueAllele1)-min(blib_scaled_traits$blueAllele1))) 
  
  
  blib_scaled_traits$blueAllele2 <- ((blib_scaled_traits$blueAllele2-min(blib_scaled_traits$blueAllele2))/(max(blib_scaled_traits$blueAllele2)-min(blib_scaled_traits$blueAllele2))) 
  
  
  blib_scaled_traits$lifeLengthAllele1 <- ((blib_scaled_traits$lifeLengthAllele1-min(blib_scaled_traits$lifeLengthAllele1))/(max(blib_scaled_traits$lifeLengthAllele1)-min(blib_scaled_traits$lifeLengthAllele1))) 
  
  blib_scaled_traits$lifeLengthAllele2 <- ((blib_scaled_traits$lifeLengthAllele2-min(blib_scaled_traits$lifeLengthAllele2))/(max(blib_scaled_traits$lifeLengthAllele2)-min(blib_scaled_traits$lifeLengthAllele2))) 
  
  
  blib_scaled_traits$turnTorqueAllele1 <- ((blib_scaled_traits$turnTorqueAllele1-min(blib_scaled_traits$turnTorqueAllele1))/(max(blib_scaled_traits$turnTorqueAllele1)-min(blib_scaled_traits$turnTorqueAllele1))) 
  
  blib_scaled_traits$turnTorqueAllele2 <- ((blib_scaled_traits$turnTorqueAllele2-min(blib_scaled_traits$turnTorqueAllele2))/(max(blib_scaled_traits$turnTorqueAllele2)-min(blib_scaled_traits$turnTorqueAllele2))) 
  
  
  blib_scaled_traits$e2repA <- ((blib_scaled_traits$e2repA-min(blib_scaled_traits$e2repA))/(max(blib_scaled_traits$e2repA)-min(blib_scaled_traits$e2repA))) 
  
  blib_scaled_traits$e2repB <- ((blib_scaled_traits$e2repB-min(blib_scaled_traits$e2repB))/(max(blib_scaled_traits$e2repB)-min(blib_scaled_traits$e2repB))) 
  
  
  
  
  
  
 # blib_scaled_traits[is.na(blib_scaled_traits)] <- 1
  
 # blib_scaled_traits[,10:27] <- blib_scaled_traits[,10:27] -1
  
}

blib_scaled_traits_g <- blib_scaled_traits%>%select(time, moveAllele1,moveAllele2,redAllele1,redAllele2,greenAllele1,greenAllele2,blueAllele1,blueAllele2,lifeLengthAllele1,lifeLengthAllele2,turnTorqueAllele1,turnTorqueAllele2,e2repA,e2repB)%>%gather(key = "locus", value = "val", -time)

traitplot <- ggplot(blib_scaled_traits_g, aes(x = time, y = val))+geom_smooth(aes(color = locus), se = FALSE, span = 0.2)+theme_bw()

traitplot

ggplot(blib_scaled_traits_g, aes(x = time, y = val))+geom_point(aes(color = locus))+theme_bw()

popplot <- ggplot(popSizes_scaled, aes(x = t, y = blibN))+geom_point()



grid.arrange(traitplot, popplot, nrow = 2)

uniq <- unique(blib_scaled_traits_g)


ggplot(uniq, aes(x = time, y = val))+geom_jitter(aes(color = locus))

nameVector <- c('sampleGroup', 'moveAllele', 'redAllele', 'greenAllele', 'blueAllele', 'lifeLengthAllele', 'turnTorqueAllele', 'e2rep')

blibPheno_haplo1 <- blibPheno%>%select(sampleGroup, moveAllele1, redAllele1, greenAllele1, blueAllele1, lifeLengthAllele1, turnTorqueAllele1, e2repA)

names(blibPheno_haplo1) <- nameVector

blibPheno_haplo2 <- blibPheno%>%select(sampleGroup, moveAllele2, redAllele2, greenAllele2, blueAllele2, lifeLengthAllele2, turnTorqueAllele2, e2repB)

names(blibPheno_haplo2) <- nameVector

blibPheno_haplo <- rbind(blibPheno_haplo1, blibPheno_haplo2)

rm(blibPheno_haplo1)
rm(blibPheno_haplo2)

blibPheno_haplo <- blibPheno_haplo %>% arrange(sampleGroup)

ggplot(blibPheno_haplo, aes(x = turnTorqueAllele))+geom_histogram()


blib_last <- blibPheno%>%filter(sampleGroup >= max(blibPheno$sampleGroup))
blib_last <- blib_last%>%mutate(id = seq(1:nrow(blib_last)))


{
blib_last$moveAllele1 <- blib_last$moveAllele1/max(blib_last$moveAllele1)
blib_last$moveAllele2 <- blib_last$moveAllele2/max(blib_last$moveAllele2)
blib_last$redAllele1 <- blib_last$redAllele1/max(blib_last$redAllele1)
blib_last$redAllele2 <- blib_last$redAllele2/max(blib_last$redAllele2)
blib_last$greenAllele1 <- blib_last$greenAllele1/max(blib_last$greenAllele1)
blib_last$greenAllele2 <- blib_last$greenAllele2/max(blib_last$greenAllele2)
blib_last$blueAllele1 <- blib_last$blueAllele1/max(blib_last$blueAllele1)
blib_last$blueAllele2 <- blib_last$blueAllele2/max(blib_last$blueAllele2)
blib_last$turnTorqueAllele1 <- blib_last$turnTorqueAllele1/max(blib_last$turnTorqueAllele1)
blib_last$turnTorqueAllele2 <- blib_last$turnTorqueAllele2/max(blib_last$turnTorqueAllele2)
blib_last$e2repA <- blib_last$e2repA/max(blib_last$e2repA)
blib_last$e2repB <- blib_last$e2repB/max(blib_last$e2repB)
}




blib_last_g <- blib_last%>%select(id, moveAllele1,moveAllele2,redAllele1,redAllele2,greenAllele1,greenAllele2,blueAllele1,blueAllele2,lifeLengthAllele1,lifeLengthAllele2,turnTorqueAllele1,turnTorqueAllele2,e2repA,e2repB)%>%gather(key = "locus", value = "val", -id)

ggplot(blib_last_g, aes(x = val, color = locus))+ geom_histogram()
ggplot(blib_last, aes(x = moveAllele1))+geom_histogram()

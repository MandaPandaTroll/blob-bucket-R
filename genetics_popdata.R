rm(list = ls())
install.packages("tidyverse")
install.packages("ggplot2")
library(tidyverse)
library(ggplot2)
library(PopGenome)

blib_pheno <- read_csv("Blib_genetics.csv")
blob_pheno <- read_csv("Blob_genetics.csv")
blyb_pheno <- read_csv("Blyb_genetics.csv")
blub_pheno <- read_csv("Blub_genetics.csv")
popSizes <- read_csv("Pop_data.csv")

blob_pheno <- blob_pheno%>%distinct()
blib_pheno <- blib_pheno%>%distinct()







popSizes_scaled <- popSizes
popSizes_scaled$blibN <- (popSizes_scaled$blibN-min(popSizes_scaled$blibN))/(max(popSizes_scaled$blibN)-min(popSizes_scaled$blibN))


popSizes_scaled$blobN <- (popSizes_scaled$blobN-min(popSizes_scaled$blobN))/(max(popSizes_scaled$blobN)-min(popSizes_scaled$blobN))

popSizes_scaled$blybN <- (popSizes_scaled$blybN-min(popSizes_scaled$blybN))/(max(popSizes_scaled$blybN)-min(popSizes_scaled$blybN))


popSizes_scaled$blubN <- (popSizes_scaled$blubN-min(popSizes_scaled$blubN))/(max(popSizes_scaled$blubN)-min(popSizes_scaled$blubN))


popSizesG <- popSizes %>%
  select(t, blibN, blobN, blybN, blubN) %>%
  gather(key = "species", value = "individuals", -t)



ggplot(popSizesG, aes(x = t, y = log10(individuals))) + 
  geom_line(aes(color = species))+
  scale_color_manual(values = c("forestgreen", "magenta","red","#FFCC66"))+theme_bw()




popSizes_scaledG <- popSizes_scaled %>%
  select(t, blibN, blobN, blybN, blubN) %>%
  gather(key = "species", value = "individuals", -t)



ggplot(popSizes_scaledG, aes(x = t, y = individuals)) + 
  geom_line(aes(color = species))+
  scale_color_manual(values = c("green", "magenta","red","#FFCC66"))+theme_bw()


ggplot(popSizes_scaledG, aes(x = t, y = individuals)) + 
  geom_smooth(aes(color = species), method = "loess", se = TRUE, span = 0.1)+
  scale_color_manual(values = c("green", "magenta","red","#FFCC66"))+theme_bw()

popcov <- cov(popSizes_scaled)










bloybdif <- anti_join(blob_pheno, blyb_pheno)



  

blob_pheno <- blob_pheno%>%mutate(sample = seq(1:nrow(blob_pheno)))
blob_pheno <- blob_pheno%>%mutate(
  moveForce = (moveAllele1+moveAllele2)/2)
blob_pheno <- blob_pheno%>%mutate(
  e2rep = (e2repAllele1+e2repAllele2)/2)


moveAlleles_blob <- blob_pheno %>%
  select(sample, moveAllele1, moveAllele2) %>%
  gather(key = "allele", value = "value", -sample)

ggplot(moveAlleles_blob, aes(x = sample, y = value)) + 
  geom_smooth(aes(color = allele), method = "loess")+
  scale_color_manual(values = c("blue", "red"))+labs(title = "Blob MoveAlleles")

moveGG <- ggplot(moveAlleles_blob)
moveGG

moveGG+ geom_point(aes(sample,value, color = allele), alpha = 0.1)


e2reps_blob <- blob_pheno %>%
  select(sample,  e2rep) %>%
  gather(key = "allele", value = "value", -sample)

e2repGG <- ggplot(e2reps_blob)
e2repGG+ geom_point(aes(sample,value, color = allele), alpha = 0.1)

ggplot(e2reps_blob, aes(x = sample, y = value)) + 
  geom_smooth(aes(color = allele), method = "loess")+
  scale_color_manual(values = c("blue", "red"))


turnAlleles_blob <- blob_pheno %>%
  select(sample, turnTorqueAllele1, turnTorqueAllele2) %>%
  gather(key = "allele", value = "value", -sample)

ggplot(turnAlleles_blob, aes(x = sample, y = value)) + 
  geom_smooth(aes(color = allele), method = "loess")+
  scale_color_manual(values = c("blue", "red"))+labs(title = "Blob turnTorqueAlleles")


sizeAlleles_blob <- blob_pheno %>%
  select(sample, sizeAllele1, sizeAllele2) %>%
  gather(key = "allele", value = "value", -sample)

ggplot(sizeAlleles_blob, aes(x = sample, y = value)) + 
  geom_smooth(aes(color = allele), method = "loess")+
  scale_color_manual(values = c("blue", "red"))+labs(title = "Blob sizeAlleles")


lifeAlleles_blob <- blob_pheno %>%
  select(sample, lifeLengthAllele1, lifeLengthAllele2) %>%
  gather(key = "allele", value = "value", -sample)

ggplot(lifeAlleles_blob, aes(x = sample, y = value)) + 
  geom_smooth(aes(color = allele), method = "loess")+
  scale_color_manual(values = c("blue", "red"))+labs(title = "Blob lifeLengthAlleles")
















blyb_pheno <- blyb_pheno%>%mutate(sample = seq(1:nrow(blyb_pheno)))
blyb_pheno <- blyb_pheno%>%mutate(
  moveForce = (moveAllele1+moveAllele2)/2)
blyb_pheno <- blyb_pheno%>%mutate(
  e2rep = (e2repAllele1+e2repAllele2)/2)


moveAlleles_blyb <- blyb_pheno %>%
  select(sample, moveAllele1, moveAllele2) %>%
  gather(key = "allele", value = "value", -sample)

ggplot(moveAlleles_blyb, aes(x = sample, y = value)) + 
  geom_smooth(aes(color = allele), method = "loess")+
  scale_color_manual(values = c("blue", "red"))+labs(title = "blyb MoveAlleles")

moveGG_blyb <- ggplot(moveAlleles_blyb)


moveGG_blyb+ geom_point(aes(sample,value, color = allele), alpha = 0.1)


e2reps_blyb <- blyb_pheno %>%
  select(sample,  e2rep) %>%
  gather(key = "allele", value = "value", -sample)

e2repGG <- ggplot(e2reps_blyb)
e2repGG+ geom_point(aes(sample,value, color = allele), alpha = 0.1)

ggplot(e2reps_blyb, aes(x = sample, y = value)) + 
  geom_smooth(aes(color = allele), method = "loess")+
  scale_color_manual(values = c("blue", "red"))


turnAlleles_blyb <- blyb_pheno %>%
  select(time, turnTorqueAllele1, turnTorqueAllele2) %>%
  gather(key = "allele", value = "value", -time)

ggplot(turnAlleles_blyb, aes(x = time, y = value)) + 
  geom_smooth(aes(color = allele), method = "loess")+
  scale_color_manual(values = c("blue", "red"))+labs(title = "blyb turnTorqueAlleles")


sizeAlleles_blyb <- blyb_pheno %>%
  select(time, sizeAllele1, sizeAllele2) %>%
  gather(key = "allele", value = "value", -time)

ggplot(sizeAlleles_blyb, aes(x = time, y = value)) + 
  geom_smooth(aes(color = allele), method = "loess")+
  scale_color_manual(values = c("blue", "red"))+labs(title = "blyb sizeAlleles")


lifeAlleles_blyb <- blyb_pheno %>%
  select(sample, lifeLengthAllele1, lifeLengthAllele2) %>%
  gather(key = "allele", value = "value", -sample)

ggplot(lifeAlleles_blyb, aes(x = sample, y = value)) + 
  geom_smooth(aes(color = allele), method = "loess")+
  scale_color_manual(values = c("blue", "red"))+labs(title = "blyb lifeLengthAlleles")












rm(list=ls())
blib_pheno <- read_csv("Blib_genetics.csv")



samples_blib <- seq(1:nrow(blib_pheno))
blib_pheno <- blib_pheno%>%mutate(sample = samples_blib)
blib_pheno <- blib_pheno%>%mutate(
  moveForce = (moveAllele1+moveAllele2)/2)

moveAlleles_blib <- blib_pheno %>%
  select(sample, moveAllele1, moveAllele2, moveForce) %>%
  gather(key = "allele", value = "value", -sample)

ggplot(moveAlleles_blib, aes(x = sample, y = value)) + 
  geom_smooth(aes(color = allele), method = "loess", se = FALSE)+
  scale_color_manual(values = c("blue", "red", "green"))+labs(title = "Blib MoveAlleles")

ggplot(moveAlleles_blib, aes(x = sample, y = value)) + 
  geom_point(aes(color = allele, alpha = 0.001))+
  scale_color_manual(values = c("blue", "red", "green"))+labs(title = "Blib MoveAlleles")

turnAlleles_blib <- blib_pheno %>%
  select(sample, turnTorqueAllele1, turnTorqueAllele2) %>%
  gather(key = "allele", value = "value", -sample)

ggplot(turnAlleles_blib, aes(x = sample, y = value)) + 
  geom_smooth(aes(color = allele), method = "loess")+
  scale_color_manual(values = c("blue", "red"))

greenAlleles_blib <- blib_pheno %>%
  select(sample, greenAllele1, greenAllele2) %>%
  gather(key = "allele", value = "value", -sample)

ggplot(greenAlleles_blib, aes(x = sample, y = value)) + 
  geom_smooth(aes(color = allele), method = "loess")+geom_point(aes(color=allele))
  scale_color_manual(values = c("blue", "red"))

blib_pheno <- blib_pheno%>%arrange(generation)


blobtimes <- blob_pheno%>%distinct(time)
poptimes <- popSizes%>%distinct(t)

fixtime <- function(x,y) if (x != y ){x = y} 

niceSeq <- data.frame(t = seq(50,11500, 50))


blob_pheno <- blob_pheno%>%rename(t = time)

blob_pheno$t <- plyr::round_any(blob_pheno$t, 50)

testjoin <- full_join(blob_pheno,popSizes)

testjoin <- testjoin%>%filter(species == "blob")



test_moveAlleles_blob <- testjoin %>%
  select(t, moveAllele1, moveAllele2,blobN) %>%
  gather(key = "allele", value = "value", -t)

ggplot(test_moveAlleles_blob, aes(x = t, y = (value))) + 
  geom_smooth(aes(color = allele), method = "loess", se = FALSE)+
  scale_color_manual(values = c("blue", "red", "magenta"))+labs(title = "Blob MoveAlleles")

ggplot(popSizes, aes(t, blobN))+geom_line()
  

summary(blob_pheno)
summary(blyb_pheno)



tripops <- popSizes%>%mutate(bloybN = blobN+blybN)
tripops <- tripops%>%select(t, blibN, bloybN, blubN)

tripops <- tripops%>%rename("Producer" = blibN)%>%rename("Primary_Consumer" = bloybN)%>%rename("Secondary_Consumer" = blubN)%>%rename(Time_Steps = t)

tripopsG <- tripops %>%
  select("Time_Steps", "Producer", "Primary_Consumer", "Secondary_Consumer") %>%
  gather(key = "species", value = "individuals", -"Time_Steps")



ggplot(tripopsG, aes(x = Time_Steps, y = log10(individuals))) + 
  geom_line(aes(color = species))+
  scale_color_manual(values = c("magenta", "forestgreen","red"))+theme_bw()

write_csv(tripops,"tripops.csv")

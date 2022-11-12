#New pheno analysis script

rm(list=ls())
gc()

library(tidyverse)
library(ggplot2)


{
blibPheno <- read.csv("Blib_genetics.csv")
  blibPheno <- blibPheno%>%distinct()
popData <- read.csv("Pop_data.csv")
}

#popData <- popData%>%mutate(blibN_scaled = blibN/max(blibN))
blibPheno <- blibPheno%>%mutate(greengene = (greenAllele1+greenAllele2)/2)
ggplot(blibPheno,aes(x = time))+geom_count(aes(y = exon_intron_ratio), colour = "blue", alpha = 0.1)
ggplot(blibPheno%>%select(time, greenAllele1, greenAllele2, greengene)%>%gather("key", "value", -time), aes(time,value, colour = key))+geom_count(alpha = 0.5)

ggplot(blibPheno%>%select(time, moveAllele1, moveAllele2)%>%gather("key", "value", -time), aes(time,value, colour = key))+geom_count(alpha = 0.5)

ggplot(blibPheno%>%select(time, greenAllele1, greenAllele2, greengene)%>%gather("key", "value", -time), aes(time,value, colour = key))+geom_smooth()

ggplot(popData%>%select(t,blibN,blybN,blobN)%>%gather(key="species",value="n",-t),aes(x = t, y = log10(n), colour = species ))+geom_line()

test_df <- merge(popData, blibPheno, by.x = "t", by.y = "time")
 test_df <- test_df%>%mutate(exon_ratio_scaled = (exon_intron_ratio-min(exon_intron_ratio))/(max(exon_intron_ratio)-min(exon_intron_ratio)) )
 
ggplot(test_df,aes(x = t))+
  geom_point(aes(y =exon_ratio_scaled))+
  geom_line(aes(y = blibN_scaled))



blibPheno <- blibPheno%>%group_by(sampleGroup)


blibPheno <- blibPheno%>%mutate( 
  move_sd = sd(c(blibPheno$moveAllele1,blibPheno$moveAllele2)), 
  move_mean = mean(c(blibPheno$moveAllele1,blibPheno$moveAllele2)), 
  green_sd = sd(c(blibPheno$greenAllele1,blibPheno$greenAllele2)),
  green_mean = mean(c(blibPheno$greenAllele1,blibPheno$greenAllele2)),
  life_sd = sd(c(blibPheno$lifeLengthAllele1,blibPheno$lifeLengthAllele2)),
  life_mean = mean(c(blibPheno$lifeLengthAllele1,blibPheno$lifeLengthAllele2)),
  turn_sd = sd(c(blibPheno$turnTorqueAllele1,blibPheno$turnTorqueAllele2)),
  turn_mean = mean(c(blibPheno$turnTorqueAllele1,blibPheno$turnTorqueAllele2)),
  e2rep_sd = sd(c(blibPheno$e2repA,blibPheno$e2repB)),
  e2rep_mean = mean(c(blibPheno$e2repA,blibPheno$e2repB)))

sds <- blibPheno%>%select(sampleGroup,move_sd,green_sd,life_sd,turn_sd,e2rep_sd)
means <- blibPheno%>%select(sampleGroup,move_mean,green_mean,life_mean,turn_mean,e2rep_mean)
    
phenostats_g <- sds%>%gather(key = "key", value = "value", -sampleGroup)
ggplot(phenostats_g,aes(x = sampleGroup, y = log10(1+value)^2, colour = key))+geom_line()

plot(sds$sampleGroup,sds$green_sd^2, type = "l")

#New pheno analysis script

rm(list=ls())
gc()
 

library(tidyverse)
library(ggplot2)


{
blibPheno <- read.csv("Blib_genetics.csv")
  blibPheno <- blibPheno%>%distinct(name, .keep_all = TRUE)
popData <- read.csv("Pop_data.csv")
}

ggplot(popData,aes(t/60, blibN))+geom_point()
#ggplot(blibPheno,aes(time, intron1))+geom_count(alpha = 0.5)


#ggplot(blibPheno%>%select(time, intron1, intron2, intron3,intron4)%>%gather("key", "value", -time), aes(time,value, colour = key))+geom_count(alpha = 0.4)
par(mfrow=c(1,2))
plot(popData$t,popData$blibN)
plot(blibPheno$time,blibPheno$Heterozygosity, type = "p")
par(mfrow=c(1,1))
ggplot(blibPheno,aes(generation,Heterozygosity))+geom_count(alpha = 0.1)+geom_smooth(method="lm")


#popData <- popData%>%mutate(blibN_scaled = blibN/max(blibN))
blibPheno <- blibPheno%>%mutate(greengene = (greenAllele1+greenAllele2)/2)
ggplot(blibPheno,aes(x = time))+geom_count(aes(y = (exon_intron_ratio)), colour = "blue", alpha = 0.1)+geom_smooth(aes(y = exon_intron_ratio),colour = "red", alpha = 0.5)


ggplot(blibPheno,aes(x = time))+geom_count(aes(y = (Heterozygosity)), colour = "blue", alpha = 0.1)+geom_smooth(aes(y = Heterozygosity),colour = "red", alpha = 0.01)
hist(blibPheno$exon_intron_ratio)

ggplot(blibPheno%>%select(time, greenAllele1, greenAllele2)%>%gather("key", "value", -time), aes(time,value, colour = key))+geom_count(alpha = 0.1)+geom_smooth()

hist((blibPheno$greenAllele1+blibPheno$greenAllele2)/2)

ggplot(blibPheno%>%select(time, lifeLengthAllele1, lifeLengthAllele2)%>%gather("key", "value", -time), aes(time,value, colour = key))+geom_jitter(alpha = 0.2)+geom_smooth(colour = "black")

hist((blibPheno$lifeLengthAllele1+blibPheno$lifeLengthAllele2)/2)

ggplot(blibPheno%>%select(time, turnTorqueAllele1, turnTorqueAllele2)%>%gather("key", "value", -time), aes(time,value, colour = key))+geom_jitter(alpha = 0.2)+geom_smooth(colour = "black")


hist((blibPheno$turnTorqueAllele1+blibPheno$turnTorqueAllele2)/2)

ggplot(blibPheno%>%select(time, blueAllele1, blueAllele2)%>%gather("key", "value", -time), aes(time,value, colour = key))+geom_count(alpha = 0.2)

ggplot(blibPheno%>%select(time, moveAllele1, moveAllele2)%>%gather("key", "value", -time), aes(time,value, colour = key))+geom_count(alpha = 0.4)+geom_smooth(aes(colour = key))

hist((blibPheno$moveAllele1+blibPheno$moveAllele2)/2)

ggplot(blibPheno%>%select(time, greenAllele1, greenAllele2)%>%gather("key", "value", -time), aes(time,value, colour = key))+geom_smooth(aes(colour = key))

ggplot(popData%>%select(t,blibN,blybN,blobN)%>%gather(key="species",value="n",-t),aes(x = t, y = log10(n), colour = species ))+scale_color_manual(values = c("green", "magenta", "yellow"))+geom_line()+theme_dark()


exonmod <- glm(exon_intron_ratio~time, data = blibPheno)
plot(exonmod)
summary(exonmod)

blibPheno <- blibPheno%>%group_by(sampleGroup)

plot(blibPheno$exon_intron_ratio, type = "l")
abline(exonmod, col = "red")

blibPheno <- blibPheno %>% group_by(time)

blibPheno<- blibPheno%>%mutate(
  moveVariance = var(c(moveAllele1,moveAllele2)),
  redVariance = var(c(redAllele1,redAllele2)),
  greenVariance = var(c(greenAllele1,greenAllele2)),
  blueVariance = var(c(blueAllele1,blueAllele2)),
  lifeLengthVariance = var(c(lifeLengthAllele1,lifeLengthAllele2)),
  turnVariance = var(c(turnTorqueAllele1,turnTorqueAllele2)),
  e2repVariance = var(c(e2repA,e2repB))
                               )
plot(popData$t,popData$blibN, type = "l")
par(mfrow = c(2,2))
plot(blibPheno$time, sqrt(blibPheno$greenVariance), type = "l", col = "darkgreen")
plot(blibPheno$time, sqrt(blibPheno$moveVariance), type = "l")
plot(blibPheno$time, sqrt(blibPheno$turnVariance), type = "l")
plot(blibPheno$time, sqrt(blibPheno$lifeLengthVariance), type = "l")
plot(blibPheno$time, sqrt(blibPheno$e2repVariance), type = "l")


test_df <- merge(popData, blibPheno, by.x = "t", by.y = "time")


mod <- lm(greenVariance~blibN, data = test_df)
summary(mod)
plot(mod)

test_df <- test_df%>%mutate(exon_ratio_scaled = (exon_intron_ratio-min(exon_intron_ratio))/(max(exon_intron_ratio)-min(exon_intron_ratio)) )


test_df <- test_df%>%mutate(
  blibN_scaled = (blibN-min(blibN))/( max(blibN)-min(blibN) ) ,
  
  greenVariance_scaled =(greenVariance-min(greenVariance))/( max(greenVariance)-min(greenVariance) ),
  
  moveVariance_scaled =(moveVariance-min(moveVariance))/( max(moveVariance)-min(moveVariance) )
  
  )

cov(test_df$blibN,test_df$moveVariance)


ggplot(test_df%>%select(t,blibN_scaled,moveVariance_scaled,greenVariance_scaled)%>%gather(key ="key", value= "value", -t), aes(x = t, y = value, colour = key))+geom_smooth()+geom_point()


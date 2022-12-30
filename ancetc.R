#blib anc etc

{
  library(ape)
  library(phytools)
  library(tidyverse)
  library(Biostrings)
  library(phangorn)
  library(factoextra)
  library(seqinr)
  library(scales)
}

#AB_aligned <- read.alignment("AB_aligned2.fasta", format = "fasta")
#Aphy <- as.phyDat(read.alignment("Abin_aligned.fasta", format = "fasta"))
#Aphy_sample <- as.phyDat(read.alignment("A_sample_aligned.fasta", format = "fasta"))
{
  rm(list=ls())
  gc()
}
#ABphy <- as.phyDat(AB_aligned)

{
  
  blib_raw <- bind_rows(lapply(list.files(pattern = "blib_sats*"), read.csv))
  blib_raw<- blib_raw%>%distinct(name, .keep_all = TRUE)
 
  
  divergence_ratios <- blib_raw%>%select(time,divergence_ratio_A, divergence_ratio_B, generation)
  divergence_ratios$meanDivergence <- (divergence_ratios$divergence_ratio_A+divergence_ratios$divergence_ratio_B)/2
  
  blib_raw<- blib_raw%>%distinct(name, .keep_all = TRUE)%>%select(-divergence_ratio_A,-divergence_ratio_B)
  
  gc()
 
}





sample_geno <- blib_raw%>%filter(sampleGroup %% 1 == 0)%>%group_by(sampleGroup)%>%slice_sample(n = 1)



plot(divergence_ratios$divergence_ratio_A,divergence_ratios$divergence_ratio_B,pch = 20, col = scales::alpha("black", 0.1))
abline(lm(divergence_ratio_B~divergence_ratio_A, data = divergence_ratios), col = "red")



ggplot(gather(divergence_ratios%>%select(-meanDivergence,-time), key = "chSet", value = "N_SNP", -generation), aes(generation,(N_SNP), colour = chSet))+geom_point(alpha=0.3)+geom_smooth(method = "lm", colour = "black")

ggplot(divergence_ratios,aes(generation, meanDivergence))+geom_smooth(method = "lm", colour = "red", alpha = 0.4)+geom_point(alpha=0.3)+ggtitle("Number of SNPs by generation")

plot(divergence_ratios$generation,((divergence_ratios$divergence_ratio_A)), col = "red", type = "p",lwd = 0.1)
points(divergence_ratios$generation,((divergence_ratios$divergence_ratio_B)), col = "blue", lwd = 0.1)
ggplot(divergence_ratios,aes(generation,( abs(divergence_ratio_A)-abs(divergence_ratio_B))))+geom_point(alpha = 0.1)+geom_smooth(method = "loess")


                             
ggplot(divergence_ratios,aes(time,((divergence_ratio_A/(9*486))+(divergence_ratio_B/(9*486)))/2 ))+geom_count(alpha = 0.1)+geom_smooth()



fit2 <- lm(meanDivergence~poly(generation,2,raw=TRUE), data=divergence_ratios)

plot(divergence_ratios$generation,divergence_ratios$meanDivergence,lwd = 0.1)
lines(divergence_ratios$generation,fit2$fitted.values, col = "blue")


#popData <- read.csv("Pop_data.csv")
#blibPheno <- read.csv("Blib_genetics.csv")
#blibPheno<- blibPheno%>%distinct(name, .keep_all = TRUE)



#blibPheno$name <- gsub("blib","",blibPheno$name)
#blibPheno$name <- as.numeric(blibPheno$name)

#first_geno <- blib_raw%>%filter(time == min(time))
#first_geno_names <- paste0("A_",first_geno$sampleGroup,"_",first_geno$sample_number)
  
#last_geno <- blib_raw[1:nrow(blib_raw%>%filter(time == min(time))),]

last_geno <- blib_raw[1,]
last_geno <- rbind(last_geno,blib_raw%>%filter(sampleGroup >= (max(sampleGroup))))


#last_geno <- blib_raw%>%filter(sampleGroup >= (max(sampleGroup)))
#last_geno <- last_geno%>%filter(str_length(testA) == max(str_length(testA)))
#last_pheno <- blibPheno%>%filter(time >= max(last_geno$time))

last_geno$name_A <- paste0("A_",last_geno$sampleGroup,"_",last_geno$sample_number)
last_geno$name_B <- paste0("B_",last_geno$sampleGroup,"_",last_geno$sample_number)

#sample_geno <- blib_raw%>%group_by(sampleGroup)%>%slice_sample(n = 5)
sample_geno$name_A <- paste0("A_",sample_geno$sampleGroup,"_",sample_geno$sample_number)
sample_geno$name_B <- paste0("B_",sample_geno$sampleGroup,"_",sample_geno$sample_number)

{
  Astring <- DNAStringSet(last_geno$testA)
  Bstring <- DNAStringSet(last_geno$testB)
  
  #names(Astring) <- paste0("A_",last_geno$sampleGroup,"_",last_geno$sample_number)
  #names(Bstring) <- paste0("B_",last_geno$sampleGroup,"_",last_geno$sample_number)
  
  names(Astring) <- last_geno$name_A
  names(Bstring) <- last_geno$name_B
  
  Astring <- unique(Astring)
  Abin <- as.DNAbin(Astring)
  Aphy <- as.phyDat(Abin)
  
  Bstring <- unique(Bstring)
  Bbin <- as.DNAbin(Bstring)
  Bphy <- as.phyDat(Bbin) 
  
  ABstring <- c(Astring,Bstring)
  ABstring <- unique(ABstring)
  ABbin <- as.DNAbin(ABstring)
  ABphy <- as.phyDat(ABbin) 
  
}




#write.csv(blib_raw,"very_long_run.csv")

write.FASTA(Abin,"Abin.fasta")
#write.FASTA(Abin,"Bbin.fasta")
#write.FASTA(ABbin,"ABbin.fasta")

#alignment <-as.phyDat( read.alignment("alignment.fas", format = "fasta"))

{
  #Create distance matrix with JC69
  DistA_JC=dist.ml(Aphy,  model = "JC69") 
  #Create NJ tree from distance matrix
  NJ_A = NJ(DistA_JC)
  
  #Calculate likelihood of NJ tree and create start tree for optimised ML tree
  ML_START_A=pml(NJ_A, Aphy)
  #Use start tree and likelihood to create optimized ML tree
  ML_A=optim.pml(ML_START_A,model="GTR",  optNni = T) 
  ML_A_root <- root(ML_A$tree,out = as.character(Astring@ranges@NAMES[1]))
}
#DECIPHER::BrowseSeqs(ABstring,highlight = 0)

tip.date <- last_geno$time[last_geno$name_A%in%ML_A$tree$tip.label]/(60)
#+rnorm(ML_A_root$tip.label,0,2)
ML_A_root <- (rtt(ML_A$tree,tip.date))

mu <- estimate.mu(ML_A_root, tip.date)

node.date <- estimate.dates(ML_A_root, tip.date, mu, nsteps = 300)
#node.date <- estimate.dates(ML_A_root, node.date, mu, nsteps = 0, lik.tol = 1e-5)

ML_A_root$edge.length <- node.date[ML_A_root$edge[, 2]] - node.date[ML_A_root$edge[, 1]]

plot(ML_A_root, cex = 0.5);axisPhylo();
plot(ladderize(root(ML_A$tree, out = "A_0_0")), cex = 0.5);axisPhylo();
elte <- ltt(drop.tip(ML_A_root, "A_0_0"))



{
  #Create distance matrix with JC69
  DistB_JC=dist.ml(Bphy,  model = "JC69") 
  #Create NJ tree from distance matrix
  NJ_B = NJ(DistB_JC)
  
  #Calculate likelihood of NJ tree and create start tree for optimised ML tree
  ML_START_B=pml(NJ_B, Bphy)
  #Use start tree and likelihood to create optimized ML tree
  ML_B=optim.pml(ML_START_B,model="GTR",  optNni = T) 
  ML_B_root <- root(ML_B$tree,out = as.character(Bstring@ranges@NAMES[1]))
}




tip.date <- last_geno$time[last_geno$name_B%in%ML_B$tree$tip.label]
#+rnorm(ML_B_root$tip.label,0,2)
ML_B_root <- (rtt(ML_B$tree,tip.date))

mu <- estimate.mu(ML_B_root, tip.date)

node.date <- estimate.dates(ML_B_root, tip.date, mu, nsteps = 300)
#node.date <- estimate.dates(ML_B_root, node.date, mu, nsteps = 0, lik.tol = 1e-4)

ML_B_root$edge.length <- node.date[ML_B_root$edge[, 2]] - node.date[ML_B_root$edge[, 1]]

plot(ML_B_root, cex = 0.5);axisPhylo();
ltt(drop.tip(ML_B_root, "B_0_0"))

mltt.plot(drop.tip(ML_A_root, "A_0_0"),drop.tip(ML_B_root, "B_0_0"), log ="y")



{
  #Create distance matrix with JC69
  DistAB_JC=dist.ml(ABphy,  model = "JC69") 
  #Create NJ tree from distance matrix
  NJ_AB = NJ(DistAB_JC)
  
  #Calculate likelihood of NJ tree and create start tree for optimised ML tree
  ML_START_AB=pml(NJ_AB, ABphy)
  #Use start tree and likelihood to create optimized ML tree
  ML_AB=optim.pml(ML_START_AB,model="GTR",  optNni = T) 
  ML_AB_root <- root(ML_AB$tree,out = as.character(ABstring@ranges@NAMES[1]))
}
gc()
ML_AB_root <- drop.tip(ML_AB_root,c("A_0_0","B_0_0"))
namevector <- c(last_geno$name_A,last_geno$name_B)
tip.date <- last_geno$time[namevector%in%ML_AB$tree$tip.label]+rnorm(ML_AB$tree$tip.label,0,10)
ML_AB_root <- (rtt(ML_AB_root,tip.date ))

mu <- estimate.mu(ML_AB_root, tip.date)

node.date <- estimate.dates(ML_AB$tree, tip.date, nsteps = 100)
#node.date <- estimate.dates(ML_AB_root, node.date, mu, nsteps = 0, lik.tol = 1e-4)

ML_AB_root$edge.length <- node.date[ML_AB_root$edge[, 2]] - node.date[ML_AB_root$edge[, 1]]

plot(ML_AB_root, cex = 0.5);axisPhylo();
ltt(ML_AB_root)





par(mfrow=c(1,2))
plot.phylo(ML_A_root, cex = 0.3, show.tip.label = T, type = "phylogram",main = "Ch A, blibs last sample group + first sample")
plot.phylo(ML_B_root, cex = 0.3, show.tip.label = T, type = "phylogram",main = "Ch B, blibs last samplegroup + first sample", direction = "leftwards")

par(mfrow=c(1,1))

plot.phylo(drop.tip(ML_AB$tree, c("A_0_0", "B_0_0")), cex = 0.3, show.tip.label = T, type = "phylogram",main = "Ch AB, blue blibs last sample group + first sample");axisPhylo()
#ML_AB_root <- drop.tip(ML_AB_root,c("A_0_0", "B_0_0"))

ML_A_root$tip.label <- gsub("A","ind",ML_A_root$tip.label)
ML_B_root$tip.label <- gsub("B","ind",ML_B_root$tip.label)

ML_B_root <- drop.tip(ML_B_root,setdiff(ML_B_root$tip.label,ML_A_root$tip.label))

ML_A_root <- drop.tip(ML_A_root,setdiff(ML_A_root$tip.label,ML_B_root$tip.label))

par(mfrow=c(1,2))
plot.phylo(ladderize(ML_A_root), cex = 0.3, show.tip.label = F, type = "unrooted")
plot.phylo(ladderize(ML_B_root), cex = 0.3, show.tip.label = F, type = "unrooted")
par(mfrow=c(1,1))


#plot.phylo(ladderize(ML_AB$tree), cex = 0.3, show.tip.label = T, type = "phylogram")
#ML_AB$tree <- root(ML_AB$tree, outgroup = c("B_57_73","B_57_7","B_57_50","B_57_6","B_57_84","B_57_214","B_57_23"))

#ML_AB$tree <- root(ML_AB$tree, outgroup = c("B_57_73"))

#plot.phylo(ladderize(ML_AB$tree), cex = 0.3, show.tip.label = T, type = "phylogram")





timetree_A <- (chronos(ML_A_root, lambda = 10, calibration = makeChronosCalib(ML_A_root, node = "root", age.min = max(last_geno$time)/(60*60),age.max = max(last_geno$time)/(60*60), interactive = FALSE, soft.bounds = FALSE)))


timetree_B <- (chronos(ML_B_root, lambda = 10, calibration = makeChronosCalib(ML_B_root, node = "root", age.min = max(last_geno$time)/(60*60),
                                                                             age.max = max(last_geno$time)/(60*60), interactive = FALSE, soft.bounds = FALSE)))




ltt(timetree_A)
ltt(timetree_B)




indtree <- speciesTree(c(timetree_A,timetree_B))



kronoviz(c(timetree_A, timetree_B, indtree), horiz = FALSE,
         type = "phylogram", cex = 0.3, font = 1)
kronoviz(c(ladderize(timetree_A), ladderize(timetree_B), ladderize(indtree)), horiz = FALSE,
         type = "phylogram", cex = 0.3, font = 1)


par(mfrow=c(1,2))
plot.phylo(ladderize(timetree_A), cex = 0.3, show.tip.label = T, type = "phylogram")
plot.phylo(ladderize(timetree_B), cex = 0.3, show.tip.label = T, type = "phylogram", direction =  "leftwards")

par(mfrow=c(1,2))
plot.phylo(ladderize(indtree), cex = 0.3, show.tip.label = T, type = "phylogram")

plot.phylo(ladderize(timetree_A), cex = 0.3, show.tip.label = T, type = "phylogram", direction =  "leftwards")


par(mfrow=c(1,1))
plot(ladderize(indtree), cex = 0.2);axisPhylo()
mltt.plot(timetree_A, timetree_B,indtree, log = "y", xlab = "Time (h)");title( main = "LTT blib" )
par(mfrow=c(1,1))
compare.chronograms(indtree,timetree_A)
yule(indtree)

densiTree(c(timetree_A,timetree_B), type = "phylogram", cex = 0.3)
mltt.plot(timetree_A,timetree_B, log = "y");title( main = "LTT plot of time trees, blue blib")







par(mfrow=c(1,1))

plot.phylo(ape::consensus(timetree_A,timetree_B), cex = 0.3, show.tip.label = T, type = "phylogram",main = "Consensus tree, blue blibs last sample group + first sample")

{
{
  Astring_sample <- DNAStringSet(sample_geno$testA)
  Bstring_sample <- DNAStringSet(sample_geno$testB)
  
 # names(Astring_sample) <- paste0("A_",sample_geno$sampleGroup,"_",sample_geno$sample_number)
 # names(Bstring_sample) <- paste0("B_",sample_geno$sampleGroup,"_",sample_geno$sample_number)
  
  names(Astring_sample) <- sample_geno$name
  names(Bstring_sample) <- sample_geno$name
  
  Astring_sample <- unique(Astring_sample)
  Abin_sample <- as.DNAbin(Astring_sample)
  Aphy_sample <- as.phyDat(Abin_sample)
  
  Bstring_sample <- unique(Bstring_sample)
  Bbin_sample <- as.DNAbin(Bstring_sample)
  Bphy_sample <- as.phyDat(Bbin_sample) 
  
  ABstring_sample <- c(Astring_sample,Bstring_sample)
  ABstring_sample <- unique(ABstring_sample)
  ABbin_sample <- as.DNAbin(ABstring_sample)
  ABphy_sample <- as.phyDat(ABbin_sample) 
}
#write.FASTA(Abin_sample,"Abin_sample.fasta")
#write.FASTA(Bbin_sample,"Bbin_sample.fasta")
#write.FASTA(ABbin_sample,"ABbin_sample.fasta")


#write.FASTA(Abin,"Abin.fas")

#alignment <-as.phyDat( read.alignment("alignment.fas", format = "fasta"))

{
  #Create distance matrix with JC69
  DistA_JC_sample=dist.ml(Aphy_sample,  model = "JC69") 
  #Create NJ tree from distance matrix
  NJ_A_sample = NJ(DistA_JC_sample)
  
  #Calculate likelihood of NJ tree and create start tree for optimised ML tree
  ML_START_A_sample=pml(NJ_A_sample, Aphy_sample)
  #Use start tree and likelihood to create optimized ML tree
  ML_A_sample=optim.pml(ML_START_A_sample,model="GTR",  optNni = T) 
  ML_A_root_sample <- root(ML_A_sample$tree,out = as.character(Astring_sample@ranges@NAMES[1]))
}
plot(ML_A_root_sample, cex = 0.3)
  

{
  #Create distance matrix with JC69
  DistB_JC_sample=dist.ml(Bphy_sample,  model = "JC69") 
  #Create NJ tree from distance matrix
  NJ_B_sample = NJ(DistB_JC_sample)
  
  #Calculate likelihood of NJ tree and create start tree for optimised ML tree
  ML_START_B_sample=pml(NJ_B_sample, Bphy_sample)
  #Use start tree and likelihood to create optimized ML tree
  ML_B_sample=optim.pml(ML_START_B_sample,model="GTR",  optNni = T) 
  ML_B_root_sample <- root(ML_B_sample$tree,out = as.character(Bstring_sample@ranges@NAMES[1]))
}

  
  {
    #Create distance matrix with JC69
    DistAB_JC_sample=dist.ml(ABphy_sample,  model = "JC69") 
    #Create NJ tree from distance matrix
    NJ_AB_sample = NJ(DistAB_JC_sample)
    
    #Calculate likelihood of NJ tree and create start tree for optimised ML tree
    ML_START_AB_sample=pml(NJ_AB_sample, ABphy_sample)
    #Use start tree and likelihood to create optimized ML tree
    ML_AB_sample=optim.pml(ML_START_AB_sample,model="GTR",  optNni = T) 
    ML_AB_root_sample <- root(ML_AB_sample$tree,out = as.character(ABstring_sample@ranges@NAMES[1]))
  }
  

  
}




tip.date <- sample_geno$time[sample_geno$name%in%ML_A_sample$tree$tip.label]
ML_A_sample$tree <- (rtt(ML_A_sample$tree,tip.date))

mu <- estimate.mu(ML_A_sample$tree, tip.date)

node.date <- estimate.dates(ML_A_sample$tree, tip.date, mu, nsteps = 100)
#node.date <- estimate.dates(ML_A_sample$tree, node.date, mu, nsteps = 0, lik.tol = 1e-4)

ML_A_sample$tree$edge.length <- node.date[ML_A_sample$tree$edge[, 2]] - node.date[ML_A_sample$tree$edge[, 1]]

ltt.plot(ML_A_sample$tree)


plot.phylo(ML_A_sample$tree, cex = 0.2, show.tip.label = T, type = "phylogram",main = "Ch A, blibs sampled over time,  timetree");axisPhylo()

plot.phylo(NJ_A_sample, cex = 0.2, show.tip.label = T, type = "phylogram",main = "Ch A, blibs sampled over time,  timetree");axisPhylo()

nodelabels(frame = "none", cex = 0.3)

DECIPHER::BrowseSeqs(Astring_sample, highlight = 0)

par(mfrow=c(1,2))
plot.phylo(ladderize(ML_A_root_sample), cex = 0.3, show.tip.label = T, type = "phylogram",main = "Ch A, blue blibs sampled over time")
plot.phylo(ladderize(ML_B_root_sample), cex = 0.3, show.tip.label = T, type = "phylogram",main = "Ch B, blue blibs sampled over time", direction = "leftwards")

par(mfrow=c(1,1))
plot.phylo(ladderize(ML_AB_root_sample), cex = 0.3, show.tip.label = T, type = "phylogram",main = "Ch AB, blue blibs sampled over time")


library(ade4)
ade4::table.paint(as.data.frame(as.matrix(dist.dna(Abin_sample))))




plot(collapse.singles(ladderize(ML_A$tree)), cex = 0.3)
#Calculate bootstrap support with 200 replicates for ML tree 1
Boot1=bootstrap.pml(ML_A, bs=100, optNni=TRUE, model="GTR")


#Add bootstrap values
BSTree1 <-  plotBS(ML_A$tree, Boot1, p=0 ,type="phylogram", main="blibbe", cex = 0.3)

sampleD.pca <- prcomp(DistB_JC, scale = TRUE)
fviz_eig(sampleD.pca)


fviz_pca_ind(sampleD.pca,
             geom = c("point"),
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             
             
)




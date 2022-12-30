#Stringcomparer



string.diff.ex<-function(a,b,exclude=c("n","N","?"),ignore.case=TRUE)
{
  if(nchar(a)!=nchar(b)) stop("Lengths of input strings differ. Please check your input.")
  if(ignore.case==TRUE)
  {
    a<-toupper(a)
    b<-toupper(b)
  }
  diff.a<-unlist(strsplit(a,split=""))
  diff.b<-unlist(strsplit(b,split=""))
  diff.d<-rbind(diff.a,diff.b)
  for(ex.loop in 1:length(exclude))
  {
    diff.d<-diff.d[,!(diff.d[1,]==exclude[ex.loop]|diff.d[2,]==exclude[ex.loop])]
  }
  differences<-sum(diff.d[1,]!=diff.d[2,])
  return(differences)
}

thediffs <- c(rep(0,nrow(blib_raw)))
for(i in 1:nrow(blib_raw)){
  thediffs[i] <- string.diff.ex(blib_raw$testA[i],blib_raw$testB[i])
}
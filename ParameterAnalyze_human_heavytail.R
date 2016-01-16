

# 
# E.s.plus.vec <- seq(0.07,0.041,-0.001)
# mu.vec <- seq(0.000025,0.0007,0.000025)

inc.array <- array(NA, c(length(E.s.plus.vec),length(mu.vec),100))

for(qq in 1:length(E.s.plus.vec)){
  for(ww in 1:length(mu.vec)){
    mypath <- file.path(getwd(),"data",paste("heavyS_",E.s.plus.vec[qq],"Mu_",mu.vec[ww],".csv", sep = ""),sep="")
    test <- read.csv(mypath)
    inc.array[qq,ww,] <- test$x
  }
}


# par(mar=c(5,6.5,4,5)+.1)
# plot(inc.array[27,18,],ylim=c(0,100),type="l",lwd=3,lty=4, xlab="Age (years)", ylab="Percent of the population with 
#      an initiated tumor",main="Population tumor incidence in humans; heavy-tailed scenario",cex.lab=1.5,cex.axis=1.3)
# lines(inc.array[27,19,],lwd=3,lty=4)
# lines(inc.array[27,20,],lwd=6,lty=4,col="red")
# lines(inc.array[27,21,],lwd=3,lty=4)
# lines(inc.array[27,22,],lwd=3,lty=4)
# # lines(inc.array[4,20,],lwd=3,lty=4)
# # lines(inc.array[4,21,],lwd=3,lty=4)
# polyp2 <- polyp[,2]
# points(polyp2,x=polyp[,1],type="o",lwd=4,col="red")
# abline(v=0);abline(h=0)


# E.s.plus.vec[27]
# mu.vec[20]


# par(mar=c(5,6.5,4,5)+.1)
# plot(inc.array[28,24,],ylim=c(0,100),type="l",lwd=3,lty=4, xlab="Age (years)", ylab="Percent of the population with 
#      an initiated tumor",main="Population tumor incidence in humans; heavy-tailed scenario",cex.lab=1.5,cex.axis=1.3)
# lines(inc.array[28,25,],lwd=3,lty=4)
# lines(inc.array[27,20,],lwd=3,lty=4,col="red")
# lines(inc.array[26,17,],lwd=3,lty=4)
# # lines(inc.array[27,22,],lwd=3,lty=4)
# # lines(inc.array[4,20,],lwd=3,lty=4)
# # lines(inc.array[4,21,],lwd=3,lty=4)
# polyp2 <- polyp[,2]
# points(polyp2,x=polyp[,1],type="o",lwd=4,col="red")


# mypath <- file.path("C:","Users","Vincent","Desktop","MS_figwork","SandMu",paste("inc.array",".RData", sep = ""))
# 
# save(inc.array,file=mypath)

# polyp <- read.csv("C:/Users/Vincent/Dropbox/Modeling tumor evolution/information/data/polyp.csv")
# polyp2 <- polyp[,2]
# points(polyp2,x=polyp[,1],type="o",lwd=4,col="red")

ls.fun <- function(incidence){
  #   incidence2 <- indidence 
  polyp.vec <- polyp$age[1:9]
  to.sum <- rep(NA, length(polyp.vec))
  for(i in 1:length(polyp.vec)){
    to.sum[i] <- abs(polyp2[i]-incidence[polyp.vec[i]])^2
  }
  least.squared <- (1/length(polyp.vec))*sum(to.sum)
  return(least.squared)
}
ls.fun(inc.array[1,4,])

ls.mat <- matrix(nrow=length(inc.array[,1,1]),ncol=length(inc.array[1,,1]))

for(jj in 1:length(ls.mat[,1])){
  for(kk in 1:length(ls.mat[1,])){
    ls.mat[jj,kk] <- ls.fun(inc.array[jj,kk,])
  }
}

# which(ls.mat==min(ls.mat))

#mypath <- file.path("C:","Users","Vincent","Desktop","MS_figwork","SandMu2",paste("ls.mat2",".RData", sep = ""))
#save(ls.mat,file=mypath)
# heatmap(ls.mat, Rowv=NA, Colv=NA, col=heat.colors(256))

ls.mat3 <- round(ls.mat,2)
# 
# require(gplots)
# par(mar=c(5,6.5,4,5)+.1)
# heatmap.2(ls.mat,Rowv=NA,Colv=NA,dendrogram='none',trace='none',cellnote=ls.mat3,notecex=0.75,labRow=E.s.plus.vec,labCol=mu.vec,xlab="Mutation Rate",ylab="E[s+]",margins=c(7,6),breaks=c(0,66,100,150,200,250,300,350,500),col=redblue,cexRow=1,key=T,density.info=c("none"))

heatmap.2(ls.mat3,Rowv=NA,
          Colv=NA, 
          dendrogram='none',
          trace='none',
          cellnote=ls.mat3,
          notecex=0.75,
          labRow=E.s.plus.vec,
          labCol=mu.vec,
          xlab="Mutation Rate",
          ylab="E[s+]",
          margins=c(9,6.5),
          breaks=c(0,66,67,75,150,200,250,300,350,500),
          col=redblue,density.info=c("none"),
          key=F,
          keysize=0.4,
          cexRow=1.3,
          cexCol=1.3)



# plot(inc.array[1,4,],ylab="Percent of the population with a tumor", xlab="Age",lty=1,type="l",lwd=5,ylim=c(0,100))
par(mar=c(5,6.5,4,5)+.1)
plot(inc.array[27,18,],ylim=c(0,100),type="l",lwd=3,lty=4, xlab="Age (years)", ylab="Percent of the population with 
     an initiated tumor",main="Population tumor incidence in humans; power-law scenario",cex.lab=1.5,cex.axis=1.3)
lines(inc.array[27,19,],lwd=3,lty=4)
lines(inc.array[27,20,],lwd=6,lty=4,col="red")
lines(inc.array[27,21,],lwd=3,lty=4)
lines(inc.array[27,22,],lwd=3,lty=4)
# lines(inc.array[4,20,],lwd=3,lty=4)
# lines(inc.array[4,21,],lwd=3,lty=4)
polyp2 <- polyp[,2]
points(polyp2,x=polyp[,1],type="o",lwd=4,col="red")
abline(v=0);abline(h=0)

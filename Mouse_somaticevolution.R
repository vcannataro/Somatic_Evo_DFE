###Author: Vincent Cannataro 



# require("caTools")
div.exp <- 1
div.pow <- 0

# lambda.0 <- 0.20  ##initial division rate of stem cells 
# N <- 6 ##number of stem cells in the stem cell niche
# N.tot <- 15
# nu.0 <- 0.333
diff.rate <- nu.0  
# crypts <- 12e5

iter <- 8 #number of mutations to calculate (and lambda densities for each mutation). This is the counter you wait for when running the code
resolut <- 0.001    #how resolute we want our mutational state space to be
x.max.distance <- 1  #mutation state space will go from 0 to this number 
mut.space <- seq(resolut, x.max.distance, by=resolut) #mutational space to calculate the density for
mut.space <- round(mut.space,nchar(resolut))

#polyp <- read.csv("C:/Users/Vincent/Dropbox/Modeling tumor evolution/information/data/polyp.csv")

###setting alpha and beta and the corresponding DFE/Pfix values
if(div.exp==1){
  alpha <- 1/((s.b+1)-1)
  beta <- 1/(1-(1-s.d))
  DFE <- function(alpha, beta, lambda, lambda.0){
    if(lambda < lambda.0){
      return((1-P.B)*(beta/lambda.0)*exp(-beta*(1-(lambda/lambda.0))))
    }
    if(lambda > lambda.0){
      return(P.B*(alpha/lambda.0)*exp(-alpha*((lambda/lambda.0)-1)))
    }
    if(lambda == lambda.0){
      return(0)
    }
  }
  
}
Pfix <- function(alpha,beta,N){
  integrand.1 <- function(u){
    ((1-(1/u))/(1-(1/u)^N))*exp(-beta*(1-u))*(1-P.B)*beta 
  }
  integrand.2 <- function(u){
    ((1-(1/u))/(1-(1/u)^N))*exp(-alpha*(u-1))*P.B*alpha
  }
  output <- integrate(integrand.1,1e-5,1)$value + integrate(integrand.2,1,Inf)$value
  return(output)
}



##A function that creates a vector containing the discretized values of the DFE
DFE.vec <- function(lambda.0, mut, alpha, beta){
  x.space <- mut
  den <- rep(0,length(x.space))
  for(i in 1:length(den)){
    den[i] <- DFE(alpha,beta,x.space[i],lambda.0)
  }
  return(den)
}


##A function that creates discretized probability of fixation for a given birth rate, lambda
prob.of.fix <- function(lambda.0,xi,N){
  out <- rep(NA,length(xi))
  for(i in 1:length(xi)){
    if(xi[i] > 0){
      out[i] <- 
        ((1 - lambda.0/(xi[i]))/(1 - (lambda.0/(xi[i]))^N))
    }else{
      out[i] <- 0
    }
  }
  out[which(is.nan(out))] <- 1/N
  return(out)
}




##Function that creates the posterior density of lambda given the DFE and probability of fixation and lambda.0

lambda.density <- function(lambda.0, xi, N, alpha, beta){
  output <- rep(NA, length(xi))
  
  prob.fix <- prob.of.fix(lambda.0=lambda.0, xi=xi, N)
  mut.dens <- DFE.vec(lambda.0,mut=xi,alpha,beta)
  
  to.sum <- rep(NA, length(xi)) 
  
  for(i in 1:length(xi)){
    to.sum[i] <- prob.fix[i]*mut.dens[i]
  }
  
  denominator<-trapz(x=xi,y=to.sum)
  
  
  
  post <- rep(NA, length(prob.fix))
  
  
  for(i in 1:length(xi)){
    post[i] <- (to.sum[i])/denominator
  }
  
  return(post)
}


##Makes a ~nice~ figure showing first lambda density against probability.fixation and DFE
plotter <- function(lambda.0, n.span, p.span, N, alpha, beta){  
  par(new=F)
  par(mar=c(5,4,4,5)+.1)
  plot(x=mut.space, y=prob.of.fix(lambda.0,mut.space, N), type="l",ylim=c(0,1),lwd=3
       ,xlab="Division Rate", main="Probability of fixation, DFE, and first lambda density", xlim=c(n.span,p.span),ylab="Probability of fixation")
  legend("topleft",paste("s.b.=",s.b, "s.d.=",s.d,"\n",
                         "Percent Ben.=", P.B*100,"\n",
                         "N=",N))
  if(div.exp==1){
    legend("topright","Exponential DFE")
  }
  if(div.pow==1){
    legend("topright","Power DFE")
  }
  abline(h=1/N)
  par(new=T)
  plot(DFE.vec(lambda.0,mut.space,alpha,beta), x=mut.space, type="l",lwd=2,xaxt="n",yaxt="n",xlab="",ylab="", col="red",xlim=c(n.span,p.span));abline(v=lambda.0, col="red"); #abline(v=0, col="blue")
  axis(4)
  mtext("Densities",side=4,line=3)
  par(new=T)
  plot(x=mut.space, y=lambda.density(lambda.0,mut.space,N,alpha,beta)
       ,lwd=4,xaxt="n",yaxt="n",xlab="",ylab="", col="green",xlim=c(n.span,p.span),lty=6,type="l")
}

##First lambda density of lambda.0
f.1 <- lambda.density(lambda.0,mut.space,N,alpha,beta)

f.1.expdiv <- f.1
##Constructs a list of the density of lambda for every possible lambda.0
##This will be useful to recall in the next section (instead of calculating it every time). 
den.list <- list()
for(j in 1:length(mut.space)){
  den.list[[j]] <- lambda.density(lambda.0=mut.space[j],mut.space,N,alpha,beta)
}

##list of the subsequent densities of lambda given another mutation
post.den <- list()
##Used to record the full density (the other list gets density > lambda cut off and redistributed
post.den.full <- list()
post.den[[1]] <- f.1
post.den.full[[1]] <- f.1

#vector recording probability of tumorigenesis for each mutation (area of lambda density
#greater than tumorigenesis threshold (nu.0)
tumor.init <- rep(0,iter)

##storing and getting rid of the first bit of tumorigenesis probability
tumor.init[1] <- trapz(x=mut.space[(which(mut.space==diff.rate)):length(mut.space)],y=post.den[[1]][(which(mut.space==diff.rate)):length(mut.space)])
##dummy density to manipulate (not needed anymore but still going to use for now because it's working)
test <- post.den[[1]]
##replacing the density after the tumorigenesis threshold with 0
test[(which(mut.space==diff.rate)+1):length(mut.space)] <- 0
##storing the new area under this density
l2.tail <- trapz(test,x=mut.space)
##renormalizing the density so the area = 1
test <- test*(1/l2.tail) 
##restoring the density to it's proper place
post.den[[1]] <- test 

###loop that creates the probability density of lambda for subsequent mutations. 
for(z in 2:iter){
  
  post.den[[z]] <- rep(0,length=length(mut.space))
  
  ##for every possible value of lambda (let's call it lambda.A)
  for(i in 1:length(mut.space)){
    #this value of lambda
    this.y <- mut.space[i]
    # a vector of 0s the length of the possible lambda values
    this.den <- rep(0,length=length(mut.space))
    ##for every other possible value of lambda (let's call it lambda.B)
    for(j in 1:length(mut.space)){
      #The value of lambda.A given you started at every lambda.B
      #for your lambda.0, multiplied by the value of the previous
      #lambda density for the lambda.B value
      this.den[j] <- den.list[[j]][i]*post.den[[z-1]][j]
      #print(j)
    }
    #Intregral of all this over possible lambda.B to get the value of lambda.A,
    #which will be the density of the new lambda for z mutations
    post.den[[z]][i] <- trapz(x=mut.space,y=this.den)
    #     print(i)
    
  }
  #same method as before for storing/replacing the area after the 
  #tumorigenesis threshold.
  test <- post.den[[z]]
  post.den.full[[z]] <- test
  tumor.init[z] <- trapz(x=mut.space[(which(mut.space==diff.rate)):length(mut.space)],y=post.den[[z]][(which(mut.space==diff.rate)):length(mut.space)])
  
  test[(which(mut.space==diff.rate)+1):length(mut.space)] <- 0
  
  l2.tail <- trapz(test,x=mut.space)
  test <- test*(1/l2.tail)
  post.den[[z]] <- test 
  print(paste("Calculation of fixed mutation number",z,"out of",iter,"Exponential Scenario"))
}



#accumulated probability of tumorigenesis in the crypt given a certain number of fixed mutations
prob.cancer <- rep(0, length=length(tumor.init))
prob.cancer[1] <- tumor.init[1]

for(i in 2:length(prob.cancer)){
  prob.cancer[i] <- tumor.init[i]*(1-prob.cancer[i-1])+prob.cancer[i-1]
}

#plot(prob.cancer)

Phi <- 1-prob.cancer


p.fix <- Pfix(alpha,beta,N) ## probability of fixation for a given DFE 
#mu <- 1e-6   ##mutation rate per division
p.fix.expdiv <- p.fix

division.rate <- lambda.0 ##division rate per stem cell per day

##days in X years (since division rate is in per days)
days <- function(weeks){
  return(weeks*7)
}


##expected fixed mutations in X years
muts <- function(weeks){
  return(p.fix*mu*N*division.rate*days(weeks))
}

##The number of mutations in the crypt to explore
n.vec <- 1:iter

##The number of years in a lifetime to look at
week.vec <- 1:(3*52)

##
incidence <- rep(0,length(week.vec))

no.cancer.vec <- rep(0,length(week.vec))

for(i in 1:length(week.vec)){
  no.cancer <- Phi[1]^(crypts*dpois(1,muts(i)))
  for(j in 2:length(n.vec)){
    no.cancer <- no.cancer*(Phi[j]^(crypts*dpois(j,muts(i))))
  }
  no.cancer.vec[i] <- no.cancer
}



tumor.init.expdiv <- tumor.init
prob.cancer.expdiv <- prob.cancer
no.cancer.vec.expdiv <- no.cancer.vec

post.den.expdiv <- post.den 
post.den.full.expdiv <- post.den.full


######################################################
######################################################
######################################################

###
##Implimenting the new DFE that incorporate percent beneficial mutation and expected values
#####HEAVY-TAILED SCENARIO FOR CHANGE IN DIVISION RATE 


div.exp <- 0
div.pow <- 1



if(div.pow==1){
  alpha <- (-1+2*(1+s.b))/((1+s.b)-1)
  beta <- 1/(1-(1-s.d))
  DFE <- function(alpha, beta, lambda, lambda.0){
    if(lambda < lambda.0){
      return((1-P.B)*(beta/lambda.0)*exp(-beta*(1-(lambda/lambda.0))))
    }
    if(lambda > lambda.0){
      return(P.B*((alpha-1)/lambda.0)*(lambda/lambda.0)^(-alpha))
    }
    if(lambda == lambda.0){
      return(0)
    }
  }
  
}
Pfix <- function(alpha,beta,N){
  integrand.1 <- function(u){
    ((1-(1/u))/(1-(1/u)^N))*exp(-beta*(1-u))*(1-P.B)*beta 
  }
  integrand.2 <- function(u){
    ((1-(1/u))/(1-(1/u)^N))*u^(-alpha)*P.B*(alpha-1)
  }
  output <- (integrate(integrand.1,1e-5,1)$value + integrate(integrand.2,1,Inf)$value)
  return(output)
}



DFE.vec <- function(lambda.0, mut, alpha, beta){
  x.space <- mut
  den <- rep(0,length(x.space))
  for(i in 1:length(den)){
    den[i] <- DFE(alpha,beta,x.space[i],lambda.0)
  }
  return(den)
}

prob.of.fix <- function(lambda.0,xi,N){
  out <- rep(NA,length(xi))
  for(i in 1:length(xi)){
    if(xi[i] > 0){
      out[i] <- 
        ((1 - lambda.0/(xi[i]))/(1 - (lambda.0/(xi[i]))^N))
    }else{
      out[i] <- 0
    }
  }
  out[which(is.nan(out))] <- 1/N
  return(out)
}

lambda.density <- function(lambda.0, xi, N, alpha, beta){
  output <- rep(NA, length(xi))
  
  prob.fix <- prob.of.fix(lambda.0=lambda.0, xi=xi, N)
  mut.dens <- DFE.vec(lambda.0,mut=xi,alpha,beta)
  
  to.sum <- rep(NA, length(xi)) 
  
  for(i in 1:length(xi)){
    to.sum[i] <- prob.fix[i]*mut.dens[i]
  }
  
  denominator<-trapz(x=xi,y=to.sum)
  
  
  
  post <- rep(NA, length(prob.fix))
  
  
  for(i in 1:length(xi)){
    post[i] <- (to.sum[i])/denominator
  }
  
  return(post)
}




f.1 <- lambda.density(lambda.0,mut.space,N,alpha,beta)


f.1.powdiv <- f.1


den.list <- list()
for(j in 1:length(mut.space)){
  den.list[[j]] <- lambda.density(lambda.0=mut.space[j],mut.space,N,alpha,beta)
}

post.den <- list()
post.den.full <- list()
post.den[[1]] <- f.1
post.den.full[[1]] <- f.1

tumor.init <- rep(0,iter)
##storing and getting rid of the first bit of tumorigenesis probability
tumor.init[1] <- trapz(x=mut.space[(which(mut.space==diff.rate)):length(mut.space)],y=post.den[[1]][(which(mut.space==diff.rate)):length(mut.space)])
test <- post.den[[1]]
test[(which(mut.space==diff.rate)+1):length(mut.space)] <- 0
l2.tail <- trapz(test,x=mut.space)
test <- test*(1/l2.tail)
post.den[[1]] <- test 


for(z in 2:iter){
  
  post.den[[z]] <- rep(0,length=length(mut.space))
  
  
  for(i in 1:length(mut.space)){
    this.y <- mut.space[i]
    this.den <- rep(0,length=length(mut.space))
    for(j in 1:length(mut.space)){
      this.den[j] <- den.list[[j]][i]*post.den[[z-1]][j]
      #print(j)
    }
    post.den[[z]][i] <- trapz(x=mut.space,y=this.den)
    #     print(i)
    
  }
  test <- post.den[[z]]
  post.den.full[[z]] <- test
  tumor.init[z] <- trapz(x=mut.space[(which(mut.space==diff.rate)):length(mut.space)],y=post.den[[z]][(which(mut.space==diff.rate)):length(mut.space)])
  
  test[(which(mut.space==diff.rate)+1):length(mut.space)] <- 0
  
  l2.tail <- trapz(test,x=mut.space)
  test <- test*(1/l2.tail)
  post.den[[z]] <- test 
  print(paste("Calculation of fixed mutation number",z,"out of",iter,"Heavy-tailed Scenario"))
}



prob.cancer <- rep(0, length=length(tumor.init))
prob.cancer[1] <- tumor.init[1]

for(i in 2:length(prob.cancer)){
  prob.cancer[i] <- tumor.init[i]*(1-prob.cancer[i-1])+prob.cancer[i-1]
}

#plot(prob.cancer)

Phi <- 1-prob.cancer


p.fix <- Pfix(alpha,beta,N) ## probability of fixation for a given DFE 
#mu <- 1e-6   ##mutation rate per division
p.fix.powerdiv <- p.fix

division.rate <- lambda.0 ##division rate per stem cell per day

##days in X years (since division rate is in per days)
days <- function(weeks){
  return(weeks*7)
}


##expected fixed mutations in X years
muts <- function(weeks){
  return(p.fix*mu*N*division.rate*days(weeks))
}

##The number of mutations in the crypt to explore
n.vec <- 1:iter

##The number of years in a lifetime to look at
week.vec <- 1:(3*52)

##
incidence <- rep(0,length(week.vec))

no.cancer.vec <- rep(0,length(week.vec))

for(i in 1:length(week.vec)){
  no.cancer <- Phi[1]^(crypts*dpois(1,muts(i)))
  for(j in 2:length(n.vec)){
    no.cancer <- no.cancer*(Phi[j]^(crypts*dpois(j,muts(i))))
  }
  no.cancer.vec[i] <- no.cancer
}

tumor.init.powerdiv <- tumor.init
prob.cancer.powerdiv <- prob.cancer
no.cancer.vec.powerdiv <- no.cancer.vec

post.den.powdiv <- post.den 
post.den.full.powdiv <- post.den.full


density.plotter <- function(expdiv,powdiv,xlim1,xlim2,ylim1,ylim2,extent){
  if(expdiv==1){
    post1.new <- c(post.den.full.expdiv[[1]][1:which(mut.space==lambda.0)-1],NaN,post.den.full.expdiv[[1]][(which(mut.space==lambda.0)+1):length(post.den.full.expdiv[[1]])])
    par(mfrow =c(1,1))
    plot(post1.new, type="l",x=mut.space, lwd=6,xlim=c(xlim1,xlim2),ylim=c(ylim1,ylim2),
         ylab="Density", xlab="Division Rate (days)",main="Probability densities of division rate given fixed mutations in mice",col="dark green",lty=2, cex.lab=1.5, cex.axis=1.5,cex.main=1)
    abline(v=lambda.0,col="black",lwd=1)
    if(extent != 1){
      for(i in 2:extent){
        lines(post.den.full.expdiv[[i]], x=mut.space, type="l", lwd=4)
      }
    }
  }
  if(powdiv==1){
    post1.new <- c(post.den.full.powdiv[[1]][1:which(mut.space==lambda.0)-1],NaN,post.den.full.powdiv[[1]][(which(mut.space==lambda.0)+1):length(post.den.full.powdiv[[1]])])
    par(mfrow =c(1,1))
    plot(post1.new, type="l",x=mut.space, lwd=6,xlim=c(xlim1,xlim2),ylim=c(ylim1,ylim2),
         ylab="Density", xlab="Division Rate (days)",main="Probability densities of division rate given fixed mutations in mice",col="dark green",lty=2, cex.lab=1.5, cex.axis=1.5,cex.main=1)
    abline(v=lambda.0,col="black",lwd=1)
    if(extent != 1){
      for(i in 2:extent){
        lines(post.den.full.powdiv[[i]], x=mut.space, type="l", lwd=4)
      }
    }
  }
}

# density.plotter(1,0,xlim1=0,xlim2=.4,ylim1=0,ylim2=35,8)
# abline(v=nu.0,lwd=3,col="red")
# abline(v=0);abline(h=0)
# 
# density.plotter(1,0,xlim1=0.32,xlim2=.34,ylim1=0,ylim2=1.5e-3,8)
# abline(v=nu.0,lwd=3,col="red")
# abline(v=0);abline(h=0)
# 

# plot(tumor.init.powerdiv,lwd=5,type="o",cex.lab=1.3,cex.axis=1.3,
#      main="Probability of tumorigenesis per fixed mutation", xlab="Fixed mutation number",
#      ylab="Probability")
# lines(tumor.init.expdiv,lwd=5,type="o",cex.lab=1.3,cex.axis=1.3, lty=2,
#       main="Probability of tumorigenesis per fixed mutation", xlab="Fixed mutation number",
#       ylab="Probability")
# abline(h=0)

# par(mar=c(4,7,4,5)+.1)
# plot((1-no.cancer.vec.powerdiv)*1e2,lty=3,lwd=5,type="l", xlab="Age (weeks)", ylab="Percent of the population with 
#      an initiated tumor",main="Population tumor incidence in mice",ylim=c(0,1e2),cex.lab=1.5,cex.axis=1.3)
# lines((1-no.cancer.vec.expdiv)*1e2,lty=2,col="black", lwd=5)
# legend("topright", c("Pareto","Exponential"), lty=c(3,2),lwd=3 )
# abline(v=0);abline(h=0)


####
#Plotting the probability of tumorigenesis for fixed mutations
#### 
# par(mar=c(4,7,4,5)+.1)
# plot(prob.cancer.powerdiv, xlab="Fixed Mutation Number",ylab="Probability of Tumorigenesis in a Crypt",type="b", lwd=4,cex.lab=1.5, cex.axis=1.5,ylim=c(0,0.0002),main="Probability of Tumorigenesis in a Mouse Crypt",lty=3)
# lines(prob.cancer.expdiv, xlab="Fixed Mutation Number",ylab="Probability of Tumorigenesis in a crypt",type="b", lwd=4,cex.lab=1.5, cex.axis=1.5,lty=2)
# legend("bottomright", c("Pareto","Exponential            "), lty=c(3,2),lwd=3 )

# plot(prob.cancer.powerdiv, xlab="Fixed Mutation Number",ylab="Probability of Tumorigenesis in a crypt",type="b", lwd=3,cex.lab=1.5, cex.axis=1.5)
# 


###
#Expected value
###
expected.finder <- function(f.x,space){
  return(trapz(y=(f.x*space),x=space))
}


expected.vec.mouseexp <- rep(NA,iter)
for(i in 1:iter){
  expected.vec.mouseexp[i] <- expected.finder(post.den.expdiv[[i]],mut.space)
}



###
#mutation profiles
###

total.probcancer <- rep(NA,iter)
mut.dist.mat <- matrix(nrow=iter,ncol=length(week.vec))
for(j in 1:length(week.vec)){
  for(i in 1:iter){
    total.probcancer[i] <- (dpois(i,muts(j))*prob.cancer[i])
  }
  total.probcancersum <- sum(total.probcancer)
  for(i in 1:iter){
    mut.dist.mat[i,j] <- (dpois(i,muts(j))*prob.cancer[i])/total.probcancersum
  }
}

# par(mar=c(5,7,4,5)+.1)
# barplot(mut.dist.mat,horiz=F,main="Mutation profile at the onset of tumorigenesis:
#         Mouse, exponential DFE on division rate model",axes=T,axisnames=T,names.arg=week.vec,xlab="Age (Weeks)",col=rainbow(iter),ylab="Probability each mutation caused the tumor",cex.axis=1.5,cex.lab=1.5)
# legend("bottomright",legend=1:3,fill=rainbow(iter))



###Author: Vincent Cannataro
##



# E.s.plus.vec <- seq(0.054,0.054,0.001)
# mu.vec <- seq(0.000525,0.0007,0.000025)

#inc.array <- array(NA, c(length(E.s.plus.vec),length(mu.vec),100))

dir.create("data")

for(qq in 1:length(E.s.plus.vec)){
  for(ww in 1:length(mu.vec)){
    
    
    require("caTools")
    div.exp <- 0
    div.pow <- 1

    
    
    #     P.B <- 0.0575 ##probability of a beneficial mutation (vs. deleterious) 
    #     
    #     lambda.0 <- round(1/7,3)  ##initial division rate of stem cells /
    #     N <- 20 ##number of stem cells in the stem cell niche
    #     N.tot <- 36
    #     nu.0 <- 0.321
    #     
    #N.tot <- round((N*(-nu.0))/(lambda.0-nu.0),1)
    
    #nu.0 <- round((lambda.0*(N + (N.tot - N)))/(N.tot-N),3)
    
    
    s.b <- E.s.plus.vec[qq] #expected value given a beneficial mutation occured
    #     s.d <- 0.217 #expected value given a deleterious mutation occured
    
    diff.rate <- nu.0  
    #     crypts <- 2e7
    
    mu <- mu.vec[ww]   ##mutation rate per division
    
    iter <- 10 #number of mutations to calculate (and lambda densities for each mutation). This is the counter you wait for when running the code
    
    
    resolut <- 0.001    #how resolute we want our mutational state space to be
    x.max.distance <- 1  #mutation state space will go from 0 to this number 
    
    
    
    mut.space <- seq(resolut, x.max.distance, by=resolut) #mutational space to calculate the density for
    mut.space <- round(mut.space,nchar(resolut))
    
    #     polyp <- read.csv("C:/Users/Vincent/Dropbox/Modeling tumor evolution/information/data/polyp.csv")
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
    }
    if(div.pow==1){
      alpha <- (-1+2*(s.b+1))/((s.b+1)-1)
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
      
      Pfix <- function(alpha,beta,N){
        integrand.1 <- function(u){
          ((1-(1/u))/(1-(1/u)^N))*exp(-beta*(1-u))*(1-P.B)*beta 
        }
        integrand.2 <- function(u){
          ((1-(1/u))/(1-(1/u)^N))*u^(-alpha)*P.B*(alpha-1)
        }
        output <- integrate(integrand.1,1e-5,1)$value + integrate(integrand.2,1,Inf)$value
        return(output)
      }
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
      #print(z)
    }
    
    
    #accumulated probability of tumorigenesis in the crypt given a certain number of fixed mutations
    prob.cancer <- rep(0, length=length(tumor.init))
    prob.cancer[1] <- tumor.init[1]
    
    for(i in 2:length(prob.cancer)){
      prob.cancer[i] <- tumor.init[i]*(1-prob.cancer[i-1])+prob.cancer[i-1]
    }
    
    Phi <- 1-prob.cancer
    
    
    p.fix <- Pfix(alpha,beta,N) ## probability of fixation for a given DFE 
    
    p.fix.expdiv <- p.fix
    
    division.rate <- lambda.0 ##division rate per stem cell per day
    
    ##days in X years (since division rate is in per days)
    days <- function(years){
      return(years*365.25)
    }
    
    
    ##expected fixed mutations in X years
    muts <- function(years){
      return(p.fix*mu*N*division.rate*days(years))
    }
    
    ##The number of mutations in the crypt to explore
    n.vec <- 1:iter
    
    ##The number of years in a lifetime to look at
    year.vec <- 1:100
    
    ##
    incidence <- rep(0,length(year.vec))
    
    no.cancer.vec <- rep(0,length(year.vec))
    
    for(i in 1:length(year.vec)){
      no.cancer <- Phi[1]^(crypts*dpois(1,muts(i)))
      for(j in 2:length(n.vec)){
        no.cancer <- no.cancer*(Phi[j]^(crypts*dpois(j,muts(i))))
      }
      no.cancer.vec[i] <- no.cancer
    }
    
    no.cancer.vec.expdiv <- no.cancer.vec
    
    mypath <- file.path(getwd(),"data",paste("heavyS_",s.b,"Mu_",mu,".csv", sep = ""),sep="")
    write.csv((1-no.cancer.vec.expdiv)*1e2,file=mypath)
    #######################################
    #################################
    print(paste(ww,qq,"out of",length(mu.vec), length(E.s.plus.vec))) 
    
  }
}

#mypath <- file.path(getwd(),"data",sep="")
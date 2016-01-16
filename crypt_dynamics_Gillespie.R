###Gillespie alg, crypt dynamics, full sim 

# mouse <- 0
# human <- 1

if(mouse == 1){
  lambda.0 <- 0.20  ##initial division rate of stem cells 
  N <- 6 ##number of stem cells in the stem cell niche
  N.tot <- 15 ##total number of stem cells in the crypt
  nu.0 <- round((lambda.0*(N + (N.tot - N)))/(N.tot-N),3)  ##determines the differentiation rate necessary to maintain N and N.tot with that lambda.0
  lifetime <- 3 #lifetime in years
  col.max <- 5e4
  
  niche.matrix <- matrix(nrow=N,ncol=col.max,data=NA)
  niche.matrix[,1] <- 1
  
  crypt.pop <- matrix(nrow=1, ncol=col.max, data=NA)
  crypt.pop[1,1] <- N.tot-N

}

if(human == 1){
  lambda.0 <- round(1/7,3)  ##initial division rate of stem cells /
  N <- 20 ##number of stem cells in the stem cell niche
  N.tot <- 36
  nu.0 <- 0.321#round((lambda.0*(N + (N.tot - N)))/(N.tot-N),3)
  lifetime  <- 100 #lifetime in years
  
  col.max <- 2e6
  
  niche.matrix <- matrix(nrow=N,ncol=col.max,data=NA)
  niche.matrix[,1] <- 1
  
  crypt.pop <- matrix(nrow=1, ncol=col.max, data=NA)
  crypt.pop[1,1] <- N.tot-N
  
  
}

mu <- 2*6.3e-5   ##mutation rate per division

P.B <- 0.0575 ##probability of a beneficial mutation (vs. deleterious) 
s.b <- 0.061 #expected value given a beneficial mutation occured
s.d <- 0.217 #expected value given a deleterious mutation occured


#niche.matrix <- matrix(nrow=N,ncol=1,data=1)

lineage.matrix <- matrix(nrow=2, ncol=1, data=c(lambda.0,nu.0))

#crypt.pop <- matrix(nrow=1, ncol=1, data=(N.tot-N))

#initialize time
t <- rep(NA, length=col.max)
t[1] <- 0
counter <- 1
#cur_t  <-  t[length(t)]
t_final  <- lifetime*365.25 #days to run sim for (unless the ppopulation gets too big)
max_N.tot <- 1e3
cur_N.tot <- N.tot

niche.replace.fun <- function(n){
  if(n == 1){
    if(runif(1) < 0.5){
      return(N)
    }else{
      return(2)
    }
  }
  if(n == N){
    if(runif(1) < 0.5){
      return(N-1)
    }else{
      return(1)
    }
  }
  if(runif(1) < 0.5){
    return(n+1)
  }else{
    return(n-1)
  }
}

mut.fun <- function(lam){
  if(runif(1) < P.B){
    return(lam+lam*rexp(1,1/s.b))
  }else{
    new_lam <- lam-lam*rexp(1,1/s.d)
    if(new_lam > 0){
      return(new_lam)
    }else{
      return(0)
    }
  }
}


niche_rates <- rep(lambda.0, N)
crypt_rates_div <- rep(lambda.0,length(crypt.pop[,1]))
crypt_rates_diff <- rep(nu.0,length(crypt.pop[,1]))


while(t[counter] < t_final & cur_N.tot < max_N.tot ){
  
  #define total new rate 
  #all division in niche, all division in crypt.pop, all differentiation in crypt pop
  
  #niche rates, 
  for(i in 1:N){
    niche_rates[i] <- lineage.matrix[1,niche.matrix[i,counter]]
  }
  
  #Pop.2 division rates
  for(i in 1:ncol(lineage.matrix)){
    crypt_rates_div[i] <- crypt.pop[i,counter]*lineage.matrix[1,i]
  }
  
  #Pop.2 differentiation rates
  for(i in 1:ncol(lineage.matrix)){
    crypt_rates_diff[i] <- crypt.pop[i,counter]*lineage.matrix[2,i]
  }
  
  new_rate <- sum(niche_rates,crypt_rates_div,crypt_rates_diff)
  
  new_wait <- rexp(1,new_rate)
  
  choice <- sample(x=seq(1, (N+2*ncol(lineage.matrix)),1),
                   replace=T, prob=c(niche_rates/new_rate,crypt_rates_div/new_rate,crypt_rates_diff/new_rate))[1]
  
  counter <- counter+1
  
  
  t[counter] <- t[counter-1]+ new_wait  #new t of this event occuring
  niche.matrix[,counter] <- niche.matrix[,counter-1]   #prep new time  size for manipulation below
  crypt.pop[,counter] <- crypt.pop[,counter-1]
  
  #cur_t <- length(t)
  
  
  
  
  #if the event occurs in the stem cell niche                 
  if(choice <= N){
    #a certain cell divides, and replaces neighbor
    neighbor <- niche.replace.fun(choice)
    #replacement cell might have a mutation that forms a new lineage
    if(runif(1) < mu){
      lineage.matrix <- cbind(lineage.matrix,lineage.matrix[,niche.matrix[choice]]) #new lineage recorded
      lineage.matrix[1,ncol(lineage.matrix)] <- mut.fun(lineage.matrix[1,ncol(lineage.matrix)]) #new lambda based off of old lambda
      
      crypt.pop <- rbind(crypt.pop, rep(0,ncol(crypt.pop))) #new lineage spot in crypt recordings
      
      crypt.pop[niche.matrix[neighbor,counter],counter] <- crypt.pop[niche.matrix[neighbor,counter],counter]+1 #pushed out neighbor joins 2nd population
      niche.matrix[neighbor,counter] <- ncol(lineage.matrix) #replacement neighbor of THE new lineage
      
      
      
    }else{  #if there is no mutation... 
      
      crypt.pop[niche.matrix[neighbor,counter],counter] <- crypt.pop[niche.matrix[neighbor,counter],counter]+1 #pushed out neighbor joins 2nd population
      
      niche.matrix[neighbor,counter] <- niche.matrix[choice,counter] #neighbor replaced by choice lineage
      
    }
    
    
  }
  
  ##If the event is an individual in pop.2 divides...
  if(choice > N & choice <= N+(length(crypt_rates_div))){
    choice <- choice-N #find the lineage number 
    
    if(runif(1) > mu){
      crypt.pop[choice, counter] <- crypt.pop[choice, counter]+1
    }else{
      lineage.matrix <- cbind(lineage.matrix,lineage.matrix[,choice]) #new lineage recorded
      lineage.matrix[1,ncol(lineage.matrix)] <- mut.fun(lineage.matrix[1,ncol(lineage.matrix)]) #new lambda based off of old lambda
      
      #adds a new spot for the lineage and an individual
      crypt.pop <- rbind(crypt.pop, rep(0,ncol(crypt.pop))) #new lineage spot in crypt recordings
      crypt.pop[nrow(crypt.pop),counter] <- 1
      
    }
  }
  
  if(choice > N+(length(crypt_rates_div))){
    choice <- choice-(N+length(crypt_rates_div))
    
    crypt.pop[choice,counter] <- crypt.pop[choice,counter]-1
    
  }
  
  cur_N.tot <- N+sum(crypt.pop[,counter])
  if(counter%%1e4==0){
    print(paste("Percent of simulation completed: ", round(((t[counter])/t_final),4)*100))
  }
  
}

niche.matrix <- niche.matrix[,1:counter]
crypt.pop <- crypt.pop[,1:counter]
t <- t[1:counter]


# which(niche.matrix==2)

# plot(crypt.pop[1,],x=t,type="l")
# plot(crypt.pop[2,],x=t,type="l")
#crypt.pop[2,1000:2000]
# summary(crypt.pop[1,])
print("Lineages in the niche at last time point:");print(niche.matrix[,counter])

counter

print(paste("Number of lineages by the end of the simulation:",length(lineage.matrix[1,])))
print("Division Rates of all lineages created");print(lineage.matrix[1,])

image(niche.matrix,col=rainbow(length(lineage.matrix[1,])))

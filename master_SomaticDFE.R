###################
#Author: Vincent Cannataro
#R code to reproduce the data and figures for the MS on DFE in somatic tissue
################### 

#Run these first
require("caTools")
require("gplots")
#rewrite this to set the working directory to where you want it 
# setwd("C:/Users/Vincent/Dropbox/Modeling tumor evolution/code/Prob_Cancer/ForPub")
polyp <- read.csv("polyp.csv")
polyp2 <- polyp[,2]

#Note: running the human posterior density / incidence codes will erase the previous mouse data stored densities in the workspace (and vice versa) 
##################################################################################
########## 

###
#When mutations affect division rate
###

#Evolution and tumor incidence in mice 

#DFE Parameters 
P.B <- 0.0575 ##probability of a beneficial mutation (vs. deleterious) 
s.b <- 0.061 #expected value given a beneficial mutation occured
s.d <- 0.217 #expected value given a deleterious mutation occured
mu <- 2*6.3e-5   ##mutation rate per division
##########

##### Mouse initial parameters 
lambda.0 <- 0.20  ##initial division rate of stem cells 
N <- 6 ##number of stem cells in the stem cell niche
N.tot <- 15
nu.0 <- 0.333 #The rate at which stem cells commit to differentiation 
crypts <- 12e5
#####

#Code to produce the posterior distributions of division rate given fixed mutation
#and incidence curves
source("Mouse_somaticevolution.R")

#####
#Plotting posterior densities and incidence curves produced in the sourced code above
#####

#Posterior Densities 
#density.plotter(1 if you want to see exponential scenario, 1 if you want to see heavy-tailed, lower x bound of the plot, upper x bound,lower y, upper y, how many densities to plot)
par(mar=c(4,6,4,5)+.1)
density.plotter(1,0,xlim1=0,xlim2=.4,ylim1=0,ylim2=35,8)
abline(v=nu.0,lwd=3,col="red")
abline(v=0);abline(h=0)
#Posterior Densities zoomed in
density.plotter(1,0,xlim1=0.315,xlim2=.34,ylim1=0,ylim2=1.5e-3,8)
abline(v=nu.0,lwd=3,col="red")
abline(v=0);abline(h=0)

#Incidence curves 
par(mar=c(5,7,4,5)+.1)
plot((1-no.cancer.vec.powerdiv)*1e2,lty=3,lwd=5,type="l", xlab="Age (weeks)", ylab="Percent of the population with 
     an initiated tumor",main="Population tumor incidence in mice",ylim=c(0,1e2),cex.lab=1.5,cex.axis=1.3)
lines((1-no.cancer.vec.expdiv)*1e2,lty=2,col="black", lwd=5)
legend("topright", c("Pareto","Exponential            "), lty=c(3,2),lwd=3 )
abline(v=0);abline(h=0)

#####"Mutational profile" plot
par(mar=c(5,7,4,5)+.1)
barplot(mut.dist.mat,horiz=F,main="Mutation profile at the onset of tumorigenesis:
        Mouse, exponential DFE on division rate model",axes=T,axisnames=T,names.arg=week.vec,xlab="Age (Weeks)",col=rainbow(iter),ylab="Probability each mutation caused the tumor",cex.axis=1.5,cex.lab=1.5)
legend("bottomright",legend=1:3,fill=rainbow(iter))


###################################################################################
########## 

#Evolution and tumor incidence in humans

#DFE Parameters 
P.B <- 0.0575 ##probability of a beneficial mutation (vs. deleterious) 
s.b <- 0.061 #expected value given a beneficial mutation occured
s.d <- 0.217 #expected value given a deleterious mutation occured
mu <- 2*6.3e-5   ##mutation rate per division
##########

##### Human initial parameters 
lambda.0 <- round(1/7,3)  ##initial division rate of stem cells /
N <- 20 ##number of stem cells in the stem cell niche
N.tot <- 36
nu.0 <- 0.321 #the rate that stem cells commit to differentiation 
crypts <- 2e7 #estimated crypts in the large intestine 
#####

#Code to produce the posterior distributions of division rate given fixed mutation
#and incidence curves
source("Human_somaticevolution.R")

#####
#Plotting posterior densities and incidence curves produced in the sourced code above
#####

#Posterior Densities 
#density.plotter(1 if you want to see exponential scenario, 1 if you want to see heavy-tailed, lower x bound of the plot, upper x bound,lower y, upper y, how many densities to plot)
par(mar=c(4,6,4,5)+.1)
density.plotter(1,0,xlim1=0,xlim2=.4,ylim1=0,ylim2=35,8)
abline(v=nu.0,lwd=3,col="red")
abline(v=0);abline(h=0)
#Posterior densities zoomed in
density.plotter(1,0,xlim1=0.315,xlim2=.34,ylim1=0,ylim2=8e-3,8)
abline(v=nu.0,lwd=3,col="red")
abline(v=0);abline(h=0)


plot(tumor.init.expdiv,type="l", xlab="Fixed Mutation Number", lwd=5,ylab="Log(Probability of Tumorigenesis 
     PER Fixed Mutation)",cex.lab=1.5,log="y",
     main="Probabilities of Tumorigenesis Per Fixed Mutation for 
     Division Rate Scenario",cex.axis=1.25)

plot(prob.cancer.expdiv, xlab="Fixed Mutation Number", cex.lab=1.5, ylab="Log(Cumulative Probability of Tumorigenesis)", type="l", lwd=5,cex.axis=1.25,log="y",main="The Cumulative Probability of Tumorigenesis
     Division Rate Scenario")
plot(prob.cancer.expdiv, xlab="Fixed Mutation Number", cex.lab=1.5, ylab="Cumulative Probability of Tumorigenesis", type="l", lwd=5,cex.axis=1.25,log="",main="The Cumulative Probability of Tumorigenesis
     Division Rate Scenario")


#Incidence curves 
par(mar=c(5,7,4,5)+.1)
plot((1-no.cancer.vec.powerdiv)*1e2,lwd=5,type="l",lty=3, xlab="Age (years)", ylab="Percent of the population with 
     an initiated tumor",main="Population tumor incidence in humans",ylim=c(0,100),cex.lab=1.5,cex.axis=1.3)
lines((1-no.cancer.vec.expdiv)*1e2,lty=2,col="black", lwd=5)
points(polyp2,x=polyp[,1],type="o",lwd=4,col="red")
abline(v=0);abline(h=0)
legend("bottomright", c("Pareto","Exponential            ","Incidence"), lty=c(3,2,1),lwd=3,col=c("black","black","red"))


####"Mutation profile" plot
par(mar=c(5,7,4,5)+.1)
barplot(mut.dist.mat,horiz=F,main="Mutation profile at the onset of tumorigenesis:
        Human, exponential DFE on division rate model",axes=T,axisnames=T,names.arg=year.vec,xlab="Age (Years)",col=rainbow(iter),ylab="Probability each mutation caused the tumor",cex.lab=1.5,cex.axis=1.5)
legend("bottomright",legend=1:7,fill=rainbow(iter))




###Expected value figs --> make sure to run the above scripts first
#############
#Here, the mouse plot has the lambda.0=0.20 value, since we 
#just ran the human scenario and lambda.0 differs in that scenario
plot(c(0.2,expected.vec.mouseexp)/0.2,x=0:8, ylab="Expected value of division rate 
     divided by the original division rate",xlab="Fixed mutation number", cex.lab=1.5,cex.axis=1.5,main="Ratio of the expected value of the division rate 
     to the original division rate given fixed mutations in mice and humans",cex=2,pch=19)
points(c(lambda.0,expected.vec.humanexp[1:8])/lambda.0,x=0:8,pch=18,cex=2,col="red")
legend("bottomleft", c("Mice","Humans      "),cex=2, pch=c(19,18),col=c("black","red"))
#############



##################################################################################
########## 

###
#When mutations affect differentiation rate
#These scripts take longer than those above,
#because the mutational space is larger (no tumorigenesis threshold
#cutting off densities for large nu)
###


#Evolution and tumor incidence in mice 

#DFE Parameters 
P.B <- 0.0575 ##probability of a beneficial mutation (vs. deleterious) 
s.b <- 0.061 #expected value given a beneficial mutation occured
s.d <- 0.217 #expected value given a deleterious mutation occured
mu <- 2*6.3e-5   ##mutation rate per division
##########

##### Mouse initial parameters 
lambda.0 <- 0.20  ##initial division rate of stem cells 
N <- 6 ##number of stem cells in the stem cell niche
N.tot <- 15
nu.0 <- 0.333 #The rate at which stem cells commit to differentiation 
crypts <- 12e5
#####

#Code to produce the posterior distributions of division rate given fixed mutation
#and incidence curves
source("Mouse_somaticevolution_DIFF.R")

#####
#Plotting posterior densities and incidence curves produced in the sourced code above
#####

#Posterior Densities 
#density.plotter(1 if you want to see exponential scenario, 1 if you want to see heavy-tailed, lower x bound of the plot, upper x bound,lower y, upper y, how many densities to plot)
par(mar=c(5,7,4,5)+.1)
density.plotter(1,xlim1=0.1,xlim2=2,ylim1=0,ylim2=15,8)
abline(v=lambda.0,lwd=3,col="red")
abline(v=nu.0)
abline(v=0);abline(h=0)
#Posterior Densities zoomed in
density.plotter(1,xlim1=0.19,xlim2=.21,ylim1=0,ylim2=0.5e-2,8)
abline(v=lambda.0,lwd=3,col="red")
abline(v=0);abline(h=0)

#Incidence curves 
par(mar=c(5,7,4,5)+.1)
plot((1-no.cancer.vec.expdiff)*1e2,lty=3,lwd=5,type="l", xlab="Age (weeks)", ylab="Percent of the population with 
     an initiated tumor",main="Population tumor incidence in mice",ylim=c(0,1e2),cex.lab=1.5,cex.axis=1.3)
# lines((1-no.cancer.vec.expdiv)*1e2,lty=2,col="black", lwd=5)
# legend("topright", c("Pareto","Exponential            "), lty=c(3,2),lwd=3 )
abline(v=0);abline(h=0)

###"mutational profile" plot
par(mar=c(5,7,4,5)+.1)
barplot(mut.dist.mat,horiz=F,main="Mutation profile at the onset of tumorigenesis:
        Mouse, exponential DFE on differentiation rate model",axes=T,axisnames=T,names.arg=week.vec,xlab="Age (Weeks)",col=rainbow(iter),ylab="Probability each mutation caused the tumor",cex.axis=1.5,cex.lab=1.5)
legend("bottomright",legend=1:3,fill=rainbow(iter))


###################################################################################
########## 

#Evolution and tumor incidence in humans

#DFE Parameters 
P.B <- 0.0575 ##probability of a beneficial mutation (vs. deleterious) 
s.b <- 0.061 #expected value given a beneficial mutation occured
s.d <- 0.217 #expected value given a deleterious mutation occured
mu <- 2*6.3e-5   ##mutation rate per division
##########

##### Human initial parameters 
lambda.0 <- round(1/7,3)  ##initial division rate of stem cells /
N <- 20 ##number of stem cells in the stem cell niche
N.tot <- 36
nu.0 <- 0.321 #the rate that stem cells commit to differentiation 
crypts <- 2e7 #estimated crypts in the large intestine 
#####

#Code to produce the posterior distributions of division rate given fixed mutation
#and incidence curves
source("Human_somaticevolution_DIFF.R")

#####
#Plotting posterior densities and incidence curves produced in the sourced code above
#####

#Posterior Densities 
#density.plotter(1 if you want to see exponential scenario, 1 if you want to see heavy-tailed, lower x bound of the plot, upper x bound,lower y, upper y, how many densities to plot)
par(mar=c(5,7,4,5)+.1)
density.plotter(1,xlim1=0.1,xlim2=2,ylim1=0,ylim2=15,8)
abline(v=nu.0,lwd=1,col="black")
abline(v=0);abline(h=0)
#Posterior densities zoomed in
density.plotter(1,xlim1=0.133,xlim2=.153,ylim1=0,ylim2=5e-4,8)
abline(v=lambda.0,lwd=3,col="red")
abline(v=0);abline(h=0)





#Incidence curves 
par(mar=c(5,7,4,5)+.1)
plot((1-no.cancer.vec.expdiff)*1e2,lwd=5,type="l",lty=3, xlab="Age (years)", ylab="Percent of the population with 
     an initiated tumor",main="Population tumor incidence in humans",ylim=c(0,100),cex.lab=1.5,cex.axis=1.3)
# lines((1-no.cancer.vec.expdiv)*1e2,lty=2,col="black", lwd=5)
points(polyp2,x=polyp[,1],type="o",lwd=4,col="red")
abline(v=0);abline(h=0)
# legend("bottomright", c("Pareto","Exponential            ","Incidence"), lty=c(3,2,1),lwd=3,col=c("black","black","red"))


###"mutation profile" plots
barplot(mut.dist.mat,horiz=F,main="Mutation profile at the onset of tumorigenesis:
        Human, exponential DFE on differentiation rate model",axes=T,axisnames=T,names.arg=year.vec,xlab="Age (Years)",col=rainbow(iter),ylab="Probability each mutation caused the tumor",cex.axis=1.5,cex.lab=1.5)
legend("bottomright",legend=1:7,fill=rainbow(iter))

#######################################################
##################
#####

####
#Expected values plotter for 
####
plot(c(0.333,mouse.expdiff.meannu[1:8])/0.333,x=0:8, ylab="Expected value of differentiation rate 
divided by the original differentiation rate",xlab="Fixed mutation number", cex.lab=1.5,cex.axis=1.5,main="Ratio of the expected value of the differentiation rate 
     to the original differentiation rate given fixed mutations in mice and humans",cex=2,pch=19)
points(c(0.321,human.expdiff.meannu[1:8])/0.321,x=0:8,pch=18,cex=2,col="red")
legend("topleft", c("Mice","Humans"),cex=2, pch=c(19,18),col=c("black","red"))




#######################################################
##################
#####

##Exploring parameter space for least squares analysis

#Generating and storing incidence curves 

##########
#Evolution and tumor incidence in humans
#Run these parameters first before sourcing ParameterExplore_....R

#DFE Parameters 
P.B <- 0.0575 ##probability of a beneficial mutation (vs. deleterious) 
s.b <- 0.061 #expected value given a beneficial mutation occured
s.d <- 0.217 #expected value given a deleterious mutation occured
mu <- 2*6.3e-5   ##mutation rate per division
##########

##### Human initial parameters 
lambda.0 <- round(1/7,3)  ##initial division rate of stem cells /
N <- 20 ##number of stem cells in the stem cell niche
N.tot <- 36
nu.0 <- 0.321 #the rate that stem cells commit to differentiation 
crypts <- 2e7 #estimated crypts in the large intestine 
#####

#Parameter Space to explore 
E.s.plus.vec <- seq(0.07,0.041,-0.001)
mu.vec <- c(seq(0.000025,0.000125,0.000025),0.000126,seq(0.00015,0.0007,0.000025))
# E.s.plus.vec <- seq(0.041,0.07,0.001) #Expected beneficial effect (s.b)
# mu.vec <- seq(0.000025,0.0007,0.000025) #Mutation rate (mu)
#########
# E.s.plus.vec <- seq(0.041,0.07,0.001) #Expected beneficial effect (s.b)
# mu.vec <- seq(0.000126,0.000126,0.000025) #Mutation rate (mu)
# #########

###
#Sourcing these files below will create a new folder in the working directory
#and store a CSV for every incidence curve for every parameter choice
#these are retrieved later for the least squares analysis 
#(saved as individual files in case code needs to be terminated and
#can be restarted to fill in the gaps)
#These scripts take a while to run, make sure you are upgraded to newest version
#of R to prevent memory allocation errors. 
#
#Every time you run these files they overwrite the previous saved incidence data files
###
source("ParameterExplore_human_exptail.R")

source("ParameterExplore_human_heavytail.R")




############################
#Analyzing the incidence curves generated above

# E.s.plus.vec <- seq(0.07,0.041,-0.001)
# mu.vec <- seq(0.000025,0.0007,0.000025)

# E.s.plus.vec <- seq(0.041,0.07,0.001) #Expected beneficial effect (s.b)
# mu.vec <- seq(0.000126,0.000126,0.000025) #Mutation rate (mu)

E.s.plus.vec <- seq(0.07,0.041,-0.001)
mu.vec <- c(seq(0.000025,0.000125,0.000025),0.000126,seq(0.00015,0.0007,0.000025))


#####Sourcing these files causes a large heatmap to be produced,
##If using RStudio you should expand the plot window 

#Outputs the 
source("ParameterAnalyze_human_exptail.R")



# E.s.plus.vec <- seq(0.07,0.041,-0.001)
# mu.vec <- seq(0.000025,0.0007,0.000025)

source("ParameterAnalyze_human_heavytail.R")

###
#For the differentiation rate case
###





#######################################################
##################
#####

##Exploring parameter space for least squares analysis

#Generating and storing incidence curves 

##########
#Evolution and tumor incidence in humans
#Run these parameters first before sourcing ParameterExplore_....R

#DFE Parameters 
P.B <- 0.0575 ##probability of a beneficial mutation (vs. deleterious) 
s.b <- 0.061 #expected value given a beneficial mutation occured
s.d <- 0.217 #expected value given a deleterious mutation occured
mu <- 2*6.3e-5   ##mutation rate per division
##########

##### Human initial parameters 
lambda.0 <- round(1/7,3)  ##initial division rate of stem cells /
N <- 20 ##number of stem cells in the stem cell niche
N.tot <- 36
nu.0 <- 0.321 #the rate that stem cells commit to differentiation 
crypts <- 2e7 #estimated crypts in the large intestine 
#####

#Parameter Space to explore 


E.s.plus.vec <- seq(0.07,0.041,-0.001)
mu.vec <- c(seq(0.00000025,0.0000025,0.00000025),seq(0.0000050,0.0000125,0.0000025),seq(0.000025,0.000125,0.000025),0.000126,seq(0.00015,0.0007,0.000025))

###
#Sourcing these files below will create a new folder in the working directory
#and store a CSV for every incidence curve for every parameter choice
#these are retrieved later for the least squares analysis 
#(saved as individual files in case code needs to be terminated and
#can be restarted to fill in the gaps)
#These scripts take a while to run, make sure you are upgraded to newest version
#of R to prevent memory allocation errors. 
#
#Every time you run these files they overwrite the previous saved incidence data files
###
source("Human_parameterexplore_DIFF.R")

# source("ParameterExplore_human_heavytail.R")




############################
#Analyzing the incidence curves generated above

E.s.plus.vec <- seq(0.07,0.041,-0.001)
mu.vec <- c(seq(0.00000025,0.0000025,0.00000025),seq(0.0000050,0.0000125,0.0000025),seq(0.000025,0.000125,0.000025),0.000126,seq(0.00015,0.0007,0.000025))
#####Sourcing these files causes a large heatmap to be produced,
##If using RStudio you should expand the plot window 

#Outputs the 
source("ParameterAnalyze_human_exptail_DIFF.R")








############################################################################
##########################

#Gillespie algorithm/method to simulate mutation accumulation and niche/crypt dynamics.


#1 for the organism to simulate a crypt in, 0 for the other
#NOTE: Human lifetime set to 100 years, and many more stem cells, so it takes much longer
mouse <- 0
human <- 1


source("crypt_dynamics_Gillespie.R")












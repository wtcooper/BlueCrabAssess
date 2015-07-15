rm(list=ls())

#SOURCE the functions
source("mcmc_funs.r")

#file name output by the model
mcmcFile="cmsa_refs.mcmc"
mcmcFile1="cmsa_refs1.mcmc"
mcmcFile2="cmsa_refs2.mcmc"
mcmcFile3="cmsa_refs3.mcmc"

library(coda)


##################################################
#run likelihood profile to get alternative starting points for multichain run (in MSY.pvl)
##################################################
#this will give estimate of how long to run for convergence
#take alternative starting points in the chains and put them in the parsX.txt files for later 
system("cmsa.exe -lprof -prsave", wait=T)




##################################################
#Single MCMC run to get burn-in, #runs, and thinning rate
##################################################
library(coda)

runMCMC("cmsa.exe",1000000,1)
Sys.time()


#############
#### East ####
#############

east=read.table("cmsa_refs_east.mcmc", sep="\t", header=T)

burn=1000
thin=1000

se=seq(burn,length(east[,1]),by=thin)
chain=mcmc(east[se,])

#plot(chain)
#raftery.diag(chain,r=0.01)

summary(chain)


east=read.table("cmsa_yearly_east.mcmc", sep="\t", header=T)

burn=1000
thin=1000

se=seq(burn,length(east[,1]),by=thin)
chain=mcmc(east[se,])

#plot(chain)
raftery.diag(chain,r=0.01)

summary(chain)

#############
#### West ####
#############

west=read.table("cmsa_refs_west.mcmc", sep="\t", header=T)

burn=1000
thin=1000

se=seq(burn,length(west[,1]),by=thin)
chain=mcmc(west[se,])

#plot(chain)
#raftery.diag(chain,r=0.01)

summary(chain)


west=read.table("cmsa_yearly_west.mcmc", sep="\t", header=T)

burn=1000
thin=1000

se=seq(burn,length(west[,1]),by=thin)
chain=mcmc(west[se,])

#plot(chain)
raftery.diag(chain,r=0.01)

summary(chain)


################
## Old 
################

runMCMC("cmsa.exe",50000000,350)
Sys.time()

#chain.1=getMCMC(mcmcFile1)
#chain.2=getMCMC(mcmcFile2)
#chain.3=getMCMC(mcmcFile3)

chain=chain.3

windows(record=T)
plot(chain)
raftery.diag(chain,r=0.01)
#gelman.diag(chain)
#gelman.plot(chain)


### This returns as a chain
singleChain=runAndGetMCMC("cmsa.exe",mcmcFile,10000000,350)
windows(recrod=T)
plot(singleChain)
raftery.diag(singleChain,r=0.01)

#Output:
#Quantile (q) = 0.025
#Accuracy (r) = +/- 0.01
#Probability (s) = 0.95 

#Burn-in  Total  Lower bound  Dependence
#(M)      (N)    (Nmin)       factor (I)
#MSY   198      50886  937           54.3     
#FMSY  896      307440 937          328.0     
#NMSY  230      56097  937           59.9     
#u0MSY 896      307440 937          328.0     
#u1MSY 896      307440 937          328.0     
#uMSY  896      307440 937          328.0     



##################################################
#Run multiple chains at different starting values for convergence check
##################################################
#NOTE: these par files give very strange results, not sure what's going on....

chain1=runAndGetMCMC_pars("cmsa.exe","pars1.txt",mcmcFile,100000,1)
chain2=runAndGetMCMC_pars("cmsa.exe","pars2.txt",mcmcFile,100000,1)
chain3=runAndGetMCMC_pars("cmsa.exe","pars3.txt",mcmcFile,100000,1)
chains=mcmc.list(chain1,chain2,chain3)
windows(record=T)
plot(chains)
gelman.diag(chains)
gelman.plot(chains)

windows(record=T)
plot(chain3)


##################################################
#Single MCMC run on final params
##################################################

singleChain=runAndGetMCMC("cmsa.exe",mcmcFile,10000,1)
summary(singleChain)
windows(record=T)
plot(singleChain)








##################################################
#Manual Single MCMC run to get burn-in, #runs, and thinning rate
##################################################
#  (1) Run ADMB manually with command options -mcmc 100000 -mcsave 1
#	(2) Run ADMB manually with command options -mceval

chain<-read.table(mcmcFile1,header=T)
plot(chain)
raftery.diag(chain,r=0.01)


#  (3) Go into MSY.pvl file, pull different starting values and put in seperate files

#  (4) Run ADMB manually with command options -mcmc 1000000 -mcsave 100 -mcpin parfilename1
#	(5) Run ADMB manually with command options -mceval

chain1<-read.table(mcmcFile,header=T)

#Repeat:
#  (4) Run ADMB manually with command options -mcmc 1000000 -mcsave 100 -mcpin parfilename2
#	(5) Run ADMB manually with command options -mceval

chain2<-read.table(mcmcFile,header=T)

#Repeat:
#  (4) Run ADMB manually with command options -mcmc 1000000 -mcsave 100 -mcpin parfilename3
#	(5) Run ADMB manually with command options -mceval

chain3<-read.table(mcmcFile,header=T)

chains=mcmc.list(chain1,chain2,chain3)
plot(chains)
gelman.diag(chains)
gelman.plot(chains)


#######################
## Direct R code that don't work.....
#######################
modelName="cmsa.exe"
numSims=10000
thin=10
mcmcArgs =paste(modelName,"-mcmc",numSims,"-mcsave",thin)
mcmcEvalArgs = paste(modelName," -mceval")

system(mcmcArgs, wait=T)
x=1
system(mcmcEvalArgs,wait=T)
x


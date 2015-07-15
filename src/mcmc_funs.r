
##################################################
##Functions for running MCMC
##################################################

library(coda)

runMCMC=function(modelName, numSims, thin) {
	
	mcmcArgs =paste(modelName,"-mcmc",numSims,"-mcsave",thin)
	mcmcEvalArgs = paste(modelName," -mceval")
	
	#run the admb model mcmc
	system(mcmcArgs, wait=T)
	
	#run -mceval
	system(mcmcEvalArgs,wait=T)
	
} #end function

getMCMC=function(mcmcFileName) {
	#get chain output
	chains<-read.table(mcmcFileName,header=T)
	return(mcmc(chains))
} #end function


runAndGetMCMC=function(modelName, mcmcFileName, numSims, thin) {
	
	mcmcArgs =paste(modelName,"-mcmc",numSims,"-mcsave",thin)
	mcmcEvalArgs = paste(modelName," -mceval")
	
	#run the admb model mcmc
	system(mcmcArgs, wait=T)
	
	#run -mceval
	system(mcmcEvalArgs,wait=T)
	
	#get chain output
	chains<-read.table(mcmcFileName,header=T)
	return(mcmc(chains))
} #end function



#runs with a parameter input file
runAndGetMCMC_pars=function(modelName, parFileName, mcmcFileName, numSims, thin) {
	
	mcmcArgs =paste(modelName,"-mcmc",numSims,"-mcsave",thin,"-mcpin",parFileName)
	mcmcEvalArgs = paste(modelName,"-mceval")
	
	#run the admb model mcmc
	system(mcmcArgs, wait=T)
	
	#run -mceval
	system(mcmcEvalArgs,wait=T)
	
	#get chain output
	chains<-read.table(mcmcFileName,header=T)
	return(mcmc(chains))
} #end function
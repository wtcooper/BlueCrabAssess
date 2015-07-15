rm(list=ls())

##SOURCE the functions
#source("mcmc_funs.r")


library(coda)


##################################################
#run likelihood profile to get alternative starting points for multichain run (in MSY.pvl)
##################################################
#this will give estimate of how long to run for convergence
#take alternative starting points in the chains and put them in the parsX.txt files for later 
#system("cmsa.exe -lprof -prsave", wait=T)
#Sys.time()


###################################################
##Single MCMC run to get burn-in, #runs, and thinning rate
###################################################

#

reps=100000
thin=1000
seed=12345
seed=15243

#Do west
setwd("./West")
#With initial parameter stating values: from likelihood profile above
#mcmcArgs =paste("cmsa.exe","-mcmc",reps,"-mcsave",thin,"-mcseed",seed,"-mcpin","par_west.txt")
#Without initial parameter starting values
mcmcArgs =paste("cmsa.exe","-mcmc",reps,"-mcsave",thin,"-mcseed",seed)
mcmcEvalArgs = paste("cmsa.exe"," -mceval")
system(mcmcArgs, wait=T)
system(mcmcEvalArgs,wait=T)


#Do east
setwd("./East")
#mcmcArgs =paste("cmsa.exe","-mcmc",reps,"-mcsave",thin,"-mcseed",seed,"-mcpin","par_east.txt")
mcmcArgs =paste("cmsa.exe","-mcmc",reps,"-mcsave",thin,"-mcseed",seed)
mcmcEvalArgs = paste("cmsa.exe"," -mceval")
system(mcmcArgs, wait=T)
system(mcmcEvalArgs,wait=T)

Sys.time()



###################################################
##Analysis of chains: plots and summary
###################################################


#############
#### East ####
#############
library(coda)

setwd("./East")
east=read.table("cmsa_refs.mcmc", sep="\t", header=T)

burn=2
thin=1

se=seq(burn,length(east[,1]),by=thin)
chain=mcmc(east[se,])

windows(record=T, width=6.5, height=8)
par(mfrow=c(5,2), mar=c(4,4,2,2), oma=c(0,0,0,0), 	
		mgp = c(1.25, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
		tck = -0.01, 		#size of tick mark; neg is outside graph
		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
		cex.lab=1, 		#size of axis title
		cex.axis=.9,		#size of axis labels
		bty="l")			#no box, just axes 


plot(chain[,c(1,2,3,5,6)], auto.layout = FALSE, ask=FALSE)  
#raftery.diag(chain,r=0.01)

summary(chain)



east=read.table("cmsa_pars.mcmc", sep="\t", header=T)
nms=names(east)
nms[1]="InitialN"; nms[2]="InitialR"; nms[3]="InitialF";
names(east)=nms

burn=2
thin=1

se=seq(burn,length(east[,1]),by=thin)
chain=mcmc(east[se,])

windows(record=T, width=6.5, height=8)
par(mfrow=c(5,2), mar=c(4,4,2,2), oma=c(0,0,0,0), 	
		mgp = c(1.25, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
		tck = -0.01, 		#size of tick mark; neg is outside graph
		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
		cex.lab=1, 		#size of axis title
		cex.axis=.9,		#size of axis labels
		bty="l")			#no box, just axes 


plot(chain, auto.layout = FALSE, ask=FALSE)  
#raftery.diag(chain,r=0.01)

summary(chain)



east=read.table("cmsa_yearly.mcmc", sep="\t", header=T)

burn=2
thin=1

se=seq(burn,length(east[,1]),by=thin)
chain=mcmc(east[se,])

#plot(chain)
raftery.diag(chain,r=0.01)

summary(chain)

rm(east,chain)


#############
#### West ####
#############

setwd("./West")

west=read.table("cmsa_refs.mcmc", sep="\t", header=T)

burn=2
thin=1

se=seq(burn,length(west[,1]),by=thin)
chain=mcmc(west[se,])

windows(record=T, width=6.5, height=8)
par(mfrow=c(5,2), mar=c(4,4,2,2), oma=c(0,0,0,0), 	
		mgp = c(1.25, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
		tck = -0.01, 		#size of tick mark; neg is outside graph
		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
		cex.lab=1, 		#size of axis title
		cex.axis=.9,		#size of axis labels
		bty="l")			#no box, just axes 


plot(chain[,c(1,2,3,5,6)], auto.layout = FALSE, ask=FALSE)  
#raftery.diag(chain,r=0.01)

summary(chain)

west=read.table("cmsa_pars.mcmc", sep="\t", header=T)
nms=names(west)
nms[1]="InitialN"; nms[2]="InitialR"; nms[3]="InitialF";
names(west)=nms
burn=2
thin=1

se=seq(burn,length(west[,1]),by=thin)
chain=mcmc(west[se,])

windows(record=T, width=6.5, height=8)
par(mfrow=c(5,2), mar=c(4,4,2,2), oma=c(0,0,0,0), 	
		mgp = c(1.25, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
		tck = -0.01, 		#size of tick mark; neg is outside graph
		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
		cex.lab=1, 		#size of axis title
		cex.axis=.9,		#size of axis labels
		bty="l")			#no box, just axes 


plot(chain, auto.layout = FALSE, ask=FALSE)  
#raftery.diag(chain,r=0.01)

summary(chain)


west=read.table("cmsa_yearly.mcmc", sep="\t", header=T)

burn=2
thin=1

se=seq(burn,length(west[,1]),by=thin)
chain=mcmc(west[se,])

#plot(chain)
raftery.diag(chain,r=0.01)

summary(chain)

rm(west,chain)



#############
#### Atlantic ####
#############
library(coda)

setwd("./Atl")
atl=read.table("cmsa_refs.mcmc", sep="\t", header=T)

burn=2
thin=1

se=seq(burn,length(atl[,1]),by=thin)
chain=mcmc(atl[se,])

windows(record=T, width=6.5, height=8)
par(mfrow=c(5,2), mar=c(4,4,2,2), oma=c(0,0,0,0), 	
		mgp = c(1.25, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
		tck = -0.01, 		#size of tick mark; neg is outside graph
		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
		cex.lab=1, 		#size of axis title
		cex.axis=.9,		#size of axis labels
		bty="l")			#no box, just axes 


plot(chain[,c(1,2,3,5,6)], auto.layout = FALSE, ask=FALSE)  
#raftery.diag(chain,r=0.01)

summary(chain)



atl=read.table("cmsa_pars.mcmc", sep="\t", header=T)
nms=names(atl)
nms[1]="InitialN"; nms[2]="InitialR"; nms[3]="InitialF";
names(atl)=nms

burn=2
thin=1

se=seq(burn,length(atl[,1]),by=thin)
chain=mcmc(atl[se,])

windows(record=T, width=6.5, height=8)
par(mfrow=c(5,2), mar=c(4,4,2,2), oma=c(0,0,0,0), 	
		mgp = c(1.25, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
		tck = -0.01, 		#size of tick mark; neg is outside graph
		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
		cex.lab=1, 		#size of axis title
		cex.axis=.9,		#size of axis labels
		bty="l")			#no box, just axes 


plot(chain, auto.layout = FALSE, ask=FALSE)  
#raftery.diag(chain,r=0.01)

summary(chain)



atl=read.table("cmsa_yearly.mcmc", sep="\t", header=T)

burn=2
thin=1

se=seq(burn,length(atl[,1]),by=thin)
chain=mcmc(atl[se,])

#plot(chain)
raftery.diag(chain,r=0.01)

summary(chain)

rm(atl,chain)




#############
#### Gulf ####
#############
library(coda)

setwd("./Glf")
glf=read.table("cmsa_refs.mcmc", sep="\t", header=T)

burn=2
thin=1

se=seq(burn,length(glf[,1]),by=thin)
chain=mcmc(glf[se,])

windows(record=T, width=6.5, height=8)
par(mfrow=c(5,2), mar=c(4,4,2,2), oma=c(0,0,0,0), 	
		mgp = c(1.25, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
		tck = -0.01, 		#size of tick mark; neg is outside graph
		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
		cex.lab=1, 		#size of axis title
		cex.axis=.9,		#size of axis labels
		bty="l")			#no box, just axes 


plot(chain[,c(1,2,3,5,6)], auto.layout = FALSE, ask=FALSE)  
#raftery.diag(chain,r=0.01)

summary(chain)



glf=read.table("cmsa_pars.mcmc", sep="\t", header=T)
nms=names(glf)
nms[1]="InitialN"; nms[2]="InitialR"; nms[3]="InitialF";
names(glf)=nms

burn=2
thin=1

se=seq(burn,length(glf[,1]),by=thin)
chain=mcmc(glf[se,])

windows(record=T, width=6.5, height=8)
par(mfrow=c(5,2), mar=c(4,4,2,2), oma=c(0,0,0,0), 	
		mgp = c(1.25, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
		tck = -0.01, 		#size of tick mark; neg is outside graph
		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
		cex.lab=1, 		#size of axis title
		cex.axis=.9,		#size of axis labels
		bty="l")			#no box, just axes 


plot(chain, auto.layout = FALSE, ask=FALSE)  
#raftery.diag(chain,r=0.01)

summary(chain)



glf=read.table("cmsa_yearly.mcmc", sep="\t", header=T)

burn=2
thin=1

se=seq(burn,length(glf[,1]),by=thin)
chain=mcmc(glf[se,])

#plot(chain)
raftery.diag(chain,r=0.01)

summary(chain)

rm(glf,chain)






###################################################
##Plot Uncertainty values
###################################################


#############
#### East ####
#############
setwd("./East")
east=read.table("MCMC_east.dat", sep="\t", header=T)

windows(record=T, width=6.5, height=8)
par(mfrow=c(3,1), mar=c(4,4,2,2), oma=c(0,0,0,0), 	
#par(mfrow=c(3,1), mar=c(3,3,1,1), oma=c(0,0,0,0), 	
				mgp = c(1.25, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
		tck = -0.01, 		#size of tick mark; neg is outside graph
		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
		cex.lab=1, 		#size of axis title
		cex.axis=1,		#size of axis labels
		bty="l")			#no box, just axes 


#Recruits

x=east$year
y1=east$RMC
y2=east$R
ylo=east$RMCLO
yhi=east$RMCHI
polyx=c(x, rev(x))
polyy=c(ylo, rev(yhi))

plot(x, y1, type="l", lwd=2, xlab=" ",  axes = TRUE, ylab=expression(Juveniles~x10^6),  ylim=c(min(c(y1,y2,ylo))*.8,max(c(y1,y2,yhi))*1.2), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
polygon(polyx, polyy, col="gray", border=NA)
lines(x, y2, lwd=2)
lines(x, y1, lty=2, lwd=2)
legend("topleft", c("Base Best Fit", "MCMC Median"), lwd=c(2,2), lty=c(1,2))

#Adults
x=east$year
y1=east$NMC
y2=east$N
ylo=east$NMCLO
yhi=east$NMCHI
polyx=c(x, rev(x))
polyy=c(ylo, rev(yhi))

plot(x, y1, type="l", lwd=2, xlab=" ",  axes = TRUE, ylab=expression(Adults~x10^6),  ylim=c(min(c(y1,y2,ylo))*.8,max(c(y1,y2,yhi))*1.2), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
polygon(polyx, polyy, col="gray", border=NA)
lines(x, y2, lwd=2)
lines(x, y1, lty=2, lwd=2)
#legend("topleft", c("Base Best Fit", "MCMC Median"), lwd=c(2,2), lty=c(1,2))

#F
x=east$year
y1=east$FMC
y2=east$F
ylo=east$FMCLO
yhi=east$FMCHI
polyx=c(x, rev(x))
polyy=c(ylo, rev(yhi))

plot(x, y1, type="l", lwd=2, xlab=" ",  axes = TRUE, ylab="F rate",  ylim=c(min(c(y1,y2,ylo))*.8,max(c(y1,y2,yhi))*1.2), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
polygon(polyx, polyy, col="gray", border=NA)
lines(x, y2, lwd=2)
lines(x, y1, lty=2, lwd=2)
#legend("topleft", c("Base Best Fit", "MCMC Median"), lwd=c(2,2), lty=c(1,2))



#############
#### West ####
#############

setwd("./West")
west=read.table("MCMC_west.dat", sep="\t", header=T)

windows(record=T, width=6.5, height=8)
par(mfrow=c(3,1), mar=c(4,4,2,2), oma=c(0,0,0,0), 	
		mgp = c(1.25, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
		tck = -0.01, 		#size of tick mark; neg is outside graph
		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
		cex.lab=1, 		#size of axis title
		cex.axis=1,		#size of axis labels
		bty="l")			#no box, just axes 


#Recruits

x=west$year
y1=west$RMC
y2=west$R
ylo=west$RMCLO
yhi=west$RMCHI
polyx=c(x, rev(x))
polyy=c(ylo, rev(yhi))

plot(x, y1, type="l", lwd=2, xlab=" ",  axes = TRUE, ylab=expression(Juveniles~x10^6),  ylim=c(min(c(y1,y2,ylo))*.8,max(c(y1,y2,yhi))*1.2), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
polygon(polyx, polyy, col="gray", border=NA)
lines(x, y2, lwd=2)
lines(x, y1, lty=2, lwd=2)
legend("topleft", c("Base Best Fit", "MCMC Median"), lwd=c(2,2), lty=c(1,2))

#Adults
x=west$year
y1=west$NMC
y2=west$N
ylo=west$NMCLO
yhi=west$NMCHI
polyx=c(x, rev(x))
polyy=c(ylo, rev(yhi))

plot(x, y1, type="l", lwd=2, xlab=" ",  axes = TRUE, ylab=expression(Adults~x10^6),  ylim=c(min(c(y1,y2,ylo))*.8,max(c(y1,y2,yhi))*1.2), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
polygon(polyx, polyy, col="gray", border=NA)
lines(x, y2, lwd=2)
lines(x, y1, lty=2, lwd=2)
#legend("topleft", c("Base Best Fit", "MCMC Median"), lwd=c(2,2), lty=c(1,2))

#F
x=west$year
y1=west$FMC
y2=west$F
ylo=west$FMCLO
yhi=west$FMCHI
polyx=c(x, rev(x))
polyy=c(ylo, rev(yhi))

plot(x, y1, type="l", lwd=2, xlab=" ",  axes = TRUE, ylab="F rate",  ylim=c(min(c(y1,y2,ylo))*.8,max(c(y1,y2,yhi))*1.2), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
polygon(polyx, polyy, col="gray", border=NA)
lines(x, y2, lwd=2)
lines(x, y1, lty=2, lwd=2)
#legend("topleft", c("Base Best Fit", "MCMC Median"), lwd=c(2,2), lty=c(1,2))





#############
#### Gulf  ####
#############
setwd("./Glf")
east=read.table("MCMC_glf.dat", sep="\t", header=T)

windows(record=T, width=6.5, height=8)
par(mfrow=c(3,1), mar=c(4,4,2,2), oma=c(0,0,0,0), 	
		mgp = c(1.25, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
		tck = -0.01, 		#size of tick mark; neg is outside graph
		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
		cex.lab=1, 		#size of axis title
		cex.axis=1,		#size of axis labels
		bty="l")			#no box, just axes 


#Recruits

x=east$year
y1=east$RMC
y2=east$R
ylo=east$RMCLO
yhi=east$RMCHI
polyx=c(x, rev(x))
polyy=c(ylo, rev(yhi))

plot(x, y1, type="l", lwd=2, xlab=" ",  axes = TRUE, ylab=expression(Juveniles~x10^6),  ylim=c(min(c(y1,y2,ylo))*.8,max(c(y1,y2,yhi))*1.2), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
polygon(polyx, polyy, col="gray", border=NA)
lines(x, y2, lwd=2)
lines(x, y1, lty=2, lwd=2)
legend("topleft", c("Base Best Fit", "MCMC Median"), lwd=c(2,2), lty=c(1,2))

#Adults
x=east$year
y1=east$NMC
y2=east$N
ylo=east$NMCLO
yhi=east$NMCHI
polyx=c(x, rev(x))
polyy=c(ylo, rev(yhi))

plot(x, y1, type="l", lwd=2, xlab=" ",  axes = TRUE, ylab=expression(Adults~x10^6),  ylim=c(min(c(y1,y2,ylo))*.8,max(c(y1,y2,yhi))*1.2), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
polygon(polyx, polyy, col="gray", border=NA)
lines(x, y2, lwd=2)
lines(x, y1, lty=2, lwd=2)
#legend("topleft", c("Base Best Fit", "MCMC Median"), lwd=c(2,2), lty=c(1,2))

#F
x=east$year
y1=east$FMC
y2=east$F
ylo=east$FMCLO
yhi=east$FMCHI
polyx=c(x, rev(x))
polyy=c(ylo, rev(yhi))

plot(x, y1, type="l", lwd=2, xlab=" ",  axes = TRUE, ylab="F rate",  ylim=c(min(c(y1,y2,ylo))*.8,max(c(y1,y2,yhi))*1.2), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
polygon(polyx, polyy, col="gray", border=NA)
lines(x, y2, lwd=2)
lines(x, y1, lty=2, lwd=2)
#legend("topleft", c("Base Best Fit", "MCMC Median"), lwd=c(2,2), lty=c(1,2))



#############
#### Atlantic ####
#############
setwd("./Atl")
east=read.table("MCMC_atl.dat", sep="\t", header=T)

windows(record=T, width=6.5, height=8)
par(mfrow=c(3,1), mar=c(4,4,2,2), oma=c(0,0,0,0), 	
		mgp = c(1.25, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
		tck = -0.01, 		#size of tick mark; neg is outside graph
		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
		cex.lab=1, 		#size of axis title
		cex.axis=1,		#size of axis labels
		bty="l")			#no box, just axes 


#Recruits

x=east$year
y1=east$RMC
y2=east$R
ylo=east$RMCLO
yhi=east$RMCHI
polyx=c(x, rev(x))
polyy=c(ylo, rev(yhi))

plot(x, y1, type="l", lwd=2, xlab=" ",  axes = TRUE, ylab=expression(Juveniles~x10^6),  ylim=c(min(c(y1,y2,ylo))*.8,max(c(y1,y2,yhi))*1.2), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
polygon(polyx, polyy, col="gray", border=NA)
lines(x, y2, lwd=2)
lines(x, y1, lty=2, lwd=2)
legend("topleft", c("Base Best Fit", "MCMC Median"), lwd=c(2,2), lty=c(1,2))

#Adults
x=east$year
y1=east$NMC
y2=east$N
ylo=east$NMCLO
yhi=east$NMCHI
polyx=c(x, rev(x))
polyy=c(ylo, rev(yhi))

plot(x, y1, type="l", lwd=2, xlab=" ",  axes = TRUE, ylab=expression(Adults~x10^6),  ylim=c(min(c(y1,y2,ylo))*.8,max(c(y1,y2,yhi))*1.2), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
polygon(polyx, polyy, col="gray", border=NA)
lines(x, y2, lwd=2)
lines(x, y1, lty=2, lwd=2)
#legend("topleft", c("Base Best Fit", "MCMC Median"), lwd=c(2,2), lty=c(1,2))

#F
x=east$year
y1=east$FMC
y2=east$F
ylo=east$FMCLO
yhi=east$FMCHI
polyx=c(x, rev(x))
polyy=c(ylo, rev(yhi))

plot(x, y1, type="l", lwd=2, xlab=" ",  axes = TRUE, ylab="F rate",  ylim=c(min(c(y1,y2,ylo))*.8,max(c(y1,y2,yhi))*1.2), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
polygon(polyx, polyy, col="gray", border=NA)
lines(x, y2, lwd=2)
lines(x, y1, lty=2, lwd=2)
#legend("topleft", c("Base Best Fit", "MCMC Median"), lwd=c(2,2), lty=c(1,2))


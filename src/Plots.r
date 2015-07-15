
rm(list = ls())

pdfFile="Graphs.pdf"
pdfOn=0  #0=off


if (pdfOn==0) windows(record=T)
if (pdfOn==1) pdf(file=pdfFile, width=7, height=10)

###########################################################
###########################################################
#Model Output Figs
###########################################################
###########################################################

#Fix the .std file to remove a space between st & dev so reads proper #columns
input_file <- readLines("cmsa.std")
input_file[1] <- " index   name           value      stddev   "
writeLines(input_file, "cmsa.std")

#Read in data files
resids<-read.table("obs_pred_results.dat",header=TRUE,sep= "")
HPD<-read.table("HPD_results.dat",header=TRUE,sep= "", fill=TRUE)
MSY<-read.table("MSY_results.dat",header=TRUE)
gen <- read.table("gen_results.dat",header=TRUE,sep= "")
SDs<- read.table("cmsa.std",header=TRUE,sep= "")
tSPR <- read.table("tSPR.dat", header=T)

  




############################################
#Observed and predicted total catch
############################################
windows(record=T,width=6.5, height=8)
par(	#mfrow=c(3,2), 
		mar=c(3,3,1,1), oma=c(0,0,0,0), 	
		mgp = c(2, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
		tck = -0.01, 		#size of tick mark; neg is outside graph
		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
		cex.lab=1, 		#size of axis title
		cex.axis=1,		#size of axis labels
		bty="l")			#no box, just axes

plotCI=0 #turn off CI's

layout(matrix(c(1,2,3,4,5,6), nrow=3, ncol=2, byrow = TRUE), 
		widths=c(1.5,1), heights=c(1,1,1))

#select total catch data
varp<-c("t")
myvars <-c("obs")
y1 <- resids[which(resids$sex==varp),myvars ]
myvars <-c("pred")
y2 <- resids[which(resids$sex==varp),myvars ]
y2_sds<-SDs[which(SDs$name==c("TC")),c("stddev")]
nyrs<-length(y2_sds)
y2_up<-y2+2*(y2_sds)
y2_lo<-y2-2*(y2_sds)
myvars <-c("year")
xc <- resids[which(resids$sex==varp),myvars ]
#par(mar=c(4,4,4,4))
plot(xc, y2, type="l", lwd=2, xlab=" ", axes = TRUE, ylab=expression(Catch~x10^6), ylim=c(min(c(y1,y2))*.8,max(c(y1,y2))*1.2),  bty="l")#,xaxs = "i", yaxs = "i",pty="s") xlim=c(min(xc)-0.5,max(xc)+0.5), 
if (plotCI>0) lines(xc,y2_up,lty=2,lwd=2)
if (plotCI>0) lines(xc,y2_lo,lty=2,lwd=2)
points(xc,y1,pch=16)
text(2000,600,labels="Total Catch")

res=((y1-y2)-mean(y1-y2))/sd(y1-y2)
plot(xc, res, type="p", lwd=2, xlab=" ", axes = TRUE, ylab="Catch Residuals",   bty="l")#,xaxs = "i", yaxs = "i",pty="s") xlim=c(min(xc)-0.5,max(xc)+0.5), 
abline(h=0, lwd=2, lty=2)



#Observed and predicted recruitment indices

#par(mfrow=c(2,1), mar=c(4,4,2,2), oma=c(0,0,0,0), 	
#		mgp = c(2, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
#		tck = -0.01, 		#size of tick mark; neg is outside graph
#		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
#		cex.lab=1, 		#size of axis title
#		cex.axis=1,		#size of axis labels
#		bty="l")			#no box, just axes 
#plotCI=0 #turn off CI's

varp<-0
varp2 <- c("r")
myvars <-c("obs")
y1 <- resids[which(resids$sex==varp & resids$a_r==varp2 & resids$obs>0.000001& resids$obs<1e100),myvars ]
myvars <-c("pred")
y2 <- resids[which(resids$sex==varp & resids$a_r==varp2),myvars ]
nyrs<-length(y2)
y2_sds<-SDs[which(SDs$name==c("re_survey_est")),c("stddev")]
y2_up<-y2+2*(y2_sds)
y2_lo<-y2-2*(y2_sds)
myvars <-c("year")
x <- resids[which(resids$sex==varp & resids$a_r==varp2 & resids$obs>0.000001& resids$obs<1e100),myvars ]
miss=nyrs-length(x)+1
ytest=y2[miss+1:length(nyrs)]
#par(mar=c(4,4,4,4))
plot(x, y2[miss:nyrs], type="l", lwd=2, xlab=" ",  axes = TRUE, ylab="Juvenile IOA", ylim=c(min(c(y1,y2))*.8,max(c(y1,y2))*1.2), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
if (plotCI>0) lines(x,y2_up[miss:nyrs],lty=2,lwd=2)
if (plotCI>0) lines(x,y2_lo[miss:nyrs],lty=2,lwd=2)
points(x,y1,pch=16)

res=((y1-y2)-mean(y1-y2))/sd(y1-y2)
plot(xc, res, type="p", lwd=2, xlab=" ", axes = TRUE, ylab="Juvenile Residuals",   bty="l")#,xaxs = "i", yaxs = "i",pty="s") xlim=c(min(xc)-0.5,max(xc)+0.5), 
abline(h=0, lwd=2, lty=2)




#Observed and predicted adult indices

#par(mfrow=c(2,1), mar=c(4,4,2,2), oma=c(0,0,0,0), 	
#		mgp = c(2, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
#		tck = -0.01, 		#size of tick mark; neg is outside graph
#		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
#		cex.lab=1, 		#size of axis title
#		cex.axis=1,		#size of axis labels
#		bty="l")			#no box, just axes plotCI=0 #turn off CI's
#
varp<-0
varp2 <- c("a")
myvars <-c("obs")
y1 <- resids[which(resids$sex==varp & resids$a_r==varp2 & resids$obs>0.000001& resids$obs<1e100),myvars ]
myvars <-c("pred")
y2 <- resids[which(resids$sex==varp & resids$a_r==varp2),myvars ]
nyrs<-length(y2)
y2_sds<-SDs[which(SDs$name==c("ad_survey_est")),c("stddev")]
y2_up<-y2+2*(y2_sds)
y2_lo<-y2-2*(y2_sds)
myvars <-c("year")
x <- resids[which(resids$sex==varp & resids$a_r==varp2 & resids$obs>0.000001& resids$obs<1e100),myvars ]
miss=nyrs-length(x)+1
ytest=y2[miss+1:length(nyrs)]
#par(mar=c(4,4,4,4))
plot(x, y2[miss:nyrs], type="l", lwd=2, xlab=" ",  axes = TRUE, ylab="Adult IOA", ylim=c(min(c(y1,y2))*1.2,max(c(y1,y2))*1.2), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
if (plotCI>0) lines(x,y2_up[miss:nyrs],lty=2,lwd=2)
if (plotCI>0) lines(x,y2_lo[miss:nyrs],lty=2,lwd=2)
points(x,y1,pch=16)

res=((y1-y2)-mean(y1-y2))/sd(y1-y2)
plot(xc, res, type="p", lwd=2, xlab=" ", axes = TRUE, ylab="Adult Residuals",   bty="l")#,xaxs = "i", yaxs = "i",pty="s") xlim=c(min(xc)-0.5,max(xc)+0.5), 
abline(h=0, lwd=2, lty=2)




############################################
#Plot stock-recruitment
############################################
windows(record=T)
par(mfrow=c(2,1), mar=c(4,4,2,2), oma=c(0,0,0,0), 	
		mgp = c(2, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
		tck = -0.01, 		#size of tick mark; neg is outside graph
		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
		cex.lab=1, 		#size of axis title
		cex.axis=1,		#size of axis labels
		bty="l")			#no box, just axes plotCI=0 #turn off CI's


#Recruits at start of year
length=length(HPD$Spawners)
SPT<-HPD$Spawners[1:length-1]
R<-HPD$Rec[2:length]
plot(SPT,R, type="p", pch=16, xlab="Spawners (t)", xlim=c(0,max(SPT)+.1), ylab="Recruits (t+1)", ylim=c(0,max(R)+1), bty="l")#,xaxs = "i", yaxs = "i",pty="s")

varp<-c("alpha")
myvars <-c("Value")
alpha <- gen[which(gen$Name==varp),myvars ]

#bias correction
varp<-c("rec_cv")
myvars <-c("Value")
rec_cv <- gen[which(gen$Name==varp),myvars ]
#alpha=alpha*exp(-0.5*(log(rec_cv*rec_cv+1)))

varp<-c("beta")
myvars <-c("Value")
beta <- gen[which(gen$Name==varp),myvars ]

varp<-c("SRType")
myvars <-c("Value")
SRType <- gen[which(gen$Name==varp),myvars ]
spawners=seq(0,max(SPT),.01)
if (SRType==1) recruits=(spawners/(alpha+beta*spawners)) #*exp(-0.5*(log(rec_cv*rec_cv+1)))
if (SRType==2) recruits=(alpha*spawners*exp(-beta*spawners)) #*exp(-0.5*(log(rec_cv*rec_cv+1)))
lines(spawners,recruits,lty=2,lwd=2)



#recruit deviations

#Observed and predicted recruitment indices
varp<-c("r")
myvars <-c("obs")
y1 <- resids[which(resids$sex==varp ),myvars ]
y1 <- y1[2:length(y1)]
myvars <-c("pred")
y2 <- resids[which(resids$sex==varp ),myvars ]
y2 <- y2[2:length(y2)]
myvars <-c("year")
x <- resids[which(resids$sex==varp ),myvars ]
x <- x[2:length(x)]

res=((y1-y2)-mean(y1-y2))/sd(y1-y2)
plot(x, res, type="p", lty=2, lwd=2, xlab=" ", axes = TRUE, ylab="Residuals", bty="l")#,xaxs = "i", yaxs = "i",pty="s")
abline(h=0, lwd=2, lty=2)


############################################
#Catch per stage
############################################
par(mfrow=c(1,1), mar=c(4,4,2,2), oma=c(0,0,0,0), 	
		mgp = c(2, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
		tck = -0.01, 		#size of tick mark; neg is outside graph
		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
		cex.lab=1, 		#size of axis title
		cex.axis=1,		#size of axis labels
		bty="l")			#no box, just axes plotCI=0 #turn off CI's

#Numbers Caught 
x <- HPD$year
rec <- HPD$Rec*HPD$u0
ad <- HPD$Adult*HPD$u1
tot <- HPD$TC

varp<-c("MSY")
myvars <-c("Value")
MSYest<- gen[which(gen$Name==varp),myvars ]

plot(x, tot, type="l", lwd=2,xlab="Year", ylab=expression("Total Catch "~x10^6), xlim=c(min(x)-.01,max(x)),ylim=c(0,max(tot)*1.2), bty="l")
lines(x, rec, type="l", lwd=2, lty=2)
lines(x, ad, type="l", lwd=2, lty=3)
legend("topleft", c("Total","Juveniles", "Adults"), lwd=c(2,2,2), lty=c(1,2,3),cex=1, bty='n')
abline(h=MSYest)

text(2009,MSYest,labels=paste("MSY =",MSYest, sep=""), pos=3,cex=.8)


############################################
#Equilibrium catch at MSY
############################################
#Observed and predicted total catch
#select total catch data
par(mfrow=c(1,1), mar=c(4,4,2,2), oma=c(0,0,0,0), 	
		mgp = c(2, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
		tck = -0.01, 		#size of tick mark; neg is outside graph
		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
		cex.lab=1, 		#size of axis title
		cex.axis=1,		#size of axis labels
		bty="l")			#no box, just axes plotCI=0 #turn off CI's

x <- MSY$Fval
y <- MSY$C_eq

varp<-c("MSY")
myvars <-c("Value")
MSYest<- gen[which(gen$Name==varp),myvars ]
varp<-c("FMSY")
myvars <-c("Value")
FMSY<- gen[which(gen$Name==varp),myvars ]
varp<-c("uMSY")
myvars <-c("Value")
uMSY<- gen[which(gen$Name==varp),myvars ]


plot(x, y, type="l", lwd=2,xlab="Fishing Rate (F)", xlim=c(0,FMSY*1.5), axes = TRUE, ylab=expression("Equilibrium Catch "~x10^6), ylim=c(0,max(y)+5), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
abline(h=MSYest)
abline(v=FMSY)

text(.6,MSYest,labels=paste("MSY =",MSYest,sep=""), pos=1,cex=.8)
text(FMSY,0,labels=paste("FMSY =",FMSY), adj=c(-0.1,-1), srt=90,cex=.8)
#points(u0MSY,MSYest,cex=5)
#text(u0MSY,MSYest,adj=c(-.3,-1),labels=paste("F-msy =",FMSY),cex=.8)




############################################
#Ref Point Status
############################################
par(mfrow=c(1,1), mar=c(4,4,2,2), oma=c(0,0,0,0), 	
		mgp = c(2, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
		tck = -0.01, 		#size of tick mark; neg is outside graph
		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
		cex.lab=1, 		#size of axis title
		cex.axis=1,		#size of axis labels
		bty="l")			#no box, just axes plotCI=0 #turn off CI's


length=length(HPD$year)
years=HPD$year
FRat=HPD$FMSYRatio
NRat=HPD$NMSYRatio
varp<-c("FMSYRatio")
myvars <-c("Value")
Fest<- gen[which(gen$Name==varp),myvars ]

varp<-c("NMSYRatio")
myvars <-c("Value")
Nest<- gen[which(gen$Name==varp),myvars ]

varp<-c("cLim")
myvars <-c("Value")
cLim<- gen[which(gen$Name==varp),myvars ]

varp<-c("projYears")
myvars <-c("Value")
projYears<- gen[which(gen$Name==varp),myvars ]

xlim=seq(0,5,by=0.1)
ylim=rep(1,length(xlim))
for (i in 1:length(xlim)){
	if (xlim[i]<cLim) ylim[i]=xlim[i]/cLim
}
x=NRat
y=FRat
plot(x[1:(length(x)-1-projYears)], y[1:(length(y)-1-projYears)], type="p", pch=20, cex=.7,xlab="N/Nmsy", ylab="F/Fmsy", 
		xlim=c(0,max(x[1:length(x)-1])*1.2),ylim=c(0,max(c(1.2,y[1:length(y)-1]*1.2))), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
text(x[1:(length(x)-1-projYears)], y[1:(length(y)-1-projYears)], paste(years[1:(length(years)-1-projYears)]), cex=0.7, pos=4, offset =0.2)
lines(xlim, ylim, lty=2, lwd=2)
points(Nest,Fest,pch=8, cex=2)
#abline(h=1,v=1,lty=2, lwd=2)
legend("topleft", c("Default control rule", "Current status"), lwd=c(2,0), lty=c(2,0), pch=c(-1,8), bty='n')



##Do limit phase plot
x=NRat
y=FRat
plot(x[1:(length(x)-projYears)], y[1:(length(y)-projYears)], type="p", pch=20, cex=.7,xlab="N/Nmsy", ylab="F/Fmsy", 
		xlim=c(0,max(x[1:length(x)])*1.2),ylim=c(0,max(c(1.2,y[1:length(y)]*1.2))), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
polygon(c(-1,1,1,-1),c(-1,-1,1,1),col=rgb(1,1,.5,1), border=NA) #yellow
polygon(c(1,max(x[1:length(x)])*1.3,max(x[1:length(x)])*1.3,1),c(1,1,max(y[1:length(y)])*1.3,max(y[1:length(y)])*1.3),col=rgb(1,1,.5,1), border=NA) #yellow
polygon(c(1,max(x[1:length(x)])*1.3,max(x[1:length(x)])*1.3,1),c(-1,-1,1,1),col=rgb(.4,.75,.5,1), border=NA) #green
polygon(c(-1,1,1,-1),c(1,1,max(y[1:length(y)])*1.3,max(y[1:length(y)])*1.3),col=rgb(1,.4,.3,1), border=NA) #red
box()
points(x[1:(length(x)-projYears)], y[1:(length(y)-projYears)], pch=20, cex=.8)
text(x[1:(length(x)-projYears)], y[1:(length(y)-projYears)], paste(years[1:(length(years)-projYears)]), cex=0.7, pos=4, offset =0.2)
abline(v=1,h=1,lty=2,lwd=2)
points(Nest,Fest,pch=8, cex=2)


##Do limit phase plot
#varp<-c("FMSY")
#myvars <-c("Value")
#FMSY<- gen[which(gen$Name==varp),myvars ]
#
#varp<-c("NMSY")
#myvars <-c("Value")
#NMSY<- gen[which(gen$Name==varp),myvars ]
#
#F=	HPD$F[1:(length(HPD$F)-1-projYears)]
#N=HPD$Adult[1:(length(HPD$Adult)-1-projYears)]
#NLim=cLim*NMSY
#FLim=rep(FMSY, length(N))
#for (i in 1:length(N)){
#	if (N[i]<NLim) FLim[i]=(FMSY*N[i])/(cLim*NMSY)
#}
#
#plot(N/NLim, F/FMSY,type="p", pch=20, cex=.7,xlab="N/NLimit", ylab="F/FLimit", 
#		xlim=c(0,max(N/NMSY)*1.2),ylim=c(0,max(F/FLim)*1.2), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
#text(N/NLim, F/FLim, paste(years[1:(length(years)-1-projYears)]), cex=0.6, pos=4, offset =0.2)
#
#abline(h=1,v=1,lty=2, lwd=2)

############
## Plot F/FLim and N/NLim
############
windows(record=T) #,width=6.5, height=3.5)
par(mfrow=c(2,1), mar=c(4,4,2,2), oma=c(0,0,0,0), 	
		mgp = c(2, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
		tck = -0.01, 		#size of tick mark; neg is outside graph
		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
		cex.lab=1, 		#size of axis title
		cex.axis=1,		#size of axis labels
		bty="l")			#no box, just axes 

varp<-c("FMSY")
myvars <-c("Value")
FMSY<- gen[which(gen$Name==varp),myvars ]

varp<-c("NMSY")
myvars <-c("Value")
NMSY<- gen[which(gen$Name==varp),myvars ]

F=	HPD$F
N=HPD$Adult
NLim=cLim*NMSY
FLim=rep(FMSY, length(N))
for (i in 1:length(N)){
	if (N[i]<NLim) FLim[i]=(FMSY*N[i])/(cLim*NMSY)
}

Year=HPD$year
FFmsy=FRat
BBmsy=NRat
FFlim=F/FLim
BBlim=N/NLim

plot(Year, FFmsy, xlab="Year", ylab="F Ratio",type="l", lwd=2,bty="l", ylim=c(0, max(c(1.2,na.omit(FFmsy),na.omit(FFlim))*1.2)))
lines(Year, FFlim, lty=2, lwd=2 )
abline(h=1, lty=3,lwd=1)
legend("topright", c("F/FMSY", "F/FLimit"), lty=c(1,2), lwd=c(2,2), bty='n', cex=.9)

plot(Year, BBmsy, xlab="Year", ylab="Biomass Ratio",type="l", lwd=2,bty="l",  ylim=c(0, max(c(1.2,na.omit(BBmsy),na.omit(BBlim))*1.2)))
lines(Year, BBlim, lty=2, lwd=2 )
abline(h=1, lty=3,lwd=1)
legend("topright", c("B/BMSY", "B/BLimit"), lty=c(1,2), lwd=c(2,2), bty='n', cex=.9)



########################################
## Retrospective Analsyes -- not part of typical output, need to cut and paste it in
windows(record=T,width=5,height=7)
par(mfrow=c(2,1), mar=c(4,4,2,2), oma=c(0,0,0,0), 	
		mgp = c(2, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
		tck = -0.01, 		#size of tick mark; neg is outside graph
		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
		cex.lab=1, 		#size of axis title
		cex.axis=1,		#size of axis labels
		bty="l")			#no box, just axes plotCI=0 #turn off CI's

###########  East ###########  

east=read.table("Retro_east.dat", sep="\t", header=T)
x=east$Year
a0=east$Adult0
a1=east$Adult1
a2=east$Adult2
a3=east$Adult3
a4=east$Adult4
a5=east$Adult5

f0=east$F0
f1=east$F1
f2=east$F2
f3=east$F3
f4=east$F4
f5=east$F5


plot(x,a0, type="l", lty=1, xlab=" ", ylab=expression(Adults~x10^6) , ylim=c(min(na.omit(c(a0,a1,a2,a3,a4,a5)))*.8,max(na.omit(c(a0,a1,a2,a3,a4,a5)))*1.2), bty="l")
lines(x,a1, type="l", lty=2)
lines(x,a2, type="l", lty=3)
lines(x,a3, type="l", lty=4)
lines(x,a4, type="l", lty=5)
lines(x,a5, type="l", lty=6)
legend("topright", c("2011 (base)", "2010","2009","2008","2007","2006"), lty=c(1,2,3,4,5,6), bty='n', cex=.8)

plot(x[1:length(x)-1],f0[1:length(x)-1], type="l", lty=1, xlab=" ", ylab="Fishing rate (F)" , xlim=c(min(x),max(x)), ylim=c(min(na.omit(c(f0,f1,f2,f3,f4,f5)))*.8,max(na.omit(c(f0,f1,f2,f3,f4,f5)))*1.2), bty="l")
lines(x[1:(length(x)-2)],f1[1:(length(x)-2)], type="l", lty=2)
lines(x[1:(length(x)-3)],f2[1:(length(x)-3)], type="l", lty=3)
lines(x[1:(length(x)-4)],f3[1:(length(x)-4)], type="l", lty=4)
lines(x[1:(length(x)-5)],f4[1:(length(x)-5)], type="l", lty=5)
lines(x[1:(length(x)-6)],f5[1:(length(x)-6)], type="l", lty=6)
#legend("topright", c("2011", "2010","2009","2008","2007","2006"), lty=c(1,2,3,4,5,6), bty='n')


###########  West ###########  

west=read.table("Retro_west.dat", sep="\t", header=T)
x=west$Year
a0=west$Adult0
a1=west$Adult1
a2=west$Adult2
a3=west$Adult3
a4=west$Adult4
a5=west$Adult5

f0=west$F0
f1=west$F1
f2=west$F2
f3=west$F3
f4=west$F4
f5=west$F5


plot(x,a0, type="l", lty=1, xlab=" ", ylab=expression(Adults~x10^6) , ylim=c(min(na.omit(c(a0,a1,a2,a3,a4,a5)))*.8,max(na.omit(c(a0,a1,a2,a3,a4,a5)))*1.2), bty="l")
lines(x,a1, type="l", lty=2)
lines(x,a2, type="l", lty=3)
lines(x,a3, type="l", lty=4)
lines(x,a4, type="l", lty=5)
lines(x,a5, type="l", lty=6)
legend("topright", c("2011 (base)", "2010","2009","2008","2007","2006"), lty=c(1,2,3,4,5,6), bty='n', cex=.8)

plot(x[1:length(x)-1],f0[1:length(x)-1], type="l", lty=1, xlab=" ", ylab="Fishing rate (F)" , xlim=c(min(x),max(x)), ylim=c(min(na.omit(c(f0,f1,f2,f3,f4,f5)))*.8,max(na.omit(c(f0,f1,f2,f3,f4,f5)))*1.2), bty="l")
lines(x[1:(length(x)-2)],f1[1:(length(x)-2)], type="l", lty=2)
lines(x[1:(length(x)-3)],f2[1:(length(x)-3)], type="l", lty=3)
lines(x[1:(length(x)-4)],f3[1:(length(x)-4)], type="l", lty=4)
lines(x[1:(length(x)-5)],f4[1:(length(x)-5)], type="l", lty=5)
lines(x[1:(length(x)-6)],f5[1:(length(x)-6)], type="l", lty=6)
#legend("topright", c("2011", "2010","2009","2008","2007","2006"), lty=c(1,2,3,4,5,6), bty='n')


par(mfrow=c(2,1), mar=c(4,4,2,2), oma=c(0,0,0,0), 	
		mgp = c(2, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
		tck = -0.01, 		#size of tick mark; neg is outside graph
		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
		cex.lab=1, 		#size of axis title
		cex.axis=1,		#size of axis labels
		bty="l")			#no box, just axes plotCI=0 #turn off CI's



###########  Atlantic ###########  
atl=read.table("Retro_atl.dat", sep="\t", header=T)
x=atl$Year
a0=atl$Adult0
a1=atl$Adult1
a2=atl$Adult2
a3=atl$Adult3
a4=atl$Adult4
a5=atl$Adult5

f0=atl$F0
f1=atl$F1
f2=atl$F2
f3=atl$F3
f4=atl$F4
f5=atl$F5


plot(x,a0, type="l", lty=1, xlab=" ", ylab=expression(Adults~x10^6) , ylim=c(min(na.omit(c(a0,a1,a2,a3,a4,a5)))*.8,max(na.omit(c(a0,a1,a2,a3,a4,a5)))*1.2), bty="l")
lines(x,a1, type="l", lty=2)
lines(x,a2, type="l", lty=3)
lines(x,a3, type="l", lty=4)
lines(x,a4, type="l", lty=5)
lines(x,a5, type="l", lty=6)
legend("topright", c("2011 (base)", "2010","2009","2008","2007","2006"), lty=c(1,2,3,4,5,6), bty='n', cex=.8)

plot(x[1:length(x)-1],f0[1:length(x)-1], type="l", lty=1, xlab=" ", ylab="Fishing rate (F)" , xlim=c(min(x),max(x)), ylim=c(min(na.omit(c(f0,f1,f2,f3,f4,f5)))*.8,max(na.omit(c(f0,f1,f2,f3,f4,f5)))*1.2), bty="l")
lines(x[1:(length(x)-2)],f1[1:(length(x)-2)], type="l", lty=2)
lines(x[1:(length(x)-3)],f2[1:(length(x)-3)], type="l", lty=3)
lines(x[1:(length(x)-4)],f3[1:(length(x)-4)], type="l", lty=4)
lines(x[1:(length(x)-5)],f4[1:(length(x)-5)], type="l", lty=5)
lines(x[1:(length(x)-6)],f5[1:(length(x)-6)], type="l", lty=6)
#legend("topright", c("2011", "2010","2009","2008","2007","2006"), lty=c(1,2,3,4,5,6), bty='n')




###########  Gulf ###########  
glf=read.table("Retro_glf.dat", sep="\t", header=T)
x=glf$Year
a0=glf$Adult0
a1=glf$Adult1
a2=glf$Adult2
a3=glf$Adult3
a4=glf$Adult4
a5=glf$Adult5

f0=glf$F0
f1=glf$F1
f2=glf$F2
f3=glf$F3
f4=glf$F4
f5=glf$F5


plot(x,a0, type="l", lty=1, xlab=" ", ylab=expression(Adults~x10^6) , ylim=c(min(na.omit(c(a0,a1,a2,a3,a4,a5)))*.8,max(na.omit(c(a0,a1,a2,a3,a4,a5)))*1.2), bty="l")
lines(x,a1, type="l", lty=2)
lines(x,a2, type="l", lty=3)
lines(x,a3, type="l", lty=4)
lines(x,a4, type="l", lty=5)
lines(x,a5, type="l", lty=6)
legend("topright", c("2011 (base)", "2010","2009","2008","2007","2006"), lty=c(1,2,3,4,5,6), bty='n', cex=.8)

plot(x[1:length(x)-1],f0[1:length(x)-1], type="l", lty=1, xlab=" ", ylab="Fishing rate (F)" , xlim=c(min(x),max(x)), ylim=c(min(na.omit(c(f0,f1,f2,f3,f4,f5)))*.8,max(na.omit(c(f0,f1,f2,f3,f4,f5)))*1.2), bty="l")
lines(x[1:(length(x)-2)],f1[1:(length(x)-2)], type="l", lty=2)
lines(x[1:(length(x)-3)],f2[1:(length(x)-3)], type="l", lty=3)
lines(x[1:(length(x)-4)],f3[1:(length(x)-4)], type="l", lty=4)
lines(x[1:(length(x)-5)],f4[1:(length(x)-5)], type="l", lty=5)
lines(x[1:(length(x)-6)],f5[1:(length(x)-6)], type="l", lty=6)
#legend("topright", c("2011", "2010","2009","2008","2007","2006"), lty=c(1,2,3,4,5,6), bty='n')











##################################################
#Plot Environment
##################################################
windows(record=T, height=4, width=7)
par(mfrow=c(1,1), mar=c(4,4,2,2), oma=c(0,0,0,0), 	
		mgp = c(2, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
		tck = -0.01, 		#size of tick mark; neg is outside graph
		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
		cex.lab=1, 		#size of axis title
		cex.axis=1,		#size of axis labels
		bty="l")			#no box, just axes plotCI=0 #turn off CI's


#min(HPD$recM)
#max(HPD$recM)
plot(HPD$year, HPD$recM, type="l",lty=1,lwd=2, ylim=c(-.5,max(c(HPD$adM, HPD$recM))*1.2),xlab="Year", ylab="Mortality", bty="l")
lines(HPD$year, HPD$adM,lty=2,lwd=2)
lines(HPD$year, HPD$MEnvRec*.1,lty=3,lwd=2)
legend("topleft", c("Recruit M", "Adult M", "USGS Streamflow"), lty=c(1,2,3), lwd=c(2,2,2),bty='n', cex=.8)








##################################################
#Plot Composite for Sensitivity
##################################################

windows(record=T, width=7, height=8.5)
par(mfrow=c(4,1), mar=c(3,3,1,1), oma=c(0,0,0,0), 	
		mgp = c(2, 0.15, 0), 	#val 1=axis title location, 2=axis labels, 3=axis line
		tck = -0.01, 		#size of tick mark; neg is outside graph
		font.lab=1, 		#font 1=normal, 2=bold, 3=italic
		cex.lab=1, 		#size of axis title
		cex.axis=1,		#size of axis labels
		bty="l")			#no box, just axes 
plotCI=0 #turn off CI's

#select total catch data
varp<-c("t")
myvars <-c("obs")
y1 <- resids[which(resids$sex==varp),myvars ]
myvars <-c("pred")
y2 <- resids[which(resids$sex==varp),myvars ]
y2_sds<-SDs[which(SDs$name==c("TC")),c("stddev")]
nyrs<-length(y2_sds)
y2_up<-y2+2*(y2_sds)
y2_lo<-y2-2*(y2_sds)
myvars <-c("year")
xc <- resids[which(resids$sex==varp),myvars ]
#par(mar=c(4,4,4,4))
plot(xc, y2, type="l", lwd=2, xlab=" ", axes = TRUE, ylab=expression(Catch~x10^6), ylim=c(min(c(y1,y2))*.8,max(c(y1,y2))*1.2),  bty="l")#,xaxs = "i", yaxs = "i",pty="s") xlim=c(min(xc)-0.5,max(xc)+0.5), 
if (plotCI>0) lines(xc,y2_up,lty=2,lwd=2)
if (plotCI>0) lines(xc,y2_lo,lty=2,lwd=2)
points(xc,y1,pch=16)
text(2000,600,labels="Total Catch")


varp<-0
varp2 <- c("r")
myvars <-c("obs")
y1 <- resids[which(resids$sex==varp & resids$a_r==varp2 & resids$obs>0.000001& resids$obs<1e100),myvars ]
myvars <-c("pred")
y2 <- resids[which(resids$sex==varp & resids$a_r==varp2),myvars ]
nyrs<-length(y2)
y2_sds<-SDs[which(SDs$name==c("re_survey_est")),c("stddev")]
y2_up<-y2+2*(y2_sds)
y2_lo<-y2-2*(y2_sds)
myvars <-c("year")
x <- resids[which(resids$sex==varp & resids$a_r==varp2 & resids$obs>0.000001& resids$obs<1e100),myvars ]
miss=nyrs-length(x)+1
ytest=y2[miss+1:length(nyrs)]
#par(mar=c(4,4,4,4))
plot(x, y2[miss:nyrs], type="l", lwd=2, xlab=" ",  axes = TRUE, ylab="Juvenile IOA", ylim=c(min(c(y1,y2))*.8,max(c(y1,y2))*1.2), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
if (plotCI>0) lines(x,y2_up[miss:nyrs],lty=2,lwd=2)
if (plotCI>0) lines(x,y2_lo[miss:nyrs],lty=2,lwd=2)
points(x,y1,pch=16)


varp<-0
varp2 <- c("a")
myvars <-c("obs")
y1 <- resids[which(resids$sex==varp & resids$a_r==varp2 & resids$obs>0.000001& resids$obs<1e100),myvars ]
myvars <-c("pred")
y2 <- resids[which(resids$sex==varp & resids$a_r==varp2),myvars ]
nyrs<-length(y2)
y2_sds<-SDs[which(SDs$name==c("ad_survey_est")),c("stddev")]
y2_up<-y2+2*(y2_sds)
y2_lo<-y2-2*(y2_sds)
myvars <-c("year")
x <- resids[which(resids$sex==varp & resids$a_r==varp2 & resids$obs>0.000001& resids$obs<1e100),myvars ]
miss=nyrs-length(x)+1
ytest=y2[miss+1:length(nyrs)]
#par(mar=c(4,4,4,4))
plot(x, y2[miss:nyrs], type="l", lwd=2, xlab=" ",  axes = TRUE, ylab="Adult IOA", ylim=c(min(c(y1,y2))*1.2,max(c(y1,y2))*1.2), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
if (plotCI>0) lines(x,y2_up[miss:nyrs],lty=2,lwd=2)
if (plotCI>0) lines(x,y2_lo[miss:nyrs],lty=2,lwd=2)
points(x,y1,pch=16)

#min(HPD$recM)
#max(HPD$recM)
plot(HPD$year, HPD$recM, type="l",lty=1,lwd=2, ylim=c(-.5,max(c(HPD$adM, HPD$recM))*1.2),xlab="Year", ylab="Natural Mortality", bty="l")
lines(HPD$year, HPD$adM,lty=2,lwd=2)
lines(HPD$year, HPD$MEnvRec*.2,lty=3,lwd=2)
legend("topleft", c("Recruit M", "Adult M", "USGS Streamflow"), lty=c(1,2,3), lwd=c(2,2,2),bty='n', cex=1)




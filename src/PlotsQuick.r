
rm(list = ls())

pdfFile="Graphs.pdf"
pdfOn=0  #0=off

twoPlots=0 #0=off

if (pdfOn==0) windows(record=T, width=4.25, height=5.75)
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




par(mfrow=c(3,1))#,mar=c(5,5,2,5),oma=c(0,0,0,0),mgp=c(2.2,1,0))
par(mfrow=c(3,1), mar=c(5.1,4.1,4.1,2.1), oma=c(0,0,0,0)) #default margins
plotCI=0 #turn off CI's


#############################################
##Observed and predicted total catch
#############################################
##select total catch data
#varp<-c("t")
#myvars <-c("obs")
#y1 <- resids[which(resids$sex==varp),myvars ]
#myvars <-c("pred")
#y2 <- resids[which(resids$sex==varp),myvars ]
#y2_sds<-SDs[which(SDs$name==c("TC")),c("stddev")]
#nyrs<-length(y2_sds)
#y2_up<-y2+2*(y2_sds)
#y2_lo<-y2-2*(y2_sds)
#myvars <-c("year")
#xc <- resids[which(resids$sex==varp),myvars ]
##par(mar=c(4,4,4,4))
#plot(xc, y2, type="l", xlab=" ", axes = TRUE, ylab=expression(Catch~x10^6), ylim=c(min(c(y1,y2))*.8,max(c(y1,y2))*1.2),  bty="l")#,xaxs = "i", yaxs = "i",pty="s") xlim=c(min(xc)-0.5,max(xc)+0.5), 
##lines(xc,y2,lty=1, lwd=2)
#if (plotCI>0) lines(xc,y2_up,lty=2,lwd=2)
#if (plotCI>0) lines(xc,y2_lo,lty=2,lwd=2)
#points(xc,y1,pch=16)
#text(2000,600,labels="Total Catch")
#
#
#
##Observed and predicted recruitment indices
#varp<-0
#varp2 <- c("r")
#myvars <-c("obs")
#y1 <- resids[which(resids$sex==varp & resids$a_r==varp2 & resids$obs>0.000001& resids$obs<1e100),myvars ]
#myvars <-c("pred")
#y2 <- resids[which(resids$sex==varp & resids$a_r==varp2),myvars ]
#nyrs<-length(y2)
#y2_sds<-SDs[which(SDs$name==c("re_survey_est")),c("stddev")]
#y2_up<-y2+2*(y2_sds)
#y2_lo<-y2-2*(y2_sds)
#myvars <-c("year")
#x <- resids[which(resids$sex==varp & resids$a_r==varp2 & resids$obs>0.000001& resids$obs<1e100),myvars ]
#miss=nyrs-length(x)+1
#ytest=y2[miss+1:length(nyrs)]
##par(mar=c(4,4,4,4))
#plot(x, y2[miss:nyrs], type="l", xlab=" ", xlim=c(min(xc)-.5,max(xc)+.5), axes = TRUE, ylab="Age-0 CPUE", ylim=c(0,max(y1)+.5), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
#lines(x,y2[miss:nyrs],lty=1, lwd=2)
#if (plotCI>0) lines(x,y2_up[miss:nyrs],lty=2,lwd=2)
#if (plotCI>0) lines(x,y2_lo[miss:nyrs],lty=2,lwd=2)
#points(x,y1,pch=16)
#
#recIOA=data.frame(year=x,obs=y1,pred=y2)
##Observed and predicted adult indices
#varp<-0
#varp2 <- c("a")
#myvars <-c("obs")
#y1 <- resids[which(resids$sex==varp & resids$a_r==varp2 & resids$obs>0.000001& resids$obs<1e100),myvars ]
#myvars <-c("pred")
#y2 <- resids[which(resids$sex==varp & resids$a_r==varp2),myvars ]
#nyrs<-length(y2)
#y2_sds<-SDs[which(SDs$name==c("ad_survey_est")),c("stddev")]
#y2_up<-y2+2*(y2_sds)
#y2_lo<-y2-2*(y2_sds)
#myvars <-c("year")
#x <- resids[which(resids$sex==varp & resids$a_r==varp2 & resids$obs>0.000001& resids$obs<1e100),myvars ]
#miss=nyrs-length(x)+1
#ytest=y2[miss+1:length(nyrs)]
##par(mar=c(4,4,4,4))
#plot(x, y2[miss:nyrs], type="l", xlab=" ", xlim=c(min(xc)-.5,max(xc)+.5), axes = TRUE, ylab="Age-1+ CPUE", ylim=c(0,max(y1)+.5), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
#lines(x,y2[miss:nyrs],lty=1, lwd=2)
#if (plotCI>0) lines(x,y2_up[miss:nyrs],lty=2,lwd=2)
#if (plotCI>0) lines(x,y2_lo[miss:nyrs],lty=2,lwd=2)
#points(x,y1,pch=16)


############################################
#Observed and predicted total catch
############################################
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
if (twoPlots==0) {
	
	par(mfrow=c(2,1))
#par(mar=c(5.1,4.1,4.1,2.1))
	
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
#Ref Points Rates relative to MSY
	############################################
	par(mfrow=c(3,1))
#Fmsy
	x <- HPD$year
	y <- HPD$F
	y_sds<-SDs[which(SDs$name==c("F")),c("stddev")]
	y_up<-y+2*(y_sds)
	y_lo<-y-2*(y_sds)
	varp<-c("FMSY")
	myvars <-c("Value")
	FMSY<- gen[which(gen$Name==varp),myvars ]
	
	plot(x[1:length(x)-1], y[1:length(y)-1], type="l", lwd=2,xlab="Year", ylab="Fishing Mortality", xlim=c(min(x)-.01,max(x)),ylim=c(min(y)-.2,max(y)+.2), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
#lines(x, y_up, lty=2, lwd=1.5)
#lines(x, y_lo, lty=2, lwd=1.5)
	lines(x[1:length(x)-1], y[1:length(y)-1], type="l", lwd=2)
	abline(h=FMSY,lty=2, lwd=2)
	text(2011,FMSY,labels="Fmsy", pos=1,cex=.8)
	
#Nmsy
	length=length(HPD$year)
	x <- HPD$year
	y <- HPD$Adult
	y_sds<-SDs[which(SDs$name==c("N")),c("stddev")]
	y_up<-y+2*(y_sds[1:length])
	y_lo<-y-2*(y_sds[1:length])
	varp<-c("NMSY")
	myvars <-c("Value")
	NMSY<- gen[which(gen$Name==varp),myvars ]
	
	plot(x[1:length(x)-1], y[1:length(y)-1], type="l", lwd=2,xlab="Year", ylab=expression("Number of Adult Crabs "~x10^6), xlim=c(min(x)-.01,max(x)),ylim=c(0,max(y)+2), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
#lines(x, y_up, lty=2, lwd=1.5)
#lines(x, y_lo, lty=2, lwd=1.5)
	lines(x[1:length(x)-1], y[1:length(y)-1], type="l", lwd=2)
	abline(h=NMSY,lty=2, lwd=2)
	text(2011,NMSY,labels="Nmsy", pos=1,cex=.8)
	
#Rmsy
	length=length(HPD$year)
	x <- HPD$year
	y <- HPD$Rec
	y_sds<-SDs[which(SDs$name==c("N")),c("stddev")]
	y_up<-y+2*(y_sds[1:length])
	y_lo<-y-2*(y_sds[1:length])
	varp<-c("RMSY")
	myvars <-c("Value")
	RMSY<- gen[which(gen$Name==varp),myvars ]
	
	plot(x[1:length(x)-1], y[1:length(y)-1], type="l", lwd=2,xlab="Year", ylab=expression("Number of Recruit Crabs "~x10^6), xlim=c(min(x)-.01,max(x)),ylim=c(0,max(y)+2), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
#lines(x, y_up, lty=2, lwd=1.5)
#lines(x, y_lo, lty=2, lwd=1.5)
	lines(x[1:length(x)-1], y[1:length(y)-1], type="l", lwd=2)
	abline(h=RMSY,lty=2, lwd=2)
	text(2011,RMSY,labels="Rmsy", pos=1,cex=.8)
	
	
	############################################
#Exploitation Rates relative to MSY
	############################################
#u0 
	par(mfrow=c(3,1))
	
	length=length(HPD$year)
	x <- HPD$year[1:length]
	y <- HPD$u0[1:length]
	varp<-c("u0MSY")
	myvars <-c("Value")
	uMSY<- gen[which(gen$Name==varp),myvars ]
	
	plot(x[1:length(x)-1], y[1:length(y)-1], type="l", lwd=2,xlab="Year", ylab="Recruit Harvest Rt (u0)", xlim=c(min(x)-.01,max(x)),ylim=c(min(y)-.05,max(y)+.05), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
	lines(x[1:length(x)-1], y[1:length(y)-1], type="l", lwd=2)
	abline(h=uMSY,lty=2, lwd=2)
	text(2011,uMSY,labels="u0 msy", pos=3,cex=1)
	
#u1
	length=length(HPD$year)
	x <- HPD$year[1:length]
	y <- HPD$u1[1:length]
	varp<-c("u1MSY")
	myvars <-c("Value")
	uMSY<- gen[which(gen$Name==varp),myvars ]
	
	plot(x[1:length(x)-1], y[1:length(y)-1], type="l", lwd=2,xlab="Year", ylab="Adult Harvest Rt (u1+)", xlim=c(min(x)-.01,max(x)),ylim=c(min(y)-.05,max(y)+.05), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
	lines(x[1:length(x)-1], y[1:length(y)-1], type="l", lwd=2)
	abline(h=uMSY,lty=2, lwd=2)
	text(2011,uMSY,labels="u1+ msy", pos=3,cex=1)
	
#uAll
	length=length(HPD$year)
	x <- HPD$year[1:length]
	y <- HPD$uAll[1:length]
	y_sds<-SDs[which(SDs$name==c("u")),c("stddev")]
	y_up<-y+2*(y_sds)
	y_lo<-y-2*(y_sds)
	varp<-c("uMSY")
	myvars <-c("Value")
	uMSY<- gen[which(gen$Name==varp),myvars ]
	
	plot(x[1:length(x)-1], y[1:length(y)-1], type="l", lwd=2,xlab="Year", ylab="Total Harvest Rt (u0+)", xlim=c(min(x)-.01,max(x)),ylim=c(min(y)-.05,max(y)+.05), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
#lines(x, y_up, lty=2, lwd=1.5)
#lines(x, y_lo, lty=2, lwd=1.5)
	lines(x[1:length(x)-1], y[1:length(y)-1], type="l", lwd=2)
	abline(h=uMSY,lty=2, lwd=2)
	text(2011,uMSY,labels="u0+ msy", pos=3,cex=1)
	
	
	
	############################################
#Equilibrium catch at MSY and Numbers Catch
	############################################
#Observed and predicted total catch
#select total catch data
	par(mfrow=c(2,1))
	x <- MSY$uAll_eq
	y <- MSY$C_eq
	plot(x, y, type="l", lwd=2,xlab="Total Harvest Rate (u0+)", xlim=c(0,.7), axes = TRUE, ylab=expression("Equilibrium Catch "~x10^6), ylim=c(0,max(y)+5), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
	varp<-c("MSY")
	myvars <-c("Value")
	MSYest<- gen[which(gen$Name==varp),myvars ]
	varp<-c("FMSY")
	myvars <-c("Value")
	FMSY<- gen[which(gen$Name==varp),myvars ]
	varp<-c("uMSY")
	myvars <-c("Value")
	uMSY<- gen[which(gen$Name==varp),myvars ]
	abline(h=MSYest)
	abline(v=uMSY)
	
	text(.1,MSYest,labels=paste("msy=",MSYest,", Fmsy=",FMSY, sep=""), pos=1,cex=.8)
	text(uMSY,0,labels=paste("u-msy =",uMSY), adj=c(-0.1,-1), srt=90,cex=.8)
#points(u0MSY,MSYest,cex=5)
#text(u0MSY,MSYest,adj=c(-.3,-1),labels=paste("F-msy =",FMSY),cex=.8)
	
#Numbers Caught 
	x <- HPD$year
	rec <- HPD$Rec*HPD$u0
	ad <- HPD$Adult*HPD$u1
	tot <- HPD$TC
	
	
	plot(x, tot, type="l", lwd=2,xlab="Year", ylab=expression("Total Catch "~x10^6), xlim=c(min(x)-.01,max(x)),ylim=c(0,max(tot)+.1), bty="l")
	lines(x, rec, type="l", lwd=2, lty=2)
	lines(x, ad, type="l", lwd=2, lty=3)
	legend("topleft", c("age-0+","age-0", "age-1+"), lwd=c(2,2,2), lty=c(1,2,3),cex=.6, bty='n')
	abline(h=MSYest)
	
	text(2009,MSYest,labels=paste("msy=",MSYest, sep=""), pos=3,cex=.8)
}

############################################
#Ref Point Status
############################################
par(mfrow=c(1,1))#, mar=c(5.1,4.5,4.1,2.1))


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
		xlim=c(0,max(x[1:length(x)-1])*1.2),ylim=c(0,max(y[1:length(y)-1])*1.2), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
text(x[1:(length(x)-1-projYears)], y[1:(length(y)-1-projYears)], paste(years[1:(length(years)-1-projYears)]), cex=0.5, pos=4, offset =0.2)
lines(xlim, ylim, lty=2, lwd=2)
points(Nest,Fest,pch=8, cex=2)
#abline(h=1,v=1,lty=2, lwd=2)
legend("topright",cex=.7, c("Default control rule", "Current status"), lwd=c(2,0), lty=c(2,0), pch=c(-1,8)) #bty='n',


#Do limit phase plot
##Do limit phase plot
x=NRat
y=FRat
plot(x[1:(length(x)-projYears)], y[1:(length(y)-projYears)], type="p", pch=20, cex=.7,xlab="N/Nmsy", ylab="F/Fmsy", 
		xlim=c(0,max(c(x[1:length(x)],1.2))*1.2),ylim=c(0,max(c(1.2,y[1:length(y)]*1.2))), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
polygon(c(-1,1,1,-1),c(-1,-1,1,1),col=rgb(1,1,.5,1), border=NA) #yellow
polygon(c(1,max(c(x[1:length(x)],1.2))*1.3,max(c(x[1:length(x)],1.2))*1.3,1),c(1,1,max(c(y[1:length(y)],1.2))*1.3,max(c(y[1:length(y)],1.2))*1.3),col=rgb(1,1,.5,1), border=NA) #yellow
polygon(c(1,max(c(x[1:length(x)],1.2))*1.3,max(c(x[1:length(x)],1.2))*1.3,1),c(-1,-1,1,1),col=rgb(.4,.75,.5,1), border=NA) #green
polygon(c(-1,1,1,-1),c(1,1,max(c(y[1:length(y)],1.2))*1.3,max(c(y[1:length(y)],1.2))*1.3),col=rgb(1,.4,.3,1), border=NA) #red
box()
points(x[1:(length(x)-projYears)], y[1:(length(y)-projYears)], pch=20, cex=.8)
text(x[1:(length(x)-projYears)], y[1:(length(y)-projYears)], paste(years[1:(length(years)-projYears)]), cex=0.5, pos=4, offset =0.2)
abline(v=1,h=1,lty=2,lwd=2)
points(Nest,Fest,pch=8, cex=2)

#Do limit phase plot
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
#plot(N/NLim, F/FLim,type="p", pch=20, cex=.7,xlab="N/NLimit", ylab="F/FLimit", 
#		xlim=c(0,max(N/NLim)*1.2),ylim=c(0,max(F/FLim)*1.2), bty="l")#,xaxs = "i", yaxs = "i",pty="s")
#text(N/NLim, F/FLim, paste(years[1:(length(years)-1-projYears)]), cex=0.6, pos=4, offset =0.2)
#
#abline(h=1,v=1,lty=2, lwd=2)

############################################
#SPR
############################################
if (twoPlots==0){
	
	
	x=tSPR$year
	y=tSPR$tSPR
	plot(x[1:length(x)-1], y[1:length(y)-1],type="l", xlim=c(min(x),max(x)),ylim=c(0, max(y)*1.1),lwd=2, xlab="Year", ylab="tSPR", bty="l")
	
}



##################################################
#Plot Environment
##################################################
par(mfrow=c(1,1))#, mar=c(5.1,4.5,4.1,2.1))


#min(HPD$recM)
#max(HPD$recM)
plot(HPD$year, HPD$recM, type="l",lty=1,lwd=2, ylim=c(-.3,3.5),xlab="Year", ylab="Mortality", bty="l")
lines(HPD$year, HPD$adM,lty=2,lwd=2)
lines(HPD$year, HPD$MEnvRec*.1,lty=3,lwd=2)
lines(HPD$year, HPD$MEnvAd*.1-.1,lty=4,lwd=2)
legend("topleft", c("Recruit M", "Adult M", "Env Index Rec", "Env Index Ad"), lty=c(1,2,3,4), lwd=c(2,2,2,2))



##################################################
#Plot Environment 2
##################################################
par(mfrow=c(1,1))#, mar=c(5.1,4.5,4.1,2.1))


#min(HPD$recM)
#max(HPD$recM)
plot(HPD$year, HPD$recM, type="l",lty=1,lwd=2, ylim=c(-.3,5),xlab="Year", ylab="Mortality", bty="l")
lines(HPD$year, HPD$adM,lty=2,lwd=2)
lines(HPD$year, HPD$MEnvRec*.1,lty=3,lwd=2)
legend("topleft", c("Recruit M", "Adult M", "Env Index"), lty=c(1,2,3), lwd=c(2,2,2), cex=.8)

##################################################
#Plot Likelihood profile for MSY
##################################################


getLikeProf <- function(filename){
	
	file=basename(filename)
	
	ifile=scan(filename,what="character",flush=T,blank.lines.skip=F, quiet=T) 
	idx=sapply(as.double(ifile),is.na) 
	indexes=list() 
	counter = 1 
	for (i in 1:length(idx)){ 
		if( (i!=length(idx)) & (idx[i]==T) & (idx[i+1]==F)) { 
			j=i+1 
			size=0 
			while ( (j<=length(idx)) & (idx[j] ==F) ) { 
				size = size+1 
				j = j +1 
			} 
			indexes[[counter ]] = c(i,size) 
			counter = counter + 1 
		} 
	} 
	newlist=list()
	for (i in 1:length(indexes)){ 
		vals <- indexes[[i]] 
		start <- vals[1] 
		rows <- vals[2] 
		if (rows==1) dum=as.double(scan(filename,skip=start,nlines=1,quiet=T,what="")) 
		else dum=as.matrix(read.table(filename,skip=start,nrow=rows,fill=T)) 
		if(is.numeric(dum))#Logical test to ensure dealing with numbers 
		{ 
			line = as.character(scan(filename,skip=start-1,nlines=1,quiet=T,what=""))
			newline=paste(line,collapse="_")
			newline=gsub("#_","",newline)
			newline=gsub("//.","",newline)
			newline=gsub("#","",newline)
			newline=gsub("//","",newline)
			newlist[[newline]]=dum
			
		} 
	}
	assign(file, newlist, envir = .GlobalEnv)
}

#getLikeProf("steep.plt")
# doesn't work  --  getLP("uMSY.plt")

#plot(steep.plt$Profile_likelihood[,1],steep.plt$Profile_likelihood[,2], type="l", lwd=2, 
#		xlab="steep Estimate", ylab="Likelihood", bty="l")

if (pdfOn==1) dev.off()
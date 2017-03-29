###############################################################################################################################################
                                             ########CP - Spectrum (mostly graphics) analyses ########
###############################################################################################################################################
rm=list(ls())
graphics.off()
library("zoo")
library("lubridate")
source("functions_global.r")

#Phytoplankton data 
sp=c("GUI", "LEP", "NIT", "PSE", "RHI", "SKE", "CHA", "AST", "GYM",  "PRP", "EUG", "CRY")
sp=sort(sp)

#Loading data from Teychan and B7
filename=paste("Teychan_base.csv",sep="")
tabbis=read.csv(filename,na.strings="NA",header=TRUE,sep=";",dec=".")
dates=as.Date(tabbis$Date)

filename2=paste("B7_base.csv",sep="")
tabbis2=read.csv(filename2,na.strings="NA",header=TRUE,sep=";",dec=".")
dates2=as.Date(tabbis2$Date)
dates_bis2=seq(dates2[1],dates2[length(dates2)],14)

#First analysis of phytoplanktonic spectra that may differ according to random values used in the reconstruction of the time series
pdf("spectrum_flore_bis_for_each_species_clean.pdf")
par(mfrow=c(3,2))
xli1=8
xli2=72
maxi_test=10

for (s in sp){
#Removing dates before 1997 for species that were not well-recognized
if(s=="CRY"|s=="EUG"|s=="GYM"){
dates_correct=dates[year(dates)>1996]
var=tabbis[year(dates)>1996,s]
}else{
dates_correct=dates
var=tabbis[,s]}
dates_bis=seq(dates_correct[1],dates_correct[length(dates_correct)],14)

#First time series replaces missing values with a fixed constant corresponding to the smallest number of planktonic cells that can be detected in a sample
var1=na.approx(var,maxgap=2,x=dates_correct,xout=dates_bis,na.rm=FALSE)
var1_bis=matrix(NA,nrow=maxi_test,ncol=length(var1))
var1_bis[1,]=var1
var1_bis[1,is.na(var1)]=10

var2=na.approx(tabbis2[,s],maxgap=2,x=dates2,xout=dates_bis2,na.rm=FALSE)
var2_bis=matrix(NA,nrow=maxi_test,ncol=length(var2))
var2_bis[1,]=var2
var2_bis[1,is.na(var2)]=10

#Spectrum computation
x=spectrum(var1_bis[1,],plot=FALSE,method="ar")
x2=spectrum(var2_bis[1,],plot=FALSE,method="ar")
xlog=spectrum(log10(var1_bis[1,]),plot=FALSE,method="ar")
x2log=spectrum(log10(var2_bis[1,]),plot=FALSE,method="ar")

#Computing periods in months to make the graph easier to understand
freq=1/x$freq*14/(365.25/12)
freq2=1/x2$freq*14/(365.25/12)

#Then, compute different spectra with different time series interpolation: a random values is drawn between 0 and the minimum value of the time series for each missing value
#Initialize matrices
x_bis=matrix(NA,nrow=maxi_test,ncol=length(x$spec))
x2_bis=matrix(NA,nrow=maxi_test,ncol=length(x2$spec))
x_bis[1,]=x$spec
x2_bis[1,]=x2$spec

xlog_bis=matrix(NA,nrow=maxi_test,ncol=length(xlog$spec))
x2log_bis=matrix(NA,nrow=maxi_test,ncol=length(x2log$spec))
xlog_bis[1,]=xlog$spec
x2log_bis[1,]=x2log$spec

#Compute spectra
for (i in 2:maxi_test){
var1_bis[i,]=var1
var1_bis[i,is.na(var1)]=runif(sum(is.na(var1)),0,min(var1,na.rm=TRUE))
var2_bis[i,]=var2
var2_bis[i,is.na(var2)]=runif(sum(is.na(var2)),0,min(var2,na.rm=TRUE))

x_bis[i,]=spectrum(var1_bis[i,],plot="FALSE",method="ar")$spec
x2_bis[i,]=spectrum(var2_bis[i,],plot="FALSE",method="ar")$spec
xlog_bis[i,]=spectrum(log10(var1_bis[i,]),plot="FALSE",method="ar")$spec
x2log_bis[i,]=spectrum(log10(var2_bis[i,]),plot="FALSE",method="ar")$spec
}

#Plot everything
yli1=min(xlog_bis)
yli2=max(xlog_bis)
plot(freq,xlog_bis[1,],"l",xlim=c(xli1,xli2),xaxt="n",xlab="Periode (months)",main=paste("log10",s),lty=1,ylim=c(yli1,yli2),cex=1.3,ylab="Spectrum",lwd=3)
for(i in 2:maxi_test){
lines(freq,xlog_bis[i,],lty=3,col="grey",cex=0.1)
}
par(new=TRUE)
yli1=min(x2log_bis)
yli2=max(x2log_bis)
plot(freq2,x2log_bis[1,],xlim=c(xli1,xli2),"l",col="blue",xaxt="n",yaxt="n",ylab="",xlab="",ylim=c(yli1,yli2),cex=1.3,lwd=3)
for(i in 2:maxi_test){
lines(freq2,x2log_bis[i,],lty=3,col="cyan",cex=0.1)
}
legend("topright",c("Teychan, fixed","Teychan, random","B7, fixed","B7, random"),col=c("black","grey","blue","cyan"),lty=c(1,2),bty="n")
labe=seq(6,12*xli2,6)
axis(1,at=labe,labels=labe)
axis(4,col="blue")
}
dev.off()

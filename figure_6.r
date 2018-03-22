###############################################################################################################################################
                                      ########CP - Plot Fig. 6 in Barraquand et al. 2017 using MARSS results ########
###############################################################################################################################################

### Initialize
rm(list=ls())
graphics.off()

library("MARSS")
library("zoo")
library("lubridate")

###Loading data
#Simulation
openfile=paste("Teychan_physics_pencen.RData")
load(openfile)

#Real data
lieu="Teychan"
sp=c("AST","NIT","PSE","SKE","CHA","GUI","LEP","RHI","GYM","PRP","CRY","EUG")
filename=paste(lieu,"_base.csv",sep="")
tabbis=read.csv(filename,na.strings="NA",header=TRUE,sep=";",dec=".")
dates=as.Date(tabbis$Date)
tab=tabbis[year(dates)>1996,]#Using data from 1997, because cryptophytes were not counted (or badly, for the first year, before)
dates=dates[year(dates)>1996]
tab_sp=tab[,sp]

###Treating missing values###
consecutif=2 #Number of missing values above which we keep the NA
timestep=14 #Regular time lapses between two observations
dates_bis=seq(dates[1],dates[length(dates)],timestep) #Regular time grid
tab_sp=na.approx(tab_sp,maxgap=consecutif,x=dates,xout=dates_bis,na.rm=FALSE) #Interpolation over regular time grid
tab_sp=log(tab_sp)
tab_sp=t(scale(tab_sp[2:(length(dates_bis)),]))
rownames(tab_sp)=sp
dates_bis=dates_bis[2:length(dates_bis)]

#Initialize graph frame
pdf("figure_6.pdf",width=20,height=20)

#Graphical parameters
par(mfrow=c(4,3),mar=c(2,2,3,1),oma=c(2,6.5,1,0.5))
fac_main=4.0
fac_axis=3.5
fac_lg=3.0
fac_lab=3.5
alwd=3.5
apc=2

#Draw the graph
tt=c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)")
for (i in 1:length(sp)){
	plot(dates_bis,fit_log$states[i,],t="l",col="black",cex.main=fac_main,cex.axis=fac_axis,cex.lab=fac_lab,lwd=alwd,cex=apc,ps=apc,main=sp[i],ylab="",xlab="",xlim=c(min(dates_bis),max(dates_bis)+50),ylim=c(-3.5,3.5),xaxt="n",yaxt="n")
	points(dates_bis,tab_sp[sp[i],],t="p",col="red",cex.main=fac_main,cex.axis=fac_axis,cex.lab=fac_lab,lwd=alwd,cex=apc,ps=apc,pch=16)
	text(min(dates_bis)+100,3.15,tt[i],cex=fac_axis)
	#x-axis
	if(i>9){
		axis(1,at=as.Date(c("2000/01/01","2005/01/01","2010/01/01","2015/01/01")),labels=rep("",4))
		mtext(c('2000','2005','2010','2015'),side=1,line=2.5,at=as.Date(c("2000/01/01","2005/01/01","2010/01/01","2015/01/01")),cex=2.)
	}
	#y-axis
	if(i%%3==1){
		axis(2,at=seq(-3,3,1),labels=c("-3","","-1","","1","","3"),cex.axis=fac_axis)
	}
	if(i==4){
		mtext("Scaled log abundance",side=2,at=-5,line=4.5,cex=fac_axis)
		legend("topright",c("Observation","Estimation"),col=c("red","black"),lty=c(NA,1),pch=c(16,NA),cex=fac_lg,lwd=alwd,pt.cex=apc,bty="n")
	}
}
dev.off()


rm(list=ls())
graphics.off()

###########Figure S5.1
path_data_post="/home/cpicoche/Documents/Plankton/data/treated/"

library("zoo")
library("lubridate")

lieu="Teychan"

###Loading data
filename=paste(path_data_post,lieu,"_base.csv",sep="")
tabbis=read.csv(filename,na.strings="NA",header=TRUE,sep=";",dec=".")
dates=as.Date(tabbis$Date)
tab=tabbis[year(dates)>1996,]#Using data from 1997, because cryptophytes were not counted (or badly, for the first year, before
dates=dates[year(dates)>1996]
colo=c("red","black")
tab_analyse=tab[,c("SAL","CHL")]
index=which((!is.na(tab_analyse[,"CHL"]))&(!is.na(tab_analyse[,"SAL"])))
tt=tab_analyse[index,]
dates_bis=dates[index]
tt=scale(tt)

pdf("comparaison_methode_clustering_Teychan_2.pdf")
cut=2
#Plot K-means
km=kmeans(tt,cut,iter.max=100,nstart=5)
par(mfrow=c(2,2))
plot(tt,col=km$cluster,pch=16,main="K-means",xlab="Salinity",ylab="Chlorophyll a")
#Plot UPGMA and Ward's method
txt=c("UPGMA","Ward's method")
ti=0
for (m in c("average","ward.D2")){
	ti=ti+1
	di=dist(tt, method="euclidean")
	ag=hclust(di,method=m)
	plou=cutree(ag,k=cut)
	plot(tt,t="n",main=txt[ti],xlab="Salinity",ylab="Chlorophyll a")
	for (i in 1:cut){
        	indi=which(plou==i)
	        points(tt[indi,"SAL"],tt[indi,"CHL"],col=colo[i],pch=16)
	}
}
#Plot per season
plot(tt,col=plou,t="n",main="Season",xlab="Salinity",ylab="Chlorophyll a")
i=0
for(d in dates_bis){
        i=i+1
        d=as.Date(d)
        if(month(d)>=1&&month(d)<=6){
                points(tt[i,"SAL"],tt[i,"CHL"],col="red",pch=16)
        }else if(month(d)>=7&&month(d)<=12){
                points(tt[i,"SAL"],tt[i,"CHL"],col="black",pch=16)
        }
}
legend("topleft",c("Jan-Jun","Jul-Dec"),col=c("red","black"),pch=16,bty="n")
dev.off()

###########Figure S5.2
fac_axis=2
cex_pch=2
fac_lab=2
alwd=4
fac_main=3
pchtest=c(16,17,18)

pdf("dynamique_Teychan_B7.pdf",width=15,height=10)
par(mfrow=c(3,3),
          oma = c(0,0,0,0) + 0.1,
          mar = c(5,5,2,1) + 0.1)
#Plotting chl and salinity for Teychan only
ag=hclust(di,method="ward.D2")
plou=cutree(ag,k=cut)
plot(tt,col=plou,pch=16,xlab="Salinity",ylab="Chlorophyll a",cex.axis=fac_axis,cex.lab=fac_lab,cex=cex_pch,lwd=alwd)
w=0
tw=c("Chlorophyll a","Salinity")
for (v in c("CHL","SAL")){
	w=w+1
	if(w==1){
		mati="Teychan"
	}else{
		mati=""
	}

	plot(dates_bis,tt[,v],t="n",lty=3,col="grey",ylab=tw[w],main=mati,xlab="",cex.axis=fac_axis,cex.lab=fac_lab,lwd=alwd,cex.main=fac_main)
	for (i in 1:cut){
        	indi=which(plou==i)
	        points(dates_bis[indi],tt[indi,v],pch=pchtest[i],col=c("black","red")[i],cex=cex_pch)
        }
}

#Plotting chl and salinity for B7 only
lieu="B7"
###Loading data
filename=paste(path_data_post,lieu,"_base.csv",sep="")
tabbis=read.csv(filename,na.strings="NA",header=TRUE,sep=";",dec=".")
dates2=as.Date(tabbis$Date)

tab_analyse=tabbis[,c("SAL","CHL")]
index=which((!is.na(tab_analyse[,"CHL"]))&(!is.na(tab_analyse[,"SAL"])))
tti=tab_analyse[index,]
dates_bis2=dates2[index]
tt2=scale(tti)
di=dist(tt2, method="euclidean")
ag=hclust(di,method="ward.D2")

plou=cutree(ag,k=cut)
plot(tt2,col=plou,pch=16,xlab="Salinity",ylab="Chlorophyll a",cex.axis=fac_axis,cex.lab=fac_lab,cex=cex_pch,lwd=alwd)
w=0
for (v in c("CHL","SAL")){
	w=w+1
	if(w==1){
		mati="B7"
	}else{
		mati=""
	}
	plot(dates_bis2,tt2[,v],t="n",lty=3,col="grey",ylab=tw[w],xlab="",main=mati,xlim=c(dates_bis[1],dates_bis[length(dates_bis)]),cex.axis=fac_axis,cex.lab=fac_lab,cex=cex_pch,lwd=alwd,cex.main=fac_main)
	for (i in 1:cut){
        	indi=which(plou==i)
        	points(dates_bis2[indi],tt2[indi,v],pch=pchtest[i],col=c("darkgrey","orange")[i],cex=cex_pch)
        }
}

#USING BOTH
tt3=rbind(tt,tt2)

di=dist(tt3, method="euclidean")
ag=hclust(di,method="ward.D2")

plou=cutree(ag,k=cut)
plot(tt3,col=plou,pch=16,xlab="Salinity",ylab="Chlorophyll a",cex.axis=fac_axis,cex.lab=fac_lab,cex=cex_pch,lwd=alwd)
w=0
for (v in c("CHL","SAL")){
	w=w+1
	if(w==1){
		mati="Teychan+B7"
	}else{
		mati=""
	}
	plot(dates_bis,tt3[1:dim(tt)[1],v],t="n",lty=1,col="black",ylab=tw[w],xlab="",main=mati,xlim=c(dates_bis[1],dates_bis[length(dates_bis)]),cex.axis=fac_axis,cex.lab=fac_lab,cex=cex_pch,lwd=alwd,cex.main=fac_main)
	for (i in 1:cut){
        	indi=which(plou==i)
	        points(dates_bis[indi[indi<=length(dates_bis)]],tt3[indi[indi<=length(dates_bis)],v],pch=16,col=c("black","red")[i],cex=cex_pch)
        	points(dates_bis2[indi[indi>length(dates_bis)]-length(dates_bis)],tt3[indi[indi>length(dates_bis)],v],pch=17,col=c("darkgrey","orange")[i],cex=cex_pch)
	}
}
dev.off()

###########Figure S5.3 and S5.4 are drawn in the file Ricker.R

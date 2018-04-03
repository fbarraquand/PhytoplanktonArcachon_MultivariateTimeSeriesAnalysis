rm(list=ls())
graphics.off()
library("zoo")
library("lubridate")
source("functions_global.r")

#Species of interest
sp=c("GUI", "LEP", "NIT", "PSE", "RHI", "SKE", "CHA", "AST", "GYM",  "PRP", "EUG", "CRY")
sp=sort(sp)

#Data
filename=paste(path_data_post,"Teychan_base.csv",sep="")
tabbis=read.csv(filename,na.strings="NA",header=TRUE,sep=";",dec=".")
dates=as.Date(tabbis$Date)
#New
filename7=paste(path_data_post,"B7_base.csv",sep="")
tabbis7=read.csv(filename7,na.strings="NA",header=TRUE,sep=";",dec=".")
dates7=as.Date(tabbis7$Date)

pdf(paste(path_graph,"PACF_ARMA_comparison_mixte.pdf",sep=""))
par(mfrow=c(3,4))
store=array(NA,c(length(sp),3,15))
for (s in 1:length(sp)){
	#CRY and EUG have not been correctly recorded for a while
	if(sp[s]=="CRY"|sp[s]=="EUG"|sp[s]=="GYM"){
		dates_correct=dates[year(dates)>1996]
		tab=tabbis[year(dates)>1996,]
		dates_correct7=dates7[year(dates7)>1996]
		tab7=tabbis7[year(dates7)>1996,]
	}else{
		dates_correct=dates
		tab=tabbis
		dates_correct7=dates7
		tab7=tabbis7}
	dates_bis=seq(dates_correct[1],dates_correct[length(dates_correct)],14)
	dates_bis7=seq(dates_correct7[1],dates_correct7[length(dates_correct7)],14)

	#Dealing with NA and zero values for plankton
	var1=log(na.approx(tab[,sp[s]],maxgap=2,x=dates_correct,xout=dates_bis,na.rm=FALSE)+1)
	var17=log(na.approx(tab7[,sp[s]],maxgap=2,x=dates_correct7,xout=dates_bis7,na.rm=FALSE)+1)

	p=find_best_pacf(var1,var17)
	aic_mat1=find_best_arma(var1,var17,type="pp")
	a1=which(aic_mat1==min(aic_mat1,na.rm=TRUE))
	title(sp[s])
	aic_mat2=find_best_arma(var1,var17,type="pp-1")
	a2=which(aic_mat2==min(aic_mat2,na.rm=TRUE))
	aic_mat3=find_best_arma(var1,var17,type="p0")
	a3=which(aic_mat3==min(aic_mat3,na.rm=TRUE))

	store[s,1,1:length(aic_mat1)]=aic_mat1
	store[s,2,1:length(aic_mat2)]=aic_mat2
	store[s,3,1:length(aic_mat3)]=aic_mat3

	sink(paste("essai_",sp[s],"_mixte.txt",sep=""))
	print(paste("lag","pp","pp-1","p0",sep=";"))
	for(l in 1:length(aic_mat1)){
		print(paste(l,aic_mat1[l],aic_mat2[l],aic_mat3[l],sep=";"))	
	}
	sink()

	print(paste("For genus",sp[s],"pacf gives",p,"arma pp gives",a1,"arma pp-1 gives",a2,"arma p0 gives",a3))
}
dev.off()
sink()

f1="ARpp_relmin.csv"
write.table(t(sp),file=f1,na="NA",sep=";",dec=".",col.names=FALSE,row.names=FALSE,append=FALSE)
f2="ARp0_relmin.csv"
write.table(t(sp),file=f2,na="NA",sep=";",dec=".",col.names=FALSE,row.names=FALSE,append=FALSE)
f1_bis="ARpp_sort.csv"
write.table(t(sp),file=f1_bis,na="NA",sep=";",dec=".",col.names=FALSE,row.names=FALSE,append=FALSE)
f2_bis="ARp0_sort.csv"
write.table(t(sp),file=f2_bis,na="NA",sep=";",dec=".",col.names=FALSE,row.names=FALSE,append=FALSE)
f3="ARpp1_relmin.csv"
write.table(t(sp),file=f3,na="NA",sep=";",dec=".",col.names=FALSE,row.names=FALSE,append=FALSE)
for (l in 1:dim(store)[3]){
	a=store[,1,l]/apply(store[,1,],1,min,na.rm=TRUE)-1.0
        write.table(t(a),file=f1,na="NA",sep=";",dec=",",col.names=FALSE,row.names=FALSE,append=TRUE)
	a=store[,3,l]/apply(store[,3,],1,min,na.rm=TRUE)-1.0
        write.table(t(a),file=f2,na="NA",sep=";",dec=",",col.names=FALSE,row.names=FALSE,append=TRUE)
	
	a=store[,2,l]/apply(store[,2,],1,min,na.rm=TRUE)-1.0
        write.table(t(a),file=f3,na="NA",sep=";",dec=",",col.names=FALSE,row.names=FALSE,append=TRUE)
}



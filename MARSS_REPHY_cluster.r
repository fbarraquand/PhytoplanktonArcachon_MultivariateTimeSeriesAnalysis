###############################################################################################################################################
                                             ########CP - Specific MAR model for CHA and AST ########
###############################################################################################################################################

###Setting path and library
rm(list=ls())
graphics.off()

library("zoo")
library("lubridate")
library("MARSS")
lieu="Teychan" #can also be B7
source("functions_global.r")


###Phytoplankton, guild definitions
sp=c("AST","CHA")
###Nutrient definition
cov3_phy=c("SAL","CumRg","MeanVent")
cov3_tot=c(cov3_phy,"TEMP","CHL") #covariates we are going to use, including temperature and chlorophyll for cluster computation

###Loading data
filename=paste(lieu,"_base.csv",sep="")
tabbis=read.csv(filename,na.strings="NA",header=TRUE,sep=";",dec=".")
dates=as.Date(tabbis$Date)
tab=tabbis[year(dates)>1996,]#Using data from 1997, because cryptophytes were not counted before (or badly, for the first year)
dates=dates[year(dates)>1996]
tab_sp=tab[,sp]

###Treating covariates to have cumulated radiation and wind energy
accum_vent=3 #We take into account wind energy 3 days before the date
tab_cov=new_covar(tab,dates,accum_vent)

###Treating missing values for both species and covariates###
consecutif=2 #Number of missing values above which we keep the NA
timestep=14 #Regular time lapses between two observations
dates_bis=seq(dates[1],dates[length(dates)],timestep) #Regular time grid
tab_sp=na.approx(tab_sp,maxgap=consecutif,x=dates,xout=dates_bis,na.rm=FALSE) #Interpolation over regular time grid
#Replace missing values
for (s in sp){
	tab_sp[is.na(tab_sp[,s]),s]=runif(sum(is.na(tab_sp[,s])),0,min(tab_sp[,s],na.rm=TRUE))
}
tab_cov_bis=matrix(NA,length(dates_bis),length(cov3_tot))
colnames(tab_cov_bis)=cov3_tot
for (c in cov3_tot){
	tab_cov_bis[,c]=approx(tab_cov[,c],x=dates,xout=dates_bis)$y
}

#Log transfo for species abundance
tab_sp=log(tab_sp)

#Defining clusters
cut=2 #Number of cluster we want
di=dist(tab_cov_bis[,c("CHL","SAL")],method="euclidean")
ag=hclust(di,method="ward.D2")
plou=cutree(ag,k=cut) #contains the indices of the different groups

#For each group (bloom and non-bloom conditions)
for (c in 1:cut){
	indi=which(plou!=c)
	tab_sp_phase=tab_sp
	tab_sp_phase[indi,]=NA

	tab_cov_phase=tab_cov_bis[,cov3_phy]
	
	tab_sp_phase=t(scale(tab_sp_phase[3:(length(dates_bis)),]))
	tab_cov_phase=t(scale(tab_cov_phase[2:(length(dates_bis)-1),]))
	rownames(tab_sp_phase)=sp
	rownames(tab_cov_phase)=cov3_phy

	#Setting null model
	B1="diagonal and unequal"
	analyse_MARSS(tab_sp_phase,tab_cov_phase,B1,paste(lieu,"_physics_null_phase_bis",c,".RData",sep=""))
	#Setting full model
	B1="unconstrained"
	analyse_MARSS(tab_sp_phase,tab_cov_phase,B1,paste(lieu,"_physics_unconstrained_phase_bis",c,".RData",sep=""))
}




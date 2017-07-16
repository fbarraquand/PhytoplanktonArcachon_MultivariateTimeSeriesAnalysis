###############################################################################################################################################
########CP - Core of MARSS analyses on real data, estimating coefficients for different scenarios of interactions between species #############
###############################################################################################################################################


###Setting path and library
rm(list=ls())
graphics.off()

library("zoo")
library("lubridate")
library("MARSS")
source("functions_global.r")
lieu="Teychan" #Or "B7"

###Phytoplankton, guild definitions
sp=c("AST","NIT","PSE","SKE","CHA","GUI","LEP","RHI","GYM","PRP","CRY","EUG")
sp_pen=c("AST","NIT","PSE") #pennate diatoms
sp_cen=c("SKE","CHA","GUI","LEP","RHI") #centric diatoms
sp_diat=c(sp_pen,sp_cen) #diatoms
sp_din=c("GYM","PRP") #dinoflagellates
###Nutrient definition
cov3_phy=c("SAL","CumRg","MeanVent") #physics-only model includes salinity, irradiance and wind energy
cov3_tot=c(cov3_phy,"TEMP") #covariates we are going to use, including temperature used for the dedicated seasonal variable

###Loading data
filename=paste(lieu,"_base.csv",sep="")
tabbis=read.csv(filename,na.strings="NA",header=TRUE,sep=";",dec=".")
dates=as.Date(tabbis$Date)
tab=tabbis[year(dates)>1996,]#Using data from 1997, because cryptophytes were not counted (or badly, for the first year), before
dates=dates[year(dates)>1996]
tab_sp=tab[,sp]

###Treating covariates to use cumulated irradiance and wind energy
accum_vent=3 #We take into account wind energy 3 days before the date
tab_cov=new_covar(tab,dates,accum_vent)

###Treating missing values for both species and covariates###
consecutif=2 #Number of missing values above which we keep NA in the planktonic data instead of interpolating them
timestep=14 #Regular time lapses between two observations
dates_bis=seq(dates[1],dates[length(dates)],timestep) #Regular time grid
tab_sp=na.approx(tab_sp,maxgap=consecutif,x=dates,xout=dates_bis,na.rm=FALSE) #Interpolation over regular time grid
#Replace missing values by random values for planktonic data
for (s in sp){
	tab_sp[is.na(tab_sp[,s]),s]=runif(sum(is.na(tab_sp[,s])),0,min(tab_sp[,s],na.rm=TRUE))
}
tab_cov_bis=matrix(NA,length(dates_bis),length(cov3_tot))
colnames(tab_cov_bis)=cov3_tot
#Interpolation for covariates
for (c in cov3_tot){
	tab_cov_bis[,c]=approx(tab_cov[,c],x=dates,xout=dates_bis)$y
}

###Seasonal component: first we project the temperature on a trigonometric function to have and seasonal component and then use only the residuals of the other covariates
t=1:length(tab_cov_bis[,'TEMP'])
omega=2*pi/(365.25/14)
sini=sin(omega*t)
cosi=cos(omega*t)
gap_temp=!is.na(tab_cov_bis[,'TEMP'])
season=lm(tab_cov_bis[,'TEMP']~sini+cosi)
for (c in cov3_tot){
	tab_cov_bis[gap_temp&!is.na(tab_cov_bis[,c]),c]=lm(tab_cov_bis[gap_temp,c]~season$fitted.values)$residuals
        tab_cov_bis[is.na(tab_cov_bis[,'TEMP']),c]=NA
                }
tab_cov_bis[,length(cov3_tot)]=season$fitted.values
cov3_tot=c(cov3_tot[1:(length(cov3_tot)-1)],"season")

#Log transformation for species abundance and scaling for all time series
tab_sp=log(tab_sp)
tab_sp=t(scale(tab_sp[2:(length(dates_bis)),]))
tab_cov=t(scale(tab_cov_bis[2:(length(dates_bis)),]))
rownames(tab_sp)=sp
rownames(tab_cov)=cov3_tot

#Estimating the null model
B1="diagonal and unequal"
analyse_MARSS(tab_sp,tab_cov,B1,paste(lieu,"_physics_null_seasonal_bis.RData",sep=""))

#Estimating the full model
B1="unconstrained"
analyse_MARSS(tab_sp,tab_cov,B1,paste(lieu,"_physics_unconstrained_seasonal_bis.RData",sep=""))

#Estimating a model allowing no interaction between pennate and centric diatoms. There is only one specific case of dinoflagellate being still able to feed on cryptophytes.
B2=matrix(list(0),nrow=length(sp),ncol=length(sp),dimnames=list(sp,sp))
for (i in 1:length(sp)){
        s=sp[i]
        for (j in 1:length(sp)){
                s2=sp[j]
                if(((s %in% sp_pen)&(s2%in% sp_pen))||((s %in% sp_cen)&(s2%in% sp_cen))||s==s2||((s=="CRY")&(s2 %in% sp_din))||((s2=="CRY")&(s %in% sp_din))||((s %in% sp_din)&(s2%in% sp_din))){
                        B2[i,j]=paste(s2,s,sep="")
                }else{
                        if(s==s2){
                                B2[i,j]=paste(s2,s,sep="")
                        }
                }
        }
}
analyse_MARSS(tab_sp,tab_cov,B2,paste(lieu,"_physics_pencen_seasonal_bis.RData",sep=""))

#Estimating a model allowing no interaction between diatoms, dinoflagellates, crypotphytes and euglenophytes. There is only one specific case of dinoflagellate being still able to feed on cryptophytes
B3=matrix(list(0),nrow=length(sp),ncol=length(sp),dimnames=list(sp,sp))
for (i in 1:length(sp)){
        s=sp[i]
        for (j in 1:length(sp)){
                s2=sp[j]
                if(((s %in% sp_diat)&(s2%in% sp_diat))||s==s2||((s=="CRY")&(s2 %in% sp_din))||((s2=="CRY")&(s %in% sp_din))||((s %in% sp_din)&(s2%in% sp_din))){
                        B3[i,j]=paste(s2,s,sep="")
                }else{
                        if(s==s2){
                                B3[i,j]=paste(s2,s,sep="")
                        }
                }
        }
}
analyse_MARSS(tab_sp,tab_cov,B3,paste(lieu,"_physics_diatdin_seasonal_bis.RData",sep=""))

#Estimating a model allowing no interactions intra-group, but interactions inter-group
B4=matrix(list(0),nrow=length(sp),ncol=length(sp),dimnames=list(sp,sp))
for (i in 1:length(sp)){
        s=sp[i]
        for (j in 1:length(sp)){
                s2=sp[j]
                if(((s %in% sp_diat)&!(s2%in% sp_diat))||(s%in% c("CRY","EUG"))||(s2%in% c("CRY","EUG"))||((s %in% sp_din)&!(s2%in% sp_din))||((s2 %in% sp_diat)&!(s%in% sp_diat))||((s2 %in% sp_din)&!(s%in% sp_din))){
                        B4[i,j]=paste(s2,s,sep="")
                }else{
                        if(s==s2){
                                B4[i,j]=paste(s2,s,sep="")
                        }
                }
        }
}
analyse_MARSS(tab_sp,tab_cov,B4,paste(lieu,"_physics_inverse_seasonal_bis.RData",sep=""))


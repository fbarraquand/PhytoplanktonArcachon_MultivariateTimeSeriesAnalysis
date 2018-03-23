##### FB - 30/11/2016 - Try other nonlinear competition models for CHA / AST, just to check there are no interactions
##### Also check pairwise Granger causality in passing

###Setting path and library
rm(list=ls())
graphics.off()

library("zoo")
library("lubridate")
library("MARSS")
source("functions_global.r")

###Phytoplankton, guild definitions
sp=c("AST","NIT","PSE","SKE","CHA","GUI","LEP","RHI","GYM","PRP","CRY","EUG")
###Nutrient definition
cov3_phy=c("SAL","CumRg","MeanVent")
cov3_nut=c('Ntot','SI','PHOS')
cov3_tot=c(cov3_phy,"TEMP") #covariates we are going to use, including temperature

###Loading data
lieu="Teychan"
filename=paste(lieu,"_base.csv",sep="")
tabbis=read.csv(filename,header=TRUE,sep=";",dec=".")#na.strings="NA"
dates=as.Date(tabbis$Date)
tab=tabbis
tab_sp=tab[,sp]

###Treating covariates to have cumulated irradiance and wind energy
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

### For a non-linear model, we can abitrarily decide that above 11 (log(N)), 
### it is bloom density, and have for each species a bloom variable = 1 or 0
### Then we use the growth rate ~ lm (log(CHA)+log(AST)+log(TEMP))
bloom_thresh=11
tab$AST_bloom=0
tab$CHA_bloom=0
tab$AST_bloom[log(tab$AST)>bloom_thresh]=1
tab$CHA_bloom[log(tab$CHA)>bloom_thresh]=1

### Each genus condition (bloom or non-bloom) affects the effects of densities on them
n=nrow(tab) # size dataset
### AST
tab$L_AST=c(log(tab$AST[1:(n-1)]),NA)
tab$L_AST_plus1 = c(log(tab$AST[2:n]),NA)
tab$r_AST=c((tab$L_AST_plus1[1:(n-1)] - tab$L_AST[1:(n-1)]),NA) #AST growth rates
### CHA
tab$L_CHA=c(log(tab$CHA[1:(n-1)]),NA)
tab$L_CHA_plus1 = c(log(tab$CHA[2:n]),NA)
tab$r_CHA=c(tab$L_CHA_plus1[1:(n-1)] - tab$L_CHA[1:(n-1)],NA) #CHA growth rates

#### Covariates
tab_cov = as.data.frame(tab_cov)

### Regression model independent of bloom periods
lm_AST = lm(tab$r_AST ~ tab$L_AST + tab$L_CHA + tab$TEMP)
print(summary(lm_AST))
lm_CHA= lm(tab$r_CHA~ tab$L_CHA + tab$L_AST + tab$TEMP)
print(summary(lm_CHA))

### Try pairwise Granger causality. 
require(lmtest)
### Now we need to center absolutely (otherwise big mistakes can be done)
tab$L_AST=log(tab$AST)-mean(log(tab$AST),na.rm=TRUE)
tab$L_CHA=log(tab$CHA)-mean(log(tab$CHA),na.rm=TRUE)
grangertest(tab$L_AST,tab$L_CHA,order=3)
# Model 1: tab$L_CHA ~ Lags(tab$L_CHA, 1:3) + Lags(tab$L_AST, 1:3)
# Model 2: tab$L_CHA ~ Lags(tab$L_CHA, 1:3)
# Res.Df Df      F Pr(>F)
# 1    220                 
# 2    223 -3 1.3903 0.2466
grangertest(tab$L_CHA,tab$L_AST,order=3)
# Granger causality test
# 
# Model 1: tab$L_AST ~ Lags(tab$L_AST, 1:3) + Lags(tab$L_CHA, 1:3)
# Model 2: tab$L_AST ~ Lags(tab$L_AST, 1:3)
# Res.Df Df      F  Pr(>F)  
# 1    220                    
# 2    223 -3 2.1281 0.09756 .
### Slight effect of CHA on AST due to confounding effects of irradiance, again? 
### Under that assumption it should disappear whenever we include CumRg into the model.
grangertest(tab$L_CHA,tab$L_AST) # slight effect

### Granger causality using 1 lag and a conditioning CumRg factor. 
### On AST dynamics
lm_full = lm(tab$r_AST[!is.na(tab$L_CHA)] ~ tab$L_AST[!is.na(tab$L_CHA)] + tab$L_CHA[!is.na(tab$L_CHA)]  + tab_cov$CumRg[!is.na(tab$L_CHA)])
lm_null = lm(tab$r_AST[!is.na(tab$L_CHA)] ~ tab$L_AST[!is.na(tab$L_CHA)] + tab_cov$CumRg[!is.na(tab$L_CHA)])
### The !is.na() yields the same number of observations which makes the AIC comparable and the Wald test possible
AIC(lm_full,lm_null)
waldtest(lm_full,lm_null) ### These are not identified as different...

### Add conditioning on TEMP (both season and irradiance taken into account)
lm_full = lm(tab$r_AST ~ tab$L_AST + tab$L_CHA +tab$TEMP + tab_cov$CumRg)
lm_null = lm(tab$r_AST ~ tab$L_AST +tab$TEMP + tab_cov$CumRg)
AIC(lm_full,lm_null) ## AST and CHA do not have the same number of NA, may bias computation. 
### Let's check that
lm_full = lm(tab$r_AST[!is.na(tab$L_CHA)] ~ tab$L_AST[!is.na(tab$L_CHA)] + tab$L_CHA[!is.na(tab$L_CHA)] +tab$TEMP[!is.na(tab$L_CHA)] + tab_cov$CumRg[!is.na(tab$L_CHA)])
lm_null = lm(tab$r_AST[!is.na(tab$L_CHA)] ~ tab$L_AST[!is.na(tab$L_CHA)] +tab$TEMP[!is.na(tab$L_CHA)] + tab_cov$CumRg[!is.na(tab$L_CHA)])
### The !is.na() yields the same number of observations which makes the AIC comparable and the Wald test possible
AIC(lm_full,lm_null)
waldtest(lm_full,lm_null) ### These are not identified as different...

### Does not change anything to the previous fit
### CHA seem to improve the fit yet the coefficient is not significantly 
### different from zero...

### Scale everyting to see what coeffs are most influential
tab$TEMP=scale(tab$TEMP)
tab_cov$CumRg=scale(tab_cov$CumRg)

### AST
tab$L_AST=scale(c(log(tab$AST[1:(n-1)]),NA))
tab$L_AST_plus1 = scale(c(log(tab$AST[2:n]),NA))
tab$r_AST=c((tab$L_AST_plus1[1:(n-1)] - tab$L_AST[1:(n-1)]),NA)
### CHA
tab$L_CHA=scale(c(log(tab$CHA[1:(n-1)]),NA))
tab$L_CHA_plus1 = scale(c(log(tab$CHA[2:n]),NA))
tab$r_CHA=c(tab$L_CHA_plus1[1:(n-1)] - tab$L_CHA[1:(n-1)],NA)
### 

### Add conditioning on TEMP (both season and irradiance taken into account)
lm_full = lm(tab$r_AST ~ tab$L_AST + tab$L_CHA +tab$TEMP + tab_cov$CumRg)
lm_null = lm(tab$r_AST ~ tab$L_AST +tab$TEMP + tab_cov$CumRg)
AIC(lm_full,lm_null) ## AST and CHA do not have the same number of NA, may bias computation. 
summary(lm_full)
summary(lm_null)

### Let's check that
lm_full = lm(tab$r_AST[!is.na(tab$L_CHA)] ~ tab$L_AST[!is.na(tab$L_CHA)] + tab$L_CHA[!is.na(tab$L_CHA)] +tab$TEMP[!is.na(tab$L_CHA)] + tab_cov$CumRg[!is.na(tab$L_CHA)])
lm_null = lm(tab$r_AST[!is.na(tab$L_CHA)] ~ tab$L_AST[!is.na(tab$L_CHA)] +tab$TEMP[!is.na(tab$L_CHA)] + tab_cov$CumRg[!is.na(tab$L_CHA)])
### The !is.na() yields the same number of observations which makes the AIC comparable and the Wald test possible
AIC(lm_full,lm_null)
waldtest(lm_full,lm_null) ### These are not identified as different...


############### Check with Ricker model 05/02/2017 #####################################
### Regression model independent of bloom periods
lm_AST_Ricker = lm(tab$r_AST ~ tab$AST + tab$CHA + tab$TEMP)
summary(lm_AST_Ricker)
lm_CHA_Ricker= lm(tab$r_CHA~ tab$CHA + tab$AST + tab$TEMP)
summary(lm_CHA_Ricker)
## Some effects there, but the model would have to fit
plot(lm_AST_Ricker)
plot(lm_CHA_Ricker)

plot(tab$AST,tab$r_AST)
plot(tab$AST,tab$r_CHA)
plot(tab$CHA,tab$r_AST)
### Difficult to see anything on those plots...
### What is clear is that the R^2 decreased (slightly). 
### Perhaps we need to simulate the Gompertz and the Ricker model
plot(tab$L_AST,tab$r_CHA)
plot(tab$AST,tab$r_CHA,xlim=c(0,1000000))
abline(a=coef(lm_CHA_Ricker)[1],b=coef(lm_CHA_Ricker)[3])
### This does not look very convincing... but still we need to check with simulations. 
### I need to extract the noise values for the simulations. 

r1=coef(lm_AST_Ricker)[1]
alpha11=coef(lm_AST_Ricker)[2]
alpha12=coef(lm_AST_Ricker)[3]
c1=coef(lm_AST_Ricker)[4]
sigma1=sigma(lm_AST_Ricker)

r2=coef(lm_CHA_Ricker)[1]
alpha22=coef(lm_CHA_Ricker)[2]
alpha21=coef(lm_CHA_Ricker)[3]
c2=coef(lm_CHA_Ricker)[4]
sigma2=sigma(lm_CHA_Ricker)

tmax=nrow(tab)
N=matrix(1,2,tmax)
y=tab$TEMP

for (t in 1:(tmax-1)){
  N[1,t+1] = N[1,t]*exp(r1+alpha11*N[1,t]+alpha12*N[2,t]+c1*y[t+1]+rnorm(1,0,sigma1))
  N[2,t+1] = N[2,t]*exp(r2+alpha21*N[1,t]+alpha22*N[2,t]+c2*y[t+1]+rnorm(1,0,sigma2))
}


#Plot parameters
fac_main=3.0
fac_axis=2.0
fac_lg=1.5
fac_lab=1.9
alwd=2.5
apc=2
awidth=15

pdf("Ricker_simulation_clean.pdf",width=awidth)
par(mar=c(5,6,1,1))
plot(1:tmax,log(N[1,]),col="green",type="o",lwd=alwd,cex.main=fac_main,cex.axis=fac_axis,cex.lab=fac_lab,cex=apc,main="",ylab="Planktonic log(abundance)",xlab="Time",xlim=c(0,tmax+50))
lines(1:tmax,log(N[2,]),col="blue",type="o",lwd=alwd,cex=apc)
legend("bottomleft",c("AST","CHA"),col=c("green","blue"),lty=1,cex=fac_main,bty="n",lwd=alwd,pch="o")
dev.off()

#### Now do the same thing with the log-linear model

r1=coef(lm_AST)[1] # Not a maximal growth rate there
alpha11=coef(lm_AST)[2]
alpha12=coef(lm_AST)[3]
c1=coef(lm_AST)[4]
sigma1=sigma(lm_AST)

r2=coef(lm_CHA)[1]
alpha22=coef(lm_CHA)[2]
alpha21=coef(lm_CHA)[3]
c2=coef(lm_CHA)[4]
sigma2=sigma(lm_CHA)

tmax=nrow(tab)
n=matrix(1,2,tmax)
y=tab$TEMP

for (t in 1:(tmax-1)){
  n[1,t+1] = n[1,t]+r1+alpha11*n[1,t]+alpha12*n[1,t]+c1*y[t+1]+rnorm(1,0,sigma1)
  n[2,t+1] = n[2,t]+r2+alpha21*n[1,t]+alpha22*n[2,t]+c2*y[t+1]+rnorm(1,0,sigma2)
}
pdf("Gompertz_simulation_clean.pdf",width=awidth)
par(mar=c(5,6,1,1))
plot(1:tmax,n[1,],col="green",type="o",lwd=alwd,cex.main=fac_main,cex.axis=fac_axis,cex.lab=fac_lab,cex=apc,main="",ylab="Planktonic log(abundance)",xlab="Time",xlim=c(0,tmax+50),ylim=c(0,20))
lines(1:tmax,n[2,],col="blue",type="o",lwd=alwd,cex=apc)
legend("topright",c("AST","CHA"),col=c("green","blue"),lty=1,cex=fac_main,lwd=alwd,bty="n",pch="o")
dev.off()
### Much more resembling dynamics. 




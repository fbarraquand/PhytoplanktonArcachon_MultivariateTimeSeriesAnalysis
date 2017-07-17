##### CP - 26/06/2017 - New model for phase-dep with environmental conditions, using script from FB
##### FB - 30/11/2016 - Try phase-dep and nonlinear competition models for CHA / AST, just to check there are no interactions

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
cov3_phy=c("CumRg") #Physics-based variables can be irradiance or temperature
cov3_tot=c(cov3_phy,"SAL","TEMP","CHL") #covariates we are going to use, including temperature and chlorophyll for cluster computation

sink("phase_dependent_growth.txt")
for(lieu in c("Teychan","B7")){

#Loading data	
filename=paste(path_data_post,lieu,"_base.csv",sep="")
tabbis=read.csv(filename,header=TRUE,sep=";",dec=".")
dates=as.Date(tabbis$Date)
tab=tabbis
tab_sp=tab[,sp]

###Treating covariates to have cumulated radiation and wind energy
accum_vent=3 #We take into account wind energy 3 days before the date
tab_cov=new_covar(tab,dates,accum_vent)

### For a non-linear model, we can set 11 in log abundance as the bloom density
### threshold. For for each species, we have a bloom variable = 1 or 0
### Then we use the growth rate ~ lm (log(CHA)+log(AST)+log(TEMP/CumRg))
### Fit a Threshold AutoRegressive model.
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
tab$r_AST=c((tab$L_AST_plus1[1:(n-1)] - tab$L_AST[1:(n-1)]),NA)
### CHA
tab$L_CHA=c(log(tab$CHA[1:(n-1)]),NA)
tab$L_CHA_plus1 = c(log(tab$CHA[2:n]),NA)
tab$r_CHA=c(tab$L_CHA_plus1[1:(n-1)] - tab$L_CHA[1:(n-1)],NA)

#### Covariates
tab_cov = as.data.frame(tab_cov)

### Regression model independent of bloom periods
print(lieu)
print("TEMPERATURE")
print("All data")
lm_AST = lm(tab$r_AST ~ tab$L_AST + tab$L_CHA + tab$TEMP)
print(summary(lm_AST))
lm_CHA= lm(tab$r_CHA~ tab$L_CHA + tab$L_AST + tab$TEMP)
print(summary(lm_CHA))

### Bloom model
print("Bloom")
lm_AST = lm(tab$r_AST[tab$AST_bloom==1] ~ tab$L_AST[tab$AST_bloom==1]  + tab$L_CHA[tab$AST_bloom==1]  + tab$TEMP[tab$AST_bloom==1] )
print(summary(lm_AST))
lm_CHA= lm(tab$r_CHA[tab$CHA_bloom==1] ~ tab$L_CHA[tab$CHA_bloom==1] + tab$L_AST[tab$CHA_bloom==1] + tab$TEMP[tab$CHA_bloom==1])
print(summary(lm_CHA))

### Non-bloom model
print("Non-bloom")
lm_AST = lm(tab$r_AST[tab$AST_bloom==0] ~ tab$L_AST[tab$AST_bloom==0]  + tab$L_CHA[tab$AST_bloom==0]  + tab$TEMP[tab$AST_bloom==0] )
print(summary(lm_AST)) ### Stronger effect of CHA there. 
lm_CHA= lm(tab$r_CHA[tab$CHA_bloom==0] ~ tab$L_CHA[tab$CHA_bloom==0] + tab$L_AST[tab$CHA_bloom==0] + tab$TEMP[tab$CHA_bloom==0])
print(summary(lm_CHA))

print("CUMRG")
print("All data")
lm_AST = lm(tab$r_AST ~ tab$L_AST + tab$L_CHA + tab_cov$CumRg)
print(summary(lm_AST))
lm_CHA= lm(tab$r_CHA~ tab$L_CHA + tab$L_AST + tab_cov$CumRg)
print(summary(lm_CHA))

### Bloom model
print("Bloom")
lm_AST = lm(tab$r_AST[tab$AST_bloom==1] ~ tab$L_AST[tab$AST_bloom==1]  + tab$L_CHA[tab$AST_bloom==1]  + tab_cov$CumRg[tab$AST_bloom==1] )
print(summary(lm_AST))
lm_CHA= lm(tab$r_CHA[tab$CHA_bloom==1] ~ tab$L_CHA[tab$CHA_bloom==1] + tab$L_AST[tab$CHA_bloom==1] + tab_cov$CumRg[tab$CHA_bloom==1])
print(summary(lm_CHA))

### Non-bloom model
print("Non-bloom")
lm_AST = lm(tab$r_AST[tab$AST_bloom==0] ~ tab$L_AST[tab$AST_bloom==0]  + tab$L_CHA[tab$AST_bloom==0]  + tab_cov$CumRg[tab$AST_bloom==0])
print(summary(lm_AST)) ### Stronger effect of CHA there. 
lm_CHA= lm(tab$r_CHA[tab$CHA_bloom==0] ~ tab$L_CHA[tab$CHA_bloom==0] + tab$L_AST[tab$CHA_bloom==0] + tab_cov$CumRg[tab$CHA_bloom==0])
print(summary(lm_CHA))

#########################CLUSTER###############
### Here, the model is the same but the threshold is based on environmental conditions
### whose multidimensional nature is summarized by clustering	
cut=2 #We want two clusters of environmental conditions
di=dist(tab_cov[,c("CHL","SAL")],method="euclidean") #We build a distance matrix
ag=hclust(di,method="ward.D2") #This distance matrix leads to clustering using ward method
plou=cutree(ag,k=cut) #We extract two clusters from this tree

#For each phase, we build a linear model for temperature and irradiance	
for (c in 1:cut){
	print(paste("Phase",c,sep=" "))
        indi=which(plou==c)
	print("TEMPERATURE")
        lm_AST = lm(tab$r_AST[indi] ~ tab$L_AST[indi]  + tab$L_CHA[indi]  + tab$TEMP[indi] )
        print(summary(lm_AST))
        lm_CHA= lm(tab$r_CHA[indi] ~ tab$L_CHA[indi] + tab$L_AST[indi] + tab$TEMP[indi])
        print(summary(lm_CHA))
	print("CUMRG")
	lm_AST = lm(tab$r_AST[indi] ~ tab$L_AST[indi]  + tab$L_CHA[indi]  + tab_cov$CumRg[indi] )
	print(summary(lm_AST))
	lm_CHA= lm(tab$r_CHA[indi] ~ tab$L_CHA[indi] + tab$L_AST[indi] + tab_cov$CumRg[indi])
	print(summary(lm_CHA))


}
}
sink()

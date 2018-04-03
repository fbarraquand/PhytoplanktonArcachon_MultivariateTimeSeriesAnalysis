###############################################################################################################################################
########################################## CP - Autoregressive model with a fixed order and covariates ########################################
###############################################################################################################################################
#This function is an interface for the ARLM_function, which computes linear models for the different species in sp, and for the covariates 
#defined in acov, for Teychan, B7, or all sites. The configuration described here corresponds to Table 3 in the manuscript (full model,
#both places taken into account, and different uses of the nutrient values

rm(list=ls())
graphics.off()
library("zoo")
library("lubridate")
library("perturb") #for collinearity index
library("sme") #for AICc
library("Hmisc") #for errbar
source("functions_global.r")
source("ARLM_function.r")

sp=c("AST","NIT","PSE","SKE","CHA","GUI","LEP","RHI","GYM","PRP","CRY","EUG")
sp=sort(sp)
place="All" #Could be "All" to use both sites, B7 or Teychan 
acov="Tot" #Could be "Phy","Nut" or "Tot"

tab_stor_AICc=matrix(NA,nrow=length(sp),ncol=6)
rownames(tab_stor_AICc)=sp
colnames(tab_stor_AICc)=c("No Season No saturation","Season No Saturation", "No Season Saturation","Season Saturation","No Season Ratio","Season Ratio")

#ARLM function (explicit seasonal component (TRUE) or raw values (FALSE), variables in the model (Tot, Nut, Phy), place (All, Teychan or B7), ratio of nutrients (TRUE) or raw values (FALSE), saturating function (TRUE) or raw values (FALSE), species to study)
#No Season No saturation
aseason=FALSE
aratio=FALSE
astress=FALSE
tab_stor_AICc[,1]=ARLM_function(aseason,acov,place,aratio,astress,sp)

#Season No saturation
aseason=TRUE
aratio=FALSE
astress=FALSE
tab_stor_AICc[,2]=ARLM_function(aseason,acov,place,aratio,astress,sp)

#No Season Saturation
aseason=FALSE
aratio=FALSE
astress=TRUE
tab_stor_AICc[,3]=ARLM_function(aseason,acov,place,aratio,astress,sp)

#Season Saturation
aseason=TRUE
aratio=FALSE
astress=TRUE
tab_stor_AICc[,4]=ARLM_function(aseason,acov,place,aratio,astress,sp)

#No Season Ratio
aseason=FALSE
aratio=TRUE
astress=FALSE
tab_stor_AICc[,5]=ARLM_function(aseason,acov,place,aratio,astress,sp)

#Season Ratio
aseason=TRUE
aratio=TRUE
astress=FALSE
tab_stor_AICc[,6]=ARLM_function(aseason,acov,place,aratio,astress,sp)

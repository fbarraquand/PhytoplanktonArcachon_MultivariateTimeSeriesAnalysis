###############################################################################################################################################
########################################## CP - Autoregressive model with a fixed order and covariates ########################################
###############################################################################################################################################

rm(list=ls())
graphics.off()
library("zoo")
library("lubridate")
library("perturb") #for collinearity index
library("sme") #for AICc
library("Hmisc") #for errbar
source("functions_global.r")

#Planktonic species we study
sp=c("AST","NIT","PSE","SKE","CHA","GUI","LEP","RHI","GYM","PRP","CRY","EUG")

#Covariates we use
cov3_phy=c("CumRg","MeanVent","SAL") #physics-only model
cov3_nut=c('Ntot','PHOS') #nutrient model
cov3_tot=c(cov3_phy) #Here, cov3_tot can aggregate several covariate vectors described above
cov3=sort(cov3_tot)
cov_saison=cov3

#Data from Teychan
filename=paste("Teychan_base.csv",sep="")
tabbis=read.csv(filename,na.strings="NA",header=TRUE,sep=";",dec=".")
dates=as.Date(tabbis$Date)

#Data from B7
filename7=paste("B7_base.csv",sep="")
tabbis7=read.csv(filename7,na.strings="NA",header=TRUE,sep=";",dec=".")
dates7=as.Date(tabbis7$Date)

#Covariate pre-processing
accum_vent=3 #We take into account wind energy 3 days before the date
tab_cov=new_covar(tabbis,dates,accum_vent,start="t") #This function transforms data to have the cumulated irradiance and wind energy over the whole time series and the sum of NH4 and NOx for nitrogen input
tab_cov7=new_covar(tabbis7,dates7,accum_vent,start="t")

nb_lag=3 #Number of lag determined in a previous analysis

centr_reduit=TRUE #Boolean for scaled variable
desaisonnalise=TRUE #Boolean for a dedicated seasonal component
use_stress=FALSE #Boolean for saturating function
doyouplot=FALSE #Boolean to plot the residuals of the linear model
figure_4=TRUE #Boolean to plot figure 4

#Options for plotting
max_ci=0.0
avg_ci=0.0
fac_main=3.0
fac_axis=2.0
fac_lg=1.5
fac_lab=1.5
alwd=2.5
apc=2
awidth=15

#We output text results in fileout
fileout="lm_test.txt"
sink(fileout)
print(paste("Simulation uses covariates",paste(cov3,collapse=" "),"with z-scored",centr_reduit,", de-season",desaisonnalise,"use_stress",use_stress,"and nb_lag",nb_lag))

for (ss in 1:length(sp)){
	s=sp[ss]
	#CRY and EUG have not been correctly recorded for a while, so we have to change the dates
        if(s=="CRY"|s=="EUG"|s=="GYM"){
                dates_correct=dates[year(dates)>1996]
                tab=tabbis[year(dates)>1996,]
                tab_covbis=tab_cov[year(dates)>1996,]

                dates_correct7=dates7[year(dates7)>1996]
                tab7=tabbis7[year(dates7)>1996,]
                tab_covbis7=tab_cov7[year(dates7)>1996,]
        }else{
                dates_correct=dates
                tab=tabbis
                tab_covbis=tab_cov

                dates_correct7=dates7
                tab7=tabbis7
                tab_covbis7=tab_cov7
        }
#Interpolation over a regular sampling time sequence
        dates_bis=seq(dates_correct[1],dates_correct[length(dates_correct)],14)
        dates_bis7=seq(dates_correct7[1],dates_correct7[length(dates_correct7)],14)

        #Dealing with NA and zero values for plankton
        var1=na.approx(tab[,s],maxgap=2,x=dates_correct,xout=dates_bis,na.rm=FALSE)
        var17=na.approx(tab7[,s],maxgap=2,x=dates_correct7,xout=dates_bis7,na.rm=FALSE)
        #Growth rates
	tx_croiss=diff(log(var1))
        tx_croiss7=diff(log(var17))

	newvar1=log(var1)
	newvar17=log(var17)

#Growth adjustment for different time lags and concatenation of observations from Teychan and B7
        tx_croissbis=tx_croiss[nb_lag:length(tx_croiss)]
        tx_croissbis7=tx_croiss7[nb_lag:length(tx_croiss7)]
#mat_ab contains the log abundance values, with a shift corresponding to the time lag for the following linear model
        mat_ab=matrix(NA,ncol=nb_lag,nrow=length(tx_croissbis))
        mat_ab7=matrix(NA,ncol=nb_lag,nrow=length(tx_croissbis7))
        for (i in 0:(nb_lag-1)){
                mat_ab[,i+1]=newvar1[(nb_lag-i):(length(newvar1)-1-i)]
                mat_ab7[,i+1]=newvar17[(nb_lag-i):(length(newvar17)-1-i)]
        }
        tx_croiss_concat=c(tx_croissbis,tx_croissbis7)
        mat_ab_concat=rbind(mat_ab,mat_ab7)

        #Same thing for covariates
        var2=na.approx(tab_covbis,maxgap=1,x=dates_correct,xout=dates_bis,na.rm=FALSE)
        var3=var2[nb_lag:(dim(var2)[1]-1),]
        var27=na.approx(tab_covbis7,maxgap=1,x=dates_correct7,xout=dates_bis7,na.rm=FALSE)
        var37=var27[nb_lag:(dim(var27)[1]-1),]

#Use of a saturating value for nitrogen. stress_function can be found in the file functions_global.r
	if(use_stress){
		if (s=="CRY"){ #cryptophytes
			threshold=1.2*2
		}else if(s=="EUG"){ #euglenophytes
			threshold=2.0*2
		}else if(s %in% c("PRP","GYM")){ #dinoflagellates
			threshold=7.0*2
		}else{ #diatoms
			threshold=1.25*2
		}
		var3[,"Ntot"]=stress_function(var3[,"Ntot"],threshold)
		var37[,"Ntot"]=stress_function(var37[,"Ntot"],threshold)
	}

#Seasonal component from temperature data
	if (desaisonnalise){
		#First, compute the dedicated seasonal variable based on temperature variation, modeled with a trigonometric function
		t=1:length(var3[,'TEMP'])
		omega=2*pi/(365.25/14)
		sini=sin(omega*t)
		cosi=cos(omega*t)
		gap_temp=!is.na(var3[,'TEMP'])
		season=lm(var3[,'TEMP']~sini+cosi)
		#Seconde, use only the residuals of the covariate when projected on the seasonal component
		for (c in cov_saison){
			var3[gap_temp&!is.na(var3[,c]),c]=lm(var3[gap_temp,c]~season$fitted.values)$residuals
			var3[is.na(var3[,'TEMP']),c]=NA
		}

		#Same process for B7 values
		t=1:length(var37[,'TEMP'])
                omega=2*pi/(365.25/14)
                sini=sin(omega*t)
                cosi=cos(omega*t)
		season7=lm(var37[,'TEMP']~sini+cosi)
		gap_temp=!is.na(var37[,'TEMP'])
		for (c in cov_saison){
			var37[gap_temp&!is.na(var37[,c]),c]=lm(var37[gap_temp,c]~season7$fitted.values)$residuals
			var37[is.na(var37[,'TEMP']),c]=NA
		}
		
		#Concatenate values for Teychan and B7
		season_concat=c(season$fitted.values,season7$fitted.values)
		if(centr_reduit){
			season_concat=scale(season_concat)
		}
	
	}
        var3_concat=rbind(var3,var37)

#Scaling
	if(centr_reduit){
		var3_concat=scale(var3_concat)
	}

#We keep only complete observations	
	id_na=apply(cbind(tx_croiss_concat,mat_ab_concat,var3_concat),2,is.na)
	id_bis=rep(TRUE,dim(id_na)[1])
	for (j in 1:dim(id_na)[2]){
		id_bis=id_bis*!id_na[,j]
	}
	id_bis=as.logical(id_bis)
	tx_croiss_concat=tx_croiss_concat[id_bis]
	mat_ab_concat=mat_ab_concat[id_bis,]
	var3_concat=var3_concat[id_bis,]
	cov3bis=cov3

#Evaluate the linear model for the species s
	if(length(cov3)>0){
#First, build the set of covariate that will be used
		com=paste("var3_concat[,'",cov3bis[1],"']",sep="")
		for (v in 2:length(cov3bis)){
				com=paste(com,"+var3_concat[,'",cov3bis[v],"']",sep="")
		}
		if(desaisonnalise){
				season_concat=season_concat[id_bis]
				com=paste(com,"+season_concat",sep="")
		}
#Then, growth rates are modelled with past abundance and covariates
		eval(parse(text=paste("ll=lm(tx_croiss_concat~mat_ab_concat+",com,")",sep="")))
		print('********COLLINEARITY TEST***********')
		print(colldiag(ll),fuzz=0.4)
		max_ci=max(max_ci,max(colldiag(ll)$condindx))
		avg_ci=avg_ci+max(colldiag(ll)$condindx)/length(sp)
		print('********COLLINEARITY TEST END***********')
		a=anova(ll)
#Finally, keep tracks of the results with anova, AIC, AICc
		print(summary(ll))
		print(a)
		print("*****************************AIC**************************")
		print(extractAIC(ll))	
		print(AICc(ll))	
		print("*****************************AIC**************************")

		if(figure_4){
	                if(ss==1){
        	        	aadd=FALSE
                	        pdf("lm_tot_lag3.pdf",width=20)
                        	par(mfrow=c(1,1),oma=c(1,1,5,1),mar=c(5,5,1,1))
	                }else{
        	        	aadd=TRUE
                	}
                	errbar(ss,coef(summary(ll))[2+nb_lag, "Estimate"],coef(summary(ll))[2+nb_lag, "Estimate"]+coef(summary(ll))[2+nb_lag, "Std. Error"],coef(summary(ll))[2+nb_lag, "Estimate"]-coef(summary(ll))[2+nb_lag, "Std. Error"],xlim=c(1,52),add=aadd,ylim=c(-0.5,0.5),xaxt="n",ylab="Covariate effect",xlab="",lwd=2.5,cex.lab=2.5,cex.axis=1.75,cex=2)
                	for (c in 2:(length(cov3)+1)){ #+1 for season
		                errbar(ss+13*(c-1),coef(summary(ll))[1+nb_lag+c, "Estimate"],coef(summary(ll))[1+nb_lag+c, "Estimate"]+coef(summary(ll))[1+nb_lag+c, "Std. Error"],coef(summary(ll))[1+nb_lag+c, "Estimate"]-coef(summary(ll))[1+nb_lag+c, "Std. Error"],add=TRUE,lwd=2.5,cex=2)
                	}
		}


#Plot the residuals
		if(doyouplot){
	               	pdf(paste(s,"_verif_residuals.pdf",sep=""))
        	       	par(mfrow=c(3,2))
	               	plot(ll,cex.main=fac_main,cex.axis=fac_axis,cex.lab=fac_lab,lwd=alwd,cex=apc,ps=apc)
			acf(residuals(ll),cex.main=fac_main,cex.axis=fac_axis,cex.lab=fac_lab,lwd=alwd,cex=apc)
			dev.off()
		}

		}
}
		if(figure_4){
	               	axis(1,at=c(seq(1,12),seq(14,25),seq(27,38),seq(40,51)),rep(c(sp),4),las=2,cex.axis=1.75)
       		       	mtext(c("Irradiance","Wind","Salinity","Season"),side=3,at=c(6,19,33,47),cex=2.5,line=3)
      	   	       	mtext(c(expression(paste("(J cm"^"-2",")",sep="")),expression(paste("(m"^"2"," s"^"-2",")",sep="")),expression(paste("(g kg"^"-1",")",sep="")),""),at=c(6,19,33,47),side=3,cex=2.5,line=0.25)
                	abline(v=13,lty=2,lwd=1.75)
	                abline(v=26,lty=2,lwd=1.75)
        	        abline(v=39,lty=2,lwd=1.75)
        	        abline(h=0,lty=2,lwd=1.75)
                	dev.off()
		}

sink()

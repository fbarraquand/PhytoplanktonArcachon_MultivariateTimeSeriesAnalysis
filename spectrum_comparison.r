###############################################################################################################################################
                                             ########CP - Spectrum (mostly graphics) analyses ########
###############################################################################################################################################
rm(list=ls())
graphics.off()

library("zoo")
library("spectral.methods")
source("functions_global.r")

#Choice of phytoplankton groups and covariates we want to use
sp=c("GUI", "LEP", "NIT", "PSE", "RHI", "SKE", "CHA", "AST", "GYM",  "PRP", "EUG", "CRY")
sp=sort(sp)
cov3=c("CumDebit","MeanNAO","Ntot","CumPrec","CumRg","SAL","TEMP","MeanVent","PHOS","SI")
cov_comparable=c("Ntot","SAL","TEMP","PHOS","SI")
cov3=sort(cov3)

#Graphics parameters
fac_main=3.0
fac_axis=2.0
fac_lg=2.0
fac_lab=2.0
alwd=2.5
ate=seq(3,12*12,3)
labe=ate
labe[seq(1,length(ate),2)]=NA
xli1=8
xli2=72

#Other parameters
maxi_test=10 #Maximum number of tests to perform on time series reconstruction
accum_vent=3 #Number of days on which we cumulate wind and irradiance
m=3 #Length of smoothing window

#Loading data
filename="Teychan_base.csv"
tab=read.csv(filename,na.strings="NA",header=TRUE,sep=";",dec=".")
dates=as.Date(tab$Date)
dates_bis=seq(dates[dates>as.Date("01/01/1997","%d/%m/%Y")][1],dates[length(dates)],14)
tab_cov=new_covar(tab,dates,accum_vent)

filename2="B7_base.csv"
tab7=read.csv(filename2,na.strings="NA",header=TRUE,sep=";",dec=".")
dates7=as.Date(tab7$Date)
dates_bis7=seq(dates7[1],dates7[length(dates7)],14)
tab_cov7=new_covar(tab7,dates7,accum_vent)

pdf("spectrum_coherence_clean.pdf",width=20)
par(mfrow=c(3,4),mar=c(4.1, 5.1, 4.1, 1.1))
        for(c in 1:length(cov3)){
                covar1=approx(tab_cov[,cov3[c]],x=dates,xout=dates_bis)$y
                id=!is.na(covar1)
		covar1_bis=covar1[id]
                if(sum(is.na(covar1))>0){ #if there are still missing values, they are at the end or at the beginning
                        info=sum(is.na(covar1))
			dd=dates_bis[id]
                }else{
                        info=0
			dd=dates_bis
                }
		#Periodogram
                x_cov=spec.pgram(covar1_bis,kernel("modified.daniell",m),taper=0,plot=FALSE)

		#Compute periods in months instead of frequency
		freq_cov=1/x_cov$freq*14/(365.25/12)
                spec_cov=x_cov$spec

		#Same analysis for B7
		covar17=approx(tab_cov7[,cov3[c]],x=dates7,xout=dates_bis7)$y
                idi7=!is.na(covar17)
		covar17_bis=covar17[idi7]
                if(sum(is.na(covar17))>0){
                        info=sum(is.na(covar17))
			dd7=dates_bis7[idi7]
                }else{
                        info=0
			dd7=dates_bis7
                }
                x_cov=spec.pgram(covar17_bis,kernel("modified.daniell",m),taper=0,plot=FALSE)

		freq_cov7=1/x_cov$freq*14/(365.25/12)
                spec_cov7=x_cov$spec
		
		#Plot periodograms
		yli2=max(c(spec_cov,spec_cov7))
                plot(freq_cov,spec_cov,"l",xaxt="n",xlim=c(xli1,xli2),xlab="Period (month)",main=cov3[c],lty=1,lwd=alwd,ylim=c(0,yli2),cex.axis=fac_axis,cex.lab=fac_lab,cex.main=fac_main,ylab="")
                lines(freq_cov7,spec_cov7,lty=1,lwd=alwd,col="blue")
        	 axis(1,at=ate,labels=labe,cex.lab=fac_lab)
		

		#Compute coherence
		 for(cbis in 1:length(cov3)){
			if (cbis!=c){
                #For Teychan
		covar2=approx(tab_cov[,cov3[cbis]],x=dates,xout=dates_bis)$y
                idbis=!is.na(covar2)
                        covar2_bis=covar2[id&idbis]
			covar1_ter=covar1[id&idbis]

                x_coherency=spec.pgram(cbind(covar1_ter,covar2_bis),kernel("modified.daniell",c(m,m)),taper=0,plot=FALSE)

		#Compute significance thresholds
                L=2*m+1 #L=2m+1
                f = qf(.995, 2, x_coherency$df-2) #alpha=5%, corrected Bonferroni because 10 variables
                C = f/(L+f)

		
                freq=1/x_coherency$freq*14/(365.25/12)
                spec=x_coherency$coh

		#For B7
                covar27=approx(tab_cov7[,cov3[cbis]],x=dates7,xout=dates_bis7)$y
                idbis7=!is.na(covar27)
                        covar2_bis7=covar27[idi7&idbis7]
                        covar1_ter7=covar17[idi7&idbis7]

                x_coherency7=spec.pgram(cbind(covar1_ter7,covar2_bis7),kernel("modified.daniell",c(m,m)),taper=0,plot=FALSE)

                freq7=1/x_coherency7$freq*14/(365.25/12)
                spec7=x_coherency7$coh


		#Plot
                yli2=max(C,max(spec),max(spec7))
                plot(freq,spec,"l",xlim=c(xli1,xli2),xaxt="n",xlab="Period (month)",main=paste("Coherence",cov3[cbis],cov3[c]),lty=1,ylim=c(0,max(C,max(spec))),cex.main=fac_main,lwd=alwd,cex.lab=fac_lab,cex.axis=fac_axis,ylab="")
                abline(h = C,col="black",lty=2,lwd=alwd)
                axis(1,at=ate,labels=labe,cex.axis=fac_axis,cex.lab=fac_lab)
                lines(freq7,spec7,"l",lty=1,col="blue",lwd=alwd)
        }
}

                plot(0,0,t="n",xaxt="n",yaxt="n",ylab="",xlab='',bty="n")
                plot(0,0,t="n",xaxt="n",yaxt="n",ylab="",xlab='',bty="n")
                legend('topleft',c('Teychan','Buoy7','Spectrum','0.5% Significance'),lty=c(1,1,1,2),col=c('black',"blue","black","black"),bty="n",cex=fac_lg)
}
dev.off()

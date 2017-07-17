
###############################################################################################################################################
########## CP           - Global functions used in the analyses ###############
###############################################################################################################################################

#Cumulate daily values, transform values to ratios when necessary, cumulate nitrogen sources
new_covar=function(tab,dates,accum_vent,start="t-1"){
#Load meteorological data
        path_data_post="./"
	load(paste(path_data_post,"vent.RData",sep=""))
	load(paste(path_data_post,"meteo.RData",sep=""))
	date_meteo=as.Date(meteo$Date)
	load(paste(path_data_post,"nao.RData",sep=""))
	
	cov_bis=c("CHL","CumDebit","MeanNAO","Ntot","PHEO","P/N","CumPrec","CumRg","SAL","SI/N","TEMP","MeanVent","PHOS","SI","MES","NH4","NOx","Rg","NAO_month","AMO_month")
	tab_cov=matrix(NA,nrow=dim(tab)[1],ncol=length(cov_bis),dimnames=list(rownames(tab),cov_bis))
#No need to change theses values
	tab_cov[,"Rg"]=tab[,"Rg"]
	tab_cov[,"NAO_month"]=tab[,"NAO_month"]
	tab_cov[,"AMO_month"]=tab[,"AMO_month"]
	tab_cov[,"CHL"]=tab[,"CHL"]
	tab_cov[,"MES"]=tab[,"MES"]
	tab_cov[,"NOx"]=tab[,"NOx"]
	tab_cov[,"NH4"]=tab[,"NH4"]
	tab_cov[,"SI"]=tab[,"SI"]
	tab_cov[,"PHOS"]=tab[,"PHOS"]
	tab_cov[,"Ntot"]=tab[,"NH4"]+tab[,"NOx"]
	tab_cov[,"PHEO"]=tab[,"PHEO"]
	tab_cov[,"P/N"]=tab[,"PHOS"]/tab_cov[,"Ntot"]
	id=is.infinite(tab_cov[,"P/N"])
	tab_cov[id,"P/N"]=NA
	tab_cov[,"SAL"]=tab[,"SAL"]
	tab_cov[,"SI/N"]=tab[,"SI"]/tab_cov[,"Ntot"]
	id=is.infinite(tab_cov[,"SI/N"])
	tab_cov[id,"SI/N"]=NA
	tab_cov[,"TEMP"]=tab[,"TEMP"]

#Start defines the period over which meteorological variables are defined: it can be from t-1 to t in the case of correlation between variables ; or it can be from t to t+1 when doing the ARLM analysis because there's a shift in data afterwards
	if(start=="t-1"){
		sta=2
		sto=length(dates)
	}else if(start=="t"){
		sta=1
		sto=length(dates)-1
	}else{
		stop("You should define the period over which you wish to integrate your physical data")
	}

	for (d in sta:sto){
		if(start=="t-1"){
        		d1=dates[d-1]
	        	d2=dates[d]
		}else if(start=="t"){
        		d1=dates[d]
	        	d2=dates[d+1]
		}
#Mean value for NAO
	        tab_cov[d,"MeanNAO"]=mean(nao$Val[nao$Date>=d1&nao$Date<d2],na.rm=TRUE)
#Cumulated precipitation, irradiance and inflow
                tab_cov[d,"CumPrec"]=tail(cumsum(meteo$Precipitation[date_meteo>=d1&date_meteo<d2]),1)
                tab_cov[d,"CumRg"]=tail(cumsum(meteo$Rg[date_meteo>=d1&date_meteo<d2]),1)
                tab_cov[d,"CumDebit"]=tail(cumsum(meteo$Debit.Eyre[date_meteo>=d1&date_meteo<d2]),1)
#Cumulated wind energy
                tab_vent=rep(NA,d2-d1)
                for (dd in 1:(d2-d1)){
                        d_ap=d1+dd
                        d_av=d_ap-accum_vent
                        id=which(vtot$date>=d_av&vtot$date<=d_ap)
                        tab_vent[dd]=sum(as.numeric(vtot$vitesse[id])^2)
                }
                tab_cov[d,"MeanVent"]=mean(tab_vent,na.rm=TRUE)
}
        return(tab_cov)
}

########## CP           - Function taking into account saturation, based on half-saturation threshold for a nutrient ###############
stress_function=function(nit,threshold){
        lambda=log(2)/threshold
        vec=1-exp(-lambda*nit)
        return(vec)
        }

########## CP           - Function embedding MARSS model to make sure that consistent options are used throughout the analyses ###############

analyse_MARSS=function(tab_sp,tab_cov,B1,filename,C1="unconstrained",amethod="kem",pi1="zero",V1="identity",aalpha=0.05,min_iter=15,max_iter=500,iter_boot=10,aconv.test.slop.tol=0.001,aabstol=0.001,boot=TRUE,Q1="diagonal and unequal",aicb=FALSE){
#Matrix        
	U1="zero" #Offset for processes
        Z1=diag(1,length(sp),length(sp)) #Identity between Y and X
        A1="zero" #Offset for observations
        R1="zero" #Covariance matrix for the observation matrix
	c1=tab_cov

        if(amethod=="kem"){
		cntl.list=list(conv.test.slope.tol=aconv.test.slop.tol,minit=min_iter,maxit=max_iter,abstol=aabstol)
	}else{
		cntl.list=list()
	}
	model.list=list(B=B1,U=U1,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=1,C=C1,c=c1)

	print(paste("Start estimation for",filename,sep=" "))
        fit_log=MARSS(tab_sp, method=amethod,model=model.list,control=cntl.list)
	print("Stop estimation")

#Boostrapping
	if(boot){
	        print("Start bootstrap")
        	cis=MARSSparamCIs(fit_log,method="parametric",nboot=iter_boot,alpha=aalpha)
        	save(fit_log,cis,file=filename)
		}
	else{
        	save(fit_log,file=filename)
	}

#AIC computation
	if(aicb){
	        print("Start AIC")
        	aic=MARSSaic(fit_log,output=c("AIC","AICc", "AICbp", "boot.params"),Options=list(nboot=iter_boot))
        	save(fit_log,aic,file=paste("aic_",filename,sep=""))

	}

        print(paste("Sauvegarde de",filename,sep=" "))
}



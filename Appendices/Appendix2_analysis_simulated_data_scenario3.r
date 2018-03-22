rm(list=ls())

##Setting path and library
#path_script="/home/cpicoche/Documents/Plankton/script/"
path_data_post="/home/frederic/Documents/Plankton/SimulatedData/Loop_MARSS/"

graphics.off()
library("zoo")
library("lubridate")
library("MARSS")

###Phytoplankton, guild definitions
sp=c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10")
iter_min=100 #200 for Griffiths ; but it seems that about 30 are enough in our case
iter_estimate=10 #not used for now
iter_scenario=10
iter_boot=1000 #Ok for Griffiths
scenario=as.character(3)
BFGS=FALSE

###Loading data
filename=paste(path_data_post,"MockPlanktonTimeSeries_Loop.csv",sep="")
tabbis=read.csv(filename,na.strings="NA",header=TRUE,sep=",",dec=".")
for (r in 1:iter_scenario){
	cov3_tot=c("Abiotic_var1","Abiotic_var2")
	tab=tabbis[tabbis$Scenario==as.character(scenario)&tabbis$Repeat==r,]
	dates=tab$Time_index
	dates_bis=dates #already interpolated
	tab_sp=tab[,sp]
	tab_cov=tab[,cov3_tot]
	tab_cov_bis=tab_cov #to keep coherence with the real data set

	#Setting model
	B1="unconstrained"
	U1="zero"
	Q1="diagonal and unequal"
	Z1=diag(1,length(sp),length(sp))
	A1="zero"
	R1="zero"
	V1=diag(1,length(sp))
	pi1="zero"
	C1="unconstrained"
	aalpha=0.05
	cntl.list=list(conv.test.slope.tol=0.001,minit=iter_min,maxit=500,abstol=0.001)#,MCInit=TRUE,silent=2)#,numInits=iter_estimate,numInitSteps=10) #Took values from Griffiths for minit,maxit,abstol. conv.test.slope is recommended in MARSS User Guide itself, numInits and numInitStep are taken as rules of thumbs (quick search on the Internet, numInits=500 elsewhere but it seems really big to me

	#Without season and with EM algorithm
	tab_sp=t(scale(tab_sp))
	tab_cov=t(scale(tab_cov_bis))
	rownames(tab_sp)=sp
	rownames(tab_cov)=cov3_tot
	c1=tab_cov
	model.list=list(B=B1,U=U1,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=1,C=C1,c=c1)
	fit_log=MARSS(tab_sp, method="kem",model=model.list,control=cntl.list)

	cis=MARSSparamCIs(fit_log,method="parametric",nboot=iter_boot,alpha=aalpha)
	save(fit_log,cis,file=paste("MockSeries_scenario",scenario,"_rep",r,"_kem_noseason.RData",sep=""))
	
	#verif BFGS
	if(BFGS){
        fit_log=MARSS(tab_sp, method="BFGS",model=model.list)
        cis=MARSSparamCIs(fit_log,method="parametric",nboot=iter_boot,alpha=aalpha)
        save(fit_log,cis,file=paste("MockSeries_scenario",scenario,"_rep",r,"_bfgs_noseason.RData",sep=""))
	}

	#With season with EM algorithm
	t=1:length(tab_cov_bis[,'Abiotic_var1'])
	omega=2*pi/(24) #Took frequency from Fred's data
	sini=sin(omega*t)
	cosi=cos(omega*t)
	season=lm(tab_cov_bis[,'Abiotic_var1']~sini+cosi)
	for (c in cov3_tot){
		tab_cov_bis[,c]=lm(tab_cov_bis[,c]~season$fitted.values)$residuals
        }
	tab_cov_bis=cbind(tab_cov_bis,season$fitted.values)
	cov3_tot=c(cov3_tot[1:(length(cov3_tot))],"season")
	tab_cov=t(scale(tab_cov_bis))
        rownames(tab_cov)=cov3_tot
        c1=tab_cov
        
        model.list=list(B=B1,U=U1,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=1,C=C1,c=c1)
        fit_log=MARSS(tab_sp, method="kem",model=model.list,control=cntl.list)
        cis=MARSSparamCIs(fit_log,method="parametric",nboot=iter_boot,alpha=aalpha)
        save(fit_log,cis,file=paste("MockSeries_scenario",scenario,"_rep",r,"_kem_season.RData",sep=""))

	if(BFGS){
        #verif BFGS
        fit_log=MARSS(tab_sp, method="BFGS",model=model.list)
        cis=MARSSparamCIs(fit_log,method="parametric",nboot=iter_boot,alpha=aalpha)
        save(fit_log,cis,file=paste("MockSeries_scenario",scenario,"_rep",r,"_bfgs_season.RData",sep=""))
	}
}

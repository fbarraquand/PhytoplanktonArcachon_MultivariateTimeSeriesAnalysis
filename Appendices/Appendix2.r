rm(list=ls())
graphics.off()

library("zoo")
library('lubridate')
source("../functions_global.r")

###########Figure S2.1
#Phytoplankton
sp=c("GUI", "LEP", "NIT", "PSE", "RHI", "SKE", "CHA", "AST", "GYM",  "PRP", "EUG", "CRY")
sp=sort(sp)
#Data
filename="../Teychan_base.csv"
tabbis=read.csv(filename,na.strings="NA",header=TRUE,sep=";",dec=".")
dates=as.Date(tabbis$Date)
#New
filename7="../B7_base.csv"
tabbis7=read.csv(filename7,na.strings="NA",header=TRUE,sep=";",dec=".")
dates7=as.Date(tabbis7$Date)

fac_main=3.7
fac_axis=2.5
fac_lg=2.5
fac_pch=2.5
fac_lab=4
alwd=2.5
awidth=20

pdf("growth_vs_abundance_ricker.pdf",height=20,width=20)
par(mfrow=c(4,3),mar=c(5.1, 6.1, 4.1, 1.1))
i=0
for (s in sp){
        i=i+1
        #CRY and EUG have not been correctly recorded for a while
        if(s=="CRY"|s=="EUG"|s=="GYM"){
                dates_correct=dates[year(dates)>1996]
                tab=tabbis[year(dates)>1996,]

                dates_correct7=dates7[year(dates7)>1996]
                tab7=tabbis7[year(dates7)>1996,]
        }else{
                dates_correct=dates
                tab=tabbis

                dates_correct7=dates7
                tab7=tabbis7
        }
        dates_bis=seq(dates_correct[1],dates_correct[length(dates_correct)],14)
        dates_bis7=seq(dates_correct7[1],dates_correct7[length(dates_correct7)],14)

        #Dealing with NA and zero values for plankton
	#Gompertz
        var1=log(na.approx(tab[,s],maxgap=2,x=dates_correct,xout=dates_bis,na.rm=FALSE))
        var17=log(na.approx(tab7[,s],maxgap=2,x=dates_correct7,xout=dates_bis7,na.rm=FALSE))
        tx_croiss=diff(var1)
        tx_croiss7=diff(var17)
	#Ricker
        var1=na.approx(tab[,s],maxgap=2,x=dates_correct,xout=dates_bis,na.rm=FALSE)
        var17=na.approx(tab7[,s],maxgap=2,x=dates_correct7,xout=dates_bis7,na.rm=FALSE)
        var1=var1[1:length(tx_croiss)]
        var17=var17[1:length(tx_croiss7)]
        if(i%%3==1){
                yylab="Growth rates (C/L)"
        }else{
                yylab=""
        }
        if(i>9){
	#Gompertz
#                xxlab="Log abundance (C/L)"
	#Ricker
                xxlab="abundance (C/L)"
        }else{
                xxlab=""
        }

        plot(var1,tx_croiss,t="p",col="black",pch=16,xlab=xxlab,ylab=yylab,cex.main=fac_main,cex.axis=fac_axis,cex.lab=fac_lab,main=s,cex=fac_pch)
        points(var17,tx_croiss7,pch=16,col="blue",cex=fac_pch)
        legend("topright",c("Teychan","Buoy 7"),col=c("black","blue"),pch=16,cex=fac_lg)
}
dev.off()

###########Figure S2.2 is included in ARLM_fixed_nb_lag_clean.r
###########Figure S2.3 and S2.4 are in Appendices2_plot_for_simulated_data.r

###########Figure S2.5
library("MARSS")

args = "../MARSS_results/B7_physics_unconstrained.RData"

load(args[1])
fit_log=fit_log
a=residuals(fit_log)

fac_main=4.0
fac_axis=2.5
fac_lg=3.5
fac_lab=4
alwd=5
apc=3
awidth=15

sp=c("AST","NIT","PSE","SKE","CHA","GUI","LEP","RHI","GYM","PRP","CRY","EUG")

graphename=paste(strsplit(args[1],".RData")[[1]][1],"_residuals_essai.pdf",sep="")
pdf(graphename,width=20)
par(mfrow=c(1,3),mar=c(7,6,3,1))
for (i in 1:length(sp)){
        plot(fit_log$states[i,],a$residuals[12+i,],cex.main=fac_main,cex.axis=fac_axis,cex.lab=fac_lab,lwd=alwd,cex=apc,ps=apc,main="",ylab="Residuals",xlab="Fitted values")
        abline(h=0,cex.main=fac_main,cex.axis=fac_axis,cex.lab=fac_lab,cex=apc,ps=apc,lty=2,lwd=alwd)
        qqnorm(a$residuals[12+i,],cex.main=fac_main,cex.axis=fac_axis,cex.lab=fac_lab,lwd=alwd,cex=apc,ps=apc,main=sp[i])
        qqline(a$residuals[12+i,],lwd=alwd)
        acf(a$residuals[12+i,],cex.main=fac_main,cex.axis=fac_axis,cex.lab=fac_lab,lwd=alwd,cex=apc,ps=apc,main="")
}
dev.off()


###########Figure S2.6 and S2.7 are the same as figure 6

###########Figure S2.8
library("stringr")
l=c("Teychan")
model=c("pencen")

for(lieu in l){
for(m in model){
rm(list=ls()[grep('fit_log',ls())])
openfile=paste("../MARSS_results/",lieu,"_physics_",m,".RData",sep="")
load(openfile)
eval(parse(text=paste("fit_log=",ls()[grep('fit_log',ls())],sep="")))

filenamebis=paste(strsplit(openfile,".RData")[[1]][1],"_cross_validation.pdf",sep="")

###Phytoplankton, guild definitions
sp=c("AST","NIT","PSE","SKE","CHA","GUI","LEP","RHI","GYM","PRP","CRY","EUG")
###Nutrient definition
cov3_phy=c("SAL","CumRg","MeanVent")
cov3_nut=c('Ntot','SI','PHOS')
cov3_tot=c(cov3_phy,"TEMP") #covariates we are going to use, including temperature

###Loading data
filename=paste("../",lieu,"_base.csv",sep="")
tabbis=read.csv(filename,na.strings="NA",header=TRUE,sep=";",dec=".")
dates=as.Date(tabbis$Date)

tab_sp=tabbis[,sp]
###Treating covariates to have cumulated radiation and wind energy
accum_vent=3 #We take into account wind energy 3 days before the date
tab_cov=new_covar(tabbis,dates,accum_vent)

###Treating missing values for both species and covariates###
consecutif=2 #Number of missing values above which we keep the NA
timestep=14 #Regular time lapses between two observations
dates_bis=seq(dates[1],dates[length(dates)],timestep) #Regular time grid
tab_sp=na.approx(tab_sp,maxgap=consecutif,x=dates,xout=dates_bis,na.rm=FALSE) #Interpolation over regular time grid
tab_cov_bis=matrix(NA,length(dates_bis),length(cov3_tot))
colnames(tab_cov_bis)=cov3_tot
for (c in cov3_tot){
        tab_cov_bis[,c]=approx(tab_cov[,c],x=dates,xout=dates_bis)$y
}

#Difference between estimation and validation
tab_sp=log(tab_sp)
tab_sp=t(scale(tab_sp[1:(length(dates_bis)),]))

tab_sp_estimation=tab_sp[,year(dates_bis[1:(length(dates_bis))])>1996]#Using data from 1997, because cryptophytes were not counted (or badly, for the first year, before
dates_bis_estimation=dates[year(dates_bis[1:(length(dates_bis))])>1996]
tab_sp_validation=tab_sp[,year(dates_bis[1:(length(dates_bis))])<=1996]#Using data from 1997, because cryptophytes were not counted (or badly, for the first year, before
dates_bis_validation=dates[year(dates_bis[1:(length(dates_bis))])<=1996]

tab_sp_validation["CRY",]=rep(0,dim(tab_sp_validation)[2])
tab_sp_validation["EUG",]=rep(0,dim(tab_sp_validation)[2])
rownames(tab_sp_estimation)=sp
rownames(tab_sp_validation)=sp

###Seasonal component
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
tab_cov_bis=t(scale(tab_cov_bis[1:(length(dates_bis)),]))

tab_cov_bis_estimation=tab_cov_bis[,year(dates_bis[1:(length(dates_bis))])>1996]
tab_cov_bis_validation=tab_cov_bis[,year(dates_bis[1:(length(dates_bis))])<=1996]
rownames(tab_cov_bis_estimation)=cov3_tot
rownames(tab_cov_bis_validation)=cov3_tot
predictee=tab_sp_validation[,3:dim(tab_sp_validation)[2]]
predicted=matrix(NA,nrow=12,ncol=dim(tab_sp_validation)[2]-2)
U=matrix(fit_log$par$U,nrow=12,ncol=4)
B=matrix(0,nrow=12,ncol=12)
nom=dimnames(fit_log$par$B)[[1]]
for (n in 1:length(nom)){
        if(grepl("^\\(",nom[n])){ #default when using unconstrained or diagonal
                i=as.numeric(strsplit(nom[n],split="[\\(,\\)]")[[1]][2])
                j=as.numeric(strsplit(nom[n],split="[\\(,\\)]")[[1]][3])
                if(is.na(i)){ #i and j might not be numeric
                        i=which(unlist(lapply(sp,function(x) grepl(x,strsplit(nom[n],split="[\\(,\\)]")[[1]][2]))))
                        j=which(unlist(lapply(sp,function(x) grepl(x,strsplit(nom[n],split="[\\(,\\)]")[[1]][3]))))
                }
        }else{
                a=str_split(nom[n],sp)
                for (abis in 1:length(a)){
                        if(length(a[[abis]])==2){
                                if((a[[abis]][1]=="")&&(a[[abis]][2]!="")){
                                        j=abis
                                }
                                if((a[[abis]][2]=="")&&(a[[abis]][1]!="")){
                                        i=abis
                                }
                        }else if(length(a[[abis]])==3){
                                if((a[[abis]][2]=="")&&(a[[abis]][1]=="")&&(a[[abis]][3]=="")){
                                        j=abis
                                        i=abis
                                }
                        }
                }
        }
        B[i,j]=fit_log$par$B[n]
}


for (data in 3:dim(tab_sp_validation)[2]){
                obstm1=tab_sp_validation[,data-1]
                covt=tab_cov_bis_validation[,data]
                obst=tab_sp_validation[,data]

                sim=rep(NA,12)
                for(i in 1:12){
                        ss=0
                        for (j in 1:12){
                                if(is.na(obstm1[j])){
                                        if(B[i,j]!=0){
                                                ss=ss+NA
                                        }else{
                                                ss=ss+0
                                        }
                                }else{
                                        ss=ss+B[i,j]*obstm1[j]
                                }
                        }
                        sim[i]=ss
                }
                sim=sim+U%*%covt

                predicted[,data-2]=sim
}

r1=0
r2=0
for (i in 1:10){ #not 12, because we don't take into account CRY and EUG
        for (j in 1:dim(predicted)[2]){
                if(!is.na(predicted[i,j]*predictee[i,j])){
                        r1=r1+predictee[i,j]^2
                        r2=r2+(predicted[i,j]-predictee[i,j])^2
                }
        }
}
r=1-r2/r1

print(paste(m,lieu))
print(r)

fac_main=3.0
fac_axis=2.5
fac_lg=2.5
fac_lab=2.5
alwd=3.5
awidth=20
dim_fin=dim(tab_sp_validation)[2]
dates_bis_validation=dates_bis_validation[3:dim_fin]
pdf(filenamebis,width=awidth)
par(mar=c(3,5,3,0.5))
for(i in 1:10){
        plot(dates_bis_validation,predictee[i,],pch=16,col="red",main=sp[i],ylab="Scaled log abundance",xlab="",cex.main=fac_main,cex.axis=fac_axis,cex.lab=fac_lab,cex=2)
        points(dates_bis_validation,predicted[i,],pch=16,col="black",cex=1.5)
        for(t in 2:(dim_fin-3)){
                arrows(dates_bis_validation[t-1],predictee[i,t-1],dates_bis_validation[t],predicted[i,t],length=0.15,lwd=alwd,col="black")
                lines(rep(dates_bis_validation[t],2),c(predicted[i,t],predictee[i,t]),col="blue",lwd=alwd+0.5,lty=2)
        }
        legend("topleft",c("Observed","Predicted","Prediction","Distance"),col=c("red","black","black","blue"),lty=c(NA,NA,1,2),pch=c(16,16,NA,NA),cex=fac_lg,lwd=alwd,bty="n")
}
dev.off()
}
}


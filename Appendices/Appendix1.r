rm(list=ls())
graphics.off()
library("lubridate")
###########Figure S1.1
#Loading data
filename="../Teychan_base.csv"
filenameB7="../B7_base.csv"
tabbis=read.csv(filename,na.strings="NA",header=TRUE,sep=";",dec=".")
tabbisB7=read.csv(filenameB7,na.strings="NA",header=TRUE,sep=";",dec=".")
tabbis$Date=as.Date(tabbis$Date)
tabbisB7$Date=as.date(tabbisB7$Date)
findnei=which(colnames(tabbis)=="NEI")

#Plotting data
pdf("boxplot_var.pdf")
var=c("TEMP","SAL","CHL")
var_jolie=c("Temperature (°C)",expression(paste("Salinity (g kg"^'-1',")",sep="")),expression(paste("Chlorophyll a (",mu,"g L"^"-1",")",sep="")),"TURB","TURBFNU")
par(mar=c(5,6,2,2))
i=0
for (v in var){
i=i+1
yli2=max(tabbis[,v],na.rm=TRUE)+0.5
yli1=min(tabbis[,v],na.rm=TRUE)-0.5
plot(0:13,rep(0,14),type="n",ylim=c(yli1,yli2),xlab="Month",ylab=var_jolie[i],xaxt="n",yaxt="n",cex.axis=2,cex.lab=2)
for (m in 1:12){
        var=tabbis[month(tabbis$Date)==m,v]
        var7=tabbisB7[month(tabbisB7$Date)==m,v]
        boxplot(var,at=m,add=TRUE,col="grey",bg="grey",pch=21,cex=2,cex.axis=1.5,range=1.5,lwd=1.5)
        boxplot(var7,at=m+0.4,add=TRUE,col="blue",pch=21,cex=2,bg="blue",yaxt="n",range=1.5,lwd=1.5)
        print(paste("Month",m,":",v,min(var7,na.rm=TRUE),mean(var7,na.rm=TRUE),max(var7,na.rm=TRUE)))
}
axis(1,at=1:12,labels=c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'),cex.axis=2)
}
dev.off()

###########Figure S1.2
pdf('stratif.pdf',width=20,height=15)
par(mfrow=c(2,1),mar=c(5.1,5.5,3,0.1))
lieu=c("Teychan","B7")
for (l in 1:length(lieu)){
###Loading data
filename=paste("../",lieu[l],"_base.csv",sep="")
tabbis=read.csv(filename,na.strings="NA",header=TRUE,sep=";",dec=".")
dates=as.Date(tabbis$Date)

plot(dates,tabbis[,"STRATIF_SAL"],col="black",lty=1,t="o",ylim=c(-2.5,5.5),xlab="",ylab="Difference bottom-surface",xlim=c(as.Date("2003-01-01"),as.Date("2016-01-01")),main=lieu[l],lwd=3,cex.axis=2.5,cex.main=4,cex.lab=2.5,pch=16,cex=1.5)
lines(dates,tabbis[,"STRATIF_TEMP"],col="blue",lty=1,lwd=3,t='o',pch=16,cex=1.5)
legend("topleft",c("Salinity","Temperature"),lty=c(1,1),col=c("black","blue"),bty="n",cex=3,lwd=3)
}
dev.off()


###########Figure S1.3
###This figure is based on the raw database we cleaned in Teychan  and B7 bases. This first database, and corresponding script, is not used afterwards and was not put in the repository. It is still available upon request.


###########Figure S1.4
library("corrplot")
library("zoo")
source("../functions_global.r")

#Data
filename="../Teychan_base.csv"
tabbis=read.csv(filename,na.strings="NA",header=TRUE,sep=";",dec=".")
dates=as.Date(tabbis$Date)
filename7="../B7_base.csv"
tabbis7=read.csv(filename7,na.strings="NA",header=TRUE,sep=";",dec=".")
dates7=as.Date(tabbis7$Date)

accum_vent=3 #We take into account wind energy 3 days before the date
tab_cov=new_covar(tabbis,dates,accum_vent,start="t-1")
tab_cov7=new_covar(tabbis7,dates7,accum_vent,start="t-1")

dates_bis=seq(dates[1],dates[length(dates)],14)
dates_bis7=seq(dates7[1],dates7[length(dates7)],14)

var2=na.approx(tab_cov,maxgap=0,x=dates,xout=dates_bis,na.rm=FALSE)
var27=na.approx(tab_cov7,maxgap=0,x=dates7,xout=dates_bis7,na.rm=FALSE)
var3_concat=rbind(var2,var27)

pdf("ARLM_corrplot.pdf",width=13,height=13)
par(mfrow=c(1,1))
cocov=colnames(var3_concat)
cocov=sort(cocov)
id=grep("CHL",cocov)
cocov=cocov[-id]
id=grep("PHEO",cocov)
cocov=cocov[-id]
M=cor(var3_concat[,cocov],method="spearman",use="pairwise.complete.obs")
corrplot(M,method="number")
dev.off()


###########Figure S1.5
fac_lab=1.0
fac_axis=1.0
fac_leg=1.0
fac_main=1.0
alwd=1.0

pdf("SiNtemp.pdf",width=20)
par(mfrow=c(1,1))
plot(dates_bis,var2[,'Ntot'],t="l",col="black",xlab="",ylab=expression(paste("Nitrogen/Silicon (",mu,"M/l)")),lwd=alwd,cex.main=fac_main,cex.axis=fac_axis,cex.lab=fac_lab)
points(dates_bis,var2[,'SI'],t="p",col="black",pch=16,lwd=alwd)
lines(dates_bis7,var27[,'Ntot'],t="l",col="blue",lwd=alwd)
points(dates_bis7,var27[,'SI'],t="p",col="blue",pch=16,lwd=alwd)
legend("topright",c("Teychan","Buoy 7","Nitrogen","Silicon"),col=c("black","blue","black","black"),pch=c(16,16,NA,16),lty=c(1,1,1,NA))
dev.off()

###########Figure S1.6
pdf("SivsN.pdf",width=13)
par(mfrow=c(1,1))
plot(var2[,'Ntot'],var2[,'SI'],t="p",col="black",pch=16,xlab=expression(paste("Total Nitrogen (",mu,"M/l)")),ylab=expression(paste("Silicon (",mu,"M/l)")),cex.main=fac_main,cex.axis=fac_axis,cex.lab=fac_lab,xlim=c(0,45))
points(var27[,'Ntot'],var27[,'SI'],t="p",col="blue",pch=16)
legend("topright",c("Teychan","Buoy 7"),col=c("black","blue"),pch=16)
dev.off()

###########Figure S1.7
library("spectral.methods")
cov3=c("CumDebit","MeanNAO","Ntot","CumPrec","CumRg","SAL","TEMP","MeanVent","PHOS","SI")
cov_comparable=c("Ntot","SAL","TEMP","PHOS","SI")
cov3=sort(cov3)

fac_main=3.0
fac_axis=2.0
fac_lg=2.0
fac_lab=2.0
alwd=2.5

maxi_test=10
accum_vent=3
ate=seq(3,12*12,3)
labe=ate
labe[seq(1,length(ate),2)]=NA
xli1=8
xli2=72
m=3

filename="../Teychan_base.csv"
tab=read.csv(filename,na.strings="NA",header=TRUE,sep=";",dec=".")
dates=as.Date(tab$Date)
dates_bis=seq(dates[dates>as.Date("01/01/1997","%d/%m/%Y")][1],dates[length(dates)],14)
tab_cov=new_covar(tab,dates,accum_vent)

filename2="../B7_base.csv"
tab7=read.csv(filename2,na.strings="NA",header=TRUE,sep=";",dec=".")
dates7=as.Date(tab7$Date)
dates_bis7=seq(dates7[1],dates7[length(dates7)],14)
tab_cov7=new_covar(tab7,dates7,accum_vent)

pdf("spectrum_and_acf.pdf",,width=20)
par(mfrow=c(3,4),mar=c(4.1, 5.1, 4.1, 1.1))
        for(c in 1:length(cov3)){
#Treating Teychan variables
                covar1=approx(tab_cov[,cov3[c]],x=dates,xout=dates_bis)$y
                if(sum(is.na(covar1))>0){ #removing NA still there, at the beginning or the end
                        id=!is.na(covar1)
                        info=sum(is.na(covar1))
                        covar1=covar1[id]
                        dd=dates_bis[id]
                }
                else{
                        info=0
                        dd=dates_bis
                }
                x_cov=spec.pgram(covar1,kernel("modified.daniell",m),taper=0,plot=FALSE)
                freq_cov=1/x_cov$freq*14/(365.25/12)
                spec_cov=x_cov$spec

#Treating B7 variables
                covar17=approx(tab_cov7[,cov3[c]],x=dates7,xout=dates_bis7)$y
                if(sum(is.na(covar17))>0){ #s'il reste encore des NA, cela signifie qu'ils sont au début ou à la fin
                        id=!is.na(covar17)
                        info=sum(is.na(covar17))
                        covar17=covar1[id]
                        dd7=dates_bis7[id]
                }
                else{
                        info=0
                        dd7=dates_bis7
                }
                x_cov=spec.pgram(covar17,kernel("modified.daniell",m),taper=0,plot=FALSE)
                freq_cov7=1/x_cov$freq*14/(365.25/12)
                spec_cov7=x_cov$spec
                
#Finally plotting the periodogram
		yli2=max(c(spec_cov,spec_cov7))
                plot(freq_cov,spec_cov,"l",xaxt="n",xlim=c(xli1,xli2),xlab="Period (month)",main=cov3[c],lty=1,lwd=alwd,ylim=c(0,yli2),cex.axis=fac_axis,cex.lab=fac_lab,cex.main=fac_main,ylab="")
                lines(freq_cov7,spec_cov7,lty=1,lwd=alwd,col="blue")
                axis(1,at=ate,labels=labe,cex.lab=fac_lab)

#Computing the acf
                for(cbis in 1:length(cov3)){
                        if (cbis!=c){
                               b=ccf(tab_cov[,cov3[c]],tab_cov[,cov3[cbis]],na.action=na.pass,plot=FALSE)
                               a=ccf(tab_cov7[,cov3[c]],tab_cov7[,cov3[cbis]],plot=FALSE,na.action=na.pass)
                               yli1=min(c(b$acf,a$acf),na.rm=TRUE)
                               yli2=max(c(b$acf,a$acf),na.rm=TRUE)
                               plot(b$lag,b$acf,t="h",col="black",main=paste(cov3[c],"vs",cov3[cbis]),ylim=c(yli1,yli2),ylab="ACF",xlab="Lag",cex.axis=fac_axis,cex.lab=fac_lab,cex.main=fac_main,lwd=alwd)
                               abline(h=0,lwd=alwd)
                               if (cov3[cbis] %in% cov_comparable){
                                       lines(a$lag+0.3,a$acf,t="h",col="blue",lwd=alwd)
                               }
                        }
               }

               plot(0,0,t="n",xaxt="n",yaxt="n",ylab="",xlab='',bty="n")
               legend("topleft",c(paste(info,"missing values"),"Teychan","Buoy 7"),bty="n",col=c(NA,"black","blue"),lty=c(NA,1,1),seg.len=0.5,cex=fac_lg)
}
dev.off()


###########Figure S1.8
fac_main=3.0
fac_axis=2.0
fac_lg=1.5
fac_lab=2.0
alwd=2.5
awidth=20

#Data
filename="../Teychan_base.csv"
tabbis=read.csv(filename,na.strings="NA",header=TRUE,sep=";",dec=".")
dates=as.Date(tabbis$Date)
filename7="../B7_base.csv"
tabbis7=read.csv(filename7,na.strings="NA",header=TRUE,sep=";",dec=".")

accum_vent=3 #We take into account wind energy 3 days before the date
tab_cov=new_covar(tabbis,dates,accum_vent,start="t")

dates_bis=seq(dates[1],dates[length(dates)],14)

t=1:length(dates_bis)
omega=2*pi/(365.25/14)
sini=sin(omega*t)
cosi=cos(omega*t)

var2=scale(na.approx(tab_cov[,"CumRg"],maxgap=1,x=dates,xout=dates_bis,na.rm=FALSE))
indi2=which(!is.na(var2))
plou2=lm(var2~sini+cosi)

var3=scale(na.approx(tab_cov[,"TEMP"],maxgap=1,x=dates,xout=dates_bis,na.rm=FALSE))
indi3=which(!is.na(var3))
plou3=lm(var3~sini+cosi)

var_exp=scale(na.approx(tab_cov[,"Ntot"],maxgap=1,x=dates,xout=dates_bis,na.rm=FALSE))
indi_exp=which(!is.na(var_exp))
pifbis=lm(var_exp~sini+cosi)

pdf("saisonnalite_comp_phase.pdf",width=awidth)
par(mfrow=c(1,1))
li1=min(min(plou2$fitted.values,na.rm=TRUE),min(plou3$fitted.values,na.rm=TRUE),min(pifbis$fitted.values,na.rm=TRUE))
li2=max(max(plou2$fitted.values,na.rm=TRUE),max(plou3$fitted.values,na.rm=TRUE),max(pifbis$fitted.values,na.rm=TRUE))
plot(dates_bis[indi2],plou2$fitted.values,col="black",t="l",xlim=c(dates_bis[1]-900,dates_bis[length(dates_bis)]),ylim=c(li1,li2),xlab="",ylab="Seasonal component",cex.main=fac_main,cex.axis=fac_axis,cex.lab=fac_lab,lwd=alwd)
lines(dates_bis[indi3],plou3$fitted.values,col="blue",t="l",lwd=alwd)
lines(dates_bis[indi2][indi_exp],pifbis$fitted.values,col="red",t="l",lwd=alwd)
legend("topleft",c("CumRg","Temp","Ntot"),col=c("black","blue","red"),lty=1,cex=fac_lg)
dev.off()

###########Figure S1.9
threshold=1.25*2
var_exp_stress=stress_function(na.approx(tab_cov[,"Ntot"],maxgap=1,x=dates,xout=dates_bis,na.rm=FALSE),threshold)
pif_stress=lm(var_exp_stress[indi3]~plou3$fitted.values)

pdf("saisonnalite_nstress.pdf",width=awidth)
li1=min(pif_stress$residuals,na.rm=TRUE)
li2=1.0
plot(dates_bis[indi3][indi_exp],pif_stress$fitted.values,ylim=c(li1,li2),t="l",xlab="",ylab="",main="N with saturation",cex.main=fac_main,cex.axis=fac_axis,cex.lab=fac_lab,lwd=alwd)
lines(dates_bis,var_exp_stress,col="red",t="p",pch=16)
points(dates_bis[indi3][indi_exp],pif_stress$residuals,col="blue",t="p",pch=16)
legend('bottomright',c("Seasonal component","Real values","Residuals"),col=c("black","red","blue"),pch=c(NA,16,16),lty=c(1,NA,NA),cex=fac_lg)
dev.off()

###########Figure S1.10--S1.12
var=c("SAL","CumRg","MeanVent")
tab_saison=matrix(NA,nrow=length(plou3$fitted.values),ncol=4)
tab_saison[,1]=plou3$fitted.values
colnames(tab_saison)=c("Season",var)
tab_sans_saison=matrix(NA,nrow=length(plou3$fitted.values),ncol=4)
tab_sans_saison[,1]=plou3$fitted.values
colnames(tab_sans_saison)=c("Season",var)
for (v in var){
        plou=scale(na.approx(tab_cov[,v],maxgap=2,x=dates,xout=dates_bis,na.rm=FALSE))
        indi_plou=which(!is.na(plou))
        plou_seas=lm(plou[indi3]~plou3$fitted.values)
        li1=min(plou_seas$residuals,na.rm=TRUE)
        li2=max(plou_seas$residuals,na.rm=TRUE)
        pdf(paste(v,"_saisonnalite.pdf",sep=""),width=awidth)
        plot(dates_bis[indi3][indi_plou],plou_seas$fitted.values,ylim=c(li1,li2),t="l",xlab="",ylab="",main=v,cex.main=fac_main,cex.axis=fac_axis,cex.lab=fac_lab,lwd=alwd)
        lines(dates_bis,plou,col="red",lty=2,pch=16,lwd=alwd)
        points(dates_bis[indi3][indi_plou],plou_seas$residuals,col="blue",t="p",pch=16)
        legend('bottomright',c("Seasonal component","Real values","Residuals"),col=c("black","red","blue"),pch=c(NA,NA,16),lty=c(1,2,NA),cex=fac_lg)
        l=0
        for (i in 1:length(plou)[1]){
                tab_sans_saison[i,v]=plou[i]
                if(is.na(plou[i])){
                        tab_saison[i,v]=NA
                }else{
                        l=l+1
                        tab_saison[i,v]=plou_seas$residuals[l]
                }
        }
        dev.off()
}

###########Figure S1.13
pdf("saisonnalite_comp_temp.pdf",width=awidth)
temp=ts(var3,start=c(year(dates_bis[1]),week(dates_bis[1])),end=c(year(tail(dates_bis,1)),week(tail(dates_bis,1))),freq=52/2) 
dd2=decompose(temp)
ddbis=stl(temp,s.window="periodic")
li1=min(min(ddbis$time.series[,"seasonal"],na.rm=TRUE),min(plou3$fitted.values,na.rm=TRUE),min(dd2$seasonal,na.rm=TRUE))
li2=max(max(ddbis$time.series[,"seasonal"],na.rm=TRUE),max(plou3$fitted.values,na.rm=TRUE),max(dd2$seasonal,na.rm=TRUE))
plot(dates_bis[indi3],plou3$fitted.values,col="black",t="l",xlim=c(dates_bis[1]-900,dates_bis[length(dates_bis)]),ylim=c(li1,li2),xlab="",ylab="Seasonal component",cex.main=fac_main,cex.axis=fac_axis,cex.lab=fac_lab,lwd=alwd)
lines(as.Date(date_decimal(c(time(dd2$seasonal)))),dd2$seasonal,col="blue",t="l",lwd=alwd)
lines(as.Date(date_decimal(c(time(ddbis$time.series)))),ddbis$time.series[,"seasonal"],col="red",t="l",lwd=alwd,lty=2)
legend("topleft",c("LM","Decompose","STL"),col=c("black","blue","red"),lty=c(1,1,2),cex=fac_lg)
dev.off()

###########Figure S1.14
cov3_phy=c("SAL","CumRg","MeanVent")
cov3_nut=c('Ntot','SI','PHOS')
cov3_tot=c(cov3_phy,"TEMP") #covariates we are going to use, including temperature

###Loading data
filename="../Teychan_base.csv"
tabbis=read.csv(filename,na.strings="NA",header=TRUE,sep=";",dec=".")
dates=as.Date(tabbis$Date)

###Treating covariates to have cumulated radiation and wind energy
accum_vent=3 #We take into account wind energy 3 days before the date
tab_cov=new_covar(tabbis,dates,accum_vent)

###Treating missing values for both species and covariates###
consecutif=2 #Number of missing values above which we keep the NA
timestep=14 #Regular time lapses between two observations
dates_bis=seq(dates[1],dates[length(dates)],timestep) #Regular time grid
tab_cov_bis=matrix(NA,length(dates_bis),length(cov3_tot))
colnames(tab_cov_bis)=cov3_tot
for (c in cov3_tot){
        tab_cov_bis[,c]=approx(tab_cov[,c],x=dates,xout=dates_bis)$y
}

fac_main=3.0
fac_axis=2.0
fac_lg=1.5
fac_lab=2.0
alwd=2.5
awidth=20
cex_pch=3

###Cosine method
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

#Average year
t=1:365
yy=unique(year(dates))
temps_journalier=c(parse_date_time("1/1987","%j/%Y"))
for (y in yy){
        for (tt in t){
                temps_journalier=c(temps_journalier,parse_date_time(paste(tt,'/',y,sep=""),"%j/%Y"))
        }
}
temps_journalier=temps_journalier[-1]

real_time=parse_date_time(dates,"%Y-%m-%d")
real_time_bis=parse_date_time(dates_bis,"%Y-%m-%d")
data=approx(tab_cov[,'TEMP'],x=real_time,xout=temps_journalier)$y

model_data=matrix(NA,nrow=365,ncol=length(yy))

for (y in 1:length(yy)){
        model_data[,y]=data[year(temps_journalier)==yy[y]]
}

mean_year=apply(model_data,1,mean,na.rm=TRUE)
plou=rep(mean_year,length(yy))

#Plot
pdf("comparaison_t.pdf",width=10,height=7)
par(mar=c(3,6,4,0.5))
plot(real_time_bis,season$fitted.values,col="red",pch=16,xlim=c(real_time_bis[1],real_time_bis[25]),ylim=c(7,23),ylab="Temperature",xlab="",main="Method comparison",cex.main=fac_main,cex.axis=fac_axis,cex.lab=fac_lab,lwd=alwd,cex=cex_pch,xaxt="n")
points(temps_journalier,plou,cex=cex_pch*0.75)
axis.POSIXct(1,at=seq(as.POSIXct("1987-01-01",tz="UTC"),as.POSIXct("1988-01-15",tz="UTC"),"1 month"),format="%d/%m",cex.axis=fac_axis*1)
legend("topleft",c("Cosine method","Average year"),col=c("red","black"),pch=c(16,1),cex=fac_lg)
dev.off()
          
###########Figure S1.15
MS=cor(tab_sans_saison,method="spearman",use="pairwise.complete.obs")
M=cor(tab_saison,method="spearman",use="pairwise.complete.obs")
pdf("ARLM_corrplot_after_season.pdf",width=10,height=10)
corrplot(MS,type="upper",method="number",col="black",cl.pos="n",tl.pos="lt")
corrplot(M,type="lower",method="number",col="black",add=TRUE,cl.pos="n",tl.pos="n")
dev.off()


###########Figure S1.16
#Phytoplankton
sp=c("GUI", "LEP", "NIT", "PSE", "RHI", "SKE", "CHA", "AST", "GYM",  "PRP", "EUG", "CRY")
sp=sort(sp)
#Data
filename="../Teychan_base.csv"
tab=read.csv(filename,na.strings="NA",header=TRUE,sep=";",dec=".")
dates=as.Date(tab$Date)
dates_correct=dates[year(dates)>1996]
#New
filename7="../B7_base.csv"
tab7=read.csv(filename7,na.strings="NA",header=TRUE,sep=";",dec=".")
dates7=as.Date(tab7$Date)
dates_correct7=dates7[year(dates7)>1996]

dates_bis=seq(dates_correct[1],dates_correct[length(dates_correct)],14)
dates_bis7=seq(dates_correct7[1],dates_correct7[length(dates_correct7)],14)

dates_mois=seq(as.Date("01/01/1997","%d/%m/%Y"),dates_correct[length(dates_correct)], by = "month")
dates_mois7=seq(as.Date("01/02/2003","%d/%m/%Y"),dates_correct7[length(dates_correct7)], by = "month")

dates_saison=seq(as.Date("21/12/1996","%d/%m/%Y"),dates_correct[length(dates_correct)], by = "quarter")
dates_saison7=seq(as.Date("21/12/2002","%d/%m/%Y"),dates_correct7[length(dates_correct7)], by = "quarter")

for(s in sp){
        #na.approx
	s1=na.approx(tab[,s],maxgap=2,x=dates,xout=dates_bis,na.rm=FALSE)
        s17=na.approx(tab7[,s],maxgap=2,x=dates7,xout=dates_bis7,na.rm=FALSE)

        #Monthly
        s2=rep(NA,length(dates_mois))
        s27=rep(NA,length(dates_mois7))
        for (d in 1:(length(dates_mois)-1)){
                d1=dates_mois[d]
                d2=dates_mois[d+1]
                s2[d]=mean(tab[(dates>=d1)&(dates<d2),s],na.rm=TRUE)
        }
        s2[length(s2)]=mean(tab[dates>=d2,s],na.rm=TRUE)
        for (d in 1:(length(dates_mois7)-1)){
                d1=dates_mois7[d]
                d2=dates_mois7[d+1]
                s27[d]=mean(tab7[(dates7>=d1)&(dates7<d2),s],na.rm=TRUE)
        }
        s27[length(s27)]=mean(tab7[dates7>=d2,s],na.rm=TRUE)

        #Quarterly
        s3=rep(NA,length(dates_saison))
        s37=rep(NA,length(dates_saison7))
        for (d in 1:(length(dates_saison)-1)){
                d1=dates_saison[d]
                d2=dates_saison[d+1]
                s3[d]=mean(tab[(dates>=d1)&(dates<d2),s],na.rm=TRUE)
        }
        s3[length(s3)]=mean(tab[dates>=d2,s],na.rm=TRUE)
        for (d in 1:(length(dates_saison7)-1)){
                d1=dates_saison7[d]
                d2=dates_saison7[d+1]
                s37[d]=mean(tab7[(dates7>=d1)&(dates7<d2),s],na.rm=TRUE)
        }
        s37[length(s37)]=mean(tab7[dates7>=d2,s],na.rm=TRUE)

        #SSA-based
        s_bis=log(na.approx(tab[dates%in%dates_correct,s],maxgap=0,x=dates_correct,xout=dates_bis,na.rm=FALSE))
        s_bis7=log(na.approx(tab7[dates7%in%dates_correct7,s],maxgap=0,x=dates_correct7,xout=dates_bis7,na.rm=FALSE))
        s4=gapfillSSA(series=s_bis,amnt.iters=c(20,20))$filled.series
        s47=gapfillSSA(series=s_bis7,amnt.iters=c(20,20))$filled.series

	#Plot
        pdf(paste("differents_groupements_temporels_",s,".pdf",sep=""),width=20)
        par(mfrow=c(2,1),mar=c(4.1, 5.1, 4.1, 1.1))
	#Teychan
        xli1=as.Date("01/11/1995","%d/%m/%Y")-700
        xli2=dates[length(dates)]
        yli1=log(min(tab[,s],na.rm=TRUE))
        yli2=log(max(tab[,s],na.rm=TRUE))
        plot(dates[dates>dates_saison[1]],log(tab[dates>dates_saison[1],s]),pch=16,xlim=c(xli1,xli2),ylim=c(yli1,yli2),main=paste("Teychan, correlation interpo/spectral :",format(cor(s1,s4,method="spearman",use="pairwise.complete.obs"),digits=3)),xlab="",ylab="Log(Abundance)",cex.main=fac_main,cex.axis=fac_axis,cex.lab=fac_lab)
        lines(dates_bis,log(s1),pch=16,lty=1,t="o",col="red")
        points(dates_mois,log(s2),pch=16,col="blue")
        points(dates_saison,log(s3),pch=16,col="green")
        points(dates_bis,s4,pch=16,col="grey")
        legend("topleft",c(paste("Real data",format(sum(is.na(tab[dates %in% dates_correct,s]))/length(tab[dates %in% dates_correct,s]),digits=3),"NA"),paste("na.approx",format(sum(is.na(s1))/length(s1),digits=3),"NA"),paste("Monthly",format(sum(is.na(s2))/length(s2),digits=3),"NA"),paste("Quarterly",format(sum(is.na(s3))/length(s3),digits=3),"NA"),paste("SSA",format(sum(is.na(s4))/length(s4),digits=3),"NA")),col=c("black","red","blue","green","grey"),pch=16,lty=1,cex=fac_lg,bty="n")
	#B7
        xli1=as.Date("01/11/2001","%d/%m/%Y")-700
        xli2=dates7[length(dates7)]
        yli1=log(min(tab7[,s],na.rm=TRUE))
        yli2=log(max(tab7[,s],na.rm=TRUE))
        plot(dates7[dates7>dates_saison7[1]],log(tab7[dates7>dates_saison7[1],s]),pch=16,cex=1.5,xlim=c(xli1,xli2),ylim=c(yli1,yli2),main=paste("Buoy7, correlation interpo/spectral :",format(cor(s17,s47,method="spearman",use="pairwise.complete.obs"),digits=3)),xlab="",ylab="Log(Abundance)",cex.main=fac_main,cex.axis=fac_axis,cex.lab=fac_lab)
        lines(dates_bis7,log(s17),pch=16,lty=1,t="o",col="red")
        points(dates_mois7,log(s27),pch=16,col="blue")
        points(dates_saison7,log(s37),pch=16,col="green")
        points(dates_bis7,s47,pch=16,col="grey")
        legend("topleft",c(paste("Real data",format(sum(is.na(tab7[dates7 %in% dates_correct7,s]))/length(tab7[dates7 %in% dates_correct7,s]),digits=3),"NA"),paste("na.approx",format(sum(is.na(s17))/length(s17),digits=3),"NA"),paste("Monthly",format(sum(is.na(s27))/length(s27),digits=3),"NA"),paste("Quarterly",format(sum(is.na(s37))/length(s37),digits=3),"NA"),paste("SSA",format(sum(is.na(s47))/length(s47),digits=3),"NA")),col=c("black","red","blue","green","grey"),pch=16,lty=1,bty="n",cex=fac_lg)

        dev.off()
}

###########Figure S1.17
library("lubridate")
pdf("spectrum_flore_bis_for_each_species.pdf")
par(mfrow=c(3,2))
xli1=8
xli2=72
maxi_test=100
for (s in sp){
	if(s=="CRY"|s=="EUG"|s=="GYM"){
		dates_correct=dates[year(dates)>1996]
		var=tabbis[year(dates)>1996,s]
	}else{
		dates_correct=dates
		var=tabbis[,s]
	}
	dates_bis=seq(dates_correct[1],dates_correct[length(dates_correct)],14)
	var1=na.approx(var,maxgap=2,x=dates_correct,xout=dates_bis,na.rm=FALSE)
	var1_bis=matrix(NA,nrow=maxi_test,ncol=length(var1))
	var1_bis[1,]=var1
	var1_bis[1,is.na(var1)]=10 #Minimum detectable abundance

	var2=na.approx(tabbis7[,s],maxgap=2,x=dates7,xout=dates_bis7,na.rm=FALSE)
	var2_bis=matrix(NA,nrow=maxi_test,ncol=length(var2))
	var2_bis[1,]=var2
	var2_bis[1,is.na(var2)]=10

	xlog=spectrum(log10(var1_bis[1,]),plot=FALSE,method="ar")
	x2log=spectrum(log10(var2_bis[1,]),plot=FALSE,method="ar")
	freq=1/xlog$freq*14/(365.25/12)
	freq2=1/x2log$freq*14/(365.25/12)
	
	xlog_bis=matrix(NA,nrow=maxi_test,ncol=length(xlog$spec))
	x2log_bis=matrix(NA,nrow=maxi_test,ncol=length(x2log$spec))
	xlog_bis[1,]=xlog$spec
	x2log_bis[1,]=x2log$spec

	for (i in 2:maxi_test){
		var1_bis[i,]=var1
		var1_bis[i,is.na(var1)]=runif(sum(is.na(var1)),0,min(var1,na.rm=TRUE))
		var2_bis[i,]=var2
		var2_bis[i,is.na(var2)]=runif(sum(is.na(var2)),0,min(var2,na.rm=TRUE))

		xlog_bis[i,]=spectrum(log10(var1_bis[i,]),plot="FALSE",method="ar")$spec
		x2log_bis[i,]=spectrum(log10(var2_bis[i,]),plot="FALSE",method="ar")$spec
	}	

	yli1=min(xlog_bis)
	yli2=max(xlog_bis)
	plot(freq,xlog_bis[1,],"l",xlim=c(xli1,xli2),xaxt="n",xlab="Periode (months)",main=paste("log10",s),lty=1,ylim=c(yli1,yli2),cex=1.3,ylab="Spectrum",lwd=3)
	for(i in 2:maxi_test){
		lines(freq,xlog_bis[i,],lty=3,col="grey",cex=0.1)
	}
	par(new=TRUE)
	yli1=min(x2log_bis)
	yli2=max(x2log_bis)
	plot(freq2,x2log_bis[1,],xlim=c(xli1,xli2),"l",col="blue",xaxt="n",yaxt="n",ylab="",xlab="",ylim=c(yli1,yli2),cex=1.3,lwd=3)
	for(i in 2:maxi_test){
		lines(freq2,x2log_bis[i,],lty=3,col="cyan",cex=0.1)
	}
	legend("topright",c("Teychan, fixed","Teychan, random","B7, fixed","B7, random"),col=c("black","grey","blue","cyan"),lty=c(1,2),bty="n")
	labe=seq(6,12*xli2,6)
	axis(1,at=labe,labels=labe)
	axis(4,col="blue")
}
dev.off()

###########Figure S1.18
sp=c("GUI", "LEP", "NIT", "PSE", "RHI", "SKE", "CHA", "AST", "GYM",  "PRP", "EUG", "CRY")
sp=sort(sp)
cov3=c("CumDebit","MeanNAO","Ntot","CumPrec","CumRg","SAL","TEMP","MeanVent","PHOS","SI")
cov_comparable=c("Ntot","SAL","TEMP","PHOS","SI")
cov3=sort(cov3)

fac_main=3.0
fac_axis=2.0
fac_lg=2.0
fac_lab=2.0
alwd=2.5

maxi_test=10
accum_vent=3
ate=seq(3,12*12,3)
labe=ate
labe[seq(1,length(ate),2)]=NA
xli1=8
xli2=72
m=3

filename=paste("../Teychan_base.csv",sep="")
tab=read.csv(filename,na.strings="NA",header=TRUE,sep=";",dec=".")
dates=as.Date(tab$Date)
dates_bis=seq(dates[dates>as.Date("01/01/1997","%d/%m/%Y")][1],dates[length(dates)],14)
tab_cov=new_covar(tab,dates,accum_vent)

filename2=paste("../B7_base.csv",sep="")
tab7=read.csv(filename2,na.strings="NA",header=TRUE,sep=";",dec=".")
dates7=as.Date(tab7$Date)
dates_bis7=seq(dates7[1],dates7[length(dates7)],14)
tab_cov7=new_covar(tab7,dates7,accum_vent)

pdf("spectrum_essai1.pdf",width=20)
par(mfrow=c(3,4),mar=c(4.1, 5.1, 4.1, 1.1))
for (s in sp){
        var1=log(na.approx(tab[,s],maxgap=2,x=dates,xout=dates_bis,na.rm=FALSE))
        var1_bis=matrix(NA,nrow=maxi_test,ncol=length(var1))
        var2=log(na.approx(tab7[,s],maxgap=2,x=dates7,xout=dates_bis7,na.rm=FALSE))
        var2_bis=matrix(NA,nrow=maxi_test,ncol=length(var2))

        s_bis=log(na.approx(tab[,s],maxgap=0,x=dates,xout=dates_bis,na.rm=FALSE))
        s_bis7=log(na.approx(tab7[,s],maxgap=0,x=dates7,xout=dates_bis7,na.rm=FALSE))
        s2=gapfillSSA(series=s_bis,amnt.iters=c(5,5),fill.margins=TRUE)$filled.series
        s27=gapfillSSA(series=s_bis7,amnt.iters=c(5,5),fill.margins=TRUE)$filled.series #Warning : normally, it should only affect a few points at the end of the series, but you never know

        x=spec.pgram(s2,kernel("modified.daniell",m),taper=0,plot=FALSE)
        freq=1/x$freq*14/(365.25/12)

        spec=x$spec
        x7=spec.pgram(s27,kernel("modified.daniell",m),taper=0,plot=FALSE)
        freq7=1/x7$freq*14/(365.25/12)
	spec7=x7$spec

        xbis=array(NA,dim=c(maxi_test,length(x$freq),2))
        xbis7=array(NA,dim=c(maxi_test,length(x7$freq),2))
        for(i in 1:maxi_test){
                var1_bis[i,]=var1
                var1_bis[i,is.na(var1)]=runif(sum(is.na(var1)),0,min(var1,na.rm=TRUE)/2)
                x=spec.pgram(var1_bis[i,],kernel("modified.daniell",m),taper=0,plot=FALSE)
                xbis[i,,1]=1/x$freq*14/(365.25/12)
                xbis[i,,2]=x$spec

                var2_bis[i,]=var2
                var2_bis[i,is.na(var2)]=runif(sum(is.na(var2)),0,min(var2,na.rm=TRUE)/2)
                x=spec.pgram(var2_bis[i,],kernel("modified.daniell",m),taper=0,plot=FALSE)
                xbis7[i,,1]=1/x$freq*14/(365.25/12)
                xbis7[i,,2]=x$spec
        }
        yli1=min(c(xbis[,,2],xbis7[,,2]))
        yli2=max(c(xbis[,,2],xbis7[,,2]))

        plot(freq,spec,xlim=c(xli1,xli2),xaxt="n",xlab="Period (month)",ylab="",main=s,lty=1,lwd=3,ylim=c(yli1,yli2),cex.main=fac_main,t="n",cex.lab=fac_lab,cex.axis=fac_axis)
	axis(1,at=ate,labels=labe,cex.axis=fac_axis,cex.lab=fac_lab)
        for(i in 1:maxi_test){
                lines(xbis[i,,1],xbis[i,,2],lty=1,lwd=alwd,col="black")
                lines(xbis7[i,,1],xbis7[i,,2],lty=1,lwd=alwd,col="blue")
        }

        for(c in 1:length(cov3)){
                covar1=approx(tab_cov[,cov3[c]],x=dates,xout=dates_bis)$y
                if(sum(is.na(covar1))>0){ 
                        id=!is.na(covar1)
                        info=sum(is.na(covar1))
                        covar1=covar1[id]
                        vv=var1_bis[1,id]
                        dd=dates_bis[id]
                }
                else{
                        info=0
                        vv=var1_bis[1,]
                        dd=dates_bis
                }
                x_coherency=spec.pgram(cbind(vv,covar1),kernel("modified.daniell",c(m,m)),taper=0,plot=FALSE)
                L=2*m+1 #L=2m+1
                f = qf(.995, 2, x_coherency$df-2) #alpha=5%, corrected Bonferroni because 10 variables
                C = f/(L+f)
                freq=1/x_coherency$freq*14/(365.25/12)
                spec=x_coherency$coh

                covar17=approx(tab_cov7[,cov3[c]],x=dates7,xout=dates_bis7)$y
                if(sum(is.na(covar17))>0){
                        id=!is.na(covar17)
                        info=sum(is.na(covar17))
                        covar17=covar17[id]
                        vv7=var2_bis[1,id]
                        dd7=dates_bis7[id]
                }
                else{
                        info=0
                        vv7=var2_bis[1,]
                        dd7=dates_bis7
                }
                x_coherency=spec.pgram(cbind(vv7,covar17),kernel("modified.daniell",c(m,m)),taper=0,plot=FALSE)
                L=2*m+1 #L=2m+1
                f = qf(.995, 2, x_coherency$df-2) #alpha=5%, corrected Bonferroni because 10 variables
                C = f/(L+f)
                freq7=1/x_coherency$freq*14/(365.25/12)
                spec7=x_coherency$coh

                yli2=max(C,max(spec),max(spec7))
                plot(freq,spec,"l",xlim=c(xli1,xli2),xaxt="n",xlab="Period (month)",main=paste("Coherence",s,cov3[c]),lty=1,ylim=c(0,max(C,max(spec),max(spec7))),cex.main=fac_main,lwd=alwd,cex.lab=fac_lab,cex.axis=fac_axis,ylab="")
                abline(h = C,col="black",lty=2,lwd=alwd)
                axis(1,at=ate,labels=labe,cex.axis=fac_axis,cex.lab=fac_lab)
                lines(freq7,spec7,"l",lty=1,col="blue",lwd=alwd)
                abline(h = C,col="blue",lty=2,lwd=alwd)
        }
                plot(0,0,t="n",xaxt="n",yaxt="n",ylab="",xlab='',bty="n")
                legend('topleft',c('Teychan','Buoy7','Spectrum','0.5% Significance'),lty=c(1,1,1,2),col=c('black',"blue","black","black"),bty="n",cex=fac_lg)
}
dev.off()


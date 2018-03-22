#This script compares the confidence intervals obtained from the Hessian computation and the bootstrap for one scenario of Simulated Data, taking into account the real value
rm(list=ls())
graphics.off()
dirin="./"

iter=1:10
alpha=0.05
ll=0.05
repet=1:10

amain=3.5
aaxis=2.5
alab=2.5
alwd=4
awidth=40
aheight=27
aleg=3.0
apch=3.0
atitle=6


id_diag=seq(1,100,11)

biais_array_global=array(NA,dim=c(3,length(repet),3,2,2)) #3 scenarios, 10 repetions, 3 coefficient types - B diag, B out of diag, C whole ; 2 versions (lagged covariate, hereafter called v2, or not); 2 models
sign_diffpm_array_signif=array(NA,dim=c(3,length(repet),3,2,2)) #significant values change from positive in real data to negative in false data ; 2 versions (lagged covariate or not) ; 2 models
sign_diffmp_array_signif=array(NA,dim=c(3,length(repet),3,2,2)) #significant values change from negative to positive ; 2 versions (lagged covariate or not) ; 2 models
proportion_val_ok=array(NA,dim=c(3,length(repet),3,2,2)) #percentage of value inside the 95% confidence intervals

for (scenario in 1:3){
dirin_estimation_v2=paste(dirin,"analyse/v2/scenario",scenario,"/",sep="")
dirin_estimation_v1=paste(dirin,"analyse/scenario",scenario,"/",sep="")

for (r in iter){
	load(paste(dirin,"data/SimulatedParameters/SimulatedParameters_repeat",r,".RData",sep=""))
	if(scenario==1){
        B_used=B1
	model=c("null_noseason","null_season")
	}else{
	B_used=B2
	model=c("noseason","season")
	}
		if(scenario==2){
			C=rep(0,length(C))
		}else{
			C=c(c1,c2)
		}

	
	fp_s2_cov_v1=0
	fp_s2_cov_v2=0
	fn_s2_inter_v1=0
	fn_s2_inter_v2=0
	fn_s3_inter_v1=0
	fn_s3_inter_v2=0


	for (m in 1:2){
	
	#First estimates from v2
		if(scenario==1){
		load(paste(dirin_estimation_v2,"MockSeries_scenario",scenario,"_rep",r,"_kem_",model[m],"_withAIC.RData",sep=""))
		}else{
		load(paste(dirin_estimation_v2,"MockSeries_scenario",scenario,"_rep",r,"_kem_",model[m],"_v2.RData",sep=""))
		}
		fit_log_v2=fit_log
		aic_v2=aic
		maxB=length(fit_log_v2$par$B)
		if(maxB<100){
			B_est_v2=diag(c(fit_log_v2$par$B))
		}else{
			B_est_v2=matrix(c(fit_log_v2$par$B),nrow=10,ncol=10)
		}
		C_est_v2=matrix(c(fit_log_v2$par$U),nrow=10)
		
        #First estimates from v1
                load(paste(dirin_estimation_v1,"MockSeries_scenario",scenario,"_rep",r,"_kem_",model[m],".RData",sep=""))
		fit_log_v1=fit_log
		cis_v1=cis
                if(maxB<100){
                        B_est_v1=diag(c(fit_log_v1$par$B))
                }else{
                        B_est_v1=matrix(c(fit_log_v1$par$B),nrow=10,ncol=10)
                }
                C_est_v1=matrix(c(fit_log_v1$par$U),nrow=10)

	if(scenario==1){	
	#Confidence intervals from boostrap for v2
		B_sup_vec=apply(aic_v2$boot.params[1:maxB,], 1, quantile, probs = 1-alpha/2)
		B_inf_vec=apply(aic_v2$boot.params[1:maxB,], 1, quantile, probs = alpha/2)
		if(maxB<100){
			B_sup_b_v2=diag(B_sup_vec)
			B_inf_b_v2=diag(B_inf_vec)
		}else{
			B_sup_b_v2=matrix(B_sup_vec,nrow=10,ncol=10)
			B_inf_b_v2=matrix(B_inf_vec,nrow=10,ncol=10)
		}
		C_sup_b_v2=matrix(apply(aic_v2$boot.params[(maxB+1):(maxB+length(fit_log_v2$par$U)),], 1, quantile, probs = 1-alpha/2),nrow=10)
		C_inf_b_v2=matrix(apply(aic_v2$boot.params[(maxB+1):(maxB+length(fit_log_v2$par$U)),], 1, quantile, probs = alpha/2),nrow=10)
	}else{
		B_sup_b_v2=matrix(NA,nrow=10,ncol=10)
		B_inf_b_v2=matrix(NA,nrow=10,ncol=10)
		C_sup_b_v2=matrix(NA,nrow=10,ncol=length(fit_log_v2$par$U)/10)
		C_inf_b_v2=matrix(NA,nrow=10,ncol=length(fit_log_v2$par$U)/10)
	}

	#Confidence intervals from hessian for v2
		if(scenario==1){
		load(paste(dirin_estimation_v2,"hessian_MockSeries_scenario",scenario,"_rep",r,"_kem_",model[m],"_withAIC.RData",sep=""))
		}else{
		load(paste(dirin_estimation_v2,"hessian_MockSeries_scenario",scenario,"_rep",r,"_kem_",model[m],"_v2.RData",sep=""))
		}
		if(maxB<100){
			B_sup_h_v2=diag(c(hi$par.upCI$B))
			B_inf_h_v2=diag(c(hi$par.lowCI$B))
		}else{
			B_sup_h_v2=matrix(c(hi$par.upCI$B),nrow=10,ncol=10)
			B_inf_h_v2=matrix(c(hi$par.lowCI$B),nrow=10,ncol=10)
		}
		C_sup_h_v2=matrix(c(hi$par.upCI$U),nrow=10)
		C_inf_h_v2=matrix(c(hi$par.lowCI$U),nrow=10)

        #Confidence intervals from bootstrap for v1
                load(paste(dirin_estimation_v1,"MockSeries_scenario",scenario,"_rep",r,"_kem_",model[m],".RData",sep=""))
                if(maxB<100){
                        B_sup_b_v1=diag(c(cis_v1$par.upCI$B))
                        B_inf_b_v1=diag(c(cis_v1$par.lowCI$B))
                }else{
                        B_sup_b_v1=matrix(c(cis_v1$par.upCI$B),nrow=10,ncol=10)
                        B_inf_b_v1=matrix(c(cis_v1$par.lowCI$B),nrow=10,ncol=10)
                }
                C_sup_b_v1=matrix(c(cis_v1$par.upCI$U),nrow=10)
                C_inf_b_v1=matrix(c(cis_v1$par.lowCI$U),nrow=10)
		
		miniB=min(c(min(c(B_inf_b_v2),na.rm=TRUE),min(c(B_inf_h_v2),na.rm=TRUE),min(B_used)))
		maxiB=max(c(max(c(B_sup_b_v2),na.rm=TRUE),max(c(B_sup_h_v2),na.rm=TRUE),max(B_used)))
		miniC=min(c(min(c(C_inf_b_v2),na.rm=TRUE),min(c(C_inf_h_v2),na.rm=TRUE),min(C)))
		maxiC=max(c(max(c(C_sup_b_v2),na.rm=TRUE),max(c(C_sup_h_v2),na.rm=TRUE),max(C)))
	
		biais_array_global[scenario,r,1,1,m]=mean(diag(B_used)-diag(B_est_v1),na.rm=TRUE)
		biais_array_global[scenario,r,1,2,m]=mean(diag(B_used)-diag(B_est_v2),na.rm=TRUE)
		biais_array_global[scenario,r,2,1,m]=mean(c(B_used)[-id_diag]-c(B_est_v1)[-id_diag],na.rm=TRUE)
		biais_array_global[scenario,r,2,2,m]=mean(c(B_used)[-id_diag]-c(B_est_v2)[-id_diag],na.rm=TRUE)
		proportion_val_ok[scenario,r,1,1,m]=(sum(diag(B_used)>=diag(B_inf_b_v1)&diag(B_used)<=diag(B_sup_b_v1),na.rm=TRUE))/10
		proportion_val_ok[scenario,r,1,2,m]=(sum(diag(B_used)>=diag(B_inf_h_v2)&diag(B_used)<=diag(B_sup_h_v2),na.rm=TRUE))/10
		proportion_val_ok[scenario,r,2,1,m]=(sum(c(B_used)[-id_diag]>=c(B_inf_b_v1)[-id_diag]&c(B_used)[-id_diag]<=c(B_sup_b_v1)[-id_diag],na.rm=TRUE))/90
		proportion_val_ok[scenario,r,2,2,m]=(sum(c(B_used)[-id_diag]>=c(B_inf_h_v2)[-id_diag]&c(B_used)[-id_diag]<=c(B_sup_h_v2)[-id_diag],na.rm=TRUE))/90
		if(m==1){
			biais_array_global[scenario,r,3,1,m]=mean(c(C)-c(C_est_v1),na.rm=TRUE)
			biais_array_global[scenario,r,3,2,m]=mean(c(C)-c(C_est_v2),na.rm=TRUE)
			proportion_val_ok[scenario,r,3,1,m]=(sum(C>=C_inf_b_v1&C<=C_sup_b_v1,na.rm=TRUE))/20
			proportion_val_ok[scenario,r,3,2,m]=(sum(C>=C_inf_h_v2&C<=C_sup_h_v2,na.rm=TRUE))/20
		}

	                
		signifB_v1=sign(B_inf_b_v1)*sign(B_sup_b_v1)
		signifC_v1=sign(C_inf_b_v1)*sign(C_sup_b_v1)
		signifB_v2=sign(B_inf_h_v2)*sign(B_sup_h_v2)
		signifC_v2=sign(C_inf_h_v2)*sign(C_sup_h_v2)

                siB1=signifB_v1
                siB1[signifB_v1<0]=NA
                siB2=signifB_v2
                siB2[signifB_v2<0]=NA
                siC1=signifC_v1
                siC1[signifC_v1<0]=NA
                siC2=signifC_v2
                siC2[signifC_v2<0]=NA

		B_est_v1_si=B_est_v1*siB1
		B_est_v2_si=B_est_v2*siB2
		C_est_v1_si=C_est_v1*siC1
		C_est_v2_si=C_est_v2*siC2
			
				if(m==1){
			if(scenario==2){
				fp_s2_cov_v1=fp_s2_cov_v1+sum(!is.na(C_est_v1_si))/10
				fp_s2_cov_v2=fp_s2_cov_v2+sum(!is.na(C_est_v2_si))/10
				fn_s2_inter_v1=fn_s2_inter_v1+sum(is.na(B_est_v1_si)[-id_diag])/10
				fn_s2_inter_v2=fn_s2_inter_v2+sum(is.na(B_est_v2_si)[-id_diag])/10
			}
			if(scenario==3){
				fn_s3_inter_v1=fn_s3_inter_v1+sum(is.na(B_est_v1_si)[-id_diag])/10
				fn_s3_inter_v2=fn_s3_inter_v2+sum(is.na(B_est_v2_si)[-id_diag])/10
			}
				}

		
		sign_diffpm_array_signif[scenario,r,1,1,m]=sum(sign(diag(B_used))>0&sign(diag(B_est_v1_si))<0,na.rm=TRUE)
		sign_diffpm_array_signif[scenario,r,1,2,m]=sum(sign(diag(B_used))>0&sign(diag(B_est_v2_si))<0,na.rm=TRUE)
		sign_diffpm_array_signif[scenario,r,2,1,m]=sum(sign(B_used)[-id_diag]>0&sign(B_est_v1_si)[-id_diag]<0,na.rm=TRUE)
		sign_diffpm_array_signif[scenario,r,2,2,m]=sum(sign(B_used)[-id_diag]>0&sign(B_est_v2_si)[-id_diag]<0,na.rm=TRUE)
		sign_diffpm_array_signif[scenario,r,3,1,m]=sum(sign(C)>0&sign(C_est_v1_si)<0,na.rm=TRUE)
		sign_diffpm_array_signif[scenario,r,3,2,m]=sum(sign(C)>0&sign(C_est_v2_si)<0,na.rm=TRUE)

		sign_diffmp_array_signif[scenario,r,1,1,m]=sum(sign(diag(B_used))<0&sign(diag(B_est_v1_si))>0,na.rm=TRUE)
		sign_diffmp_array_signif[scenario,r,1,2,m]=sum(sign(diag(B_used))<0&sign(diag(B_est_v2_si))>0,na.rm=TRUE)
		sign_diffmp_array_signif[scenario,r,2,1,m]=sum(sign(B_used)[-id_diag]<0&sign(B_est_v1_si)[-id_diag]>0,na.rm=TRUE)
		sign_diffmp_array_signif[scenario,r,2,2,m]=sum(sign(B_used)[-id_diag]<0&sign(B_est_v2_si)[-id_diag]>0,na.rm=TRUE)
		sign_diffmp_array_signif[scenario,r,3,1,m]=sum(sign(C)<0&sign(C_est_v1_si)>0,na.rm=TRUE)
		sign_diffmp_array_signif[scenario,r,3,2,m]=sum(sign(C)<0&sign(C_est_v2_si)>0,na.rm=TRUE)

}
}
	if(scenario==2){
	print(paste("Pour scenario",scenario,"il y a",fp_s2_cov_v1,"faux positifs dans v1 et",fp_s2_cov_v2,"dans v2"))
	print(paste("Pour scenario",scenario,"il y a",fn_s2_inter_v1,"faux negatifs dans v1 et",fn_s2_inter_v2,"dans v2"))
	}
	if(scenario==3){
	print(paste("Pour scenario",scenario,"il y a",fn_s3_inter_v1,"faux negatifs dans v1 et",fn_s3_inter_v2,"dans v2"))
	}
}

pdf("/home/cpicoche/Documents/Plankton/essai_graphe_bias_pour_article.pdf",width=awidth,height=aheight/1.5)
par(mfrow=c(1,3),mar=c(8,10,10,4))
name_scenar=c("Evt","Inter seul","Evt+Inter")
name_scenar=c("S1","S2","S3")
tipi=c("Intra","Inter","Covariates")
tipi=c("Intra-group","Inter-group","Covariates")
for (tp in 1:3){
mini=min(biais_array_global[,,tp,,],na.rm=TRUE)
maxi=max(biais_array_global[,,tp,,],na.rm=TRUE)
mini=min(biais_array_global[,,tp,2,],na.rm=TRUE)
maxi=max(biais_array_global[,,tp,2,],na.rm=TRUE)
mini=-0.015
maxi=0.4
plot(1:3,apply(biais_array_global[,,tp,2,1],1,mean,na.rm=TRUE),xlim=c(0.7,4),ylim=c(mini,maxi),xlab="",ylab="",cex.main=amain*2,cex.axis=aaxis*2,cex.lab=alab*2,pch=16,cex=apch*2,xaxt="n",col="black")
arrows(1:3,apply(biais_array_global[,,tp,2,1],1,min,na.rm=TRUE),1:3,apply(biais_array_global[,,tp,2,1],1,max,na.rm=TRUE),code=3,length=0.2,angle=90,col='black',lwd=alwd*2.5)

abline(h=0.0,lty=3,lwd=alwd/2)
axis(1,at=1:3,cex.axis=aaxis*2,labels=FALSE,cex.axis=aaxis*2)
mtext(name_scenar,side=1,line=5.5,at=1:3,cex=amain)
title(tipi[tp],line=2,cex.main=amain*2)
if(tp==1){
mtext("% Bias",side=2,line=5.5,at=0.2,cex=amain)
}
}
dev.off()

pdf("/home/cpicoche/Documents/Plankton/essai_graphe_propok_pour_article.pdf",width=awidth,height=aheight/1.5)
par(mfrow=c(1,3),mar=c(8,15,10,4))
name_scenar=c("Evt seul","Inter seul","Evt+Inter")
name_scenar=c("S1","S2","S3")
tipi=c("Intra","Inter","Covariates")
tipi=c("Intra-group","Inter-group","Covariates")
for (tp in 1:3){
mini=min(proportion_val_ok[,,tp,,],na.rm=TRUE)
maxi=max(proportion_val_ok[,,tp,,],na.rm=TRUE)
mini=min(proportion_val_ok[,,tp,2,],na.rm=TRUE)*100
maxi=max(proportion_val_ok[,,tp,2,],na.rm=TRUE)*100
mini=0
maxi=100
plot(1:3,100*apply(proportion_val_ok[,,tp,2,1],1,mean,na.rm=TRUE),xlim=c(0.7,4),ylim=c(mini,maxi),xlab="",ylab="",cex.main=amain*2,cex.axis=aaxis*2,cex.lab=alab*2,pch=16,cex=apch*2,xaxt="n",col="black")
arrows(1:3,100*apply(proportion_val_ok[,,tp,2,1],1,min,na.rm=TRUE),1:3,100*apply(proportion_val_ok[,,tp,2,1],1,max,na.rm=TRUE),code=3,length=0.2,angle=90,col='black',lwd=alwd*2.5)


axis(1,at=1:3,cex.axis=aaxis*2,labels=FALSE,cex.axis=aaxis*2)
mtext(name_scenar,side=1,line=4.5,at=1:3,cex=amain)
title(tipi[tp],line=2,cex.main=amain*2)
if(tp==1){
mtext("% in the confidence interval",side=2,line=6,at=100*0.5,cex=amain*1.1)
}
}
dev.off()

pdf("/home/cpicoche/Documents/Plankton/SimulatedData/Loop_MARSS/analyse/v2/essai_graphe_signdiff.pdf",width=awidth,height=aheight)
par(mfrow=c(1,3),mar=c(8,15,10,4))
name_scenar=c("S1","S2","S3")
tipi=c("Diag","OffDiag","Covariates")
for (tp in 1:3){
mini=min(sign_diffmp_array_signif[,,tp,,],na.rm=TRUE)
maxi=max(sign_diffpm_array_signif[,,tp,,],na.rm=TRUE)
plot(1:3,apply(sign_diffpm_array_signif[,,tp,2,1],1,mean,na.rm=TRUE),xlim=c(0.7,4),ylim=c(mini,maxi),xlab="",ylab="",cex.main=amain*2,cex.axis=aaxis*2,cex.lab=alab*2,pch=16,cex=apch*2,xaxt="n",col="black")
arrows(1:3,apply(sign_diffpm_array_signif[,,tp,2,1],1,min,na.rm=TRUE),1:3,apply(sign_diffpm_array_signif[,,tp,1,1],1,max,na.rm=TRUE),code=3,length=0.2,angle=90,col='black',lwd=alwd*2.5)
points(c(1.1,2.1,3.1),apply(sign_diffpm_array_signif[,,tp,2,2],1,mean,na.rm=TRUE),col="black",pch=16,cex=apch*2)
arrows(c(1.1,2.1,3.1),apply(sign_diffpm_array_signif[,,tp,2,2],1,min,na.rm=TRUE),c(1.1,2.1,3.1),apply(sign_diffpm_array_signif[,,tp,1,2],1,max,na.rm=TRUE),code=3,length=0.2,angle=90,col='black',lwd=alwd*2.5,lty=2)

points(c(1.2,2.2,3.2),apply(sign_diffmp_array_signif[,,tp,2,1],1,mean,na.rm=TRUE),col="grey",pch=16,cex=apch*2)
arrows(c(1.2,2.2,3.2),apply(sign_diffmp_array_signif[,,tp,2,1],1,min,na.rm=TRUE),c(1.2,2.2,3.2),apply(sign_diffmp_array_signif[,,tp,2,1],1,max,na.rm=TRUE),code=3,length=0.2,angle=90,col='grey',lwd=alwd*2.5,lty=1)

points(c(1.3,2.3,3.3),apply(sign_diffmp_array_signif[,,tp,2,2],1,mean,na.rm=TRUE),col="grey",pch=16,cex=apch*2)
arrows(c(1.3,2.3,3.3),apply(sign_diffmp_array_signif[,,tp,2,2],1,min,na.rm=TRUE),c(1.3,2.3,3.3),apply(sign_diffmp_array_signif[,,tp,2,2],1,max,na.rm=TRUE),code=3,length=0.2,angle=90,col='grey',lwd=alwd*2.5,lty=2)
axis(1,at=1:3,cex.axis=aaxis*2,labels=FALSE,cex.axis=aaxis*2)
mtext(name_scenar,side=1,line=4.5,at=1:3,cex=amain*1.9)
title(tipi[tp],line=2,cex.main=amain*2)
if(tp==1){
mtext("Number",side=2,line=5.5,at=0.6,cex=amain*1.9)
legend("topleft",c("No season + to -","Season + to -","No season - to +"),col=c("black","black","grey"),lty=c(1,2,1),lwd=alwd*2,cex=aleg*1.5)
}
}
dev.off()
 

###############################################################################################################################################
                                  ########CP - Code used to produce Fig. 3 in Barraquand et al. 2017 ########
###############################################################################################################################################

rm(list=ls())
graphics.off()

tab=read.table("arlm_aicc.csv",header=TRUE,sep=",") #Reads the table containing all AICc from ARX(3)
sp=tab$SP

#Limits for the y-axis
yli1=-4
yli2=25

pdf("figure_3.pdf",width=15)
par(mar=c(3,6,1,1))
plot(0,0,xlim=c(0,length(sp)*3),ylim=c(yli1,yli2),t="n",xlab="",ylab=expression(paste(Delta,"AICc=AICc-AICc(full model)",sep="")),xaxt="n",lwd=2.5,cex.lab=2.5,cex.axis=1.75,cex=2)
axis(1,at=seq(2,37,3),lab=sp,cex=2.5,cex.axis=1.75)
abline(h=0,lty=1,lwd=2.5)
for (s in 1:length(sp)){
	physics=tab[s,"P_S_X"]-tab[s,"F_S_N"] #AICc_physics - AICc_full model
	nutrients=tab[s,"N_S_N"]-tab[s,"F_S_N"] #AICc_nutrients - AICc_full model

	rect((s-1)*3+1,0,(s-1)*3+2,physics,col="blue") #Difference between physics and full model
	rect((s-1)*3+2,0,(s-1)*3+3,nutrients,col="orange") #Difference between nutrients and full model
}
	legend("topright",c("Physics only","Nutrients only"),fill=c("blue","orange"),col=c("blue","orange"),bty="n",cex=2) 
dev.off()

###############################################################################################################################################
                                      ########CP - Plot Fig. 5 in Barraquand et al. 2017 using MARSS results ########
###############################################################################################################################################

### Initialize
graphics.off()
rm(list=ls())
library("stringr")

#Graphical parameters
pm=0.1
ab=0.1
fact=1.5
acex=2.5

fac_main=3.0
fac_axis=1.7
fac_lg=1.5
fac_lab=1.5
alwd=2.5
apc=2

pdf("figure_5.pdf",width=17,height=15)

#Load results for Teychan
load(paste("Teychan_physics_pencen.RData",sep=""))
cisv1=eval(parse(text=paste("cisv1=",ls(pattern='cis'),sep="")))
cis=cisv1

sp=dimnames(cis$model$data)[[1]] #Studied species
var=dimnames(cis$call$model$c)[[1]] #Studied covariates
nom=dimnames(cis$par$B)[[1]] #Names of the interactions

#Draw the graph frame
plot(0,0,t='n',xlim=c(0.75,(length(var)+length(sp))),ylim=c(0.75,length(sp)),yaxt="n",xaxt="n",xlab="",ylab="")
axis(2,at=1:length(sp),lab=rev(sp),cex.axis=fac_axis)
axis(1,at=1:(length(sp)+length(var)),lab=c(sp,"Sal.","Irr.","Wind","Seas."),cex.axis=fac_axis)
for (i in 1:length(sp)){
	abline(h=i,lty=2)
}
abline(v=12.75,lty=2)

#get the value from the MARSS object, according to the names of the coefficients
for (n in 1:length(nom)){
	if(grepl("^\\(",nom[n])){ #default when using unconstrained or diagonal, names are of the form (1,2) for (i,j)
		i=as.numeric(strsplit(nom[n],split="[\\(,\\)]")[[1]][2])
		j=as.numeric(strsplit(nom[n],split="[\\(,\\)]")[[1]][3])
		if(is.na(i)){ #i and j might not be numeric
			i=which(unlist(lapply(sp,function(x) grepl(x,strsplit(nom[n],split="[\\(,\\)]")[[1]][2]))))
			j=which(unlist(lapply(sp,function(x) grepl(x,strsplit(nom[n],split="[\\(,\\)]")[[1]][3]))))
		}
	}else{
		a=str_split(nom[n],sp) #The coefficient is XY where X eats Y, X and Y being in the species list
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

	baseline=length(sp)-i+1 #Compute the position of the species on the y-axis
        if(i==j){ #For intragroup interactions, withdraw 1 to the value, to make them comparable to intergroup interactions
		mini=min(0,(cis$par$B[n]-1)*fact) 
		maxi=max(0,(cis$par$B[n]-1)*fact)
		rect(j-pm,baseline+mini,j+pm,baseline+maxi,col="black") #Draw the bar corresponding to the coefficient value
        }else{
		mini=min(0,(cis$par$B[n])*fact)
		maxi=max(0,(cis$par$B[n])*fact)
		rect(j-pm,baseline+mini,j+pm,baseline+maxi,col="black")
        }
        if(!is.na(cis$par.se$B[n])){
                if(cis$par.upCI$B[n]*cis$par.lowCI$B[n]>0){ #If upper and lower values of the confidence intervals have the same sign, the coefficient is deemed significant
			points(j,baseline+maxi+ab,pch='*',cex=acex)
                }
	}
}

#Covariate effects
l=0
for (j in 1:length(var)){
        for (i in 1:length(sp)){
		l=l+1
		ibis=(i-1)*2+1 #Position on the x-axis
		baseline=length(sp)-i+1 #Position on the y-axis
		mini=min(0,(cis$par$U[l])*fact)
		maxi=max(0,(cis$par$U[l])*fact)
		rect(length(sp)+j-pm,baseline+mini,length(sp)+j+pm,baseline+maxi,col="black") #Draw a line
                if(!is.na(cis$par.se$U[l])){
                        if(cis$par.upCI$U[l]*cis$par.lowCI$U[l]>0){
				points(j+length(sp),baseline+maxi+ab,pch='*',col="black",cex=acex)
                        }
                }
        }
}

#Load B7 results and do exactly the same thing, with a small shift on the x-axis
rm(list=ls(pattern='cis'))
load(paste("B7_physics_pencen.RData",sep=""))
cisv1=eval(parse(text=paste("cisv1=",ls(pattern='cis'),sep=""))) 
cis=cisv1

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
	ibis=i*2
	baseline=length(sp)-i+1
        if(i==j){
		mini=min(0,(cis$par$B[n]-1)*fact)
		maxi=max(0,(cis$par$B[n]-1)*fact)
		rect(j+2*pm,baseline+mini,j+4*pm,baseline+maxi,col="blue")
       	}else{
		mini=min(0,(cis$par$B[n])*fact)
		maxi=max(0,(cis$par$B[n])*fact)
		rect(j+2*pm,baseline+mini,j+4*pm,baseline+maxi,col="blue")
       	}
        if(!is.na(cis$par.se$B[n])){
                if(cis$par.upCI$B[n]*cis$par.lowCI$B[n]>0){
			points(3*pm+j,baseline+maxi+ab,pch='*',col="blue",cex=acex)
                }
        }
}

l=0
for (j in 1:length(var)){
        for (i in 1:length(sp)){
                l=l+1
		ibis=i*2
		baseline=length(sp)-i+1
		mini=min(0,(cis$par$U[l])*fact)
		maxi=max(0,(cis$par$U[l])*fact)
		rect(length(sp)+j+2*pm,baseline+mini,length(sp)+j+4*pm,baseline+maxi,col="blue")
                if(!is.na(cis$par.se$U[l])){
                        if(cis$par.upCI$U[l]*cis$par.lowCI$U[l]>0){
				points(3*pm+j+length(sp),baseline+maxi+ab,pch='*',col="blue",cex=acex)
                        }
                }

        }
}
dev.off()

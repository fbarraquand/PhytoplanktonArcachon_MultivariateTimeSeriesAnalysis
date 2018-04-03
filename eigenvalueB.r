#This function computes the eigenvalue of the interaction matrices
rm(list=ls())
graphics.off()
library("stringr")


model=c("pencen","diatdin","unconstrained","inverse")
site=c("Teychan","B7")

for (s in site){
	for (m in model){
		file=paste("MARSS_results/",s,"physics",m,'seasonal_bis_NEW_NOMISTAKE.RData',sep="_")
		print(file)
		load(file)
sp=dimnames(fit_log$model$data)[[1]]
B=matrix(0,length(sp),length(sp))
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
        if(i==j){
                B[i,j]=fit_log$par$B[n]-1
        }else{
                B[i,j]=fit_log$par$B[n]

	}
}



		eig_B=eigen(B)$values
		print(max(Mod(eig_B)))
	}
}
		

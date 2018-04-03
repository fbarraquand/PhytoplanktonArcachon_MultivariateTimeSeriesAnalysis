#This function should take whatever number of estimated values (and sometimes real values for simulated datasets) and return the values of AICc and BIC
rm(list=ls())
graphics.off()

#Define datasets in RData
#For simulated data
name_model=c("null","unconstrained","pencen","diatdin","inverse")
#name_model=c("pencen")
lieu=c("Teychan")

dirin="./MARSS_results/"

for(ll in 1:length(lieu)){
	print(lieu[ll])
	for (model in 1:length(name_model)){
				print(name_model[model])
				load(paste(dirin,lieu[ll],"_physics_",name_model[model],".RData",sep=""))
				bibic=-2*fit_log$logLik+log(fit_log$samp.size/12)*fit_log$num.params
				print("AICc")
				print(fit_log$AICc)
				print("BIC")
				print(bibic)
		}
	}
			

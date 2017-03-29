###############################################################################################################################################
### FB 11/06/2016 - Plankton-like MAR models, with seasonal covariate (24 units/year), for which we try to infer ecological interactions ######
### FB 05/07/2016 - Introducing a loop, setting the effects of variable 2 to both positive and negative values, saving outputs regularly ######
###############################################################################################################################################

### Basic stuff
rm(list=ls())
graphics.off()
set.seed(42) 

### I assume two environmental variables, both related to seasonality + some autocorrelated noise
### Several population dynamics datasets generated for the 10 species in response to their numbers and the environment. 

################### Simulated datasets #####################################################################################################

#1. Null model, no interactions, only response to the environment
### All species respond positively to variable 1
### Then half the other species respond negatively to variable 2 and half positively

#2. Only interactions, mostly negative, some positive, with no effect of the environment. Strong dd intragenus. 

#3. Combination of #1 and #2. 

#############################################################################################################################################

################### SIMULATION ##############################################################################################################
### General parameters
n.time<-500  # Number of years or more generally temporal units
index_time<- 1:n.time
n.species<-10 # Number of species in the community

### Loop on repeats
n.repeats<-10 

for (k.repeats in 1:n.repeats)
{

x<-matrix(0, nrow=n.time,n.species,byrow=TRUE) # Matrix of log-abundances or log-densities
epsilon<-matrix(0, nrow=n.time-1,n.species,byrow=TRUE) # Noise on growth rates
x1<-matrix(1,n.species,byrow=TRUE) # Vector of initial abundances - could be changed
x[1,]<-x1
a0<-runif(n.species,2,3) # Zero-log-density growth rates
seasonality<-2*sin(2*pi*index_time/24)	# must be enough to affect the growth rates
sigma<-matrix(0.5,n.species,byrow=TRUE) # standard error

### Environmental variables
y1noise<-arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=n.time,sd=sqrt(0.5) )
y2noise<-arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=n.time,sd=sqrt(1) ) # more noisy

y1<-seasonality+y1noise
y2<-seasonality+y2noise

cor(y1,y2)

### Interaction matrix 

#1. No interactions 

B1<-matrix(0, ncol=n.species, nrow=n.species) 
for (i in 1:n.species){B1[i,i]<-runif(1,-0.2,0.75)}

B<-B1

#2. effect of environmental covariates
c1=runif(n.species,0.8,1)
c2=runif(n.species,-0.5,0.5)

### Simulation of data

for (t in 1:(n.time-1))
	{
    		epsilon[t,]<-rnorm(n.species,mean = 0, sd = sigma)# Noise
		x[t+1,]<-a0 + y1[t]*c1 + y2[t]*c2 + B %*% x[t,] + epsilon[t,] ### B is defined with different diagonal coefficients
}

### Plotting B matrix
image(B)
image(B-diag(B))

### Equilibrium values - some algebra shows that we have 
x_star<- - solve(B) %*% a0 ### 
x_star
eigen(B)$values ### 

matplot(x)

### Zoom in on the second and third year
matplot(x[25:75,])
matlines(x[25:75,])

DataPlankton1=data.frame(index_time,as.numeric(y1),as.numeric(y2),x)
names(DataPlankton1)
DataPlankton1$Scenario=1
DataPlankton1$Repeat=k.repeats

####################################################

#2. Only interactions
B2<-matrix(runif(n.species*n.species,-0.3,0.1), ncol=n.species) 
for (i in 1:n.species){B2[i,i]<-runif(1,-0.05,0.5)}
B<-B2

### Simulation of data

for (t in 1:(n.time-1))
	{
    		epsilon[t,]<-rnorm(n.species,mean = 0, sd = sigma)# Noise
		x[t+1,]<-a0 + B %*% x[t,] + epsilon[t,] ### B is defined with different diagonal coefficients
}

### Plotting B matrix
image(B)
image(B-diag(B))

### Equilibrium values - some algebra shows that we have 
x_star<- - solve(B) %*% a0 ### 
x_star
eigen(B)$values ### 

matplot(x)

### Zoom in on the second and third year
matplot(x[25:75,])
matlines(x[25:75,])


DataPlankton2=data.frame(index_time,as.numeric(y1),as.numeric(y2),x)
names(DataPlankton2)
DataPlankton2$Scenario=2
DataPlankton2$Repeat=k.repeats

######################################################################################################################################

#3. Both 
B<-B2

for (t in 1:(n.time-1))
	{
    		epsilon[t,]<-rnorm(n.species,mean = 0, sd = sigma)# Noise
		x[t+1,]<-a0 + y1[t]*c1 + y2[t]*c2 + B %*% x[t,] + epsilon[t,] ### B is defined with different diagonal coefficients
}


matplot(x)

### Zoom in on the second and third year
matplot(x[25:75,])
matlines(x[25:75,])

DataPlankton3=data.frame(index_time,as.numeric(y1),as.numeric(y2),x)
names(DataPlankton3)
DataPlankton3$Scenario=3
DataPlankton3$Repeat=k.repeats

######### Save covariables and parameters ##############
save(a0,y1,y2,c1,c2,B1,B2,file=paste("SimulatedParameters_repeat",as.character(k.repeats),"_noseason.RData",sep=""))

########################################
DataPlankton=rbind(DataPlankton1,DataPlankton2,DataPlankton3)
names(DataPlankton)[1:3]=c("Time_index","Abiotic_var1","Abiotic_var2")

### Adding to the overall data structure before saving
if (k.repeats==1){DataPlankton_Loop=DataPlankton #initializing
} else {
DataPlankton_Loop=rbind(DataPlankton_Loop,DataPlankton) #stacking with previous repeat
}

}#end of loop on repeats

### Write to file 

write.csv(DataPlankton_Loop,file="MockPlanktonTimeSeries_Loop_noseason.csv")





# *******************************************************************
# *****	JAGS scripts accompanying 'Bayesian multivariate        *****
# *****	meta-regression: A tutorial'
# *****	Author: Amy Mulick, LSHTM							                  *****
# ***** February 2024                                           *****
# *****                                                         *****
# *****	jcode.1: Pareto hyperprior on precision			            *****
# *****	jcode.2: Gamma hyperprior on precision			            *****
# *****	jcode.3: Uniform hyperprior on log variance		        	*****
# ***** jcode.4: Uniform hyperprior on variance			            *****
# ***** jcode.5: Uniform hyperprior on standard deviation			  *****
# ***** jcode.6: Half normal hyperprior on standard deviation		*****
# *******************************************************************

jcode.1 <- "model{

############## LIKELIHOOD

for(i in 1:N){               #number of studies
  for(j in 1:M){               #number of exposures 
    for(k in 1:P){               #number of outcomes
      
      #within-study outcome distributions
      logor[i,j,k]   ~ dnorm(mu[i,j,k], pow(se[i,j,k],-2))
      
      #between-study mean of data-generating mechanisms
      mu[i,j,k] <- u[i,k] + b[k]*exposure[i,j,k]
      
      #impute missing SEs
      se[i,j,k] ~ dunif(0, 1)
      
    }
  }
}

############## PRIOR DISTRIBUTIONS

for (k in 1:P) {               #number of outcomes 
  b[k] ~ dnorm(0, 1/100000)       #FIXED outcome-specific slopes
  
  for(i in 1:N){               #number of studies
    u[i,k] ~ dnorm(0,inv.tau2[k])    #RANDOM studyXoutcome slope
  }
  
}

############## HYPERPRIORS

#Prior 1 (precision)
for(k in 1:P){              
  inv.tau2[k]     ~ dpar(1, a.tau)
  tau2[k]   <- 1/inv.tau2[k]
}

}"



jcode.2 <- "model{

############## LIKELIHOOD

for(i in 1:N){               #number of studies
  for(j in 1:M){               #number of exposures 
    for(k in 1:P){               #number of outcomes
      
      logor[i,j,k]   ~ dnorm(mu[i,j,k], pow(se[i,j,k],-2))
      mu[i,j,k] <- u[i,k] + b[k]*exposure[i,j,k]
      se[i,j,k] ~ dunif(0, 1)
      
    }
  }
}

############## PRIOR DISTRIBUTIONS

for (k in 1:P) {               #number of outcomes 
  b[k] ~ dnorm(0, 1/100000)       #FIXED outcome-specific slopes
  
  for(i in 1:N){               #number of studies
    u[i,k] ~ dnorm(0,inv.tau2[k])    #RANDOM studyXoutcome effect, does not differ by exposure
  }
  
}

############## HYPERPRIORS

#Prior 2 (precision)
for(k in 1:P){
  inv.tau2[k]     ~ dgamma(a.tau, a.tau)
  tau2[k]   <- 1/inv.tau2[k]
}

}"




jcode.3 <- "model{

############## LIKELIHOOD

for(i in 1:N){               #number of studies
  for(j in 1:M){               #number of exposures 
    for(k in 1:P){               #number of outcomes
      
      logor[i,j,k]   ~ dnorm(mu[i,j,k], pow(se[i,j,k],-2))
      mu[i,j,k] <- u[i,k] + b[k]*exposure[i,j,k]
      se[i,j,k] ~ dunif(0, 1)
      
    }
  }
}

############## PRIOR DISTRIBUTIONS

for (k in 1:P) {               #number of outcomes 
  b[k] ~ dnorm(0, 1/100000)       #FIXED outcome-specific slopes
  
  for(i in 1:N){               #number of studies
    u[i,k] ~ dnorm(0,inv.tau2[k])    #RANDOM studyXoutcome effect, does not differ by exposure
  }
  
}

############## HYPERPRIORS

#Prior 3 (log variance)
for(k in 1:P){
  inv.tau2[k]   <- 1/tau2[k]
  tau2[k] <- exp(log.tau2[k])
  log.tau2[k]     ~ dunif(-10, a.tau)
}

}"



jcode.4 <- "model{

############## LIKELIHOOD

for(i in 1:N){               #number of studies
  for(j in 1:M){               #number of exposures 
    for(k in 1:P){               #number of outcomes
      
      logor[i,j,k]   ~ dnorm(mu[i,j,k], pow(se[i,j,k],-2))
      mu[i,j,k] <- u[i,k] + b[k]*exposure[i,j,k]
      se[i,j,k] ~ dunif(0, 1)
      
    }
  }
}

############## PRIOR DISTRIBUTIONS

for (k in 1:P) {               #number of outcomes 
  b[k] ~ dnorm(0, 1/100000)       #FIXED outcome-specific slopes
  
  for(i in 1:N){               #number of studies
    u[i,k] ~ dnorm(0,inv.tau2[k])    #RANDOM studyXoutcome effect, does not differ by exposure
  }
  
}

############## HYPERPRIORS

#Prior 4 (variance)
for(k in 1:P){
  inv.tau2[k]   <- 1/tau2[k]
  tau2[k]         ~ dunif(0, a.tau)
}

}"



jcode.5 <- "model{

############## LIKELIHOOD

for(i in 1:N){               #number of studies
  for(j in 1:M){               #number of exposures 
    for(k in 1:P){               #number of outcomes
      
      logor[i,j,k]   ~ dnorm(mu[i,j,k], pow(se[i,j,k],-2))
      mu[i,j,k] <- u[i,k] + b[k]*exposure[i,j,k]
      se[i,j,k] ~ dunif(0, 1)
      
    }
  }
}

############## PRIOR DISTRIBUTIONS

for (k in 1:P) {               #number of outcomes 
  b[k] ~ dnorm(0, 1/100000)       #FIXED outcome-specific slopes
  
  for(i in 1:N){               #number of studies
    u[i,k] ~ dnorm(0,inv.tau2[k])    #RANDOM studyXoutcome effect, does not differ by exposure
  }
  
}

############## HYPERPRIORS

#Prior 5 (standard deviation)
for(k in 1:P){
  inv.tau2[k]   <- 1/tau2[k]
  tau2[k] <- tau[k]^2
  tau[k]          ~ dunif(0, a.tau)
}

}"




jcode.6 <- "model{

############## LIKELIHOOD

for(i in 1:N){               #number of studies
  for(j in 1:M){               #number of exposures 
    for(k in 1:P){               #number of outcomes
      
      logor[i,j,k]   ~ dnorm(mu[i,j,k], pow(se[i,j,k],-2))
      mu[i,j,k] <- u[i,k] + b[k]*exposure[i,j,k]
      se[i,j,k] ~ dunif(0, 1)
      
    }
  }
}

############## PRIOR DISTRIBUTIONS

for (k in 1:P) {               #number of outcomes 
  b[k] ~ dnorm(0, 1/100000)       #FIXED outcome-specific slopes
  
  for(i in 1:N){               #number of studies
    u[i,k] ~ dnorm(0,inv.tau2[k])    #RANDOM studyXoutcome effect, does not differ by exposure
  }
  
}

############## HYPERPRIORS

#Prior 6 (standard deviation)
for(k in 1:P){
  inv.tau2[k]   <- 1/tau2[k]
  tau2[k] <- tau[k]^2
  tau[k]          ~ dnorm(0, a.tau) T(0,)
}

}"


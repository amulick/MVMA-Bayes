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
      
      logor[i,j,k]   ~ dnorm(mu[i,j,k], pow(se[i,j,k],-2))          #outcome distribution
      mu[i,j,k] <- u[i,k] + b1[k] + b2[k]*exposure[i,j,k]      #data generating model
      
      se[i,j,k] ~ dunif(0, imp.se)        #impute missing SEs (does not overwrite observed SEs)
      
    }
  }
}

############## PRIOR DISTRIBUTIONS

# fixed effects
for (k in 1:P) {               #number of outcomes 
  b1[k] ~ dnorm(0, 1/100000)       #outcome-specific intercept
  b2[k] ~ dnorm(0, 1/100000)       #outcome-specific slope 
}

# random effects
for(i in 1:N){               #number of studies
  for(k in 1:P){               #number of outcomes
    u[i,k] ~ dnorm(0,inv.tau2[k])    #random studyXoutcome effect, does not differ by exposure
  }
}

############## HYPERPRIORS

#Prior 1 (precision)
for(j in 1:P){              
  inv.tau2[j]     ~ dpar(1, a.tau)
  tau2[j]   <- 1/inv.tau2[j]
}

}"



jcode.2 <- "model{

############## LIKELIHOOD

for(i in 1:N){               #number of studies
  for(j in 1:M){               #number of exposures 
    for(k in 1:P){               #number of outcomes
      
      logor[i,j,k]   ~ dnorm(mu[i,j,k], pow(se[i,j,k],-2))          #outcome distribution
      mu[i,j,k] <- u[i,k] + b1[k] + b2[k]*exposure[i,j,k]      #data generating model
      
      se[i,j,k] ~ dunif(0, imp.se)        #impute missing SEs (does not overwrite observed SEs)
      
    }
  }
}

############## PRIOR DISTRIBUTIONS

# fixed effects
for (k in 1:P) {               #number of outcomes 
  b1[k] ~ dnorm(0, 1/100000)       #outcome-specific intercept
  b2[k] ~ dnorm(0, 1/100000)       #outcome-specific slope 
}

# random effects
for(i in 1:N){               #number of studies
  for(k in 1:P){               #number of outcomes
    u[i,k] ~ dnorm(0,inv.tau2[k])    #random studyXoutcome effect
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
      
      logor[i,j,k]   ~ dnorm(mu[i,j,k], pow(se[i,j,k],-2))          #outcome distribution
      mu[i,j,k] <- u[i,k] + b1[k] + b2[k]*exposure[i,j,k]      #data generating model
      
      se[i,j,k] ~ dunif(0, imp.se)        #impute missing SEs (does not overwrite observed SEs)
      
    }
  }
}

############## PRIOR DISTRIBUTIONS

# fixed effects
for (k in 1:P) {               #number of outcomes 
  b1[k] ~ dnorm(0, 1/100000)       #outcome-specific intercept
  b2[k] ~ dnorm(0, 1/100000)       #outcome-specific slope 
}

# random effects
for(i in 1:N){               #number of studies
  for(k in 1:P){               #number of outcomes
    u[i,k] ~ dnorm(0,inv.tau2[k])    #random studyXoutcome effect
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
      
      logor[i,j,k]   ~ dnorm(mu[i,j,k], pow(se[i,j,k],-2))          #outcome distribution
      mu[i,j,k] <- u[i,k] + b1[k] + b2[k]*exposure[i,j,k]      #data generating model
      
      se[i,j,k] ~ dunif(0, imp.se)        #impute missing SEs (does not overwrite observed SEs)
      
    }
  }
}

############## PRIOR DISTRIBUTIONS

# fixed effects
for (k in 1:P) {               #number of outcomes 
  b1[k] ~ dnorm(0, 1/100000)       #outcome-specific intercept
  b2[k] ~ dnorm(0, 1/100000)       #outcome-specific slope 
}

# random effects
for(i in 1:N){               #number of studies
  for(k in 1:P){               #number of outcomes
    u[i,k] ~ dnorm(0,inv.tau2[k])    #random studyXoutcome effect
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
      
      logor[i,j,k]   ~ dnorm(mu[i,j,k], pow(se[i,j,k],-2))          #outcome distribution
      mu[i,j,k] <- u[i,k] + b1[k] + b2[k]*exposure[i,j,k]      #data generating model
      
      se[i,j,k] ~ dunif(0, imp.se)        #impute missing SEs (does not overwrite observed SEs)
      
    }
  }
}

############## PRIOR DISTRIBUTIONS

# fixed effects
for (k in 1:P) {               #number of outcomes 
  b1[k] ~ dnorm(0, 1/100000)       #outcome-specific intercept
  b2[k] ~ dnorm(0, 1/100000)       #outcome-specific slope 
}

# random effects
for(i in 1:N){               #number of studies
  for(k in 1:P){               #number of outcomes
    u[i,k] ~ dnorm(0,inv.tau2[k])    #random studyXoutcome effect
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
      
      logor[i,j,k]   ~ dnorm(mu[i,j,k], pow(se[i,j,k],-2))          #outcome distribution
      mu[i,j,k] <- u[i,k] + b1[k] + b2[k]*exposure[i,j,k]      #data generating model
      
      se[i,j,k] ~ dunif(0, imp.se)        #impute missing SEs (does not overwrite observed SEs)
      
    }
  }
}

############## PRIOR DISTRIBUTIONS

# fixed effects
for (k in 1:P) {               #number of outcomes 
  b1[k] ~ dnorm(0, 1/100000)       #outcome-specific intercept
  b2[k] ~ dnorm(0, 1/100000)       #outcome-specific slope 
}

# random effects
for(i in 1:N){               #number of studies
  for(k in 1:P){               #number of outcomes
    u[i,k] ~ dnorm(0,inv.tau2[k])    #random studyXoutcome effect
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


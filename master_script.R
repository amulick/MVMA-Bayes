# *******************************************************************
# *****	Analysis scripts accompanying 'Bayesian multivariate    *****
# *****	meta-regression: A tutorial'
# *****	Author: Amy Mulick, LSHTM							                  *****
# ***** February 2024                                           *****
# *****                                                         *****
# *******************************************************************

rm(list=ls())

#### Load packages, data & scripts
library("rjags")
library("mixmeta")

load(file="../data/ecz_sev_MVMA_data.RData")
d <- ecz_sev_MVMA_data

source("JAGS_scripts.R") 



###########################################
######  prepare data for Bayesian analysis ######
###########################################

# ensure correct order of dataset
d <- d[order(d$outcome,d$exposure,d$study),]

# JAGS input parameters
jdat <- list(N=length(table(d$study)), 
             M=length(table(d$exposure)), 
             P=length(table(d$outcome)),
             imp.se=0.4, a.tau=0.005)
jdat$logor=array(d$logor, dim = c(jdat$N, jdat$M, jdat$P))
jdat$se=array(d$se, dim = c(jdat$N, jdat$M, jdat$P))
jdat$exposure=array(as.numeric(as.factor(d$exposure)), dim=c(jdat$N, jdat$M, jdat$P))

# matrix of hyperprior limits
hp.scale <- matrix(c(0.001, 0.25, 0.001, 0.1, 10, 2, 1000, 4, 100, 2, 100, 2),
                   nrow = 6, byrow = T)

# parameters to collect in output
to.collect <- c("b1", "b2","tau2")



###########################################
######  Bayesian analyses  ######
###########################################

# models
for (y in 1:dim(hp.scale)[1]) {        #loop over different hyperpriors
  for (z in 1:dim(hp.scale)[2]){       #loop over different limits
    
    jdat$a.tau <- hp.scale[y, z]  
    set.seed(1000399)
    jmod <- jags.model(textConnection(get(paste0("jcode.", y))), data=jdat, n.chains=2, n.adapt=1000)
    update(jmod,5000)
    
    nam <- paste("jpos", y, z, sep=".")
    assign(nam, coda.samples(jmod, to.collect, n.iter=100000, thin=5))
    
  }
}

rm(nam)



###########################################
######  graphs    ######
###########################################

# diagnostic trace plots: Figure XX
for (y in 1:dim(hp.scale)[1]) {        #loop over different hyperpriors
  for (z in 1:dim(hp.scale)[2]){       #loop over different limits
    
    par(mfrow=c(5,4), mar=c(1,1,4,1), oma=c(0,0,2,0))
    traceplot(get(paste("jpos", y, z, sep=".")))
    title(main=paste("Bayesian model, hyperprior", y, 
                     ifelse(z==1, "(vague)", "(lightly informative)")),
          outer=T)
    
  }
}

# raw data: Figure XX
{
  plot(NULL, axes = F, xlab = "", ylab = "Relative Risk", ylim = log(c(0.5, 5)), 
       xlim = c(0,3.1), main = "Observed effects of eczema severity on CVD")
  abline(h=0)
  axis(1, at = 0:3, labels = c("None", levels(d$exposure)))
  axis(2, at = log(c(0.5,0.8,1,2,5)), labels = c(0.5,0.8,1,2,5), las = 2)
  
  c <- 2
  for (i in levels(d$outcome)){
    for (j in levels(d$study)) {
      
      d.t <- d[d$outcome==i & d$study==j, ]
      d.t$exposure <- as.numeric(d.t$exposure)
      lines(c(0, d.t$exposure + (c-3)/80), c(0, d.t$logor), 
            lwd=2, type = "b", col=c)
      segments(d.t$exposure + (c-3)/60, d.t$logor - 1.96*d.t$se,
               d.t$exposure + (c-3)/60, d.t$logor + 1.96*d.t$se, 
               lwd=2, col=c)
    }
    c <- c+1
  }
  
  legend("topleft",levels(d$outcome),col=2:7,lty=1,ncol=2,pch=1, lwd=2)
  rm(d.t, c)
}

# forest plots: Figure XX
{
  d.t <- data.frame()
  for (y in 1:dim(hp.scale)[1]) {        #loop over different hyperpriors
    for (z in 1:dim(hp.scale)[2]){       #loop over different limits
      
      m <- summary(get(paste("jpos", y, z, sep=".")))
      dd <- data.frame(cbind(m[["statistics"]][grep("b2", rownames(m$statistics)),],
                             m[["quantiles"]][grep("b2", rownames(m$quantiles)),]))
      dd$pr <- y
      dd$v.i <- z
      d.t <- rbind(d.t, dd)
      
    }
  }
  rm(dd, m)
  d.t$CrI <- paste0("(",
                    format(exp(d.t$X2.5.),digits=2, nsmall=2), ", ",
                    format(exp(d.t$X97.5.),digits=2, nsmall=2), ")")  
  d.t$desc <- paste0("HP ", d.t$pr, ifelse(d.t$v.i==1, ": vague", ": lightly inf."))
  
  titles <- levels(d$outcome)
  
  for (y in 1:length(titles)) {
    e <- paste0("b2.", y, ".")
    dd <- d.t[grep(e, rownames(d.t)), ]
    
    plot(NA, axes=F, xlim=log(c(0.75, 2)), ylim=c(1,12), 
         main= titles[y], xlab = "Hazard ratio", ylab = "")
    axis(1, at=log(c(0.75, 1, 1.5, 2)), labels=paste(c(0.75, 1, 1.5, 2)))
    abline(v=0)
    text(log(2), 12:1, paste(format(exp(dd[, "Mean"]), digits=3), dd[, "CrI"]), adj=1)
    mtext(dd[, "desc"] , side=2, at=12:1, las=2, adj=0, line=2)
    mtext("--> Eczema severity increases risk", side=1, at=log(1.01), cex=0.75, adj=0)
    mtext("Eczema severity decreases risk <--", side=1, at=log(0.99), cex=0.75, adj=1, line=1.5)
    
    points(dd[, "Mean"], 12:1, pch = 15)
    segments(dd[, "X2.5."], 12:1, dd[, "X97.5."], 12:1, )
  }
  b.res <- d.t[d.t$pr==6 & d.t$v.i==1, c("Mean", "X2.5.", "X97.5.", "CrI")]  #save best estimates
  rm(dd, d.t)
}

# frequentist analysis for comparison: Figure XX
{
  m <- list()
  j=1
  for (i in levels(d$outcome)) {
    m[[j]] <- mixmeta(logor ~ as.numeric(exposure), se, random = ~ 1|study,
                      data=d[d$outcome == i,])[c("coefficients", "vcov")]
    j <- j+1
  }
  
  d.t <- data.frame(effect=unlist(m)[grep("exposure", names(unlist(m)))], 
                    se=sqrt(unlist(m)[grep("vcov4", names(unlist(m)))]))
  d.t$lb <- d.t$effect - 1.96*d.t$se
  d.t$ub <- d.t$effect + 1.96*d.t$se
  d.t$CI <- paste0(" (", 
                   format(exp(d.t$lb), digits = 2), ", ",
                   format(exp(d.t$ub), digits = 3), ")")
  d.t$outcome <- levels(d$outcome)
  
  plot(NA, axes=F, xlim=log(c(0.75, 2)), ylim=c(0,6), 
       main= "Univariate frequentist meta-regression", xlab = "Hazard ratio", ylab = "")
  axis(1, at=log(c(0.75, 1, 1.5, 2)), labels=paste(c(0.75, 1, 1.5, 2)))
  abline(v=0)
  
  text(log(2), 6:1, paste(format(exp(d.t[, "effect"]), digits=3), d.t[, "CI"]), adj=1)
  mtext(d.t[, "outcome"] , side=2, at=6:1, las=2, adj=0, line=2)
  mtext("--> Eczema severity increases risk", side=1, at=log(1.01), cex=0.75, adj=0)
  mtext("Eczema severity decreases risk <--", side=1, at=log(0.99), cex=0.75, adj=1, line=1.5)
  
  points(d.t[, "effect"], 6:1, pch = 15)
  segments(d.t[, "lb"], 6:1, d.t[, "ub"], 6:1, )
  points(b.res[, "Mean"], 6:1-0.25, pch = 15, col="gray")
  segments(b.res[, "X2.5."], 6:1-0.25, b.res[, "X97.5."], 6:1-0.25, col="gray")
  text(log(2), 6:1-0.25, paste(format(exp(b.res[, "Mean"]), digits=3), b.res[, "CrI"]),
       adj=1, col="gray")
  
  rm(d.t, m)
}

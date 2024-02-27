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

source("../scripts/JAGS_scripts.R") 



###########################################
######  prepare data for Bayesian analysis ######
###########################################

# ensure correct order of dataset
d <- d[order(d$outcome,d$exposure,d$study),]

# JAGS input parameters
jdat <- list(N=length(table(d$study)), 
             M=length(table(d$exposure)), 
             P=length(table(d$outcome)),
             a.tau=0.005)
jdat$logor=array(d$logor, dim = c(jdat$N, jdat$M, jdat$P))
jdat$se=array(d$se, dim = c(jdat$N, jdat$M, jdat$P))
jdat$exposure=array(as.numeric(as.factor(d$exposure)), dim=c(jdat$N, jdat$M, jdat$P))

# matrix of hyperprior limits
hp.scale <- matrix(c(0.001, 0.25, 0.001, 0.1, 10, 2, 1000, 4, 100, 2, 1/100, 1/2),
                   nrow = 6, byrow = T)

# parameters to collect in output
to.collect <- c("b","tau2")



###########################################
######  Bayesian analyses  ######
###########################################

# models with random seeds set for each chain for reproducible results
for (y in 1:dim(hp.scale)[1]) {        #loop over different hyperpriors
  for (z in 1:dim(hp.scale)[2]){       #loop over different limits
    
    jdat$a.tau <- hp.scale[y, z]  
    set.seed(1000399)
    jmod <- jags.model(textConnection(get(paste0("jcode.", y))), data=jdat, 
                       n.chains=2, n.adapt=1000, 
                       inits=list(list(.RNG.name = "base::Wichmann-Hill",
                                       .RNG.seed = 123456),
                                  list(.RNG.name = "base::Wichmann-Hill",
                                       .RNG.seed = 7890)))
    update(jmod,5000)
    
    nam <- paste("jpos", y, z, sep=".")
    assign(nam, coda.samples(jmod, to.collect, n.iter=100000, thin=5))
    
  }
}

rm(nam)



###########################################
######  graphs    ######
###########################################

# diagnostic plots
for (y in 1:dim(hp.scale)[1]) {        #loop over different hyperpriors
  for (z in 1:dim(hp.scale)[2]){       #loop over different limits
    
    jpos <- get(paste("jpos", y, z, sep="."))
    tit <- paste("Bayesian model, hyperprior", y, 
                 ifelse(z==1, "(vague)", "(lightly informative)"))
    
    # transform to analysis scale if relevant. HP1&2 were precision, 
    # 3&4 variance (no transformation needed), 5&6 SD
    for (i in 1:length(jpos)) {
      if (y==1 | y==2) jpos[[i]][, grep("tau2\\[", colnames(jpos[[i]]))] <- 
          apply(jpos[[i]][, grep("tau2\\[", colnames(jpos[[i]]))], 1:2, function(x) 1/x)
      if (y==5 | y==6) jpos[[i]][, grep("tau2\\[", colnames(jpos[[i]]))] <- 
          apply(jpos[[i]][, grep("tau2\\[", colnames(jpos[[i]]))], 1:2, function(x) sqrt(x))
    }
    
    par(mfrow=c(4,3), mar=c(1,2,4,1), oma=c(0,0,2,0))
    traceplot(jpos)
    title(main=tit, outer=T)
    gelman.plot(jpos, auto.layout = F, ylim=c(1, 1.5))
    title(main=tit, outer=T)
    
  }
}

# Figure 1
{
  par(mfrow=c(1,1), mar=c(4,2,3,1))
  plot(NULL, axes = F, xlab = "", ylab = "Relative Risk", ylim = log(c(0.5, 5)), 
       xlim = c(0,3.1), main = "Observed effects of eczema on CVD, by severity")
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

# Figure 2
{  # create dataset
  d.t <- data.frame()
  for (y in 1:dim(hp.scale)[1]) {        #loop over different hyperpriors
    for (z in 1:dim(hp.scale)[2]){       #loop over different limits
      
      m <- summary(get(paste("jpos", y, z, sep=".")))
      dd <- data.frame(cbind(m[["statistics"]], m[["quantiles"]]))
      dd$pr <- y
      dd$v.i <- z
      d.t <- rbind(d.t, dd)
      
    }
  }
  rm(dd, m)
  d.t$desc <- paste0("HP ", d.t$pr, ifelse(d.t$v.i==1, ": vague", ": lightly inf."))
  titles <- levels(d$outcome)
  
  # graph labels
  par(mfcol = c(2,7), mar=c(4,2,2,1), oma=c(1,1,1,3))
  
  plot(NA, axes=F, xlim=log(c(0.75, 2)), ylim=c(1,13), xlab = "", ylab = "")
  mtext(d.t$desc[grep("b\\[1.", rownames(d.t))], side=2, 
        at=0.2+12:1, las=2, adj=0, line=-5, cex=0.9)
  
  plot(NA, axes=F, xlim=log(c(0.75, 2)), ylim=c(1,13), xlab = "", ylab = "")
  mtext(d.t$desc[grep("b\\[1.", rownames(d.t))], side=2, 
        at=0.2+12:1, las=2, adj=0, line=-5, cex=0.9)
  
  #graph results
  for (y in 1:length(titles)) {
    e <- paste0("b\\[", y, ".")
    f <- paste0("tau2\\[", y, ".")
    b <- d.t[grep(e, rownames(d.t)), ]
    t <- d.t[grep(f, rownames(d.t)), ]
    b$CrI <- paste0(format(exp(b$X50.),   digits=3), " (",
                    format(exp(b$X2.5.),  digits=2, nsmall=2),", ",
                    format(exp(b$X97.5.), digits=2, nsmall=2),")")  
    
    t$CrI <- paste0(format(sqrt(t$X50.),   digits=3), " (",
                    format(sqrt(t$X2.5.),  digits=2, nsmall=2),", ",
                    format(sqrt(t$X97.5.), digits=2, nsmall=2),")")  
    # betas
    plot(NA, axes=F, xlim=log(c(0.75, 2)), ylim=c(1,13), main=titles[y],
         xlab = "Hazard ratio", ylab = "")
    axis(1, at=log(c(0.75, 1, 1.5, 2)), labels=paste(c(0.75, 1, 1.5, 2)))
    abline(v=0)
    
    text(log(2.6), 0.2+c(12.6,12:1), c("Median (95% CrI)", b[, "CrI"]), adj=1,
         font=c(2, rep(1,12)), cex=c(1, rep(0.95,12)), xpd=NA)
    points(b[, "X50."], 12:1, pch = 15)
    segments(b[, "X2.5."], 12:1, b[, "X97.5."], 12:1)
    
    # taus
    plot(NA, axes=F, xlim=sqrt(c(0, 4)), ylim=c(1,13), 
         xlab = "Between-study SD \n(units in log(HR))", ylab = "")
    axis(1, at=sqrt(c(0, 0.25, 1, 4)), labels=paste(sqrt(c(0, 0.25, 1, 4))))
    
    text(sqrt(6), 0.2+c(12.6,12:1), c("Median (95% CrI)", t[, "CrI"]), adj=1,
         font=c(2, rep(1,12)), cex=c(1, rep(0.95,12)), xpd=NA)
    points(sqrt(t[, "X50."]), 12:1, pch = 15)
    segments(sqrt(t[, "X2.5."]), 12:1, sqrt(t[, "X97.5."]), 12:1)
  }

  #save best estimates  & delete unnecessary objects
  b.res <- d.t[d.t$pr==3 & d.t$v.i==2, c("X50.", "X2.5.", "X97.5.")]
  b.res <- b.res[grep("b\\[", rownames(b.res)),]
  rm(d.t, b, t, e, f)
}

# Figure 3
{
  m <- list()
  j=1
  for (i in levels(d$outcome)) {
    m[[j]] <- mixmeta(logor ~ as.numeric(exposure) - 1, se, random = ~ 1|study,
                      data=d[d$outcome == i,])[c("coefficients", "vcov")]
    j <- j+1
  }
  
  d.t <- data.frame(effect=unlist(m)[grep("exposure", names(unlist(m)))], 
                    se=sqrt(unlist(m)[grep("vcov", names(unlist(m)))]))
  d.t$lb <- d.t$effect - 1.96*d.t$se
  d.t$ub <- d.t$effect + 1.96*d.t$se
  d.t$CI <- paste0(format(exp(d.t$effect), digits = 3), " (", 
                   format(exp(d.t$lb),     digits = 2), ", ",
                   format(exp(d.t$ub),     digits = 3), ")")
  d.t$outcome <- levels(d$outcome)
  b.res$CrI <- paste0(format(exp(b.res[, "X50."]), digits=3), " (", 
                      format(exp(b.res[, "X2.5."]), digits=3), ", ",
                      format(exp(b.res[, "X97.5."]), digits=3), ")")
  
  par(mfrow=c(1,1), mar=c(4,1,2,1), oma=c(1,1,1,1))
  plot(NA, axes=F, xlim=log(c(0.75, 2)), ylim=c(0,6), 
       main= "Frequentist vs. Bayesian comparison", xlab = "Hazard ratio", ylab = "")
  axis(1, at=log(c(0.75, 1, 1.5, 2)), labels=paste(c(0.75, 1, 1.5, 2)))
  abline(v=0)
  
  text(log(2), 6:1, d.t[, "CI"], adj=1)
  mtext(d.t[, "outcome"] , side=2, at=6:1, las=2, adj=0, line=2)
  mtext("--> More severe eczema increases risk", side=1, at=log(1.01), cex=0.75, adj=0)
  
  points(d.t[, "effect"], 6:1, pch = 15)
  segments(d.t[, "lb"], 6:1, d.t[, "ub"], 6:1, )
  points(b.res[, "X50."], 6:1-0.25, pch = 15, col="gray")
  segments(b.res[, "X2.5."], 6:1-0.25, b.res[, "X97.5."], 6:1-0.25, col="gray")
  text(log(2), 6:1-0.25, b.res$CrI, adj=1, col="gray")
  
  rm(d.t, m)
}



###########################################
######  multivariate frequentist attempt    ######
###########################################


d <- d[order(d$outcome, d$study, d$exposure),]
dd <- reshape(d[, c("study", "outcome", "exposure", "logor", "se")], 
              idvar = c("study","exposure"), timevar = "outcome", direction = "wide") 

model <- mixmeta(as.matrix(dd[, grep("logor", names(dd))]) ~ as.numeric(dd$exposure) -1, 
                 S=as.matrix(dd[, grep("se[.]", names(dd))]), 
                 random= ~ 1|dd$study, bscov = "diag")

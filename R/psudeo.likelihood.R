pseudo.neg.log.like <- function(par, lengths, k, delta, V, data, weights, outcome.model=c("PO","ACAT")){

  in.par <- par
  if(outcome.model == "PO" & k > 2){
    in.par[1:(k-1)] <- cumsum(par[1:(k-1)])
    in.par[lengths[1] + lengths[2] + (1:(k-1))] <- cumsum(par[lengths[1] + lengths[2] + (1:(k-1))])
    in.par[lengths[1] + lengths[2] + lengths[3] + (1:(k-1))] <- cumsum(par[lengths[1] + lengths[2] + lengths[3] + (1:(k-1))])
  }
  
  par.list <- list(a = in.par[1:lengths[1]],
                   e = in.par[(lengths[1]+1):(lengths[1]+lengths[2])],
                   b = in.par[(lengths[1]+lengths[2]+1):(lengths[1]+lengths[2]+lengths[3])],
                   c = in.par[(lengths[1]+lengths[2]+lengths[3]+1):sum(lengths)])
  
  if (lengths[2]==0) par.list$e <- NULL

  v.probs <- estimate.outcome.probabilities(par=par.list, data=data, delta=delta, k=k, outcome.model=outcome.model)

  #get the probability of the specific outcome
  p.alpha <- v.probs$D0[cbind(1:nrow(data), data[,V])]
  p.beta.gamma <- v.probs$D1[cbind(1:nrow(data), data[,V])]
  
  part.A.1 <- sum(log(p.beta.gamma[data[delta]==1]))
  #print(part.A.1)
  part.A.2 <- sum(log(((weights$w1 * p.beta.gamma) + (weights$w0 * p.alpha))[data[delta]==0]))
  #print(part.A.2)
  return(-(part.A.1 + part.A.2))
}

old.pseudo.neg.log.like <- function(par, dat, cox.weights, rangeZ, rangeT, rangeR, rangeTR, rangeV, V.model)
{
  v.alpha0 <- par[c("V.alpha1", "V.alpha2")] 
  v.beta0 <- par[c("V.beta1", "V.beta2")]
  v.gamma0 <- par[c("V.gamma1", "V.gamma2")]
  v.alphaZ <- par[c("V.alphaZ1", "V.alphaZ2")]
  v.gammaZ <- par[c("V.gammaZ1", "V.gammaZ2")]
  v.alphaR <- par["V.alphaR"] 
  v.betaR <- par["V.betaR"]
  v.gammaR <- par["V.gammaR"]
  v.gammaT = par["V.gammaT"]
  v.gammaTR = par["V.gammaTR"]
  
  if(V.model=="PO"){
    v.alpha0["V.alpha2"] <- v.alpha0["V.alpha2"]+v.alpha0["V.alpha1"]
    v.beta0["V.beta2"] <- v.beta0["V.beta2"]+v.beta0["V.beta1"]
    v.gamma0["V.gamma2"] <- v.gamma0["V.gamma2"]+v.gamma0["V.gamma1"]
    v.probs <- probs.po(alpha0 = v.alpha0, alphaZ = v.alphaZ, alphaR = v.alphaR, 
                      beta0 = v.beta0, gammaZ = v.gammaZ,
                      gamma0 = v.gamma0, gammaR = v.gammaR,
                      gammaT = v.gammaT, gammaTR = v.gammaTR,
                      rangeZ=rangeZ, rangeT=rangeT, rangeR=rangeR,
                      rangeTR=rangeTR, dat = as.matrix(dat[-rangeV]))
  }
  if(V.model=="ACAT"){
    v.probs <- probs.acat(alpha0 = v.alpha0, alphaZ = v.alphaZ, alphaR = v.alphaR, 
                        beta0 = v.beta0, gammaZ = v.gammaZ,
                        gamma0 = v.gamma0, gammaR = v.gammaR,
                        gammaT = v.gammaT, gammaTR = v.gammaTR,
                        rangeZ=rangeZ, rangeT=rangeT, rangeR=rangeR,
                        rangeTR=rangeTR, dat = as.matrix(dat[-rangeV]))
  }
  
  #get the probability of the specific outcome
  p.alpha <- v.probs$D0[cbind(1:nrow(dat), dat$V)]
  p.beta.gamma <- v.probs$D1[cbind(1:nrow(dat), dat$V)]
  
  part.A.1 <- sum(log(p.beta.gamma[dat$delta==1]))
  part.A.2 <- sum(log(((cox.weights$w1 * p.beta.gamma) + (cox.weights$w0 * p.alpha))[dat$delta==0]))
  
  return(-(part.A.1 + part.A.2))
}

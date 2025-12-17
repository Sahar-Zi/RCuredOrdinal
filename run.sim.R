run.sim <- function(seed = 2212308, n, k, par, outcome.model = c("PO","ACAT"),
                    formula.c, formula.a, equal.effect, Tau, R, delta,
                    replications = 1, B = 0, alpha = 0.05, save.iter = FALSE){
  
  timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  
  summary.sim <- function(est.values, true.par, est.se.values, CI.between){
    avg <- colMeans(est.values)
    true.d <- true.par$d
    true.tau <- c(true.par$Tau)
    true.a <- true.par$a
    true.b <- true.par$b
    true.c <- true.par$c
    bias <- c(true.d, true.tau, true.a, true.b, true.c) - avg
    SD <- colStdevs(est.values)
    names(bias) <- names(avg)
    if(!is.na(est.se.values[1,1])){
      SE <- colMeans(est.se.values)
      CP <- colMeans(CI.between)
      return(data.frame(avg, bias, SD, SE, CP))
    }
    
    return(data.frame(avg, bias, SD))
  }
  
  # ---- Check the formulas only for the simulation, there are the same checks in the ordercuredfit
  if(is.null(formula.a) || is.null(formula.c) 
     || identical(as.character(formula.a), "") || identical(as.character(formula.c), "")){stop("one or both of the models is empty")}
  
  if(all.vars(formula.a)[1] != all.vars(formula.c)[1]){stop("response must be the same for both models")}
  
  if(delta %in% c(all.vars(formula.a), all.vars(formula.c))){stop("delta cannot be in the models")}
  
  if(Tau %in% all.vars(formula.a)){stop("Tau cannot be in model.a, it is not defined for the uncured fraction")}
  # ----
  
  p <- length(unlist(par)) - 2 # Would change to -1 after using the weibull cure model
  seeds <- data.frame(matrix(NA, nrow = replications, ncol = B + 1))
  replication.errors <- 0
  
  est.values <- matrix(NA, nrow = replications, ncol = p)
  naive.est.values <- matrix(NA, nrow = replications, ncol = k + 2)
  cc.est.values <- matrix(NA, nrow = replications, ncol = k + 4)
  
  est.se.values <- matrix(NA, nrow = replications, ncol = p)
  naive.est.se.values <- matrix(NA, nrow = replications, ncol = k + 2)
  cc.est.se.values <- matrix(NA, nrow = replications, ncol = k + 4)
  
  CIL <- matrix(NA, nrow = replications, ncol = p)
  CIU <- matrix(NA, nrow = replications, ncol = p)
  CI.between <- matrix(NA, nrow = replications, ncol = p)
  
  rep.i <- 1
  seed.c <- seed
  
  if (save.iter) {
    log_file <- paste0("LOG-Simulation-", "-", timestamp, "-rep=", replications,
                       "-n=", n, "-B=", B, "-v.model=", outcome.model, ".txt")
  }
  
  cat("\n----------------------------------------------\n")
  while(rep.i < replications + 1){
    
    log_output <- capture.output({
      cat("\n----------------------------------------------\n")
      cat(paste("Replication number", rep.i, "\n"))
    })
    
    bootstrap.error.counts <- integer(replications)
    bootstrap.errors <- 0
    
    seeds[rep.i,1] <- seed.c
    set.seed(seed.c)
    seed.c <- seed.c + 1
    cat(paste("Replication number", rep.i, "is Running...\n"))
    
    dat <- gen.demo.data(n = n, k = k, par = par, outcome.model = outcome.model)
    
    est <- tryCatch({
      ordcuredfit(formula.c = formula.c, formula.a = formula.a, equal.effect = equal.effect,
                  Tau=Tau, R=R, delta=delta, data=dat, outcome.model=outcome.model)},
      error = function(e) {
        message("Error at replication ", rep.i, ": ", e$message) 
        NULL  # Return NA or another placeholder value
      })
    if (is.null(est)) next
    
    est.values[rep.i,] <- unlist(est$par)
    naive.est.values[rep.i,] <- unlist(est$naive[,"Estimate"])
    cc.est.values[rep.i,] <- unlist(est$cc[,"Estimate"])
    naive.est.se.values[rep.i,] <- unlist(est$naive[,"Std. Error"])
    cc.est.se.values[rep.i,] <- unlist(est$cc[,"Std. Error"])
    
    colnames(est.values) <- colnames(est.se.values) <- colnames(CIL) <- colnames(CIU) <- colnames(CI.between) <- names(unlist(est$par))
    colnames(naive.est.values) <- colnames(naive.est.se.values) <- names(unlist(est$naive[,"Estimate"]))
    colnames(cc.est.values) <- colnames(cc.est.se.values) <- names(unlist(est$cc[,"Estimate"]))
    
    if (!is.null(est$trace)) {
      log_output <- c(log_output,
                      "Optimization Trace:", est$trace)
    }
    
    log_output <- c(log_output,
                    "Estimates:", capture.output(print(est.values[rep.i,])),
                    "Naive Estimates:", capture.output(print(naive.est.values[rep.i,])),
                    "Complete Case Estimates:", capture.output(print(cc.est.values[rep.i,])),"\n")
    
    if (save.iter) {
      cat(log_output, file = log_file, sep = "\n", append = TRUE)
    } else {
      cat(log_output, sep = "\n")
    }
    
    if (B>0){
      est.values.BS <- matrix(NA, nrow = B, ncol = p)
      rep.b <- 1
      
      while(rep.b < B + 1){
        # Set bootstrap seed
        set.seed(seed.c)
        seeds[rep.i, rep.b + 1] <- seed.c
        seed.c <- seed.c + 1
        
        # Generate bootstrap sample
        boot.dat <- dat[sample(n, size=n, replace=TRUE),]
        
        # Analyze bootstrap sample
        boot.est <- tryCatch({
          ordcuredfit(formula.c = formula.c, formula.a = formula.a, equal.effect = equal.effect,
                      Tau=Tau, R=R, delta=delta, data=boot.dat, outcome.model="PO")},
          error = function(e) {
            message("Error at replication ", i, " Bootstrap ",rep.b, ": ", e$message) 
            bootstrap.errors <<- bootstrap.errors + 1
            NULL
          })
        if (is.null(boot.est)) {
          replication.errors <- replication.errors + 1
          next
        }
        
        ## LOGGING THE BOOTSTRAP TRACE IF save.iter == TRUE
        if (save.iter && !is.null(boot.est$trace)) {
          bs_header <- paste0("---- BOOTSTRAP ", rep.b, " OF REPLICATION ", rep.i, " ----")
          cat(c(bs_header, boot.est$trace, ""), file = log_file, sep = "\n", append = TRUE)
        }
        
        est.values.BS[rep.b,] <- unlist(boot.est$par)
        rep.b <- rep.b + 1
      }
      
      est.se.values[rep.i,] <- colStdevs(est.values.BS)
      
      for (j in 1:ncol(est.values.BS)){
        ci_bs <- quantile(est.values.BS[,j],c(alpha/2,1-(alpha/2)))
        CIL[rep.i,j] <- ci_bs[1]
        CIU[rep.i,j] <- ci_bs[2]
        CI.between[rep.i,j] <- between(unlist(par)[-c(4,5)][j], ci_bs[1], ci_bs[2]) # minus scale and shape
      }
      
      if (save.iter) {
        cat(paste0("Bootstrap errors in replication ", rep.i, ": ", bootstrap.errors, "\n"),
            file = log_file, append = TRUE)
      }
      bootstrap.error.counts[rep.i] <- bootstrap.errors
    }
    
    rep.i <- rep.i + 1
    
    cat("\n----------------------------------------------\n")
  }
  
  summary.table <- summary.sim(est.values, par, est.se.values, CI.between)
  
  if (save.iter) {
    cat("\n========== Simulation Summary ==========\n", file = log_file, append = TRUE)
    cat(paste("Total failed replications:", replication.errors, "\n"), file = log_file, append = TRUE)
    cat(paste("Average bootstrap errors per replication:", mean(bootstrap.error.counts), "\n"),
        file = log_file, append = TRUE)
    cat("========================================\n", file = log_file, append = TRUE)
    
    write_out <- function(timestamp, data, name) {
      file_conn <- file(paste0("Simulation-", name, "-", timestamp, "-rep=", replications,
                               "-n=", n, "-B=", B, "-v.model=", outcome.model, ".txt"), "w")
      suppressWarnings(
        write.table(data, file = file_conn, sep = "\t", na = "", row.names = FALSE, col.names = TRUE, append = TRUE)
      )
      close(file_conn)
    }
    
    write_out(timestamp, est.values, "est.vals")
    write_out(timestamp, est.se.values, "est.se.vals")
    write_out(timestamp, naive.est.values, "naive.est.vals")
    write_out(timestamp, cc.est.values, "cc.est.vals")
    write_out(timestamp, naive.est.se.values, "naive.est.se.vals")
    write_out(timestamp, cc.est.se.values, "cc.est.se.vals")
    write_out(timestamp, CIL, "CIL")
    write_out(timestamp, CIU, "CIU")
    write_out(timestamp, CI.between, "CI.between")
    write_out(timestamp, summary.table, "summary.table")
  }
  
  return(list(summary.table = summary.table,
              true.params = par, est.values = est.values, est.se.values = est.se.values,
              CIL = CIL, CIU = CIU, CI.between = CI.between,
              naive.est.values = naive.est.values, cc.est.values = cc.est.values,
              naive.est.se.values = naive.est.se.values, cc.est.se.values = cc.est.se.values, 
              seed = seeds))
}

run.sim.param <- function(seed = 2212308, n, k, par, outcome.model = c("PO","ACAT"),
                    survform, cureform, formula.a, formula.e, formula.b, formula.c, 
                    equal.effect, Tau, R, delta, var = var,
                    replications = 1, alpha = 0.05, save.iter = FALSE){
  
  timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  
  CI <- function(est.values, est.se.values, alpha, true.par){
    l <- est.values - qnorm(1 - alpha / 2) * est.se.values
    u <- est.values + qnorm(1 - alpha / 2) * est.se.values
    return((true.par >= l) & (true.par <= u))
  }
  
  summary.sim <- function(est.values, true.par, est.se.values, CI.est, est.se.inv.values = NULL, CI.est.inv = NULL){
    names(true.par) <- names(est.values)
    avg <- colMeans(est.values)
    bias <- true.par - avg
    SD <- colStdevs(est.values)
    names(bias) <- names(avg)
    if(!is.na(est.se.values[1,1])){
      SE <- colMeans(est.se.values)
      CP <- colMeans(CI.est)
      if(!is.null(est.se.inv.values)){
        SE.inv <- colMeans(est.se.inv.values)
        CP.inv <- colMeans(CI.est.inv)
        return(data.frame(true = true.par, avg, bias, SD, 
                          "SE-sandwich" = SE, "SE-inv" = SE.inv, "CP-sandwich" = CP, "CP-inv" = CP.inv))
      }
      return(data.frame(true = true.par, avg, bias, SD, SE, CP))
    }
    
    return(data.frame(true = true.par, avg, bias, SD))
  }
  
  # ---- Check the formulas only for the simulation, there are the same checks in the ordercuredfit
  if(is.null(formula.a) || is.null(formula.c) 
     || identical(as.character(formula.a), "") || identical(as.character(formula.c), "")){stop("one or both of the models is empty")}
  
  if(all.vars(formula.a)[1] != all.vars(formula.c)[1]){stop("response must be the same for both models")}
  
  if(delta %in% c(all.vars(formula.a), all.vars(formula.c))){stop("delta cannot be in the models")}
  
  if(Tau %in% all.vars(formula.a)){stop("Tau cannot be in model.a, it is not defined for the uncured fraction")}
  # ----
  
  p <- length(unlist(par)) # Would change to -1 after using the weibull cure model
  seeds <<- data.frame(matrix(NA, nrow = replications, ncol = 1))
  colnames(seeds) <<- c("Seed")
  replication.errors <- 0
  
  est.values <<- est.val <- matrix(NA, nrow = replications, ncol = p)
  naive.est.values <<- naive.est.val <- matrix(NA, nrow = replications, ncol = k + 2)
  cc.est.values <<- cc.est.val <- matrix(NA, nrow = replications, ncol = k + 4)
  
  est.cov.mat <<- array(NA, dim = c(p, p, reps))
  est.cov.inv.mat <<- array(NA, dim = c(p, p, reps))
  est.se.values <<- est.se.val <- matrix(NA, nrow = replications, ncol = p)
  est.se.inv.values <<- est.se.inv.val <- matrix(NA, nrow = replications, ncol = p)
  naive.est.se.values <<- naive.est.se.val <- matrix(NA, nrow = replications, ncol = k + 2)
  cc.est.se.values <<- cc.est.se.val <- matrix(NA, nrow = replications, ncol = k + 4)

  CI.est <<- ci.est <- matrix(NA, nrow = replications, ncol = p)
  CI.est.inv <<- ci.est.inv <- matrix(NA, nrow = replications, ncol = p)
  CI.naive <<- ci.naive <- matrix(NA, nrow = replications, ncol = k+2)
  CI.cc <<- ci.cc <- matrix(NA, nrow = replications, ncol = k+4)
  
  rep.i <- 1
  seed.c <- seed
  
  if (save.iter) {
    log_file <- paste0("LOG-Simulation-", "-", timestamp, "-rep=", replications,
                       "-n=", n, "-v.model=", outcome.model, ".txt")
  }
  
  cat("\n----------------------------------------------\n")
  while(rep.i < replications + 1){
    
    log_output <- capture.output({
      cat("\n----------------------------------------------\n")
      cat(paste("Replication number", rep.i, "\n"))
    })
    
    seeds[rep.i,1] <<- seed.c
    set.seed(seed.c)
    seed.c <- seed.c + 1
    
    cat(paste("Replication number", rep.i, "is Running...\n"))

    dat <- gen.demo.data(n = n, k = k, par = par, outcome.model = outcome.model)
    
    est <- tryCatch({
      ordcuredfit.param(survform = survform, cureform = cureform, formula.a = formula.a, formula.e = formula.e,
                        formula.b = formula.b, formula.c = formula.c, equal.effect = equal.effect,
                  Tau = Tau, R = R, delta = delta, data = dat, outcome.model = outcome.model, var = var)},
      error = function(e) {
        message("Error at replication ", rep.i, ": ", e$message) 
        NULL  # Return NA or another placeholder value
      })
    if (is.null(est)) next
    if (!is.null(est$trace)) {
      log_output <- c(log_output,
                      "Optimization Trace:", est$trace)
    }

    est.values[rep.i,] <<- est.val[rep.i,] <- unlist(est$par)
    colnames(est.values) <<- colnames(est.val) <- names(unlist(est$par))
    naive.est.values[rep.i,] <<- naive.est.val[rep.i,] <- unlist(est$naive[,"Estimate"])
    colnames(naive.est.values) <<- colnames(naive.est.val) <- names(unlist(est$naive[,"Estimate"]))
    cc.est.values[rep.i,] <<- cc.est.val[rep.i,] <- unlist(est$cc[,"Estimate"])
    colnames(cc.est.values) <<- colnames(cc.est.val) <- names(unlist(est$cc[,"Estimate"]))
    
    log_output <- c(log_output,
                    "Estimates:", capture.output(print(est.val[rep.i,])),
                    "Naive Estimates:", capture.output(print(naive.est.val[rep.i,])),
                    "Complete Case Estimates:", capture.output(print(cc.est.val[rep.i,])),"\n")
    
    if (save.iter) {
      cat(log_output, file = log_file, sep = "\n", append = TRUE)
    } else {
      cat(log_output, sep = "\n")
    }
    
    if(var){
      naive.est.se.values[rep.i,] <<- naive.est.se.val[rep.i,] <- unlist(est$naive[,"Std. Error"])
      cc.est.se.values[rep.i,] <<- cc.est.se.val[rep.i,] <- unlist(est$cc[,"Std. Error"])
      
      colnames(est.se.values) <<- colnames(est.se.val) <- names(unlist(est$par))
      colnames(est.se.inv.values) <<- colnames(est.se.inv.val) <- names(unlist(est$par))
      colnames(naive.est.se.values) <<- colnames(naive.est.se.val) <- names(unlist(est$naive[,"Estimate"]))
      colnames(cc.est.se.values) <<- colnames(cc.est.se.val) <- names(unlist(est$cc[,"Estimate"]))
      
      #est.se.values[rep.i,] <<- est.se.val[rep.i,] <- runif(21,0,1)
      #est.se.inv.values[rep.i,] <<- est.se.inv.val[rep.i,] <- runif(21,0,2)
      
      est.se.values[rep.i,] <<- est.se.val[rep.i,] <- sqrt(diag(est.cov.mat[,,rep.i] <- est$variance$stacked.v.est))
      est.se.inv.values[rep.i,] <<- est.se.inv.val[rep.i,] <- sqrt(diag(est.cov.inv.mat[,,rep.i] <- -est$variance$G.tilde.inv))
      
      #v.est <- variance.est(est = est, outcome.model = outcome.model, data = as.data.frame(dat))
      #est.se.values[rep.i,] <<- est.se.val[rep.i,] <- sqrt(diag(est.cov.matrix[,,rep.i] <<- v.est$stacked.v.est))
      
      CI.est[rep.i,] <<- ci.est[rep.i,] <- CI(est.val[rep.i,], est.se.val[rep.i,], alpha, unlist(par))
      CI.est.inv[rep.i,] <<- ci.est.inv[rep.i,] <- CI(est.val[rep.i,], est.se.inv.val[rep.i,], alpha, unlist(par))
      CI.naive[rep.i,] <<- ci.naive[rep.i,] <- CI(naive.est.val[rep.i,], naive.est.se.val[rep.i,], alpha, unlist(par$a))
      CI.cc[rep.i,] <<- ci.cc[rep.i,] <- CI(cc.est.val[rep.i,], cc.est.se.val[rep.i,], alpha, 
                                            c(unlist(par$c)[1:(k-1)],unlist(par$e),unlist(par$c)[k:length(par$c)]))
    }
    rep.i <- rep.i + 1
    
    save.image(paste("Simulation-Environment-", timestamp, "-rep=", replications,
                     "-n=", n, "-v.model=", outcome.model, ".RData", sep=""))
    cat("\n----------------------------------------------\n")
  }
  
  summary.table.est <- summary.sim(est.val, unlist(par), est.se.val, ci.est, est.se.inv.val, ci.est.inv)
  summary.table.naive <- summary.sim(naive.est.val, unlist(par$a), naive.est.se.val, ci.naive)
  summary.table.cc <- summary.sim(cc.est.val, c(unlist(par$c)[1:(k-1)],unlist(par$e),unlist(par$c)[k:length(par$c)]), cc.est.se.val, ci.cc)
  
  if (save.iter) {
    cat("\n========== Simulation Summary ==========\n", file = log_file, append = TRUE)
    cat(paste("Total failed replications:", replication.errors, "\n"), file = log_file, append = TRUE)
    cat("========================================\n", file = log_file, append = TRUE)
    
    write_out <- function(timestamp, data, name) {
      file_conn <- file(paste0("Simulation-", name, "-", timestamp, "-rep=", replications,
                               "-n=", n, "-v.model=", outcome.model, ".txt"), "w")
      suppressWarnings(
        write.table(data, file = file_conn, sep = "\t", na = "", row.names = TRUE, col.names = TRUE, append = TRUE)
      )
      close(file_conn)
    }
    
    write_out(timestamp, est.val, "est.vals")
    write_out(timestamp, naive.est.val, "naive.est.vals")
    write_out(timestamp, cc.est.val, "cc.est.vals")
    if(var){
      write_out(timestamp, est.se.val, "est.se.vals")
      write_out(timestamp, est.se.inv.val, "est.se.inv.val")
      write_out(timestamp, ci.est, "ci.est")
    }
    write_out(timestamp, summary.table.est, "summary.table.est")
    write_out(timestamp, summary.table.naive, "summary.table.naive")
    write_out(timestamp, summary.table.cc, "summary.table.cc")
    
    save.image(paste("Simulation-Environment-", timestamp, "-rep=", replications,
                     "-n=", n, "-v.model=", outcome.model, ".RData", sep=""))
  }
  
  if(var){
    return(list(summary.table.est = summary.table.est, summary.table.naive = summary.table.naive, summary.table.cc = summary.table.cc,
                true.params = par, est.values = est.val, 
                est.se.values = est.se.val, est.se.inv.values = est.se.inv.val,
                est.cov.mat = est.cov.mat, est.cov.inv.mat = est.cov.inv.mat,
                CI.est = ci.est, CI.est.inv = ci.est.inv,
                naive.est.values = naive.est.val, cc.est.values = cc.est.val,
                naive.est.se.values = naive.est.se.val, cc.est.se.values = cc.est.se.val,
                CI.naive = ci.naive, CI.cc = ci.cc,
                seed = seeds))
  }
  return(list(summary.table.est = summary.table.est, summary.table.naive = summary.table.naive, summary.table.cc = summary.table.cc,
              true.params = par, est.values = est.val,
              naive.est.values = naive.est.val, cc.est.values = cc.est.val,
              seed = seeds))
}


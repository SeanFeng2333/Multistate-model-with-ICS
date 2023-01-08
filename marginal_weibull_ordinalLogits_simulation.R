#### load packages ####
library(here)
library(tidyverse)
library(tinytex)
library(knitr)
library(kableExtra)
library(msm)
library(minqa)
library(nlme)
library(car)
library(contrast)
library(lme4)
library(data.table)
library(distr)
library(survival)
library(parallel)


#### batch simulation arguments ####

args <- (commandArgs(TRUE))
print(args)

if (length(args)==0){
  print("No arguments supplied.")
  
  ## supply default values
  baseseed <- 2345
  B = 1000 # number of resamples 
  m = 500 # number of clusters 
  baseline_haz12 <- 0.3 # baseline hazard of transition 1-2
  baseline_haz23 <- 0.3 # baseline hazard of transition 2-3
  beta2_12 <- 0.5 # effect of MetS on transition 1-2
  beta2_23 <- 0.5 # effect of MetS on transition 2-3
  ita1 <- -2 # intercept for first cluster size group
  ita2 <- -1 # intercept for second cluster size group
  ita_x <- 0.8 # effect of MetS on cluster size group
  ita_u <- 1 # degree of the ICS
  visittimes <- 5 # number of visits 
  gamma <- 1 # shape parameter of the Weibull distribution
  tau <- 1 # 1/tau is the variance of the frailty
  
  ncores <- 8 # number of cores used for simulation
  numsim <- 2 # number of simulation for each core
} else{
  for (i in (1:length(args))) eval(parse(text=args[[i]]))
}


simtotal <- ncores * numsim # total number of simulation


task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
simgroup_min <- (task_id-1) * simtotal 
simgroup_max <- task_id * simtotal

# simulation function
simulation.func <- function(x){
  from <- seq(1, simtotal, by = numsim)[x]
  to <- seq(numsim, simtotal, by = numsim)[x]
  
  # create empty matrices for results
  naive.exp.est.mat <- matrix(NA, numsim, 6)
  naive.exp.se.mat <- matrix(NA, numsim, 6)
  naive.weibull.est.mat <- matrix(NA, numsim, 6)
  naive.weibull.se.mat <- matrix(NA, numsim, 6)
  gee.weibull.est.mat <- matrix(NA, numsim, 6)
  gee.weibull.se.mat <- matrix(NA, numsim, 6)
  wcr.weibull.est.mat <- matrix(NA, numsim, 6)
  wcr.weibull.se.mat <- matrix(NA, numsim, 6)
  wsf.weibull.est.mat <- matrix(NA, numsim, 6)
  wsf.weibull.se.mat <- matrix(NA, numsim, 6)
  wsf.weibull.pw.est.mat <- matrix(NA, numsim, 11)
  wsf.weibull.pw.se.mat <- matrix(NA, numsim, 11)
  
  for (s in 1:numsim){
    
    # set unique seed for each iteration
    sim <- from + s-1 + simgroup_min
    simseed <- from + s-1 + simgroup_min
    set.seed(simseed)
    
    out <- tryCatch({  
      #### Generate ICS dataset ####
      # generate frailty
      w <- rgamma(m, tau, tau)
      frailty <- w
      # generate clusters with ICS
      cluster_size_grp <- c()
      cluster_size <- c()
      x2 <- rbinom(m, size=1, prob=0.5)
      for (i in 1:m) {
        p1 <- exp(ita1+ita_x*x2[i]+ita_u*log(frailty[i]))/(1+exp(ita1+ita_x*x2[i]+ita_u*log(frailty[i])))
        p2 <- exp(ita2+ita_x*x2[i]+ita_u*log(frailty[i]))/(1+exp(ita2+ita_x*x2[i]+ita_u*log(frailty[i]))) - 
          exp(ita1+ita_x*x2[i]+ita_u*log(frailty[i]))/(1+exp(ita1+ita_x*x2[i]+ita_u*log(frailty[i])))
        p3 <- 1 - exp(ita2+ita_x*x2[i]+ita_u*log(frailty[i]))/(1+exp(ita2+ita_x*x2[i]+ita_u*log(frailty[i])))
        random.grp <- rmultinom(1,1,c(p1, p2, p3))
        cluster_size_grp[i] <- 1*random.grp[1] + 2*random.grp[2] + 3*random.grp[3]
        if (cluster_size_grp[i]==1) {cluster_size[i]=ceiling(runif(1, 5, 15))}
        else if (cluster_size_grp[i]==2) {cluster_size[i]=ceiling(runif(1, 15, 25))}
        else if (cluster_size_grp[i]==3) {cluster_size[i]=ceiling(runif(1, 25, 28))}
      }
      
      cluster_size_vec <- rep(cluster_size, cluster_size)
      cluster_size_grp_vec <- rep(cluster_size_grp, cluster_size)
      cluster_id <- rep(1:m, cluster_size)
      obs_num <- length(cluster_id)
      
      # create a vector of individual based on the sameple size of cluster
      individual <- c()
      for (i in 1:m){
        individual <- c(individual, 1:cluster_size[i])
      }
      cluster_individual <- cluster_id + individual * 0.01
      
      # create a frailty vector
      frailty_vec <- rep(frailty, cluster_size)
      
      ## generate time to event from exponential distribution 
      # The survival times base on drawing from a U~unif(0,1), U is the survival probability
      u12 <- runif(obs_num, 0, 1)
      u23 <- runif(obs_num, 0, 1)
      # generate covariates in cluster level
      x0 <- rep(1, obs_num)
      x2 <- rbinom(m, size=1, prob=0.5)
      x2 <- rep(x2, cluster_size)
      
      # exact transition times
      exact.t.12 <- ( -log(1-u12) / (baseline_haz12 * frailty_vec * exp(beta2_12 * x2 )))^(1/gamma)
      exact.t.23 <- ( -log(1-u23) / (baseline_haz23 * frailty_vec * exp(beta2_23 * x2 )) + exact.t.12^gamma)^(1/gamma) 
      
      #### visit times ####
      # vtime
      vtime.1 <- runif(m, .01, .01)
      for(i in 2:visittimes) { 
        last_nam <- paste("vtime.", i-1, sep = "")
        nam <- paste("vtime.", i, sep = "")
        assign(nam, runif(m, 2, 2)+ get(last_nam))
      }
      
      # repeat visit times for all tooth
      for(i in 1:visittimes) { 
        nam <- paste("vtime.", i, sep = "")
        assign(nam, rep(get(nam), cluster_size))
      }
      
      # visit states
      for(i in 1:visittimes) { 
        nam <- paste("vstate.", i, sep = "")
        assign(nam, case_when(exact.t.12 > get(paste("vtime.",i ,sep = "")) ~ 1,
                              exact.t.12 <= get(paste("vtime.",i ,sep = "")) & exact.t.23 > get(paste("vtime.",i ,sep = "")) ~ 2,
                              exact.t.12 < get(paste("vtime.",i ,sep = "")) & exact.t.23 <= get(paste("vtime.",i ,sep = "")) ~ 3))
      }
      
      # create visit data
      visit.time.dat<-c()
      visit.state.dat<-c()
      time.name <- c()
      state.name <- c()
      for(i in 1:visittimes) { 
        visit.time.dat <- cbind(visit.time.dat, get(paste("vtime.", i, sep = "")))
        visit.state.dat <- cbind(visit.state.dat, get(paste("vstate.", i, sep = "")))
        time.name <- c(time.name, paste("vtime.", i, sep = ""))
        state.name <- c(state.name, paste("vstate.", i, sep = ""))
      }
      colnames(visit.time.dat) <- time.name
      colnames(visit.state.dat) <- state.name
      
      # create dataset
      sim.data <- data.frame(cluster_id, individual, cluster_individual, frailty_vec, cluster_size_vec,
                             exact.t.12, exact.t.23, x0, x1, x2, visit.time.dat, visit.state.dat)

      # create long format dataset:
      state <- sim.data %>% 
        pivot_longer(
          cols = starts_with("vstate."),
          values_to = "state"
        ) %>%
        dplyr::select(cluster_individual, state)
      
      sim.data.long <- sim.data %>% 
        pivot_longer(
          cols = starts_with("vtime."),
          names_to = "visit",
          values_to = "visit.time"
        ) %>%
        cbind(state = state$state) %>% 
        dplyr::select(-starts_with("vstate.")) 
      
      sim.data.long.1 <- sim.data.long %>% filter(!(state == 3))
      
      # keep the first teeth loss
      sim.data.long.2 <- sim.data.long %>% filter(state == 3) %>% 
        group_by(cluster_id, individual) %>% 
        mutate(state = ifelse(row_number() == 1, state, NA)) %>% 
        na.omit()
      
      sim.data.long <- rbind(sim.data.long.1, sim.data.long.2) %>% 
        arrange(cluster_id, individual) %>% 
        group_by(cluster_id, individual) 
      
      # remove individuals with only one visit 
      sim.data.long.1visit <- sim.data.long %>% 
        group_by(cluster_id, individual) %>% 
        mutate(obsnum =dplyr::n()) %>% 
        filter((obsnum==1)) 
      n_distinct(sim.data.long.1visit$cluster_id)
      
      sim.data.long <- sim.data.long %>% 
        mutate(obsnum=dplyr::n())%>% 
        filter(!(obsnum==1)) 
      
      sim.data.long <- sim.data.long %>% 
        group_by(cluster_id, individual) %>% 
        mutate(transit=c(state[1], state[1:((dplyr::n())-1)] + 0.1*state[2:(dplyr::n())]),
               time.diff=visit.time - lag(visit.time),
               t1=lag(visit.time),
               t2=visit.time,
               transit2.t2=case_when(transit >=2.0 ~ time.diff))
      
      sim.data.long[is.na(sim.data.long)] <- 0
      
      sim.data.wsf <- sim.data.long %>% filter(!(transit == 1 | transit==2))
      
      # two datasets for piecewise WSF 
      sim.data.wsf1 <- sim.data.wsf %>% filter(t1 < (visittimes-1)*2/2)
      sim.data.wsf2 <- sim.data.wsf %>% filter(t1 >= (visittimes-1)*2/2)
      
      #### naive ####
      # specify initial intensity matrix 
      twoway3.q <- rbind(c(0, 0.1, 0), c(0, 0, 0.3), c(0, 0, 0))
      rownames(twoway3.q) <- colnames(twoway3.q) <- c("state1", "state2","state3")
      
      naive.model <- msm(state ~ visit.time, 
                         covariates=~x2, 
                         center = F,
                         subject = cluster_individual, data = sim.data.long, 
                         opt.method="optim",method="BFGS",  
                         hessian=T,
                         qmatrix=twoway3.q)
      
      naive.est <- naive.model$estimates
      naive.cov <- naive.model$covmat
      exp.naive.est <- naive.est
      exp.naive.se <- sqrt(diag(naive.cov))
      
      #### naive weibull ####
      inital.para <- c(log(1), as.vector(naive.est))
      
      logliklihood <- function(inital.para, data){
        k <- inital.para[1]
        beta120 <- inital.para[2]
        beta121 <- inital.para[4]
        beta230 <- inital.para[3]
        beta231 <- inital.para[5]
        q12 <- with(data, exp(beta120 * x0 + beta121 * x2))
        q23 <- with(data, exp(beta230 * x0 + beta231 * x2))
        p11 <- with(data, exp(-q12*(t2^(exp(k)) - t1^(exp(k)))))
        p22 <- with(data, exp(-q23*(t2^(exp(k)) - t1^(exp(k)))))
        p12 <- with(data, q12/(q12-q23) * (p22-p11))
        p23 <- with(data, 1- p22)
        p13 <- with(data, 1 - p11 - p12)
        p13[p13<=0] = 10^-10
        llk <- with(data, -2*sum((1*(transit==1.1)*log(p11)  + 
                                    1*(transit==1.2)*log(p12) + 
                                    1*(transit==1.3)*log(p13) +
                                    1*(transit==2.2)*log(p22) + 
                                    1*(transit==2.3)*log(p23))))
        return(llk)
      }
      naive.weibull.opt <- optim(inital.para, logliklihood, data = sim.data.wsf, method = "BFGS", hessian = T)
      weibull.naive.est <-  naive.weibull.opt$par
      naive.weibull.cov <- solve(0.5*naive.weibull.opt$hessian)
      naive.weibull.cor <- cor(naive.weibull.cov)    
      weibull.naive.se <- sqrt(diag(naive.weibull.cov))
      
      #### regular GEE ####
      set.seed(simseed)
      inital.para <- c(log(1), as.vector(naive.est))
      logliklihood <- function(inital.para, data){
        gamma <- exp(inital.para[1])
        beta120 <- (inital.para[2])
        beta121 <- inital.para[4]
        beta230 <- (inital.para[3])
        beta231 <- inital.para[5]
        q12 <- with(data, exp(beta120+beta121*x2))
        q23 <- with(data, exp(beta230+beta231*x2))
        p11 <- with(data, exp(-q12*(t2^gamma - t1^gamma)))
        p22 <- with(data, exp(-q23*(t2^gamma - t1^gamma)))
        p12 <- with(data, q12/(q12-q23) * (p22-p11))
        p23 <- with(data, 1- p22)
        p13 <- with(data, 1 - p11 - p12)
        p13[p13<=0] = 10^-10
        llk <- with(data, -2*sum((1*(transit==1.1)*log(p11)  + 
                                    1*(transit==1.2)*log(p12) + 
                                    1*(transit==1.3)*log(p13) +
                                    1*(transit==2.2)*log(p22) + 
                                    1*(transit==2.3)*log(p23))))
        return(llk)
      }
      
      weibull.optim.est <- optim(inital.para, logliklihood, data = sim.data.wsf, method = "BFGS", hessian = T, control = list(trace=1))
      weibull.gee.est <-  weibull.optim.est$par
      H <- weibull.optim.est$hessian
      sqrt(diag(solve(H)))
      # covariance matrix 
      k <- weibull.optim.est$par[1]
      beta120 <- weibull.optim.est$par[2]
      beta230 <- weibull.optim.est$par[3]
      beta121 <- weibull.optim.est$par[4]
      beta231 <- weibull.optim.est$par[5]
      
      score_k <- c()
      score_beta120 <- c()
      score_beta230 <- c()
      score_beta121 <- c()
      score_beta231 <- c()
      for (i in 1:nrow(sim.data.wsf)){
        x0 <- sim.data.wsf$x0[i]
        x <- sim.data.wsf$x2[i]
        t1 <- sim.data.wsf$t1[i]
        t2 <- sim.data.wsf$t2[i]
        transit <- sim.data.wsf$transit[i]
        cluster.size <- sim.data.wsf$cluster_size_vec[i]
        
        q12 <- exp(beta120 * x0 + beta121 * x)
        q23 <- exp(beta230 * x0 + beta231 * x)
        y <- q12/(q12-q23)
        
        ##### transition probability #####
        p11t <- exp(-q12 * (t2^exp(k) - t1^exp(k)))
        p22t <- exp(-q23 * (t2^exp(k) - t1^exp(k)))
        p12t <- (q12/(q12-q23)) * (p22t - p11t)
        p13t <- 1-p11t-p12t
        p23t <- 1-p22t
        
        ##### first derivative of p #####
        dp11t_dk <- -p11t*q12*exp(k)*(t2^(exp(k))*log(t2) - t1^(exp(k))*log(t1))
        dp11t_dbeta120 <- -p11t*q12*(t2^exp(k) - t1^exp(k))*x0
        dp11t_dbeta230 <- 0
        dp11t_dbeta121 <- -p11t*q12*(t2^exp(k) - t1^exp(k))*x
        dp11t_dbeta231 <- 0
        
        dp22t_dk <- -p22t*q23*exp(k)*(t2^(exp(k))*log(t2) - t1^(exp(k))*log(t1))
        dp22t_dbeta120 <- 0
        dp22t_dbeta230 <- -p22t*q23*(t2^exp(k) - t1^exp(k))*x0
        dp22t_dbeta121 <- 0
        dp22t_dbeta231 <- -p22t*q23*(t2^exp(k) - t1^exp(k))*x
        
        dy_dk <- 0
        dy_dbeta120 <- -q12*q23*x0/(q12-q23)^2
        dy_dbeta230 <- q12*q23*x0/(q12-q23)^2
        dy_dbeta121 <- -q12*q23*x/(q12-q23)^2
        dy_dbeta231 <- q12*q23*x/(q12-q23)^2
        
        dp12t_dk <- y*(dp22t_dk - dp11t_dk)
        dp12t_dbeta120 <- dy_dbeta120*(p22t-p11t) + y*(dp22t_dbeta120-dp11t_dbeta120)
        dp12t_dbeta230 <- dy_dbeta230*(p22t-p11t) + y*(dp22t_dbeta230-dp11t_dbeta230)
        dp12t_dbeta121 <- dy_dbeta121*(p22t-p11t) + y*(dp22t_dbeta121-dp11t_dbeta121)
        dp12t_dbeta231 <- dy_dbeta231*(p22t-p11t) + y*(dp22t_dbeta231-dp11t_dbeta231)
        
        dp13t_dk <- -dp11t_dk - dp12t_dk
        dp13t_dbeta120 <- -dp11t_dbeta120 - dp12t_dbeta120
        dp13t_dbeta230 <- -dp11t_dbeta230 - dp12t_dbeta230
        dp13t_dbeta121 <- -dp11t_dbeta121 - dp12t_dbeta121
        dp13t_dbeta231 <- -dp11t_dbeta231 - dp12t_dbeta231
        
        dp23t_dk <- -dp22t_dk
        dp23t_dbeta120 <- -dp22t_dbeta120
        dp23t_dbeta230 <- -dp22t_dbeta230
        dp23t_dbeta121 <- -dp22t_dbeta121
        dp23t_dbeta231 <- -dp22t_dbeta231
        
        ##### 5*1 Vector of score function #####
        score_k <- c(score_k, (1*(transit==1.1) * dp11t_dk/p11t + 
                                 1*(transit==1.2) * dp12t_dk/p12t + 
                                 1*(transit==1.3) * dp13t_dk/p13t + 
                                 1*(transit==2.2) * dp22t_dk/p22t + 
                                 1*(transit==2.3) * dp23t_dk/p23t) )
        
        score_beta120 <- c(score_beta120, (1*(transit==1.1) * dp11t_dbeta120/p11t + 
                                             1*(transit==1.2) * dp12t_dbeta120/p12t + 
                                             1*(transit==1.3) * dp13t_dbeta120/p13t + 
                                             1*(transit==2.2) * dp22t_dbeta120/p22t + 
                                             1*(transit==2.3) * dp23t_dbeta120/p23t))
        
        score_beta230 <- c(score_beta230, (1*(transit==1.1) * dp11t_dbeta230/p11t + 
                                             1*(transit==1.2) * dp12t_dbeta230/p12t + 
                                             1*(transit==1.3) * dp13t_dbeta230/p13t + 
                                             1*(transit==2.2) * dp22t_dbeta230/p22t + 
                                             1*(transit==2.3) * dp23t_dbeta230/p23t) )
        
        score_beta121 <- c(score_beta121, (1*(transit==1.1) * dp11t_dbeta121/p11t + 
                                             1*(transit==1.2) * dp12t_dbeta121/p12t + 
                                             1*(transit==1.3) * dp13t_dbeta121/p13t + 
                                             1*(transit==2.2) * dp22t_dbeta121/p22t + 
                                             1*(transit==2.3) * dp23t_dbeta121/p23t) )
        
        score_beta231 <- c(score_beta231, (1*(transit==1.1) * dp11t_dbeta231/p11t + 
                                             1*(transit==1.2) * dp12t_dbeta231/p12t + 
                                             1*(transit==1.3) * dp13t_dbeta231/p13t + 
                                             1*(transit==2.2) * dp22t_dbeta231/p22t + 
                                             1*(transit==2.3) * dp23t_dbeta231/p23t))
      }
      sim.data.wsf$score_k <- score_k
      sim.data.wsf$score_beta120 <- score_beta120
      sim.data.wsf$score_beta230 <- score_beta230
      sim.data.wsf$score_beta121 <- score_beta121
      sim.data.wsf$score_beta231 <- score_beta231
      
      score_k_cluster <- aggregate(x = sim.data.wsf$score_k,
                                   by = list(sim.data.wsf$cluster_id), 
                                   FUN = sum)[,2]
      
      score_beta120_cluster <- aggregate(x = sim.data.wsf$score_beta120,  # Specify data column
                                         by = list(sim.data.wsf$cluster_id), # Specify group indicator
                                         FUN = sum)[,2]
      
      score_beta230_cluster <- aggregate(x = sim.data.wsf$score_beta230,  
                                         by = list(sim.data.wsf$cluster_id), 
                                         FUN = sum)[,2]
      
      score_beta121_cluster <- aggregate(x = sim.data.wsf$score_beta121,
                                         by = list(sim.data.wsf$cluster_id), 
                                         FUN = sum)[,2]
      
      score_beta231_cluster <- aggregate(x = sim.data.wsf$score_beta231,
                                         by = list(sim.data.wsf$cluster_id), 
                                         FUN = sum)[,2]
      
      
      U1_U1 <- sum(score_k_cluster*score_k_cluster)
      U1_U2 <- sum(score_k_cluster*score_beta120_cluster)
      U1_U3 <- sum(score_k_cluster*score_beta230_cluster)
      U1_U4 <- sum(score_k_cluster*score_beta121_cluster)
      U1_U5 <- sum(score_k_cluster*score_beta231_cluster)
      U2_U2 <- sum(score_beta120_cluster*score_beta120_cluster)
      U2_U3 <- sum(score_beta120_cluster*score_beta230_cluster)
      U2_U4 <- sum(score_beta120_cluster*score_beta121_cluster)
      U2_U5 <- sum(score_beta120_cluster*score_beta231_cluster)
      U3_U3 <- sum(score_beta230_cluster*score_beta230_cluster)
      U3_U4 <- sum(score_beta230_cluster*score_beta121_cluster)
      U3_U5 <- sum(score_beta230_cluster*score_beta231_cluster)
      U4_U4 <- sum(score_beta121_cluster*score_beta121_cluster)
      U4_U5 <- sum(score_beta121_cluster*score_beta231_cluster)
      U5_U5 <- sum(score_beta231_cluster*score_beta231_cluster)
      
      V <- matrix(c(U1_U1, U1_U2, U1_U3, U1_U4, U1_U5,
                    U1_U2, U2_U2, U2_U3, U2_U4, U2_U5,
                    U1_U3, U2_U3, U3_U3, U3_U4, U3_U5,
                    U1_U4, U2_U4, U3_U4, U4_U4, U4_U5,
                    U1_U5, U2_U5, U3_U5, U4_U5, U5_U5), nrow = 5)
      
      weibull.gee.cov <- solve(0.5*H) %*% V %*% solve(0.5*H)
      gee.weibull.cor <- cor(weibull.gee.cov)
      weibull.gee.se <- sqrt(diag(weibull.gee.cov))
      
      #### wcr weibull ####
      
      beta.mat <- c()
      covmat.sum <- matrix(rep(0, 25), nrow = 5)
      
      logliklihood <- function(inital.para, data){
        k <- inital.para[1]
        beta120 <- inital.para[2]
        beta121 <- inital.para[4]
        beta230 <- inital.para[3]
        beta231 <- inital.para[5]
        q12 <- with(data, exp(beta120 * x0 + beta121 * (x2)))
        q23 <- with(data, exp(beta230 * x0 + beta231 * (x2)))
        p11 <- with(data, exp(-q12*(t2^(exp(k)) - t1^(exp(k)))))
        p22 <- with(data, exp(-q23*(t2^(exp(k)) - t1^(exp(k)))))
        p12 <- with(data, q12/(q12-q23) * (p22-p11))
        p23 <- with(data, 1- p22)
        p13 <- with(data, 1 - p11 - p12)
        p13[p13<=0] = 10^-10
        
        llk <- with(data, -2*sum((1*(transit==1.1)*log(p11)  + 
                                    1*(transit==1.2)*log(p12) + 
                                    1*(transit==1.3)*log(p13) +
                                    1*(transit==2.2)*log(p22) + 
                                    1*(transit==2.3)*log(p23))))
        return(llk)
      }
      
      
      
      inital.para <- c(log(1), as.vector(naive.est))
      
      for (b in 1:B){
        set.seed(simseed*10000+b)
        ##resampled dataset
        randomindividual <- sim.data.wsf %>% group_by(cluster_id) %>% distinct(cluster_individual) %>% sample_n(1)
        sim.resample <- sim.data.wsf %>% 
          filter(cluster_individual %in% randomindividual$cluster_individual) 
        
        weibull.optim.est <- optim(inital.para, logliklihood, data = sim.resample, method = "BFGS", hessian = T)
        beta.q <- weibull.optim.est$par
        H <- weibull.optim.est$hessian
        covmat.q <- solve(0.5*H)
        beta.mat <- cbind(beta.mat, beta.q)
        covmat.sum <- covmat.sum + covmat.q
      }
      
      beta.mean <- rowSums(beta.mat)/B
      beta.mat.dev <- beta.mat - matrix(rep(beta.mean, B), ncol = B)
      beta.var <-(covmat.sum/B - (beta.mat.dev %*% t(beta.mat.dev))/(B))
      weibull.wcr.est <- beta.mean
      wcr.cov <- beta.var
      wcr.weibull.cor <- cor(beta.var)
      weibull.wcr.se <- sqrt(diag(wcr.cov))
      
      
      #### Weibull WSF ####
      set.seed(simseed)
      inital.para <- c(log(1), as.vector(naive.est))
      logliklihood <- function(inital.para, data){
        gamma <- exp(inital.para[1])
        beta120 <- (inital.para[2])
        beta121 <- inital.para[4]
        beta230 <- (inital.para[3])
        beta231 <- inital.para[5]
        q12 <- with(data, exp(beta120+beta121*x2))
        q23 <- with(data, exp(beta230+beta231*x2))
        p11 <- with(data, exp(-q12*(t2^gamma - t1^gamma)))
        p22 <- with(data, exp(-q23*(t2^gamma - t1^gamma)))
        p12 <- with(data, q12/(q12-q23) * (p22-p11))
        p23 <- with(data, 1- p22)
        p13 <- with(data, 1 - p11 - p12)
        p13[p13<=0] = 10^-10
        llk <- with(data, -2*sum((1/cluster_size_vec)*(1*(transit==1.1)*log(p11)  + 
                                                         1*(transit==1.2)*log(p12) + 
                                                         1*(transit==1.3)*log(p13) +
                                                         1*(transit==2.2)*log(p22) + 
                                                         1*(transit==2.3)*log(p23))))
        return(llk)
      }
      
      weibull.optim.est <- optim(inital.para, logliklihood, data = sim.data.wsf, method = "BFGS", hessian = T, control = list(trace=1))
      weibull.wsf.est <-  weibull.optim.est$par
      H <- weibull.optim.est$hessian
      sqrt(diag(solve(H)))
      # covariance matrix 
      k <- weibull.optim.est$par[1]
      beta120 <- weibull.optim.est$par[2]
      beta230 <- weibull.optim.est$par[3]
      beta121 <- weibull.optim.est$par[4]
      beta231 <- weibull.optim.est$par[5]
      
      score_k <- c()
      score_beta120 <- c()
      score_beta230 <- c()
      score_beta121 <- c()
      score_beta231 <- c()
      for (i in 1:nrow(sim.data.wsf)){
        x0 <- sim.data.wsf$x0[i]
        x <- sim.data.wsf$x2[i]
        t1 <- sim.data.wsf$t1[i]
        t2 <- sim.data.wsf$t2[i]
        transit <- sim.data.wsf$transit[i]
        cluster.size <- sim.data.wsf$cluster_size_vec[i]
        
        q12 <- exp(beta120 * x0 + beta121 * x)
        q23 <- exp(beta230 * x0 + beta231 * x)
        y <- q12/(q12-q23)
        
        ##### transition probability #####
        p11t <- exp(-q12 * (t2^exp(k) - t1^exp(k)))
        p22t <- exp(-q23 * (t2^exp(k) - t1^exp(k)))
        p12t <- (q12/(q12-q23)) * (p22t - p11t)
        p13t <- 1-p11t-p12t
        p23t <- 1-p22t
        
        ##### first derivative of p #####
        dp11t_dk <- -p11t*q12*exp(k)*(t2^(exp(k))*log(t2) - t1^(exp(k))*log(t1))
        dp11t_dbeta120 <- -p11t*q12*(t2^exp(k) - t1^exp(k))*x0
        dp11t_dbeta230 <- 0
        dp11t_dbeta121 <- -p11t*q12*(t2^exp(k) - t1^exp(k))*x
        dp11t_dbeta231 <- 0
        
        dp22t_dk <- -p22t*q23*exp(k)*(t2^(exp(k))*log(t2) - t1^(exp(k))*log(t1))
        dp22t_dbeta120 <- 0
        dp22t_dbeta230 <- -p22t*q23*(t2^exp(k) - t1^exp(k))*x0
        dp22t_dbeta121 <- 0
        dp22t_dbeta231 <- -p22t*q23*(t2^exp(k) - t1^exp(k))*x
        
        dy_dk <- 0
        dy_dbeta120 <- -q12*q23*x0/(q12-q23)^2
        dy_dbeta230 <- q12*q23*x0/(q12-q23)^2
        dy_dbeta121 <- -q12*q23*x/(q12-q23)^2
        dy_dbeta231 <- q12*q23*x/(q12-q23)^2
        
        dp12t_dk <- y*(dp22t_dk - dp11t_dk)
        dp12t_dbeta120 <- dy_dbeta120*(p22t-p11t) + y*(dp22t_dbeta120-dp11t_dbeta120)
        dp12t_dbeta230 <- dy_dbeta230*(p22t-p11t) + y*(dp22t_dbeta230-dp11t_dbeta230)
        dp12t_dbeta121 <- dy_dbeta121*(p22t-p11t) + y*(dp22t_dbeta121-dp11t_dbeta121)
        dp12t_dbeta231 <- dy_dbeta231*(p22t-p11t) + y*(dp22t_dbeta231-dp11t_dbeta231)
        
        dp13t_dk <- -dp11t_dk - dp12t_dk
        dp13t_dbeta120 <- -dp11t_dbeta120 - dp12t_dbeta120
        dp13t_dbeta230 <- -dp11t_dbeta230 - dp12t_dbeta230
        dp13t_dbeta121 <- -dp11t_dbeta121 - dp12t_dbeta121
        dp13t_dbeta231 <- -dp11t_dbeta231 - dp12t_dbeta231
        
        dp23t_dk <- -dp22t_dk
        dp23t_dbeta120 <- -dp22t_dbeta120
        dp23t_dbeta230 <- -dp22t_dbeta230
        dp23t_dbeta121 <- -dp22t_dbeta121
        dp23t_dbeta231 <- -dp22t_dbeta231
        
        ##### 5*1 Vector of score function #####
        score_k <- c(score_k, (1*(transit==1.1) * dp11t_dk/p11t + 
                                 1*(transit==1.2) * dp12t_dk/p12t + 
                                 1*(transit==1.3) * dp13t_dk/p13t + 
                                 1*(transit==2.2) * dp22t_dk/p22t + 
                                 1*(transit==2.3) * dp23t_dk/p23t) / cluster.size)
        
        score_beta120 <- c(score_beta120, (1*(transit==1.1) * dp11t_dbeta120/p11t + 
                                             1*(transit==1.2) * dp12t_dbeta120/p12t + 
                                             1*(transit==1.3) * dp13t_dbeta120/p13t + 
                                             1*(transit==2.2) * dp22t_dbeta120/p22t + 
                                             1*(transit==2.3) * dp23t_dbeta120/p23t) / cluster.size)
        
        score_beta230 <- c(score_beta230, (1*(transit==1.1) * dp11t_dbeta230/p11t + 
                                             1*(transit==1.2) * dp12t_dbeta230/p12t + 
                                             1*(transit==1.3) * dp13t_dbeta230/p13t + 
                                             1*(transit==2.2) * dp22t_dbeta230/p22t + 
                                             1*(transit==2.3) * dp23t_dbeta230/p23t) / cluster.size)
        
        score_beta121 <- c(score_beta121, (1*(transit==1.1) * dp11t_dbeta121/p11t + 
                                             1*(transit==1.2) * dp12t_dbeta121/p12t + 
                                             1*(transit==1.3) * dp13t_dbeta121/p13t + 
                                             1*(transit==2.2) * dp22t_dbeta121/p22t + 
                                             1*(transit==2.3) * dp23t_dbeta121/p23t) / cluster.size)
        
        score_beta231 <- c(score_beta231, (1*(transit==1.1) * dp11t_dbeta231/p11t + 
                                             1*(transit==1.2) * dp12t_dbeta231/p12t + 
                                             1*(transit==1.3) * dp13t_dbeta231/p13t + 
                                             1*(transit==2.2) * dp22t_dbeta231/p22t + 
                                             1*(transit==2.3) * dp23t_dbeta231/p23t) / cluster.size)
      }
      sim.data.wsf$score_k <- score_k
      sim.data.wsf$score_beta120 <- score_beta120
      sim.data.wsf$score_beta230 <- score_beta230
      sim.data.wsf$score_beta121 <- score_beta121
      sim.data.wsf$score_beta231 <- score_beta231
      
      score_k_cluster <- aggregate(x = sim.data.wsf$score_k,
                                   by = list(sim.data.wsf$cluster_id), 
                                   FUN = sum)[,2]
      
      score_beta120_cluster <- aggregate(x = sim.data.wsf$score_beta120,  # Specify data column
                                         by = list(sim.data.wsf$cluster_id), # Specify group indicator
                                         FUN = sum)[,2]
      
      score_beta230_cluster <- aggregate(x = sim.data.wsf$score_beta230,  
                                         by = list(sim.data.wsf$cluster_id), 
                                         FUN = sum)[,2]
      
      score_beta121_cluster <- aggregate(x = sim.data.wsf$score_beta121,
                                         by = list(sim.data.wsf$cluster_id), 
                                         FUN = sum)[,2]
      
      score_beta231_cluster <- aggregate(x = sim.data.wsf$score_beta231,
                                         by = list(sim.data.wsf$cluster_id), 
                                         FUN = sum)[,2]
      
      
      U1_U1 <- sum(score_k_cluster*score_k_cluster)
      U1_U2 <- sum(score_k_cluster*score_beta120_cluster)
      U1_U3 <- sum(score_k_cluster*score_beta230_cluster)
      U1_U4 <- sum(score_k_cluster*score_beta121_cluster)
      U1_U5 <- sum(score_k_cluster*score_beta231_cluster)
      U2_U2 <- sum(score_beta120_cluster*score_beta120_cluster)
      U2_U3 <- sum(score_beta120_cluster*score_beta230_cluster)
      U2_U4 <- sum(score_beta120_cluster*score_beta121_cluster)
      U2_U5 <- sum(score_beta120_cluster*score_beta231_cluster)
      U3_U3 <- sum(score_beta230_cluster*score_beta230_cluster)
      U3_U4 <- sum(score_beta230_cluster*score_beta121_cluster)
      U3_U5 <- sum(score_beta230_cluster*score_beta231_cluster)
      U4_U4 <- sum(score_beta121_cluster*score_beta121_cluster)
      U4_U5 <- sum(score_beta121_cluster*score_beta231_cluster)
      U5_U5 <- sum(score_beta231_cluster*score_beta231_cluster)
      
      V <- matrix(c(U1_U1, U1_U2, U1_U3, U1_U4, U1_U5,
                    U1_U2, U2_U2, U2_U3, U2_U4, U2_U5,
                    U1_U3, U2_U3, U3_U3, U3_U4, U3_U5,
                    U1_U4, U2_U4, U3_U4, U4_U4, U4_U5,
                    U1_U5, U2_U5, U3_U5, U4_U5, U5_U5), nrow = 5)
      
      weibull.wsf.cov <- solve(0.5*H) %*% V %*% solve(0.5*H)
      wsf.weibull.cor <- cor(weibull.wsf.cov)
      weibull.wsf.se <- sqrt(diag(weibull.wsf.cov))
      #### Weibull WSF of data1 ####
      set.seed(simseed)
      inital.para <- c(log(1), as.vector(naive.est))
      logliklihood <- function(inital.para, data){
        gamma <- exp(inital.para[1])
        beta120 <- (inital.para[2])
        beta121 <- inital.para[4]
        beta230 <- (inital.para[3])
        beta231 <- inital.para[5]
        q12 <- with(data, exp(beta120+beta121*x2))
        q23 <- with(data, exp(beta230+beta231*x2))
        p11 <- with(data, exp(-q12*(t2^gamma - t1^gamma)))
        p22 <- with(data, exp(-q23*(t2^gamma - t1^gamma)))
        p12 <- with(data, q12/(q12-q23) * (p22-p11))
        p23 <- with(data, 1- p22)
        p13 <- with(data, 1 - p11 - p12)
        p13[p13<=0] = 10^-10
        llk <- with(data, -2*sum((1/cluster_size_vec)*(1*(transit==1.1)*log(p11)  + 
                                                         1*(transit==1.2)*log(p12) + 
                                                         1*(transit==1.3)*log(p13) +
                                                         1*(transit==2.2)*log(p22) + 
                                                         1*(transit==2.3)*log(p23))))
        return(llk)
      }
      
      weibull.optim.est1 <- optim(inital.para, logliklihood, data = sim.data.wsf1, method = "BFGS", hessian = T, control = list(trace=1))
      weibull.wsf.est1 <-  weibull.optim.est1$par
      H1 <- weibull.optim.est1$hessian
      
      # covariance matrix 
      k <- weibull.optim.est1$par[1]
      beta120 <- weibull.optim.est1$par[2]
      beta230 <- weibull.optim.est1$par[3]
      beta121 <- weibull.optim.est1$par[4]
      beta231 <- weibull.optim.est1$par[5]
      
      score_k <- c()
      score_beta120 <- c()
      score_beta230 <- c()
      score_beta121 <- c()
      score_beta231 <- c()
      for (i in 1:nrow(sim.data.wsf1)){
        x0 <- sim.data.wsf1$x0[i]
        x <- sim.data.wsf1$x2[i]
        t1 <- sim.data.wsf1$t1[i]
        t2 <- sim.data.wsf1$t2[i]
        transit <- sim.data.wsf1$transit[i]
        cluster.size <- sim.data.wsf1$cluster_size_vec[i]
        
        q12 <- exp(beta120 * x0 + beta121 * x)
        q23 <- exp(beta230 * x0 + beta231 * x)
        y <- q12/(q12-q23)
        
        ##### transition probability #####
        p11t <- exp(-q12 * (t2^exp(k) - t1^exp(k)))
        p22t <- exp(-q23 * (t2^exp(k) - t1^exp(k)))
        p12t <- (q12/(q12-q23)) * (p22t - p11t)
        p13t <- 1-p11t-p12t
        p23t <- 1-p22t
        
        ##### first derivative of p #####
        dp11t_dk <- -p11t*q12*exp(k)*(t2^(exp(k))*log(t2) - t1^(exp(k))*log(t1))
        dp11t_dbeta120 <- -p11t*q12*(t2^exp(k) - t1^exp(k))*x0
        dp11t_dbeta230 <- 0
        dp11t_dbeta121 <- -p11t*q12*(t2^exp(k) - t1^exp(k))*x
        dp11t_dbeta231 <- 0
        
        dp22t_dk <- -p22t*q23*exp(k)*(t2^(exp(k))*log(t2) - t1^(exp(k))*log(t1))
        dp22t_dbeta120 <- 0
        dp22t_dbeta230 <- -p22t*q23*(t2^exp(k) - t1^exp(k))*x0
        dp22t_dbeta121 <- 0
        dp22t_dbeta231 <- -p22t*q23*(t2^exp(k) - t1^exp(k))*x
        
        dy_dk <- 0
        dy_dbeta120 <- -q12*q23*x0/(q12-q23)^2
        dy_dbeta230 <- q12*q23*x0/(q12-q23)^2
        dy_dbeta121 <- -q12*q23*x/(q12-q23)^2
        dy_dbeta231 <- q12*q23*x/(q12-q23)^2
        
        dp12t_dk <- y*(dp22t_dk - dp11t_dk)
        dp12t_dbeta120 <- dy_dbeta120*(p22t-p11t) + y*(dp22t_dbeta120-dp11t_dbeta120)
        dp12t_dbeta230 <- dy_dbeta230*(p22t-p11t) + y*(dp22t_dbeta230-dp11t_dbeta230)
        dp12t_dbeta121 <- dy_dbeta121*(p22t-p11t) + y*(dp22t_dbeta121-dp11t_dbeta121)
        dp12t_dbeta231 <- dy_dbeta231*(p22t-p11t) + y*(dp22t_dbeta231-dp11t_dbeta231)
        
        dp13t_dk <- -dp11t_dk - dp12t_dk
        dp13t_dbeta120 <- -dp11t_dbeta120 - dp12t_dbeta120
        dp13t_dbeta230 <- -dp11t_dbeta230 - dp12t_dbeta230
        dp13t_dbeta121 <- -dp11t_dbeta121 - dp12t_dbeta121
        dp13t_dbeta231 <- -dp11t_dbeta231 - dp12t_dbeta231
        
        dp23t_dk <- -dp22t_dk
        dp23t_dbeta120 <- -dp22t_dbeta120
        dp23t_dbeta230 <- -dp22t_dbeta230
        dp23t_dbeta121 <- -dp22t_dbeta121
        dp23t_dbeta231 <- -dp22t_dbeta231
        
        ##### 5*1 Vector of score function #####
        score_k <- c(score_k, (1*(transit==1.1) * dp11t_dk/p11t + 
                                 1*(transit==1.2) * dp12t_dk/p12t + 
                                 1*(transit==1.3) * dp13t_dk/p13t + 
                                 1*(transit==2.2) * dp22t_dk/p22t + 
                                 1*(transit==2.3) * dp23t_dk/p23t) / cluster.size)
        
        score_beta120 <- c(score_beta120, (1*(transit==1.1) * dp11t_dbeta120/p11t + 
                                             1*(transit==1.2) * dp12t_dbeta120/p12t + 
                                             1*(transit==1.3) * dp13t_dbeta120/p13t + 
                                             1*(transit==2.2) * dp22t_dbeta120/p22t + 
                                             1*(transit==2.3) * dp23t_dbeta120/p23t) / cluster.size)
        
        score_beta230 <- c(score_beta230, (1*(transit==1.1) * dp11t_dbeta230/p11t + 
                                             1*(transit==1.2) * dp12t_dbeta230/p12t + 
                                             1*(transit==1.3) * dp13t_dbeta230/p13t + 
                                             1*(transit==2.2) * dp22t_dbeta230/p22t + 
                                             1*(transit==2.3) * dp23t_dbeta230/p23t) / cluster.size)
        
        score_beta121 <- c(score_beta121, (1*(transit==1.1) * dp11t_dbeta121/p11t + 
                                             1*(transit==1.2) * dp12t_dbeta121/p12t + 
                                             1*(transit==1.3) * dp13t_dbeta121/p13t + 
                                             1*(transit==2.2) * dp22t_dbeta121/p22t + 
                                             1*(transit==2.3) * dp23t_dbeta121/p23t) / cluster.size)
        
        score_beta231 <- c(score_beta231, (1*(transit==1.1) * dp11t_dbeta231/p11t + 
                                             1*(transit==1.2) * dp12t_dbeta231/p12t + 
                                             1*(transit==1.3) * dp13t_dbeta231/p13t + 
                                             1*(transit==2.2) * dp22t_dbeta231/p22t + 
                                             1*(transit==2.3) * dp23t_dbeta231/p23t) / cluster.size)
      }
      sim.data.wsf1$score_k <- score_k
      sim.data.wsf1$score_beta120 <- score_beta120
      sim.data.wsf1$score_beta230 <- score_beta230
      sim.data.wsf1$score_beta121 <- score_beta121
      sim.data.wsf1$score_beta231 <- score_beta231
      
      score_k_cluster <- aggregate(x = sim.data.wsf1$score_k,
                                   by = list(sim.data.wsf1$cluster_id), 
                                   FUN = sum)[,2]
      
      score_beta120_cluster <- aggregate(x = sim.data.wsf1$score_beta120,  # Specify data column
                                         by = list(sim.data.wsf1$cluster_id), # Specify group indicator
                                         FUN = sum)[,2]
      
      score_beta230_cluster <- aggregate(x = sim.data.wsf1$score_beta230,  
                                         by = list(sim.data.wsf1$cluster_id), 
                                         FUN = sum)[,2]
      
      score_beta121_cluster <- aggregate(x = sim.data.wsf1$score_beta121,
                                         by = list(sim.data.wsf1$cluster_id), 
                                         FUN = sum)[,2]
      
      score_beta231_cluster <- aggregate(x = sim.data.wsf1$score_beta231,
                                         by = list(sim.data.wsf1$cluster_id), 
                                         FUN = sum)[,2]
      
      
      U1_U1 <- sum(score_k_cluster*score_k_cluster)
      U1_U2 <- sum(score_k_cluster*score_beta120_cluster)
      U1_U3 <- sum(score_k_cluster*score_beta230_cluster)
      U1_U4 <- sum(score_k_cluster*score_beta121_cluster)
      U1_U5 <- sum(score_k_cluster*score_beta231_cluster)
      U2_U2 <- sum(score_beta120_cluster*score_beta120_cluster)
      U2_U3 <- sum(score_beta120_cluster*score_beta230_cluster)
      U2_U4 <- sum(score_beta120_cluster*score_beta121_cluster)
      U2_U5 <- sum(score_beta120_cluster*score_beta231_cluster)
      U3_U3 <- sum(score_beta230_cluster*score_beta230_cluster)
      U3_U4 <- sum(score_beta230_cluster*score_beta121_cluster)
      U3_U5 <- sum(score_beta230_cluster*score_beta231_cluster)
      U4_U4 <- sum(score_beta121_cluster*score_beta121_cluster)
      U4_U5 <- sum(score_beta121_cluster*score_beta231_cluster)
      U5_U5 <- sum(score_beta231_cluster*score_beta231_cluster)
      
      V1 <- matrix(c(U1_U1, U1_U2, U1_U3, U1_U4, U1_U5,
                     U1_U2, U2_U2, U2_U3, U2_U4, U2_U5,
                     U1_U3, U2_U3, U3_U3, U3_U4, U3_U5,
                     U1_U4, U2_U4, U3_U4, U4_U4, U4_U5,
                     U1_U5, U2_U5, U3_U5, U4_U5, U5_U5), nrow = 5)
      
      weibull.wsf.cov1 <- solve(0.5*H1) %*% V1 %*% solve(0.5*H1)
      wsf.weibull.cor1 <- cor(weibull.wsf.cov1)
      weibull.wsf.se1 <- sqrt(diag(weibull.wsf.cov1))
      
      
      
      
      #### Weibull WSF of data2 ####
      inital.para <- c(log(1), as.vector(naive.est))
      logliklihood <- function(inital.para, data){
        gamma <- exp(inital.para[1])
        beta120 <- (inital.para[2])
        beta121 <- inital.para[4]
        beta230 <- (inital.para[3])
        beta231 <- inital.para[5]
        q12 <- with(data, exp(beta120+beta121*x2))
        q23 <- with(data, exp(beta230+beta231*x2))
        p11 <- with(data, exp(-q12*(t2^gamma - t1^gamma)))
        p22 <- with(data, exp(-q23*(t2^gamma - t1^gamma)))
        p12 <- with(data, q12/(q12-q23) * (p22-p11))
        p23 <- with(data, 1- p22)
        p13 <- with(data, 1 - p11 - p12)
        p13[p13<=0] = 10^-10
        llk <- with(data, -2*sum((1/cluster_size_vec)*(1*(transit==1.1)*log(p11)  + 
                                                         1*(transit==1.2)*log(p12) + 
                                                         1*(transit==1.3)*log(p13) +
                                                         1*(transit==2.2)*log(p22) + 
                                                         1*(transit==2.3)*log(p23))))
        return(llk)
      }
      
      weibull.optim.est2 <- optim(inital.para, logliklihood, data = sim.data.wsf2, method = "BFGS", hessian = T, control = list(trace=1))
      weibull.wsf.est2 <-  weibull.optim.est2$par
      H2 <- weibull.optim.est2$hessian
      
      # covariance matrix 
      k <- weibull.optim.est2$par[1]
      beta120 <- weibull.optim.est1$par[2]
      beta230 <- weibull.optim.est1$par[3]
      beta121 <- weibull.optim.est1$par[4]
      beta231 <- weibull.optim.est1$par[5]
      
      score_k <- c()
      score_beta120 <- c()
      score_beta230 <- c()
      score_beta121 <- c()
      score_beta231 <- c()
      for (i in 1:nrow(sim.data.wsf2)){
        x0 <- sim.data.wsf2$x0[i]
        x <- sim.data.wsf2$x2[i]
        t1 <- sim.data.wsf2$t1[i]
        t2 <- sim.data.wsf2$t2[i]
        transit <- sim.data.wsf2$transit[i]
        cluster.size <- sim.data.wsf2$cluster_size_vec[i]
        
        q12 <- exp(beta120 * x0 + beta121 * x)
        q23 <- exp(beta230 * x0 + beta231 * x)
        y <- q12/(q12-q23)
        
        ##### transition probability #####
        p11t <- exp(-q12 * (t2^exp(k) - t1^exp(k)))
        p22t <- exp(-q23 * (t2^exp(k) - t1^exp(k)))
        p12t <- (q12/(q12-q23)) * (p22t - p11t)
        p13t <- 1-p11t-p12t
        p23t <- 1-p22t
        
        ##### first derivative of p #####
        dp11t_dk <- -p11t*q12*exp(k)*(t2^(exp(k))*log(t2) - t1^(exp(k))*log(t1))
        dp11t_dbeta120 <- -p11t*q12*(t2^exp(k) - t1^exp(k))*x0
        dp11t_dbeta230 <- 0
        dp11t_dbeta121 <- -p11t*q12*(t2^exp(k) - t1^exp(k))*x
        dp11t_dbeta231 <- 0
        
        dp22t_dk <- -p22t*q23*exp(k)*(t2^(exp(k))*log(t2) - t1^(exp(k))*log(t1))
        dp22t_dbeta120 <- 0
        dp22t_dbeta230 <- -p22t*q23*(t2^exp(k) - t1^exp(k))*x0
        dp22t_dbeta121 <- 0
        dp22t_dbeta231 <- -p22t*q23*(t2^exp(k) - t1^exp(k))*x
        
        dy_dk <- 0
        dy_dbeta120 <- -q12*q23*x0/(q12-q23)^2
        dy_dbeta230 <- q12*q23*x0/(q12-q23)^2
        dy_dbeta121 <- -q12*q23*x/(q12-q23)^2
        dy_dbeta231 <- q12*q23*x/(q12-q23)^2
        
        dp12t_dk <- y*(dp22t_dk - dp11t_dk)
        dp12t_dbeta120 <- dy_dbeta120*(p22t-p11t) + y*(dp22t_dbeta120-dp11t_dbeta120)
        dp12t_dbeta230 <- dy_dbeta230*(p22t-p11t) + y*(dp22t_dbeta230-dp11t_dbeta230)
        dp12t_dbeta121 <- dy_dbeta121*(p22t-p11t) + y*(dp22t_dbeta121-dp11t_dbeta121)
        dp12t_dbeta231 <- dy_dbeta231*(p22t-p11t) + y*(dp22t_dbeta231-dp11t_dbeta231)
        
        dp13t_dk <- -dp11t_dk - dp12t_dk
        dp13t_dbeta120 <- -dp11t_dbeta120 - dp12t_dbeta120
        dp13t_dbeta230 <- -dp11t_dbeta230 - dp12t_dbeta230
        dp13t_dbeta121 <- -dp11t_dbeta121 - dp12t_dbeta121
        dp13t_dbeta231 <- -dp11t_dbeta231 - dp12t_dbeta231
        
        dp23t_dk <- -dp22t_dk
        dp23t_dbeta120 <- -dp22t_dbeta120
        dp23t_dbeta230 <- -dp22t_dbeta230
        dp23t_dbeta121 <- -dp22t_dbeta121
        dp23t_dbeta231 <- -dp22t_dbeta231
        
        ##### 5*1 Vector of score function #####
        score_k <- c(score_k, (1*(transit==1.1) * dp11t_dk/p11t + 
                                 1*(transit==1.2) * dp12t_dk/p12t + 
                                 1*(transit==1.3) * dp13t_dk/p13t + 
                                 1*(transit==2.2) * dp22t_dk/p22t + 
                                 1*(transit==2.3) * dp23t_dk/p23t) / cluster.size)
        
        score_beta120 <- c(score_beta120, (1*(transit==1.1) * dp11t_dbeta120/p11t + 
                                             1*(transit==1.2) * dp12t_dbeta120/p12t + 
                                             1*(transit==1.3) * dp13t_dbeta120/p13t + 
                                             1*(transit==2.2) * dp22t_dbeta120/p22t + 
                                             1*(transit==2.3) * dp23t_dbeta120/p23t) / cluster.size)
        
        score_beta230 <- c(score_beta230, (1*(transit==1.1) * dp11t_dbeta230/p11t + 
                                             1*(transit==1.2) * dp12t_dbeta230/p12t + 
                                             1*(transit==1.3) * dp13t_dbeta230/p13t + 
                                             1*(transit==2.2) * dp22t_dbeta230/p22t + 
                                             1*(transit==2.3) * dp23t_dbeta230/p23t) / cluster.size)
        
        score_beta121 <- c(score_beta121, (1*(transit==1.1) * dp11t_dbeta121/p11t + 
                                             1*(transit==1.2) * dp12t_dbeta121/p12t + 
                                             1*(transit==1.3) * dp13t_dbeta121/p13t + 
                                             1*(transit==2.2) * dp22t_dbeta121/p22t + 
                                             1*(transit==2.3) * dp23t_dbeta121/p23t) / cluster.size)
        
        score_beta231 <- c(score_beta231, (1*(transit==1.1) * dp11t_dbeta231/p11t + 
                                             1*(transit==1.2) * dp12t_dbeta231/p12t + 
                                             1*(transit==1.3) * dp13t_dbeta231/p13t + 
                                             1*(transit==2.2) * dp22t_dbeta231/p22t + 
                                             1*(transit==2.3) * dp23t_dbeta231/p23t) / cluster.size)
      }
      sim.data.wsf2$score_k <- score_k
      sim.data.wsf2$score_beta120 <- score_beta120
      sim.data.wsf2$score_beta230 <- score_beta230
      sim.data.wsf2$score_beta121 <- score_beta121
      sim.data.wsf2$score_beta231 <- score_beta231
      
      score_k_cluster <- aggregate(x = sim.data.wsf2$score_k,
                                   by = list(sim.data.wsf2$cluster_id), 
                                   FUN = sum)[,2]
      
      score_beta120_cluster <- aggregate(x = sim.data.wsf2$score_beta120,  # Specify data column
                                         by = list(sim.data.wsf2$cluster_id), # Specify group indicator
                                         FUN = sum)[,2]
      
      score_beta230_cluster <- aggregate(x = sim.data.wsf2$score_beta230,  
                                         by = list(sim.data.wsf2$cluster_id), 
                                         FUN = sum)[,2]
      
      score_beta121_cluster <- aggregate(x = sim.data.wsf2$score_beta121,
                                         by = list(sim.data.wsf2$cluster_id), 
                                         FUN = sum)[,2]
      
      score_beta231_cluster <- aggregate(x = sim.data.wsf2$score_beta231,
                                         by = list(sim.data.wsf2$cluster_id), 
                                         FUN = sum)[,2]
      
      
      U1_U1 <- sum(score_k_cluster*score_k_cluster)
      U1_U2 <- sum(score_k_cluster*score_beta120_cluster)
      U1_U3 <- sum(score_k_cluster*score_beta230_cluster)
      U1_U4 <- sum(score_k_cluster*score_beta121_cluster)
      U1_U5 <- sum(score_k_cluster*score_beta231_cluster)
      U2_U2 <- sum(score_beta120_cluster*score_beta120_cluster)
      U2_U3 <- sum(score_beta120_cluster*score_beta230_cluster)
      U2_U4 <- sum(score_beta120_cluster*score_beta121_cluster)
      U2_U5 <- sum(score_beta120_cluster*score_beta231_cluster)
      U3_U3 <- sum(score_beta230_cluster*score_beta230_cluster)
      U3_U4 <- sum(score_beta230_cluster*score_beta121_cluster)
      U3_U5 <- sum(score_beta230_cluster*score_beta231_cluster)
      U4_U4 <- sum(score_beta121_cluster*score_beta121_cluster)
      U4_U5 <- sum(score_beta121_cluster*score_beta231_cluster)
      U5_U5 <- sum(score_beta231_cluster*score_beta231_cluster)
      
      V2 <- matrix(c(U1_U1, U1_U2, U1_U3, U1_U4, U1_U5,
                     U1_U2, U2_U2, U2_U3, U2_U4, U2_U5,
                     U1_U3, U2_U3, U3_U3, U3_U4, U3_U5,
                     U1_U4, U2_U4, U3_U4, U4_U4, U4_U5,
                     U1_U5, U2_U5, U3_U5, U4_U5, U5_U5), nrow = 5)
      
      weibull.wsf.cov2 <- solve(0.5*H2) %*% V2 %*% solve(0.5*H2)
      wsf.weibull.cor2 <- cor(weibull.wsf.cov2)
      weibull.wsf.se2 <- sqrt(diag(weibull.wsf.cov2))
      
    }, error=function(e) e)
    #### result ####
    if (inherits(out, "error")) {
      naive.exp.est.mat[s,] <- c(sim, rep(NA, 5))
      naive.exp.se.mat[s,] <- c(sim, rep(NA, 5))
      naive.weibull.est.mat[s,] <- c(sim, rep(NA, 5))
      naive.weibull.se.mat[s,] <- c(sim, rep(NA, 5))
      gee.weibull.est.mat[s,] <- c(sim, rep(NA, 5))
      gee.weibull.se.mat[s,] <- c(sim, rep(NA, 5))
      wcr.weibull.est.mat[s,] <- c(sim, rep(NA, 5))
      wcr.weibull.se.mat[s,] <- c(sim, rep(NA, 5))
      wsf.weibull.est.mat[s,] <- c(sim, rep(NA, 5))
      wsf.weibull.se.mat[s,] <- c(sim, rep(NA, 5))
      wsf.weibull.pw.est.mat[s,] <- c(sim, rep(NA, 10))
      wsf.weibull.pw.se.mat[s,] <- c(sim, rep(NA, 10))
      
      
    } else {
      naive.exp.est.mat[s,] <- c(sim, NA, as.vector(exp.naive.est))
      naive.exp.se.mat[s,] <- c(sim, NA, as.vector(exp.naive.se))
      naive.weibull.est.mat[s,] <- c(sim, as.vector(weibull.naive.est))
      naive.weibull.se.mat[s,] <- c(sim, as.vector(weibull.naive.se))
      gee.weibull.est.mat[s,] <- c(sim, as.vector(weibull.gee.est))
      gee.weibull.se.mat[s,] <- c(sim, as.vector(weibull.gee.se))
      wcr.weibull.est.mat[s,] <- c(sim, as.vector(weibull.wcr.est))
      wcr.weibull.se.mat[s,] <- c(sim, as.vector(weibull.wcr.se))
      wsf.weibull.est.mat[s,] <- c(sim, as.vector(weibull.wsf.est))
      wsf.weibull.se.mat[s,] <- c(sim, as.vector(weibull.wsf.se))
      wsf.weibull.pw.est.mat[s,] <- c(sim, as.vector(weibull.wsf.est1), as.vector(weibull.wsf.est2))
      wsf.weibull.pw.se.mat[s,] <- c(sim, weibull.wsf.se1, weibull.wsf.se2)
      
    }
  }
  return(list(naive.exp.est.mat, naive.exp.se.mat, naive.weibull.est.mat, naive.weibull.se.mat, 
              gee.weibull.est.mat, gee.weibull.se.mat,
              wcr.weibull.est.mat, wcr.weibull.se.mat, wsf.weibull.est.mat, wsf.weibull.se.mat,
              wsf.weibull.pw.est.mat, wsf.weibull.pw.se.mat))
}



system.time(simulation.out <- mclapply(1:ncores, simulation.func, mc.cores = ncores, mc.silent = F))
simout.df <- as.data.frame(do.call(rbind, simulation.out))

#### organize the simulation resutls ####

naive.exp.est.result <- as.data.frame(do.call(rbind, simout.df$V1))
naive.exp.se.result <- as.data.frame(do.call(rbind, simout.df$V2))
naive.weibull.est.result <- as.data.frame(do.call(rbind, simout.df$V3))
naive.weibull.se.result <- as.data.frame(do.call(rbind, simout.df$V4))
wcr.weibull.est.result <- as.data.frame(do.call(rbind, simout.df$V5))
wcr.weibull.se.result <- as.data.frame(do.call(rbind, simout.df$V6))
wsf.weibull.est.result <- as.data.frame(do.call(rbind, simout.df$V7))
wsf.weibull.se.result <- as.data.frame(do.call(rbind, simout.df$V8))
wsf.weibull.pw.est.result <- as.data.frame(do.call(rbind, simout.df$V9))
wsf.weibull.pw.se.result <- as.data.frame(do.call(rbind, simout.df$V10))

wcr.weibull.est.result <- wcr.weibull.est.result[complete.cases(wcr.weibull.se.result),]
wcr.weibull.se.result <- wcr.weibull.se.result[complete.cases(wcr.weibull.se.result),]
dim(wcr.weibull.est.result)
dim(wcr.weibull.se.result)

colnames(naive.exp.est.result) <- c("sim", "log.gamma", "log.lambda120", "log.lambda230", "beta12", "beta23")
colnames(naive.exp.se.result) <- c("sim", "log.gamma", "log.lambda120", "log.lambda230", "beta12", "beta23")
colnames(naive.weibull.est.result) <- c("sim", "log.gamma", "log.lambda120", "log.lambda230", "beta12", "beta23")
colnames(naive.weibull.se.result) <- c("sim", "log.gamma", "log.lambda120", "log.lambda230", "beta12", "beta23")
colnames(wcr.weibull.est.result) <- c("sim", "log.gamma", "log.lambda120", "log.lambda230", "beta12", "beta23")
colnames(wcr.weibull.se.result) <- c("sim", "log.gamma", "log.lambda120", "log.lambda230", "beta12", "beta23")
colnames(wsf.weibull.est.result) <- c("sim", "log.gamma", "log.lambda120", "log.lambda230", "beta12", "beta23")
colnames(wsf.weibull.se.result) <- c("sim", "log.gamma", "log.lambda120", "log.lambda230", "beta12", "beta23")
colnames(wsf.weibull.pw.est.result) <- c("sim", "log.gamma.t1", "log.lambda120.t1", "log.lambda230.t1", "beta12.t1", "beta23.t1", "log.gamma.t2", "log.lambda120.t2", "log.lambda230.t2", "beta12.t2", "beta23.t2")
colnames(wsf.weibull.pw.se.result) <- c("sim", "log.gamma.t1", "log.lambda120.t1", "log.lambda230.t1", "beta12.t1", "beta23.t1", "log.gamma.t2", "log.lambda120.t2", "log.lambda230.t2", "beta12.t2", "beta23.t2")

if (task_id==1){
  write.csv(naive.exp.est.result, paste0("naive.exp.est_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"), row.names =F, col.names = T)
  write.csv(naive.exp.se.result, paste0("naive.exp.se_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"), row.names =F, col.names = T)
  write.csv(naive.weibull.est.result, paste0("naive.weibull.est_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"), row.names =F, col.names = T)
  write.csv(naive.weibull.se.result, paste0("naive.weibull.se_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"), row.names =F, col.names = T)
  write.csv(wcr.weibull.est.result, paste0("wcr.weibull.est_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"), row.names =F, col.names = T)
  write.csv(wcr.weibull.se.result, paste0("wcr.weibull.se_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"), row.names =F, col.names = T)
  write.csv(wsf.weibull.est.result, paste0("wsf.weibull.est_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"), row.names =F, col.names = T)
  write.csv(wsf.weibull.se.result, paste0("wsf.weibull.se_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"), row.names =F, col.names = T)
  write.csv(wsf.weibull.pw.est.result, paste0("wsf.weibull.pw.est_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"), row.names =F, col.names = T)
  write.csv(wsf.weibull.pw.se.result, paste0("wsf.weibull.pw.se_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"), row.names =F, col.names = T)
}else {
  write.csv(naive.exp.est.result, paste0("naive.exp.est_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"), row.names =F, col.names = F)
  write.csv(naive.exp.se.result, paste0("naive.exp.se_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"), row.names =F, col.names = F)
  write.csv(naive.weibull.est.result, paste0("naive.weibull.est_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"), row.names =F, col.names = F)
  write.csv(naive.weibull.se.result, paste0("naive.weibull.se_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"), row.names =F, col.names = F)
  write.csv(wcr.weibull.est.result, paste0("wcr.weibull.est_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"), row.names =F, col.names = F)
  write.csv(wcr.weibull.se.result, paste0("wcr.weibull.se_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"), row.names =F, col.names = F)
  write.csv(wsf.weibull.est.result, paste0("wsf.weibull.est_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"), row.names =F, col.names = F)
  write.csv(wsf.weibull.se.result, paste0("wsf.weibull.se_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"), row.names =F, col.names = F)
  write.csv(wsf.weibull.pw.est.result, paste0("wsf.weibull.pw.est_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"), row.names =F, col.names = F)
  write.csv(wsf.weibull.pw.se.result, paste0("wsf.weibull.pw.se_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"), row.names =F, col.names = F)
}


naive.exp.est.mean <- unlist(lapply(naive.exp.est.result, mean))
naive.exp.est.sd <- unlist(lapply(naive.exp.est.result, sd))
naive.exp_est <- c( naive.exp.est.mean[2:6])
naive.exp_sd <- c( naive.exp.est.sd[2:6])
naive.exp_mse <- c(unlist(lapply(naive.exp.se.result, mean))[2:6])

naive.weibull.est.mean <- unlist(lapply(naive.weibull.est.result, mean))
naive.weibull.est.sd <- unlist(lapply(naive.weibull.est.result, sd))
naive.weibull_est <- c( naive.weibull.est.mean[2:6])
naive.weibull_sd <- c( naive.weibull.est.sd[2:6])
naive.weibull_mse <- c(unlist(lapply(naive.weibull.se.result, mean))[2:6])

wcr.weibull.est.mean <- unlist(lapply(wcr.weibull.est.result, mean))
wcr.weibull.est.sd <- unlist(lapply(wcr.weibull.est.result, sd))
wcr.weibull_est <- c( wcr.weibull.est.mean[2:6])
wcr.weibull_sd <- c( wcr.weibull.est.sd[2:6])
wcr.weibull_mse <- c(unlist(lapply(wcr.weibull.se.result, mean))[2:6])

wsf.weibull.est.mean <- unlist(lapply(wsf.weibull.est.result, mean))
wsf.weibull.est.sd <- unlist(lapply(wsf.weibull.est.result, sd))
wsf.weibull_est <- c( wsf.weibull.est.mean[2:6])
wsf.weibull_sd <- c( wsf.weibull.est.sd[2:6])
wsf.weibull_mse <- c(unlist(lapply(wsf.weibull.se.result, mean))[2:6])

wsf.weibull.pw.est.mean <- unlist(lapply(wsf.weibull.pw.est.result, mean))
wsf.weibull.pw.est.sd <- unlist(lapply(wsf.weibull.pw.est.result, sd))
wsf.weibull.pw_est <- c( wsf.weibull.pw.est.mean[2:11])
wsf.weibull.pw_sd <- c( wsf.weibull.pw.est.sd[2:11])
wsf.weibull.pw_mse <- c(unlist(lapply(wsf.weibull.pw.se.result, mean))[2:11])


parameter.true <- c(log(gamma), log(baseline_haz12), log(baseline_haz23),  beta2_12 , beta2_23)


sim_final_result <- round(rbind(parameter.true,
                                naive.exp_est, naive.exp_sd, naive.exp_mse,
                                naive.weibull_est, naive.weibull_sd, naive.weibull_mse,
                                wcr.weibull_est, wcr.weibull_sd, wcr.weibull_mse,
                                wsf.weibull_est, wsf.weibull_sd, wsf.weibull_mse), 4)
colnames(sim_final_result) <- c("log.gamma", "log.lambda120", "log.lambda230", "beta12", "beta23")

sim_final_result_wsf.pw <- round(rbind(parameter.true,wsf.weibull.pw_est, wsf.weibull.pw_sd, wsf.weibull.pw_mse ),4)
colnames(sim_final_result_wsf.pw) <- c("log.gamma.t1", "log.lambda120.t1", "log.lambda230.t1", "beta12.t1", "beta23.t1", "log.gamma.t2", "log.lambda120.t2", "log.lambda230.t2", "beta12.t2", "beta23.t2")


write.csv(sim_final_result, paste0("marginal.result_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, ".csv"),  row.names =T)
write.csv(sim_final_result_wsf.pw, paste0("wsf.pw.result_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, ".csv"),  row.names =T)

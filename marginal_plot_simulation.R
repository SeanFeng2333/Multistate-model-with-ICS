#### load packages ####
library(ggplot2)
library(grid)
library(gridExtra)
library(lemon)

#### Ordinal logits m500 tau 1 #### 
m = 500
tau=1
baseline_haz12 <- 0.3
baseline_haz23 <- 0.3
beta2_12 <- 0.5
beta2_23 <- 0.5
visittimes <- 5
gamma=1
task_id=1


wsf.weibull.est <- read.csv(paste0("wsf.weibull.est_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"))[,2:6]
wsf.weibull.se <- read.csv(paste0("wsf.weibull.se_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"))[,2:6]
sim_final_result <- read.csv(paste0("marginal.result_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, ".csv"))


wsf.weibull.pw.est <- read.csv(paste0("wsf.weibull.pw.est_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"))[2:11]
wsf.weibull.pw.se <- read.csv(paste0("wsf.weibull.pw.se_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"))[2:11]
sim_final_result_wsf.pw <- read.csv(paste0("wsf.pw.result_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, ".csv"))


HR12.wsf.weibull.lci <- mean(exp(wsf.weibull.est$beta12+qnorm(0.025)*wsf.weibull.se$beta12))
HR12.wsf.weibull.uci <- mean(exp(wsf.weibull.est$beta12+qnorm(0.975)*wsf.weibull.se$beta12))

HR23.wsf.weibull.lci <- mean(exp(wsf.weibull.est$beta23+qnorm(0.025)*wsf.weibull.se$beta23))
HR23.wsf.weibull.uci <- mean(exp(wsf.weibull.est$beta23+qnorm(0.975)*wsf.weibull.se$beta23))


wsf.weibull.lambda120.lci <- mean(exp(wsf.weibull.est$log.lambda120) + 
                                    qnorm(0.025)* exp(wsf.weibull.est$log.lambda120)*(wsf.weibull.se$log.lambda120))
wsf.weibull.lambda120.uci <- mean(exp(wsf.weibull.est$log.lambda120) + 
                                    qnorm(0.975)* exp(wsf.weibull.est$log.lambda120)*(wsf.weibull.se$log.lambda120))
wsf.weibull.lambda230.lci <- mean(exp(wsf.weibull.est$log.lambda230) + 
                                    qnorm(0.025)* exp(wsf.weibull.est$log.lambda230)*(wsf.weibull.se$log.lambda230))
wsf.weibull.lambda230.uci <- mean(exp(wsf.weibull.est$log.lambda230) + 
                                    qnorm(0.975)* exp(wsf.weibull.est$log.lambda230)*(wsf.weibull.se$log.lambda230))
wsf.weibull.gamma.lci <- mean(exp(wsf.weibull.est$log.gamma) + 
                                qnorm(0.025)* exp(wsf.weibull.est$log.gamma)*(wsf.weibull.se$log.gamma))
wsf.weibull.gamma.uci <- mean(exp(wsf.weibull.est$log.gamma) + 
                                qnorm(0.975)* exp(wsf.weibull.est$log.gamma)*(wsf.weibull.se$log.gamma))


#wsf.weibull.beta23.se <- exp(wsf.weibull.est$beta23)*wsf.weibull.se$beta23
#HR23.wsf.weibull.lci <- mean(exp(wsf.weibull.est$beta23)+qnorm(0.025)*wsf.weibull.beta23.se)
#HR23.wsf.weibull.uci <- mean(exp(wsf.weibull.est$beta23)+qnorm(0.975)*wsf.weibull.beta23.se)


wsf.weibull.pw.beta12.t1.lci <- wsf.weibull.pw.est$beta12.t1+qnorm(0.025)*wsf.weibull.pw.se$beta12.t1
wsf.weibull.pw.beta12.t1.uci <- wsf.weibull.pw.est$beta12.t1+qnorm(0.975)*wsf.weibull.pw.se$beta12.t1
HR12.wsf.weibull.pw.t1.lci <- mean(exp(wsf.weibull.pw.beta12.t1.lci))
HR12.wsf.weibull.pw.t1.uci <- mean(exp(wsf.weibull.pw.beta12.t1.uci))

wsf.weibull.pw.beta12.t2.lci <- wsf.weibull.pw.est$beta12.t2+qnorm(0.025)*wsf.weibull.pw.se$beta12.t2
wsf.weibull.pw.beta12.t2.uci <- wsf.weibull.pw.est$beta12.t2+qnorm(0.975)*wsf.weibull.pw.se$beta12.t2
HR12.wsf.weibull.pw.t2.lci <- mean(exp(wsf.weibull.pw.beta12.t2.lci))
HR12.wsf.weibull.pw.t2.uci <- mean(exp(wsf.weibull.pw.beta12.t2.uci))

wsf.weibull.pw.beta23.t1.lci <- wsf.weibull.pw.est$beta23.t1+qnorm(0.025)*wsf.weibull.pw.se$beta23.t1
wsf.weibull.pw.beta23.t1.uci <- wsf.weibull.pw.est$beta23.t1+qnorm(0.975)*wsf.weibull.pw.se$beta23.t1
HR23.wsf.weibull.pw.t1.lci <- mean(exp(wsf.weibull.pw.beta23.t1.lci))
HR23.wsf.weibull.pw.t1.uci <- mean(exp(wsf.weibull.pw.beta23.t1.uci))

wsf.weibull.pw.beta23.t2.lci <- wsf.weibull.pw.est$beta23.t2+qnorm(0.025)*wsf.weibull.pw.se$beta23.t2
wsf.weibull.pw.beta23.t2.uci <- wsf.weibull.pw.est$beta23.t2+qnorm(0.975)*wsf.weibull.pw.se$beta23.t2
HR23.wsf.weibull.pw.t2.lci <- mean(exp(wsf.weibull.pw.beta23.t2.lci))
HR23.wsf.weibull.pw.t2.uci <- mean(exp(wsf.weibull.pw.beta23.t2.uci))


##### plot of marginal hazard ratios #####
t <- 0:10
t1 <- 0:4
t2 <- 4:10

HR12.true = function(t){
  cum.basehaz12.true <- baseline_haz12 * t^gamma
  (tau + cum.basehaz12.true)/(tau + cum.basehaz12.true*exp(beta2_12))*exp(beta2_12)
}
HR12.wsf.weibull <- exp(sim_final_result[8,5])*t^0
HR12.wsf.weibull.t1 <- exp(sim_final_result_wsf.pw[2,5])*t1^0
HR12.wsf.weibull.t2 <- exp( sim_final_result_wsf.pw[2,7])*t2^0

HR23.true = function(t){
  cum.basehaz23.true <- baseline_haz23 * t^gamma
  (tau + cum.basehaz23.true)/(tau + cum.basehaz23.true*exp(beta2_23))*exp(beta2_23)
}
HR23.wsf.weibull <- exp(sim_final_result[8,6])*t^0
HR23.wsf.weibull.t1 <- exp(sim_final_result_wsf.pw[2,6])*t1^0
HR23.wsf.weibull.t2 <- exp( sim_final_result_wsf.pw[2,8])*t2^0


#par(mar = c(4.1, 4.4, 4.1, 1.9), xaxs="i", yaxs="i")
HR12plot_tau1 <- ggplot(data.frame(t=c(0, 10)), aes(x=t)) + 
  stat_function(fun=HR12.true, aes(color="TRUE")) + 
  geom_line(aes(y = HR12.wsf.weibull[1], color = "WSF")) +
  geom_line(aes(x=c(0,4), y = HR12.wsf.weibull.t1[1], color = "2-piece WSF")) +
  geom_line(aes(x=c(4,10), y = HR12.wsf.weibull.t2[1], color = "2-piece WSF")) +
  geom_ribbon(aes(ymin= HR12.wsf.weibull.lci, ymax=HR12.wsf.weibull.uci), linetype=2, alpha=0.1, fill="red")+
  geom_ribbon(aes(x=c(0,4), ymin= HR12.wsf.weibull.pw.t1.lci, ymax=HR12.wsf.weibull.pw.t1.uci), linetype=2, alpha=0.1, fill="blue")+
  geom_ribbon(aes(x=c(4,10), ymin= HR12.wsf.weibull.pw.t2.lci, ymax=HR12.wsf.weibull.pw.t2.uci), linetype=2, alpha=0.1, fill="blue")+
  scale_colour_manual("", 
                      breaks = c("TRUE", "WSF", "2-piece WSF"),
                      values = c("black","red", "blue"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("solid", "solid", "solid")))) + 
  ylim(0.8,1.7)+
  labs(y = "Marginal hazard ratio", x="Years")+
  theme(legend.position="top", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 


HR23plot_tau1 <- ggplot(data.frame(t=c(0, 10)), aes(x=t)) + 
  stat_function(fun=HR23.true, aes(color="TRUE")) + 
  geom_line(aes(y = HR23.wsf.weibull[1], color = "WSF")) +
  geom_line(aes(x=c(0,4), y = HR23.wsf.weibull.t1[1], color = "2-piece WSF")) +
  geom_line(aes(x=c(4,10), y = HR23.wsf.weibull.t2[1], color = "2-piece WSF")) +
  geom_ribbon(aes(ymin= HR23.wsf.weibull.lci, ymax=HR23.wsf.weibull.uci), linetype=2, alpha=0.1, fill="red")+
  geom_ribbon(aes(x=c(0,4), ymin= HR23.wsf.weibull.pw.t1.lci, ymax=HR23.wsf.weibull.pw.t1.uci), linetype=2, alpha=0.1, fill="blue")+
  geom_ribbon(aes(x=c(4,10), ymin= HR23.wsf.weibull.pw.t2.lci, ymax=HR23.wsf.weibull.pw.t2.uci), linetype=2, alpha=0.1, fill="blue")+
  scale_colour_manual("", 
                      breaks = c("TRUE", "WSF", "2-piece WSF"),
                      values = c("black","red", "blue"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("solid", "solid", "solid")))) + 
  ylim(0.8,1.7)+
  #labs(y = "Hazard ratio", title = "Marginal hazard ratio of transition 2-to-3") + 
  labs(y = "Marginal hazard ratio", x="Years")+
  theme(legend.position="top", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 




##### Cumulative baseline hazard #####


t <- (0:100)*0.1
cumhaz12.true = function(t){
  (tau)*(log(tau+baseline_haz12*t) - log(tau))
}
cumhaz12.wsf.weibull = function(t){
  sim_final_result[8,8]*t^sim_final_result[8,7]
}

cumhaz12.wsf.weibull.pw = function(t){
  (sim_final_result_wsf.pw[2,10])*t^sim_final_result_wsf.pw[2,9]
}


cumhaz23.true = function(t){
  (1+tau)*(log(tau+baseline_haz12*t)-log(tau))
}
cumhaz23.wsf.weibull = function(t){
  sim_final_result[8,9]*t^sim_final_result[8,7]
}
cumhaz23.wsf.weibull.pw = function(t){
  (sim_final_result_wsf.pw[2,11])*t^sim_final_result_wsf.pw[2,9]
}


cumhaz12plot_tau1 <-  ggplot(data.frame(t=c(0, 10)), aes(x=t)) + 
  stat_function(fun=cumhaz12.true, aes(color="TRUE")) + 
  stat_function(fun=cumhaz12.wsf.weibull, aes(color="WSF")) + 
  stat_function(fun=cumhaz12.wsf.weibull.pw, aes(color="2-piece WSF")) + 
  scale_colour_manual("", 
                      breaks = c("TRUE", "WSF", "2-piece WSF"),
                      values = c("black","red", "blue"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("solid", "solid", "solid")))) + 
  ylim(0,1.5)+
  labs(y = "Marginal cumulative baseline hazard", x= "Years")+
  theme(legend.position="top", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 


cumhaz23plot_tau1 <- ggplot(data.frame(t=c(0, 10)), aes(x=t)) + 
  stat_function(fun=cumhaz23.true, aes(color="TRUE")) + 
  stat_function(fun=cumhaz23.wsf.weibull, aes(color="WSF")) + 
  stat_function(fun=cumhaz23.wsf.weibull.pw, aes(color="2-piece WSF")) + 
  scale_colour_manual("", 
                      breaks = c("TRUE", "WSF", "2-piece WSF"),
                      values = c("black","red", "blue"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("solid", "solid", "solid")))) + 
  ylim(0,3)+
  labs(y = "Marginal cumulative baseline hazard", x= "Years")+
  theme(legend.position="top", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 



grid_arrange_shared_legend(HR12plot_tau1, HR23plot_tau1, cumhaz12plot_tau1, cumhaz23plot_tau1, ncol=2,nrow = 2, 
                           position="right")





#### Ordinal logits m500 tau 10 ####
m = 500
tau=10
baseline_haz12 <- 0.3
baseline_haz23 <- 0.3
beta2_12 <- 0.5
beta2_23 <- 0.5
visittimes <- 5
gamma=1
task_id=1

wsf.weibull.est <- read.csv(paste0("wsf.weibull.est_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"))[,2:6]
wsf.weibull.se <- read.csv(paste0("wsf.weibull.se_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"))[,2:6]
sim_final_result <- read.csv(paste0("marginal.result_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, ".csv"))

wsf.weibull.pw.est <- read.csv(paste0("wsf.weibull.pw.est_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"))[,2:11]
wsf.weibull.pw.se <- read.csv(paste0("wsf.weibull.pw.se_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, "task_id", task_id, ".csv"))[,2:11]
sim_final_result_wsf.pw <- read.csv(paste0("wsf.pw.result_m", m, "ordinalLogits", "tau", tau, "gamma", gamma, "visit", visittimes, ".csv"))

HR12.wsf.weibull.lci <- mean(exp(wsf.weibull.est$beta12+qnorm(0.025)*wsf.weibull.se$beta12))
HR12.wsf.weibull.uci <- mean(exp(wsf.weibull.est$beta12+qnorm(0.975)*wsf.weibull.se$beta12))

HR23.wsf.weibull.lci <- mean(exp(wsf.weibull.est$beta23+qnorm(0.025)*wsf.weibull.se$beta23))
HR23.wsf.weibull.uci <- mean(exp(wsf.weibull.est$beta23+qnorm(0.975)*wsf.weibull.se$beta23))


wsf.weibull.lambda120.lci <- mean(exp(wsf.weibull.est$log.lambda120) + 
                                    qnorm(0.025)* exp(wsf.weibull.est$log.lambda120)*(wsf.weibull.se$log.lambda120))
wsf.weibull.lambda120.uci <- mean(exp(wsf.weibull.est$log.lambda120) + 
                                    qnorm(0.975)* exp(wsf.weibull.est$log.lambda120)*(wsf.weibull.se$log.lambda120))
wsf.weibull.lambda230.lci <- mean(exp(wsf.weibull.est$log.lambda230) + 
                                    qnorm(0.025)* exp(wsf.weibull.est$log.lambda230)*(wsf.weibull.se$log.lambda230))
wsf.weibull.lambda230.uci <- mean(exp(wsf.weibull.est$log.lambda230) + 
                                    qnorm(0.975)* exp(wsf.weibull.est$log.lambda230)*(wsf.weibull.se$log.lambda230))
wsf.weibull.gamma.lci <- mean(exp(wsf.weibull.est$log.gamma) + 
                                qnorm(0.025)* exp(wsf.weibull.est$log.gamma)*(wsf.weibull.se$log.gamma))
wsf.weibull.gamma.uci <- mean(exp(wsf.weibull.est$log.gamma) + 
                                qnorm(0.975)* exp(wsf.weibull.est$log.gamma)*(wsf.weibull.se$log.gamma))



wsf.weibull.pw.beta12.t1.lci <- wsf.weibull.pw.est$beta12.t1+qnorm(0.025)*wsf.weibull.pw.se$beta12.t1
wsf.weibull.pw.beta12.t1.uci <- wsf.weibull.pw.est$beta12.t1+qnorm(0.975)*wsf.weibull.pw.se$beta12.t1
HR12.wsf.weibull.pw.t1.lci <- mean(exp(wsf.weibull.pw.beta12.t1.lci))
HR12.wsf.weibull.pw.t1.uci <- mean(exp(wsf.weibull.pw.beta12.t1.uci))

wsf.weibull.pw.beta12.t2.lci <- wsf.weibull.pw.est$beta12.t2+qnorm(0.025)*wsf.weibull.pw.se$beta12.t2
wsf.weibull.pw.beta12.t2.uci <- wsf.weibull.pw.est$beta12.t2+qnorm(0.975)*wsf.weibull.pw.se$beta12.t2
HR12.wsf.weibull.pw.t2.lci <- mean(exp(wsf.weibull.pw.beta12.t2.lci))
HR12.wsf.weibull.pw.t2.uci <- mean(exp(wsf.weibull.pw.beta12.t2.uci))

wsf.weibull.pw.beta23.t1.lci <- wsf.weibull.pw.est$beta23.t1+qnorm(0.025)*wsf.weibull.pw.se$beta23.t1
wsf.weibull.pw.beta23.t1.uci <- wsf.weibull.pw.est$beta23.t1+qnorm(0.975)*wsf.weibull.pw.se$beta23.t1
HR23.wsf.weibull.pw.t1.lci <- mean(exp(wsf.weibull.pw.beta23.t1.lci))
HR23.wsf.weibull.pw.t1.uci <- mean(exp(wsf.weibull.pw.beta23.t1.uci))

wsf.weibull.pw.beta23.t2.lci <- wsf.weibull.pw.est$beta23.t2+qnorm(0.025)*wsf.weibull.pw.se$beta23.t2
wsf.weibull.pw.beta23.t2.uci <- wsf.weibull.pw.est$beta23.t2+qnorm(0.975)*wsf.weibull.pw.se$beta23.t2
HR23.wsf.weibull.pw.t2.lci <- mean(exp(wsf.weibull.pw.beta23.t2.lci))
HR23.wsf.weibull.pw.t2.uci <- mean(exp(wsf.weibull.pw.beta23.t2.uci))


##### plot of marginal hazard ratios #####

t <- 0:10
t1 <- 0:4
t2 <- 4:10
HR12.true_tau10 = function(t){
  cum.basehaz12.true <- baseline_haz12 * t^gamma
  (tau + cum.basehaz12.true)/(tau + cum.basehaz12.true*exp(beta2_12))*exp(beta2_12)
}
HR12.wsf.weibull_tau10 <- exp(sim_final_result[8,5])*t^0
HR12.wsf.weibull.t1_tau10 <- exp(sim_final_result_wsf.pw[2,5])*t1^0
HR12.wsf.weibull.t2_tau10 <- exp( sim_final_result_wsf.pw[2,7])*t2^0

HR23.true_tau10 = function(t){
  cum.basehaz23.true <- baseline_haz23 * t^gamma
  (tau + cum.basehaz23.true)/(tau + cum.basehaz23.true*exp(beta2_23))*exp(beta2_23)
}
HR23.wsf.weibull_tau10 <- exp(sim_final_result[8,6])*t^0
HR23.wsf.weibull.t1_tau10 <- exp(sim_final_result_wsf.pw[2,6])*t1^0
HR23.wsf.weibull.t2_tau10 <- exp( sim_final_result_wsf.pw[2,8])*t2^0


HR12plot_tau10 <- ggplot(data.frame(t=c(0, 10)), aes(x=t)) + 
  stat_function(fun=HR12.true_tau10, aes(color="TRUE")) + 
  geom_line(aes(y = HR12.wsf.weibull_tau10[1], color = "WSF")) +
  geom_line(aes(x=c(0,4), y = HR12.wsf.weibull.t1_tau10[1], color = "2-piece WSF")) +
  geom_line(aes(x=c(4,10), y = HR12.wsf.weibull.t2_tau10[1], color = "2-piece WSF")) +
  geom_ribbon(aes(ymin= HR12.wsf.weibull.lci, ymax=HR12.wsf.weibull.uci), linetype=2, alpha=0.1, fill="red")+
  geom_ribbon(aes(x=c(0,4), ymin= HR12.wsf.weibull.pw.t1.lci, ymax=HR12.wsf.weibull.pw.t1.uci), linetype=2, alpha=0.1, fill="blue")+
  geom_ribbon(aes(x=c(4,10), ymin= HR12.wsf.weibull.pw.t2.lci, ymax=HR12.wsf.weibull.pw.t2.uci), linetype=2, alpha=0.1, fill="blue")+
  scale_colour_manual("", 
                      breaks = c("TRUE", "WSF", "2-piece WSF"),
                      values = c("black","red", "blue"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("solid", "solid", "solid")))) + 
  ylim(1,1.8)+
  labs(y = "Marginal hazard ratio", x="Years", title = "Transition 1-to-2")+
  theme(legend.position="top", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 





HR23plot_tau10 <- ggplot(data.frame(t=c(0, 10)), aes(x=t)) + 
  stat_function(fun=HR23.true_tau10, aes(color="TRUE")) + 
  geom_line(aes(y = HR23.wsf.weibull_tau10[1], color = "WSF")) +
  geom_line(aes(x=c(0,4), y = HR23.wsf.weibull.t1_tau10[1], color = "2-piece WSF")) +
  geom_line(aes(x=c(4,10), y = HR23.wsf.weibull.t2_tau10[1], color = "2-piece WSF")) +
  geom_ribbon(aes(ymin= HR23.wsf.weibull.lci, ymax=HR23.wsf.weibull.uci), linetype=2, alpha=0.1, fill="red")+
  geom_ribbon(aes(x=c(0,4), ymin= HR23.wsf.weibull.pw.t1.lci, ymax=HR23.wsf.weibull.pw.t1.uci), linetype=2, alpha=0.1, fill="blue")+
  geom_ribbon(aes(x=c(4,10), ymin= HR23.wsf.weibull.pw.t2.lci, ymax=HR23.wsf.weibull.pw.t2.uci), linetype=2, alpha=0.1, fill="blue")+
  scale_colour_manual("", 
                      breaks = c("TRUE", "WSF", "2-piece WSF"),
                      values = c("black","red", "blue"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("solid", "solid", "solid")))) + 
  ylim(1,1.8)+
  labs(y = "Marginal hazard ratio", x="Years", title = "Transition 2-to-3")+
  theme(legend.position="top", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 




##### Cumulative baseline hazard #####


t <- (0:100)*0.1
cumhaz12.true_tau10 = function(t){
  (tau)*(log(tau+baseline_haz12*t) - log(tau))
}
cumhaz12.wsf.weibull_tau10 = function(t){
  sim_final_result[8,8]*t^sim_final_result[8,7]
}

cumhaz12.wsf.weibull.pw_tau10 = function(t){
  (sim_final_result_wsf.pw[2,10])*t^sim_final_result_wsf.pw[2,9]
}


cumhaz23.true_tau10 = function(t){
  (1+tau)*(log(tau+baseline_haz12*t)-log(tau))
}
cumhaz23.wsf.weibull_tau10 = function(t){
  sim_final_result[8,9]*t^sim_final_result[8,7]
}
cumhaz23.wsf.weibull.pw_tau10 = function(t){
  (sim_final_result_wsf.pw[2,11])*t^sim_final_result_wsf.pw[2,9]
}



cumhaz12plot_tau10 <- ggplot(data.frame(t=c(0, 10)), aes(x=t)) + 
  stat_function(fun=cumhaz12.true_tau10, aes(color="TRUE")) + 
  stat_function(fun=cumhaz12.wsf.weibull_tau10, aes(color="WSF")) + 
  stat_function(fun=cumhaz12.wsf.weibull.pw_tau10, aes(color="2-piece WSF")) + 
  scale_colour_manual("", 
                      breaks = c("TRUE", "WSF", "2-piece WSF"),
                      values = c("black","red", "blue"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("solid", "solid", "solid")))) + 
  ylim(0,3)+
  labs(y = "Marginal cumulative baseline hazard", x= "Years")+
  theme(legend.position="top", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 



cumhaz23plot_tau10 <- ggplot(data.frame(t=c(0, 10)), aes(x=t)) + 
  stat_function(fun=cumhaz23.true_tau10, aes(color="TRUE")) + 
  stat_function(fun=cumhaz23.wsf.weibull_tau10, aes(color="WSF")) + 
  stat_function(fun=cumhaz23.wsf.weibull.pw_tau10, aes(color="2-piece WSF")) + 
  scale_colour_manual("", 
                      breaks = c("TRUE", "WSF", "2-piece WSF"),
                      values = c("black","red", "blue"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("solid", "solid", "solid")))) + 
  ylim(0,3)+
  labs(y = "Marginal cumulative baseline hazard", x= "Years")+
  theme(legend.position="top", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 


grid_arrange_shared_legend(HR12plot_tau10, HR23plot_tau10, cumhaz12plot_tau10, cumhaz23plot_tau10, ncol=2,nrow = 2, 
                           position="right")



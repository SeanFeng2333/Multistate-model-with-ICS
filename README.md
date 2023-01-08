# Multistate model with ICS 

R codes for simulation study in "Marginal Clustered Multistate Models for Longitudinal Progressive Processes with Informative Cluster Size"


## marginal_weibull_ordinalLogits_simulation.R

R code for running simulation with ICS generated from the ordinal logistic regression model

The code includes: 

1. Code for data generation

2. Fitting simulated data with
   i. naive multistate model that assumes an exponential baseline hazard and independent subjects 
   ii. naive multistate model that assumes a Weibull baseline hazard and independent subjects 
   iii. GEE with independent correlaiton structures
   iv. Within cluster resampling method
   v. Weighted score function method
   vi. 2-piece weighted score function method

3. Organize and save the results



## marginal_weibull_quantileICS_simulation.R

R code for running simulation with ICS generated from the quantiles of the frailties

The code includes: 

1. Code for data generation

2. Fitting simulated data with
   i. naive multistate model that assumes an exponential baseline hazard and independent subjects 
   ii. naive multistate model that assumes a Weibull baseline hazard and independent subjects 
   iii. GEE with independent correlaiton structures
   iv. Within cluster resampling method
   v. Weighted score function method
   vi. 2-piece weighted score function method

3. Organize and save the results



## marginal_plot_simulation.R

R code for creating marginal hazard ratios and marginal cumulative baseline hazards for simulation results 

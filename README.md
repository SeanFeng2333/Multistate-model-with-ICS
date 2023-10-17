# Multistate model with ICS 

R codes for simulation study in "Marginal Clustered Multistate Models for Longitudinal Progressive Processes with Informative Cluster Size".

Note: the authors used computing clusters to run the computationally intensive simulation. 



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

3. Organize and save the results. The results were saved in the folder "Marginal.inference.ordinalLogits"



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

3. Organize and save the results. The results were saved in the folder "Marginal.inference.3quantileICS"



##  Marginal.inference.ordinalLogits

This folder cantains intermediate results from marginal_weibull_ordinalLogits_simulation.R. 

For each of the four simulation scenarios (m = 100, 500; tau = 1, 10), we ran the simulation 1000 times and saved the output.

We also summarized the simulation results including mean estiamtes, mean standard deviations and empirical standard errors for each sceanrio. 



##  Marginal.inference.3quantileICS

This folder cantains intermediate results from marginal_weibull_quantileICS_simulation.R. 

For each of the four simulation scenarios (m = 100, 500; tau = 1, 10), we ran the simulation 1000 times and saved the output.

We also summarized the simulation results including mean estiamtes, mean standard deviations and empirical standard errors for each sceanrio. 



## marginal_plot_simulation.R

R code for creating marginal hazard ratios and marginal cumulative baseline hazards for simulation results. 

Users can use the results from the folders "Marginal.inference.ordinalLogits" and "Marginal.inference.3quantileICS" to reporduce the plots in the article. 





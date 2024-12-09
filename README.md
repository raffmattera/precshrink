Here you can find the files to replicate the results of the paper "Improved precision matrix estimation for mean-variance portfolio
 selection". You can find two main folders, one for simulations and another one for the empirical analyses. The important codes
are the following:

A) Simulations code.R: this files implements the simulation study;
A.1) covestim.R: this files includes the functions used to estimate covariance matrix estimators in the paper. These
functions are used in the simulation study;

B) Applications code.R: this file implements the empirical anlysis, computes the Sharpe ration and the pvalues;
B.1) mvp.R: this code includes the function adopted for building portfolios with different estimators.
The results are equivalent of applying the estimators in covestim.R and then plug-in the estimates
in the MeanVariancePortfolio formula. In this codes we however use the results of Remark 1 in the paper;

C) Results evaluation code.R: this provides a faster way of proceeding the results obtained in the codes above.

Moreover, for the empirical analyses, users can also load the results of the out-of-sample portfolios to access
the final results faster. These are the files "results_m180.RDS", "results_m240.RDS" and "results_m360.RDS" for different
values of M.

The R version used is the 4.4.1. The packages, with their versions, are the following.

 PeerPerformance_2.2.5, lmtest_0.9-40, zoo_1.8-12 , sandwich_3.1-0,      
 data.table_1.15.4, doParallel_1.0.17, iterators_1.0.14, foreach_1.5.2,        
 MASS_7.3-60.2, nlshrink_1.0.1, RiskPortfolios_2.1.7,  rio_1.1.1         

* Files:
  - Files with names XXX_ngp.R are used when n > p.
  - Files with names XXX_nlp.R are used when n <= p.
  - If XXX = 'competitors', then it implements factor analysis using FANC, SPC, Rockova and George's method (without and with varimax rotation); see the main paper for citations.
  - If XXX = 'xfa', then it implements factor analysis using xfa and returns the values of Lambda and sigma^2 across the alpha-eta grid; see the main paper for the definition of alpha-eta grid.
  - If XXX = 'mdl_wts', then the results of xfa_YYY.R are imported from '../result/YYY/' directory, where YYY = {ngp, nlp}, and posterior weights of Lambda matrices are calculated across the alpha-eta grid using EBIC; see the main paper for the exact form of EBIC.
  - File simulation.R simulates the loadings matrices and data using the setup described in Section 5 of the main paper. 
  - File figures.R plots the results shown in Section 5 of the main paper.
  - File fit_xfa.R includes 'xfa' function that fits xfa across the alpha-eta grid. One can use this function to fit and select a factor model using xfa. 

* Use the files as follows:
  0. Skip steps 1-4 and plot the results using step 5. Steps 1-4 are time consuming, so we recommend using them on a cluster.
  1. Run simulation.R file to store all the data sets in data directory.
  2. Choose the combination of (mtd, id) that corresponds to the right simulation combination of (n, p, replication).
  3. Run XXX_YYY.R file as 

     R CMD BATCH --no-save --no-restore "--args mtd id" XXX_YYY.R XXX_YYY_mtd_id.rout, 

     where XXX = {competitors, xfa}, YYY = {ngp, nlp}, and (mtd, id) are replaced by the choices in step 2. This will store the results for the competitors and the partial result for xfa in the 'results/YYY' directory.
  4. If XXX = 'xfa', then run mld_wts_YYY.R file as 
     
     R CMD BATCH --no-save --no-restore "--args mtd id" mdl_wts_YYY.R mdl_wts_YYY_mtd_id.rout, 

     where (mtd, id) are replaced by the choices in step 2. This will store the final result for xfa in the 'results/YYY' directory.
  5. Run

     R CMD BATCH --no-save --no-restore figures.R figs.rout, 

     to import the results, plot them, and save the files in 'result/img' directory.

* Comments:
  - Running all the replications is time consuming. We have summarized the results for rank (rnks), mean square error (mse), false discovery rate (fdr), and true positive rate (tpr) in files YYY_rnks.rds, YYY_mse.rds, YYY_fdr.rds, YYY_tpr.rds, where YYY = {ngp, nlp}.
  - File figures.R uses the files YYY_rnks.rds, YYY_mse.rds, YYY_fdr.rds, YYY_tpr.rds to plot the results.
  - The files XXX_ngp.R and XXX_nlp.R are coded so that they can be used on a cluster. For example, one can run all the replications on an SGE cluster for n = 5000 using a qsub file that includes the commands listed below:

  #!/bin/bash
  #$ -N xfa_ngp_5000
  #$ -t 1-50

  module load R/3.2.1
  
  R CMD BATCH --no-save --no-restore "--args 5 $SGE_TASK_ID" xfa_ngp.R xfa_ngp_5k_$SGE_TASK_ID.rout
 
* Contact:
  Please email Sanvesh Srivastava <sanvesh-srivastava@uiowa.edu> if you have any questions related to the code.

   
  
  

library(BSFG)  # devtools::install_github('deruncie/SparseFactorMixedModel',subdir='BSFG',ref='re-write-code')
library(GridLMM) #devtools::install_github('deruncie/GridLMM', ref = "develop")

#load my data 
rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS")
myPhenos <- read.table("data/GEMMA_eachAt_Bc/02_GEMMA/binMAF20NA10.fam")
mySNP <- read.table("data/GEMMA_eachAt_Bc/01_PLINK/dpbinMAF20NA10.ped")
myKin <- read.table("data/GEMMA_eachAt_Bc/03_kmat/binMAF20NA10_kmat1_pheno1.cXX.txt")
myRun <- 1.1
myPath <- "data/Runcie_BSFG"

Y = myPhenos # nxp gene expression matrix, already normalized and ready to go!
  # Should it be quantile-normalized?
  # According to Matthew Stephens, it's not necessary to normalize it, because a factor or two will become the normalization factor anyway
  data = mySNP # data.frame with n rows giving the sample information, in the same order as Y
  # In particular, a column "ID" of genotype ID's that match with the rownames of the Kinship matrix
  # also contains columns for any fixed effect factors that should be modeled.
  K_geno = myKin # n x n kinship matrix with rownames corresponding to genotype IDs
  # cis_genotypes = # list with p elements, each a n x b matrix with the local cis-genotypes for each of the p genes

  runID = myRun # unique identifier for this chain. Should identify the model and the chain. Probably want to run at least 2 chains to compare.

  results_folder = myPath #path to folder to store the results


  # initialize priors
  run_parameters = BSFG_control(
    num_NA_groups = 12,
    scale_Y = FALSE,
    simulation = FALSE,
    h2_divisions = 20,
    h2_step_size = .3,
    burn = 0000,
    thin = 10,
    K = 70,   # I'd start with 200 factors. The idea is to have more factors than are needed
    verbose=T
  )


priors = BSFG_priors(
  tot_Y_var = list(V = 0.5,   nu = 5),      # Prior variance of trait residuals after accounting for fixed effects and factors
  tot_F_var = list(V = 1, nu = 20),     # Prior variance of factor traits. This is included to improve MCMC mixing, but can be turned off by setting nu very large
  Lambda_prior = list(
    sampler = sample_Lambda_prec_horseshoe,
    prop_0 = 0.5,                      # Prior guess of number of genes on each factor
    delta_l = list(shape = 3, rate = 1), # prior on decrease in probability of inclusion of a gene on a factor. shape/rate is the expected decrease in odds of inclusion from factor j to j+1
    delta_iterations_factor = 100
  ),
  B2_prior = list(
    sampler = sample_B2_prec_horseshoe,
    prop_0 = 0.1
  ),
  cis_effects_prior = list(
    prec = 1e-6    # prior on cis-eQTL effect sizes. May need to be larger if there are multiple cis-genotypes per gene
  ),
  h2_priors_resids_fun = function(h2s,n) 1,  # Function that returns the prior density for any value of the h2s vector (ie the vector of random effect proportional variances across all random effects. 1 means constant prior. Alternative: pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
  h2_priors_factors_fun = function(h2s,n) 1 # See above. Another choice is one that gives 50% weight to h2==0: ifelse(h2s == 0,n,n/(n-1))
)

BSFG_state = setup_model_BSFG(Y,~ (1|ID),
                              data=data,
                              relmat = list(ID = K_geno),  # To run model #5, exclude this line
                              # cis_genotypes = cis_genotypes, # To run model #6, exclude this line
                              run_parameters=run_parameters,
                              run_ID = runID
)


BSFG_state = set_priors_BSFG(BSFG_state,priors)
BSFG_state = initialize_variables_BSFG(BSFG_state)
BSFG_state = initialize_BSFG(BSFG_state)
BSFG_state$Posterior$posteriorSample_params = c(BSFG_state$Posterior$posteriorSample_params,'Lambda_c2')
BSFG_state = clear_Posterior(BSFG_state)

n_samples = 100;  # how many samples to collect at once?
for(i  in 1:200) {
  print(sprintf('Run %d',i))
  BSFG_state = sample_BSFG(BSFG_state,n_samples,grainSize=1)  # run MCMC chain n_samples iterations. grainSize is a paramter for parallelization (smaller = more parallelization)

  BSFG_state = save_posterior_chunk(BSFG_state)  # save any accumulated posterior samples in the database to release memory
  print(BSFG_state) # print status of current chain
  plot(BSFG_state) # make some diagnostic plots. These are saved in a pdf booklet: diagnostic_plots.pdf

  # set of commands to run during burn-in period to help chain converge
  if(BSFG_state$current_state$nrun < BSFG_state$run_parameters$burn || i < 70) {
    # Factor order doesn't "mix" well in the MCMC. We can help it by manually re-ordering from biggest to smallest, and by removing factors that are too correlated
    BSFG_state = reorder_factors(BSFG_state,drop_cor_threshold = 0.6)
    BSFG_state = clear_Posterior(BSFG_state)
  }
  if(i < 100){
    BSFG_state = clear_Posterior(BSFG_state)
  }
}

# extract posterior samples of factors, and calculate posterior means
F_ind = get_posterior_mean(load_posterior_param(BSFG_state,'F'))
saveRDS(F_ind,file = sprintf('%s/%s_F_ind.rds',results_folder,BSFG_state$run_ID))  # the factor scores for GWAS
cis_effects = load_posterior_param(BSFG_state,'cis_effects')
saveRDS(cis_effects,file = sprintf('%s/%s_cis_effects.rds',results_folder,BSFG_state$run_ID))  # the posterior samples of the cis_effects (comment if not included in model)

Lambda_m_eff = load_posterior_param(BSFG_state,'Lambda_m_eff')
saveRDS(Lambda_m_eff,file = sprintf('%s/%s_Lambda_m_eff.rds',results_folder,BSFG_state$run_ID)) # the posterior samples of the effective number of genes/factor
Lambda_c2 = load_posterior_param(BSFG_state,'Lambda_c2')
saveRDS(Lambda_c2,file = sprintf('%s/%s_Lambda_c2.rds',results_folder,BSFG_state$run_ID))
Lambda = load_posterior_param(BSFG_state,'Lambda')
pm_Lambda = get_posterior_mean(Lambda)
saveRDS(pm_Lambda,file = sprintf('%s/%s_pm_Lambda.rds',results_folder,BSFG_state$run_ID))  # the posterior mean of the factor loadings

# system(sprintf('rm -rf %s',BSFG_state$run_ID)) # optional, to clean up at the end



# run GWAS with GridLMM. Alternatively, this could probably be scripted for GEMMA since I think there is only 1 random effects

X = # n x b matrix of genotypes. Needs rownames that correspond to data$ID

data$y = F_ind[,1]
m_null = GridLMM_ML(y~(1|ID),data,relmat = list(ID=K_geno),tolerance = .01,diagonalize = F)
V_setup = m_null$setup

p_value_matrix = matrix(0,ncol(X),ncol(F_ind))
for(i in 1:ncol(F_ind)) {
  print(sprintf('Column %d',i))
  data$y = F_ind[,i]
  res = GridLMM_GWAS(y~(1|ID),~1,~0,data = data,X = X,X_ID = 'ID',V_setup = V_setup,h2_step = 0.01,
                     centerX = FALSE, scaleX = FALSE,
                     verbose=F)
  p_value_matrix[,i] = res$results$p_value_REML
}
saveRDS(p_value_matrix,sprintf('%s/%s_p_value_matrix.rds',results_folder,run_ID))

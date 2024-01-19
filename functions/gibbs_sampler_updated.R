#' Load libraies
load_libraries = function(){
  shhh = function(lib_name){ # It's a library, so shhh!
    suppressWarnings(suppressMessages(require(lib_name, character.only = TRUE)))
  }
  
  shhh("tidyverse")
  shhh("ACutils")
  shhh("mvtnorm")
  shhh("salso")
  shhh("FGM")
  shhh("gmp")
  shhh("mcclust")
  shhh("mcclust.ext")
  shhh("logr")
  shhh("tidygraph")
  shhh("ggraph")
  shhh("igraph")
  shhh("Rcpp")
  shhh("RcppArmadillo")
  shhh("RcppEigen")
  
  source("utility_functions.R");
  source("bulky_functions.R");
  source("data_generation.R")
  sourceCpp("wade.cpp")
  Rcpp::sourceCpp('UpdateParamsGSL.cpp')
}


#' Gibbs sampler modified
#'
#' @param niter desired number of effective iterations.
#' @param seed seed for reproducibility.
#' @param initialization_value_list list containing values for initialization.
#' @param sigma integer corresponding to the value of sigma_prior_0 in set_options.
#' @param theta integer corresponding to the value of theta_prior_0 in set_options.
#' @param p integer corresponing to the number of basis.
#' @param graph_density graph density.
#' @param beta_sig2 variance of the Beta
#' @param BaseMat TODO
#' @param data TODO
#' @param a_tau_eps,
#' @param b_tau_eps,
#' @param sigma_mu,
#' @param r,
#' @param Update_Beta,
#' @param Update_Mu,
#' @param Update_Tau
#'
#' @return TODO
#' @export TODO
#'
#' @examples
new_gibbs_sampler = function(
    niter,
    seed,
    initialization_value_list,
    sigma, 
    theta, 
    p, 
    graph_density, 
    beta_sig2, 
    BaseMat, 
    data
    a_tau_eps,
    b_tau_eps,
    sigma_mu,
    r,
    Update_Beta,
    Update_Mu,
    Update_Tau
  ){
  
  # Load libraries
  load_libraries()
  
  # Set seed
  set.seed(seed)
  
  # Create a list to store the chains
  chains <- list(
    Beta = vector("list", length = niter),
    mu = vector("list", length = niter),
    tau_eps = vector("list", length = niter),
    K = vector("list", length = niter),
    G = vector("list", length = niter),
    z = vector("list", length = niter),
    rho = vector("list", length = niter),
    time = vector("list", length = niter)
  )
  
  # Initialization of the chains
  chains$Beta[[1]] <- initialization_value_list$Beta
  chains$mu[[1]] <- initialization_value_list$mu
  chains$tau_eps[[1]] <- initialization_value_list$tau_eps 
  chains$K[[1]] <- initialization_value_list$K
  chains$G[[1]] <- initialization_value_list$G
  chains$z[[1]] <- initialization_value_list$z
  chains$time <- initialization_value_list$time
  chains$rho[[1]] <- initialization_value_list$rho
  weights_a <- initialization_value_list$weights_a
  weights_d <- initialization_value_list$weights_d
  total_weights <- initialization_value_list$total_weights
  total_K = initialization_value_list$total_K
  total_graphs = initialization_value_list$total_graphs
  
  # Set the starting graph_start to NULL
  graph_start = NULL
  
  # TODO
  tbase_base = t(BaseMat)%*%BaseMat
  tbase_data = t(BaseMat)%*%t(data)
  Sdata = sum(diag(data%*%t(data)))
  
  # Start running the gibb sampler
  for(s in 2:niter){
    
    fit = UpdateParamsGSL(
      chains$Beta[[s-1]],
      chains$mu[[s-1]],
      chains$tau_eps[[s-1]],
      chains$K[[s-1]],
      tbase_base,
      tbase_data,
      Sdata,
      a_tau_eps,
      b_tau_eps,
      sigma_mu,
      r,
      Update_Beta,
      Update_Mu,
      Update_Tau
    )
    
    # Beta
    chains$Beta[[s]] <- fit$Beta
    
    # mu
    chains$mu[[s]] <- fit$mu
    
    # tau
    chains$tau_eps[[s]] <- fit$tau_eps
    
    
    options = set_options(
      sigma_prior_0 = sigma,
      sigma_prior_parameters = list("a"=1,"b"=1,"c"=1,"d"=1),
      theta_prior_0 = theta,
      theta_prior_parameters = list("c"=1,"d"=1),
      rho0 = chains$rho[[s-1]], # start with one group
      weights_a0 = weights_a,
      weights_d0 = weights_d,
      total_weights0 = total_weights,
      total_K0 = total_K,
      total_graphs0 = total_graphs,
      graph = graph_start,
      alpha_target = 0.234,
      beta_mu = graph_density,   # da modificare (expected value beta distr of the graph)
      beta_sig2 = beta_sig2,     # da modificare (var beta distr del grafo, fra 0 e 0.25)
      d = 3,                     # param della G wishart (default 3)
      alpha_add = 0.5,
      adaptation_step = 1/(p*1000),
      update_sigma_prior = TRUE,
      update_theta_prior = TRUE,
      update_weights = TRUE,
      update_partition = TRUE,
      update_graph = TRUE,
      perform_shuffle = TRUE
      )
    
    # Single iteration of the old Gibbs_sampler
    res <- Gibbs_sampler(
      data = t(fit$Beta - fit$mu),  # da provare
      niter = 1, # niter finali, già tolto il burn in
      nburn = 0,
      thin = 1,
      options = options,
      seed = 22111996,
      print = FALSE
    )
    
    # Update computed parameters
    z = do.call(rbind, lapply(res$rho, rho_to_z))
    
    chains$rho[[s]] <- res$rho[[1]]
    chains$K[[s]] <- res$K [[1]]
    chains$G[[s]] <- res$G [[1]]
    chains$z[[s]] <- z
    chains$time[[s]] <- res$execution_time
    
    weights_a <- res$weights_a
    weights_d <- res$weights_d
    
    total_weights <- res$total_weights
    total_K <- res$total_K [[1]]
    total_graphs <- res$total_graphs [[1]]
    
    graph_start = res$bdgraph_start
    
  }
  
  return(chains)
}
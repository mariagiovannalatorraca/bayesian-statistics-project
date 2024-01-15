## Load libraries
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

## Load custom functions
source("utility_functions.R");source("bulky_functions.R");source("data_generation.R")
sourceCpp("wade.cpp")

Rcpp::sourceCpp('UpdateParamsGSL.cpp')

library('RcppGSL')
library(fda)

# Load Data
# setwd("G:/.shortcut-targets-by-id/1Ck2MctcmCBWOueSMeW4Zr8IOvOaOzdqz/BicoccaDrive/TutoratoPolimi_2324/Materiale/FunctionalRandomGraph_shared/Codici")
load("purees.Rdat")
curves = purees$data
wavelengths = purees$wavelengths
strawberry <- curves[which(curves$Group == "Strawberry"), ]
data = strawberry[,-1]*100
data = as.matrix(data) # n x r

#Generate basis
p =  40
n = dim(data)[1]
r = dim(data)[2]
range_x = range(wavelengths)

data_W <- t(as.matrix(data))
basis <- create.bspline.basis(rangeval=range(wavelengths), nbasis=40, norder=3) 
data.fd <- Data2fd(y = data_W, argvals = wavelengths, basisobj = basis)
plot.fd(data.fd, main="B-splines")

BaseMat <- eval.basis(wavelengths,basis) # matrix r x p 
#plot(basis)



# Initialization 
seed = 22111996
set.seed(seed)
Beta = matrix(rnorm(n=p*n), nrow = p, ncol = n)   # matrice iniziale con errore
mu = rnorm(n=p)         # mu iniziale con errore
tau_eps = 100      # lo considera già al quadrato
K = rWishart(n=1,df = p+10, Sigma = diag(p) )
K = K[,,1]
# ACutils::ACheatmap(K,center_value = NULL, remove_diag = T)
tbase_base = t(BaseMat)%*%BaseMat # p x p (phi_t * phi)
tbase_data = t(BaseMat)%*%t(data) # p x n (phi_t * Y_t)  mettiamo insieme tutti i beta, verranno etrsatte le colonne nel file C
Sdata = sum(diag(data%*%t(data))) # inefficient calculation ((2b + Sdata) è il b di tau_eps a posteriori) qui mette phi*beta = 0 come inizializz.
Update_Beta <- TRUE
Update_Mu <- TRUE
Update_Tau <- TRUE
a_tau_eps <- b_tau_eps <- 2
sigma_mu = 100

# Define variance of the Beta
beta_sig2 = 0.2

# Compute graph density
graph_density = 0.3  # come si sceglie?

niter <- 10000     #almeno 10000

# Create a list for chains
chains <- list(Beta = vector("list", length = niter),
               mu = vector("list", length = niter),
               tau_eps = vector("list", length = niter),
               K = vector("list", length = niter),
               G = vector("list", length = niter),
               z = vector("list", length = niter),
               rho = vector("list", length = niter),
               time = vector("list", length = niter)
               ) 

# Initialization of the chains
chains$Beta[[1]] <- Beta
chains$mu[[1]] <- mu
chains$tau_eps[[1]] <- tau_eps 
chains$K[[1]] <- K
chains$G[[1]] <- diag(1,p,p)
chains$z[[1]] <- rep(1,p)

chains$time <- 0

chains$rho[[1]] <- p

sigma <-0.5
theta<-1

weights_a <- rep(1,p-1)
weights_d <- rep(1,p-1)
total_weights <- 0
total_K = matrix(0,p,p)
total_graphs = matrix(0,p,p)
graph_start = NULL

#dim(fit$U)
for(s in 2:niter) {

  fit = UpdateParamsGSL( chains$Beta[[s-1]], chains$mu[[s-1]], chains$tau_eps[[s-1]], 
                         chains$K[[s-1]], tbase_base, tbase_data, Sdata,   
                      a_tau_eps, b_tau_eps, sigma_mu, r,
                      Update_Beta, Update_Mu, Update_Tau )
  
  # Beta
  chains$Beta[[s]] <- fit$Beta
  
  # mu
  chains$mu[[s]] <- fit$mu
  
  # tau
  chains$tau_eps[[s]] <- fit$tau_eps
  
  
  options = set_options(  sigma_prior_0=sigma,
                          sigma_prior_parameters=list("a"=1,"b"=1,"c"=1,"d"=1),
                          theta_prior_0=theta,
                          theta_prior_parameters=list("c"=1,"d"=1),
                          rho0=chains$rho[[s-1]], # start with one group
                          weights_a0=weights_a,
                          weights_d0=weights_d,
                          total_weights0=total_weights,
                          total_K0 = total_K,
                          total_graphs0 = total_graphs,
                          graph = graph_start,
                          alpha_target=0.234,
                          beta_mu=graph_density,   # da modificare (expected value beta distr of the graph)
                          beta_sig2=beta_sig2,     # da modificare (var beta distr del grafo, fra 0 e 0.25)
                          d=3,                     # param della G wishart (default 3)
                          alpha_add=0.5,
                          adaptation_step=1/(p*1000),
                          update_sigma_prior=TRUE,
                          update_theta_prior=TRUE,
                          update_weights=TRUE,
                          update_partition=TRUE,
                          update_graph=TRUE,
                          perform_shuffle=TRUE  )
  
  
  ## Run 
  
  res <- Gibbs_sampler(
    data = t(fit$Beta - fit$mu),  # da provare
    niter = 1, # niter finali, già tolto il burn in
    nburn = 0,
    thin = 1,
    options = options,
    seed = 22111996,
    print = FALSE
  )
  
  
  z = do.call(rbind, lapply(res$rho, rho_to_z))
  
  chains$rho[[s]] <- res$rho[[1]]
  
  chains$K[[s]] <- res$K [[1]]
  
  chains$G[[s]] <- res$G [[1]]   #primo valore NULL o identità 
  
  chains$z[[s]] <- z      #primo valore NULL o mettiamo tutti 1
  
  #salvare qui dentro i tempi per ogni K
  chains$time[[s]] <- res$execution_time
  
  weights_a <- res$weights_a
  weights_d <- res$weights_d
  
  total_weights <- res$total_weights
  total_K <- res$total_K [[1]]
  total_graphs <- res$total_graphs [[1]]
  
  # update the initial graph
  graph_start = res$bdgraph_start
  
}




# # burn in 
burn_in <- 1000
# chain_beta <- chain_beta[(burn_in + 1):niter]
# chain_mu <- chain_mu[(burn_in + 1):niter]
# chain_tau <- chain_tau[(burn_in + 1):niter]



## Graph

# DA CAMBIARE: MEDIA CON I TEMPI OTTENUTI DA RES

# Extract last plinks
last_plinks = tail(res$G, n=1)[[1]]

# Criterion 1 to select the threshold (should not work very well) and assign final graph
threshold = 0.5
G_est <- matrix(0,p,p)
G_est[which(last_plinks>threshold)] = 1

#Criterion 2 to select the threshold
bfdr_select = BFDR_selection(last_plinks, tol = seq(0.1, 1, by = 0.001))

# Inspect the threshold and assign final graph
bfdr_select$best_treshold
G_est = bfdr_select$best_truncated_graph



ACutils::ACheatmap(
  tail(res$K,n=1)[[1]],  
  use_x11_device = F,
  horizontal = F,
  main = "Estimated Precision matrix",
  center_value = NULL,
  col.upper = "black",
  col.center = "grey50",
  col.lower = "white"
)


# K aggiornata
# somma dei tempi*K/tot_tempo
sum <- 0
for(i in (burn_in+1):niter){
  sum <- sum + chains$time[[i]]*chains$K[[i]]
}

K_fin <- sum/sum(chains$time)

ACutils::ACheatmap(
  K_fin,  
  use_x11_device = F,
  horizontal = F,
  main = "Estimated Precision matrix",
  center_value = NULL,
  col.upper = "black",
  col.center = "grey50",
  col.lower = "white"
)

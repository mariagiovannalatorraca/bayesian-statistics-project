# Libraries ---------------------------------------------------------------
library(BDgraph)
library(fda)


# Custom function ---------------------------------------------------------
#' Generate_BlockDiagonal
#' This function generate data by constructing a block diagonal matrix. 
#' The idea is taken from Sun et at. (2014) -  Adaptive Variable Clustering in Gaussian Graphical Models
#' We add an edge between two variables with probability p_in if they are in the same cluster, otherwise with
#' probability p_out. The corresponding value of the precision matrix is set equal to elem_out. To ensure positive
#' definiteness, we increase its eigenvalues summing the absolute value of the minimum eigenvalue plus min_eigenval_correction
#' Given the precision matrix, the covariance matrix, the data and the graph are generate accordingly. 
#' Note, the precision matrix is not drawn from a G-Wishart distribution
#' @param n [scalar] sample size
#' @param z_true [vector] cluster allocation vector
#' @param p_in [scalar] probability of having an edge between variables in the same cluster
#' @param p_out [scalar] probability of having an edge between variables in the different cluster
#' @param elem_out [scalar] value of elements of the precision matrix outside the diagonal
#' @param min_eigenval_correction [scalar] value to be added to the minimum eigenvalue
#' @param seed [integer] random seed
#'
#' @export
Generate_BlockDiagonal = function(n,
                                  z_true,
                                  p_in = 1,
                                  p_out = 0,
                                  elem_out = 5,
                                  min_eigenval_correction = 3,
                                  seed = 25254) {
  
  set.seed(seed)
  p = length(z_true) #set p
  # Generate precision and covariance matrices
  Omega = matrix(0,nrow = p,ncol = p)
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      
      if(z_true[i]==z_true[j]){
        # draw with prob p_in
        if(runif(n=1)<p_in)
          Omega[i,j] = elem_out
      }else{
        # draw with prob p_out
        if(runif(n=1)<p_out)
          Omega[i,j] = elem_out
      }
      
    }
  }
  Omega = Omega + t(Omega) # define lower diagonal part
  min_eig = eigen(Omega)$values[p] # compute eigenvalues, may be negative or 0
  Omega_true = Omega + diag(x=abs(min_eig) + min_eigenval_correction,nrow = p)  # ensure positive definiteness
  Sigma_true = solve(Omega_true) # compute covariance matrix
  
  # Generate data
  data = matrix(NA,nrow = n, ncol = p)
  data = t(apply(data, 1, FUN = function(x){ mvtnorm::rmvnorm(n = 1, sigma = Sigma_true)  }))
  U = t(data)%*%data
  
  # Compute underlying graph
  G = Omega_true
  G[abs(G)>0] = 1
  
  return(list("data"=data,"U"=U,"Prec"=Omega_true,"Cov"=Sigma_true,"Graph"=G))
}


# Data simulation ---------------------------------------------------------

# Set the seed
seed = 22111996
# Simulate a block graph
set.seed(seed)
p <- 40
Nclust = 3 # number of clusters
n_j = rep(floor(p/Nclust),Nclust) # cluster cardinalities
n_j[Nclust] = n_j[Nclust] + (p-sum(n_j)) # make sure n_j sums to p

z_true = rep(0,p)
inizio = c(1,cumsum(n_j)[1:(Nclust-1)]+1)
fine   = cumsum(n_j)
for(j in 1:Nclust){
  z_true[ inizio[j]:fine[j] ] = j
}
sim_KG = Generate_BlockDiagonal(n = p, z_true = z_true, p_in = 1, p_out = 0, seed = seed,
                                elem_out = 100, min_eigenval_correction = 100)
K_true = sim_KG$Prec
ACutils::ACheatmap(
  K_true,  
  use_x11_device = F,
  horizontal = F,
  main = "True Precision matrix",
  center_value = NULL,
  col.upper = "black",
  col.center = "grey50",
  col.lower = "white"
)
ACutils::ACheatmap(
  sim_KG$Cov,  
  use_x11_device = F,
  horizontal = F,
  main = "True Precision matrix",
  center_value = NULL,
  col.upper = "black",
  col.center = "grey50",
  col.lower = "white"
)


# Simulate the beta coefficients given the precision matrix
n <- 300
x_mu <- seq(0, 1, length.out=p)
mu_true = ACutils::dmix(x = x_mu, w_j = c(0.9,0.1), mu_vec = c(0.3,0.8), sigma_vec = c(0.1,0.05))
plot(x_mu, mu_true, type='l')
beta_true <- rmvnorm(n = n, mean =  mu_true , sigma=sim_KG$Cov)
plot(beta_true[1,], type='l')


# Build the phi matrix
r <- 200
x <- seq(0, 1, length.out=r)
basis <- create.bspline.basis(rangeval=range(x), nbasis=40, norder=3) 
BaseMat <- eval.basis(x, basis) # matrix r x p 

# Simulate the functional data 
tau_eps = 100
cov_y <- (1/tau_eps) * diag(1, r, r)
y_hat_true = matrix(nrow = n, ncol = r)

for (i in c(1:n)){
  mean_y <- as.vector(BaseMat%*%beta_true[i, ])
  y_hat_true[i, ] <- rmvnorm(n=1, mean=mean_y, sigma=cov_y)
}

plot(x, y_hat_true[1, ], type='l')

# Save data in a cvs
BaseMat_path = "BaseMat.csv"
y_hat_true_path = "y_hat_true.csv"
beta_true_path = "beta_true.csv"
mu_true_path = "C:/Users/volor/Desktop/Bayesian Statistics/Project/R/mu_true.csv"

saveRDS(sim_KG, file ="simKG.rds")

write.csv(BaseMat, file = BaseMat_path, row.names = FALSE)
write.csv(y_hat_true, file = y_hat_true_path, row.names = FALSE)
write.csv(beta_true, file = beta_true_path, row.names = FALSE)
write.csv(mu_true, file = mu_true_path, row.names = FALSE)

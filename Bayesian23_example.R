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

## Generate data: define settings

# Define true clustering
rho_true = c(8,4,8,5)

# Set seed for data generation
seed_data_gen = 22111996

# Define number of observations
n = 500

# Define variance of the Beta
beta_sig2 = 1/16

z_true = rho_to_z(rho_true)
r_true = z_to_r(z_true)
p = length(z_true)

# Generate data
sim = Generate_BlockDiagonal(n = n, z_true = z_true)


ACutils::ACheatmap(
  sim$Graph,
  use_x11_device = F,
  horizontal = F,
  main = "Graph",
  center_value = NULL,
  col.upper = "black",
  col.center = "grey50",
  col.lower = "white"
)

ACutils::ACheatmap(
  sim$Prec,
  use_x11_device = F,
  horizontal = F,
  main = "Precision matrix",
  center_value = NULL,
  col.upper = "black",
  col.center = "grey50",
  col.lower = "white"
)

ACutils::ACheatmap(
  sim$Cov,
  use_x11_device = F,
  horizontal = F,
  main = "Covariance",
  center_value = NULL,
  col.upper = "black",
  col.center = "grey50",
  col.lower = "white"
)

## Compute graph density
graph_density = sum(sim$Graph) / (p*(p-1))



## Set options

options = set_options(  sigma_prior_0=0.5,
                        sigma_prior_parameters=list("a"=1,"b"=1,"c"=1,"d"=1),
                        theta_prior_0=1,
                        theta_prior_parameters=list("c"=1,"d"=1),
                        rho0=p, # start with one group
                        weights_a0=rep(1,p-1),
                        weights_d0=rep(1,p-1),
                        alpha_target=0.234,
                        beta_mu=graph_density,
                        beta_sig2=beta_sig2,
                        d=3,
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
  data = sim$data,
  niter = 10,
  nburn = 2,
  thin = 1,
  options = options,
  seed = 123456,
  print = TRUE
)

## Save

# Before saving, also append true data
res$true_rho = rho_true
res$true_precision = sim$Prec
res$true_graph = sim$Graph
# remove self loops
res$true_graph[col(res$true_graph)==row(res$true_graph)] = 0

# save an object to a file
# saveRDS(res, file = filename_data)

# Posterior analysis

## Recomputing the partition in other forms and the number of groups

rho_true = res$true_rho
r_true = rho_to_r(rho_true)
z_true = rho_to_z(rho_true)
p = length(z_true)
num_clusters_true = length(rho_true)
rho = res$rho
r = do.call(rbind, lapply(res$rho, rho_to_r))
z = do.call(rbind, lapply(res$rho, rho_to_z))
num_clusters = do.call(rbind, lapply(res$rho, length))
num_clusters = as.vector(num_clusters)


### Acceptance frequency
mean(res$accepted)

### Barplot of changepoints
bar_heights = colSums(r)
cp_true = which(r_true==1)
color <- ifelse(seq_along(bar_heights) %in% c(cp_true), "red", "gray")

barplot(
  bar_heights,
  names = seq_along(bar_heights),
  border = "NA",
  space = 0,
  yaxt = "n",
  main="Changepoint frequency distribution",
  #col = color,
  cex.names=.6,
  las=2
)

abline(v=cp_true-0.5, col="red", lwd=2)
legend("topright", legend=c("True"), col=c("red"),
       bty = "n",
       lty = 1,
       cex = 0.6)

### Evolution of the number of clusters
plot(
  x = seq_along(num_clusters),
  y = num_clusters,
  type = "n",
  xlab = "Iterations",
  ylab = "Number of groups",
  main = "Number of groups - Traceplot"
)
lines(x = seq_along(num_clusters), y = num_clusters)
abline(h = length(z_to_rho(z_true)),
       col = "red",
       lwd = 4)
legend("topleft", legend=c("True"), col=c("red"),
       lty = 1,
       cex = 1)

barplot(
  prop.table(table(num_clusters)),
  xlab = "Number of groups",
  ylab = "Relative Frequency",
  main = paste(
    "Number of groups - Relative Frequency\n",
    "Last:",
    tail(num_clusters, n = 1),
    "- Mean:",
    round(mean(num_clusters), 2),
    "- True:",
    num_clusters_true
  )
)

### Evolution of the Rand Index

# computing rand index for each iteration
rand_index = apply(z, 1, mcclust::arandi, z_true)

# plotting the traceplot of the index
plot(
  x = seq_along(rand_index),
  y = rand_index,
  type = "n",
  xlab = "Iterations",
  ylab = "Rand Index",
  main = paste(
    "Rand Index - Traceplot\n",
    "Last:",
    round(tail(rand_index, n=1), 3),
    "- Mean:",
    round(mean(rand_index), 2)
  )
)
lines(x = seq_along(rand_index), y = rand_index)
abline(h = 1, col = "red", lwd = 4)

### Retrieving best partition using VI on visited ones (order is guaranteed here)

# Here we are satisfied with finding the optimal partition only in the set of those visited, not in all the possible ones. I expect it could work even worse. But at least it guarantees to find an admissible one.
# I would say that it is the implementation of formula (13) of the Corradin-Danese paper (https://doi.org/10.1016/j.ijar.2021.12.019).

# compute VI
sim_matrix <- salso::psm(z)
dists <- VI_LB(z, psm_mat = sim_matrix)

# select best partition (among the visited ones)
best_partition_index = which.min(dists)
rho_est = rho[[best_partition_index]]
z_est = z[best_partition_index,]

# VI loss
dists[best_partition_index]

# select best partition
unname(z_est)

# compute Rand Index
mcclust::arandi(z_est, z_true)

## Graph

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

### Standardized Hamming distance
SHD = sum(abs(sim$true_graph - G_est)) / (p^2 - p)
SHD

### Plot estimated matrices
ACutils::ACheatmap(
  last_plinks,
  use_x11_device = F,
  horizontal = F,
  main = "Estimated plinks matrix",
  center_value = NULL,
  col.upper = "black",
  col.center = "grey50",
  col.lower = "white"
)

ACutils::ACheatmap(
  G_est,
  use_x11_device = F,
  horizontal = F,
  main = "Estimated Graph",
  center_value = NULL,
  col.upper = "black",
  col.center = "grey50",
  col.lower = "white"
)

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

### Evolution of the Kullback-Leibler
kl_dist = do.call(rbind, lapply(res$K, function(k) {
  ACutils::KL_dist(res$true_precision, k)
}))

last = round(tail(kl_dist, n=1), 3)
plot(
  x = seq_along(kl_dist),
  y = kl_dist,
  type = "n",
  xlab = "Iterations",
  ylab = "K-L distance",
  main = paste("Kullback-Leibler distance\nLast value:", last)
)
lines(x = seq_along(kl_dist), y = kl_dist)



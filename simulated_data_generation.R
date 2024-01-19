library(BDgraph)
library(fda)

# Set the seed
seed = 22111996

# Simulate a block graph
set.seed(seed)
p <- 40
adj <- graph.sim(p = p, graph = "cluster", class = 3 , vis = T, prob = 1)
adj

# Simulate the precision matrix given the graph
set.seed(seed)
K_true <- rgwish(n = 1, adj = adj, D = 0.001*diag(p))

# Simulate the beta coefficients given the precision matrix
set.seed(seed)
n <- 300
x_mu <- seq(0, 1, length.out=p)
mu_true = dgamma(x = x_mu, shape = 2, rate = 10) 
plot(x_mu, mu_true, type='l')
beta_true <- rmvnorm(n = n, mean =  mu_true , sigma=solve(K_true))
plot(beta_true[1,], type='l')

# Build the phi matrix
r <- 200
x <- seq(0, 1, length.out=r)
basis <- create.bspline.basis(rangeval=range(x), nbasis=40, norder=3) 
BaseMat <- eval.basis(x, basis) # matrix r x p 

# Simulate the functional data 
tau_eps = 100
cov_y <- 1/tau_eps * diag(1, r, r)
y_hat_true = matrix(nrow = n, ncol = r)

for (i in c(1:n)){
  mean_y <- as.vector(BaseMat%*%beta_true[i, ])
  y_hat_true[i, ] <- rmvnorm(n=1, mean=mean_y, sigma=cov_y)
}

x11()
plot(x, y_hat_true[1, ], type='l')

# Save data in a cvs
BaseMat_path = "BaseMat.csv"
y_hat_true_path = "y_hat_true.csv"

write.csv(BaseMat, file = BaseMat_path, row.names = FALSE)
write.csv(y_hat_true, file = y_hat_true_path, row.names = FALSE)


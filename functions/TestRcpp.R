# Load Data
setwd("G:/.shortcut-targets-by-id/1Ck2MctcmCBWOueSMeW4Zr8IOvOaOzdqz/BicoccaDrive/TutoratoPolimi_2324/Materiale/FunctionalRandomGraph_shared/Codici")
load("purees.Rdat")
curves = purees$data
wavelengths = purees$wavelengths
strawberry <- curves[which(curves$Group == "Strawberry"), ]
data = strawberry[,-1]*100
data = as.matrix(data) # n x r

# Generate Basis
p =  40
n = dim(data)[1]
r = dim(data)[2]
range_x = range(wavelengths)
basis = BGSL::Generate_Basis(n_basis = p, range = range_x, grid_points = wavelengths, order =3)
BaseMat = basis$BaseMat # r x p matrix

# Perform 1 iteration
seed = 1234
set.seed(seed)
Beta = matrix(rnorm(n=p*n), nrow = p, ncol = n)
mu = rnorm(n=p)
tau_eps = 100
K = rWishart(n=1,df = p+10, Sigma = diag(p) )
K = K[,,1]
# ACutils::ACheatmap(K,center_value = NULL, remove_diag = T)
tbase_base = t(BaseMat)%*%BaseMat # p x p
tbase_data = t(BaseMat)%*%t(data) # p x n
Sdata = sum(diag(data%*%t(data))) # inefficient calculation
Update_Beta <- TRUE
Update_Mu <- TRUE
Update_Tau <- TRUE
a_tau_eps <- b_tau_eps <- 2
sigma_mu = 100

Rcpp::sourceCpp('UpdateParams.cpp')
fit = UpdateParams( Beta, mu, tau_eps, K,
                    tbase_base, tbase_data, Sdata,
                    a_tau_eps, b_tau_eps, sigma_mu, r,
                    Update_Beta, Update_Mu, Update_Tau )

fit$Beta[,1]
fit$mu
fit$tau_eps
dim(fit$U)

Rcpp::sourceCpp('UpdateParamsGSL.cpp')
fit = UpdateParamsGSL( Beta, mu, tau_eps, K,
                        tbase_base, tbase_data, Sdata, 
                        a_tau_eps, b_tau_eps, sigma_mu, r,
                        Update_Beta, Update_Mu, Update_Tau )

fit$Beta[,1]
fit$mu
fit$tau_eps
dim(fit$U)

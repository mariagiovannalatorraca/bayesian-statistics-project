// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


//' Update of regression parameters
//'
//' @param data [matrix], r x n matrix (not used in this function)
//' @param Beta [matrix], p x n matrix, where p is the number of nodes and n the number of curves
//' @param mu [vector], length p
//' @param tau_eps [double]
//' @param K [matrix], p x p. This is a precision matrix, but be symmetric and positive definite
//' @param tbase_base [matrix], p x p. This contains \Phi^T * \Phi, where \Phi is a r x p matrix of spline coefficients. r is the number of grid points (the same grid for all splines is assumed)
//' @param tbase_data [matrix], p x n. This is  \Phi^T * data, such that tbase_data.col(i) = \Phi^T * y_i
//' @param Sdata [double] that contains the sum of the inner products for all data y_i, \sum_{i=1}^n (y_i,y_i) = \sum_{i=1}^n y_i^T * y_i
//' @param a_tau_eps [double] prior shape hyperparameter for tau_eps
//' @param b_tau_eps [double] prior rate hyperparameter for tau_eps
//' @param sigma_mu [double] variance hyperparameter for mu
//' @param r [int] r is the number of grid points (the same grid for all splines is assumed)
//' @param UpdateBeta [bool] set TRUE to update Beta matrix
//' @param UpdateMu [bool] set TRUE to update mu vector matrix
//' @param UpdateTau [bool] set TRUE to update tau_eps
// [[Rcpp::export]]
Rcpp::List UpdateParams(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& Beta, 
                        Eigen::VectorXd& mu, // p-length vector
                        double& tau_eps, // scalar
                        const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& K, 
                        const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& tbase_base, 
                        const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& tbase_data, 
                        const double& Sdata, const double& a_tau_eps, const double& b_tau_eps, const double& sigma_mu,
                        const unsigned int& r,
                        const bool& UpdateBeta, const bool& UpdateMu, const bool& UpdateTau
                        ){

  // typedef
  using IdxType         = std::size_t;
  using MatRow          = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>; 
  using MatCol          = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>; 
  using VecRow          = Eigen::RowVectorXd;
  using VecCol          = Eigen::VectorXd;
  using CholTypeRow     = Eigen::LLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Lower>;
  using CholTypeCol     = Eigen::LLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>, Eigen::Lower>;
  
  // initilizations and basic computations
  const unsigned int p = Beta.rows(); // number of nodes
  const unsigned int n = Beta.cols(); // number of curves
  const MatRow Irow( MatRow::Identity(p,p) ); // Identity matrix
  const MatRow one_over_sigma_mu((1/sigma_mu)*Irow); 

  Function rnorm_R("rnorm"); // rnorm function from R
  Function rgamma_R("rgamma");  // rgamma function from R
  
  //1) Update mu
  VecCol S_beta( Beta.rowwise().sum() ); // \sum_{i=1}^n \beta_i
  CholTypeRow chol_invM( one_over_sigma_mu + n * K ); // chol_invM is such that chol_invM*t(chol_invM) = invM
  MatRow M( chol_invM.solve(Irow) ); // M = (invM)^-1, full conditional Covariance matrix
  VecCol m( M*(K*S_beta) ); // full conditional mean

  if(UpdateMu){
    NumericVector z1( rnorm_R( Named("n",p), Named("mean",0), Named("sd")=1 ) ); // p iid draws from N(0,1)
    Eigen::Map<Eigen::VectorXd> z1_eigen(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(z1)); // cast from Rcpp::NumericVector to Eigen VectorXd

    mu = m + M.llt().matrixL().transpose()*z1_eigen; // posterior sample for \mu
  }


  //2) Update Beta
  CholTypeRow chol_invBn( tau_eps * tbase_base + K ); // \tau_epsilon * \Phi^T * \Phi + K
  MatRow Bn( chol_invBn.solve(Irow) ); // full conditional Covariance matrix
  MatRow chol_Bn( Bn.llt().matrixL() ); // lower trinagular cholesky decomposition of covariance matrix
  VecCol Kmu( K*mu ); // K*\mu
  MatCol U(MatCol::Zero(p,p)); // p x p matrix that will store sum_{i=1}^n( (\beta_i - \mu) * (\beta_i - \mu)^T  ) 
  double b_tau_eps_post(Sdata + 2.0 * b_tau_eps); // inizialize posterior rate parameter for tau_eps 

  for(unsigned int i = 0; i < n; i++){
    VecCol bn_i = Bn*(tau_eps*tbase_data.col(i) + Kmu); // full conditional mean

    VecCol beta_i; // must be defined outside if-else scope
    if(UpdateBeta){
      NumericVector z2( rnorm_R( Named("n",p), Named("mean",0), Named("sd")=1) ); // p iid draws from N(0,1)
      Eigen::Map<Eigen::VectorXd> z2_eigen(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(z2)); // cast from Rcpp::NumericVector to Eigen VectorXd
      beta_i = bn_i + chol_Bn.transpose()*z2_eigen;  // posterior sample for \beta_i  
      Beta.col(i) = beta_i; // save result in Beta matrix    
    }
    else{
      beta_i = Beta.col(i); // take the old value 
    }

    U += (beta_i - mu)*(beta_i - mu).transpose(); // Update U matrix
    b_tau_eps_post += beta_i.dot(tbase_base*beta_i) - 2.0*beta_i.dot(tbase_data.col(i));  // update posterior rate parameter for tau_eps
  }

  //Precision tau_epsilon
  double a_tau_eps_post = (n*r + a_tau_eps)*0.5; // compute shape full conditional
  b_tau_eps_post /= 2.0; 

  if(UpdateTau){
    NumericVector tau_eps_rcpp = rgamma_R( Named("n",1), Named("shape",a_tau_eps_post), Named("rate")=b_tau_eps_post ); // 1 posterior sample from tau_epsilon
    tau_eps = tau_eps_rcpp[0]; // cast as a double    
  }


  // Create list to be returned
  return Rcpp::List::create(  Rcpp::Named("Beta") = Beta,
                              Rcpp::Named("mu") = mu,
                              Rcpp::Named("tau_eps") = tau_eps,
                              Rcpp::Named("U") =  U
                            );
  
}



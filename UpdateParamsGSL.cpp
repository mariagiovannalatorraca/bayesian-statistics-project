// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppGSL)]]
#include <Rcpp.h>
#include <RcppEigen.h>

// GSL
#include <RcppGSL.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>     //For random number generators
#include <gsl/gsl_randist.h> //For random variates and probability density functions
#include <gsl/gsl_cdf.h>   //For cumulative density functions
#include <gsl/gsl_linalg.h> //For cholesky decomposition
#include <gsl/gsl_sf.h>     //For special functions such as factorials


//Containers
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <tuple>

//Generic
#define _USE_MATH_DEFINES
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <limits>
#include <memory>
#include <exception>
#include <numeric>
#include <iterator> //for std::inserter
#include <utility>  //for std::forward
#include <tuple>
#include <type_traits>
#include <functional>

//Eigen
#include <Eigen/Dense>
#include <Eigen/Cholesky>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector colNorm(const RcppGSL::Matrix & G) {
  int k = G.ncol();
  Rcpp::NumericVector n(k);           // to store results
  for (int j = 0; j < k; j++) {
    RcppGSL::VectorView colview = gsl_matrix_const_column (G, j);
    n[j] = gsl_blas_dnrm2(colview);
  }
  return n;                           // return vector
}



namespace sample{ //use the sample:: namespace to avoid clashes with R or other packages

  /*--------------------------------------------------------------
    Random number generator wrapper
  ----------------------------------------------------------------*/

  //This class simply wraps in c++ code the construction and desctruction of a gsl_rng object.
  //I had to remove std::random_device because there is a bug when compiling in window (returs always the same value).
  // reference for bug -> https://en.cppreference.com/w/cpp/numeric/random/random_device, https://sourceforge.net/p/mingw-w64/bugs/338/
  class GSL_RNG{
    public:

      //constructor 1: takes one seed. If seed is 0, generates random seed
      GSL_RNG(unsigned int const & _seed){
        if(_seed == 0){
          seed = static_cast<unsigned int>(std::chrono::steady_clock::now().time_since_epoch().count());
          std::seed_seq seq = {seed}; //seed provived here has to be random. Than std::seed_seq adds entropy becasuse steady_clock is not sufficientyl widespread
          std::vector<unsigned int> seeds(1);
          seq.generate(seeds.begin(), seeds.end());
          seed = seeds[0];
          //std::cout<<"seed = "<<seed<<std::endl;
        }
        else{
          seed = _seed;
        }
        gsl_rng_env_setup();
        r = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_set(r,seed);
      }

      //constructor 0: default constructor. It is equvalent to the previous one using seed=0.
      GSL_RNG(){
        gsl_rng_env_setup();
        r = gsl_rng_alloc(gsl_rng_default);
        seed = static_cast<unsigned int>(std::chrono::steady_clock::now().time_since_epoch().count());
        std::seed_seq seq = {seed}; //seed provived here has to be random. Than std::seed_seq adds entropy becasuse steady_clock is not sufficientyl widespread
        std::vector<unsigned int> seeds(1);
        seq.generate(seeds.begin(), seeds.end());
        seed = seeds[0];
        gsl_rng_set(r,seed);
      }

      //descructor, can not be removed because gsl_rng_alloc() is allocating memory
      ~GSL_RNG(){
        gsl_rng_free(r);
      }

      //print some information
      void print_info()const{
        printf ("generator type: %s\n", gsl_rng_name(r));
        std::cout<<"seed = "<<seed<<std::endl;
      }

      //call operator, return the engine ready to generate a number
      gsl_rng* operator()()const{
        return r;
      }

      //set seed, not tested much
      inline void set_seed(unsigned int const & s){
        seed = s;
        gsl_rng_set(r,seed);
      }

      //getter
      inline unsigned int get_seed() const{
        return seed;
      }
    private:
      gsl_rng * r; //the random number generator
      unsigned int seed; //the seed
  };

  /*--------------------------------------------------------------
    Random number distribution wrappers
  ----------------------------------------------------------------*/

  //Callable object to draw a sample from sampling from Unif([0,1])
  struct runif
  {
    double operator()(GSL_RNG const & engine)const{
      return gsl_rng_uniform(engine()); //gsl_rng_uniform is a function, nothing has to be de-allocated
    }
    double operator()()const{
      return this->operator()(GSL_RNG ());
      //return runif()(GSL_RNG ())

      /*GSL_RNG () creates a GSL_RNG obj calling the default constructor and than calls it call operator.
      In other words, it generates a random number generator, generates a number and destroys the generator.
      It is equivalent to return gsl_rng_uniform( GSL_RNG ()() ).
      In this second case, one () is for the default constructor, the second () is for the call operator.  */
    }
  };

  //Callable object to draw a sample from Unif({0,...,N-1})
  struct runif_int
  {
    unsigned int operator()(GSL_RNG const & engine, unsigned int const & N)const{
      return gsl_rng_uniform_int(engine(), N); //gsl_rng_uniform_int is a function, nothing has to be de-allocated.
    }
    unsigned int operator()(unsigned int const & N)const{
      return this->operator()(GSL_RNG (), N);
      /*GSL_RNG () creates a GSL_RNG obj calling the default constructor and than calls it call operator.
       In other words, it generates a random number generator, generates a number and destroys the generator*/
    }
  };
  

  //Callable object to draw a sample from N(mean,sd).
  // --> NB  it takes the standard deviation as input! <--
  struct rnorm
  {
    //Gets the engine
    //N(mean,sd)
    double operator()(GSL_RNG const & engine, double const & mean, double const & sd)const{
      return gsl_ran_gaussian_ziggurat(engine(),sd) + mean;
    }

    //Gets the engine
    //N(0,1)
    double operator()(GSL_RNG const & engine)const{
      return gsl_ran_gaussian_ziggurat(engine(), 1.0);
    }

    //Engine defaulted
    //N(mean,sd)
    double operator()(double const & mean, double const & sd)const{
      return this->operator()(GSL_RNG (), mean, sd);
    }

    //Engine defaulted
    //N(0,1)
    double operator()()const{
      return gsl_ran_gaussian_ziggurat(GSL_RNG ()(),1.0); //the first () is for the constructor, the second il for the call operator
    }
  };

  //Callable object to draw a sample from Gamma(shape,scale).
  // --> NB  Watch out notation! gsl uses the scale parameter, in Europe we are used to the rate parameter (rate = 1/scale) <--
  struct rgamma{

    //Gets the engine
    //Gamma(shape,scale)
    double operator()(GSL_RNG const & engine, double const & shape, double const & scale)const{
      return gsl_ran_gamma(engine(),shape,scale);
    }

    //Engine defaulted
    //Gamma(shape,scale)
    double operator()(double const & shape, double const & scale)const{
      return gsl_ran_gamma(GSL_RNG ()(),shape,scale);
    }
  };

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
Rcpp::List UpdateParamsGSL(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& Beta, 
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

  sample::GSL_RNG engine(0); // initialize random engine with default random seed
  sample::rnorm rnorm; // define callable object to generate random samples from a 1d Gaussian
  sample::rgamma rgamma; // watch out, this is parametrized as shape and scale! scale is 1/rate.
  
  //1) Update mu
  VecCol S_beta( Beta.rowwise().sum() ); // \sum_{i=1}^n \beta_i
  CholTypeRow chol_invM( one_over_sigma_mu + n * K ); // chol_invM is such that chol_invM*t(chol_invM) = invM
  MatRow M( chol_invM.solve(Irow) ); // M = (invM)^-1, full conditional Covariance matrix
  VecCol m( M*(K*S_beta) ); // full conditional mean

  if(UpdateMu){
    VecCol z1_eigen{VecCol::Constant(p,0.0)}; 

    // rnorm only draws 1 variable at the time. A for loop is needed to generate p variables
    for(unsigned int jj = 0; jj < p; jj++)
      z1_eigen(jj) = rnorm(engine, 0.0, 1.0);

    //NumericVector z1( rnorm_R( Named("n",p), Named("mean",0), Named("sd")=1 ) ); // p iid draws from N(0,1)
    //Eigen::Map<Eigen::VectorXd> z1_eigen(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(z1)); // cast from Rcpp::NumericVector to Eigen VectorXd

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
      VecCol z2_eigen{VecCol::Constant(p,0.0)}; 

      // rnorm only draws 1 variable at the time. A for loop is needed to generate p variables
      for(unsigned int jj = 0; jj < p; jj++)
        z2_eigen(jj) = rnorm(engine, 0.0, 1.0);

      //NumericVector z2( rnorm_R( Named("n",p), Named("mean",0), Named("sd")=1) ); // p iid draws from N(0,1)
      //Eigen::Map<Eigen::VectorXd> z2_eigen(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(z2)); // cast from Rcpp::NumericVector to Eigen VectorXd
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
    //NumericVector tau_eps_rcpp = rgamma_R( Named("n",1), Named("shape",a_tau_eps_post), Named("rate")=b_tau_eps_post ); // 1 posterior sample from tau_epsilon
    //tau_eps = tau_eps_rcpp[0]; // cast as a double    

    tau_eps = rgamma(engine, a_tau_eps_post, 1.0/b_tau_eps_post);

  }


  // Create list to be returned
  return Rcpp::List::create(  Rcpp::Named("Beta") = Beta,
                              Rcpp::Named("mu") = mu,
                              Rcpp::Named("tau_eps") = tau_eps,
                              Rcpp::Named("U") =  U
                            );
  
}

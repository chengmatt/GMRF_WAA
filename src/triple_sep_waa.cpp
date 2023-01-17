// Purpose: To implement triple separability as a weight-at-age matrix
// for EBS pollock
// Creator: Matthew LH. Cheng (UAF-CFOS)
// Date: 1/15/23

#include<TMB.hpp>

template<class Type> 
Type objective_function<Type>::operator() () 
  {
  
  using namespace density; // Define namespace to use multivariate distributions
  using namespace Eigen; // Define namespace for Eigen functions (i.e., sparse matrix)
  
  // Define TMB dimensions
  DATA_VECTOR(years); // vector of number of years
  DATA_VECTOR(ages); // vector of ages
  int n_years = years.size(); // integer of number of years
  int n_ages = ages.size(); // integer for number of ages
  int total_n = n_years * n_ages; // integer for nyears * nages
  
  // Define matrices for weight-at-age 
  DATA_MATRIX(X_at); // n_years * n_ages
  DATA_MATRIX(Xse_at); // n_years * n_ages
  
  // Index matrix to loop through to consruct precision matrix
  DATA_MATRIX(ay_Index);  // (n_years * n_ages) * 2 
  
  // Define parameters and precision matrix here
  PARAMETER(rho_a); // Correlation by age
  PARAMETER(rho_y); // Correlation by year
  PARAMETER(rho_c); // Correlation by cohort
  PARAMETER(log_sigma2); // Variance of the GMRF process
  PARAMETER( ln_L0 );
  PARAMETER( ln_Linf );
  PARAMETER( ln_k );
  PARAMETER( ln_alpha );
  PARAMETER( ln_beta );
  PARAMETER_ARRAY(Y_at); // Random effects for weight-at-age vector
  
  // Construct matrix for B (neighborhood)
  matrix<Type> B(total_n,total_n); 
  // Construct matrix for Identity Matrix
  matrix<Type> I(total_n,total_n); 
  // Construct matrix for Omega (Variances)
  matrix<Type> Omega(total_n,total_n); 
  // Construct dense precision matrix here
  matrix<Type> Q_dense(total_n,total_n);
  
  // CREATE PRECISION MATRIX ---------------------------------------
  for(int n = 0; n < total_n; n++) {
    
    // Define year and age objects
    Type age = ay_Index(n,0);
    Type year = ay_Index(n,1); 

    if(age > 1) {
      
      // Get column index for years
      for(int n1 = 0; n1 < total_n; n1++) {
        if(ay_Index(n1, 0) == age - 1 && ay_Index(n1, 1) == year) 
        B(n, n1) = rho_y;
      } // figure out index where this condition is met and input correlation parameter
      
    } // end age > 1 
    
    if(year > 1) {
      
      // Get column index for years
      for(int n1 = 0; n1 < total_n; n1++) {
        if(ay_Index(n1, 0) == age && ay_Index(n1, 1) == year - 1) 
          B(n, n1) = rho_a;
      } // figure out index where this condition is met and input correlation parameter
      
    } // if year > 1 
    
    if(year > 1 && age > 1) {
      
      // Get column index for years
      for(int n1 = 0; n1 < total_n; n1++) {
        if(ay_Index(n1, 0) == age - 1 && ay_Index(n1, 1) == year - 1) 
          B(n,n1) = rho_c; // correlation by cohort
      } // figure out index where this condition is met and input correlation parameter
      
    } // if both year and age > 1
    
  } // end n loop
  
  // Fill in identity matrix diagonals with 1s
  for(int i = 0; i < total_n; i++) {
    for(int j = 0; j < total_n; j++) {
      if(i == j) I(i,j) = Type(1.0);
        else I(i,j) = Type(0.0);
    } // j loop
  } // i loop
  
  // Fill in identity matrix diagonals with variance (not sure this is correct...)
  for(int i = 0; i < total_n; i++) {
    for(int j = 0; j < total_n; j++) {
      if(i == j) Omega(i,j) = Type(1/exp(log_sigma2));
      else Omega(i,j) = Type(0.0);
    } // j loop
  } // i loop
  
  // Now, do calculations to construct (Q = (I - t(B)) %*% Omega %*% (I-B))
  Q_dense = (I - B.transpose()) * Omega * (I-B);
  // Coerce to sparse precision matrix
  Eigen::SparseMatrix< Type > Q_sparse = asSparseMatrix(Q_dense);
  
  // Define joint negative likelihood here
  Type jnLL = 0;

  // LIKELIHOOD SECTION ------------------------
  array<Type> mu_at(Y_at.rows(), Y_at.cols());
  Type L0 = exp(ln_L0);
  Type Linf = exp(ln_Linf);
  Type k = exp(ln_k);
  Type alpha = exp(ln_alpha);
  Type beta = exp(ln_beta);
  for(int a = 0; a < X_at.rows(); a++) {
  for(int t = 0; t < X_at.cols(); t++) {
    mu_at(a,t) = Linf - (Linf-L0)*exp(-k*Type(a));
    mu_at(a,t) = alpha * pow( mu_at(a,t), beta );
  }}

  // Evaluate WAA data likelihood
  for(int a = 0; a < X_at.rows(); a++) {
  for(int t = 0; t < X_at.cols(); t++) {
    jnLL -= dnorm(Y_at(a,t), X_at(a,t), Xse_at(a,t), true);
  }}

  // Evaluate GMRF with precision matrix estimating cohort, year, and age correlations
  array<Type> eps_at(Y_at.rows(), Y_at.cols());
  eps_at = Y_at - mu_at;
  jnLL += GMRF(Q_sparse)( eps_at );

  // REPORT SECTION ------------------------
  REPORT(jnLL);
  REPORT(Q_sparse);
  REPORT(I);
  REPORT(B);
  REPORT(Omega);
  REPORT( mu_at );

  return(jnLL);
  
}
  
  

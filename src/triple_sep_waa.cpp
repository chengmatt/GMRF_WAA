// Purpose: To implement triple separability as a weight-at-age matrix
// for EBS pollock
// Creator: Matthew LH. Cheng (UAF-CFOS)
// Date: 1/15/23

#include<TMB.hpp>

template<class Type>
// @description: Function that constructs a precision matrix, separable along the
// year, age, and cohort axis
Eigen::SparseMatrix<Type> construct_Q(int n_years, // Integer of years
                                      int n_ages, // Integer of ages
                                      matrix<Type> ay_Index, // Index matrix to construct
                                      Type rho_y, // Correlation by years
                                      Type rho_a, // Correlation by ages
                                      Type rho_c, // Correlation by cohort
                                      Type log_sigma2, // Variance parameter governing GMRF
                                      int Var_Param // Parameterization of Variance ==0 (Conditional), == 1(Marginal)
                                        ) {
  
  // Dimension to construct matrices
  int total_n = n_years * n_ages; 
  
  // Construct matrices for precision matrix
  Eigen::SparseMatrix<Type> B(total_n,total_n); // B matrix
  Eigen::SparseMatrix<Type> I(total_n,total_n); // Identity matrix
  Eigen::SparseMatrix<Type> Omega(total_n,total_n); // Omega matrix (variances)
  Eigen::SparseMatrix<Type> Q_sparse(total_n, total_n); // Precision matrix
  
  for(int n = 0; n < total_n; n++) {
    
    // Define year and age objects
    Type age = ay_Index(n,0);
    Type year = ay_Index(n,1); 
    
    // Constructing B matrix to determine where the correlation pars should go
    if(age > 1) {
      
      // Get column index for years
      for(int n1 = 0; n1 < total_n; n1++) {
        if(ay_Index(n1, 0) == age - 1 && ay_Index(n1, 1) == year) 
          B.coeffRef(n, n1) = rho_y;
      } // n1 loop
      
    } // end age > 1 
    
    if(year > 1) {
      
      // Get column index for years
      for(int n1 = 0; n1 < total_n; n1++) {
        if(ay_Index(n1, 0) == age && ay_Index(n1, 1) == year - 1) 
          B.coeffRef(n, n1) = rho_a;
      } // n1 loop
      
    } // if year > 1 
    
    if(year > 1 && age > 1) {
      
      // Get column index for years
      for(int n1 = 0; n1 < total_n; n1++) {
        if(ay_Index(n1, 0) == age - 1 && ay_Index(n1, 1) == year - 1) 
          B.coeffRef(n,n1) = rho_c; // correlation by cohort
      } // n1 loop 
      
    } // if both year and age > 1
    
  } // end n loop
  
  // Fill in identity matrix diagonals with 1s
  for(int i = 0; i < total_n; i++) {
    for(int j = 0; j < total_n; j++) {
      if(i == j) I.coeffRef(i,j) = Type(1.0);
      else I.coeffRef(i,j) = Type(0.0);
    } // j loop
  } // i loop
  
  // Fill in Omega matrix here (variances)
  if(Var_Param == 0) { // Conditional variance
    
    for(int i = 0; i < total_n; i++) {
      for(int j = 0; j < total_n; j++) {
        if(i == j) Omega.coeffRef(i,j) = 1/exp(log_sigma2);
        else Omega.coeffRef(i,j) = Type(0.0);
      } // j loop
    } // i loop
    
  } // end if conditional variance
  
  // Need to solve Omega recursively to derive stationary (marginal) variance
  if(Var_Param == 1) { // Marginal Variance
    
    // Construct container objects
    matrix<Type> L(total_n, total_n); // L Matrix
    matrix<Type> tmp_I = I; // Temporary Identity Matrix
    matrix<Type> tmp_B = B; // Temporary B Matrix
    vector<Type> d(total_n); // Store variance calculations
    L = (tmp_I-tmp_B).inverse(); // Solve for L (not sure how to invert 
    // a sparse matrix, which is why I coerced them to regular matrices)

    for(int n = 0; n < total_n; n++) {
      
      if(n == 0) {
        d(n) = exp(log_sigma2); // Our marginal variance parameter
      } else{
        
        Type cumvar; // Cumulative Variance container
        int tmp_n_idx = n; // temporary container to do seq_len() like indexing

        for(int n1 = 0; n1 < tmp_n_idx; n1++) { // calculate cumulative variance
          cumvar += L(n, n1) * d(n1) * L(n, n1); // sum across these
        } // n1 loop
        
        // Calculate diagonal values of Omega
        d(n) = (exp(log_sigma2) - cumvar ) / pow(L(n, n), 2);
          
      } // end else statement
      
    } // end n loop
    
    // Now fill in our diagonals for Omega
    for(int i = 0; i < total_n; i++) {
      for(int j = 0; j < total_n; j++) {
        if(i == j) Omega.coeffRef(i,j) = 1/d(i);
        else Omega.coeffRef(i,j) = Type(0.0);
      } // j loop
    } // i loop
    
  } // end if marginal variance
  
  // Now, do calculations to construct (Q = (I - t(B)) %*% Omega %*% (I-B))
  Eigen::SparseMatrix<Type> B_transpose = B.transpose(); // transpose B matrix 
  
  // Calculate Precision Matrix
  Q_sparse = (I - B_transpose) * Omega * (I-B);
  
  return(Q_sparse);
} // end construct_Q function


template<class Type> 
Type objective_function<Type>::operator() () 
  {
  
  using namespace density; // Define namespace to use multivariate distributions
  using namespace Eigen; // Define namespace for Eigen functions (i.e., sparse matrix)
  
  // DATA SECTION ---------------------------------------------
  
  // Define TMB dimensions
  DATA_VECTOR(years); // vector of number of years
  DATA_VECTOR(ages); // vector of ages
  int n_years = years.size(); // integer of number of years
  int n_ages = ages.size(); // integer for number of ages
  int total_n = n_years * n_ages; // integer for nyears * nages
  
  // Define matrices for weight-at-age 
  DATA_MATRIX(X_at); // n_years * n_ages (Observered WAA)
  DATA_MATRIX(Xse_at); // n_years * n_ages (Observed sd WAA)
  
  // Index matrix to loop through to consruct precision matrix
  DATA_MATRIX(ay_Index);  // (n_years * n_ages) * 2 
  DATA_INTEGER(Var_Param); // Variance parameterization of Precision Matrix
  // == 0 (Conditional), == 1(Marginal)
  
  // PARAMETER SECTION ----------------------------------------
  PARAMETER(rho_a); // Correlation by age
  PARAMETER(rho_y); // Correlation by year
  PARAMETER(rho_c); // Correlation by cohort
  PARAMETER(log_sigma2); // Variance of the GMRF process
  PARAMETER(ln_L0); // vonB length at age 0
  PARAMETER(ln_Linf); // vonB asymptotic length
  PARAMETER(ln_k); // vonB brody growth coefficient
  PARAMETER(ln_alpha); // WL relationship conversion factor
  PARAMETER(ln_beta); // WL relationship allometric scaling
  PARAMETER_ARRAY(Y_at); // Random Effects Weight-at-age array
  
  // Define precision matrix for GMRF
  Eigen::SparseMatrix<Type> Q_sparse(total_n, total_n); // Precision matrix
  
  // Construct precision matrix here
  Q_sparse = construct_Q(n_years, n_ages, ay_Index, 
                         rho_y, rho_a, rho_c, 
                         log_sigma2, Var_Param);
  
  // Define joint negative likelihood here
  Type jnLL = 0;

  // LIKELIHOOD SECTION ------------------------
  array<Type> mu_at(Y_at.rows(), Y_at.cols()); // Mean weight at age across time
  
  // Transform vonB and WL parameters
  Type L0 = exp(ln_L0);
  Type Linf = exp(ln_Linf);
  Type k = exp(ln_k);
  Type alpha = exp(ln_alpha);
  Type beta = exp(ln_beta);
  
  // vonB and allometric weight-length relationship
  for(int a = 0; a < X_at.rows(); a++) {
  for(int t = 0; t < X_at.cols(); t++) {
    mu_at(a,t) = Linf - (Linf-L0)*exp(-k*Type(a));
    mu_at(a,t) = alpha * pow( mu_at(a,t), beta );
  } // t loop
} // a loop

  // Evaluate WAA data likelihood
  for(int a = 0; a < X_at.rows(); a++) {
  for(int t = 0; t < X_at.cols(); t++) {
    jnLL -= dnorm(X_at(a,t), Y_at(a,t), Xse_at(a,t), true);
  } // t loop
} // a loop

  // Evaluate GMRF with precision matrix estimating cohort, year, and age correlations
  array<Type> eps_at(Y_at.rows(), Y_at.cols()); // matrix of process errors
  eps_at = Y_at - mu_at; // process errors relative to the mean across age and time
  jnLL += GMRF(Q_sparse)(eps_at); 

  // REPORT SECTION ------------------------
  REPORT(jnLL);
  REPORT(Q_sparse);
  // REPORT(I);
  // REPORT(B);
  // REPORT(Omega);
  REPORT( mu_at );

  return(jnLL);
  
}
  
  

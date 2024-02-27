// Purpose: To implement a GMRF as a weight-at-age matrix for EBS pollock
// Creator: Matthew LH. Cheng (UAF-CFOS)
// Date: 1/15/23

#include<TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Function to assemble sparse precision matrix
template<class Type>
// @description: Function that constructs a precision matrix, separable along the
// year, age, and cohort axis. Var_Param allows users to switch between conditional
// variance, and marginal variance.
Eigen::SparseMatrix<Type> construct_Q(int n_years, // Integer of years
                                      int n_ages, // Integer of ages
                                      matrix<Type> ay_Index, // Index matrix to construct
                                      Type rho_y, // Partial correlation by years
                                      Type rho_a, // Partial correlation by ages
                                      Type rho_c, // Partial correlation by cohort
                                      Type log_sigma2, // Variance parameter governing GMRF
                                      int Var_Param // Parameterization of Variance ==0 (Conditional), == 1(Marginal)
                                        ) {
  
  // Dimension to construct matrices
  int total_n = n_years * n_ages; 
  
  // Construct matrices for precision matrix
  Eigen::SparseMatrix<Type> B(total_n,total_n); // B matrix
  Eigen::SparseMatrix<Type> I(total_n,total_n); // Identity matrix
  I.setIdentity(); // Set I to identity matrix
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
  
  // Fill in Omega matrix here (variances)
  if(Var_Param == 0) { // Conditional variance
    
    for(int i = 0; i < total_n; i++) {
      for(int j = 0; j < total_n; j++) {
        if(i == j) Omega.coeffRef(i,j) = 1/exp(log_sigma2);
        else Omega.coeffRef(i,j) = Type(0.0);
      } // j loop
    } // i loop
    
  } // end if conditional variance
  
  if(Var_Param == 1) { // Marginal Variance
    
    // Construct container objects
    matrix<Type> L(total_n, total_n); // L Matrix
    matrix<Type> tmp_I_B = I-B; // Temporary Matrix to store I-B
    L =  tmp_I_B.inverse(); // Invert to get L
    vector<Type> d(total_n); // Store variance calculations
    
    for(int n = 0; n < total_n; n++) {
      if(n == 0) {
        d(n) = exp(log_sigma2); // marginal variance parameter
      } else{
        
        Type cumvar = 0; // Cumulative Variance Container
        
        for(int n1 = 0; n1 < n; n1++) {
          cumvar += L(n, n1) * d(n1) * L(n, n1);
        } // n1 loop
        
        // Calculate diagonal values for omega
        d(n) = (exp(log_sigma2) - cumvar) / pow(L(n, n), 2);
        
      } // else loop
    } // n loop
    
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
  DATA_INTEGER(n_proj_years); // integer for projection years
  int n_years = years.size() + n_proj_years; // integer of number of years
  int n_ages = ages.size(); // integer for number of ages
  int total_n = (n_years + n_proj_years) * n_ages; // integer for nyears * nages
  
  // Define matrices for weight-at-age 
  DATA_MATRIX(X_at); // n_years * n_ages (Observered WAA)
  DATA_MATRIX(Xsd_at); // n_years * n_ages (Observed CV WAA)
  
  // Index matrix to loop through to consruct precision matrix
  DATA_MATRIX(ay_Index);  // (n_years * n_ages) * 2 
  DATA_INTEGER(Var_Param); // Variance parameterization of Precision Matrix
  // == 0 (Conditional), == 1(Marginal)
  
  // PARAMETER SECTION ----------------------------------------
  PARAMETER(rho_a); // Partial correlation by age
  PARAMETER(rho_y); // Partial correlation by year
  PARAMETER(rho_c); // Partial correlation by cohort
  PARAMETER(log_sigma2); // Variance of the GMRF process
  PARAMETER(ln_L0); // vonB length at age 0
  PARAMETER(ln_Linf); // vonB asymptotic length
  PARAMETER(ln_k); // vonB brody growth coefficient
  PARAMETER(ln_alpha); // WL relationship conversion factor
  PARAMETER(ln_beta); // WL relationship allometric scaling
  PARAMETER_ARRAY(ln_Y_at); // Random Effects Weight-at-age array
  
  // Define precision matrix for GMRF
  Eigen::SparseMatrix<Type> Q_sparse(total_n, total_n); // Precision matrix
  
  // Construct precision matrix here
  Q_sparse = construct_Q(n_years, n_ages, ay_Index, 
                         rho_y, rho_a, rho_c, 
                         log_sigma2, Var_Param);
  
  // Define joint negative likelihood here
  Type jnLL = 0;

  // LIKELIHOOD SECTION ------------------------
  array<Type> mu_at(ln_Y_at.rows(), ln_Y_at.cols()); // Mean weight at age across time
  
  // Transform vonB and WL parameters
  Type L0 = exp(ln_L0);
  Type Linf = exp(ln_Linf);
  Type k = exp(ln_k);
  Type alpha = exp(ln_alpha);
  Type beta = exp(ln_beta);
  
  // vonB and allometric weight-length relationship
  for(int a = 0; a < X_at.rows(); a++) {
  for(int t = 0; t < X_at.cols(); t++) {
    mu_at(a,t) = Linf - (Linf-L0)*exp(-k*Type(a)); // LA calculations
    mu_at(a,t) = alpha * pow( mu_at(a,t), beta ); // WL calculations
  } // t loop
} // a loop

  // Evaluate WAA data likelihood
  for(int a = 0; a < X_at.rows(); a++) {
  for(int t = 0; t < X_at.cols(); t++) {
    if( !isNA(X_at(a,t)) ){
      if(Xsd_at(a,t) > 0) jnLL -= dnorm(log(X_at(a,t)), ln_Y_at(a,t), Xsd_at(a,t), true); // only evaluate if sd > 0
    } // if we are not doing projections
  } // t loop
} // a loop
  
  // Evaluate GMRF with precision matrix estimating cohort, year, and age correlations
  array<Type> eps_at(ln_Y_at.rows(), ln_Y_at.cols()); // array of process errors
  eps_at = ln_Y_at - log(mu_at); // process errors relative to the mean across age and time
  jnLL += GMRF(Q_sparse)(eps_at); 

  // REPORT SECTION ------------------------
  REPORT(jnLL);
  REPORT(Q_sparse);
  // REPORT(I);
  // REPORT(B);
  // REPORT(Omega);
  // REPORT(L);
  REPORT(mu_at);
  ADREPORT(ln_Y_at); 

  return(jnLL);
  
} // end objective function
  
  

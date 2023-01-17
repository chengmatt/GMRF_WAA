# purpose: To check dimensions and calcualtions of matrices from TMB
# Creator: Matthew LH. Cheng (UAF-CFOS)


# Check B matrix ----------------------------------------------------------

n_a = 13
n_t = 31
total_n <- n_a * n_t
ay_Index <- as.matrix(expand.grid("age" = seq_len(n_a), 
                                  "year" = seq_len(n_t) ))

rho_a = 0.26
rho_y = 0.835
rho_c = 0.05
margvar = exp(-6.31168892)

B_mat <- matrix(NA, nrow = total_n, ncol = total_n)

for(n in 1:total_n) {
  
  age = ay_Index[n,1];
  year = ay_Index[n,2]; 
  
  if(age > 1) {
    
    for(n1 in 1:total_n) {
      if(ay_Index[n1, 1] == age - 1 & ay_Index[n1, 2] == year) 
        B_mat[n,n1] = rho_y
    } 
  } 
  
  if(year > 1) {
    
    for(n1 in 1:total_n) {
      if(ay_Index[n1, 1] == age & ay_Index[n1, 2] == year - 1) 
        B_mat[n,n1] = rho_a
    } 
  } 
  
  if(year > 1 & age > 1) {
    
    for(n1 in 1:total_n) {
      if(ay_Index[n1, 1] == age - 1 & ay_Index[n1, 2] == year - 1) 
        B_mat[n,n1] = rho_c
    } 
    
  } 
  
} 

# Fill in w/ zeros to check against construct precision fxn
B_mat[is.na(B_mat)] <- 0
I <- matrix(0, total_n, total_n)
diag(I) <- 1

I - t(B_mat)


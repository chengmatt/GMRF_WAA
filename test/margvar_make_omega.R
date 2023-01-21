# To derive omega here and check TMB code
n_a = 13
n_t = 31
total_n = n_a * n_t
log_sigma2 = log(0.1)
d1 <- rep(NA, total_n)
L = solve(I-B)

for(n in 1:total_n) {
  if(n == 1) {
    d1[n] <- exp(log_sigma2)
  } else{
    # Temp index
    tmp_n_indx <- n
    cumvar <- vector()
    for(n1 in 1:tmp_n_indx - 1) {
      cumvar[n1] <- L[n, n1] * d1[n1] * L[n, n1]
    } # n1 loop
    d1[n] = (margvar-sum(cumvar)) / L[n,n]^2
  } # else loop
} # n loop

# Solve Omega recursively for stationary variance
d = rep(NA, nrow(index))
for( n in 1:nrow(index) ){
  if(n==1){
    d[n] = margvar
  }else{
    cumvar = sum(L[n,seq_len(n-1)] * d[seq_len(n-1)] * L[n,seq_len(n-1)])
    d[n] = (margvar-cumvar) / L[n,n]^2
  }
}

# Checking to make sure the 1st loop == 2nd loop 
d1 == d

Omega = sparseMatrix( i=seq_len(n_a*n_t), j=seq_len(n_a*n_t), x=d, dims=rep(n_a*n_t,2) )
Omega_inv = sparseMatrix( i=seq_len(n_a*n_t), j=seq_len(n_a*n_t), x=1/d1, dims=rep(n_a*n_t,2) )

# Eq. 2 from Ver Hoef et al. 2018 "On the relationship between conditional (CAR) and simultaneous (SAR) autoregressive models"
#Q = tcrossprod(I-t(B))
#Q = (I-t(B)) %*% (I-B)
Q = (I-t(B)) %*% Omega_inv %*% (I-B)
diag(solve(Q))
Var = L %*% Omega %*% t(L)

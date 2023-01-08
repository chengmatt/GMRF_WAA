
library(Matrix)
library(mvtnorm)

#'@param n_a Number of ages
#'@param n_t Number of years
#'@param pcorr_age Rho correlation by age
#'@param pcorr_year Rho correlation by year
#'@param pcorr_cohort Rho correlation by cohort

# Precision matrix for a GMRF
make_precision = function( n_a, n_t, pcorr_age, pcorr_year, pcorr_cohort ){

  # Creates a grid of age and year - one age class for every year
  index = expand.grid( "age"=seq_len(n_a), "year"=seq_len(n_t) )

  # I, J, and X are all NULL.
  i = j = x = NULL
  
  for(n in 1:nrow(index)){  # Loop here specifies which age,year combination
    # are "neighbors" - analogous to a spatial weighting matrix/neighborhood
    # structure.
    
    # TESTING
    # n <- 2
    age = index[n,'age'] 
    year = index[n,'year'] 
    
    if( age>1 ){
      i = c( i, n ) 
      j = c( j, which(index[,'age']==(age-1) & index[,'year']==year) )
      x = c(x, pcorr_year) 
      # print(paste(i,j,x))
    }
    
    if( year>1 ){
      i = c( i, n ) # Fill in the nth position
      j = c( j, which(index[,'age']==age & index[,'year']==(year-1)) )
      x = c(x, pcorr_age) 
    }
    
    if( age>1 & year>1 ){
      i = c( i, n ) 
      j = c( j, which(index[,'age']==(age-1) & index[,'year']==(year-1)) )
      x = c(x, pcorr_cohort) 
    }
  } 
  
  # Assemble SAR precision
  B = sparseMatrix( i=i, j=j, x=x, dims=rep(n_a*n_t,2) )  # Fill in sparse matrix with correaltion values
  I = sparseMatrix( i=seq_len(n_a*n_t), j=seq_len(n_a*n_t), x=rep(1,n_a*n_t) ) # Identity matrix here
  
  # Eq. 2 from Ver Hoef et al. 2018 "On the relationship between conditional (CAR) and simultaneous (SAR) autoregressive models"
  Q = tcrossprod(I-t(B)) # Get CAR precision matrix here
  # Q = (I-t(B)) %*% (I-B)
  
  # i = sequential ages i.e., 1,2,3
  # j = repeating years i.e., 1,1,1
  # Each i,j combination = age year combination
  # i,j and j,i with values are "neighbors"

}

# Explore
n_a = 5 # Number of ages
n_t = 20 # Number of years
pcorr_age = 0 # values w/i same age are similar
pcorr_year = 0 # values w/i same year are similar
pcorr_cohort = 1 # values w/i each cohort are similar

Q = make_precision(n_a, n_t, pcorr_age, pcorr_year, pcorr_cohort) # precision matrix
V = solve( Q )  # invert precision matrix to get covariance matrix
Vdense = as.matrix(V)

# Draw from MVN w/ mean vector of 0
Y_at = matrix(rmvnorm(n=1, mean=rep(0,n_a*n_t), sigma=Vdense), 
              nrow=n_a, ncol=n_t)

# Coerce to matrix
Y_mat <- reshape::melt(Y_at) %>% rename(Age = X1, Year = X2) 

# Plot!
ggplot(Y_mat, aes(factor(Year), factor(Age), fill = value )) +
geom_tile() +
scale_fill_gradient2()


# image( y=seq_len(n_a), x=seq_len(n_t), z=t(Y_at), xlab="Year", ylab="Age" )


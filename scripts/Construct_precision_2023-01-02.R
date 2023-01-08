library(here)
library(Matrix)
library(mvtnorm)
library(tidyverse)

# pcorr_year ->  Correlation for ages within a year
# pcorr_age ->  Correlation for years within an age
make_precision <-
function( n_a,
          n_t,
          pcorr_age,
          pcorr_year,
          pcorr_cohort,
          margvar = 1,
          what = "Q" ){

  index = expand.grid( "age"=seq_len(n_a), "year"=seq_len(n_t) )

  i = j = x = NULL
  for( n in 1:nrow(index) ){
    age = index[n,'age']
    year = index[n,'year']
    if( age>1 ){
      i = c( i, n )
      j = c( j, which(index[,'age']==(age-1) & index[,'year']==year) )
      x = c(x, pcorr_year)
    }
    if( year>1 ){
      i = c( i, n )
      j = c( j, which(index[,'age']==age & index[,'year']==(year-1)) )
      x = c(x, pcorr_age)
    }
    if( age>1 & year>1 ){
      i = c( i, n )
      j = c( j, which(index[,'age']==(age-1) & index[,'year']==(year-1)) )
      x = c(x, pcorr_cohort)
    }
  }

  #marg_pcorr_cohort = pcorr_cohort + 2*pcorr_age*pcorr_year
  #condvar = (1 - pcorr_age^2 - pcorr_year^2 - marg_pcorr_cohort^2) * margvar
  condvar = (1 - pcorr_age^2 - pcorr_year^2 - pcorr_cohort^2) * margvar
  #d = NULL
  #for( n in 1:nrow(index) ){
  #  age = index[n,'age']
  #  year = index[n,'year']
  #  if( age==1 & year==1 ){
  #    d = c(d, margvar)
  #  }else if( age==1 & year>1 ){
  #    d = c(d, margvar * (1-pcorr_age^2))
  #  }else if( age>1 & year==1 ){
  #    d = c(d, margvar * (1-pcorr_year^2))
  #  }else{
  #    d = c(d, condvar)
  #  }
  #}

  # Assemble SAR precision
  B = sparseMatrix( i=i, j=j, x=x, dims=rep(n_a*n_t,2) )
  I = sparseMatrix( i=seq_len(n_a*n_t), j=seq_len(n_a*n_t), x=rep(1,n_a*n_t) )
  L = solve(I-B)

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
  if(any(d<0)) stop("Check d")

  #
  Omega = sparseMatrix( i=seq_len(n_a*n_t), j=seq_len(n_a*n_t), x=d, dims=rep(n_a*n_t,2) )
  Omega_inv = sparseMatrix( i=seq_len(n_a*n_t), j=seq_len(n_a*n_t), x=1/d, dims=rep(n_a*n_t,2) )

  # Eq. 2 from Ver Hoef et al. 2018 "On the relationship between conditional (CAR) and simultaneous (SAR) autoregressive models"
  #Q = tcrossprod(I-t(B))
  #Q = (I-t(B)) %*% (I-B)
  Q = (I-t(B)) %*% Omega_inv %*% (I-B)
  Var = L %*% Omega %*% t(L)

  if(what=="Q") return(Q)
  if(what=="Var") return(Var)
  if(what=="Omega") return(Omega)
  if(what=="dmat") return(matrix(d,nrow=n_a,ncol=n_t))
}

# Explore
n_a = 8
n_t = 25
pcorr_age = 0.2
pcorr_year = 0.7
pcorr_cohort = 0.1
#marg_var = condvar / (1 - pcorr_age^2 - pcorr_year^2 - pcorr_cohort^2)
margvar = 0.1
# margvar * (1 - pcorr_age^2 - pcorr_year^2)

Q = make_precision(n_a, n_t, pcorr_age, pcorr_year, pcorr_cohort, margvar)
Omega = make_precision(n_a, n_t, pcorr_age, pcorr_year, pcorr_cohort, margvar, what="Omega")
D = make_precision(n_a, n_t, pcorr_age, pcorr_year, pcorr_cohort, margvar, what="dmat")
V = solve( Q )
Vdense = as.matrix(V)
matrix( diag(Vdense), nrow=n_a, ncol=n_t ) # Check marginal variance

Y_at = matrix( rmvnorm(n=1, mean=rep(0,n_a*n_t), sigma=Vdense), nrow=n_a, ncol=n_t )
image( y=seq_len(n_a), x=seq_len(n_t), z=t(Y_at), xlab="Year", ylab="Age" )

# Check variance
mean(Y_at^2)




# Try difference parameterizations and plot -----------------------------------------

pcorr_age = c(0.1, 0.2, 0.7)
pcorr_year = c(0.2, 0.7, 0.1)
pcorr_cohort = c(0.7, 0.1, 0.2)
margvar = 0.3

# Get label types
type = c("pa = 0.1 py = 0.2 pc = 0.7",
         "pa = 0.2 py = 0.7 pc = 0.1",
         "pa = 0.7 pc = 0.1 pc = 0.2")

# Make labels bquote to plot greek letters
pc_high = paste(bquote("~rho[a] == 0.1~~"), bquote("~rho[y] == 0.2~~"), bquote("~rho[c] == 0.7" )  )
py_high = paste(bquote("~rho[a] == 0.2~~"), bquote("~rho[y] == 0.7~~"), bquote("~rho[c] == 0.1" )  )
pa_high = paste(bquote("~rho[a] == 0.7~~"), bquote("~rho[y] == 0.1~~"), bquote("~rho[c] == 0.2" )  )

# Store values
corr_all <- data.frame()

# Loop through to create Y_at matrix dataframe for different corr values
set.seed(3124)

for(i in 1:length(pcorr_age)) {
  
  # Make precision matrix
  Q = make_precision(n_a, n_t, pcorr_age[i], pcorr_year[i], pcorr_cohort[i], margvar)
  # Get dense covariance matrix
  V = solve( Q )
  Vdense = as.matrix(V)
  
  # Simulate random draws here
  Y_at = matrix( rmvnorm(n=1, mean=rep(0,n_a*n_t), sigma=Vdense), nrow=n_a, ncol=n_t )
  
  # Munge dataframe
  Y_mat <- reshape::melt(Y_at) %>% 
    rename(Age = X1, Year = X2) %>% 
    mutate(type = type[i])
  
  # now bind into one df
  corr_all <- rbind(Y_mat, corr_all)

} # end i loop

# Get greek letters 
corr_all <- corr_all %>% 
  mutate(type = factor(type, levels = c("pa = 0.7 pc = 0.1 pc = 0.2",
                                        "pa = 0.2 py = 0.7 pc = 0.1",
                                        "pa = 0.1 py = 0.2 pc = 0.7"),
                              labels = c(pa_high, py_high, pc_high)))
 
png(here("figs", "fig2_rmvnorm_panel.png"), height = 750, width = 1850)
# Plot!
print(
  ggplot(corr_all, aes(factor(Year), factor(Age), fill = value )) +
    geom_tile(alpha = 1.5) +
    facet_wrap(~type, ncol = 3, labeller = label_parsed) +
    scale_fill_gradient2() +
    scale_y_discrete(breaks = seq(0, 8, 2)) +
    scale_x_discrete(breaks = seq(0, 25, 5)) +
    labs(x = "Year", y = "Age", fill = "Value") +
    theme_bw() +
    theme(legend.position = "top",
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 28),
          legend.key.width = unit(1.5, "cm"),
          panel.spacing = unit(3, "lines"), # facet wrap spacing
          strip.text = element_text(size = 28),
          axis.title = element_text(size = 28),
          axis.text = element_text(size = 25, color = "black")) 
)
dev.off()

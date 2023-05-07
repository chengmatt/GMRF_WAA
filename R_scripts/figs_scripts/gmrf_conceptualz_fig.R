
# setwd(R'(C:\Users\James.Thorson\Desktop\Work files\Collaborations\2023 -- Cheng triple-separable smoother)')

# Visualize graph
library(diagram)
library(here)
n_a = 3
n_t = 4
index = expand.grid( "age"=seq_len(n_a), "year"=seq_len(n_t) )
names <- paste0( rep(1:n_a, n_t), ", ", rep(1:n_t, each=n_a) )
box.type = rep("square",n_a*n_t)

# Construct M
M <- array(0, dim=c(length(names),length(names)), dimnames=list(names,names))
i = j = x = NULL
pcorr_year = 0.2 # "rho_1"
pcorr_age = 0.1 # "rho_2"
pcorr_cohort = 0.7 # "rho_3"
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
M[cbind(i,j)] = x

# Plot
png( here("figs", "fig1_gmrf_conceptual.png"), width=6, height=3, res=200, un="in")
  par( mar=c(0,0,0,0))
  plotmat( M, pos=cbind((index[,2]-0.5)/n_t,(rev(index[,1])-0.5)/n_a), 
           curve = 0, name = names, lwd = 2, arr.pos = 0.55,
           box.lwd = 2, box.type = box.type, box.size=0.05,
           box.prop = 0.5, cex.txt=0.9, arr.length = 0.3)
  mtext( side=1, outer=FALSE, text="Year", line=-2, font=2, cex=1.5)
  mtext( side=2, outer=FALSE, text="Age", line=-2, font=2, cex=1.5)
dev.off()


#library(fda)
#library(lme4)

# Functional ICC computation
# input: 
# test, retest: matrices of dimension (n * T) where n is the number of individuals and T is the number of evaluation points
#               note that the dimensions of the two matrices must coincide
#               and functional data must be evaluated on an equally spaced grid of points
#
# output:
# list of two elements
# ICC: integrated ICC
# pointwise_ICC: vector of length T with the pointwise evaluation of the ICC
functional_ICC = function(test,retest,alpha=0.05){
  t = dim(test)[2]
  n = dim(test)[1]
  N = 2*n
  
  
  nik = matrix(nrow=n,ncol=2,data=1)
  Ni = rowSums(nik)
  Nk = colSums(nik)
  K = length(Nk)
  ks = (N - sum(Ni^2)/N) / (n-1)
  kt = (N - sum(Nk^2)/N) / (K-1)
  
  
  if(dim(retest)[1] != n){
    stop("Number of individuals in test and retest must coincide")
  }
  if(dim(retest)[2] != t){
    stop("Number of evaluation points in test and retest must coincide")
  }
  
  pointwise_ICC =
    #pointwise_ICC.2= 
    low_ICC = up_ICC = numeric(t)
  for(point in 1:t){
    data.all = data.frame(response = c(test[,point],retest[,point]), individual = factor(c(1:n,1:n)), measurement = factor(c(rep(1,n),rep(2,n))))
    anova = aov(response ~ individual + measurement, data = data.all)
    MS = summary(anova)[[1]]$'Mean Sq'
    #ICC.A1.Pointwise.2[point] = (MS[1]-MS[3])/(MS[1]+MS[3]+2/n*(MS[2]-MS[3]))
    MSs = MS[1]
    MSt = MS[2]
    MSe = MS[3]
    pointwise_ICC[point] = (kt*(MSs-MSe)) / (ks*MSt + kt*MSs + (ks*kt-ks-kt)*MSe) 
    ####
    #pointwise_ICC.2[point]=(MS[1]-MS[3])/(MS[1]+MS[3]+2/n*(MS[2]-MS[3]))
    a = ks*pointwise_ICC[point]/(kt*(1-pointwise_ICC[point]))
    b = 1 + ks*pointwise_ICC[point]*(kt-1)/(kt*(1-pointwise_ICC[point]))
    nu = (a*MSt+b*MSe)^2 / ((a*MSt)^2/(K-1) + (b*MSe)^2/(N-n-K+1 ))
    Fstar1 = qf(1-alpha/2,n-1,nu)
    Fstar2 = qf(1-alpha/2,nu,n-1)
    low_ICC[point] = kt*(MSs-Fstar1*MSe)/(kt*MSs + Fstar1*(K*MSt+(kt*ks-kt-ks)*MSe) )
    up_ICC[point] = kt*(Fstar2*MSs-MSe)/((ks*MSt+(kt*ks-kt-ks)*MSe) + Fstar2*kt*MSs)
  }
  
  return(list(ICC=mean(pointwise_ICC), 
              #pointwise_ICC.2=pointwise_ICC.2,
              pointwise_ICC=pointwise_ICC,
              low=low_ICC,up=up_ICC))
}

######################################
boot_function=function(data,number_of_resample){
  matrix_data=matrix(NA,length(data),number_of_resample)
  for (i in 1:number_of_resample) {
    data_sample=sample(data,size = length(data),replace = T)
    matrix_data[,i]=data_sample
  }
  return(matrix_data)
}
#######################################
phi_basis=function(Domain,Break_points,N_order){
  phi=bsplineS(x=Domain,breaks = Break_points,norder=N_order)
  return(phi)
}
#######################################
Real_ICC=function(sigma_B,sigma_E,p,ncoef,domain,breaks_points,n_order,c_i){
  sigma_B_matrix=matrix(sigma_B,nrow = p,ncol = ncoef,byrow = T)
  phi=phi_basis(domain,breaks_points,n_order)
  sigmaB_Y=(sigma_B_matrix)^2*(phi)^2
  sigma_E_matrix=matrix(sigma_E,nrow = p,ncol = ncoef,byrow = T)
  sigmaE_Y=(sigma_E_matrix)^2*(phi)^2
  theta_2=2*(c_i^2)
  Real_ICC=rowSums(sigmaB_Y)/(rowSums(sigmaB_Y)+(theta_2)+rowSums(sigmaE_Y))
  return(Real_ICC)
}
#######################################
#Calculation of area
area_approx = function(curve1,curve2,domain,p){
  f1 <- approxfun(domain, curve1-curve2)     # piecewise linear function
  f2 <- function(x) abs(f1(x)) # take the positive value
  out=integrate(f2, domain[1], domain[p], stop.on.error = FALSE)
  return(out$value)
}
##################################

boot_conf = function(sample,estimate,alpha=0.05){
  n = length(sample)
  ave = mean(sample)
  se = sd(sample)
  z = qnorm(1-alpha/2)
  norm = c(estimate-z*se, estimate+z*se)
  q_l = quantile(sample, alpha/2, type = 4)
  q_u = quantile(sample, 1-alpha/2, type = 4)
  perc = c(q_l, q_u)
  basic = c(2*estimate-q_u, 2*estimate-q_l)
  list(normal=norm,percentile = perc,basic = basic)
}

###############################

BCa <- function(y1,y2,sample,estimate,alpha = 0.05) {
  
  n_ind = dim(y1)[1]
  Int_ICC = numeric(n_ind)
  
  for (j in 1:n_ind) {
    
    Int_ICC[j] = functional_ICC(y1[-j,],y2[-j,])$ICC
  }
  
  #Desired quantiles
  u <- c(alpha/2, 1-alpha/2) 
  
  #Compute constants
  z0 <- qnorm(mean(sample < estimate))
  zu <- qnorm(u)
  
  Int_ICC_i_mean = mean(Int_ICC)
  
  a = sum((Int_ICC-Int_ICC_i_mean)^3)/(6*(sum((Int_ICC-Int_ICC_i_mean)^2))^(3/2))
  #Adjusted quantiles
  u_adjusted <- pnorm(z0 + (z0+zu)/(1-a*(z0+zu))) 
  
  #Accelerated Bootstrap CI
  bca = quantile(sample, u_adjusted,type=4)
  names(bca)=c("BCa Low", "BCa Up")
  list(interval=bca)
}







###################################
# functional_ICC_lmer = function(test,retest,alpha=0.05){
#   t = dim(test)[2]
#   pointwise_ICC_lmer = numeric(t)
#   for(point in 1:t){
#     data.all = data.frame(response = c(test[,point],retest[,point]), individual = factor(c(1:n,1:n)), measurement = factor(c(rep(1,n),rep(2,n))))
#     lmer_model=lmer(response ~ (1 | measurement) + (1 | individual), data=data.all)
#     var_lmer=as.data.frame(VarCorr(lmer_model),comp="Variance")
#     sig_r=var_lmer[1,"vcov"]
#     sig_c=var_lmer[2,"vcov"]
#     sig_e=var_lmer[3,"vcov"]
#     pointwise_ICC_lmer[point]=sig_r/(sig_r+sig_c+sig_e)
#   }
#   return(list(ICC_lmer=pointwise_ICC_lmer))
# }
# 
# pointwise_ICC-pointwise_ICC_lmer

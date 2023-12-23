library(mvtnorm)
library(fda)
library(GET)
library(doParallel)


#set the source path file for functions_ICC.R code like "/pfs/stor10/users/home/s/seydi/functions.R"
source()
#sigma between (row effect) like sigma_B=c(2.5,5,5,5,2.5) 
#C_i = the column effect -> C_i = 0.3
#n = number of individuals (number of curves)
#S= #number of MC simulations
ICC_band <- function(S_B,C_i,n_indiv,S_MC,id_code){
  

sigma_B=S_B
c_i=c(C_i,-C_i) #the column effect
sigma_E=c(1,1,1,1,1) #sigma error
n=n_indiv
ncoef = 5 # number of basis coefficients
p = 100 # time points 
domain = seq(0,100,len=p) # domain
breaks_points=c(0,50,100) # for calculation of basis function in real ICC function
sigma_B_diag = diag(sigma_B) # covariance matrix between basis coefficients
sigma_E_diag = diag(sigma_E) # covariance matrix error term basis coefficients
n_order=4 #order of bsplines
basis = create.bspline.basis(range(domain),nbasis=ncoef,norder=n_order) # b spline basis
mucoef = rep(0,ncoef) # mean vector of basis coefficients 
mucoef[2] = 5
mucoef[ncoef-1] = -5
mue=rep(0,ncoef) #mean of error

number_cores=detectCores()
registerDoParallel(number_cores-1) 

S=S_MC
n_boot=1000 # number of bootstrapping
# out_MC=list() 
# out_original=matrix(NA,nrow = p ,ncol = S)
# up_band=matrix(NA,nrow =p ,ncol = S)
# low_band=matrix(NA,nrow =p ,ncol = S)
# li is a list of parallelization output
li = foreach (k = (1:S),.combine = cbind,.packages = c("fda","mvtnorm"))%dopar% {
    
  
    c_e_i  = rmvnorm(n,mucoef,(sigma_B_diag)^2) #random basis coef for between individuals (row effect)
    c_e_1  = rmvnorm(n,mue,(sigma_E_diag)^2) #random basis coef for error in test
    c_e_2  = rmvnorm(n,mue,(sigma_E_diag)^2) #random basis coef for error in retest

    # functional format
    f_e_i = fd(coef=t(c_e_i),basisobj=basis) 
    f_e1 = fd(coef=t(c_e_1),basisobj=basis)
    f_e2 = fd(coef=t(c_e_2),basisobj=basis)
    
    #discretize curves
    e_i  = t(eval.fd(domain,f_e_i))
    e1   = t(eval.fd(domain,f_e1))
    e2   = t(eval.fd(domain,f_e2))
    
    
    y1 = e_i + c_i[1] + e1 #test
    y2 = e_i + c_i[2] + e2 #retest

    
    out_function_ICC = functional_ICC(y1,y2)
    ICC.Pointwise = out_function_ICC$pointwise_ICC #pointwise ICC
    up = out_function_ICC$up #upper pointwise band based on F
    low = out_function_ICC$low #lower pointwise band based on F
    Integrated_ICC = out_function_ICC$ICC
    ##########
    
    # Resampling of curves
    
    #boot_data_list1=boot_data_list2=list(n_boot)
    #m1 = m2 = matrix(NA,n,p)
    s = boot_function(data=1:n,number_of_resample = n_boot) #boot_function is defined in functions part
    #View(s)
    ICC.Pointwise.boot = matrix(NA,p,n_boot)
    ICC.boot=numeric()
      for (i in 1:n_boot) {
        #    print(paste("iter MC",k,"iter boot",i))
        # for (j in 1:n) {
        #   m1[j,] = y1[s[j,i],]
        #   m2[j,] = y2[s[j,i],]
        # }

           m1 = y1[s[,i],]
           m2 = y2[s[,i],]

      
      #boot_data_list1[[i]]=m1
      #boot_data_list2[[i]]=m2
      m_function_ICC=  functional_ICC(m1,m2)   
      ICC.Pointwise.boot[,i] = m_function_ICC$pointwise_ICC #ICC curves based on bootstrap samples
      ICC.boot[i]=m_function_ICC$ICC
      }
    
    bca_CI=BCa(y1,y2,ICC.boot,Integrated_ICC)$interval
    
    list('out_MC' = ICC.Pointwise.boot
         ,'out_original' = ICC.Pointwise
         ,'up_band' = up
         ,'low_band' = low
         ,'Integrated_ICC' = Integrated_ICC
         ,'ICC_boot_integ' = ICC.boot
         ,'BCa' = bca_CI)
}


#####################################
complete_cover_norm = complete_cover_perc =
  complete_cover_basic = complete_cover_bca = 0
width_norm = width_perc = width_basic = width_bca = numeric(S)

#####################################

#check validity of bands for three differents bands
#1- rank
#2- erl
#3- pointwise F

complete_cover_rank = 0
complete_cover_erl = 0
complete_cover_F = 0
complete_cover_z = 0
area_rank=numeric(S)
area_erl=numeric(S)
area_F=numeric(S)
area_z=numeric(S)
real_ICC=Real_ICC(sigma_B,sigma_E,p,ncoef,domain,breaks_points,n_order,c_i)
for (q in 1:S) {
  
  observations = li["out_MC",paste0("result.",q)][[1]]
  
  l = list(r = domain, obs = observations)
  boot.cset = create_curve_set(l)
  cr_rank = central_region(boot.cset, coverage=0.95, type="rank", central = "mean")
  control_rank = sum(real_ICC >= cr_rank$lo &
                       real_ICC <= cr_rank$hi)
  area_rank[q] = area_approx(cr_rank$hi, cr_rank$lo,domain,p)
  cr_erl = central_region(boot.cset, coverage=0.95, type = "erl", central = "mean")
  control_erl = sum(real_ICC >= cr_erl$lo &
                      real_ICC <= cr_erl$hi)
  area_erl[q] = area_approx(cr_erl$hi, cr_erl$lo,domain,p)
  control_F = sum(real_ICC >= li["low_band",paste0("result.",q)][[1]] &
                    real_ICC <= li["up_band",paste0("result.",q)][[1]])
  area_F[q] = area_approx(li["up_band",paste0("result.",q)][[1]],
                          li["low_band",paste0("result.",q)][[1]],domain,p)
  estim=li["out_original",paste0("result.",q)][[1]]
  #mean=apply(observations, 1, mean)
  #ICC_hat=li["out_original", 
  #paste0("result.",plot_sample)][[1]]
  sd=apply(observations, 1, sd)
  up_z=estim+qnorm(.975)*sd
  low_z=estim+qnorm(.025)*sd
  control_z=sum(real_ICC >= low_z &
                  real_ICC <= up_z)
  area_z[q] = area_approx(up_z, low_z,domain,p)


  ################
  mean_true = mean(real_ICC)
  samp = li["ICC_boot_integ",paste0("result.",q)][[1]]
  est = li["Integrated_ICC",paste0("result.",q)][[1]]
  normal = boot_conf(sample =  samp, estimate = est)$normal
  width_norm[q] = normal[2]-normal[1]
  percentile = boot_conf(sample =  samp, estimate = est)$percentile
  width_perc[q] = percentile[2]-percentile[1]
  basic = boot_conf(sample =  samp, estimate = est)$basic
  width_basic[q] = basic[2]-basic[1]
  bca = li["BCa",paste0("result.",q)][[1]]
  width_bca[q] = bca[2]-bca[1]
  
  if (normal[1] <= mean_true && normal[2] >= mean_true){
    complete_cover_norm = complete_cover_norm + 1
  }
  if (percentile[1] <= mean_true && percentile[2] >= mean_true){
    complete_cover_perc = complete_cover_perc + 1
  }
  if (basic[1] <= mean_true && basic[2] >= mean_true){
    complete_cover_basic = complete_cover_basic + 1
  }
  if (bca[1] <= mean_true && bca[2] >= mean_true){
    complete_cover_bca = complete_cover_bca + 1
  }
  
  ###############
  
  if(control_rank == p){
    complete_cover_rank = complete_cover_rank + 1
  }
  if (control_erl == p) {
    complete_cover_erl = complete_cover_erl + 1
  }
  if (control_F == p) {
    complete_cover_F = complete_cover_F + 1
  }
  if (control_z==p) {
      complete_cover_z=complete_cover_z+1
  }
    
  
}
out <- data.frame(ccr=complete_cover_rank,cce=complete_cover_erl,
                  ccf=complete_cover_F,
                  ccz=complete_cover_z,   ar=mean(area_rank),   ae=mean(area_erl),
                  aF=mean(area_F), az=mean(area_z),
                  Int_norm=complete_cover_norm , Int_perc=complete_cover_perc ,
                    Int_basic=complete_cover_basic , Int_bca=complete_cover_bca,
                  w_n=mean(width_norm) , w_p=mean(width_perc) ,
                  w_basic= mean(width_basic) , w_bca=mean(width_bca)
                  )

write.csv2(out,file = paste0(id_code,".csv"))
}




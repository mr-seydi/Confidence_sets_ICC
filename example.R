
#set source path of main "/pfs/stor10/users/home/s/seydi/main.R"
source()


#Set combination of parameters that you want to see the output

#sigma between (row effect) like sigma_B=c(2.5,5,5,5,2.5) 
#C_i = the column (test-retest) effect like C_i = 0.3
#n = number of individuals (number of curves)
#S= #number of MC simulations
#id_code is arbitrary code name of the output file that you get



#ICC_band <- function(S_B,C_i,n_indiv,S_MC,id_code)

ICC_band(S_B = c(2.5,5,5,5,2.5), C_i = 0.1, n_indiv = 10,
         S_MC = 1000, id_code = 1)

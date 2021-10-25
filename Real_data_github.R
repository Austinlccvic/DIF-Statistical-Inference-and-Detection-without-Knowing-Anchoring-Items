library(fastGHQuad)
library(quantreg)
library(MASS)
library(ltm)

### Read in data and pre-processing ###
female=read.table('females.dat', header=TRUE, skip=0)
male=read.table('males.dat', header = TRUE, skip = 0)
female=as.matrix(female)
male=as.matrix(male)

#First 32 are p items
#Next 23 are e items
#Remaining are n items (24)

p_idx=1:32
e_idx=33:55
n_idx=56:ncol(female)
dat_p=rbind(female[, p_idx], male[, p_idx])
dat_e=rbind(female[, e_idx], male[, e_idx])
dat_n=rbind(female[, n_idx], male[, n_idx])

n_female=dim(female)[1] #823
n_male=dim(male)[1] #609
x.vec=c(rep(0, n_female), rep(1, n_male))



gauss=gaussHermiteData(n=31)  
w=gauss$w
z=gauss$x


###Train using data set p###
N=dim(dat_p)[1]
J=dim(dat_p)[2]

set.seed(2)
#Train model
r_p=EM_2PL(a=runif(J), b=runif(J), g=runif(J), x=x.vec, dat_p, z, mu=0.5, w)
info=information(a=r_p$alpha.vec, b=r_p$beta.vec, g=r_p$gamma.vec, x=x.vec, dat_p, z, mu=r_p$mu, q=r_p$post)
varian=solve(info)
mu_hat=r_p$mu
a_hat=r_p$alpha.vec
b_hat=r_p$beta.vec
g_hat=r_p$gamma.vec
M=50000
sim=mvrnorm(n=M, mu=c(mu_hat, a_hat, b_hat, g_hat[2:J]), Sigma=varian)
sim_a=sim[, 2:(J+1)]
sim_gam=cbind(rep(0, M),sim[, (2*J+2):(3*J)])
rqfit = rq(sim_gam[1, ] ~ -1+a_hat)
gamm_transform=sim_gam[1, ]-rqfit$coefficients[1]*a_hat
for (k in 2:M){
  rqfit <- rq(sim_gam[k, ] ~ -1+a_hat)
  gamm_transform_c=sim_gam[k, ]-rqfit$coefficients[1]*a_hat
  gamm_transform=cbind(gamm_transform, gamm_transform_c)
}
upper_p=c()
lower_p=c()
for (j in 1:J){
  upper_p=c(upper_p, quantile(gamm_transform[j, ], 0.975))
  lower_p=c(lower_p, quantile(gamm_transform[j, ], 0.025))
}

#gamma values after transformation
rqfit <- rq(g_hat ~ -1+a_hat)
transform_p=g_hat-rqfit$coefficients[1]*a_hat
#0.080199111  0.244312915 -0.760644500  0.752013290  0.519911269
#-0.034662248  0.021428379 -0.001644909  0.344869622  0.561173585
#0.000000000  0.367542226  0.618374443 -0.536699971  0.130527843
#0.773998636  0.351661765 -0.248711891 -0.589529133 -0.124219135
#-0.145990681 -0.242751422 -0.006196621  0.004901160 -0.107805552
#-0.210992003  0.179947844 -0.557039139  0.084361490 -0.327670348
#-0.139850302  0.280599547


#item numbers
#p_25 p_29 p_30 p_34 p_37 p_42 p_48 p_50 p_56 p_73 p_75 p_91 p_95 p_2 p_5
#p_7 p_9 p_12 p_14 p_18 p_21 p_41 p_54 p_59 p_64 p_68 p_79 p_81 p_85 p_88
#p_96 p_99






#To find a list of p-values
ga_trans=gamm_transform-transform_p%*%t(rep(1, M))
p_value_p=c()
for (j in 1:J){
  stat=transform_p[j]
  dist=as.vector(ga_trans[j, ])
  p_value_c=length(c(which(dist>abs(stat)), which(dist<   -abs(stat))))/M
  p_value_p=c(p_value_p, p_value_c)
}

order_p=order(p_value_p)
p_value_p[order_p]
#List of p-values
#0.00000 0.00002 0.00004 0.00008 0.00014 0.00050 0.00534 0.00646
#0.00736 0.00750 0.01800 0.05182 0.06756 0.11788 0.16000 0.20586
#0.28952 0.29974 0.30520 0.31474 0.32266 0.32318 0.37490 0.44788
#0.48100 0.53366 0.72276 0.81756 0.90470 0.92302 0.93352 0.94764

#Item order
item_order_p=c(25,29,30,34,37,42,48,50,56,73,75,91,95,2,5,7,9,12,14, 18, 21, 41, 54, 59, 64,68, 79, 81, 85, 88, 96, 99)
item_order_p[order_p]
#7 14 34 95 81  2 73 30 37  9 88 91 29 56 41 99  5 79 12 68 21 18 96 64 85 25 42 48 59 54 50 75




item_reverse=14:length(item_order_p)
transform_p[item_reverse]=-transform_p[item_reverse]
upp_reverse=upper_p[item_reverse]
upper_p[item_reverse]=-lower_p[item_reverse]
lower_p[item_reverse]=-upp_reverse

##Benjamini Hochberg
ordered_p=p_value_p[order_p]
rank_p=1:J
#rank_p=rank(p_value_p, ties.method='random')
critical_value=c()
for (i in 1:J){
  critical_value=c(critical_value, 0.05*rank_p[i]/J)
}

p_v=ordered_p[1]
i=1
while (p_v<critical_value[i]){
  i=i+1
  p_v=ordered_p[i]
}
threshold=ordered_p[i-1]
#0.0075
#So 10 items are selected
#7 14 34 95 81  2 73 30 37  9





####Hard thresholding method
gamm_transform=transform_p
#To remember the original index of each element
gamm_o=list()
for (i in 1:length(gamm_transform)){
  gamm_o[[i]]=c(i, gamm_transform[i])
}
ord=order(abs(gamm_transform))
gamm_ord=list()
for (i in 1:length(gamm_transform)){
  idx=ord[i]
  gamm_ord[[i]]=gamm_o[[idx]]
}

g_temp=NULL
BIC_temp=c()
AIC_temp=c()
for (num0 in c(17:25)){
  gamm_trial=gamm_ord
  for (i in 1:num0){
    gamm_trial[[i]][2]=0
  }
  gam_final_trail=c()
  for (i in 1:length(gamm_transform)){
    gam_final_trail[gamm_trial[[i]][1]]=gamm_trial[[i]][2]
  }
  prop_refit=EM_2PL_non_zero(a=a_hat, b=b_hat, g=gam_final_trail, x=x.vec, dat_p, z, mu=mu_hat+rqfit$coefficients[1], w, tol=0.001)
  mll_trail = prop_refit$MLL
  k=length(which(gam_final_trail!=0))+length(which(c(a_hat, b_hat, mu_hat+rqfit$coefficients[1])!=0))
  AIC_temp=c(AIC_temp, 2*k+2*mll_trail)
  BIC_temp=c(BIC_temp, k*log(N)+2*mll_trail)
  g_temp=rbind(g_temp, gam_final_trail)
}
BIC_temp
#41905.09 41899.91 41895.82 41895.24 41893.24 41900.48 41898.95
#41904.47 41921.99
idx_opt_threshold=which.min(BIC_temp)
BIC_temp[idx_opt_threshold]#41893.24
gam_final_BIC=g_temp[idx_opt_threshold, ]
#11 DIF items
#0.0000000  0.0000000 -0.7606445  0.7520133  0.5199113  0.0000000
#0.0000000  0.0000000  0.0000000  0.5611736  0.0000000  0.3675422
#0.6183744  0.5367000  0.0000000 -0.7739986 -0.3516618  0.0000000
#0.5895291  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
#0.0000000  0.0000000  0.0000000  0.5570391  0.0000000  0.0000000
#0.0000000  0.0000000


id_p=c(3,4,5,10,12,13,14,16,17,19,28)
item_order_p[id_p]
#30 34 37 73 91 95  2  7  9 14 81


AIC_temp
#41483.74 41483.83 41485.01 41489.70 41492.96 41505.47 41509.20
#41519.99 41542.78
idx_opt_threshold=which.min(AIC_temp)
AIC_temp[idx_opt_threshold]#41483.74
gam_final_AIC=g_temp[idx_opt_threshold, ]
#15 DIF items
#0.0000000  0.0000000 -0.7606445  0.7520133  0.5199113  0.0000000
#0.0000000  0.0000000  0.3448696  0.5611736  0.0000000  0.3675422
#0.6183744  0.5367000  0.0000000 -0.7739986 -0.3516618  0.2487119
#0.5895291  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
#0.0000000  0.0000000  0.0000000  0.5570391  0.0000000  0.3276703
#0.0000000 -0.2805995






###Train using data set e###
N=dim(dat_e)[1]
J=dim(dat_e)[2]

ind_female=which(x.vec==0)
ind_male=which(x.vec==1)
dat_e_fem=dat_e[ind_female,]
dat_e_male=dat_e[ind_male,]

N_fem=dim(dat_e_fem)[1]
J_fem=dim(dat_e_fem)[2]
N_ma=dim(dat_e_male)[1]
J_ma=dim(dat_e_male)[2]
res_female <- ltm(dat_e_fem ~ z1)
alpha=as.vector(res_female$coefficients[,2])
beta=as.vector(res_female$coefficients[,1])

#To estimate theta
#quad_female=quad(a=alpha, b=beta, g=rep(0, J_fem), x=rep(0, N_fem), dat=dat_e_fem, z, mu=0, w)
#theta_fem=apply(quad_female*(rep(1, N_fem)%*%t(z)), 1, sum)


#res_male <- ltm(dat_e_male ~ z1)
#alpha_ma=as.vector(res_male$coefficients[,2])
#beta_ma=as.vector(res_male$coefficients[,1])
#gamm=alpha_ma-alpha


#Switch item 1 and item 22
col=dat_e[,(J-1)]
dat_e[,(J-1)]=dat_e[,1]
dat_e[,1]=col



set.seed(3)
#Train model
r_e=EM_2PL(a=alpha, b=beta, g=rep(0, J), x=x.vec, dat_e, z, mu=0.5, w)
info=information(a=r_e$alpha.vec, b=r_e$beta.vec, g=r_e$gamma.vec, x=x.vec, dat_e, z, mu=r_e$mu, q=r_e$post)
varian=solve(info)
mu_hat=r_e$mu
a_hat=r_e$alpha.vec
b_hat=r_e$beta.vec
g_hat=r_e$gamma.vec
M=10000
sim=mvrnorm(n=M, mu=c(mu_hat, a_hat, b_hat, g_hat[2:J]), Sigma=varian)
sim_a=sim[, 2:(J+1)]
sim_gam=cbind(rep(0, M), sim[, (2*J+2):(3*J)])
rqfit <- rq(sim_gam[1, ] ~  -1+a_hat)
gamm_transform=sim_gam[1, ]-rqfit$coefficients[1]*a_hat
for (k in 2:M){
  rqfit <- rq(sim_gam[k, ] ~  -1+a_hat)
  gamm_transform_c=sim_gam[k, ]-rqfit$coefficients[1]*a_hat
  gamm_transform=cbind(gamm_transform, gamm_transform_c)
}
upper_e=c()
lower_e=c()
for (j in 1:J){
  upper_e=c(upper_e, quantile(gamm_transform[j, ], 0.975))
  lower_e=c(lower_e, quantile(gamm_transform[j, ], 0.025))
}


rqfit <- rq(g_hat ~  -1 + a_hat)
transform_e=g_hat-rqfit$coefficients[1]*a_hat

#Item numbers (with order of first and J-1 changed)
#e_33 e_6 e_11 e_16 e_20 e_28 e_36 e_40 e_45 e_51 e_55 e_58 e_61
#e_63 e_67 e_69 e_72 e_78 e_90 e_94 e_24 e_1 e_47




#To return a list of p-values
ga_trans=gamm_transform-transform_e%*%t(rep(1, M))
p_value_e=c()
for (j in 1:J){
  stat=transform_e[j]
  dist=as.vector(ga_trans[j, ])
  p_value_c=length(c(which(dist>abs(stat)), which(dist<   -abs(stat))))/M
  p_value_e=c(p_value_e, p_value_c)
}


order_e=order(p_value_e)
p_value_e[order_e]
#0.00014 0.00108 0.00114 0.00194 0.00262 0.00340 0.00578 0.00734
#0.01206 0.01260 0.02264 0.04164 0.10804 0.12628 0.12706 0.41380
#0.62640 0.64666 0.76164 0.77082 0.83650 0.84300 0.92682

#Item order
#e_1 e_6 e_11 e_16 e_20 e_28 e_36 e_40 e_45 e_51 e_55 e_58 e_61
#e_63 e_67 e_69 e_72 e_78 e_90 e_94 e_24 e_33 e_47
item_order_e=c(33, 6, 11, 16, 20, 28, 36, 40, 45, 51, 55, 58, 61, 63, 67, 69, 72, 78, 90, 94, 24, 1, 47)
item_order_e[order_e]
#63 90 36  6 67 33 61 94 78 51 58 11 28 20 55  1 40 16 24 69 47 45 72





##Benjamini Hochberg
ordered_e=p_value_e[order_e]
rank_e=1:J
#rank_p=rank(p_value_p, ties.method='random')
critical_value=c()
for (i in 1:J){
  critical_value=c(critical_value, 0.05*rank_e[i]/J)
}

p_v=ordered_e[1]
i=1
while (p_v<critical_value[i]){
  i=i+1
  p_v=ordered_e[i]
}
threshold=ordered_e[i-1]
#0.02264
#So 11 items are selected





####Hard thresholding method
gamm_transform=transform_e
#To remember the original index of each element
gamm_o=list()
for (i in 1:length(gamm_transform)){
  gamm_o[[i]]=c(i, gamm_transform[i])
}
ord=order(abs(gamm_transform))
gamm_ord=list()
for (i in 1:length(gamm_transform)){
  idx=ord[i]
  gamm_ord[[i]]=gamm_o[[idx]]
}

g_temp=NULL
BIC_temp=c()
AIC_temp=c()
for (num0 in c(8:16)){
  gamm_trial=gamm_ord
  for (i in 1:num0){
    gamm_trial[[i]][2]=0
  }
  gam_final_trail=c()
  for (i in 1:length(gamm_transform)){
    gam_final_trail[gamm_trial[[i]][1]]=gamm_trial[[i]][2]
  }
  prop_refit=EM_2PL_non_zero(a=a_hat, b=b_hat, g=gam_final_trail, x=x.vec, dat_e, z, mu=mu_hat+rqfit$coefficients[1], w, tol=0.001)
  mll_trail = prop_refit$MLL
  k=length(which(gam_final_trail!=0))+length(which(c(a_hat, b_hat, mu_hat+rqfit$coefficients[1])!=0))
  AIC_temp=c(AIC_temp, 2*k+2*mll_trail)
  BIC_temp=c(BIC_temp, k*log(N)+2*mll_trail)
  g_temp=rbind(g_temp, gam_final_trail)
}
BIC_temp
#41589.46 41585.97 41582.43 41581.96 41588.64 41588.63 41598.85
#41616.72 41624.27
idx_opt_threshold=which.min(BIC_temp)
BIC_temp[idx_opt_threshold]#41581.96
gam_final_BIC=g_temp[idx_opt_threshold, ]
#12 DIF items
#0.6343775 -0.5043473 -0.4019456  0.0000000  0.0000000  0.0000000
#-0.5484168  0.0000000  0.0000000  0.5857987  0.0000000 -0.4706421
# 0.3474469  0.4662340  0.4401900  0.0000000  0.0000000  0.5741121
#-0.5222093  0.4973941  0.0000000  0.0000000  0.0000000

id_E=c(1,2,3,7,10,12,13,14,15,18,19, 20)
item_order_e[id_E]
#33  6 11 36 51 58 61 63 67 78 90 94

AIC_temp
#41262.91 41264.69 41266.42 41271.22 41283.16 41288.42 41303.90
#41327.05 41339.86
idx_opt_threshold=which.min(AIC_temp)
AIC_temp[idx_opt_threshold]#41262.91
gam_final_AIC=g_temp[idx_opt_threshold, ]
#15 DIF items
#0.6343775 -0.5043473 -0.4019456  0.0000000 -0.3153472  0.2747863
#-0.5484168  0.0000000  0.0000000  0.5857987  0.2368193 -0.4706421
# 0.3474469  0.4662340  0.4401900  0.0000000  0.0000000  0.5741121
#-0.5222093  0.4973941  0.0000000  0.0000000  0.0000000


rqfit <- rq(g_hat ~  -1 + a_hat)
transform_e=g_hat-rqfit$coefficients[1]*a_hat
item_reverse=c(1, length(item_order_e)-2,length(item_order_e))
transform_e[item_reverse]=-transform_e[item_reverse]
upp_reverse=upper_e[item_reverse]
upper_e[item_reverse]=-lower_e[item_reverse]
lower_e[item_reverse]=-upp_reverse
#-0.634377547 -0.504347278 -0.401945612  0.065622170 -0.315347247
#0.274786269 -0.548416835 -0.061561324  0.003893889  0.585798731
#0.236819290 -0.470642101  0.347446930  0.466234025  0.440189989
#0.033556994 -0.007647829  0.574112115 -0.522209322 -0.497394124
#0.015580704  0.100324777  0.000000000










###Train using data set n###

N=dim(dat_n)[1]
J=dim(dat_n)[2]

set.seed(1)
#Train model
r_n=EM_2PL(a=runif(J), b=runif(J), g=runif(J), x=x.vec, dat_n, z, mu=0.5, w)
info=information(a=r_n$alpha.vec, b=r_n$beta.vec, g=r_n$gamma.vec, x=x.vec, dat_n, z, mu=r_n$mu, q=r_n$post)
varian=solve(info)
mu_hat=r_n$mu
a_hat=r_n$alpha.vec
b_hat=r_n$beta.vec
g_hat=r_n$gamma.vec
sim=mvrnorm(n=M, mu=c(mu_hat, a_hat, b_hat, g_hat[2:J]), Sigma=varian)
sim_a=sim[, 2:(J+1)]
sim_gam=cbind(rep(0, M),sim[, (2*J+2):(3*J)])
rqfit <- rq(sim_gam[1, ] ~ -1+a_hat)
gamm_transform=sim_gam[1, ]-rqfit$coefficients[1]*a_hat
for (k in 2:M){
  rqfit <- rq(sim_gam[k, ] ~ -1+a_hat)
  gamm_transform_c=sim_gam[k, ]-rqfit$coefficients[1]*a_hat
  gamm_transform=cbind(gamm_transform, gamm_transform_c)
}
upper_n=c()
lower_n=c()
for (j in 1:J){
  upper_n=c(upper_n, quantile(gamm_transform[j, ], 0.975))
  lower_n=c(lower_n, quantile(gamm_transform[j, ], 0.025))
}

rqfit <- rq(g_hat ~ -1+a_hat)
transform_n=g_hat-rqfit$coefficients[1]*a_hat
#0.000000000 -0.571950625 -0.205264584  0.509491300 -0.541899906
#0.305902131  0.044495002 -0.010286251 -0.252174596 -0.114791005
#0.237440575  0.317913256  0.296092993  0.049413108  0.520565597
#-0.431399977  0.004186488 -0.112683847  0.433713698  0.564449913
#-0.525509547 -0.401094283 -0.008144173 -0.131023440


#n_3 n_8 n_13 n_17 n_22

#n_26 n_31 n_35 n_38 n_43 n_46 n_52 n_60 n_65 n_70 n_74 n_76 n_80 n_83 n_84

#n_87 n_92 n_97 n_100





#To return a list of p-values
ga_trans=gamm_transform-transform_n%*%t(rep(1, M))
p_value_n=c()
for (j in 1:J){
  stat=transform_n[j]
  dist=as.vector(ga_trans[j, ])
  p_value_c=length(c(which(dist>abs(stat)), which(dist<   -abs(stat))))/M
  p_value_n=c(p_value_n, p_value_c)
}



order_n=order(p_value_n)

p_value_n[order_n]
#0.00002 0.00010 0.00014 0.00014 0.00022 0.00026 0.00042 0.00376
#0.00450 0.00828 0.02028 0.03514 0.07220 0.11180 0.14794 0.25704
#0.32968 0.33684 0.63190 0.63214 0.77634 0.82176 0.88542 0.93842

#Item order
item_order_n=c(3, 8, 13, 17, 22,26, 31, 35, 38, 43, 46, 52, 60, 65, 70, 74, 76, 80, 83, 84,87, 92, 97, 100)
item_order_n[order_n]
#8  84  70  87  22  17  74  83  92  52  60  26  38  46  13 100  43 80  65  31   3  35  97  76



##Benjamini Hochberg
ordered_n=p_value_n[order_n]
rank_n=1:J
#rank_p=rank(p_value_p, ties.method='random')
critical_value=c()
for (i in 1:J){
  critical_value=c(critical_value, 0.05*rank_n[i]/J)
}

p_v=ordered_n[1]
i=1
while (p_v<critical_value[i]){
  i=i+1
  p_v=ordered_n[i]
}
critical_value
threshold=ordered_n[i-1] #0.02028
#So only 11 item is selected
#8  84  70  87  22  17  74  83  92  52  60





####Hard thresholding method
gamm_transform=transform_n
#To remember the original index of each element
gamm_o=list()
for (i in 1:length(gamm_transform)){
  gamm_o[[i]]=c(i, gamm_transform[i])
}
ord=order(abs(gamm_transform))
gamm_ord=list()
for (i in 1:length(gamm_transform)){
  idx=ord[i]
  gamm_ord[[i]]=gamm_o[[idx]]
}

g_temp=NULL
BIC_temp=c()
AIC_temp=c()
for (num0 in c(8:17)){
  gamm_trial=gamm_ord
  for (i in 1:num0){
    gamm_trial[[i]][2]=0
  }
  gam_final_trail=c()
  for (i in 1:length(gamm_transform)){
    gam_final_trail[gamm_trial[[i]][1]]=gamm_trial[[i]][2]
  }
  prop_refit=EM_2PL_non_zero(a=a_hat, b=b_hat, g=gam_final_trail, x=x.vec, dat_n, z, mu=mu_hat+rqfit$coefficients[1], w, tol=0.001)
  mll_trail = prop_refit$MLL
  k=length(which(gam_final_trail!=0))+length(which(c(a_hat, b_hat, mu_hat+rqfit$coefficients[1])!=0))
  AIC_temp=c(AIC_temp, 2*k+2*mll_trail)
  BIC_temp=c(BIC_temp, k*log(N)+2*mll_trail)
  g_temp=rbind(g_temp, gam_final_trail)
}
BIC_temp
#44170.54 44164.77 44159.64 44157.57 44152.62 44155.63 44155.24
#44159.32 44167.54 44178.27
idx_opt_threshold=which.min(BIC_temp)
BIC_temp[idx_opt_threshold]#44152.62
gam_final_BIC=g_temp[idx_opt_threshold, ]
#12 DIF items
#0.0000000 -0.5719506  0.0000000  0.5094913 -0.5418999  0.3059021
#0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.3179133
#0.2960930  0.0000000  0.5205656 -0.4314000  0.0000000  0.0000000
#0.4337137  0.5644499 -0.5255095 -0.4010943  0.0000000  0.0000000

id_N=c(2, 4, 5, 6, 12,13,15,16,19,20,21,22)
item_order_n[id_N]
#8 17 22 26 52 60 70 74 83 84 87 92


AIC_temp
#43828.20 43827.69 43827.83 43831.03 43831.34 43839.62 43844.50
#43853.85 43867.33 43883.33
idx_opt_threshold=which.min(AIC_temp)
AIC_temp[idx_opt_threshold]#43827.69
gam_final_AIC=g_temp[idx_opt_threshold, ]
#14 DIF items
#0.0000000 -0.5719506 -0.2052646  0.5094913 -0.5418999  0.3059021
#0.0000000  0.0000000 -0.2521746  0.0000000  0.2374406  0.3179133
#0.2960930  0.0000000  0.5205656 -0.4314000  0.0000000  0.0000000
#0.4337137  0.5644499 -0.5255095 -0.4010943  0.0000000  0.0000000









#### Plot of confidence intervals ####

J_p=dim(dat_p)[2]
J_e=dim(dat_e)[2]
J_n=dim(dat_n)[2]
pdf("CI_all_real_data.pdf", width=12, height=4)
par(mfrow=c(1,3), mai=c(0.6,0.6,0.6,0.6))

plotCI(x=1:J_p,y=transform_p[order_p],ui=upper_p[order_p], li=lower_p[order_p],
       xlab="Items",ylab='DIF', main='Item Set P', type = 'p',axes=F, gap=0, col='blue')#, slty=peopleMat[,5],pt.bg=peopleMat[,7], pch=peopleMat[,6],xlab="",ylab='People parameters');
axis(1,  at=1:J_p,labels =item_order_p[order_p], cex.axis=.7, las=2);
axis(2, at=seq(-1,2.5,0.5), labels =T)
abline(a=0, b=0, col='red')

plotCI(x=1:J_e,y=transform_e[order_e],ui=upper_e[order_e], li=lower_e[order_e],
       xlab="Items",ylab='DIF', main='Item Set E', type = 'p',axes=F, gap=0, col='blue')#, slty=peopleMat[,5],pt.bg=peopleMat[,7], pch=peopleMat[,6],xlab="",ylab='People parameters');
axis(1,  at=1:J_e,labels =item_order_e[order_e], cex.axis=.75, las=2);
axis(2, at=seq(-1,2.5,0.5), labels =T)
abline(a=0, b=0, col='red')

plotCI(x=1:J_n,y=transform_n[order_n],ui=upper_n[order_n], li=lower_n[order_n],
       xlab="Items",ylab='DIF', main='Item Set N', type = 'p',axes=F, gap=0, col='blue')#, slty=peopleMat[,5],pt.bg=peopleMat[,7], pch=peopleMat[,6],xlab="",ylab='People parameters');
axis(1,  at=1:J_n,labels =item_order_n[order_n], cex.axis=.75, las=2);
axis(2, at=seq(-1,2.5,0.5), labels =T)
abline(a=0, b=0, col='red')

dev.off()













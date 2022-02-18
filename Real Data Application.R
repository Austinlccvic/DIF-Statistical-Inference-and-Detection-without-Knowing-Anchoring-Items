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
M=10000
sim=mvrnorm(n=M, mu=c(mu_hat, a_hat, b_hat, rep(0, length(g_hat[2:J]))), Sigma=varian)
sim_a=sim[, 2:(J+1)]
sim_gam=cbind(rep(0, M), sim[, (2*J+2):(3*J)])
rqfit = rq(sim_gam[1, ] ~ -1+sim_a[1,])
gamm_transform=sim_gam[1, ]-rqfit$coefficients[1]*sim_a[1,]
for (k in 2:M){
  rqfit <- rq(sim_gam[k, ] ~ -1+sim_a[k,])
  gamm_transform_c=sim_gam[k, ]-rqfit$coefficients[1]*sim_a[k,]
  gamm_transform=cbind(gamm_transform, gamm_transform_c)
}

#gamma values after transformation
rqfit <- rq(g_hat ~ -1+a_hat)
transform_p=g_hat-rqfit$coefficients[1]*a_hat
round(transform_p, 3)
#-0.046  0.256 -0.769  0.757  0.539 -0.028  0.025  0.004  0.358
#0.560  0.003  0.377  0.614 -0.536  0.128  0.780  0.350 -0.251
#-0.583 -0.102 -0.139 -0.241 -0.008  0.000 -0.105 -0.216  0.186
#-0.552  0.085 -0.318 -0.138  0.289


upper_p=c()
lower_p=c()
for (j in 1:J){
  upper_p=c(upper_p, quantile(gamm_transform[j, ], 0.975))
  lower_p=c(lower_p, quantile(gamm_transform[j, ], 0.025))
}

upper_p=upper_p+transform_p
lower_p=lower_p+transform_p




#item numbers
#p_25 p_29 p_30 p_34 p_37 p_42 p_48 p_50 p_56 p_73 p_75 p_91 p_95 p_2 p_5
#p_7 p_9 p_12 p_14 p_18 p_21 p_41 p_54 p_59 p_64 p_68 p_79 p_81 p_85 p_88
#p_96 p_99






#To find a list of p-values
p_value_p=c()
for (j in 1:J){
  stat=transform_p[j]
  dist=as.vector(gamm_transform[j, ])
  p_value_c=length(c(which(dist>abs(stat)), which(dist<   -abs(stat))))/M
  p_value_p=c(p_value_p, p_value_c)
}

order_p=order(p_value_p)
#List of p-values
p_value_p[order_p]
#List of p-values
#  0.0014 0.0015 0.0057 0.0061 0.0104 0.0140 0.0364 0.0619 0.0625 0.0681 0.1235 0.2217 0.2304 0.2442
#  0.3389 0.3780 0.4389 0.4557 0.4567 0.5187 0.5515 0.5529 0.5819 0.5888 0.6080 0.7527 0.8441 0.8787
#  0.9447 0.9528 0.9559 0.9616

#Item order
item_order_p=c(25,29,30,34,37,42,48,50,56,73,75,91,95,2,5,7,9,12,14, 18, 21, 41, 54, 59, 64,68, 79, 81, 85, 88, 96, 99)
item_order_p[order_p]
#14  7 34 81 95  2 30 73  9 37 88 91 29 56 99 41 12 68  5 79 96 21 64 18 85 25 42 48 54 50 75 59




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
  critical_value=c(critical_value, 0.10*rank_p[i]/J)
}

p_v=ordered_p[1]
i=1
while (p_v<critical_value[i]){
  i=i+1
  p_v=ordered_p[i]
}
threshold=ordered_p[i-1]
#0.0015
#So 2 items are selected
#14, 7





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
for (num0 in c(10:25)){
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

idx_opt_threshold=which.min(BIC_temp)
BIC_temp[idx_opt_threshold]
gam_final_BIC=g_temp[idx_opt_threshold, ]

idx_opt_threshold=which.min(AIC_temp)
AIC_temp[idx_opt_threshold]#41483.74
gam_final_AIC=g_temp[idx_opt_threshold, ]







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
sim=mvrnorm(n=M, mu=c(mu_hat, a_hat, b_hat, rep(0, length(g_hat[2:J]))), Sigma=varian)
sim_a=sim[, 2:(J+1)]
sim_gam=cbind(rep(0, M), sim[, (2*J+2):(3*J)])
rqfit = rq(sim_gam[1, ] ~ -1+sim_a[1,])
gamm_transform=sim_gam[1, ]-rqfit$coefficients[1]*sim_a[1,]
for (k in 2:M){
  rqfit <- rq(sim_gam[k, ] ~ -1+sim_a[k,])
  gamm_transform_c=sim_gam[k, ]-rqfit$coefficients[1]*sim_a[k,]
  gamm_transform=cbind(gamm_transform, gamm_transform_c)
}


rqfit <- rq(g_hat ~  -1 + a_hat)
transform_e=g_hat-rqfit$coefficients[1]*a_hat

#Item numbers (with order of first and J-1 changed)
#e_33 e_6 e_11 e_16 e_20 e_28 e_36 e_40 e_45 e_51 e_55 e_58 e_61
#e_63 e_67 e_69 e_72 e_78 e_90 e_94 e_24 e_1 e_47


upper_e=c()
lower_e=c()
for (j in 1:J){
  upper_e=c(upper_e, quantile(gamm_transform[j, ], 0.975))
  lower_e=c(lower_e, quantile(gamm_transform[j, ], 0.025))
}


upper_e=upper_e+transform_e
lower_e=lower_e+transform_e




#To return a list of p-values
p_value_e=c()
for (j in 1:J){
  stat=transform_e[j]
  dist=as.vector(gamm_transform[j, ])
  p_value_c=length(c(which(dist>abs(stat)), which(dist<   -abs(stat))))/M
  p_value_e=c(p_value_e, p_value_c)
}


order_e=order(p_value_e)
p_value_e[order_e]
# 0.0000 0.0004 0.0006 0.0011 0.0013 0.0013 0.0016 0.0019 0.0031 0.0051 0.0199 0.0310 0.0644 0.0958
# 0.1278 0.4073 0.6185 0.6439 0.7819 0.8371 0.9291 0.9364 0.9391

#Item order
#e_1 e_6 e_11 e_16 e_20 e_28 e_36 e_40 e_45 e_51 e_55 e_58 e_61
#e_63 e_67 e_69 e_72 e_78 e_90 e_94 e_24 e_33 e_47
item_order_e=c(33, 6, 11, 16, 20, 28, 36, 40, 45, 51, 55, 58, 61, 63, 67, 69, 72, 78, 90, 94, 24, 1, 47)
item_order_e[order_e]
#63 36 90  6 33 67 51 78 94 61 58 11 28 55 20  1 40 16 69 24 45 47 72





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
threshold_e=ordered_e[i-1]
threshold_e
#0.0199
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
for (num0 in c(8:20)){
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

idx_opt_threshold=which.min(BIC_temp)
BIC_temp[idx_opt_threshold]
gam_final_BIC=g_temp[idx_opt_threshold, ]

idx_opt_threshold=which.min(AIC_temp)
AIC_temp[idx_opt_threshold]
gam_final_AIC=g_temp[idx_opt_threshold, ]



rqfit <- rq(g_hat ~  -1 + a_hat)
transform_e=g_hat-rqfit$coefficients[1]*a_hat
item_reverse=c(1, length(item_order_e)-2,length(item_order_e))
transform_e[item_reverse]=-transform_e[item_reverse]
upp_reverse=upper_e[item_reverse]
upper_e[item_reverse]=-lower_e[item_reverse]
lower_e[item_reverse]=-upp_reverse











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

M=10000
sim=mvrnorm(n=M, mu=c(mu_hat, a_hat, b_hat, rep(0, length(g_hat[2:J]))), Sigma=varian)
sim_a=sim[, 2:(J+1)]
sim_gam=cbind(rep(0, M), sim[, (2*J+2):(3*J)])
rqfit = rq(sim_gam[1, ] ~ -1+sim_a[1,])
gamm_transform=sim_gam[1, ]-rqfit$coefficients[1]*sim_a[1,]
for (k in 2:M){
  rqfit <- rq(sim_gam[k, ] ~ -1+sim_a[k,])
  gamm_transform_c=sim_gam[k, ]-rqfit$coefficients[1]*sim_a[k,]
  gamm_transform=cbind(gamm_transform, gamm_transform_c)
}

rqfit <- rq(g_hat ~ -1+a_hat)
transform_n=g_hat-rqfit$coefficients[1]*a_hat
#-0.523827953 -0.571442467 -0.205188135  0.509468810 -0.541677503
#  0.298886864  0.048818464  0.000000000 -0.248975474 -0.109883946
#  0.249829864  0.320615302  0.298054545  0.046934922  0.522431909
# -0.429414166  0.005705671 -0.106797787  0.444818477  0.568188766
# -0.524766526 -0.405899875 -0.011189998 -0.130705057



upper_n=c()
lower_n=c()
for (j in 1:J){
  upper_n=c(upper_n, quantile(gamm_transform[j, ], 0.975))
  lower_n=c(lower_n, quantile(gamm_transform[j, ], 0.025))
}


upper_n=upper_n+transform_n
lower_n=lower_n+transform_n




#To return a list of p-values
p_value_n=c()
for (j in 1:J){
  stat=transform_n[j]
  dist=as.vector(gamm_transform[j, ])
  p_value_c=length(c(which(dist>abs(stat)), which(dist<   -abs(stat))))/M
  p_value_n=c(p_value_n, p_value_c)
}



order_n=order(p_value_n)

p_value_n[order_n]
# 0.0004 0.0006 0.0007 0.0014 0.0016 0.0026 0.0026 0.0037 0.0130 0.0152 0.0264 0.0487 0.0994 0.1553
# 0.1856 0.2337 0.3365 0.4417 0.4694 0.7116 0.7376 0.9220 0.9531 0.9550

#Item order
item_order_n=c(3, 8, 13, 17, 22,26, 31, 35, 38, 43, 46, 52, 60, 65, 70, 74, 76, 80, 83, 84,87, 92, 97, 100)
item_order_n[order_n]
# 8  22  87  84  70   3  74  17  92  83  52  60  26  38  46  13 100  43  80  31  65  97  76  35




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

threshold_n=ordered_n[i-1] #0.0152
threshold_n
#So only 10 item is selected
#8  22  87  84  70   3  74  17  92  83






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
for (num0 in c(8:20)){
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

idx_opt_threshold=which.min(BIC_temp)
BIC_temp[idx_opt_threshold]#44152.62
gam_final_BIC=g_temp[idx_opt_threshold, ]



idx_opt_threshold=which.min(AIC_temp)
AIC_temp[idx_opt_threshold]#43827.69
gam_final_AIC=g_temp[idx_opt_threshold, ]






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













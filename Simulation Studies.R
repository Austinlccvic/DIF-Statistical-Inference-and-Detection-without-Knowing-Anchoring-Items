library(fastGHQuad)
library(quantreg)
library(MASS)
library(zoo)
source("Functions.R")

sigma <- 0.5
d_range <- 1 # 1 = small d_j, 2 = large d_j
g_abs <- 1 # 1 = small DIF, 2 = large DIF
g_prop <- 1 # 1 = high DIF, 2 = medium DIF 3 = low DIF  
N <- 500 # numeric

if (d_range==1)
	beta.vec = rep(c(0.8, 0.2, -0.4, -1, 1), 5) # small d_j
if (d_range==2)
	beta.vec = rep(c(0.8, -0.4, -1.2, -2, 2), 5) # large d_j

if (g_prop==1)
	gamma.vec=c(rep(0, 10), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 3)) # high DIF
if (g_prop==2)
	gamma.vec=c(rep(0, 15), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2)) # medium DIF
if (g_prop==3)
	gamma.vec=c(rep(0, 20), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 1)) # low DIF

alpha.vec=rep(c(1.3,1.4,1.5,1.7,1.6), 5)
n0=N/2
n1=N/2
x.vec=c(rep(0, n0), rep(1, n1))
x.mx=matrix(x.vec)

### Proposed method
r=coverage_proposed_inference(seed=120, N=N, sig=0.05, gamma.vec=g_abs*gamma.vec, beta.vec=beta.vec,  alpha.vec=alpha.vec, mu=0.5, sigma=sigma, x=x.mx, ite=10,  sig_list=seq(0,1,0.02))
write.table(r$cover, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('coverage_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'.csv'))
write.table(r$alpha, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('alpha_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'.csv'))
write.table(r$beta, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('beta_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'.csv'))
write.table(r$gamma, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('gamma_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'.csv'))
write.table(r$mu, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('mu_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'.csv'))
write.table(r$sigma, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('sigma_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'.csv'))
write.table(r$TPR_gam, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('TPR_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'.csv'))
write.table(r$FPR_gam, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('FPR_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'.csv'))
write.table(r$p_value, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('pvalue_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'.csv'))



### Anchor 1
anchor = c(1)
r=coverage_LRT_method(N=N, sig=0.05, gamma.vec=g_abs*gamma.vec, beta.vec=beta.vec, alpha.vec=alpha.vec, mu=0.5, sigma=sigma, x=x.mx, anchor=anchor, ite=10, rep=rep, sig_list=seq(0,1,0.02))
write.table(r$cover, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor1_coverage_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$alpha, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor1_alpha_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$beta, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor1_beta_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$gamma, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor1_gamma_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$mu, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor1_mu_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$sigma, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor1_sigma_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$TPR, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor1_TPR_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$FPR, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor1_FPR_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$p_value, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor1_pvalue_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))

### Anchor 5
anchor=1:5
r=coverage_LRT_method(N=N, sig=0.05, gamma.vec=g_abs*gamma.vec, beta.vec=beta.vec,  alpha.vec=alpha.vec, mu=0.5, sigma=sigma, x=x.mx, anchor=anchor, ite=10, rep=rep, sig_list=seq(0,1,0.02))
write.table(r$cover, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor5_coverage_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$alpha, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor5_alpha_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$beta, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor5_beta_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$gamma, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor5_gamma_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$mu, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor5_mu_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$sigma, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor5_sigma_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$TPR, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor5_TPR_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$FPR, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor5_FPR_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$p_value, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor5_pvalue_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))

### Anchor 10
anchor=seq(1,10, 1)
r=coverage_LRT_method(N=N, sig=0.05, gamma.vec=g_abs*gamma.vec, beta.vec=beta.vec,  alpha.vec=alpha.vec, mu=0.5, sigma=sigma, x=x.mx, anchor=anchor, ite=10, rep=rep, sig_list=seq(0,1,0.02))
write.table(r$cover, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor10_coverage_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$alpha, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor10_alpha_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$beta, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor10_beta_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$gamma, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor10_gamma_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$mu, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor10_mu_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$sigma, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor10_sigma_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$TPR, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor10_TPR_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$FPR, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor10_FPR_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))
write.table(r$p_value, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('anchor10_pvalue_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop,'_rep_', rep, '.csv'))

### Lasso
gauss=gaussHermiteData(n=35)#We use 35 Gaussian quadrature
w=gauss$w
z=gauss$x

r=compare_method(N=N, penalty=penalty, g=g_abs*gamma.vec, b=beta.vec,  a=alpha.vec, mu=0.5, sigma=sigma, w=w, z=z, k=1, rep=rep, lr=lr, tol=0.005)

write.table(r$gam_L1_opt, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('Lasso_opt_gamma_lr',lr,'_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop, '_rep_', rep, '.csv'))
write.table(r$TPR_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('Lasso_TPR_lr',lr,'_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop, '_rep_', rep, '.csv'))
write.table(r$FPR_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('Lasso_FPR_lr',lr,'_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop, '_rep_', rep, '.csv'))
write.table(r$gam_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('Lasso_gamma_lr',lr,'_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop, '_rep_', rep, '.csv'))
write.table(r$AIC_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('Lasso_AIC_lr',lr,'_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop, '_rep_', rep, '.csv'))
write.table(r$opt_penalty_AIC, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('Lasso_opt_AIC_lr',lr,'_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop, '_rep_', rep, '.csv'))
write.table(r$BIC_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('Lasso_BIC_lr',lr,'_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop, '_rep_', rep, '.csv'))
write.table(r$opt_penalty_BIC, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('Lasso_opt_BIC_lr',lr,'_N',N,'_sigma',sigma,'_d_range',d_range,'_g_abs',g_abs,'_g_prop',g_prop, '_rep_', rep, '.csv'))



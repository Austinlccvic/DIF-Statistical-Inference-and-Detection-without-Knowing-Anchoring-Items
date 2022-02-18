library(fastGHQuad)
library(quantreg)
library(MASS)
library(zoo)

### For the proposed inference method ####
#Large DIF
gamma.vec=c(rep(0, 15), rep(c(-1, 1.3, -0.9, 1.2, 1), 2))
beta.vec=rep(c(0.8, 0.2, -0.4, -1, 1), 5)
alpha.vec=rep(c(1.3,1.4,1.5,1.7,1.6), 5)

N=500
r=coverage_proposed_inference(N=N, sig=0.05, gamma.vec=gamma.vec, beta.vec=beta.vec, alpha.vec=alpha.vec, sig_list=seq(0, 1, 0.02))
write.table(r$cover, sep=",",  col.names=FALSE, row.names=FALSE, file = 'coverage_N500_large.csv')
write.table(r$alpha, sep=",",  col.names=FALSE, row.names=FALSE, file = 'alpha_N500_large.csv')
write.table(r$beta, sep=",",  col.names=FALSE, row.names=FALSE, file = 'beta_N500_large.csv')
write.table(r$gamma, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamma_N500_large.csv')
write.table(r$mu, sep=",",  col.names=FALSE, row.names=FALSE, file = 'mu_N500_large.csv')
write.table(r$TPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPR_N500_large.csv')
write.table(r$FPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPR_N500_large.csv')
write.table(r$p_value, sep=",",  col.names=FALSE, row.names=FALSE, file = 'pvalue_N500_large.csv')


N=1000
r=coverage_proposed_inference(N=N, sig=0.05, gamma.vec=gamma.vec, beta.vec=beta.vec, alpha.vec=alpha.vec, sig_list=seq(0, 1, 0.02))
write.table(r$cover, sep=",",  col.names=FALSE, row.names=FALSE, file = 'coverage_N1000_large.csv')
write.table(r$alpha, sep=",",  col.names=FALSE, row.names=FALSE, file = 'alpha_N1000_large.csv')
write.table(r$beta, sep=",",  col.names=FALSE, row.names=FALSE, file = 'beta_N1000_large.csv')
write.table(r$gamma, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamma_N1000_large.csv')
write.table(r$mu, sep=",",  col.names=FALSE, row.names=FALSE, file = 'mu_N1000_large.csv')
write.table(r$TPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPR_N1000_large.csv')
write.table(r$FPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPR_N1000_large.csv')
write.table(r$p_value, sep=",",  col.names=FALSE, row.names=FALSE, file = 'pvalue_N1000_large.csv')

#Small DIF
gamma.vec=c(rep(0, 15), rep(c(-0.4, 0.4, -0.5, 0.4, 0.4), 2))
beta.vec=rep(c(0.8, 0.2, -0.4, -1, 1), 5)
alpha.vec=rep(c(1.3,1.4,1.5,1.7,1.6), 5)
N=500

r=coverage_proposed_inference(N=N, sig=0.05, gamma.vec=gamma.vec, beta.vec=beta.vec, alpha.vec=alpha.vec, sig_list=seq(0, 1, 0.02))
write.table(r$cover, sep=",",  col.names=FALSE, row.names=FALSE, file = 'coverage_N500_small.csv')
write.table(r$alpha, sep=",",  col.names=FALSE, row.names=FALSE, file = 'alpha_N500_small.csv')
write.table(r$beta, sep=",",  col.names=FALSE, row.names=FALSE, file = 'beta_N500_small.csv')
write.table(r$gamma, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamma_N500_small.csv')
write.table(r$mu, sep=",",  col.names=FALSE, row.names=FALSE, file = 'mu_N500_small.csv')
write.table(r$TPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPR_N500_small.csv')
write.table(r$FPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPR_N500_small.csv')
write.table(r$p_value, sep=",",  col.names=FALSE, row.names=FALSE, file = 'pvalue_N500_small.csv')


N=1000
r=coverage_proposed_inference(N=N, sig=0.05, gamma.vec=gamma.vec, beta.vec=beta.vec, alpha.vec=alpha.vec, sig_list=seq(0, 1, 0.02))
write.table(r$cover, sep=",",  col.names=FALSE, row.names=FALSE, file = 'coverage_N1000_small.csv')
write.table(r$alpha, sep=",",  col.names=FALSE, row.names=FALSE, file = 'alpha_N1000_small.csv')
write.table(r$beta, sep=",",  col.names=FALSE, row.names=FALSE, file = 'beta_N1000_small.csv')
write.table(r$gamma, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamma_N1000_small.csv')
write.table(r$mu, sep=",",  col.names=FALSE, row.names=FALSE, file = 'mu_N1000_small.csv')
write.table(r$TPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPR_N1000_small.csv')
write.table(r$FPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPR_N1000__small.csv')
write.table(r$p_value, sep=",",  col.names=FALSE, row.names=FALSE, file = 'pvalue_N1000_small.csv')

### For LRT method ###

#Small DIF
gamma.vec=c(rep(0, 15), rep(c(-0.4, 0.4, -0.5, 0.4, 0.4), 2))
beta.vec=rep(c(0.8, 0.2, -0.4, -1, 1), 5)
alpha.vec=rep(c(1.3,1.4,1.5,1.7,1.6), 5)
N=500

anchor=1:5
r=coverage_LRT(N=N,sig=0.05, gamma.vec=gamma.vec, beta.vec=beta.vec, alpha.vec=alpha.vec,anchor=anchor,sig_list=seq(0, 1, 0.02))
write.table(r$cover, sep=",",  col.names=FALSE, row.names=FALSE, file = 'coverage_N500_small_anchor5.csv')
write.table(r$alpha, sep=",",  col.names=FALSE, row.names=FALSE, file = 'alpha_N500_small_anchor5.csv')
write.table(r$beta, sep=",",  col.names=FALSE, row.names=FALSE, file = 'beta_N500_small_anchor5.csv')
write.table(r$gamma, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamma_N500_small_anchor5.csv')
write.table(r$mu, sep=",",  col.names=FALSE, row.names=FALSE, file = 'mu_N500_small_anchor5.csv')
write.table(r$TPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPR_N500_small_anchor5.csv')
write.table(r$FPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPR_N500_small_anchor5.csv')
write.table(r$p_value, sep=",",  col.names=FALSE, row.names=FALSE, file = 'pvalue_N500_small_anchor5.csv')



anchor=c(1)
r=coverage_LRT(N=N,sig=0.05, gamma.vec=gamma.vec, beta.vec=beta.vec, alpha.vec=alpha.vec,anchor=anchor,sig_list=seq(0, 1, 0.02))
write.table(r$cover, sep=",",  col.names=FALSE, row.names=FALSE, file = 'coverage_N500_small_anchor1.csv')
write.table(r$alpha, sep=",",  col.names=FALSE, row.names=FALSE, file = 'alpha_N500_small_anchor1.csv')
write.table(r$beta, sep=",",  col.names=FALSE, row.names=FALSE, file = 'beta_N500_small_anchor1.csv')
write.table(r$gamma, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamma_N500_small_anchor1.csv')
write.table(r$mu, sep=",",  col.names=FALSE, row.names=FALSE, file = 'mu_N500_small_anchor1.csv')
write.table(r$TPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPR_N500_small_anchor1.csv')
write.table(r$FPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPR_N500_small_anchor1.csv')
write.table(r$p_value, sep=",",  col.names=FALSE, row.names=FALSE, file = 'pvalue_N500_small_anchor1.csv')


anchor=seq(1,10, 1)
r=coverage_LRT(N=N,sig=0.05, gamma.vec=gamma.vec, beta.vec=beta.vec, alpha.vec=alpha.vec,anchor=anchor,sig_list=seq(0, 1, 0.02))
write.table(r$cover, sep=",",  col.names=FALSE, row.names=FALSE, file = 'coverage_N500_small_anchor10.csv')
write.table(r$alpha, sep=",",  col.names=FALSE, row.names=FALSE, file = 'alpha_N500_small_anchor10.csv')
write.table(r$beta, sep=",",  col.names=FALSE, row.names=FALSE, file = 'beta_N500_small_anchor10.csv')
write.table(r$gamma, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamma_N500_small_anchor10.csv')
write.table(r$mu, sep=",",  col.names=FALSE, row.names=FALSE, file = 'mu_N500_small_anchor10.csv')
write.table(r$TPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPR_N500_small_anchor10.csv')
write.table(r$FPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPR_N500_small_anchor10.csv')
write.table(r$p_value, sep=",",  col.names=FALSE, row.names=FALSE, file = 'pvalue_N500_small_anchor10.csv')


N=1000
anchor=1:5
r=coverage_LRT(N=N,sig=0.05, gamma.vec=gamma.vec, beta.vec=beta.vec, alpha.vec=alpha.vec,anchor=anchor,sig_list=seq(0, 1, 0.02))
write.table(r$cover, sep=",",  col.names=FALSE, row.names=FALSE, file = 'coverage_N1000_small_anchor5.csv')
write.table(r$alpha, sep=",",  col.names=FALSE, row.names=FALSE, file = 'alpha_N1000_small_anchor5.csv')
write.table(r$beta, sep=",",  col.names=FALSE, row.names=FALSE, file = 'beta_N1000_small_anchor5.csv')
write.table(r$gamma, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamma_N1000_small_anchor5.csv')
write.table(r$mu, sep=",",  col.names=FALSE, row.names=FALSE, file = 'mu_N1000_small_anchor5.csv')
write.table(r$TPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPR_N1000_small_anchor5.csv')
write.table(r$FPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPR_N1000_small_anchor5.csv')
write.table(r$p_value, sep=",",  col.names=FALSE, row.names=FALSE, file = 'pvalue_N1000_small_anchor5.csv')



anchor=c(1)
r=coverage_LRT(N=N,sig=0.05, gamma.vec=gamma.vec, beta.vec=beta.vec, alpha.vec=alpha.vec,anchor=anchor,sig_list=seq(0, 1, 0.02))
write.table(r$cover, sep=",",  col.names=FALSE, row.names=FALSE, file = 'coverage_N1000_small_anchor1.csv')
write.table(r$alpha, sep=",",  col.names=FALSE, row.names=FALSE, file = 'alpha_N1000_small_anchor1.csv')
write.table(r$beta, sep=",",  col.names=FALSE, row.names=FALSE, file = 'beta_N1000_small_anchor1.csv')
write.table(r$gamma, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamma_N1000_small_anchor1.csv')
write.table(r$mu, sep=",",  col.names=FALSE, row.names=FALSE, file = 'mu_N1000_small_anchor1.csv')
write.table(r$TPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPR_N1000_small_anchor1.csv')
write.table(r$FPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPR_N1000_small_anchor1.csv')
write.table(r$p_value, sep=",",  col.names=FALSE, row.names=FALSE, file = 'pvalue_N1000_small_anchor1.csv')


anchor=seq(1,10, 1)
r=coverage_LRT(N=N,sig=0.05, gamma.vec=gamma.vec, beta.vec=beta.vec, alpha.vec=alpha.vec,anchor=anchor,sig_list=seq(0, 1, 0.02))
write.table(r$cover, sep=",",  col.names=FALSE, row.names=FALSE, file = 'coverage_N1000_small_anchor10.csv')
write.table(r$alpha, sep=",",  col.names=FALSE, row.names=FALSE, file = 'alpha_N1000_small_anchor10.csv')
write.table(r$beta, sep=",",  col.names=FALSE, row.names=FALSE, file = 'beta_N1000_small_anchor10.csv')
write.table(r$gamma, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamma_N1000_small_anchor10.csv')
write.table(r$mu, sep=",",  col.names=FALSE, row.names=FALSE, file = 'mu_N1000_small_anchor10.csv')
write.table(r$TPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPR_N1000_small_anchor10.csv')
write.table(r$FPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPR_N1000_small_anchor10.csv')
write.table(r$p_value, sep=",",  col.names=FALSE, row.names=FALSE, file = 'pvalue_N1000_small_anchor10.csv')



#Large DIF
gamma.vec=c(rep(0, 15), rep(c(-1, 1.3, -0.9, 1.2, 1), 2))
beta.vec=rep(c(0.8, 0.2, -0.4, -1, 1), 5)
alpha.vec=rep(c(1.3,1.4,1.5,1.7,1.6), 5)
N=500

anchor=1:5
r=coverage_LRT(N=N, sig=0.05, gamma.vec=gamma.vec, beta.vec=beta.vec, alpha.vec=alpha.vec,anchor=anchor,sig_list=seq(0, 1, 0.02))
write.table(r$cover, sep=",",  col.names=FALSE, row.names=FALSE, file = 'coverage_N500_large_anchor5.csv')
write.table(r$alpha, sep=",",  col.names=FALSE, row.names=FALSE, file = 'alpha_N500_large_anchor5.csv')
write.table(r$beta, sep=",",  col.names=FALSE, row.names=FALSE, file = 'beta_N500_large_anchor5.csv')
write.table(r$gamma, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamma_N500_large_anchor5.csv')
write.table(r$mu, sep=",",  col.names=FALSE, row.names=FALSE, file = 'mu_N500_large_anchor5.csv')
write.table(r$TPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPR_N500_large_anchor5.csv')
write.table(r$FPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPR_N500_large_anchor5.csv')
write.table(r$p_value, sep=",",  col.names=FALSE, row.names=FALSE, file = 'pvalue_N500_large_anchor5.csv')

anchor=c(1)
r=coverage_LRT(N=N,sig=0.05, gamma.vec=gamma.vec, beta.vec=beta.vec, alpha.vec=alpha.vec,anchor=anchor,sig_list=seq(0, 1, 0.02))
write.table(r$cover, sep=",",  col.names=FALSE, row.names=FALSE, file = 'coverage_N500_large_anchor1.csv')
write.table(r$alpha, sep=",",  col.names=FALSE, row.names=FALSE, file = 'alpha_N500_large_anchor1.csv')
write.table(r$beta, sep=",",  col.names=FALSE, row.names=FALSE, file = 'beta_N500_large_anchor1.csv')
write.table(r$gamma, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamma_N500_large_anchor1.csv')
write.table(r$mu, sep=",",  col.names=FALSE, row.names=FALSE, file = 'mu_N500_large_anchor1.csv')
write.table(r$TPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPR_N500_large_anchor1.csv')
write.table(r$FPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPR_N500_large_anchor1.csv')
write.table(r$p_value, sep=",",  col.names=FALSE, row.names=FALSE, file = 'pvalue_N500_large_anchor1.csv')

anchor=seq(1,10, 1)
r=coverage_LRT(N=N,sig=0.05, gamma.vec=gamma.vec, beta.vec=beta.vec, alpha.vec=alpha.vec,anchor=anchor,sig_list=seq(0, 1, 0.02))
write.table(r$cover, sep=",",  col.names=FALSE, row.names=FALSE, file = 'coverage_N500_large_anchor10.csv')
write.table(r$alpha, sep=",",  col.names=FALSE, row.names=FALSE, file = 'alpha_N500_large_anchor10.csv')
write.table(r$beta, sep=",",  col.names=FALSE, row.names=FALSE, file = 'beta_N500_large_anchor10.csv')
write.table(r$gamma, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamma_N500_large_anchor10.csv')
write.table(r$mu, sep=",",  col.names=FALSE, row.names=FALSE, file = 'mu_N500_large_anchor10.csv')
write.table(r$TPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPR_N500_large_anchor10.csv')
write.table(r$FPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPR_N500_large_anchor10.csv')
write.table(r$p_value, sep=",",  col.names=FALSE, row.names=FALSE, file = 'pvalue_N500_large_anchor10.csv')


N=1000
anchor=1:5
r=coverage_LRT(N=N, sig=0.05, gamma.vec=gamma.vec, beta.vec=beta.vec, alpha.vec=alpha.vec,anchor=anchor,sig_list=seq(0, 1, 0.02))
write.table(r$cover, sep=",",  col.names=FALSE, row.names=FALSE, file = 'coverage_N1000_large_anchor5.csv')
write.table(r$alpha, sep=",",  col.names=FALSE, row.names=FALSE, file = 'alpha_N1000_large_anchor5.csv')
write.table(r$beta, sep=",",  col.names=FALSE, row.names=FALSE, file = 'beta_N1000_large_anchor5.csv')
write.table(r$gamma, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamma_N1000_large_anchor5.csv')
write.table(r$mu, sep=",",  col.names=FALSE, row.names=FALSE, file = 'mu_N1000_large_anchor5.csv')
write.table(r$TPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPR_N1000_large_anchor5.csv')
write.table(r$FPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPR_N1000_large_anchor5.csv')
write.table(r$p_value, sep=",",  col.names=FALSE, row.names=FALSE, file = 'pvalue_N1000_large_anchor5.csv')



anchor=c(1)
r=coverage_LRT(N=N,sig=0.05, gamma.vec=gamma.vec, beta.vec=beta.vec, alpha.vec=alpha.vec,anchor=anchor,sig_list=seq(0, 1, 0.02))
write.table(r$cover, sep=",",  col.names=FALSE, row.names=FALSE, file = 'coverage_N1000_large_anchor1.csv')
write.table(r$alpha, sep=",",  col.names=FALSE, row.names=FALSE, file = 'alpha_N1000_large_anchor1.csv')
write.table(r$beta, sep=",",  col.names=FALSE, row.names=FALSE, file = 'beta_N1000_large_anchor1.csv')
write.table(r$gamma, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamma_N1000_large_anchor1.csv')
write.table(r$mu, sep=",",  col.names=FALSE, row.names=FALSE, file = 'mu_N1000_large_anchor1.csv')
write.table(r$TPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPR_N1000_large_anchor1.csv')
write.table(r$FPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPR_N1000_large_anchor1.csv')
write.table(r$p_value, sep=",",  col.names=FALSE, row.names=FALSE, file = 'pvalue_N1000_large_anchor1.csv')


anchor=seq(1,10, 1)
r=coverage_LRT(N=N,sig=0.05, gamma.vec=gamma.vec, beta.vec=beta.vec, alpha.vec=alpha.vec,anchor=anchor,sig_list=seq(0, 1, 0.02))
write.table(r$cover, sep=",",  col.names=FALSE, row.names=FALSE, file = 'coverage_N1000_large_anchor10.csv')
write.table(r$alpha, sep=",",  col.names=FALSE, row.names=FALSE, file = 'alpha_N1000_large_anchor10.csv')
write.table(r$beta, sep=",",  col.names=FALSE, row.names=FALSE, file = 'beta_N1000_large_anchor10.csv')
write.table(r$gamma, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamma_N1000_large_anchor10.csv')
write.table(r$mu, sep=",",  col.names=FALSE, row.names=FALSE, file = 'mu_N1000_large_anchor10.csv')
write.table(r$TPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPR_N1000_large_anchor10.csv')
write.table(r$FPR, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPR_N1000_large_anchor10.csv')
write.table(r$p_value, sep=",",  col.names=FALSE, row.names=FALSE, file = 'pvalue_N1000_large_anchor10.csv')






#### Compare the hard-thresholding and the LASSO method ####
#Large DIF
gamma.vec=c(rep(0, 15), rep(c(-1, 1.3, -0.9, 1.2, 1), 2))
beta.vec=rep(c(0.8, 0.2, -0.4, -1, 1), 5)
alpha.vec=rep(c(1.3,1.4,1.5,1.7,1.6), 5)
N=500
seed=1

penalty=seq(0.01, 0.15, (0.15-0.01)/19)
threshold=seq(0, 1.8, 1.8/19)
r=compare_method(N=N, penalty=penalty, thresholds=threshold, a=alpha.vec, b=beta.vec, g=gamma.vec, z, mu=0.5, w, k=1, seed, tol = 0.001)
write.table(r$TPR_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPRN500_L1_large.csv')
write.table(r$FPR_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPRN500_L1_large.csv')  
write.table(r$TPR_prop, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPRN500_hard_large.csv')
write.table(r$FPR_prop, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPRN500_hard_large.csv')  

write.table(r$AIC_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = 'AICN500_L1_large.csv')
write.table(r$BIC_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = 'BICN500_L1_large.csv')  
write.table(r$gam_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamallN500_L1_large.csv')
write.table(r$gam_L1_opt, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamoptN500_L1_large.csv')
write.table(r$opt_penalty_BIC, sep=",",  col.names=FALSE, row.names=FALSE, file = 'BICpenaltyoptN500_L1_large.csv')
write.table(r$AIC_prop, sep=",",  col.names=FALSE, row.names=FALSE, file = 'AICN500_hard_large.csv')
write.table(r$BIC_prop, sep=",",  col.names=FALSE, row.names=FALSE, file = 'BICN500_hard_large.csv')  
write.table(r$gam_prop, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamallN500_hard_large.csv')
write.table(r$gam_prop_opt_BIC, sep=",",  col.names=FALSE, row.names=FALSE, file = 'BICgamoptN500_hard_large.csv')  
write.table(r$opt_threshold, sep=",",  col.names=FALSE, row.names=FALSE, file = 'thresholdoptN500_hard_large.csv')


N=1000
penalty=seq(0.01, 0.15, (0.15-0.01)/19)
threshold=seq(0, 1.8, 1.8/19)
seed=1
r=compare_method(N=N, penalty=penalty, thresholds=threshold, a=alpha.vec, b=beta.vec, g=gamma.vec, z, mu=0.5, w, k=1, seed, tol = 0.001)
write.table(r$TPR_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPRN1000_L1_large.csv')
write.table(r$FPR_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPRN1000_L1_large.csv')
write.table(r$TPR_prop, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPRN1000_hard_large.csv')
write.table(r$FPR_prop, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPRN1000_hard_large.csv')

write.table(r$AIC_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = 'AICN1000_L1_large.csv')
write.table(r$BIC_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = 'BICN1000_L1_large.csv')
write.table(r$gam_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamallN1000_L1_large.csv')
write.table(r$gam_L1_opt, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamoptN1000_L1_large.csv')
write.table(r$opt_penalty_BIC, sep=",",  col.names=FALSE, row.names=FALSE, file = 'BICpenaltyoptN1000_L1_large.csv')
write.table(r$AIC_prop, sep=",",  col.names=FALSE, row.names=FALSE, file = 'AICN1000_hard_large.csv')
write.table(r$BIC_prop, sep=",",  col.names=FALSE, row.names=FALSE, file = 'BICN1000_hard_large.csv')
write.table(r$gam_prop, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamallN1000_hard_large.csv')
write.table(r$gam_prop_opt_BIC, sep=",",  col.names=FALSE, row.names=FALSE, file = 'BICgamoptN1000_hard_large.csv')
write.table(r$opt_threshold, sep=",",  col.names=FALSE, row.names=FALSE, file = 'thresholdoptN1000_hard_large.csv')



## Small DIF ##
N=500
gamma.vec=c(rep(0, 15), rep(c(-0.4, 0.4, -0.5, 0.4, 0.4), 2))
beta.vec=rep(c(0.8, 0.2, -0.4, -1, 1), 5)
alpha.vec=rep(c(1.3,1.4,1.5,1.7,1.6), 5)
penalty=seq(0.01, 0.09, (0.09-0.01)/19)
threshold=seq(0, 0.9, 0.9/19)
seed=1
r=compare_method(N=N, penalty=penalty, thresholds=threshold, a=alpha.vec, b=beta.vec, g=gamma.vec, z, mu=0.5, w, k=1, seed, tol = 0.001)
write.table(r$TPR_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPRN500_L1_small.csv')
write.table(r$FPR_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPRN500_L1_small.csv')  
write.table(r$TPR_prop, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPRN500_hard_small.csv')
write.table(r$FPR_prop, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPRN500_hard_small.csv')  

write.table(r$AIC_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = 'AICN500_L1_small.csv')
write.table(r$BIC_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = 'BICN500_L1_small.csv')  
write.table(r$gam_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamallN500_L1_small.csv')
write.table(r$gam_L1_opt, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamoptN500_L1_small.csv')
write.table(r$opt_penalty_BIC, sep=",",  col.names=FALSE, row.names=FALSE, file = 'BICpenaltyoptN500_L1_small.csv')
write.table(r$AIC_prop, sep=",",  col.names=FALSE, row.names=FALSE, file = 'AICN500_hard_small.csv')
write.table(r$BIC_prop, sep=",",  col.names=FALSE, row.names=FALSE, file = 'BICN500_hard_small.csv')  
write.table(r$gam_prop, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamallN500_hard_small.csv')
write.table(r$gam_prop_opt_BIC, sep=",",  col.names=FALSE, row.names=FALSE, file = 'BICgamoptN500_hard_small.csv')  
write.table(r$opt_threshold, sep=",",  col.names=FALSE, row.names=FALSE, file = 'thresholdoptN500_hard_small.csv')





N=1000
penalty=seq(0.01, 0.09, (0.09-0.01)/19)
threshold=seq(0, 0.9, 0.9/19)
seed=1
r=compare_method(N=N, penalty=penalty, thresholds=threshold, a=alpha.vec, b=beta.vec, g=gamma.vec, z, mu=0.5, w, k=1, seed, tol = 0.001)
write.table(r$TPR_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPRN1000_L1_small.csv')
write.table(r$FPR_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPRN1000_L1_small.csv')
write.table(r$TPR_prop, sep=",",  col.names=FALSE, row.names=FALSE, file = 'TPRN1000_hard_small.csv')
write.table(r$FPR_prop, sep=",",  col.names=FALSE, row.names=FALSE, file = 'FPRN1000_hard_small.csv')

write.table(r$AIC_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = 'AICN1000_L1_small.csv')
write.table(r$BIC_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = 'BICN1000_L1_small.csv')
write.table(r$gam_L1, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamallN1000_L1_small.csv')
write.table(r$gam_L1_opt, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamoptN1000_L1_small.csv')
write.table(r$opt_penalty_BIC, sep=",",  col.names=FALSE, row.names=FALSE, file = 'BICpenaltyoptN1000_L1_small.csv')
write.table(r$AIC_prop, sep=",",  col.names=FALSE, row.names=FALSE, file = 'AICN1000_hard_small.csv')
write.table(r$BIC_prop, sep=",",  col.names=FALSE, row.names=FALSE, file = 'BICN1000_hard_small.csv')
write.table(r$gam_prop, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamallN1000_hard_small.csv')
write.table(r$gam_prop_opt_BIC, sep=",",  col.names=FALSE, row.names=FALSE, file = 'BICgamoptN1000_hard_small.csv')
write.table(r$opt_threshold, sep=",",  col.names=FALSE, row.names=FALSE, file = 'thresholdoptN1000_hard_small.csv')








##### Compare the FDR and TPR of the proposed inference method and the LRT method #####
##Proposed method
p_N500_4=read.csv(file = 'pvalue_N500_small.csv', header = F)
p_N500_1=read.csv(file = 'pvalue_N500_large.csv', header = F)
p_N1000_4=read.csv(file = 'pvalue_N1000_small.csv', header = F)
p_N1000_1=read.csv(file = 'pvalue_N1000_large.csv', header = F)


fdr_500_4=FDR(pmat=p_N500_4, gamma.vec=c(rep(0, 15), rep(c(-0.4, 0.4, -0.5, 0.4, 0.4), 2)), fdr=0.05)
fdr_500_1=FDR(pmat=p_N500_1, gamma.vec=c(rep(0, 15), rep(c(-1, 1.3, -0.9, 1.2, 1), 2)), fdr=0.05)
fdr_1000_4=FDR(pmat=p_N1000_4, gamma.vec=c(rep(0, 15), rep(c(-0.4, 0.4, -0.5, 0.4, 0.4), 2)), fdr=0.05)
fdr_1000_1=FDR(pmat=p_N1000_1, gamma.vec=c(rep(0, 15), rep(c(-1, 1.3, -0.9, 1.2, 1), 2)), fdr=0.05)

fdr_proposed=c(fdr_500_4,fdr_500_1,fdr_1000_4,fdr_1000_1)



##LRT method
#anchor set 1
p_N500_4=read.csv(file = 'pvalue_N500_small_anchor1.csv', header = F)
p_N500_1=read.csv(file = 'pvalue_N500_large_anchor1.csv', header = F)
p_N1000_4=read.csv(file = 'pvalue_N1000_small_anchor1.csv', header = F)
p_N1000_1=read.csv(file = 'pvalue_N1000_large_anchor1.csv', header = F)
fdr_500_4=FDR(pmat=p_N500_4, gamma.vec=c(rep(0, 14), rep(c(-0.4, 0.4, -0.5, 0.4, 0.4), 2)), fdr=0.05)
fdr_500_1=FDR(pmat=p_N500_1, gamma.vec=c(rep(0, 14), rep(c(-1, 1.3, -0.9, 1.2, 1), 2)), fdr=0.05)
fdr_1000_4=FDR(pmat=p_N1000_4, gamma.vec=c(rep(0, 14), rep(c(-0.4, 0.4, -0.5, 0.4, 0.4), 2)), fdr=0.05)
fdr_1000_1=FDR(pmat=p_N1000_1, gamma.vec=c(rep(0, 14), rep(c(-1, 1.3, -0.9, 1.2, 1), 2)), fdr=0.05)

fdr_anchor1=c(fdr_500_4,fdr_500_1,fdr_1000_4,fdr_1000_1)


#anchor set 5
p_N500_4=read.csv(file = 'pvalue_N500_small_anchor5.csv', header = F)
p_N500_1=read.csv(file = 'pvalue_N500_large_anchor5.csv', header = F)
p_N1000_4=read.csv(file = 'pvalue_N1000_small_anchor5.csv', header = F)
p_N1000_1=read.csv(file = 'pvalue_N1000_large_anchor5.csv', header = F)
fdr_500_4=FDR(pmat=p_N500_4, gamma.vec=c(rep(0, 10), rep(c(-0.4, 0.4, -0.5, 0.4, 0.4), 2)), fdr=0.05)
fdr_500_1=FDR(pmat=p_N500_1, gamma.vec=c(rep(0, 10), rep(c(-1, 1.3, -0.9, 1.2, 1), 2)), fdr=0.05)
fdr_1000_4=FDR(pmat=p_N1000_4, gamma.vec=c(rep(0, 10), rep(c(-0.4, 0.4, -0.5, 0.4, 0.4), 2)), fdr=0.05)
fdr_1000_1=FDR(pmat=p_N1000_1, gamma.vec=c(rep(0, 10), rep(c(-1, 1.3, -0.9, 1.2, 1), 2)), fdr=0.05)

fdr_anchor5=c(fdr_500_4,fdr_500_1,fdr_1000_4,fdr_1000_1)



#anchor set 10
p_N500_4=read.csv(file = 'pvalue_N500_small_anchor10.csv', header = F)
p_N500_1=read.csv(file = 'pvalue_N500_large_anchor10.csv', header = F)
p_N1000_4=read.csv(file = 'pvalue_N1000_small_anchor10.csv', header = F)
p_N1000_1=read.csv(file = 'pvalue_N1000_large_anchor10.csv', header = F)
fdr_500_4=FDR(pmat=p_N500_4, gamma.vec=c(rep(0, 5), rep(c(-0.4, 0.4, -0.5, 0.4, 0.4), 2)), fdr=0.05)
fdr_500_1=FDR(pmat=p_N500_1, gamma.vec=c(rep(0, 5), rep(c(-1, 1.3, -0.9, 1.2, 1), 2)), fdr=0.05)
fdr_1000_4=FDR(pmat=p_N1000_4, gamma.vec=c(rep(0, 5), rep(c(-0.4, 0.4, -0.5, 0.4, 0.4), 2)), fdr=0.05)
fdr_1000_1=FDR(pmat=p_N1000_1, gamma.vec=c(rep(0, 5), rep(c(-1, 1.3, -0.9, 1.2, 1), 2)), fdr=0.05)

fdr_anchor10=c(fdr_500_4,fdr_500_1,fdr_1000_4,fdr_1000_1)





##### Compare area under the ROC curves of the proposed inference method and the LRT method #####
#proposed methed
TPR_N500_4=read.csv(file = 'TPR_N500_small.csv', header = F)
FPR_N500_4=read.csv(file = 'FPR_N500_small.csv', header = F)
TPR_N500_1=read.csv(file = 'TPR_N500_large.csv', header = F)
FPR_N500_1=read.csv(file = 'FPR_N500_large.csv', header = F)
TPR_N1000_4=read.csv(file = 'TPR_N1000_small.csv', header = F)
FPR_N1000_4=read.csv(file = 'FPR_N1000__small.csv', header = F)
TPR_N1000_1=read.csv(file = 'TPR_N1000_large.csv', header = F)
FPR_N1000_1=read.csv(file = 'FPR_N1000_large.csv', header = F)
prop_N500_4=eva_auc(tpr=TPR_N500_4, fpr=FPR_N500_4)
prop_N500_1=eva_auc(tpr=TPR_N500_1, fpr=FPR_N500_1)
prop_N1000_4=eva_auc(tpr=TPR_N1000_4, fpr=FPR_N1000_4)
prop_N1000_1=eva_auc(tpr=TPR_N1000_1, fpr=FPR_N1000_1)
auc_prop=c(prop_N500_4,prop_N500_1,prop_N1000_4,prop_N1000_1)
 

#LRT
#anchor1
TPR_N500_4=read.csv(file = 'TPR_N500_small_anchor1.csv', header = F)
FPR_N500_4=read.csv(file = 'FPR_N500_small_anchor1.csv', header = F)
TPR_N500_1=read.csv(file = 'TPR_N500_large_anchor1.csv', header = F)
FPR_N500_1=read.csv(file = 'FPR_N500_large_anchor1.csv', header = F)
TPR_N1000_4=read.csv(file = 'TPR_N1000_small_anchor1.csv', header = F)
FPR_N1000_4=read.csv(file = 'FPR_N1000_small_anchor1.csv', header = F)
TPR_N1000_1=read.csv(file = 'TPR_N1000_large_anchor1.csv', header = F)
FPR_N1000_1=read.csv(file = 'FPR_N1000_large_anchor1.csv', header = F)
lrt_N500_4=eva_auc(tpr=TPR_N500_4, fpr=FPR_N500_4)
lrt_N500_1=eva_auc(tpr=TPR_N500_1, fpr=FPR_N500_1)
lrt_N1000_4=eva_auc(tpr=TPR_N1000_4, fpr=FPR_N1000_4)
lrt_N1000_1=eva_auc(tpr=TPR_N1000_1, fpr=FPR_N1000_1)
auc_lrt1=c(lrt_N500_4, lrt_N500_1, lrt_N1000_4, lrt_N1000_1)



#anchor5
TPR_N500_4=read.csv(file = 'TPR_N500_small_anchor5.csv', header = F)
FPR_N500_4=read.csv(file = 'FPR_N500_small_anchor5.csv', header = F)
TPR_N500_1=read.csv(file = 'TPR_N500_large_anchor5.csv', header = F)
FPR_N500_1=read.csv(file = 'FPR_N500_large_anchor5.csv', header = F)
TPR_N1000_4=read.csv(file = 'TPR_N1000_small_anchor5.csv', header = F)
FPR_N1000_4=read.csv(file = 'FPR_N1000_small_anchor5.csv', header = F)
TPR_N1000_1=read.csv(file = 'TPR_N1000_large_anchor5.csv', header = F)
FPR_N1000_1=read.csv(file = 'FPR_N1000_large_anchor5.csv', header = F)
lrt_N500_4=eva_auc(tpr=TPR_N500_4, fpr=FPR_N500_4)
lrt_N500_1=eva_auc(tpr=TPR_N500_1, fpr=FPR_N500_1)
lrt_N1000_4=eva_auc(tpr=TPR_N1000_4, fpr=FPR_N1000_4)
lrt_N1000_1=eva_auc(tpr=TPR_N1000_1, fpr=FPR_N1000_1)
auc_lrt5=c(lrt_N500_4, lrt_N500_1, lrt_N1000_4, lrt_N1000_1)


#anchor10
TPR_N500_4=read.csv(file = 'TPR_N500_small_anchor10.csv', header = F)
FPR_N500_4=read.csv(file = 'FPR_N500_small_anchor10.csv', header = F)
TPR_N500_1=read.csv(file = 'TPR_N500_large_anchor10.csv', header = F)
FPR_N500_1=read.csv(file = 'FPR_N500_large_anchor10.csv', header = F)
TPR_N1000_4=read.csv(file = 'TPR_N1000_small_anchor10.csv', header = F)
FPR_N1000_4=read.csv(file = 'FPR_N1000_small_anchor10.csv', header = F)
TPR_N1000_1=read.csv(file = 'TPR_N1000_large_anchor10.csv', header = F)
FPR_N1000_1=read.csv(file = 'FPR_N1000_large_anchor10.csv', header = F)
lrt_N500_4=eva_auc(tpr=TPR_N500_4, fpr=FPR_N500_4)
lrt_N500_1=eva_auc(tpr=TPR_N500_1, fpr=FPR_N500_1)
lrt_N1000_4=eva_auc(tpr=TPR_N1000_4, fpr=FPR_N1000_4)
lrt_N1000_1=eva_auc(tpr=TPR_N1000_1, fpr=FPR_N1000_1)
auc_lrt10=c(lrt_N500_4, lrt_N500_1, lrt_N1000_4, lrt_N1000_1)



### Explore the coverage rates of the confidence intervals constructed using the proposed inference method and the LRT method ###
#Anchor method
coverage_anchorN500_4=read.csv(file = 'coverage_N500_small_anchor5.csv', header = F)
dim(coverage_anchorN500_4)
gamma.vec=c(rep(0, 15), rep(c(-0.4, 0.4, -0.5, 0.4, 0.4), 2))
ind1=which(gamma.vec!=0)-5
ind0=1:10

cov_rate_N500_anchor4=apply(coverage_anchorN500_4, 2, mean)
cov_anchorN500_4_pos=cov_rate_N500_anchor4[ind1]
cov_anchorN500_4_neg=cov_rate_N500_anchor4[ind0]


coverage_anchorN500_1=read.csv(file = 'coverage_N500_large_anchor5.csv', header = F)
dim(coverage_anchorN500_1)
gamma.vec=c(rep(0, 15), rep(c(-1, 1.3, -0.9, 1.2, 1), 2))
ind1=which(gamma.vec!=0)-5
ind0=1:10

cov_rate_anchorN500_1=apply(coverage_anchorN500_1, 2, mean)
cov_anchorN500_1_pos=cov_rate_anchorN500_1[ind1]
cov_anchorN500_1_neg=cov_rate_anchorN500_1[ind0]


coverage_anchorN1000_4=read.csv(file = 'coverage_N1000_small_anchor5.csv', header = F)
dim(coverage_anchorN1000_4)
gamma.vec=c(rep(0, 15), rep(c(-0.4, 0.4, -0.5, 0.4, 0.4), 2))
ind1=which(gamma.vec!=0)-5
ind0=1:10

cov_rate_N1000_anchor4=apply(coverage_anchorN1000_4, 2, mean)
cov_anchorN1000_4_pos=cov_rate_N1000_anchor4[ind1]
cov_anchorN1000_4_neg=cov_rate_N1000_anchor4[ind0]


coverage_anchorN1000_1=read.csv(file = 'coverage_N1000_large_anchor5.csv', header = F)
dim(coverage_anchorN1000_1)
gamma.vec=c(rep(0, 15), rep(c(-1, 1.3, -0.9, 1.2, 1), 2))
ind1=which(gamma.vec!=0)-5
ind0=1:10

cov_rate_anchorN1000_1=apply(coverage_anchorN1000_1, 2, mean)
cov_anchorN1000_1_pos=cov_rate_anchorN1000_1[ind1]
cov_anchorN1000_1_neg=cov_rate_anchorN1000_1[ind0]


#proposed method
coverage_N500_4=read.csv(file = 'coverage_N500_small.csv', header = F)
dim(coverage_N500_4)
gamma.vec=c(rep(0, 15), rep(c(-0.4, 0.4, -0.5, 0.4, 0.4), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N500_4=apply(coverage_N500_4, 2, mean)
cov_4_N500_pos=cov_rate_N500_4[ind1]
cov_4_N500_neg=cov_rate_N500_4[ind0]


coverage_N500_1=read.csv(file = 'coverage_N500_large.csv', header = F)
dim(coverage_N500_1)
gamma.vec=c(rep(0, 15), rep(c(-1, 1.3, -0.9, 1.2, 1), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N500_1=apply(coverage_N500_1, 2, mean)
cov_1_N500_pos=cov_rate_N500_1[ind1]
cov_1_N500_neg=cov_rate_N500_1[ind0]

coverage_N1000_4=read.csv(file = 'coverage_N1000_small.csv', header = F)
dim(coverage_N1000_4)
gamma.vec=c(rep(0, 15), rep(c(-0.4, 0.4, -0.5, 0.4, 0.4), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N1000_4=apply(coverage_N1000_4, 2, mean)
cov_4_N1000_pos=cov_rate_N1000_4[ind1]
cov_4_N1000_neg=cov_rate_N1000_4[ind0]


coverage_N1000_1=read.csv(file = 'coverage_N1000_large.csv', header = F)
dim(coverage_N1000_1)
gamma.vec=c(rep(0, 15), rep(c(-1, 1.3, -0.9, 1.2, 1), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N1000_1=apply(coverage_N1000_1, 2, mean)
cov_1_N1000_pos=cov_rate_N1000_1[ind1]
cov_1_N1000_neg=cov_rate_N1000_1[ind0]



#scatter-plots
pdf("coverage_rates_scatter.pdf", width=12, height=8)
par(mfrow=c(2,4),cex=0.65, mai=c(0.4,0.4,0.4,0.4))
plot(1:25, cov_rate_N500_4, axes=F, 
     xlab='', ylab='', main='(a) N=500, small DIF', ylim=c(0.9,1.0), col = 'cadetblue')
axis(1, at=1:25, labels =1:25, cex.axis=.6, las=2)
axis(2, at=seq(0.9,1.00,0.01), labels = T, cex.axis=0.8, las=1)
abline(a=0.95, b=0, col='red')

plot(1:25, cov_rate_N500_1, axes=F, 
     xlab='', ylab='', main='(b)  N=500, large DIF', ylim=c(0.9,1.0), col = 'cadetblue')
axis(1, at=1:25, labels =1:25, cex.axis=.6, las=2)
axis(2, at=seq(0.9, 1.00, 0.01), labels = T, cex.axis=0.8, las=1)
abline(a=0.95, b=0, col='red')

plot(1:25, cov_rate_N1000_4, axes=F, 
     xlab='', ylab='', main='(c)  N=1000, small DIF', ylim=c(0.9,1.0), col = 'cadetblue')
axis(1, at=1:25, labels =1:25, cex.axis=.6, las=2)
axis(2, at=seq(0.9,1.00,0.01), labels = T, cex.axis=0.8, las=1)
abline(a=0.95, b=0, col='red')

plot(1:25, cov_rate_N1000_1, axes=F, 
     xlab='', ylab='', main='(d) N=1000, large DIF', ylim=c(0.9,1.0), col = 'cadetblue')
axis(1, at=1:25, labels =1:25, cex.axis=.6, las=2)
axis(2, at=seq(0.9,1.00,0.01), labels = T, cex.axis=0.8, las=1)
abline(a=0.95, b=0, col='red')

plot(6:25, cov_rate_N500_anchor4, axes=F,
     xlab='', ylab='', main='(e) N=500, small DIF', ylim=c(0.9,1.0), col = 'cadetblue')
axis(1, at=1:25, labels =1:25, cex.axis=.6, las=2)
axis(2, at=seq(0.9,1.00,0.01), labels = T, cex.axis=0.8, las=1)
abline(a=0.95, b=0, col='red')

plot(6:25, cov_rate_anchorN500_1, axes=F,
     xlab='', ylab='', main='(f) N=500, large DIF', ylim=c(0.9,1.0), col = 'cadetblue')
axis(1, at=1:25, labels =1:25, cex.axis=.6, las=2)
axis(2, at=seq(0.9,1.00,0.01), labels = T, cex.axis=0.8, las=1)
abline(a=0.95, b=0, col='red')

plot(6:25, cov_rate_N1000_anchor4, axes=F,
     xlab='', ylab='', main='(g) N=1000, small DIF', ylim=c(0.9,1.0), col = 'cadetblue')
axis(1, at=1:25, labels =1:25, cex.axis=.6, las=2)
axis(2, at=seq(0.9,1.00,0.01), labels = T, cex.axis=0.8, las=1)
abline(a=0.95, b=0, col='red')

plot(6:25, cov_rate_anchorN1000_1, axes=F,
     xlab='', ylab='', main='(h) N=1000, large DIF', ylim=c(0.9,1.0), col = 'cadetblue')
axis(1, at=1:25, labels =1:25, cex.axis=.6, las=2)
axis(2, at=seq(0.9,1.00,0.01), labels = T, cex.axis=0.8, las=1)
abline(a=0.95, b=0, col='red')
dev.off()


#### Compare TPR and FPR of the proposed hard-thresholding and the LASSO method using BIC as tuning criteria ####
#LASSO
TPR_N500_4=read.csv(file = 'TPRN500_L1_small.csv', header = F)
FPR_N500_4=read.csv(file = 'FPRN500_L1_small.csv', header = F)
TPR_N500_1=read.csv(file = 'TPRN500_L1_large.csv', header = F)
FPR_N500_1=read.csv(file = 'FPRN500_L1_large.csv', header = F)
TPR_N1000_4=read.csv(file = 'TPRN1000_L1_small.csv', header = F)
FPR_N1000_4=read.csv(file = 'FPRN1000_L1_small.csv', header = F)
TPR_N1000_1=read.csv(file = 'TPRN1000_L1_large.csv', header = F)
FPR_N1000_1=read.csv(file = 'FPRN1000_L1_large.csv', header = F)
TPR_L1=round(c(mean(TPR_N500_4[,1]),mean(TPR_N500_1[,1]),mean(TPR_N1000_4[,1]),mean(TPR_N1000_1[,1])),3)
FPR_L1=round(c(mean(FPR_N500_4[,1]),mean(FPR_N500_1[,1]),mean(FPR_N1000_4[,1]),mean(FPR_N1000_1[,1])),3)


#Proposed
TPR_N500_4=read.csv(file = 'TPRN500_hard_small.csv', header = F)
FPR_N500_4=read.csv(file = 'FPRN500_hard_small.csv', header = F)
TPR_N500_1=read.csv(file = 'TPRN500_hard_large.csv', header = F)
FPR_N500_1=read.csv(file = 'FPRN500_hard_large.csv', header = F)
TPR_N1000_4=read.csv(file = 'TPRN1000_hard_small.csv', header = F)
FPR_N1000_4=read.csv(file = 'FPRN1000_hard_small.csv', header = F)
TPR_N1000_1=read.csv(file = 'TPRN1000_hard_large.csv', header = F)
FPR_N1000_1=read.csv(file = 'FPRN1000_hard_large.csv', header = F)
TPR_prop=round(c(mean(TPR_N500_4[,1]),mean(TPR_N500_1[,1]),mean(TPR_N1000_4[,1]),mean(TPR_N1000_1[,1])),3)
FPR_prop=round(c(mean(FPR_N500_4[,1]),mean(FPR_N500_1[,1]),mean(FPR_N1000_4[,1]),mean(FPR_N1000_1[,1])),3)



#### Compare area under the ROC curves of the proposed hard-thresholding and the LASSO method using BIC as tuning criteria ####
#LASSO
gam_N500_4=read.csv(file = 'gamallN500_L1_small.csv', header = F)
gam_N500_1=read.csv(file = 'gamallN500_L1_large.csv', header = F)
gam_N1000_4=read.csv(file = 'gamallN1000_L1_small.csv', header = F)
gam_N1000_1=read.csv(file = 'gamallN1000_L1_large.csv', header = F)
res_N500_4=eva_tpr_fpr(gam=gam_N500_4, gamma.vec=c(rep(0, 15), rep(c(-0.4, 0.4, -0.5, 0.4, 0.4), 2)), rep=100)
auc_N500_4=eva_auc(tpr=res_N500_4$TPR, fpr=res_N500_4$FPR)
res_N500_1=eva_tpr_fpr(gam=gam_N500_1, gamma.vec=c(rep(0, 15), rep(c(-1, 1.3, -0.9, 1.2, 1), 2)), rep=100)
auc_N500_1=eva_auc(tpr=res_N500_1$TPR, fpr=res_N500_1$FPR)
res_N1000_4=eva_tpr_fpr(gam=gam_N1000_4, gamma.vec=c(rep(0, 15), rep(c(-0.4, 0.4, -0.5, 0.4, 0.4), 2)), rep=100)
auc_N1000_4=eva_auc(tpr=res_N1000_4$TPR, fpr=res_N1000_4$FPR)
res_N1000_1=eva_tpr_fpr(gam=gam_N1000_1, gamma.vec=c(rep(0, 15), rep(c(-1, 1.3, -0.9, 1.2, 1), 2)), rep=100)
auc_N1000_1=eva_auc(tpr=res_N1000_1$TPR, fpr=res_N1000_1$FPR)
auc_l1=round(c(auc_N500_4,auc_N500_1,auc_N1000_4,auc_N1000_1),3)



#Proposed hard-thresholding method
gam_N500_4=read.csv(file = 'gamallN500_hard_small.csv', header = F)
gam_N500_1=read.csv(file = 'gamallN500_hard_large.csv', header = F)
gam_N1000_4=read.csv(file = 'gamallN1000_hard_small.csv', header = F)
gam_N1000_1=read.csv(file = 'gamallN1000_hard_large.csv', header = F)
res_N500_4=eva_tpr_fpr(gam=gam_N500_4, gamma.vec=c(rep(0, 15), rep(c(-0.4, 0.4, -0.5, 0.4, 0.4), 2)), rep=100)
auc_N500_4=eva_auc(tpr=res_N500_4$TPR, fpr=res_N500_4$FPR)
res_N500_1=eva_tpr_fpr(gam=gam_N500_1, gamma.vec=c(rep(0, 15), rep(c(-1, 1.3, -0.9, 1.2, 1), 2)), rep=100)
auc_N500_1=eva_auc(tpr=res_N500_1$TPR, fpr=res_N500_1$FPR)
res_N1000_4=eva_tpr_fpr(gam=gam_N1000_4, gamma.vec=c(rep(0, 15), rep(c(-0.4, 0.4, -0.5, 0.4, 0.4), 2)), rep=100)
auc_N1000_4=eva_auc(tpr=res_N1000_4$TPR, fpr=res_N1000_4$FPR)
res_N1000_1=eva_tpr_fpr(gam=gam_N1000_1, gamma.vec=c(rep(0, 15), rep(c(-1, 1.3, -0.9, 1.2, 1), 2)), rep=100)
auc_N1000_1=eva_auc(tpr=res_N1000_1$TPR, fpr=res_N1000_1$FPR)
auc_prop=round(c(auc_N500_4,auc_N500_1,auc_N1000_4,auc_N1000_1),3)




























































































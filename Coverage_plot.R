
### proposed method

## N = 500, small DIF, sigma 0.5
# d1 prop 1
coverage_N500_sigma0.5_d_range1_g_abs1_g_prop1=read.csv(file = paste0(proposed_results,'coverage_N500_sigma0.5_d_range1_g_abs1_g_prop1.csv'), header = F)
dim(coverage_N500_sigma0.5_d_range1_g_abs1_g_prop1)
gamma.vec=c(rep(0, 9), c(-0.6, 0.6, -0.65, 0.7), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop1=apply(coverage_N500_sigma0.5_d_range1_g_abs1_g_prop1, 2, mean)
cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop1_pos=cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop1[ind1]
cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop1_neg=cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop1[ind0]

# d1 prop 2
coverage_N500_sigma0.5_d_range1_g_abs1_g_prop2=read.csv(file = paste0(proposed_results,'coverage_N500_sigma0.5_d_range1_g_abs1_g_prop2.csv'), header = F)
dim(coverage_N500_sigma0.5_d_range1_g_abs1_g_prop2)
gamma.vec=c(rep(0, 15), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop2=apply(coverage_N500_sigma0.5_d_range1_g_abs1_g_prop2, 2, mean)
cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop2_pos=cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop2[ind1]
cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop2_neg=cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop2[ind0]

# d1 prop 3
coverage_N500_sigma0.5_d_range1_g_abs1_g_prop3=read.csv(file = paste0(proposed_results,'coverage_N500_sigma0.5_d_range1_g_abs1_g_prop3.csv'), header = F)
dim(coverage_N500_sigma0.5_d_range1_g_abs1_g_prop3)
gamma.vec=c(rep(0, 20), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 1))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop3=apply(coverage_N500_sigma0.5_d_range1_g_abs1_g_prop3, 2, mean)
cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop3_pos=cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop3[ind1]
cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop3_neg=cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop3[ind0]

# d2 prop 1
coverage_N500_sigma0.5_d_range2_g_abs1_g_prop1=read.csv(file = paste0(proposed_results,'coverage_N500_sigma0.5_d_range2_g_abs1_g_prop1.csv'), header = F)
dim(coverage_N500_sigma0.5_d_range2_g_abs1_g_prop1)
gamma.vec=c(rep(0, 9), c(-0.6, 0.6, -0.65, 0.7), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop1=apply(coverage_N500_sigma0.5_d_range2_g_abs1_g_prop1, 2, mean)
cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop1_pos=cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop1[ind1]
cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop1_neg=cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop1[ind0]

# d1 prop 2
coverage_N500_sigma0.5_d_range2_g_abs1_g_prop2=read.csv(file = paste0(proposed_results,'coverage_N500_sigma0.5_d_range2_g_abs1_g_prop2.csv'), header = F)
dim(coverage_N500_sigma0.5_d_range2_g_abs1_g_prop2)
gamma.vec=c(rep(0, 15), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop2=apply(coverage_N500_sigma0.5_d_range2_g_abs1_g_prop2, 2, mean)
cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop2_pos=cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop2[ind1]
cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop2_neg=cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop2[ind0]

# d1 prop 3
coverage_N500_sigma0.5_d_range2_g_abs1_g_prop3=read.csv(file = paste0(proposed_results,'coverage_N500_sigma0.5_d_range2_g_abs1_g_prop3.csv'), header = F)
dim(coverage_N500_sigma0.5_d_range2_g_abs1_g_prop3)
gamma.vec=c(rep(0, 20), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 1))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop3=apply(coverage_N500_sigma0.5_d_range2_g_abs1_g_prop3, 2, mean)
cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop3_pos=cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop3[ind1]
cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop3_neg=cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop3[ind0]

## N = 500, Large DIF, sigma 0.5
# d1 prop 1
coverage_N500_sigma0.5_d_range1_g_abs2_g_prop1=read.csv(file = paste0(proposed_results,'coverage_N500_sigma0.5_d_range1_g_abs2_g_prop1.csv'), header = F)
dim(coverage_N500_sigma0.5_d_range1_g_abs2_g_prop1)
gamma.vec=c(rep(0, 9), c(-0.6, 0.6, -0.65, 0.7), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop1=apply(coverage_N500_sigma0.5_d_range1_g_abs2_g_prop1, 2, mean)
cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop1_pos=cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop1[ind1]
cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop1_neg=cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop1[ind0]

# d1 prop 2
coverage_N500_sigma0.5_d_range1_g_abs2_g_prop2=read.csv(file = paste0(proposed_results,'coverage_N500_sigma0.5_d_range1_g_abs2_g_prop2.csv'), header = F)
dim(coverage_N500_sigma0.5_d_range1_g_abs2_g_prop2)
gamma.vec=c(rep(0, 15), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop2=apply(coverage_N500_sigma0.5_d_range1_g_abs2_g_prop2, 2, mean)
cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop2_pos=cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop2[ind1]
cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop2_neg=cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop2[ind0]

# d1 prop 3
coverage_N500_sigma0.5_d_range1_g_abs2_g_prop3=read.csv(file = paste0(proposed_results,'coverage_N500_sigma0.5_d_range1_g_abs2_g_prop3.csv'), header = F)
dim(coverage_N500_sigma0.5_d_range1_g_abs2_g_prop3)
gamma.vec=c(rep(0, 20), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 1))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop3=apply(coverage_N500_sigma0.5_d_range1_g_abs2_g_prop3, 2, mean)
cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop3_pos=cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop3[ind1]
cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop3_neg=cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop3[ind0]

# d2 prop 1
coverage_N500_sigma0.5_d_range2_g_abs2_g_prop1=read.csv(file = paste0(proposed_results,'coverage_N500_sigma0.5_d_range2_g_abs2_g_prop1.csv'), header = F)
dim(coverage_N500_sigma0.5_d_range2_g_abs2_g_prop1)
gamma.vec=c(rep(0, 9), c(-0.6, 0.6, -0.65, 0.7), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop1=apply(coverage_N500_sigma0.5_d_range2_g_abs2_g_prop1, 2, mean)
cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop1_pos=cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop1[ind1]
cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop1_neg=cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop1[ind0]

# d1 prop 2
coverage_N500_sigma0.5_d_range2_g_abs2_g_prop2=read.csv(file = paste0(proposed_results,'coverage_N500_sigma0.5_d_range2_g_abs2_g_prop2.csv'), header = F)
dim(coverage_N500_sigma0.5_d_range2_g_abs2_g_prop2)
gamma.vec=c(rep(0, 15), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop2=apply(coverage_N500_sigma0.5_d_range2_g_abs2_g_prop2, 2, mean)
cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop2_pos=cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop2[ind1]
cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop2_neg=cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop2[ind0]

# d1 prop 3
coverage_N500_sigma0.5_d_range2_g_abs2_g_prop3=read.csv(file = paste0(proposed_results,'coverage_N500_sigma0.5_d_range2_g_abs2_g_prop3.csv'), header = F)
dim(coverage_N500_sigma0.5_d_range2_g_abs2_g_prop3)
gamma.vec=c(rep(0, 20), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 1))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop3=apply(coverage_N500_sigma0.5_d_range2_g_abs2_g_prop3, 2, mean)
cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop3_pos=cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop3[ind1]
cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop3_neg=cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop3[ind0]

## N = 1000, small DIF, sigma 0.5
# d1 prop 1
coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop1=read.csv(file = paste0(proposed_results,'coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop1.csv'), header = F)
dim(coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop1)
gamma.vec=c(rep(0, 9), c(-0.6, 0.6, -0.65, 0.7), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop1=apply(coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop1, 2, mean)
cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop1_pos=cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop1[ind1]
cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop1_neg=cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop1[ind0]

# d1 prop 2
coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop2=read.csv(file = paste0(proposed_results,'coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop2.csv'), header = F)
dim(coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop2)
gamma.vec=c(rep(0, 15), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop2=apply(coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop2, 2, mean)
cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop2_pos=cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop2[ind1]
cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop2_neg=cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop2[ind0]

# d1 prop 3
coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop3=read.csv(file = paste0(proposed_results,'coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop3.csv'), header = F)
dim(coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop3)
gamma.vec=c(rep(0, 20), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 1))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop3=apply(coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop3, 2, mean)
cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop3_pos=cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop3[ind1]
cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop3_neg=cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop3[ind0]

# d2 prop 1
coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop1=read.csv(file = paste0(proposed_results,'coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop1.csv'), header = F)
dim(coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop1)
gamma.vec=c(rep(0, 9), c(-0.6, 0.6, -0.65, 0.7), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop1=apply(coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop1, 2, mean)
cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop1_pos=cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop1[ind1]
cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop1_neg=cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop1[ind0]

# d1 prop 2
coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop2=read.csv(file = paste0(proposed_results,'coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop2.csv'), header = F)
dim(coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop2)
gamma.vec=c(rep(0, 15), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop2=apply(coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop2, 2, mean)
cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop2_pos=cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop2[ind1]
cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop2_neg=cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop2[ind0]

# d1 prop 3
coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop3=read.csv(file = paste0(proposed_results,'coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop3.csv'), header = F)
dim(coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop3)
gamma.vec=c(rep(0, 20), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 1))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop3=apply(coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop3, 2, mean)
cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop3_pos=cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop3[ind1]
cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop3_neg=cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop3[ind0]

## N = 1000, Large DIF, sigma 0.5
# d1 prop 1
coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop1=read.csv(file = paste0(proposed_results,'coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop1.csv'), header = F)
dim(coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop1)
gamma.vec=c(rep(0, 9), c(-0.6, 0.6, -0.65, 0.7), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop1=apply(coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop1, 2, mean)
cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop1_pos=cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop1[ind1]
cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop1_neg=cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop1[ind0]

# d1 prop 2
coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop2=read.csv(file = paste0(proposed_results,'coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop2.csv'), header = F)
dim(coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop2)
gamma.vec=c(rep(0, 15), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop2=apply(coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop2, 2, mean)
cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop2_pos=cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop2[ind1]
cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop2_neg=cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop2[ind0]

# d1 prop 3
coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop3=read.csv(file = paste0(proposed_results,'coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop3.csv'), header = F)
dim(coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop3)
gamma.vec=c(rep(0, 20), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 1))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop3=apply(coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop3, 2, mean)
cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop3_pos=cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop3[ind1]
cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop3_neg=cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop3[ind0]

# d2 prop 1
coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop1=read.csv(file = paste0(proposed_results,'coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop1.csv'), header = F)
dim(coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop1)
gamma.vec=c(rep(0, 9), c(-0.6, 0.6, -0.65, 0.7), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop1=apply(coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop1, 2, mean)
cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop1_pos=cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop1[ind1]
cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop1_neg=cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop1[ind0]

# d1 prop 2
coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop2=read.csv(file = paste0(proposed_results,'coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop2.csv'), header = F)
dim(coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop2)
gamma.vec=c(rep(0, 15), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop2=apply(coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop2, 2, mean)
cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop2_pos=cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop2[ind1]
cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop2_neg=cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop2[ind0]

# d1 prop 3
coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop3=read.csv(file = paste0(proposed_results,'coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop3.csv'), header = F)
dim(coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop3)
gamma.vec=c(rep(0, 20), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 1))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop3=apply(coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop3, 2, mean)
cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop3_pos=cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop3[ind1]
cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop3_neg=cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop3[ind0]

## anchor 5, N = 500, small DIF, sigma 0.5

# d1 prop 1
anchor5_coverage_N500_sigma0.5_d_range1_g_abs1_g_prop1=NULL
for (rep in seq(0, 90, 10)){
  cover_temp = read.csv(file = paste0(LRT_high, 'anchor5_coverage_N500_sigma0.5_d_range1_g_abs1_g_prop1_rep_', rep,'.csv'), header = F)
  cover_temp = as.matrix(cover_temp)
  anchor5_coverage_N500_sigma0.5_d_range1_g_abs1_g_prop1=rbind(anchor5_coverage_N500_sigma0.5_d_range1_g_abs1_g_prop1, cover_temp)
}
dim(anchor5_coverage_N500_sigma0.5_d_range1_g_abs1_g_prop1)
gamma.vec=c(rep(0, 9), c(-0.6, 0.6, -0.65, 0.7), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop1=apply(anchor5_coverage_N500_sigma0.5_d_range1_g_abs1_g_prop1, 2, mean)
anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop1_pos=anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop1[ind1]
anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop1_neg=anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop1[ind0]

# d1 prop 2
anchor5_coverage_N500_sigma0.5_d_range1_g_abs1_g_prop2=read.csv(file = paste0(LRT_results,'anchor5_coverage_N500_sigma0.5_d_range1_g_abs1_g_prop2.csv'), header = F)
dim(anchor5_coverage_N500_sigma0.5_d_range1_g_abs1_g_prop2)
gamma.vec=c(rep(0, 15), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop2=apply(anchor5_coverage_N500_sigma0.5_d_range1_g_abs1_g_prop2, 2, mean)
anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop2_pos=anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop2[ind1]
anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop2_neg=anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop2[ind0]
anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop2

# d1 prop 3
anchor5_coverage_N500_sigma0.5_d_range1_g_abs1_g_prop3=read.csv(file = paste0(LRT_results,'anchor5_coverage_N500_sigma0.5_d_range1_g_abs1_g_prop3.csv'), header = F)
dim(anchor5_coverage_N500_sigma0.5_d_range1_g_abs1_g_prop3)
gamma.vec=c(rep(0, 20), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 1))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop3=apply(anchor5_coverage_N500_sigma0.5_d_range1_g_abs1_g_prop3, 2, mean)
anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop3_pos=anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop3[ind1]
anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop3_neg=anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop3[ind0]
anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop3

# d2 prop 1
anchor5_coverage_N500_sigma0.5_d_range2_g_abs1_g_prop1=NULL
for (rep in seq(0, 90, 10)){
  cover_temp = read.csv(file = paste0(LRT_high, 'anchor5_coverage_N500_sigma0.5_d_range2_g_abs1_g_prop1_rep_', rep,'.csv'), header = F)
  cover_temp = as.matrix(cover_temp)
  anchor5_coverage_N500_sigma0.5_d_range2_g_abs1_g_prop1=rbind(anchor5_coverage_N500_sigma0.5_d_range2_g_abs1_g_prop1, cover_temp)
}
dim(anchor5_coverage_N500_sigma0.5_d_range2_g_abs1_g_prop1)
gamma.vec=c(rep(0, 9), c(-0.6, 0.6, -0.65, 0.7), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop1=apply(anchor5_coverage_N500_sigma0.5_d_range2_g_abs1_g_prop1, 2, mean)
anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop1_pos=anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop1[ind1]
anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop1_neg=anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop1[ind0]

# d1 prop 2
anchor5_coverage_N500_sigma0.5_d_range2_g_abs1_g_prop2=read.csv(file = paste0(LRT_results,'anchor5_coverage_N500_sigma0.5_d_range2_g_abs1_g_prop2.csv'), header = F)
dim(anchor5_coverage_N500_sigma0.5_d_range2_g_abs1_g_prop2)
gamma.vec=c(rep(0, 15), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop2=apply(anchor5_coverage_N500_sigma0.5_d_range2_g_abs1_g_prop2, 2, mean)
anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop2_pos=anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop2[ind1]
anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop2_neg=anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop2[ind0]

# d1 prop 3
anchor5_coverage_N500_sigma0.5_d_range2_g_abs1_g_prop3=read.csv(file = paste0(LRT_results,'anchor5_coverage_N500_sigma0.5_d_range2_g_abs1_g_prop3.csv'), header = F)
dim(anchor5_coverage_N500_sigma0.5_d_range2_g_abs1_g_prop3)
gamma.vec=c(rep(0, 20), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 1))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop3=apply(anchor5_coverage_N500_sigma0.5_d_range2_g_abs1_g_prop3, 2, mean)
anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop3_pos=anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop3[ind1]
anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop3_neg=anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop3[ind0]

## N = 500, Large DIF, sigma 0.5
# d1 prop 1
anchor5_coverage_N500_sigma0.5_d_range1_g_abs2_g_prop1=NULL
for (rep in seq(0, 90, 10)){
  cover_temp = read.csv(file = paste0(LRT_high, 'anchor5_coverage_N500_sigma0.5_d_range1_g_abs2_g_prop1_rep_', rep,'.csv'), header = F)
  cover_temp = as.matrix(cover_temp)
  anchor5_coverage_N500_sigma0.5_d_range1_g_abs2_g_prop1=rbind(anchor5_coverage_N500_sigma0.5_d_range1_g_abs2_g_prop1, cover_temp)
}
dim(anchor5_coverage_N500_sigma0.5_d_range1_g_abs2_g_prop1)
gamma.vec=c(rep(0, 9), c(-0.6, 0.6, -0.65, 0.7), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop1=apply(anchor5_coverage_N500_sigma0.5_d_range1_g_abs2_g_prop1, 2, mean)
anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop1_pos=anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop1[ind1]
anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop1_neg=anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop1[ind0]

# d1 prop 2
anchor5_coverage_N500_sigma0.5_d_range1_g_abs2_g_prop2=read.csv(file = paste0(LRT_results,'anchor5_coverage_N500_sigma0.5_d_range1_g_abs2_g_prop2.csv'), header = F)
dim(anchor5_coverage_N500_sigma0.5_d_range1_g_abs2_g_prop2)
gamma.vec=c(rep(0, 15), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop2=apply(anchor5_coverage_N500_sigma0.5_d_range1_g_abs2_g_prop2, 2, mean)
anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop2_pos=anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop2[ind1]
anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop2_neg=anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop2[ind0]

# d1 prop 3
anchor5_coverage_N500_sigma0.5_d_range1_g_abs2_g_prop3=read.csv(file = paste0(LRT_results,'anchor5_coverage_N500_sigma0.5_d_range1_g_abs2_g_prop3.csv'), header = F)
dim(anchor5_coverage_N500_sigma0.5_d_range1_g_abs2_g_prop3)
gamma.vec=c(rep(0, 20), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 1))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop3=apply(anchor5_coverage_N500_sigma0.5_d_range1_g_abs2_g_prop3, 2, mean)
anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop3_pos=anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop3[ind1]
anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop3_neg=anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop3[ind0]

# d2 prop 1
anchor5_coverage_N500_sigma0.5_d_range2_g_abs2_g_prop1=NULL
for (rep in seq(0, 90, 10)){
  cover_temp = read.csv(file = paste0(LRT_high, 'anchor5_coverage_N500_sigma0.5_d_range2_g_abs2_g_prop1_rep_', rep,'.csv'), header = F)
  cover_temp = as.matrix(cover_temp)
  anchor5_coverage_N500_sigma0.5_d_range2_g_abs2_g_prop1=rbind(anchor5_coverage_N500_sigma0.5_d_range2_g_abs2_g_prop1, cover_temp)
}
dim(anchor5_coverage_N500_sigma0.5_d_range2_g_abs2_g_prop1)
gamma.vec=c(rep(0, 9), c(-0.6, 0.6, -0.65, 0.7), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop1=apply(anchor5_coverage_N500_sigma0.5_d_range2_g_abs2_g_prop1, 2, mean)
anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop1_pos=anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop1[ind1]
anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop1_neg=anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop1[ind0]

# d1 prop 2
anchor5_coverage_N500_sigma0.5_d_range2_g_abs2_g_prop2=read.csv(file = paste0(LRT_results,'anchor5_coverage_N500_sigma0.5_d_range2_g_abs2_g_prop2.csv'), header = F)
dim(anchor5_coverage_N500_sigma0.5_d_range2_g_abs2_g_prop2)
gamma.vec=c(rep(0, 15), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop2=apply(anchor5_coverage_N500_sigma0.5_d_range2_g_abs2_g_prop2, 2, mean)
anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop2_pos=anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop2[ind1]
anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop2_neg=anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop2[ind0]

# d1 prop 3
anchor5_coverage_N500_sigma0.5_d_range2_g_abs2_g_prop3=read.csv(file = paste0(LRT_results,'anchor5_coverage_N500_sigma0.5_d_range2_g_abs2_g_prop3.csv'), header = F)
dim(anchor5_coverage_N500_sigma0.5_d_range2_g_abs2_g_prop3)
gamma.vec=c(rep(0, 20), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 1))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop3=apply(anchor5_coverage_N500_sigma0.5_d_range2_g_abs2_g_prop3, 2, mean)
anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop3_pos=anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop3[ind1]
anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop3_neg=anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop3[ind0]

## N = 1000, small DIF, sigma 0.5
# d1 prop 1
anchor5_coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop1=NULL
for (rep in seq(0, 90, 10)){
  cover_temp = read.csv(file = paste0(LRT_high, 'anchor5_coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop1_rep_', rep,'.csv'), header = F)
  cover_temp = as.matrix(cover_temp)
  anchor5_coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop1=rbind(anchor5_coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop1, cover_temp)
}
dim(anchor5_coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop1)
gamma.vec=c(rep(0, 9), c(-0.6, 0.6, -0.65, 0.7), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop1=apply(anchor5_coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop1, 2, mean)
anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop1_pos=anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop1[ind1]
anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop1_neg=anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop1[ind0]

# d1 prop 2
anchor5_coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop2=read.csv(file = paste0(LRT_results,'anchor5_coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop2.csv'), header = F)
dim(anchor5_coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop2)
gamma.vec=c(rep(0, 15), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop2=apply(anchor5_coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop2, 2, mean)
anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop2_pos=anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop2[ind1]
anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop2_neg=anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop2[ind0]

# d1 prop 3
anchor5_coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop3=read.csv(file = paste0(LRT_results,'anchor5_coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop3.csv'), header = F)
dim(anchor5_coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop3)
gamma.vec=c(rep(0, 20), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 1))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop3=apply(anchor5_coverage_N1000_sigma0.5_d_range1_g_abs1_g_prop3, 2, mean)
anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop3_pos=anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop3[ind1]
anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop3_neg=anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop3[ind0]

# d2 prop 1
anchor5_coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop1=NULL
for (rep in seq(0, 90, 10)){
  cover_temp = read.csv(file = paste0(LRT_high, 'anchor5_coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop1_rep_', rep,'.csv'), header = F)
  cover_temp = as.matrix(cover_temp)
  anchor5_coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop1=rbind(anchor5_coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop1, cover_temp)
}
dim(anchor5_coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop1)
gamma.vec=c(rep(0, 9), c(-0.6, 0.6, -0.65, 0.7), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop1=apply(anchor5_coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop1, 2, mean)
anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop1_pos=anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop1[ind1]
anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop1_neg=anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop1[ind0]

# d1 prop 2
anchor5_coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop2=read.csv(file = paste0(LRT_results,'anchor5_coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop2.csv'), header = F)
dim(anchor5_coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop2)
gamma.vec=c(rep(0, 15), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop2=apply(anchor5_coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop2, 2, mean)
anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop2_pos=anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop2[ind1]
anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop2_neg=anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop2[ind0]

# d1 prop 3
anchor5_coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop3=read.csv(file = paste0(LRT_results,'anchor5_coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop3.csv'), header = F)
dim(anchor5_coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop3)
gamma.vec=c(rep(0, 20), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 1))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop3=apply(anchor5_coverage_N1000_sigma0.5_d_range2_g_abs1_g_prop3, 2, mean)
anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop3_pos=anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop3[ind1]
anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop3_neg=anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop3[ind0]

## N = 1000, Large DIF, sigma 0.5
# d1 prop 1
anchor5_coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop1=NULL
for (rep in seq(0, 90, 10)){
  cover_temp = read.csv(file = paste0(LRT_high, 'anchor5_coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop1_rep_', rep,'.csv'), header = F)
  cover_temp = as.matrix(cover_temp)
  anchor5_coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop1=rbind(anchor5_coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop1, cover_temp)
}
dim(anchor5_coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop1)
gamma.vec=c(rep(0, 9), c(-0.6, 0.6, -0.65, 0.7), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop1=apply(anchor5_coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop1, 2, mean)
anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop1_pos=anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop1[ind1]
anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop1_neg=anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop1[ind0]

# d1 prop 2
anchor5_coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop2=read.csv(file = paste0(LRT_results,'anchor5_coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop2.csv'), header = F)
dim(anchor5_coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop2)
gamma.vec=c(rep(0, 15), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop2=apply(anchor5_coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop2, 2, mean)
anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop2_pos=anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop2[ind1]
anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop2_neg=anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop2[ind0]

# d1 prop 3
anchor5_coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop3=read.csv(file = paste0(LRT_results,'anchor5_coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop3.csv'), header = F)
dim(anchor5_coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop3)
gamma.vec=c(rep(0, 20), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 1))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop3=apply(anchor5_coverage_N1000_sigma0.5_d_range1_g_abs2_g_prop3, 2, mean)
anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop3_pos=anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop3[ind1]
anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop3_neg=anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop3[ind0]

# d2 prop 1
anchor5_coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop1=NULL
for (rep in seq(0, 90, 10)){
  cover_temp = read.csv(file = paste0(LRT_high, 'anchor5_coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop1_rep_', rep,'.csv'), header = F)
  cover_temp = as.matrix(cover_temp)
  anchor5_coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop1=rbind(anchor5_coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop1, cover_temp)
}
dim(anchor5_coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop1)
gamma.vec=c(rep(0, 9), c(-0.6, 0.6, -0.65, 0.7), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop1=apply(anchor5_coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop1, 2, mean)
anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop1_pos=anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop1[ind1]
anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop1_neg=anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop1[ind0]

# d1 prop 2
anchor5_coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop2=read.csv(file = paste0(LRT_results,'anchor5_coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop2.csv'), header = F)
dim(anchor5_coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop2)
gamma.vec=c(rep(0, 15), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop2=apply(anchor5_coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop2, 2, mean)
anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop2_pos=anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop2[ind1]
anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop2_neg=anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop2[ind0]

# d1 prop 3
anchor5_coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop3=read.csv(file = paste0(LRT_results,'anchor5_coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop3.csv'), header = F)
dim(anchor5_coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop3)
gamma.vec=c(rep(0, 20), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 1))
ind1=which(gamma.vec!=0)
ind0=which(gamma.vec==0)

anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop3=apply(anchor5_coverage_N1000_sigma0.5_d_range2_g_abs2_g_prop3, 2, mean)
anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop3_pos=anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop3[ind1]
anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop3_neg=anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop3[ind0]

#scatter-plots

pdf("coverage_rates_scatter.pdf", width=12, height=8)
par(mfrow=c(2,4),cex=0.65, mai=c(0.4,0.4,0.4,0.4))
plot(1:25, cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop1, axes=F, 
     xlab='', ylab='', main='(a) N=500, small DIF', ylim=c(0.85,1.0), col = 'cadetblue', pch=19)
axis(1, at=1:25, labels =1:25, cex.axis=.6, las=2)
axis(2, at=seq(0.8,1.00,0.01), labels = T, cex.axis=0.8, las=1)
abline(a=0.95, b=0, col='red')
points(1:25, cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop2, col = 'purple', pch=17)
points(1:25, cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop3, col = 'brown', pch=15)

points(1:25, cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop1, col = 'cadetblue', pch=7)
points(1:25, cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop2, col = 'purple', pch=9)
points(1:25, cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop3, col = 'brown', pch=10)


plot(1:25, cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop1, axes=F, 
     xlab='', ylab='', main='(b) N=500, large DIF', ylim=c(0.85,1.0), col = 'cadetblue', pch=19)
axis(1, at=1:25, labels =1:25, cex.axis=.6, las=2)
axis(2, at=seq(0.8,1.00,0.01), labels = T, cex.axis=0.8, las=1)
abline(a=0.95, b=0, col='red')
points(1:25, cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop2, col = 'purple', pch=17)
points(1:25, cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop3, col = 'brown', pch=15)

points(1:25, cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop1, col = 'cadetblue', pch=7)
points(1:25, cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop2, col = 'purple', pch=9)
points(1:25, cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop3, col = 'brown', pch=10)


plot(1:25, cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop1, axes=F, 
     xlab='', ylab='', main='(c) N=1000, small DIF', ylim=c(0.85,1.0), col = 'cadetblue', pch=19)
axis(1, at=1:25, labels =1:25, cex.axis=.6, las=2)
axis(2, at=seq(0.8,1.00,0.01), labels = T, cex.axis=0.8, las=1)
abline(a=0.95, b=0, col='red')
points(1:25, cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop2, col = 'purple', pch=17)
points(1:25, cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop3, col = 'brown', pch=15)

points(1:25, cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop1, col = 'cadetblue', pch=7)
points(1:25, cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop2, col = 'purple', pch=9)
points(1:25, cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop3, col = 'brown', pch=10)

# 7
plot(1:25, cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop1, axes=F, 
     xlab='', ylab='', main='(d) N=1000, large DIF', ylim=c(0.85,1.0), col = 'cadetblue', pch=19)
axis(1, at=1:25, labels =1:25, cex.axis=.6, las=2)
axis(2, at=seq(0.8,1.00,0.01), labels = T, cex.axis=0.8, las=1)
abline(a=0.95, b=0, col='red')
points(1:25, cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop2, col = 'purple', pch=17)
points(1:25, cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop3, col = 'brown', pch=15)

points(1:25, cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop1, col = 'cadetblue', pch=7)
points(1:25, cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop2, col = 'purple', pch=9)
points(1:25, cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop3, col = 'brown', pch=10)


plot(6:25, anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop1, axes=F, 
     xlab='', ylab='', main='(e) N=500, small DIF', ylim=c(0.85,1.0), col = 'cadetblue', pch=19)
axis(1, at=6:25, labels =6:25, cex.axis=.6, las=2)
axis(2, at=seq(0.8,1.00,0.01), labels = T, cex.axis=0.8, las=1)
abline(a=0.95, b=0, col='red')
points(6:25, anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop2, col = 'purple', pch=17)
points(6:25, anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs1_g_prop3, col = 'brown', pch=15)

points(6:25, anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop1, col = 'cadetblue', pch=7)
points(6:25, anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop2, col = 'purple', pch=9)
points(6:25, anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs1_g_prop3, col = 'brown', pch=10)

plot(6:25, anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop1, axes=F, 
     xlab='', ylab='', main='(f) N=500, large DIF', ylim=c(0.85,1.0), col = 'cadetblue', pch=19)
axis(1, at=6:25, labels =6:25, cex.axis=.6, las=2)
axis(2, at=seq(0.8,1.00,0.01), labels = T, cex.axis=0.8, las=1)
abline(a=0.95, b=0, col='red')
points(6:25, anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop2, col = 'purple', pch=17)
points(6:25, anchor5_cov_rate_N500_sigma0.5_d_range1_g_abs2_g_prop3, col = 'brown', pch=15)

points(6:25, anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop1, col = 'cadetblue', pch=7)
points(6:25, anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop2, col = 'purple', pch=9)
points(6:25, anchor5_cov_rate_N500_sigma0.5_d_range2_g_abs2_g_prop3, col = 'brown', pch=10)


plot(6:25, anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop1, axes=F, 
     xlab='', ylab='', main='(g) N=500, small DIF', ylim=c(0.85,1.0), col = 'cadetblue', pch=19)
axis(1, at=6:25, labels =6:25, cex.axis=.6, las=2)
axis(2, at=seq(0.8,1.00,0.01), labels = T, cex.axis=0.8, las=1)
abline(a=0.95, b=0, col='red')
points(6:25, anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop2, col = 'purple', pch=17)
points(6:25, anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs1_g_prop3, col = 'brown', pch=15)

points(6:25, anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop1, col = 'cadetblue', pch=7)
points(6:25, anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop2, col = 'purple', pch=9)
points(6:25, anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs1_g_prop3, col = 'brown', pch=10)

plot(6:25, anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop1, axes=F, 
     xlab='', ylab='', main='(h) N=500, large DIF', ylim=c(0.85,1.0), col = 'cadetblue', pch=19)
axis(1, at=6:25, labels =6:25, cex.axis=.6, las=2)
axis(2, at=seq(0.8,1.00,0.01), labels = T, cex.axis=0.8, las=1)
abline(a=0.95, b=0, col='red')
points(6:25, anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop2, col = 'purple', pch=17)
points(6:25, anchor5_cov_rate_N1000_sigma0.5_d_range1_g_abs2_g_prop3, col = 'brown', pch=15)

points(6:25, anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop1, col = 'cadetblue', pch=7)
points(6:25, anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop2, col = 'purple', pch=9)
points(6:25, anchor5_cov_rate_N1000_sigma0.5_d_range2_g_abs2_g_prop3, col = 'brown', pch=10)


dev.off()

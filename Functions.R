library(fastGHQuad)
library(quantreg)
library(MASS)


gauss=gaussHermiteData(n=35)#We use 35 Gaussian quadrature
w=gauss$w
z=gauss$x

if_zero=function(x.vec){
  if (all(x.vec==0)){
    return(1)
  }
  else{return(0)}
}


### Estimation ###
#Transform Gauss Hermite Quadrature z into th, th_k=sqrt(2)*sigma*z_k+mu for k=1,...,n for all i=1,...,N.
#Now x is a N by p binary matrix
#mu and sigma are a p-dim vectors
#Return matrix Theta of dim N by n.
tran=function(z, x, mu, sigma){
  n=length(z)
  N=dim(x)[1]
  p=dim(x)[2]
  mu.vec = c(mu)
  sigma.vec = c(sigma)
  
  mean.vec = (rep(1, N)%*%t(mu.vec))*x
  var.vec = (rep(1, N)%*%t(sigma.vec^2))*x
  zero.vec=apply(x, 1, if_zero)
  var.vec=var.vec+zero.vec
  sd.vec=var.vec^0.5
  Theta=sqrt(2)*(sd.vec%*%t(rep(1, n)))*(rep(1, N)%*%t(z))+mean.vec%*%t(rep(1, n)) #N by n 
  return(Theta)
}

evaluate_prod=function(v){# for evaluating posterior
  L=length(v)
  return(exp(sum(v, na.rm = TRUE)))
}



#Evaluate the posterior quadrature values, return a matrix N by n.
#a, b are alpha and beta of length J. dat is of dim N by J. g is gamma of dim J by p.
#Now x is a N by p binary matrix
#mu and sigma are a p-dim vectors
#z, w are the Gauss Hermite quadratures and weights of length n, 
quad=function(a, b, g, x, dat, z, mu, sigma, w){ #z is quadrature of dim n
  n=length(z)
  N=dim(dat)[1]
  J=dim(dat)[2]
  p=dim(x)[2]
  Theta=tran(z, x, mu, sigma) #Transform Gauss Hermite Quadrature z into th, th_i=sqrt(2)*z_i+mu for i=1,...,n
  th=Theta[, 1]
  ta=rep(1, N)%*%t(a)
  tb=rep(1, N)%*%t(b)
  tg=x%*%t(g)
  
  temp=(th%*%t(rep(1, J)))*ta+ tb + tg #N by J
  prob=exp(temp)/(1+exp(temp)) #dim N by J
  q=apply(log((prob^dat)*((1-prob)^(1-dat))), 1, evaluate_prod) #length N vector
  for (k in 2:n){ #matrix evaluate
    th=Theta[, k]
    temp=(th%*%t(rep(1, J)))*ta+tb+tg
    prob=exp(temp)/(1+exp(temp))
    q_c=apply(log((prob^dat)*((1-prob)^(1-dat))), 1, evaluate_prod)
    q=cbind(q, q_c)
  }
  Z=rowSums(q*(rep(1, N)%*%t(w)))
  return((q*(rep(1, N)%*%t(w)))/(Z%*%t(rep(1, n))))
}




# these target functions are for optim function later

#Evaluate target function at j=1 (i.e. anchor item), phi_j=(a_j, b_j, g_j1,...,g_jp), i.e. negative log-likelihood at phi_j
#Y_j=Y[,j]
# for a_j, b_j at fixed g at j=1 anchor item
target_function1=function(phi_j, x, Y_j, z,  mu_pres, sigma_pres, q){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  a_j=phi_j[1]
  b_j=phi_j[2]
  g_j=rep(0, p)
  
  #Evaluate parameter dependent log-likelihood
  Theta=tran(z, x, mu_pres, sigma_pres)
  temp=Theta*a_j+b_j*(rep(1,N)%*%t(rep(1, n)))+(x%*%g_j)%*%t(rep(1, n)) #N by n
  t2=Y_j%*%t(rep(1, n))*temp-log(1+exp(temp)) #The parameter dependent term in log-likelihood
  res=-sum(q*t2, na.rm = TRUE)
  return(res/N)
}

#Evaluate target function at phi_j=(a_j, b_j) given (g_j1, ... g_jp) at previous step values.
#Y_j=Y[,j]
#
# simlar to target_function at no anchoring item
target_function_ab=function(ab, g, x, Y_j, z,  mu_pres, sigma_pres, q){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  a_j=ab[1]
  b_j=ab[2]
  g_j=g
  
  #Evaluate log-likelihood
  Theta=tran(z, x, mu_pres, sigma_pres)
  temp=Theta*a_j+b_j*(rep(1,N)%*%t(rep(1, n)))+(x%*%g_j)%*%t(rep(1, n)) #N by n
  t2=Y_j%*%t(rep(1, n))*temp-log(1+exp(temp)) #The remaining term in log-likelihood
  res=-sum(q*t2, na.rm = TRUE)
  return(res/N)
}

#Evaluate target function at g_jk, given (a_j, b_j, g_j1,...g_j(k-1),g_j(k+1),...g_jp) at previous step values,
#Y_k=Y[[i in group k],j]
#ab=(a_j, b_j)
#gp=(g_j1,...g_j(k-1)) ###
#ga=(g_j(k+1),...g_jp)
#gk=g[k]
#ind_k is the indices of individuals belonging to group k
target_function_g=function(ab, k, gk, gp, x_k, Y_k, z, ind_k, mu_pres, sigma_pres, q){
  N_k=dim(x_k)[1]
  p=dim(x_k)[2]
  n=length(z)
  a_j=ab[1]
  b_j=ab[2]
  gp[k]=gk
  g_j=gp
  q_k=q[ind_k, ]
  
  #Evaluate parameter dependent log-likelihood
  Theta=tran(z, x_k, mu_pres, sigma_pres)
  temp=Theta*a_j+b_j*(rep(1,N_k)%*%t(rep(1, n)))+(x_k%*%g_j)%*%t(rep(1, n)) #N_k by n
  t2=Y_k%*%t(rep(1, n))*temp-log(1+exp(temp)) #The parameter dependent term in log-likelihood
  res=-sum(q_k*t2, na.rm = TRUE)
  return(res/N_k)
}



target_mean=function(para_k, k, x, z, mu_pres, sigma_pres, q){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  me=mu_pres
  me[k]=para_k
  
  #Evaluate log-likelihood
  Theta=tran(z, x, mu_pres, sigma_pres)
  mean.v=apply(x*(rep(1, N)%*%t(me)), 1, sum)
  var.v=apply(x*(rep(1, N)%*%t(sigma_pres^2)), 1, sum)
  zero.vec=apply(x, 1, if_zero)
  var.v=var.v+zero.vec
  sd.v=var.v^0.5
  t1=-0.5*((Theta-mean.v%*%t(rep(1,n)))/(sd.v%*%t(rep(1,n))))^2 #normal term in log-likelihood  #remove prior term
  log_likelihood=t1 #N by n
  res=-sum(q*log_likelihood)
  return(res/N)
}

target_var=function(para_k, k, x, z, mu_pres, sigma_pres, q){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  sig=sigma_pres
  sig[k]=para_k
  
  #Evaluate log-likelihood
  Theta=tran(z, x, mu_pres, sigma_pres)
  mean.v=apply(x*(rep(1, N)%*%t(mu_pres)), 1, sum)
  var.v=apply(x*(rep(1, N)%*%t(sig^2)), 1, sum)
  zero.vec=apply(x, 1, if_zero)
  var.v=var.v+zero.vec
  sd.v=var.v^0.5
  t1=-0.5*((Theta-mean.v%*%t(rep(1,n)))/(sd.v%*%t(rep(1,n))))^2 #normal term in log-likelihood  #remove prior term
  t2=sum(-0.5*log(2*pi*(var.v)))
  log_likelihood=t1 #N by n
  res=-sum(q*t1)-t2
  return(res/N)
}





#Evaluate gradient with respect to phi
grad_phi1=function(phi_j, x, Y_j, z, mu_pres, sigma_pres, q){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  a_j=phi_j[1]
  b_j=phi_j[2]
  g_j=rep(0, p)
  
  #Evaluate gradient
  Theta=tran(z, x, mu_pres, sigma_pres)
  temp=Theta*a_j+b_j*(rep(1,N)%*%t(rep(1, n)))+(x%*%g_j)%*%t(rep(1, n)) #N by n
  prob=exp(temp)/(1+exp(temp))
  
  da=sum(Theta*(Y_j%*%t(rep(1, n))-prob)*q, na.rm = TRUE)
  db=sum((Y_j%*%t(rep(1, n))-prob)*q, na.rm = TRUE)
  return(-c(da, db)/N)
}


#Evaluate gradient with respect to a_j, b_j for j not anchor items.
grad_phi_ab=function(ab, g, x, Y_j, z, mu_pres, sigma_pres, q){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  a_j=ab[1]
  b_j=ab[2]
  g_j=g
  
  #Evaluate gradient
  Theta=tran(z, x, mu_pres, sigma_pres)
  temp=Theta*a_j+b_j*(rep(1,N)%*%t(rep(1, n)))+(x%*%g_j)%*%t(rep(1, n)) #N by n
  prob=exp(temp)/(1+exp(temp))
  
  da=sum(Theta*(Y_j%*%t(rep(1, n))-prob)*q, na.rm = TRUE)
  db=sum((Y_j%*%t(rep(1, n))-prob)*q, na.rm = TRUE)
  return(-c(da, db)/N)
}


#Evaluate gradient with respect to g_k
grad_phi_g=function(ab, k, gk, gp, x_k, Y_k, z, ind_k,  mu_pres, sigma_pres, q){
  N_k=dim(x_k)[1]
  p=dim(x_k)[2]
  n=length(z)
  a_j=ab[1]
  b_j=ab[2]
  gp[k]=gk
  g_j=gp
  q_k=q[ind_k, ]
  
  #Evaluate gradient
  Theta=tran(z, x_k, mu_pres, sigma_pres)
  temp=Theta*a_j+b_j*(rep(1,N_k)%*%t(rep(1, n)))+(x_k%*%g_j)%*%t(rep(1, n)) #N by n
  prob=exp(temp)/(1+exp(temp))
  dg=sum((Y_k%*%t(rep(1, n))-prob)*q_k, na.rm = TRUE)
  return(-dg/N_k)
}




grad_mean=function(para_k, k, x, z, mu_pres, sigma_pres, q){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  me=mu_pres
  me[k]=para_k
  
  #Evaluate log-likelihood
  Theta=tran(z, x, mu_pres, sigma_pres)
  mean.v=apply(x*(rep(1, N)%*%t(me)), 1, sum)
  var.v=apply(x*(rep(1, N)%*%t(sigma_pres^2)), 1, sum)
  zero.vec=apply(x, 1, if_zero)
  var.v=var.v+zero.vec
  dmu_k=sum(((x[,k]/var.v)%*%t(rep(1,n)))*(Theta-mean.v%*%t(rep(1, n)))*q)
  return(-dmu_k/N)
}


grad_var=function(para_k, k, x, z, mu_pres, sigma_pres, q){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  sig=sigma_pres
  sig[k]=para_k
  
  #Evaluate log-likelihood
  Theta=tran(z, x, mu_pres, sigma_pres)
  mean.v=apply(x*(rep(1, N)%*%t(mu_pres)), 1, sum)
  var.v=apply(x*(rep(1, N)%*%t(sig^2)), 1, sum)
  zero.vec=apply(x, 1, if_zero)
  var.v=var.v+zero.vec
  t1=sum(-x[,k]/2/var.v)
  t2=sum(((x[,k]%*%t(rep(1,n)))*q*(Theta-mean.v%*%t(rep(1,n)))^2)/2/((var.v%*%t(rep(1,n)))^2))
  dsig_k=t1+t2
  return(-dsig_k/N)
}


#Evaluate marginal log-likelihood
mml=function(a, b, g, x, dat, z, mu, sigma, w){ #q is quadrature
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  J=dim(dat)[2]
  
  #Evaluate log-likelihood
  Theta=tran(z, x, mu, sigma) #N by n
  th=Theta[, 1]
  temp=(th%*%t(a))+rep(1, N)%*%t(b)+x%*%t(g) #N by J
  prob=exp(temp)/(1+exp(temp))
  res=apply((prob^dat)*((1-prob)^(1-dat)), 1, prod, na.rm=TRUE)*w[1]
  for (k in 2:n){
    th=Theta[, k]
    temp=(th%*%t(a))+rep(1, N)%*%t(b)+x%*%t(g) #N by J
    prob=exp(temp)/(1+exp(temp))
    res_r=apply((prob^dat)*((1-prob)^(1-dat)), 1, prod, na.rm=TRUE)*w[k]
    res=res+res_r
  }
  res=-sum(log(res*(pi)^(-0.5)))
  return(res)
}


EM_2PL_inference <- function(a, b, g, x, dat, z, mu, sigma, w, ite, tol = 0.001){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  J=dim(dat)[2]
  phi0=matrix(0, J, 2+p)
  phi1=cbind(a, b, g)
  mu0=rep(0, p)
  mu1=mu
  sigma0=rep(0, p)
  sigma1=sigma
  ind_group=list()
  dat_group=list()
  x_group=list()
  for (k in 1:p){
    ind_c=which(x[,k]==1)
    ind_group[[k]]=ind_c
    dat_group[[k]]=dat[ind_c,]
    x_group[[k]]=x[ind_c,]
  }
  while(max(c(abs(phi0-phi1)), abs(mu0-mu1), abs(sigma0-sigma1))>tol){#max
    phi0=phi1
    mu0=mu1
    sigma0=sigma1
    
    #E-step
    quadrature=quad(a=phi0[,1], b=phi0[,2], g=phi0[,c(3:(2+p))], x, dat, z, mu0, sigma0, w)
    
    #M-step
    for(j in 1:J){#update phi_j one by one
      if(j==1){ # anchor item
        par_updates = optim(par=phi1[1, c(1,2)], fn=target_function1, gr=grad_phi1, method = "L-BFGS-B", x=x, Y_j=dat[,1], z=z, mu_pres=mu0,  sigma_pres=sigma0, q=quadrature, lower=-5, upper=5)
        phi1[1, c(1,2)]=par_updates$par
        
      }
      else{
        for (l in 1:ite){
          par_updates_ab = optim(par=phi1[j, c(1,2)], fn=target_function_ab, gr=grad_phi_ab, method = "L-BFGS-B", g=phi1[j, 3:(2+p)], x=x, Y_j=dat[,j], z=z, mu_pres=mu0,  sigma_pres=sigma0, q=quadrature, lower=-5, upper=5)
          phi1[j, c(1,2)]=par_updates_ab$par
          for (k in 1:p){
            par_updates_g = optim(par=phi1[j, 2+k], fn=target_function_g, gr=grad_phi_g, method = "L-BFGS-B", ab=phi1[j, c(1,2)], k=k, gp=phi1[j, 3:(2+p)], x_k=matrix(x_group[[k]]), Y_k=dat_group[[k]][,j],ind_k=ind_group[[k]], z=z, mu_pres=mu0,  sigma_pres=sigma0, q=quadrature, lower=-5, upper=5)
            phi1[j, 2+k]=par_updates_g$par
          }
        }
      }
    }
    #ab, k, gk, gp, x_k, Y_k, z, ind_k,  mu_pres, sigma_pres, q
    #To update mu one by one
    for (k in 1:p){
      par_updates = optim(par=mu0[k], fn=target_mean, gr=grad_mean, method = "L-BFGS-B", k=k,x=x, z=z, mu_pres=mu0,  sigma_pres=sigma0, q=quadrature, lower=-5, upper=5)
      mu1[k]=par_updates$par
    }
    
    #To update sigma one by one
    for (k in 1:p){
      par_updates = optim(par=sigma0[k], fn=target_var, gr=grad_var, method = "L-BFGS-B", k=k,x=x, z=z, mu_pres=mu0,  sigma_pres=sigma0, q=quadrature, lower=0, upper=5)
      sigma1[k]=par_updates$par
    }
    
    
    MLL = mml(a=phi1[,1], b=phi1[,2], g=phi1[, c(3:(2+p))], x, dat, z, mu=mu1, sigma=sigma1, w)
    print(MLL) 
  }
  list(mu=mu1, sigma=sigma1, alpha.vec=phi1[, 1], beta.vec=phi1[, 2], gamma.vec = phi1[, c(3:(p+2))], post=quadrature);
}

## main difference of below lrt with respect to proposed:
## 1. anchor set
## 2. M step deal with anchor set
## all other steps remain the same

EM_2PL_lrt <- function(a, b, g, x, dat, z, mu, sigma, w, anchor, ite, tol = 0.001){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  J=dim(dat)[2]
  phi0=matrix(0, J, 2+p)
  phi1=cbind(a, b, g)
  mu0=rep(0, p)
  mu1=mu
  sigma0=rep(0, p)
  sigma1=sigma
  ind_group=list()
  dat_group=list()
  x_group=list()
  for (k in 1:p){
    ind_c=which(x[,k]==1)
    ind_group[[k]]=ind_c
    dat_group[[k]]=dat[ind_c,]
    x_group[[k]]=x[ind_c,]
  }
  while(max(c(abs(phi0-phi1)), abs(mu0-mu1), abs(sigma0-sigma1))>tol){#max
    phi0=phi1
    mu0=mu1
    sigma0=sigma1
    
    #E-step
    quadrature=quad(a=phi0[,1], b=phi0[,2], g=phi0[,c(3:(2+p))], x, dat, z, mu0, sigma0, w)
    
    #M-step
    for(j in 1:J){#update phi_j one by one
      if(j %in% anchor){ # anchor item
        par_updates = optim(par=phi1[1, c(1,2)], fn=target_function1, gr=grad_phi1, method = "L-BFGS-B", x=x, Y_j=dat[,1], z=z, mu_pres=mu0,  sigma_pres=sigma0, q=quadrature, lower=-5, upper=5)
        phi1[1, c(1,2)]=par_updates$par
        
      }
      else{
        for (l in 1:ite){
          par_updates_ab = optim(par=phi1[j, c(1,2)], fn=target_function_ab, gr=grad_phi_ab, method = "L-BFGS-B", g=phi1[j, 3:(2+p)], x=x, Y_j=dat[,j], z=z, mu_pres=mu0,  sigma_pres=sigma0, q=quadrature, lower=-5, upper=5)
          phi1[j, c(1,2)]=par_updates_ab$par
          for (k in 1:p){
            par_updates_g = optim(par=phi1[j, 2+k], fn=target_function_g, gr=grad_phi_g, method = "L-BFGS-B", ab=phi1[j, c(1,2)], k=k, gp=phi1[j, 3:(2+p)], x_k=matrix(x_group[[k]]), Y_k=dat_group[[k]][,j],ind_k=ind_group[[k]], z=z, mu_pres=mu0,  sigma_pres=sigma0, q=quadrature, lower=-5, upper=5)
            phi1[j, 2+k]=par_updates_g$par
          }
        }
      }
    }
    #ab, k, gk, gp, x_k, Y_k, z, ind_k,  mu_pres, sigma_pres, q
    #To update mu one by one
    for (k in 1:p){
      par_updates = optim(par=mu0[k], fn=target_mean, gr=grad_mean, method = "L-BFGS-B", k=k,x=x, z=z, mu_pres=mu0,  sigma_pres=sigma0, q=quadrature, lower=-5, upper=5)
      mu1[k]=par_updates$par
    }
    
    #To update sigma one by one
    for (k in 1:p){
      par_updates = optim(par=sigma0[k], fn=target_var, gr=grad_var, method = "L-BFGS-B", k=k,x=x, z=z, mu_pres=mu0,  sigma_pres=sigma0, q=quadrature, lower=0, upper=5)
      sigma1[k]=par_updates$par
    }
    
    
    MLL = mml(a=phi1[,1], b=phi1[,2], g=phi1[, c(3:(2+p))], x, dat, z, mu=mu1, sigma=sigma1, w)
    print(MLL) 
  }
  list(mu=mu1, sigma=sigma1, alpha.vec=phi1[, 1], beta.vec=phi1[, 2], gamma.vec = phi1[, c(3:(p+2))], post=quadrature);
}




#Return the expected second derivative of individual i.
#Here x_i=x[i,], t_i=Theta[i, ], f_i=Posterior[i, ]
#Return (2+J)p+2J by (2+J)p+2J matrix
#
# Louis Identity
first_term=function(a, b, g, x_i, t_i, mu, sigma, f_i){
  J=length(a)
  n=length(t_i)
  p=length(x_i)
  if (all(x_i==0)){
    va=1
  }
  else{va=sum(x_i*(sigma^2))}
  me=sum(x_i*mu)
  dmu=-x_i%*%t(x_i)/va #p by p
  dmusigma=-x_i%*%t(x_i)*(sum((t_i-me)*f_i)/(va^2)) #p by p
  dsigma=(x_i%*%t(x_i))/2/(va^2)-x_i%*%t(x_i)*(sum((t_i-me)^2*f_i))/(va^4) #p by p
  
  
  temp=a%*%t(t_i)+b%*%t(rep(1, n))+(g*x_i)%*%t(rep(1,n)) #J by n ##
  prob=exp(temp)/(1+exp(temp)) ##J by n
  v=prob*(1-prob) #variance matrix, J by n
  
  da=rowSums(-(rep(1, J)%*%t(t_i^2))*v*(rep(1, J)%*%t(f_i))) #J dim vector
  dab=rowSums(-(rep(1, J)%*%t(t_i))*v*(rep(1, J)%*%t(f_i))) #J dim vector
  dag=rowSums(-(rep(1, J)%*%t(t_i))*v*(rep(1, J)%*%t(f_i)))%*%t(x_i) #J by p  each row corresponds to dajdgjk, k=1,...,p, others 0.
  
  
  db=rowSums(-v*(rep(1, J)%*%t(f_i))) #J dim vector
  dbg=rowSums(-v*(rep(1, J)%*%t(f_i)))%*%t(x_i) #J by p, each row corresponds to dbjdgjk, others 0.
  
  Eva=rowSums(-v*(rep(1, J)%*%t(f_i)))
  dg=matrix(0, J*p, J*p)
  xexpansion=x_i%*%t(x_i)
  for (j in 1:J){
    dg[((j-1)*p+1):(j*p),((j-1)*p+1):(j*p)]=Eva[j]*xexpansion
  }
  
  
  #Fill into the second derivative matrix, fill in upper triangle first, then add the transpose
  res=matrix(0, (p+2)*J+2*p, (p+2)*J+2*p)
  res[1:p, 1:p]=dmu
  res[1:p, (p+1):(2*p)]=dmusigma
  res[(p+1):(2*p), 1:p]=t(dmusigma)
  #res[(p+1):(2*p), 1:p]=t(dmusigma)
  res[(p+1):(2*p), (p+1):(2*p)]=dsigma
  res[(2*p+1):(2*p+2*J),(2*p+1):(2*p+2*J)]=diag(c(da, db))
  res[(2*p+1):(2*p+J),(2*p+J+1):(2*p+2*J)]=diag(dab)
  res[(2*p+J+1):(2*p+2*J), (2*p+1):(2*p+J)]=diag(dab)
  resag=matrix(0, J, J*p)
  resbg=matrix(0, J, J*p)
  for (j in 1:J){
    resag[j,((j-1)*p+1):(j*p)]=dag[j,]
    resbg[j,((j-1)*p+1):(j*p)]=dbg[j,]
  }
  res[(2*p+1):(2*p+J),(2*p+2*J+1):(2*p+2*J+J*p)]=resag
  res[(2*p+2*J+1):(2*p+2*J+J*p), (2*p+1):(2*p+J)]=t(resag)
  res[(2*p+J+1):(2*p+2*J),(2*p+2*J+1):(2*p+2*J+J*p)]=resbg
  res[(2*p+2*J+1):(2*p+2*J+J*p), (2*p+J+1):(2*p+2*J)]=t(resbg)
  res[(2*p+2*J+1):(2*p+2*J+J*p),(2*p+2*J+1):(2*p+2*J+J*p)]=dg
  return(res)
}




#Return the expected (dl/dphi)(dl/dphi)^t  of individual i.
#Here x_i=x.vec[i,], t_i=Theta[i, ], f_i=Posterior[i, ], y_i=dat[i, ]
#Return (2+J)p+2J by (2+J)p+2J matrix
second_term=function(a, b, g, x_i, y_i, t_i, mu, sigma, f_i){
  J=length(a)
  n=length(t_i)
  p=length(x_i)
  
  if (all(x_i==0)){
    va=1
  }
  else{va=sum(x_i*(sigma^2))}
  me=sum(x_i*mu)
  
  temp=a%*%t(t_i)+b%*%t(rep(1, n))+(g*x_i)%*%t(rep(1,n)) #J by n ##
  prob=exp(temp)/(1+exp(temp)) ##J by n
  
  da1=(rep(1, J)%*%t(t_i))*(y_i%*%t(rep(1, n))-prob) #dl/da  J by n
  db1=(y_i%*%t(rep(1, n))-prob) #dl/db  J by n
  dg1=NULL #J*p by n
  dg2=0 #J*p by J*p
  dif=y_i%*%t(rep(1,n))-prob
  for (i in 1:n){
    #pro=prob[, i]#J dim
    diff=dif[,i]#J dim
    dg_c=t(diff%*%t(x_i)) #p by J of one dimension of numerical expectation
    dg_c=c(dg_c) # convert to a J*p dim vector column-wise
    dg1=cbind(dg1, dg_c)
    dg_c=(dg_c%*%t(dg_c))*f_i[i] #J*p by J*p result of one dim of numerical expectation
    dg2=dg2+dg_c
  }
  
  dsigma1=-(x_i/2/va)%*%t(rep(1, n))+(x_i/2/(va^2))%*%t(((t_i-me)^2)) #p by n
  dmu1=(x_i/va)%*%t((t_i-me)) #p by n
  
  dmu2=(dmu1*(rep(1, p)%*%t(f_i)))%*%t(dmu1) #store E[(dl/dmu)^2], p by p
  dsigma2=(dsigma1*(rep(1, p)%*%t(f_i)))%*%t(dsigma1) #store E[(dl/dsigma^2)^2], p by p
  
  dmusigma=(dmu1*(rep(1, p)%*%t(f_i)))%*%t(dsigma1) #store E[(dl/dmu)(dl/dsigma)], p by p
  dmua=(dmu1*(rep(1, p)%*%t(f_i)))%*%t(da1) #store E[(dl/dmu)(dl/da)^t], p by J 
  dmub=(dmu1*(rep(1, p)%*%t(f_i)))%*%t(db1) #store E[(dl/dmu)(dl/db)^t], p by J 
  dmug=(dmu1*(rep(1, p)%*%t(f_i)))%*%t(dg1) #store E[(dl/dmu)(dl/dg)^t] p by J*p
  dsigmaa=(dsigma1*(rep(1, p)%*%t(f_i)))%*%t(da1) #store E[(dl/dsigma^2)(dl/da0)^t], p by J
  dsigmab=(dsigma1*(rep(1, p)%*%t(f_i)))%*%t(db1) #store E[(dl/dsigma^2)(dl/db)^t], p by J
  dsigmag=(dsigma1*(rep(1, p)%*%t(f_i)))%*%t(dg1) #store E[(dl/dsigma^2)(dl/dg)^t], p by J
  
  
  da2=(da1*(rep(1, J)%*%t(f_i)))%*%t(da1) #store E[(dl/da)(dl/da)^t], J by J
  dab=(da1*(rep(1, J)%*%t(f_i)))%*%t(db1) #store E[(dl/da)(dl/db)^t], J by J
  dag=(da1*(rep(1, J)%*%t(f_i)))%*%t(dg1) #store E[(dl/da)(dl/dg)^t], J by J*p
  
  db2=(db1*(rep(1, J)%*%t(f_i)))%*%t(db1)#store E[(dl/db)(dl/db)^t], J by J
  dbg=(db1*(rep(1, J)%*%t(f_i)))%*%t(dg1)#store E[(dl/db)(dl/dg)^t], J by J*p
  
  
  #Fill into the second derivative matrix
  res=matrix(0, (p+2)*J+2*p, (p+2)*J+2*p)
  res[1:p, 1:p]=dmu2
  res[1:p, (p+1):(2*p)]=dmusigma
  res[(p+1):(2*p), 1:p]=t(dmusigma)
  res[1:p, (2*p+1):(2*p+J)]=dmua
  res[(2*p+1):(2*p+J), 1:p]=t(dmua)
  res[1:p, (2*p+J+1):(2*p+2*J)]=dmub
  res[(2*p+J+1):(2*p+2*J), 1:p]=t(dmub)
  res[1:p, (2*p+2*J+1):(2*p+2*J+J*p)]=dmug
  res[(2*p+2*J+1):(2*p+2*J+J*p), 1:p]=t(dmug)
  res[(p+1):(2*p), (p+1):(2*p)]=dsigma2
  res[(p+1):(2*p), (2*p+1):(2*p+J)]=dsigmaa
  res[(2*p+1):(2*p+J), (p+1):(2*p)]=t(dsigmaa)
  res[(p+1):(2*p), (2*p+J+1):(2*p+2*J)]=dsigmab
  res[(2*p+J+1):(2*p+2*J), (p+1):(2*p)]=t(dsigmab)
  res[(p+1):(2*p), (2*p+2*J+1):(2*p+2*J+p*J)]=dsigmag
  res[(2*p+2*J+1):(2*p+2*J+p*J), (p+1):(2*p)]=t(dsigmag)
  res[(2*p+1):(2*p+J),(2*p+1):(2*p+J)]=da2
  res[(2*p+1):(2*p+J),(2*p+J+1):(2*p+2*J)]=dab
  res[(2*p+J+1):(2*p+2*J), (2*p+1):(2*p+J)]=t(dab)
  res[(2*p+1):(2*p+J),(2*p+2*J+1):(2*p+2*J+p*J)]=dag
  res[(2*p+2*J+1):(2*p+2*J+p*J), (2*p+1):(2*p+J)]=t(dag)
  res[(2*p+J+1):(2*p+2*J),(2*p+J+1):(2*p+2*J)]=db2
  res[(2*p+J+1):(2*p+2*J),(2*p+2*J+1):(2*p+2*J+J*p)]=dbg
  res[(2*p+2*J+1):(2*p+2*J+J*p), (2*p+J+1):(2*p+2*J)]=t(dbg)
  res[(2*p+2*J+1):(2*p+2*J+J*p), (2*p+2*J+1):(2*p+2*J+J*p)]=dg2
  return(res)
}



#Return the E(dl/dphi)E(dl/dphi)^t  of individual i.
#Here x_i=x.vec[i, ], t_i=Theta[i, ], f_i=Posterior[i, ], y_i=dat[i, ]
#Return (2+J)p+2J by (2+J)p+2J matrix
third_term=function(a, b, g, x_i, y_i, t_i, mu, sigma, f_i){
  J=length(a)
  n=length(t_i)
  p=length(x_i)
  if (all(x_i==0)){
    va=1
  }
  else{va=sum(x_i*(sigma^2))}
  me=sum(x_i*mu)
  
  temp=a%*%t(t_i)+b%*%t(rep(1, n))+(g*x_i)%*%t(rep(1,n)) #J by n ##
  prob=exp(temp)/(1+exp(temp)) ##J by n
  dsigma=apply(-(x_i/2/va)%*%t(rep(1, n))+(x_i/2/(va^2))%*%t(((t_i-me)^2)*f_i), 1, sum) #p dim
  dmu=apply((x_i/va)%*%t((t_i-me)*f_i), 1, sum) #p dim
  
  da=apply((rep(1, J)%*%t(f_i))*(rep(1, J)%*%t(t_i))*(y_i%*%t(rep(1, n))-prob),1,sum) #dl/da  J dim
  db=apply((rep(1, J)%*%t(f_i))*(y_i%*%t(rep(1, n))-prob),1,sum) #dl/db  J dim
  dg=0 #J*p by n
  dif=y_i%*%t(rep(1,n))-prob
  for (i in 1:n){
    diff=dif[,i] #J dim
    dg_c=t(diff%*%t(x_i)*f_i[i]) #p by J of one dimension of numerical expectation
    dg_c=c(dg_c) #convert to a J*p dim vector column-wise
    dg=dg+dg_c
  }
  
  dl1=c(dmu, dsigma, da, db, dg)
  return(dl1%*%t(dl1))
}




#return information matrix using Louis identity
#q is the posterior, N by n
information=function(a, b, g, x, dat, z, mu, sigma, q, anchor=c(0)){
  J=length(a)
  n=dim(q)[2]
  N=dim(dat)[1]
  p=dim(x)[2]
  Theta=tran(z, x, mu, sigma)
  res=matrix(0, (p+2)*J+2*p, (p+2)*J+2*p)
  if (any(anchor == 0)){
    anchor_idx = 1
  }
  else{
    anchor_idx=1:(p*anchor[length(anchor)])
  }
  for (i in 1:N){
    x_i=x[i, ]
    y_i=dat[i, ]
    t_i=Theta[i,]
    f_i=q[i, ]
    res=res+first_term(a, b, g, x_i, t_i, mu, sigma, f_i) + second_term(a, b, g, x_i, y_i, t_i, mu, sigma, f_i) - third_term(a, b, g, x_i, y_i, t_i, mu, sigma, f_i)
  }
  #res=res[,-(2*p+2*J+l_anchor+1)] # changed with possible anchor accordingly
  #res=res[-(2*p+2*J+l_anchor+1), ]
  
  res=res[ , -(2*p+2*J+anchor_idx)]
  res=res[-(2*p+2*J+anchor_idx), ]
  
  return(-res)
}


#Training LASSO method

#Evaluate gradient with respect to phi
grad_phi=function(phi_j, x, Y_j, z, mu_pres, sigma_pres, q){
  N=dim(x)[1]
  n=length(z)
  a_j=phi_j[1]
  b_j=phi_j[2]
  g_j=phi_j[3]
  
  #Evaluate gradient
  Theta=tran(z, x, mu_pres, sigma_pres)
  temp=Theta*a_j-b_j*(rep(1,N)%*%t(rep(1, n)))+(x%*%t(rep(1, n)))*g_j #N by n
  p=exp(temp)/(1+exp(temp))
  da=sum(Theta*(Y_j%*%t(rep(1, n))-p)*q)/N
  db=sum((-Y_j%*%t(rep(1, n))+p)*q)/N
  dg=sum((x%*%t(rep(1, n)))*(Y_j%*%t(rep(1, n))-p)*q)/N
  return(-c(da, db, dg))
}

##Evaluate Hessian for better initialization of learning rate
hessian=function(a, b, g, x, dat, z, mu, sigma, q){
  N=dim(dat)[1]
  J=length(a)
  Theta=tran(z, x, mu, sigma)
  res=matrix(0, 3*J+2, 3*J+2)
  #res=matrix(0, 3*J, 3*J)
  for (i in 1:N){
    x_i=x[i]
    y_i=dat[i, ]
    t_i=Theta[i,]
    f_i=q[i, ]
    res=res+first_term(a, b, g, x_i, t_i, mu, sigma, f_i)
  }
  
  return(-diag(res)/N)
}



EM_2PL_L1 <- function(a, b, g, x, dat, z, mu, sigma, w, lmd, tau, k, tol){
  J=length(a)
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  phi0=matrix(0, J, 2+p)
  phi1=cbind(a, b, g)
  mu0=rep(0, p)
  mu1=mu
  sigma0=rep(0, p)
  sigma1=sigma
  
  ind1=which(x==1)
  n1=length(ind1)
  ind_group=list()
  dat_group=list()
  x_group=list()
  for (s in 1:p){
    ind_c=which(x[,s]==1)
    ind_group[[s]]=ind_c
    dat_group[[s]]=dat[ind_c,]
    x_group[[s]]=x[ind_c,]
  }
  alpha = rep(c(1.3,1.4,1.5,1.7,1.6), 5)
  gamma = c(rep(0, 15), rep(c(-0.6, 0.6, -0.65, 0.7, 0.65), 2))
  while(max(c(max(abs(phi0-phi1)), abs(mu0-mu1), abs(sigma0-sigma1)))>tol){
    phi0=phi1
    mu0=mu1
    sigma0=sigma1
    #E-step
    quadrature=quad(a=phi0[,1], b=phi0[,2], g=phi0[,3], x, dat, z, mu0, sigma0, w)
    
    #M-step
    for(j in 1:J){#update phi_j one by one
      a_c=phi1[j, 1]
      b_c=phi1[j, 2]
      g_c=phi1[j, 3]
      for(i in 1:k){
        grad_ab=grad_phi_ab(c(a_c, b_c), g_c, x, Y_j=dat[, j], z, mu_pres=mu0, sigma_pres=sigma0, q=quadrature) 
        grad_g=grad_phi_g(c(a_c, b_c), k=1, gk=g_c, gp=g_c, x_k=matrix(x_group[[i]]), Y_k=dat_group[[i]][,j], z, ind_k=ind_group[[i]], mu_pres=mu0, sigma_pres=sigma0, q=quadrature)

        a_c= a_c - 0.1*grad_ab[1]
        b_c= b_c - tau[J+j]*grad_ab[2]
        g_c= g_c - tau[2*J+j]*grad_g
        if (abs(g_c)< tau[2*J+j]*lmd/2){
          g_c=0
        }
        else if(g_c<  - tau[2*J+j]*lmd/2){
          g_c= g_c + tau[2*J+j]*lmd/2
        }
        else if (g_c >  tau[2*J+j]*lmd/2) {
          g_c= g_c - tau[2*J+j]*lmd/2
        }
      }
      phi1[j, ]=c(a_c, b_c, g_c)
    }
    
    for (k in 1:p){
      par_updates = optim(par=mu0[k], fn=target_mean, gr=grad_mean, method = "L-BFGS-B", k=k, x=x, z=z, mu_pres=mu0,  sigma_pres=sigma0, q=quadrature, lower=-5, upper=5)
      mu1[k]=par_updates$par
    }
    
    #To update sigma one by one
    for (k in 1:p){
      par_updates = optim(par=sigma0[k], fn=target_var, gr=grad_var, method = "L-BFGS-B", k=k, x=x, z=z, mu_pres=mu0,  sigma_pres=sigma0, q=quadrature, lower=0, upper=5)
      sigma1[k]=par_updates$par
    }
    
    #MLL = mml_with_penalty(a=phi1[,1], b=phi1[,2], g=phi1[,3], x, dat, z, mu, w)
    #print(MLL)
  }
  list(mu=mu, sigma=sigma, alpha.vec=phi1[,1], beta.vec = phi1[,2], gamma.vec = phi1[,3], post=quadrature);
}

negative_log_likelihood=function(a, b, g, x, dat, z, mu, mu_pres, sigma, sigma_pres, w, q){ #q is quadrature
  N=length(x)
  n=length(z)
  J=dim(dat)[2]
  
  #Evaluate log-likelihood
  Theta=tran(z, x, mu_pres, sigma_pres) #N by n
  t1=log((rep(1, N)%*%t(w)))
  th=Theta[, 1]
  temp=(th%*%t(rep(1, J)))*(rep(1, N)%*%t(a))-rep(1, N)%*%t(b)+(x%*%t(rep(1, J)))*(rep(1, N)%*%t(g)) #N by J
  t2=rowSums(dat*temp-log(1+exp(temp))) #The remaining term in log-likelihood
  for (k in 2:n){
    th=Theta[, k]
    temp=(th%*%t(rep(1, J)))*(rep(1, N)%*%t(a))-rep(1, N)%*%t(b)+(x%*%t(rep(1, J)))*(rep(1, N)%*%t(g))
    t2_c=rowSums(dat*temp-log(1+exp(temp)))
    t2=cbind(t2, t2_c)
  }
  log_likelihood=t1+t2 #N by n
  res=-sum(q*log_likelihood)
  return(res)
  
}

##de-bias the non-zero entries obtained from LASSO method/hard-thresholding method so as to compute the MLL for AIC/BIC
EM_2PL_non_zero <- function(a, b, g, x, dat, z, mu, sigma, w, tol){
  J=length(a)
  p=dim(x)[2]
  
  ind1=which(x==1)
  n1=length(ind1)
  quadrature=quad(a, b, g, x, dat, z, mu, sigma, w)
  idx0=which(g==0)
  MLL = negative_log_likelihood(a, b, g, x, dat, z, mu=mu, mu_pres=mu, sigma=sigma, sigma_pres=sigma, w, q=quadrature)#, q=quadrature)
  #mml=function(a, b, g, x, dat, z, mu, sigma, w)
  #print(MLL)
  ind_group=list()
  dat_group=list()
  x_group=list()
  for (k in 1:p){
    ind_c=which(x[,k]==1)
    ind_group[[k]]=ind_c
    dat_group[[k]]=dat[ind_c,]
    x_group[[k]]=x[ind_c,]
  }
  
  phi0=matrix(0, J, 3)
  phi1=cbind(a, b, g)
  mu0=c(0)
  sigma0=c(0)
  while(max(c(max(abs(phi0-phi1)), abs(mu0-mu), abs(sigma0-sigma)))>tol){#max
    phi0=phi1
    mu0=mu
    sigma0=sigma
    #E-step
    quadrature=quad(a=phi0[,1], b=phi0[,2], g=phi0[,3], x, dat, z, mu0, sigma0, w)
    
    #M-step
    for(j in 1:J){#update phi_j one by one
      if (j %in% idx0){
        par_updates = optim(par=phi1[j,c(1,2)], fn=target_function1, gr=grad_phi1, method = "L-BFGS-B", x=x, Y_j=dat[,j], z, mu_pres=mu0, sigma_pres=sigma0, q=quadrature, lower=-3, upper=3)
        phi1[j,c(1,2)]=par_updates$par
      }
      
      else{
        for (l in 1:10){#ite=10
          par_updates_ab = optim(par=phi1[j, c(1,2)], fn=target_function_ab, gr=grad_phi_ab, method = "L-BFGS-B", g=phi1[j, 3:(2+p)], x=x, Y_j=dat[,j], z=z, mu_pres=mu0,  sigma_pres=sigma0, q=quadrature, lower=-5, upper=5)
          phi1[j, c(1,2)]=par_updates_ab$par
          for (k in 1:p){
            par_updates_g = optim(par=phi1[j, 2+k], fn=target_function_g, gr=grad_phi_g, method = "L-BFGS-B", ab=phi1[j, c(1,2)], k=k, gp=phi1[j, 3:(2+p)], x_k=matrix(x_group[[k]]), Y_k=dat_group[[k]][,j],ind_k=ind_group[[k]], z=z, mu_pres=mu0,  sigma_pres=sigma0, q=quadrature, lower=-5, upper=5)
            phi1[j, 2+k]=par_updates_g$par
          }
        }
      }
    }
    #ab, k, gk, gp, x_k, Y_k, z, ind_k,  mu_pres, sigma_pres, q
    #To update mu one by one
    for (k in 1:p){
      par_updates = optim(par=mu0[k], fn=target_mean, gr=grad_mean, method = "L-BFGS-B", k=k,x=x, z=z, mu_pres=mu0,  sigma_pres=sigma0, q=quadrature, lower=-5, upper=5)
      mu[k]=par_updates$par
    }
    
    #To update sigma one by one
    for (k in 1:p){
      par_updates = optim(par=sigma0[k], fn=target_var, gr=grad_var, method = "L-BFGS-B", k=k,x=x, z=z, mu_pres=mu0,  sigma_pres=sigma0, q=quadrature, lower=0, upper=5)
      sigma[k]=par_updates$par
    }

    MLL = mml(a=phi1[,1], b=phi1[,2], g=phi1[,3], x, dat, z, mu, sigma, w)
    print(MLL) 
  }
  list(MLL=MLL);
}


## For proposed p-value based method

coverage_proposed_inference=function(seed, N, sig, gamma.vec, beta.vec,  alpha.vec, mu, sigma, x, ite,  sig_list){
  set.seed(seed)
  J=length(alpha.vec)
  p=dim(x)[2]
  n_sig=length(sig_list)
  gamm=c(t(gamma.vec))
  ind1_gam=which(gamm!=0)
  ind0_gam=which(gamm==0)
  n_p_gam=length(ind1_gam)
  n_n_gam=length(ind0_gam)
  
  cover=NULL
  cover_gam=NULL
  p_value=NULL
  p_value_gam=NULL
  alpha=NULL
  beta=NULL
  gamma=NULL
  rgamma=NULL
  TPR_gam=NULL
  FPR_gam=NULL
  r_mu=NULL
  r_sigma=NULL
  r_info=NULL
  p_val_hotelling=NULL
  for (l in 1:100){
    
    #Generate data
    mu.vec = c(mu)
    sigma.vec = c(sigma)
    mu.v=x*rep(1, N)%*%t(mu.vec)
    var.v=x*(rep(1, N)%*%t(sigma.vec^2))
    zero.vec=apply(x, 1, if_zero)
    var.v=var.v+zero.vec
    sigma.v=var.v^0.5
    th=rnorm(n=N,mean=mu.v, sd=sigma.v)
    temp=(th%*%t(rep(1, J)))*(rep(1, N)%*%t(alpha.vec))+rep(1, N)%*%t(beta.vec)+x%*%t(gamma.vec) #N by J
    prob=exp(temp)/(1+exp(temp)) #dim N by J
    dat = matrix(0, N, J);
    dat[] = rbinom(N*J, 1, prob);
    
    #Train model
    r=EM_2PL_inference(a=alpha.vec,  b=beta.vec, g=gamma.vec, x=x, dat=dat, z, mu=mu, sigma=sigma, w,ite,  tol = 0.005)
    g_hat=r$gamma.vec
    mu_hat=r$mu
    sigma_hat=r$sigma
    a_hat=r$alpha.vec
    b_hat=r$beta.vec
    
    #Evaluate variance
    info=information(a=a_hat, b=b_hat, g=g_hat, x=x, dat, z, mu=mu_hat, sigma=sigma_hat, q=r$post)
    varian=solve(info)
    eig=eigen(varian)
    eig=eig$values
    
    while (any(eig<0)){
      
      #Generate data
      mu.v=apply(x*rep(1, N)%*%t(mu), 1, sum)
      var.v=apply(x*(rep(1, N)%*%t(sigma^2)), 1, sum)
      zero.vec=apply(x, 1, if_zero)
      var.v=var.v+zero.vec
      sigma.v=var.v^0.5
      th=rnorm(n=N,mean=mu.v, sd=sigma.v)
      temp=(th%*%t(rep(1, J)))*(rep(1, N)%*%t(alpha.vec))+rep(1, N)%*%t(beta.vec)+x%*%t(gamma.vec) #N by J
      prob=exp(temp)/(1+exp(temp)) #dim N by J
      dat = matrix(0, N, J);
      dat[] = rbinom(N*J, 1, prob);
      
      
      
      #Train model
      r=EM_2PL_inference(a=alpha.vec,  b=beta.vec, g=gamma.vec, x=x, dat=dat, z, mu=mu, sigma=sigma, w, ite, tol = 0.005)
      g_hat=r$gamma.vec
      mu_hat=r$mu
      sigma_hat=r$sigma
      a_hat=r$alpha.vec
      b_hat=r$beta.vec
      
      
      #Evaluate variance
      info=information(a=a_hat, b=b_hat, g=g_hat, x=x, dat, z, mu=mu_hat, sigma=sigma_hat, q=r$post)
      varian=solve(info)
      eig=eigen(varian)
      eig=eig$values
    }
    
    r_info=rbind(r_info, info)
    r_mu=c(r_mu, mu_hat)
    r_sigma=c(r_sigma, sigma_hat)
    alpha=rbind(alpha, a_hat)
    beta=rbind(beta, b_hat)
    gamma=rbind(gamma, g_hat)
    
    
    #simulate new gammas to get p-values and CIs
    transform_p=NULL
    rqfit=rq(g_hat ~ -1 + r$alpha.vec)
    transform_p=g_hat-rqfit$coefficients[1]*r$alpha.vec
    
    sim=mvrnorm(n=10000, mu=c(mu_hat, sigma_hat, a_hat,  b_hat, rep(0, (J-1)*p)), Sigma=varian)
    sim_a=sim[, (2*p+1):(2*p+J)]
    sim_gam=cbind(rep(0, 10000), sim[ ,(2*p+2*J+1):(2*p+2*J+(J-1)*p)])
    
    gamm_transform=NULL
    for (k in 1:10000){
      rqfit = rq(sim_gam[k, ] ~ -1+ sim_a[k,])
      gamm_transform_c=sim_gam[k, ]-rqfit$coefficients[1]*sim_a[k,]
      gamm_transform=cbind(gamm_transform, gamm_transform_c)
    }
  
    #To get p-values
    p_value_curr=c()
    for (j in 1:J){
      stat=transform_p[j]
      dist=as.vector(gamm_transform[j, ])
      p_value_c=length(c(which(dist>abs(stat)), which(dist<   -abs(stat))))/10000
      p_value_curr=c(p_value_curr, p_value_c)
    }
    p_value=rbind(p_value, p_value_curr)
    
    #Coverage rate
    cover_c=c()
    for (j in 1:J){
      upp_c=transform_p[j]+quantile(gamm_transform[j, ], 1-sig/2)
      low_c=transform_p[j]+quantile(gamm_transform[j, ], sig/2)
      if (gamma.vec[j]>=low_c & gamma.vec[j]<=upp_c){
        cover_c=c(cover_c, 1)
      }
      else{cover_c=c(cover_c, 0)}
    }
    cover=rbind(cover, cover_c)
    
    #To get ROC curve
    TPR_c=c()
    FPR_c=c()
    for (s in sig_list){
      upper=c()
      lower=c()
      for (j in 1:(J*p)){
        upper=c(upper, transform_p[j]+quantile(gamm_transform[j, ], 1-s/2))
        lower=c(lower, transform_p[j]+quantile(gamm_transform[j, ], s/2))
      }
      
      upper_p=upper[ind1_gam]
      lower_p=lower[ind1_gam]
      upper_n=upper[ind0_gam]
      lower_n=lower[ind0_gam]
      TPR_c=c(TPR_c, length(c(which(upper_p<0), which(lower_p>0)))/n_p_gam)
      FPR_c=c(FPR_c, length(c(which(upper_n<0), which(lower_n>0)))/n_n_gam)
    }
    TPR_gam=rbind(TPR_gam, TPR_c)
    FPR_gam=rbind(FPR_gam, FPR_c)
    
  
    
  }
  
  
  return(list(p_val_hotelling=p_val_hotelling,  alpha=alpha, beta=beta, gamma=gamma, mu=r_mu, sigma=r_sigma,
              cover=cover, p_value=p_value, TPR_gam=TPR_gam, FPR_gam=FPR_gam,r_info=r_info))
}

#For LRT method
coverage_LRT=function(seed, N, sig, gamma.vec, beta.vec,  alpha.vec, mu, sigma, x, anchor, ite,  sig_list){
  set.seed(seed)
  J=length(alpha.vec)
  p=dim(x)[2]
  n_sig=length(sig_list)
  gamm=c(t(gamma.vec))
  ind1_gam=which(gamm!=0)
  ind0_gam=which(gamm==0)
  n_p_gam=length(ind1_gam)
  n_n_gam=length(ind0_gam)
  
  cover=NULL
  cover_gam=NULL
  p_value=NULL
  p_value_gam=NULL
  alpha=NULL
  beta=NULL
  gamma=NULL
  rgamma=NULL
  TPR=NULL
  FPR=NULL
  r_mu=NULL
  r_sigma=NULL
  r_info=NULL
  p_val_hotelling=NULL
  for (l in 1:100){
    
    #Generate data
    mu.vec = c(mu)
    sigma.vec = c(sigma)
    mu.v=x*rep(1, N)%*%t(mu.vec)
    var.v=x*(rep(1, N)%*%t(sigma.vec^2))
    zero.vec=apply(x, 1, if_zero)
    var.v=var.v+zero.vec
    sigma.v=var.v^0.5
    th=rnorm(n=N,mean=mu.v, sd=sigma.v)
    temp=(th%*%t(rep(1, J)))*(rep(1, N)%*%t(alpha.vec))+rep(1, N)%*%t(beta.vec)+x%*%t(gamma.vec) #N by J
    prob=exp(temp)/(1+exp(temp)) #dim N by J
    dat = matrix(0, N, J);
    dat[] = rbinom(N*J, 1, prob);
    
    #Train model
    r=EM_2PL_lrt(a=alpha.vec,  b=beta.vec, g=gamma.vec, x=x, dat=dat, z, mu=mu, sigma=sigma, w, anchor=anchor, ite,  tol = 0.005)
    g_hat=r$gamma.vec
    mu_hat=r$mu
    sigma_hat=r$sigma
    a_hat=r$alpha.vec
    b_hat=r$beta.vec
    
    #Evaluate variance
    info=information(a=a_hat, b=b_hat, g=g_hat, x=x, dat, z, mu=mu_hat, sigma=sigma_hat, q=r$post, anchor=anchor)
    varian=solve(info)
    eig=eigen(varian)
    eig=eig$values
    
    while (any(eig<0)){
      
      #Generate data
      mu.v=apply(x*rep(1, N)%*%t(mu), 1, sum)
      var.v=apply(x*(rep(1, N)%*%t(sigma^2)), 1, sum)
      zero.vec=apply(x, 1, if_zero)
      var.v=var.v+zero.vec
      sigma.v=var.v^0.5
      th=rnorm(n=N,mean=mu.v, sd=sigma.v)
      temp=(th%*%t(rep(1, J)))*(rep(1, N)%*%t(alpha.vec))+rep(1, N)%*%t(beta.vec)+x%*%t(gamma.vec) #N by J
      prob=exp(temp)/(1+exp(temp)) #dim N by J
      dat = matrix(0, N, J);
      dat[] = rbinom(N*J, 1, prob);
      
      
      
      #Train model
      r=EM_2PL_lrt(a=alpha.vec,  b=beta.vec, g=gamma.vec, x=x, dat=dat, z, mu=mu, sigma=sigma, w, anchor=anchor, ite, tol = 0.005)
      g_hat=r$gamma.vec
      mu_hat=r$mu
      sigma_hat=r$sigma
      a_hat=r$alpha.vec
      b_hat=r$beta.vec
      
      #Evaluate variance
      info=information(a=a_hat, b=b_hat, g=g_hat, x=x, dat, z, mu=mu_hat, sigma=sigma_hat, q=r$post, anchor=anchor)
      varian=solve(info)
      eig=eigen(varian)
      eig=eig$values
    }
    
    r_info=rbind(r_info, info)
    r_mu=c(r_mu, mu_hat)
    r_sigma=c(r_sigma, sigma_hat)
    alpha=rbind(alpha, a_hat)
    beta=rbind(beta, b_hat)
    gamma=rbind(gamma, g_hat)

    var_all = diag(varian)
    sd_g=(var_all[(2*p+2*J+1):length(var_all)])^0.5
    transform_p=g_hat[-anchor]
    
    #Get coverages
    cover_c=c()
    gam_no_anchor=gamma.vec[-c(anchor)]
    for (j in 1:length(transform_p)){
      upp_c=transform_p[j]+qnorm(1-sig/2)*sd_g[j] ### assignment for sd_g ?
      low_c=transform_p[j]-qnorm(1-sig/2)*sd_g[j]
      if (gam_no_anchor[j] >=low_c & gam_no_anchor[j]<=upp_c){
        cover_c=c(cover_c, 1)
      }
      else{cover_c=c(cover_c, 0)}
    }
    cover=rbind(cover, cover_c)
    
    #To get p-values
    
    p_value_curr=c()
    for (j in 1:length(transform_p)){
      stat=transform_p[j]
      p_value_c=2*(1-pnorm(abs(stat), 0, sd_g[j]))
      p_value_curr=c(p_value_curr, p_value_c)
    }
    p_value=rbind(p_value, p_value_curr)
    
    #To get ROC curve
    TPR_c=c()
    FPR_c=c()
    for (s in sig_list){
      upper=c()
      lower=c()
      for (j in 1:length(transform_p)){
        upper=c(upper, transform_p[j]+qnorm(1-s/2)*sd_g[j])
        lower=c(lower, transform_p[j]-qnorm(1-s/2)*sd_g[j])
      }
      
      upper_p=upper[ind1_gam]
      lower_p=lower[ind1_gam]
      upper_n=upper[ind0_gam]
      lower_n=lower[ind0_gam]
      TPR_c=c(TPR_c, length(c(which(upper_p<0), which(lower_p>0)))/n_p_gam)
      FPR_c=c(FPR_c, length(c(which(upper_n<0), which(lower_n>0)))/n_n_gam)
    }
    TPR=rbind(TPR, TPR_c)
    FPR=rbind(FPR, FPR_c)
  }
  
  
  
  return(list(alpha=alpha, beta=beta, gamma=gamma, mu=r_mu, sigma=r_sigma,
              cover=cover, p_value=p_value, TPR=TPR, FPR=FPR,r_info=r_info))
}


#### Get performances e.g. FPR, TPR, AIC, BIC and etc for the LASSO method ####
compare_method=function(N, penalty, a, b, g, z, mu, sigma, w, k, seed, tol = 0.005){
  set.seed(seed)
  J=length(g)
  n0=N/2
  n1=N/2
  x.vec=c(rep(0, n0), rep(1, n1))
  x = matrix(x.vec)
  ind1=which(g!=0)
  ind0=which(g==0)
  n_p=length(ind1)
  n_n=length(ind0)
  BIC_L1=c()
  BIC_prop=c()
  AIC_L1=c()
  AIC_prop=c()
  gam_L1_opt=NULL
  gam_prop_opt=NULL
  gam_L1=NULL
  gam_prop=NULL
  TPR_L1=c()
  TPR_prop=c()
  FPR_L1=c()
  FPR_prop=c()
  opt_penalty_AIC=c()
  opt_penalty_BIC=c()
  opt_threshold_AIC=c()
  opt_threshold_BIC=c()
  for (q in 1:20){
  
    #Generate data
    mu.vec = c(mu)
    sigma.vec = c(sigma)
    mu.v=x*rep(1, N)%*%t(mu.vec)
    var.v=x*(rep(1, N)%*%t(sigma.vec^2))
    zero.vec=apply(x, 1, if_zero)
    var.v=var.v+zero.vec
    sigma.v=var.v^0.5
    th=rnorm(n=N,mean=mu.v, sd=sigma.v)
    temp=(th%*%t(rep(1, J)))*(rep(1, N)%*%t(a))+rep(1, N)%*%t(b)+x%*%t(g) #N by J
    prob=exp(temp)/(1+exp(temp)) #dim N by J
    dat = matrix(0, N, J);
    dat[] = rbinom(N*J, 1, prob);
    
    #Perform BIC selection of penalty
    quadrature=quad(a, b, g, x, dat, z, mu=mu, sigma=sigma, w)
    info=hessian(a, b, g, x=x, dat, z, mu=mu, sigma=sigma, q=quadrature)
    
    g_temp=NULL
    BIC_temp=c()
    AIC_temp=c()
    for (l in penalty){
      r_l1=EM_2PL_L1(runif(J), runif(J), runif(J), x=x, dat, z, mu=mu, sigma=sigma, w, lmd=l, tau=info, k=1, tol)
      g_l1=r_l1$gamma.vec
      g_temp=rbind(g_temp, g_l1)
      print(g_l1)
      r_l1_mll=EM_2PL_non_zero(a=r_l1$alpha.vec, b=r_l1$beta.vec, g=g_l1, x=x, dat, z, mu=r_l1$mu, sigma=r_l1$sigma, w, tol)
      k=length(which(g_l1!=0))+length(which(c(r_l1$alpha.vec, r_l1$beta.vec, r_l1$mu, r_l1$sigma)!=0))
      mll=r_l1_mll$MLL
      BIC_temp=c(BIC_temp, k*log(N)+2*mll)
      AIC_temp=c(AIC_temp, 2*k+2*mll)
    }
    gam_L1=rbind(gam_L1, g_temp)
    idx_opt_lmd=which.min(BIC_temp)
    opt_penalty_BIC=c(opt_penalty_BIC, penalty[idx_opt_lmd])
    opt_penalty_AIC=c(opt_penalty_AIC, penalty[which.min(AIC_temp)])
    BIC_L1=c(BIC_L1, min(BIC_temp))
    AIC_L1=c(AIC_L1,  AIC_temp[idx_opt_lmd])
    opt_g=g_temp[idx_opt_lmd, ]
    gam_L1_opt=rbind(gam_L1_opt, opt_g)
    g_hat_p= opt_g[ind1]
    g_hat_n= opt_g[ind0]
    
    TPR_L1=c(TPR_L1, length(which(g_hat_p!=0))/n_p)
    FPR_L1=c(FPR_L1, length(which(g_hat_n!=0))/n_n)
    
  }  
  
  return(list(TPR_L1=TPR_L1, FPR_L1=FPR_L1, AIC_L1=AIC_L1, BIC_L1=BIC_L1, gam_L1=gam_L1, gam_L1_opt=gam_L1_opt, opt_penalty_AIC=opt_penalty_AIC, opt_penalty_BIC=opt_penalty_BIC))
}



###Benjamini&Hochberg procedures

bh=function(p_value, fdr){#return a list of indices that are predicted to be positive
  J=length(p_value)
  ord=order(p_value)
  ordered_pvalue=p_value[ord]
  rank=1:J
  critical_value=c()
  for (i in 1:J){
    critical_value=c(critical_value, fdr*rank[i]/J)
  }
  

  if(all(critical_value < ordered_pvalue)){
    return(c(0))
  }
  else{
    i = max(which(critical_value >= ordered_pvalue))
    threshold=ordered_pvalue[i]
    print(which(p_value<=threshold))
    return(which(p_value<=threshold))
  }
}

FDR=function(pmat, gamma.vec, fdr){
  ind0=which(gamma.vec==0)
  J=length(gamma.vec)
  False_dis_rate=c()
  for (i in 1:dim(pmat)[1]){
    prd_ind=bh(p_value=pmat[i,], fdr=fdr)
    if (prd_ind[1]==0){
      False_dis_rate=c(False_dis_rate, 0)
    }
    
    else{
      wrong=length(intersect(prd_ind, ind0))
      False_dis_rate=c(False_dis_rate, round(wrong/length(prd_ind), 3))
    }
  }
  return(mean(False_dis_rate))
}



###Evaluate area under ROC curve ###
eva_tpr_fpr=function(gam, gamma.vec, rep){
  ind1=which(gamma.vec!=0)
  ind0=which(gamma.vec==0)
  n_p=length(ind1)
  n_n=length(ind0)
  TPR=NULL
  FPR=NULL
  for (i in 1:rep){
    tpr_c=c()
    fpr_c=c()
    for (j in 1:19){# 19 penalties
      ga=as.numeric(gam[(19*(i-1)+j), ])
      g_hat_p=ga[ind1]
      g_hat_n=ga[ind0]
      tpr_c=c(tpr_c, length(which(g_hat_p!=0))/n_p)
      fpr_c=c(fpr_c, length(which(g_hat_n!=0))/n_n)
    }
    TPR=rbind(TPR, tpr_c)
    FPR=rbind(FPR, fpr_c)
  }
  return(list(TPR=TPR, FPR=FPR))
}

eva_auc=function(tpr, fpr){
  tpr_c=apply(tpr, 2, mean)
  fpr_c=apply(fpr, 2, mean)
  tpr_c=c(tpr_c, 1)
  fpr_c=c(fpr_c, 1)
  id = order(fpr_c)
  auc=sum(diff(fpr_c[id])*rollmean(tpr_c[id],2))
  return(mean(auc))
}



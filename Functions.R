library(fastGHQuad)
library(quantreg)
library(MASS)
library(zoo)


gauss=gaussHermiteData(n=31)#We use 31 Gaussian quadrature
w=gauss$w
z=gauss$x


### Estimation ###
#Transform Gauss Hermite Quadrature z into th, th_k=sqrt(2)*z_k+mu for k=1,...,n for all i=1,...,N.
#Return matrix Theta of dim N by n.
tran=function(z, x, mu){
  n=length(z)
  N=length(x)
  Theta=sqrt(2)*(rep(1, N)%*%t(z))+(mu*x) #N by n 
  return(Theta)
}



evaluate_prod=function(v){
  L=length(v)
  return(exp(sum(v)))
}


#Evaluate the posterior quadrature values, return a matrix N by n.
#a, b, g are alpha, beta and gamma of length J. dat is of dim N by J.
#x is a length N vector.
#z, w are the Gauss Hermite quadratures and weights of length n, 
quad=function(a, b, g, x, dat, z, mu, w){ #z is quadrature of dim n
  n=length(z)
  N=dim(dat)[1]
  J=dim(dat)[2]
  Theta=tran(z, x, mu) #Transform Gauss Hermite Quadrature z into th, th_i=sqrt(2)*z_i+mu for i=1,...,n
  th=Theta[, 1]
  temp=(th%*%t(rep(1, J)))*(rep(1, N)%*%t(a))-rep(1, N)%*%t(b)+(x%*%t(rep(1, J)))*(rep(1, N)%*%t(g)) #N by J
  p=exp(temp)/(1+exp(temp)) #dim N by J
  q=apply(log((p^dat)*((1-p)^(1-dat))), 1, evaluate_prod) #length N vector
  for (k in 2:n){ #matrix evaluate
    th=Theta[, k]
    temp=(th%*%t(rep(1, J)))*(rep(1, N)%*%t(a))-rep(1, N)%*%t(b)+(x%*%t(rep(1, J)))*(rep(1, N)%*%t(g))
    p=exp(temp)/(1+exp(temp))
    q_c=apply(log((p^dat)*((1-p)^(1-dat))), 1, evaluate_prod)
    q=cbind(q, q_c)
  }
  Z=rowSums(q*(rep(1, N)%*%t(w)))
  return((q*(rep(1, N)%*%t(w)))/(Z%*%t(rep(1, n))))
}



negative_log_likelihood=function(a, b, g, x, dat, z, mu, mu_pres, w, q){ #q is quadrature
  N=length(x)
  n=length(z)
  J=dim(dat)[2]
  
  #Evaluate log-likelihood
  Theta=tran(z, x, mu_pres) #N by n
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


#Evaluate target function at phi_j=(a_j, b_j, g_j), i.e. negative log-likelihood at phi_j
#Y_j=Y[,j]
target_function=function(phi_j, x, Y_j, z,  mu_pres, q){
  N=length(x)
  n=length(z)
  n1=length(which(x==1))
  n0=length(which(x==0))
  a_j=phi_j[1]
  b_j=phi_j[2]
  g_j=phi_j[3]
  
  #Evaluate log-likelihood
  Theta=tran(z, x, mu_pres)
  t1=-0.5*(Theta-c(rep(0, n0), rep(mu_pres, n1))%*%t(rep(1, n)))^2 #normal term in log-likelihood  #remove prior term
  temp=Theta*a_j-b_j*(rep(1,N)%*%t(rep(1, n)))+(x%*%t(rep(1, n)))*g_j #N by n
  t2=Y_j%*%t(rep(1, n))*temp-log(1+exp(temp)) #The remaining term in log-likelihood
  log_likelihood=t1+t2 #N by n
  
  res=-sum(q*log_likelihood)
  return(res)
}

#Evaluate gradient with respect to phi
grad_phi=function(phi_j, x, Y_j, z, mu_pres, q){
  N=length(x)
  n=length(z)
  a_j=phi_j[1]
  b_j=phi_j[2]
  g_j=phi_j[3]
  
  #Evaluate gradient
  Theta=tran(z, x, mu_pres)
  temp=Theta*a_j-b_j*(rep(1,N)%*%t(rep(1, n)))+(x%*%t(rep(1, n)))*g_j #N by n
  p=exp(temp)/(1+exp(temp))
  da=sum(Theta*(Y_j%*%t(rep(1, n))-p)*q)/N
  db=sum((-Y_j%*%t(rep(1, n))+p)*q)/N
  dg=sum((x%*%t(rep(1, n)))*(Y_j%*%t(rep(1, n))-p)*q)/N
  return(-c(da, db, dg))
}


#Evaluate marginal log-likelihood
mml=function(a, b, g, x, dat, z, mu, w){ #q is quadrature
  N=length(x)
  n=length(z)
  J=dim(dat)[2]
  
  #Evaluate log-likelihood
  Theta=tran(z, x, mu) #N by n
  th=Theta[, 1]
  temp=(th%*%t(rep(1, J)))*(rep(1, N)%*%t(a))-rep(1, N)%*%t(b)+(x%*%t(rep(1, J)))*(rep(1, N)%*%t(g)) #N by J
  p=exp(temp)/(1+exp(temp))
  res=(p^dat)*((1-p)^(1-dat))*w[1]*(pi)^(-0.5)
  for (k in 2:n){
    th=Theta[, k]
    temp=(th%*%t(rep(1, J)))*(rep(1, N)%*%t(a))-rep(1, N)%*%t(b)+(x%*%t(rep(1, J)))*(rep(1, N)%*%t(g))
    p=exp(temp)/(1+exp(temp))
    res_r=(p^dat)*((1-p)^(1-dat))*w[k]*(pi)^(-0.5)
    res=res+res_r
  }
  res=-sum(log(res))
  return(res)
}

#EM algorithm for the proposed inference method
EM_2PL_inference <- function(a, b, g, x, dat, z, mu, w, tol = 0.001){
  J=length(a)
  ind1=which(x==1)
  n1=length(ind1)
  quadrature=quad(a, b, g, x, dat, z, mu, w)
  phi0=matrix(0, J, 3)
  phi1=cbind(a, b, g)
  mu0=0
  while(max(c(max(abs(phi0-phi1)), abs(mu0-mu)))>tol){#max
    phi0=phi1
    mu0=mu
    
    #E-step
    quadrature=quad(a=phi0[,1], b=phi0[,2], g=phi0[,3], x, dat, z, mu0, w)
    
    #M-step
    for(j in 1:J){#update phi_j one by one
      par_updates = optim(par=phi1[j, ], fn=target_function, gr=grad_phi, method = "L-BFGS-B", x=x, Y_j=dat[,j], z, mu_pres=mu0, q=quadrature, lower=-3, upper=3)
      phi1[j, ]=par_updates$par
      phi1[1, 3]=0
    }
    Theta=tran(z, x, mu0)
    mu=sum((Theta*quadrature)[ind1, ])/n1
    MLL = mml(a=phi1[,1], b=phi1[,2], g=phi1[,3], x, dat, z, mu, w)
    print(MLL) 
  }
  list(mu=mu, alpha.vec=phi1[,1], beta.vec = phi1[,2], gamma.vec = phi1[,3], post=quadrature);
}



#EM algorithm for likelihood ratio method that requires knowing anchor set of items
EM_2PL_lrt <- function(a, b, g, x, dat, z, mu, w, anchor, tol = 0.001){#anchor is a list of indices of item
  J=length(a)
  ind1=which(x==1)
  n1=length(ind1)
  quadrature=quad(a, b, g, x, dat, z, mu, w)
  phi0=matrix(0, J, 3)
  phi1=cbind(a, b, g)
  mu0=0
  while(max(c(max(abs(phi0-phi1)), abs(mu0-mu)))>tol){#max
    phi0=phi1
    mu0=mu
    #E-step
    quadrature=quad(a=phi0[,1], b=phi0[,2], g=phi0[,3], x, dat, z, mu0, w)
    
    #M-step
    for(j in 1:J){#update phi_j one by one
      par_updates = optim(par=phi1[j, ], fn=target_function, gr=grad_phi, method = "L-BFGS-B", x=x, Y_j=dat[,j], z, mu_pres=mu0, q=quadrature, lower=-3, upper=3)
      phi1[j, ]=par_updates$par
      #phi1[1, 3]=0
    }
    phi1[anchor, 3]=0
    Theta=tran(z, x, mu0)
    mu=sum((Theta*quadrature)[ind1, ])/n1
    MLL = mml(a=phi1[,1], b=phi1[,2], g=phi1[,3], x, dat, z, mu, w)
    print(MLL) 
  }
  list(mu=mu, alpha.vec=phi1[,1], beta.vec = phi1[,2], gamma.vec = phi1[,3], post=quadrature);
}




#Return the expected second derivative of individual i.
#Here x_i=x.vec[i], t_i=Theta[i, ], f_i=Posterior[i, ], y_i=Y[i, ]
#Return (3J+1) by (3J+1) matrix
first_term=function(a, b, g, x_i, y_i, t_i, mu, f_i){
  J=length(a)
  n=length(t_i)
  dmu=-x_i
  temp=a%*%t(t_i)-b%*%t(rep(1, n))+x_i*(g%*%t(rep(1, n))) #J by n
  p=exp(temp)/(1+exp(temp)) ##J by n
  v=p*(1-p) #variance matrix, J by n
  da=rowSums(-(rep(1, J)%*%t(t_i^2))*v*(rep(1, J)%*%t(f_i)))
  db=rowSums(-v*(rep(1, J)%*%t(f_i)))
  dg=rowSums(-(rep(1, J)%*%t(rep(x_i^2, n)))*v*(rep(1, J)%*%t(f_i)))
  dab=rowSums((rep(1, J)%*%t(t_i))*v*(rep(1, J)%*%t(f_i)))
  dag=rowSums(-(rep(1, J)%*%t(x_i*t_i))*v*(rep(1, J)%*%t(f_i)))
  dbg=rowSums((rep(1, J)%*%t(rep(x_i, n)))*v*(rep(1, J)%*%t(f_i)))
  #Fill into the second derivative matrix
  res=matrix(0, 3*J, 3*J)
  diag(res)=c(da, db, dg)
  res[1:J, (J+1):(2*J)]=diag(dab)
  res[1:J, (2*J+1):(3*J)]=diag(dag)
  res[(J+1):(2*J), 1:J]=diag(dab)
  res[(J+1):(2*J), (2*J+1):(3*J)]=diag(dbg)
  res[(2*J+1):(3*J), 1:J]=diag(dag)
  res[(2*J+1):(3*J), (J+1):(2*J)]=diag(dbg)
  res=rbind(rep(0, 3*J), res)
  res=cbind(c(dmu, rep(0, 3*J)), res)
  return(res)
}


#Return the expected (dl/dphi)(dl/dphi)^t  of individual i.
#Here x_i=x.vec[i], t_i=Theta[i, ], f_i=Posterior[i, ], y_i=dat[i, ]
#Return (3J+1) by (3J+1) matrix
second_term=function(a, b, g, x_i, y_i, t_i, mu, f_i){
  J=length(a)
  n=length(t_i)
  temp=a%*%t(t_i)-b%*%t(rep(1, n))+x_i*(g%*%t(rep(1, n))) #J by n
  p=exp(temp)/(1+exp(temp)) ##J by n
  da1=(rep(1, J)%*%t(t_i))*(y_i%*%t(rep(1, n))-p) #dl/da  J by n
  db1=(-y_i%*%t(rep(1, n))+p) #dl/db  J by n
  dg1=x_i*(y_i%*%t(rep(1, n))-p) #dl/dg  J by n
  dmu2=sum(x_i*((t_i-mu)^2)*f_i) #store E[(dl/dmu)^2]
  dmua=rowSums((rep(1, J)%*%t(x_i*(t_i-mu)))*da1*(rep(1, J)%*%t(f_i))) #store E[(dl/dmu)(dl/da)^t]
  dmub=rowSums((rep(1, J)%*%t(x_i*(t_i-mu)))*db1*(rep(1, J)%*%t(f_i))) #store E[(dl/dmu)(dl/db)^t]
  dmug=rowSums((rep(1, J)%*%t(x_i*(t_i-mu)))*dg1*(rep(1, J)%*%t(f_i))) #store E[(dl/dmu)(dl/dg)^t]
  da2=matrix(0, J, J) #store E[(dl/da)(dl/da)^t]
  db2=matrix(0, J, J)#store E[(dl/db)(dl/db)^t]
  dg2=matrix(0, J, J)#store E[(dl/dg)(dl/dg)^t]
  dab=matrix(0, J, J)#store E[(dl/da)(dl/db)^t]
  dag=matrix(0, J, J)#store E[(dl/da)(dl/dg)^t]
  dbg=matrix(0, J, J)#store E[(dl/db)(dl/dg)^t]
  
  for (k in 1:n){#taking the expectation
    da2=da2+da1[, k]%*%t(da1[, k])*f_i[k]
    db2=db2+db1[, k]%*%t(db1[, k])*f_i[k]
    dg2=dg2+dg1[, k]%*%t(dg1[, k])*f_i[k]
    dab=dab+da1[, k]%*%t(db1[, k])*f_i[k]
    dag=dag+da1[, k]%*%t(dg1[, k])*f_i[k]
    dbg=dbg+db1[, k]%*%t(dg1[, k])*f_i[k]
  }
  result=matrix(0, 3*J, 3*J)
  result[1:J, 1:J]=da2
  result[1:J, (J+1):(2*J)]=dab
  result[(J+1):(2*J), 1:J]=t(dab)
  result[1:J, (2*J+1):(3*J)]=dag
  result[(2*J+1):(3*J), 1:J]=t(dag)
  result[(J+1):(2*J), (J+1):(2*J)]=db2
  result[(J+1):(2*J), (2*J+1):(3*J)]=dbg
  result[(2*J+1):(3*J),(J+1):(2*J)]=t(dbg)
  result[(2*J+1):(3*J),(2*J+1):(3*J)]=dg2
  result=rbind(c(dmua, dmub, dmug), result)
  result=cbind(c(dmu2, dmua, dmub, dmug), result)
  return(result)
}


#Return the expected (dl/dphi)(dl/dphi)^t  of individual i.
#Here x_i=x.vec[i], t_i=Theta[i, ], f_i=Posterior[i, ], y_i=dat[i, ]
#Return (3J+1) by (3J+1) matrix
third_term=function(a, b, g, x_i, y_i, t_i, mu, f_i){
  J=length(a)
  n=length(t_i)
  temp=a%*%t(t_i)-b%*%t(rep(1, n))+x_i*(g%*%t(rep(1, n))) #J by n
  p=exp(temp)/(1+exp(temp)) ##J by n
  dmu=sum(x_i*(t_i-mu)*f_i) #E[dl/dmu]
  da=rowSums((rep(1, J)%*%t(t_i))*(y_i%*%t(rep(1, n))-p)*(rep(1, J)%*%t(f_i)))#E[dl/da]
  db=rowSums((-y_i%*%t(rep(1, n))+p)*(rep(1, J)%*%t(f_i)))#E[dl/db]
  dg=rowSums(x_i*(y_i%*%t(rep(1, n))-p)*(rep(1, J)%*%t(f_i)))#E[dl/dg]
  dl1=c(dmu, da, db, dg)
  #dl1=c(da, db, dg)
  return(dl1%*%t(dl1))
}

#return information matrix
#q is the posterior, N by n
information=function(a, b, g, x, dat, z, mu, q, anchor){
  N=length(x)
  J=length(a)
  Theta=tran(z, x, mu)
  res=matrix(0, 3*J+1, 3*J+1)
  for (i in 1:N){
    x_i=x[i]
    y_i=dat[i, ]
    t_i=Theta[i,]
    f_i=q[i, ]
    res=res+first_term(a, b, g, x_i, y_i, t_i, mu, f_i) + second_term(a, b, g, x_i, y_i, t_i, mu, f_i) - third_term(a, b, g, x_i, y_i, t_i, mu, f_i)
  }
  res=res[,-(2*J+anchor+1)]
  res=res[-(2*J+anchor+1), ]
  return(-res)
}



###Estimation for the LASSO method###
##Evaluate Hessian for better initialization of learning rate
hessian=function(a, b, g, x, dat, z, mu, q){
  N=length(x)
  J=length(a)
  Theta=tran(z, x, mu)
  res=matrix(0, 3*J+1, 3*J+1)
  #res=matrix(0, 3*J, 3*J)
  for (i in 1:N){
    x_i=x[i]
    y_i=dat[i, ]
    t_i=Theta[i,]
    f_i=q[i, ]
    res=res+first_term(a, b, g, x_i, y_i, t_i, mu, f_i)
  }
  
  return(-diag(res)/N)
}


#Training LASSO method
EM_2PL_L1 <- function(a, b, g, x, dat, z, mu, w, lmd, tau, k, tol){
  J=length(a)
  ind1=which(x==1)
  n1=length(ind1)
  phi0=matrix(0, J, 3)
  phi1=cbind(a, b, g)
  mu0=0
  while(max(c(max(abs(phi0-phi1)), abs(mu0-mu)))>tol){
    phi0=phi1
    mu0=mu
    #E-step
    quadrature=quad(a=phi0[,1], b=phi0[,2], g=phi0[,3], x, dat, z, mu0, w)
    
    #M-step
    for(j in 1:J){#update phi_j one by one
      a_c=phi1[j, 1]
      b_c=phi1[j, 2]
      g_c=phi1[j, 3]
      for(i in 1:k){
        grad=grad_phi(c(a_c, b_c, g_c), x, Y_j=dat[, j], z, mu_pres=mu0, q=quadrature)
        a_c= a_c - tau[j]*grad[1]
        b_c= b_c - tau[J+j]*grad[2]
        g_c= g_c - tau[2*J+j]*grad[3]
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
    
    Theta=tran(z, x, mu0)
    mu=sum((Theta*quadrature)[ind1, ])/n1
    #MLL = mml_with_penalty(a=phi1[,1], b=phi1[,2], g=phi1[,3], x, dat, z, mu, w)
    #print(MLL) 
    #print(phi1[,3])
  }
  list(mu=mu, alpha.vec=phi1[,1], beta.vec = phi1[,2], gamma.vec = phi1[,3], post=quadrature);
}

##de-bias the non-zero entries obtained from LASSO method so as to compute the BIC
EM_2PL_non_zero <- function(a, b, g, x, dat, z, mu, w, tol){
  J=length(a)
  ind1=which(x==1)
  n1=length(ind1)
  quadrature=quad(a, b, g, x, dat, z, mu, w)
  idx0=which(g==0)
  MLL = negative_log_likelihood(a, b, g, x, dat, z, mu=mu, mu_pres=mu, w, q=quadrature)
  #print(MLL)
  phi0=matrix(0, J, 3)
  phi1=cbind(a, b, g)
  mu0=0
  while(max(c(max(abs(phi0-phi1)), abs(mu0-mu)))>tol){#max
    phi0=phi1
    mu0=mu
    #E-step
    quadrature=quad(a=phi0[,1], b=phi0[,2], g=phi0[,3], x, dat, z, mu0, w)
    
    #M-step
    for(j in 1:J){#update phi_j one by one
      par_updates = optim(par=phi1[j, ], fn=target_function, gr=grad_phi, method = "L-BFGS-B", x=x, Y_j=dat[,j], z, mu_pres=mu0, q=quadrature, lower=-3, upper=3)
      phi1[j, ]=par_updates$par
      #phi1[1, 3]=0
    }
    phi1[idx0, 3]=0
    Theta=tran(z, x, mu0)
    mu=sum((Theta*quadrature)[ind1, ])/n1
    MLL = mml(a=phi1[,1], b=phi1[,2], g=phi1[,3], x, dat, z, mu, w)
    #print(MLL) 
  }
  list(MLL=MLL);
}



### Coverage rate ####

#For LRT method
coverage_LRT=function(N, sig, gamma.vec, beta.vec, alpha.vec, anchor, sig_list){#sig_list is used to construct ROC
  J=length(gamma.vec)
  n0=N/2
  n1=N/2
  x.vec=c(rep(0, n0), rep(1, n1))
  ind1=which(gamma.vec!=0) 
  ind0=which(gamma.vec==0)
  n_p=length(ind1)
  n_n=length(ind0)
  ind0_c=ind0[-anchor]
  ind0_c=ind0_c-length(anchor)
  ind1_c=ind1-length(anchor)
  nn_c=n_n-length(anchor)
  n_sig=length(sig_list)
  set.seed(5)
  cover=NULL
  p_value=NULL
  alpha=NULL
  beta=NULL
  gamma=NULL
  TPR=NULL
  FPR=NULL
  r_mu=c()
  for (l in 1:100){
    #Generate data
    theta.vec=c(rnorm(n0, mean = 0, sd=1), rnorm(n1, mean = 0.5, sd=1))
    temp = (theta.vec %*% t(rep(1, J))*(rep(1, N)%*%t(alpha.vec))) - rep(1, N) %*% t(beta.vec) + (x.vec%*%t(rep(1, J)))*(rep(1, N)%*%t(gamma.vec))
    prob = 1/(1+exp(-temp));
    dat = matrix(0, N, J);
    dat[] = rbinom(N*J, 1, prob); 
    
    #Train model
    r=EM_2PL(a=alpha.vec, b=beta.vec, g=gamma.vec, x=x.vec, dat, z, mu=0.5, w, anchor = anchor)
    g_hat=r$gamma.vec
    mu_hat=r$mu
    a_hat=r$alpha.vec
    b_hat=r$beta.vec
    
    #Evaluate variance
    info=information(a=r$alpha.vec, b=r$beta.vec, g=r$gamma.vec, x=x.vec, dat, z, mu=r$mu, q=r$post, anchor=anchor)
    varian=solve(info)
    eig=eigen(varian)
    eig=eig$values
    
    if (any(eig<0)){
      theta.vec=c(rnorm(n0, mean = 0, sd=1), rnorm(n1, mean = 0.5, sd=1))
      temp = (theta.vec %*% t(rep(1, J))*(rep(1, N)%*%t(alpha.vec))) - rep(1, N) %*% t(beta.vec) + (x.vec%*%t(rep(1, J)))*(rep(1, N)%*%t(gamma.vec))
      prob = 1/(1+exp(-temp));
      dat = matrix(0, N, J);
      dat[] = rbinom(N*J, 1, prob); 
      
      #Train model
      r=EM_2PL(a=alpha.vec, b=beta.vec, g=gamma.vec, x=x.vec, dat, z, mu=0.5, w, anchor = anchor)
      g_hat=r$gamma.vec
      mu_hat=r$mu
      a_hat=r$alpha.vec
      b_hat=r$beta.vec
      #Evaluate variance
      info=information(a=r$alpha.vec, b=r$beta.vec, g=r$gamma.vec, x=x.vec, dat, z, mu=r$mu, q=r$post,anchor=anchor)
      varian=solve(info)
      eig=eigen(varian)
      eig=eig$values
    }
    r_mu=c(r_mu, mu_hat)
    alpha=rbind(alpha, a_hat)
    beta=rbind(beta, b_hat)
    gamma=rbind(gamma, g_hat)
    
    #simulate new gammas to get p-values and CIs
    transform_p=g_hat[-anchor]
    sim=mvrnorm(n=10000, mu=c(mu_hat, a_hat, b_hat, transform_p), Sigma=varian)
    sim_gam=sim[, (2*J+2):(dim(sim)[2])]
    
    #Get coverages
    cover_c=c()
    gam_no_anchor=gamma.vec[-c(anchor)]
    for (j in 1:dim(sim_gam)[2]){
      upp_c=quantile(sim_gam[, j], 1-sig/2)
      low_c=quantile(sim_gam[, j], sig/2)
      if (gam_no_anchor[j]>=low_c & gam_no_anchor[j]<=upp_c){
        cover_c=c(cover_c, 1)
      }
      else{cover_c=c(cover_c, 0)}
    }
    cover=rbind(cover, cover_c)
    
    #To get p-values
    ga_trans=sim_gam-rep(1, 10000)%*%t(transform_p)
    p_value_curr=c()
    for (j in 1:dim(sim_gam)[2]){
      stat=transform_p[j]
      dist=as.vector(ga_trans[, j])
      p_value_c=length(c(which(dist>abs(stat)), which(dist<   -abs(stat))))/10000
      p_value_curr=c(p_value_curr, p_value_c)
    }
    p_value=rbind(p_value, p_value_curr)
    
    #To get ROC curve
    TPR_c=c()
    FPR_c=c()
    for (s in sig_list){
      upper=c()
      lower=c()
      for (j in 1:dim(sim_gam)[2]){
        upper=c(upper, quantile(sim_gam[, j], 1-s/2))
        lower=c(lower, quantile(sim_gam[, j], s/2))
      }
      
      upper_p=upper[ind1_c]
      lower_p=lower[ind1_c]
      upper_n=upper[ind0_c]
      lower_n=lower[ind0_c]
      TPR_c=c(TPR_c, length(c(which(upper_p<0), which(lower_p>0)))/n_p)
      FPR_c=c(FPR_c, length(c(which(upper_n<0), which(lower_n>0)))/nn_c)
    }
    TPR=rbind(TPR, TPR_c)
    FPR=rbind(FPR, FPR_c)
  }
  
  return(list(alpha=alpha, beta=beta, gamma=gamma, mu=r_mu, TPR=TPR, FPR=FPR, cover=cover, p_value=p_value))
}


### For proposed inference method ###

coverage_proposed_inference=function(N, sig, gamma.vec, beta.vec, alpha.vec, sig_list){
  J=length(gamma.vec)
  n0=N/2
  n1=N/2
  x.vec=c(rep(0, n0), rep(1, n1))
  ind1=which(gamma.vec!=0)
  ind0=which(gamma.vec==0)
  n_p=length(ind1)
  n_n=length(ind0)
  n_sig=length(sig_list)
  
  set.seed(5)
  cover=NULL
  p_value=NULL
  alpha=NULL
  beta=NULL
  gamma=NULL
  TPR=NULL
  FPR=NULL
  r_mu=c()
  for (l in 1:100){
    #Generate data
    theta.vec=c(rnorm(n0, mean = 0, sd=1), rnorm(n1, mean = 0.5, sd=1))
    temp = (theta.vec %*% t(rep(1, J))*(rep(1, N)%*%t(alpha.vec))) - rep(1, N) %*% t(beta.vec) + (x.vec%*%t(rep(1, J)))*(rep(1, N)%*%t(gamma.vec))
    prob = 1/(1+exp(-temp));
    dat = matrix(0, N, J);
    dat[] = rbinom(N*J, 1, prob); 
    
    #Train model
    r=EM_2PL(a=alpha.vec, b=beta.vec, g=gamma.vec, x=x.vec, dat, z, mu=0.5, w)
    g_hat=r$gamma.vec
    mu_hat=r$mu
    a_hat=r$alpha.vec
    b_hat=r$beta.vec
    #Evaluate variance
    info=information(a=r$alpha.vec, b=r$beta.vec, g=r$gamma.vec, x=x.vec, dat, z, mu=r$mu, q=r$post)
    varian=solve(info)
    eig=eigen(varian)
    eig=eig$values
    
    while (any(eig<0)){
      theta.vec=c(rnorm(n0, mean = 0, sd=1), rnorm(n1, mean = 0.5, sd=1))
      temp = (theta.vec %*% t(rep(1, J))*(rep(1, N)%*%t(alpha.vec))) - rep(1, N) %*% t(beta.vec) + (x.vec%*%t(rep(1, J)))*(rep(1, N)%*%t(gamma.vec))
      prob = 1/(1+exp(-temp));
      dat = matrix(0, N, J);
      dat[] = rbinom(N*J, 1, prob); 
      
      #Train model
      r=EM_2PL(a=alpha.vec, b=beta.vec, g=gamma.vec, x=x.vec, dat, z, mu=0.5, w)
      g_hat=r$gamma.vec
      mu_hat=r$mu
      a_hat=r$alpha.vec
      b_hat=r$beta.vec
      #Evaluate variance
      info=information(a=r$alpha.vec, b=r$beta.vec, g=r$gamma.vec, x=x.vec, dat, z, mu=r$mu, q=r$post)
      varian=solve(info)
      eig=eigen(varian)
      eig=eig$values
      
    }
    r_mu=c(r_mu, mu_hat)
    alpha=rbind(alpha, a_hat)
    beta=rbind(beta, b_hat)
    gamma=rbind(gamma, g_hat)
    
    #simulate new gammas to get p-values and CIs
    rqfit = rq(r$gamma.vec ~ -1 + r$alpha.vec)
    transform_p=g_hat-rqfit$coefficients[1]*a_hat
    sim=mvrnorm(n=10000, mu=c(mu_hat, a_hat, b_hat, g_hat[2:J]), Sigma=varian)
    sim_a=sim[, 2:(J+1)]
    sim_gam=cbind(rep(0, 10000), sim[, (2*J+2):(3*J)])
    
    gamm_transform=NULL
    for (k in 1:10000){
      rqfit = rq(sim_gam[k, ] ~ -1+ sim_a[k,])
      gamm_transform_c=sim_gam[k, ]-rqfit$coefficients[1]*sim_a[k,]
      gamm_transform=cbind(gamm_transform, gamm_transform_c)
    }
    
    #To get p-values
    ga_trans=gamm_transform-transform_p%*%t(rep(1, 10000))
    p_value_curr=c()
    for (j in 1:J){
      stat=transform_p[j]
      dist=as.vector(ga_trans[j, ])
      p_value_c=length(c(which(dist>abs(stat)), which(dist<   -abs(stat))))/10000
      p_value_curr=c(p_value_curr, p_value_c)
    }
    p_value=rbind(p_value, p_value_curr)
    
    
    cover_c=c()
    for (j in 1:J){
      upp_c=quantile(gamm_transform[j, ], 1-sig/2)
      low_c=quantile(gamm_transform[j, ], sig/2)
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
      for (j in 1:J){
        upper=c(upper, quantile(gamm_transform[j, ], 1-s/2))
        lower=c(lower, quantile(gamm_transform[j, ], s/2))
      }
      
      upper_p=upper[ind1]
      lower_p=lower[ind1]
      upper_n=upper[ind0]
      lower_n=lower[ind0]
      TPR_c=c(TPR_c, length(c(which(upper_p<0), which(lower_p>0)))/n_p)
      FPR_c=c(FPR_c, length(c(which(upper_n<0), which(lower_n>0)))/n_n)
    }
    TPR=rbind(TPR, TPR_c)
    FPR=rbind(FPR, FPR_c)
    
    
  }
  
  
  return(list(alpha=alpha, beta=beta, gamma=gamma, mu=r_mu, TPR=TPR, FPR=FPR, cover=cover, p_value=p_value))
}



#### Compare performances e.g. FPR, TPR, AIC, BIC and etc for the proposed hard-thresholding method and the LASSO method ####
compare_method=function(N, penalty, thresholds, a, b, g, z, mu, w, k, tol = 0.001){
  J=length(g)
  n0=N/2
  n1=N/2
  x.vec=c(rep(0, n0), rep(1, n1))
  ind1=which(g!=0)
  ind0=which(g==0)
  n_p=length(ind1)
  n_n=length(ind0)
  set.seed(5)
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
  opt_penalty=c()
  opt_threshold=c()
  for (q in 1:100){
    
    #Generate data
    theta.vec=c(rnorm(n0, mean = 0, sd=1), rnorm(n1, mean = 0.5, sd=1))
    temp = (theta.vec %*% t(rep(1, J))*(rep(1, N)%*%t(a))) - rep(1, N) %*% t(b) + (x.vec%*%t(rep(1, J)))*(rep(1, N)%*%t(g))
    prob = 1/(1+exp(-temp));
    dat = matrix(0, N, J);
    dat[] = rbinom(N*J, 1, prob); 
    
    
    #Perform BIC selection of penalty
    quadrature=quad(a, b, g, x.vec, dat, z, mu=0.5, w)
    info=hessian(a, b, g, x=x.vec, dat, z, mu=0.5, q=quadrature)
    g_temp=NULL
    BIC_temp=c()
    AIC_temp=c()
    for (l in penalty){
      r_l1=EM_2PL_L1(runif(J), runif(J), runif(J), x=x.vec, dat, z, mu=0.5, w, lmd=l, tau=info, k=1, tol)
      g_l1=r_l1$gamma.vec
      g_temp=rbind(g_temp, g_l1)
      r_l1_mll=EM_2PL_non_zero(a=r_l1$alpha.vec, b=r_l1$beta.vec, g=g_l1, x=x.vec, dat, z, mu=r_l1$mu, w, tol)
      k=length(which(g_l1!=0))+length(which(c(r_l1$alpha.vec, r_l1$beta.vec, r_l1$mu)!=0))
      mll=r_l1_mll$MLL
      BIC_temp=c(BIC_temp, k*log(N)+2*mll)
      AIC_temp=c(AIC_temp, 2*k+2*mll)
    }
    gam_L1=rbind(gam_L1, g_temp)
    idx_opt_lmd=which.min(BIC_temp)
    opt_penalty=c(opt_penalty, penalty[idx_opt_lmd])
    BIC_L1=c(BIC_L1, min(BIC_temp))
    AIC_L1=c(AIC_L1,  AIC_temp[idx_opt_lmd])
    opt_g=g_temp[idx_opt_lmd, ]
    gam_L1_opt=rbind(gam_L1_opt, opt_g)
    g_hat_p= opt_g[ind1]
    g_hat_n= opt_g[ind0]
    
    TPR_L1=c(TPR_L1, length(which(g_hat_p!=0))/n_p)
    FPR_L1=c(FPR_L1, length(which(g_hat_n!=0))/n_n)
    
    
    
    #Proposed method
    r=EM_2PL(runif(J), runif(J), runif(J), x=x.vec, dat, z, mu=0.5, w, tol=tol)
    g_hat=r$gamma.vec
    mu_hat=r$mu
    a_hat=r$alpha.vec
    b_hat=r$beta.vec
    
    rqfit <- rq(g_hat ~ -1+a_hat)
    gamm_transform=g_hat-rqfit$coefficients[1]*a_hat
    
    #To remember the original index of each element
    g_temp=NULL
    BIC_temp=c()
    AIC_temp=c()
    for (l in thresholds){
      gam_t=gamm_transform
      gam_t[abs(gam_t)<l]=0
      g_temp=rbind(g_temp, gam_t)
      prop_refit=EM_2PL_non_zero(a=a_hat, b=b_hat, g=gam_t, x=x.vec, dat, z, mu=mu_hat+rqfit$coefficients[1], w, tol)
      k=length(which(gam_t!=0)) + length(which(c(a_hat, b_hat, mu_hat+rqfit$coefficients[1])!=0))
      mll=prop_refit$MLL
      BIC_temp=c(BIC_temp, k*log(N)+2*mll)
      AIC_temp=c(AIC_temp, 2*k+2*mll)
    }
    gam_prop=rbind(gam_prop, g_temp)
    idx_opt_threshold=which.min(BIC_temp)
    opt_threshold=c(opt_threshold, thresholds[idx_opt_threshold])
    BIC_prop=c(BIC_prop, min(BIC_temp))
    AIC_prop=c(AIC_prop,  AIC_temp[idx_opt_threshold])
    gam_final=g_temp[idx_opt_threshold, ]
    gam_prop_opt=rbind(gam_prop_opt, gam_final)
    g_hat_p=gam_final[ind1]
    g_hat_n=gam_final[ind0]
    TPR_prop=c(TPR_prop, length(which(g_hat_p!=0))/n_p)
    FPR_prop=c(FPR_prop, length(which(g_hat_n!=0))/n_n)
    
    
  }  
  
  return(list(TPR_L1=TPR_L1, FPR_L1=FPR_L1, AIC_L1=AIC_L1, BIC_L1=BIC_L1, gam_L1=gam_L1, gam_L1_opt=gam_L1_opt, opt_penalty=opt_penalty, 
              TPR_prop=TPR_prop, FPR_prop=FPR_prop, AIC_prop=AIC_prop, BIC_prop=BIC_prop, gam_prop=gam_prop, gam_prop_opt=gam_prop_opt, opt_threshold=opt_threshold))
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
  p_v=ordered_pvalue[1]
  i=1
  while (p_v<critical_value[i]){
    i=i+1
    p_v=ordered_pvalue[i]
  }
  if(i==1){
    return(c(0))
  }
  else{
    threshold=ordered_pvalue[1, i-1]
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
    for (j in 1:20){
      ga=as.numeric(gam[(20*(i-1)+j), ])
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

















































































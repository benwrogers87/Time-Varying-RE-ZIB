#Change so that random intercept is centered at fixed intercept rather than mean 0
#This way intercept has closed form posterior, given all else
#
#pseudo data approach



#rm(list = ls(all.names = TRUE))


library(MCMCpack) #For Wishart sampling
library(mvtnorm) #Multivariate normal
library(parallel)

expit<-function(x){exp(x)/(1+exp(x))}

log_bern<-function(x,prob){
  x*log(prob) + (1-x)*log(1-prob)
}


log_bin <- function(x, size, prob){
  lchoose(size,x) + x*log(prob) + (size-x)*log(1-prob) - log(1-size*log(1-prob))
}


log_norm_pdf<-function(x, mu, sig){
  
  -.5*t(x-mu)%*%solve(sig)%*%(x-mu) - .5*log(det(sig))
  
}




#######################################################################



# data_aug <- function(y, X, beta1, beta2, gamma1, gamma2){
#   clust <- rep(1,length(y))
#   V <- y
# 
#   pie <- expit(X%*%beta1 + rep(gamma1, each=4))
# 
#   clust <-   rbinom(length(y),1, pie)
# 
#   clust[y>0] <- 1
# 
#   theta <- expit(gamma2)
# 
#   V[clust==0] <- rbinom(sum(na.omit(clust)==0), 90, theta[clust==0])
#   V[is.na(y)] <- rbinom(sum(is.na(y)), 90, theta[is.na(y)])
# 
#   return(list(V=V, clust=clust))
# }



samp.clust <- function(alpha1, alpha2, gamma1, gamma2, psi){
  p <- rep(NA, length(y))
  theta_0 <- rep(NA, length(y))
  clust <- rep(NA, length(y))
  V <- rep(NA, length(y))
    linmod2 <- X%*%alpha2 + gamma2[ll] + psi*gamma1[ll]
    linmod1 <- X%*%alpha1 + gamma1[ll]
    p <- expit(linmod1)
    theta_0 <- (1-expit(linmod2))^(90)
    prob <- theta_0/(p + theta_0)
    clust <- suppressWarnings(rbinom(length(prob), 1, prob))
    clust[y > 0] <- 1
    V <- y*clust
    #zclust <- which(clust[,k]==0)
    
    #V[zclust,k] <- rpois(length(zclust),lambda[zclust,k])
  return(list(clust=clust, V=V))
}
#data_aug(y, X, beta1, beta2, gamma1, gamma2)

#Binomial Regression model



log_lik_zero <- function(clust, alpha1, gamma1){
  
  linmod <- X%*%alpha1 + gamma1[ll]
  sum(clust*linmod - log(1+exp(linmod)), na.rm=TRUE)
}



log_lik_hurdle <- function(V, alpha2, psi, gamma1, gamma2){
  linmod <- X%*%alpha2 + gamma2[ll] + psi*gamma1[ll]
  sum(V*linmod - K*log(1+exp(linmod)), na.rm=TRUE)
} 

#log_lik_zero(alpha1, gamma1)

log_post_alpha1<-function(clust, alpha1, gamma1){
  
  if((sum(abs(alpha1)>10)>0) | (sum(is.na(alpha1))>0)){
    return(-Inf)
  }else{
    like <- log_lik_zero(clust, alpha1, gamma1)
    prior <- sum(dnorm(alpha1,prior.mn1,prior.sd1, log=TRUE))
    
    post <- sum(like, na.rm=TRUE) +sum(prior)
    if(is.nan(post)){
      return(-Inf)
    }
    else{
      return(post)
    }
  }
}


gradient_alpha1 <- function(clust, alpha1, gamma1){
  d <- length(alpha1)
  diff <- rep(NA, d)
  
  q <- (1/prior.sd1^2)*diag(d)
    
  expon <- exp(X%*%alpha1 + gamma1[ll])
    
  lik_dot <- clust - expon/(1+expon)
    
  for(m in 1:d){
      
    diff[m] <- sum(X[,m]*lik_dot, na.rm=TRUE)
  }
    
  diff <- as.numeric(diff) - q%*%alpha1
    
  #diff[abs(alpha1)>10] <- 0
  
  return(diff)
}


# gradient_alpha1_numerical <- function(alpha1){
#   N <- length(alpha1)
#   diff <- rep(NA, N)
#   e <- .00001
#     for(i in 1:N){
#       th_hi <- alpha1
#       th_lo <- alpha1
#       th_hi[i] <- alpha1[i] + e
#       th_lo[i] <- alpha1[i] - e
#       #((((1/deltat)*cbind(Z,Zt)%*%alpha1 + sigma1*th_hi[ll] + log(deltat)))-((1/deltat)*cbind(Z,Zt)%*%alpha1 + sigma1*th_lo[ll] + log(deltat)))/(2*e)
#       diff[i]<- (log_post_alpha1(th_hi) - log_post_alpha1(th_lo))/(2*e)
#   }
#   return (diff)
# }

# diff1_a1 <- gradient_alpha1(alpha1)
# diff_a1_num <- gradient_alpha1_numerical(alpha1)
# diff1_a1
# diff_a1_num
#gradient_alpha1(alpha1)



gradient2_alpha1 <- function(alpha1, gamma1){
  d <- length(alpha1)
  diff <- rep(NA, d)
  
  q <- (1/prior.sd1^2)
  
  expon <- exp(X%*%alpha1 + gamma1[ll])
  
  lik_dot <- - expon/((1+expon)^2)
  
  for(m in 1:d){
    
    diff[m] <- sum(X[,m]*lik_dot, na.rm=TRUE)
  }
  
  diff <- as.numeric(diff) - q
  
  #diff[abs(alpha1)>10] <- 0
  
  return(diff)
}



# gradient2_alpha1_numerical <- function(alpha1){
#   N <- length(alpha1)
#   diff <- rep(NA, N)
#   e <- .00001
#   for(i in 1:N){
#     th_hi <- alpha1
#     th_lo <- alpha1
#     th_hi[i] <- alpha1[i] + e
#     th_lo[i] <- alpha1[i] - e
#     #((((1/deltat)*cbind(Z,Zt)%*%alpha1 + sigma1*th_hi[ll] + log(deltat)))-((1/deltat)*cbind(Z,Zt)%*%alpha1 + sigma1*th_lo[ll] + log(deltat)))/(2*e)
#     diff[i]<- (gradient_alpha1_numerical(th_hi)[i] - gradient_alpha1_numerical(th_lo)[i])/(2*e)
#   }
#   return (diff)
# }

# diff2_a1 <- gradient2_alpha1(alpha1, gamma1)
# diff2_a1_num <- gradient2_alpha1_numerical(alpha1)
# diff2_a1
# diff2_a1_num



alpha1_norm_approx <- function(clust, alpha1, gamma1){
  alp1 <- alpha1
  
  
  for(step in 1:steps.a1){
    grad1 <- gradient_alpha1(clust, alp1, gamma1)
    grad2 <- gradient2_alpha1(alp1, gamma1)
    
    Q_approx <- - diag(grad2) 
    mode_approx <- c(alp1) - diag(1/grad2)%*%c(grad1)
    alp1 <- mode_approx
  }
  return(list(mode_approx=c(mode_approx), Q_approx=Q_approx))
}

#Good starting value for alpha1:
#alpha1 <- c(  0.50338437, -0.81447167, -1.01128499, -1.48048324, -0.08876805, 0.07802240, -0.01134941)
# steps.a1 <- 1
# alpha1_norm_approx(alpha1)





alpha1_sample_norm <- function(clust, alpha1, gamma1){
  a1_approx <- alpha1_norm_approx(clust, alpha1, gamma1)
  

  alpha1_star <- as.numeric(rmvnorm(1, a1_approx$mode_approx, solve(a1_approx$Q_approx)))

  
  a1_approx_star <- alpha1_norm_approx(clust, alpha1_star, gamma1)
  
  log_post_old <-  log_post_alpha1(clust, alpha1, gamma1)
  log_post_star <- log_post_alpha1(clust, alpha1_star, gamma1)
  
   J_old <- sum(dnorm(alpha1, a1_approx_star$mode_approx, sqrt(diag(1/a1_approx_star$Q_approx)), log=TRUE))
   J_star <- sum(dnorm(alpha1_star, a1_approx$mode_approx, sqrt(diag(1/a1_approx$Q_approx)), log=TRUE))
  
  r <- exp(log_post_star + sum(J_old) - log_post_old - sum(J_star))
  if(!is.finite(r)){r <- 0}
  if(runif(1) < r){
    th <- alpha1_star
  }else{
    th <- alpha1
  }
  p_jump <- min(1,r)
  return(list(alpha1=th, p_jump=p_jump))
}

#alpha1_sample_norm(clust, alpha1, gamma1)

#log_lik_hurdle(alpha2, psi, gamma1, gamma2)


log_post_gamma1 <- function(clust, V, alpha1, alpha2, gamma1, gamma2, psi, sigma1){
  if((sum(abs(gamma1)>10)>0) | (sum(is.na(gamma1))>0)){
    return(-Inf)
  }else{
    like <- log_lik_zero(clust, alpha1, gamma1)
    count_like <- log_lik_hurdle(V, alpha2, psi, gamma1, gamma2)
    prior <- sum(dnorm(gamma1,0,sigma1, log=TRUE))
    
    post <- sum(like, na.rm=TRUE) + sum(prior) + sum(count_like)
    if(is.nan(post)){
      return(-Inf)
    }
    else{
      return(post)
    }
  }
  
}


gradient_gamma1 <- function(clust, V, alpha1, alpha2, gamma1, gamma2, psi, sigma1){
  d <- length(gamma1)
  diff <- rep(NA, d)
  
  q <- (1/sigma1^2)
  
  expon <- exp(X%*%alpha1 + gamma1[ll])
  expon2 <- exp(X%*%alpha2 + gamma2[ll] + psi*gamma1[ll])
  
  lik_dot <- clust - expon/(1+expon)
  
  lik_dot2 <- psi*V - psi*K*expon2/(1 + expon2)
  
  for(m in 1:d){
    
    diff[m] <- sum(lik_dot[ll==m], na.rm=TRUE) + sum(lik_dot2[ll==m], na.rm=TRUE)
  }
  
  diff <- as.numeric(diff) - q*gamma1
  
  #diff[abs(alpha1)>10] <- 0
  
  return(diff)
}


# gradient_gamma1_numerical <- function(gamma1){
#   N <- length(gamma1)
#   diff <- rep(NA, N)
#   e <- .00001
#   for(i in 1:N){
#     th_hi <- gamma1
#     th_lo <- gamma1
#     th_hi[i] <- gamma1[i] + e
#     th_lo[i] <- gamma1[i] - e
#     #((((1/deltat)*cbind(Z,Zt)%*%alpha1 + sigma1*th_hi[ll] + log(deltat)))-((1/deltat)*cbind(Z,Zt)%*%alpha1 + sigma1*th_lo[ll] + log(deltat)))/(2*e)
#     diff[i]<- (log_post_gamma1(th_hi) - log_post_gamma1(th_lo))/(2*e)
#   }
#   return (diff)
# }
# 
# diff1_g1 <- gradient_gamma1(gamma1)
# diff_g1_num <- gradient_gamma1_numerical(gamma1)
# diff1_g1 - diff_g1_num



gradient2_gamma1 <- function(clust, V, alpha1, alpha2, gamma1, gamma2, psi, sigma1){
  d <- length(gamma1)
  diff <- rep(NA, d)
  
  q <- (1/sigma1^2)
  
  expon <- exp(X%*%alpha1 + gamma1[ll])
  expon2 <- exp(X%*%alpha2 + gamma2[ll] + psi*gamma1[ll])
  
  lik_dot <- - expon/((1+expon)^2) *!is.na(clust)
  lik_dot2 <- - psi*psi*K*expon2/((1+expon2)^2)*!is.na(V)
  
  for(m in 1:d){
    
    diff[m] <- sum(lik_dot[ll==m], na.rm=TRUE) + sum(lik_dot2[ll==m], na.rm=TRUE)
  }
  
  diff <- as.numeric(diff) - q
  
  #diff[abs(alpha1)>10] <- 0
  
  return(diff)
}


# 
# gradient2_gamma1_numerical <- function(gamma1){
#   N <- length(gamma1)
#   diff <- rep(NA, N)
#   e <- .00001
#   for(i in 1:N){
#     th_hi <- gamma1
#     th_lo <- gamma1
#     th_hi[i] <- gamma1[i] + e
#     th_lo[i] <- gamma1[i] - e
#     #((((1/deltat)*cbind(Z,Zt)%*%alpha1 + sigma1*th_hi[ll] + log(deltat)))-((1/deltat)*cbind(Z,Zt)%*%alpha1 + sigma1*th_lo[ll] + log(deltat)))/(2*e)
#     diff[i]<- (gradient_gamma1(th_hi)[i] - gradient_gamma1(th_lo)[i])/(2*e)
#   }
#   return (diff)
# }

#  diff2_g1 <- gradient2_gamma1(gamma1)
#  diff2_g1_num <- gradient2_gamma1_numerical(gamma1)
#  head(diff2_g1)
#  head(diff2_g1_num)
#
# diff2_g1 - diff2_g1_num

gamma1_norm_approx <- function(clust, V, alpha1, alpha2, gamma1, gamma2, psi, sigma1){
  alp1 <- gamma1
  
  
  for(step in 1:steps.g1){
    grad1 <- gradient_gamma1(clust, V, alpha1, alpha2, alp1, gamma2, psi, sigma1)
    grad2 <- gradient2_gamma1(clust, V, alpha1, alpha2, alp1, gamma2, psi, sigma1)
    
    Q_approx <- - diag(grad2) 
    mode_approx <- c(alp1) - (1/grad2)*c(grad1)
    alp1 <- mode_approx
  }
  return(list(mode_approx=c(mode_approx), Q_approx=Q_approx))
}


 gamma1_sample_norm <- function(clust, V, alpha1, alpha2, gamma1, gamma2, psi, sigma1){
   a1_approx <-  gamma1_norm_approx(clust, V, alpha1, alpha2, gamma1, gamma2, psi, sigma1)
   
   
   gamma1_star <- as.numeric(rmvnorm(1, a1_approx$mode_approx, solve(a1_approx$Q_approx)))
   
   
   a1_approx_star <-  gamma1_norm_approx(clust, V, alpha1, alpha2, gamma1_star, gamma2, psi, sigma1)
   
   log_post_old <-  log_post_gamma1(clust, V, alpha1, alpha2, gamma1, gamma2, psi, sigma1)
   log_post_star <- log_post_gamma1(clust, V, alpha1, alpha2, gamma1_star, gamma2, psi, sigma1)
   
   J_old <- sum(dnorm(gamma1, a1_approx_star$mode_approx, sqrt(diag(1/a1_approx_star$Q_approx)), log=TRUE))
   J_star <- sum(dnorm(gamma1_star, a1_approx$mode_approx, sqrt(diag(1/a1_approx$Q_approx)), log=TRUE))
   
   r <- exp(log_post_star + sum(J_old) - log_post_old - sum(J_star))
   if(!is.finite(r)){r <- 0}
   if(runif(1) < r){
     th <- gamma1_star
   }else{
     th <- gamma1
   }
   p_jump <- min(1,r)
   return(list(gamma1=th, p_jump=p_jump))
 }

 
 # 
 # steps.g1 <- 2
 # 
 # gamma1_sample_norm(clust, V, alpha1, alpha2, gamma1, gamma2, psi, sigma1)

 
 log_post_sigma1 <- function(gamma1, sigma1){
    
   like <-  sum(dnorm(gamma1, 0, sigma1, log=TRUE))
   
   prior <- dnorm(sigma1, prior.sigma1.mn, prior.sigma1.sd, log=TRUE )
   
   like+prior
 }
 
 
 #prior.sigma1.mn <- 0
 #prior.sigma1.sd <- 1
 #log_post_sigma1(sigma1)
 
 
 
 sigma1_update <- function(gamma1, sigma1){
   sigma_star <- rnorm(1, sigma1, sigma1_jump)
   
   
   r=0
   p_jump <- 0
   if(sigma_star > 0){
     #Continue for positive components tryCatch(if(runif(1) < r){sigma <- sigma_star}, error=function(w) print(r))
     
     log_post_old <- log_post_sigma1 (gamma1, sigma1)
     log_post_star <- log_post_sigma1 (gamma1, sigma_star)
     r <- log_post_star - log_post_old
     
     if(log(runif(1)) < r){sigma1 <- sigma_star}
     p_jump <- exp(min(r,log(1)))
     
   }
   
   return (list(sigma1=sigma1, p_jump=p_jump))
 }
 
 #sigma1_update(sigma1)
 
 
 ##########Count Model
 
 log_binom <- function(x, theta){
   loglik <- rep(NA, length(theta))
   loglik <- sum(x*log(theta) + (K-x) * log((1-theta)), na.rm=TRUE)
   loglik
 }
 

 
 #log_lik_hurdle(alpha2, psi, gamma1, gamma2)
 
 log_post_alpha2 <- function(V, alpha2, psi, gamma1, gamma2){
   if((sum(abs(alpha2)>10)>0) | (sum(is.na(alpha2))>0)){
     return(-Inf)
   }else{
     like <- log_lik_hurdle(V, alpha2, psi, gamma1, gamma2)
     prior <- sum(dnorm(alpha2,prior.mn2,prior.sd2, log=TRUE))
     
     post <- sum(like, na.rm=TRUE) +sum(prior)
     if(is.nan(post)){
       return(-Inf)
     }
     else{
       return(post)
     }
   }
   
 }
 
 #log_post_alpha2(alpha2)
 
 
 
 
 
 gradient_alpha2 <- function(V, alpha2, psi, gamma1, gamma2){
   d <- length(alpha2)
   diff <- rep(NA, d)
   
   q <- (1/prior.sd2^2)*diag(d)
   
   expon <- exp(X%*%alpha2 + gamma2[ll] + psi*gamma1[ll])
   
   lik_dot <- V - K*expon/(1+expon)
   
   for(m in 1:d){
     
     diff[m] <- sum(X[,m]*lik_dot, na.rm=TRUE)
   }
   
   diff <- as.numeric(diff) - q%*%alpha2
   
   #diff[abs(alpha1)>10] <- 0
   
   return(diff)
 }
 
 # 
 # gradient_alpha2_numerical <- function(alpha2){
 #   N <- length(alpha2)
 #   diff <- rep(NA, N)
 #   e <- .00001
 #   for(i in 1:N){
 #     th_hi <- alpha2
 #     th_lo <- alpha2
 #     th_hi[i] <- alpha2[i] + e
 #     th_lo[i] <- alpha2[i] - e
 #     #((((1/deltat)*cbind(Z,Zt)%*%alpha1 + sigma1*th_hi[ll] + log(deltat)))-((1/deltat)*cbind(Z,Zt)%*%alpha1 + sigma1*th_lo[ll] + log(deltat)))/(2*e)
 #     diff[i]<- (log_post_alpha2(th_hi) - log_post_alpha2(th_lo))/(2*e)
 #   }
 #   return (diff)
 # }
 # 
 # diff1_a2 <- gradient_alpha2(alpha2)
 # diff_a2_num <- gradient_alpha2_numerical(alpha2)
 # diff1_a2
 # diff_a2_num
 #gradient_alpha1(alpha1)
 
 
 
 gradient2_alpha2 <- function(alpha2, psi, gamma1, gamma2){
   d <- length(alpha2)
   diff <- rep(NA, d)
   
   q <- (1/prior.sd2^2)
   
   expon <- exp(X%*%alpha2 + gamma2[ll] + psi*gamma1[ll])
   
   lik_dot <- - K*expon/((1+expon)^2)
   
   for(m in 1:d){
     
     diff[m] <- sum(X[,m]*lik_dot, na.rm=TRUE)
   }
   
   diff <- as.numeric(diff) - q
   
   #diff[abs(alpha1)>10] <- 0
   
   return(diff)
 }
 
 
 
 # gradient2_alpha2_numerical <- function(alpha2){
 #   N <- length(alpha2)
 #   diff <- rep(NA, N)
 #   e <- .00001
 #   for(i in 1:N){
 #     th_hi <- alpha2
 #     th_lo <- alpha2
 #     th_hi[i] <- alpha2[i] + e
 #     th_lo[i] <- alpha2[i] - e
 #     #((((1/deltat)*cbind(Z,Zt)%*%alpha1 + sigma1*th_hi[ll] + log(deltat)))-((1/deltat)*cbind(Z,Zt)%*%alpha1 + sigma1*th_lo[ll] + log(deltat)))/(2*e)
 #     diff[i]<- (gradient_alpha2_numerical(th_hi)[i] - gradient_alpha2_numerical(th_lo)[i])/(2*e)
 #   }
 #   return (diff)
 # }
 # 
 # diff2_a2 <- gradient2_alpha2(alpha2)
 # diff2_a2_num <- gradient2_alpha2_numerical(alpha2)
 # diff2_a2
 # diff2_a2_num
 
 
 
 alpha2_norm_approx <- function(V, alpha2, psi, gamma1, gamma2){
   alp1 <- alpha2
   
   
   for(step in 1:steps.a2){
     grad1 <- gradient_alpha2(V, alp1, psi, gamma1, gamma2)
     grad2 <- gradient2_alpha2(alp1, psi, gamma1, gamma2)
     
     Q_approx <- - diag(grad2) 
     mode_approx <- c(alp1) - diag(1/grad2)%*%c(grad1)
     alp1 <- mode_approx
   }
   return(list(mode_approx=c(mode_approx), Q_approx=Q_approx))
 }
 
 #Good starting value for alpha1:
 #alpha1 <- c(  0.50338437, -0.81447167, -1.01128499, -1.48048324, -0.08876805, 0.07802240, -0.01134941)
 # steps.a1 <- 1
 # alpha2_norm_approx(alpha2)
 
 
 
 
 
 alpha2_sample_norm <- function(V, alpha2, psi, gamma1, gamma2){
   a2_approx <- alpha2_norm_approx(V, alpha2, psi, gamma1, gamma2)
   
   alpha2_star <- as.numeric(rmvnorm(1, a2_approx$mode_approx, solve(a2_approx$Q_approx)))
   
   
   a2_approx_star <- alpha2_norm_approx(V, alpha2_star, psi, gamma1, gamma2)
   
   log_post_old <-  log_post_alpha2(V, alpha2, psi, gamma1, gamma2)
   log_post_star <- log_post_alpha2(V, alpha2_star, psi, gamma1, gamma2)
   
   J_old <- sum(dnorm(alpha2, a2_approx_star$mode_approx, sqrt(diag(1/a2_approx_star$Q_approx)), log=TRUE))
   J_star <- sum(dnorm(alpha2_star, a2_approx$mode_approx, sqrt(diag(1/a2_approx$Q_approx)), log=TRUE))
   
   r <- exp(log_post_star + sum(J_old) - log_post_old - sum(J_star))
   if(!is.finite(r)){r <- 0}
   if(runif(1) < r){
     th <- alpha2_star
   }else{
     th <- alpha2
   }
   p_jump <- min(1,r)
   return(list(alpha2=th, p_jump=p_jump))
 }
 
 
 #alpha2_sample_norm(V, alpha2, psi, gamma1, gamma2)
 
 
 log_post_gamma2 <- function(V, alpha2, psi, gamma1, gamma2, sigma2){
   n <- length(gamma2)
   linmod <- X%*%alpha2 + gamma2[ll] + psi*gamma1[ll]
   like <- rep(NA, n)
   prior <- rep(NA, n)
   for(i in 1:n){
     like[i] <- sum(V[ll==i]*linmod[ll==i] - K*log(1+exp(linmod[ll==i])), na.rm=TRUE)
     prior[i] <- dnorm(gamma2[i],0,sigma2, log=TRUE)
   }
   post <- like + prior
   post[is.nan(post)] <- -Inf
  return(post)
   
 }
 
 #log_post_gamma2(V, alpha2, psi, gamma1, gamma2, sigma2)
 
 gradient_gamma2 <- function(V, alpha2, psi, gamma1, gamma2, sigma2){
   d <- length(gamma2)
   diff <- rep(NA, d)
   
   q <- (1/sigma2^2)
   
   expon <- exp(X%*%alpha2 + gamma2[ll] + psi*gamma1[ll])
   
   lik_dot <- V - K*expon/(1+expon)
   
   for(m in 1:d){
     
     diff[m] <- sum(lik_dot[ll==m], na.rm=TRUE)
   }
   
   diff <- as.numeric(diff) - q*gamma2
   
   #diff[abs(alpha1)>10] <- 0
   
   return(diff)
 }
 
 
 # gradient_gamma2_numerical <- function(gamma2){
 #   N <- length(gamma2)
 #   diff <- rep(NA, N)
 #   e <- .00001
 #   for(i in 1:N){
 #     th_hi <- gamma2
 #     th_lo <- gamma2
 #     th_hi[i] <- gamma2[i] + e
 #     th_lo[i] <- gamma2[i] - e
 #     #((((1/deltat)*cbind(Z,Zt)%*%alpha1 + sigma1*th_hi[ll] + log(deltat)))-((1/deltat)*cbind(Z,Zt)%*%alpha1 + sigma1*th_lo[ll] + log(deltat)))/(2*e)
 #     diff[i]<- (log_post_gamma2(th_hi) - log_post_gamma2(th_lo))/(2*e)
 #   }
 #   return (diff)
 # }
 # 
 #  diff1_g2 <- gradient_gamma2(gamma2)
 #  diff_g2_num <- gradient_gamma2_numerical(gamma2)
 #  diff1_g2 - diff_g2_num

 
 
 gradient2_gamma2 <- function(V, alpha2, psi, gamma1, gamma2, sigma2){
   d <- length(gamma2)
   diff <- rep(NA, d)
   
   q <- (1/sigma2^2)
   
   expon <- exp(X%*%alpha2 + gamma2[ll] + psi*gamma1[ll])
   
   lik_dot <- (- K*expon/((1+expon)^2))*!is.na(V)
   
   for(m in 1:d){
     
     diff[m] <- sum(lik_dot[ll==m], na.rm=TRUE)
   }
   
   diff <- as.numeric(diff) - q
   
   #diff[abs(alpha1)>10] <- 0
   
   return(diff)
 }
 
 

 # gradient2_gamma2_numerical <- function(alpha2, psi, gamma1, gamma2, sigma2){
 #   N <- length(gamma2)
 #   diff <- rep(NA, N)
 #   e <- .00001
 #   for(i in 1:N){
 #     th_hi <- gamma2
 #     th_lo <- gamma2
 #     th_hi[i] <- gamma2[i] + e
 #     th_lo[i] <- gamma2[i] - e
 #     #((((1/deltat)*cbind(Z,Zt)%*%alpha1 + sigma1*th_hi[ll] + log(deltat)))-((1/deltat)*cbind(Z,Zt)%*%alpha1 + sigma1*th_lo[ll] + log(deltat)))/(2*e)
 #     diff[i]<- (gradient_gamma2(V, alpha2, psi, gamma1, th_hi, sigma2)[i] - gradient_gamma2(V, alpha2, psi, gamma1, th_lo, sigma2)[i])/(2*e)
 #   }
 #   return (diff)
 # }
 # 
 # diff2_g2 <- gradient2_gamma2(V, alpha2, psi, gamma1, gamma2, sigma2)
 # diff2_g2_num <- gradient2_gamma2_numerical(alpha2, psi, gamma1, gamma2, sigma2)
 # head(diff2_g2)
 # head(diff2_g2_num)
 # 
 # diff2_g2 - diff2_g2_num
 
 gamma2_norm_approx <- function(V, alpha2, psi, gamma1, gamma2, sigma2){
   alp1 <- gamma2
   
   
   for(step in 1:steps.g2){
     grad1 <- gradient_gamma2(V, alpha2, psi, gamma1, alp1, sigma2)
     grad2 <- gradient2_gamma2(V, alpha2, psi, gamma1, alp1, sigma2)
     
     Q_approx <- - grad2
     mode_approx <- c(alp1) - (1/grad2)*c(grad1)
     alp1 <- mode_approx
   }
   return(list(mode_approx=c(mode_approx), Q_approx=Q_approx))
 }
 
 
 # gamma2_sample_norm <- function(V, alpha2, psi, gamma1, gamma2, sigma2){
 #   a1_approx <-  gamma2_norm_approx(V, alpha2, psi, gamma1, gamma2, sigma2)
 #   
 #   
 #   gamma2_star <- as.numeric(rmvnorm(1, a1_approx$mode_approx, solve(a1_approx$Q_approx)))
 #   
 #   
 #   a1_approx_star <-  gamma2_norm_approx(V, alpha2, psi, gamma1, gamma2_star, sigma2)
 #   
 #   log_post_old <-  log_post_gamma2(V, alpha2, psi, gamma1, gamma2, sigma2)
 #   log_post_star <- log_post_gamma2(V, alpha2, psi, gamma1, gamma2_star, sigma2)
 #   
 #   J_old <- sum(dnorm(gamma2, a1_approx_star$mode_approx, sqrt(diag(1/a1_approx_star$Q_approx)), log=TRUE))
 #   J_star <- sum(dnorm(gamma2_star, a1_approx$mode_approx, sqrt(diag(1/a1_approx$Q_approx)), log=TRUE))
 #   
 #   r <- exp(log_post_star + sum(J_old) - log_post_old - sum(J_star))
 #   if(!is.finite(r)){r <- 0}
 #   if(runif(1) < r){
 #     th <- gamma2_star
 #   }else{
 #     th <- gamma2
 #   }
 #   p_jump <- min(1,r)
 #   return(list(gamma2=th, p_jump=p_jump))
 # }
 
 
 
 gamma2_sample_norm <- function(V, alpha2, psi, gamma1, gamma2, sigma2){
   n <- length(gamma2)
   
   a1_approx <-  gamma2_norm_approx(V, alpha2, psi, gamma1, gamma2, sigma2)
   
   
   gamma2_star <- as.numeric(rnorm(n, a1_approx$mode_approx, 1/sqrt(a1_approx$Q_approx)))
   
   
   a1_approx_star <-  gamma2_norm_approx(V, alpha2, psi, gamma1, gamma2_star, sigma2)
   
   J_old <- rep(NA, n)
   J_star <- rep(NA, n)
   #not_pos_semi <- rep(0,n)
   
   log_post_old <-  log_post_gamma2(V, alpha2, psi, gamma1, gamma2, sigma2)
   log_post_star <- log_post_gamma2(V, alpha2, psi, gamma1, gamma2_star, sigma2)
   
   J_old <- dnorm(gamma2, a1_approx_star$mode_approx, 1/sqrt(a1_approx_star$Q_approx), log=TRUE)
   J_star <- dnorm(gamma2_star, a1_approx$mode_approx, 1/sqrt(a1_approx$Q_approx), log=TRUE)
   
   r <- exp(log_post_star + J_old - log_post_old - J_star)
   r[r > 1] <- 1
   r[is.nan(r)] <- 0
   
   jump <- which(runif(n) < r)
   th <- gamma2
   th[jump] <- gamma2_star[jump]
   p_jump <- r
   return(list(gamma2=th, p_jump=mean(p_jump, na.rm=TRUE)))
 }
 
 
 #gamma2_sample_norm(V, alpha2, psi, gamma1, gamma2, sigma2)
 
 
 
 # inits <- list(alpha1=alpha1, alpha2=alpha2, gamma1=gamma1, gamma2=gamma2, psi=psi, sigma1=sigma1, sigma2=sigma2)
 # saveRDS(inits, file="InitsRI.rds")
 # 

 
 
 
 
 
 
 log_post_sigma2<-function(gamma2, sigma2){
   
   prior_gamma <-  sum(dnorm(gamma2, 0, sigma2, log=TRUE))
   
   prior <- dnorm(sigma2, prior.sigma2.mn, prior.sigma2.sd, log=TRUE )
   
   prior_gamma+prior
 }
 
 
 
 
 
 
 
 sigma2_update <- function(gamma2, sigma2){
   sigma_star <- rnorm(1, sigma2, sigma2_jump)
   
   
   r=0
   p_jump <- 0
   if(sigma_star > 0){
     #Continue for positive components tryCatch(if(runif(1) < r){sigma <- sigma_star}, error=function(w) print(r))
     
     log_post_old <- log_post_sigma2 (gamma2, sigma2)
     log_post_star <- log_post_sigma2 (gamma2, sigma_star)
     r <- log_post_star - log_post_old
     
     if(log(runif(1)) < r){sigma2 <- sigma_star}
     p_jump <- exp(min(r,log(1)))
     
   }
   
   return (list(sigma2=sigma2, p_jump=p_jump))
 }
 
 
 
 
 
 
 #sigma2_update(gamma2, sigma2)
 
 log_post_psi <- function(V, alpha2, psi, gamma1, gamma2){
   like <- log_lik_hurdle(V, alpha2, psi, gamma1, gamma2)
   prior <- dnorm(psi, prior.psi.mn, prior.psi.sd, log=TRUE)
   
   like+prior
 }
 
 
 
 psi_update <- function(V, alpha2, psi, gamma1, gamma2){
   
   psi_star <- rnorm(1, psi, psi_jump)
   
   
   r=0
   p_jump <- 0
   if(psi_star>0){
     #Continue for positive components
     log_post_old <- log_post_psi (V, alpha2, psi, gamma1, gamma2)
     log_post_star <- log_post_psi (V, alpha2, psi_star, gamma1, gamma2)
     r <- log_post_star - log_post_old
     
     if(log(runif(1)) < r){psi <- psi_star}
     p_jump <- exp(min(r,log(1)))
     
   }
   
   return (list(psi=psi, p_jump=p_jump))
   
 }
 
 
 
 #psi_update(alpha2, psi, gamma1, gamma2)
 
 



count_step <- function(V, alpha2, gamma1, gamma2, sigma2, psi){

  alpha2_temp <- alpha2_sample_norm(V, alpha2, psi, gamma1, gamma2)
  alpha2 <- alpha2_temp$alpha2
  pj_alpha2 <- alpha2_temp$p_jump
  
  gamma2_temp <- gamma2_sample_norm(V, alpha2, psi, gamma1, gamma2, sigma2)
  gamma2 <- gamma2_temp$gamma2
  pj_gamma2 <- gamma2_temp$p_jump
  
  sigma2_temp <- sigma2_update(gamma2, sigma2)
  sigma2<-sigma2_temp$sigma2
  pj_sigma2 <- sigma2_temp$p_jump

  psi_temp <- psi_update(V, alpha2, psi, gamma1, gamma2)
  psi <- psi_temp$psi
  pj_psi <- psi_temp$p_jump
  
  
  return(list(alpha2 = alpha2, psi=psi, gamma2 = gamma2, sigma2 = sigma2, pj_alpha2 = pj_alpha2, pj_psi = pj_psi, pj_gamma2 = pj_gamma2, pj_sigma2 = pj_sigma2))
}


######## Zero model update:


zero_step <- function(clust, V, alpha1, alpha2, gamma1, gamma2, sigma1, sigma2, psi){
  
  alpha1_temp <- alpha1_sample_norm(clust, alpha1, gamma1)
  alpha1 <- alpha1_temp$alpha1
  pj_alpha1 <- alpha1_temp$p_jump
  
  gamma1_temp <- gamma1_sample_norm(clust, V, alpha1, alpha2, gamma1, gamma2, psi, sigma1)
  gamma1 <- gamma1_temp$gamma1
  pj_gamma1 <- gamma1_temp$p_jump
  
  sigma1_temp <- sigma1_update(gamma1, sigma1)
  sigma1 <- sigma1_temp$sigma1
  pj_sigma1 <- sigma1_temp$p_jump
  
  return(list(alpha1 = alpha1, gamma1 = gamma1, sigma1 = sigma1, pj_alpha1 = pj_alpha1, pj_gamma1 = pj_gamma1, pj_sigma1 = pj_sigma1))
}

#zero_step(clust, V, alpha1, alpha2, gamma1, gamma2, sigma1, sigma2, psi)

#LLPD Function
pd <- function(alpha1, alpha2, gamma1, gamma2, psi){
  lik <- rep(NA,length(y))
  
  theta <- expit(X%*%alpha2 + gamma2[ll] + psi*gamma1[ll])
  p <- expit(X%*%alpha1 + gamma1[ll])
  
  pos_clus <- which((y>0)&(!is.na(y)))
  lik[pos_clus] = p[pos_clus] * dbinom(y[pos_clus], 90, theta[pos_clus])
  
  zero_clus <- which((y==0)&(!is.na(y)))
  lik[zero_clus] <-  1 - p[zero_clus] + p[zero_clus]*(theta[zero_clus]^90)
  
  lik
  
}

#pd(alpha1, alpha2, gamma1, gamma2)




chain_run <- function(inits,  iter=2000, this.chain, thin=1)
{
  
  ##Initialize everything!
  
  iter.keep <- floor(iter/thin)
  
  #Initialize values for each chain (would be better to give them different starting points)
  alpha1 <- inits[[this.chain]]$alpha1
  sigma1 <- inits[[this.chain]]$sigma1
  gamma1<- inits[[this.chain]]$gamma1
  alpha2 <- inits[[this.chain]]$alpha2
  psi <- inits[[this.chain]]$psi
  sigma2 <- inits[[this.chain]]$sigma2
  gamma2<- inits[[this.chain]]$gamma2
  
  p_jump1 <- array (NA, c(iter.keep, 3))
  sims1 <- array (NA, c(iter.keep, 8))
  
  #Adjust this for gamma2 matrix
  p_jump2 <- array (NA, c(iter.keep, 4))
  sims2 <- array (NA, c(iter.keep, P+2))
  gamma2_sims<- array(NA, c(iter.keep, length(gamma2)))
  gamma1_sims<- array(NA, c(iter.keep, length(gamma2)))
  lik <- array(NA, c(iter.keep, length(y)))

  for(i in 1:iter){
    
    #Data Augmentation step
    dat <- samp.clust(alpha1, alpha2, gamma1, gamma2, psi)
    V <- dat$V
    clust <- dat$clust
    
    
    update1 <- zero_step(clust, V, alpha1, alpha2, gamma1, gamma2, sigma1, sigma2, psi)
    
    alpha1 <- update1$alpha1
    gamma1 <- update1$gamma1
    sigma1 <- update1$sigma1
    
    
    update2 <- count_step(V, alpha2, gamma1, gamma2, sigma2, psi)
    
    alpha2 <- update2$alpha2
    gamma2 <- update2$gamma2
    sigma2 <- update2$sigma2
    psi <- update2$psi
    if(i%%100==0){
      cat(i, "\t")
    }
    
    if(i%%thin==0)
    {
      ii=floor(i/thin)  
      sims1[ii,] <- c (alpha1, sigma1)
      gamma1_sims[ii,]<-gamma1
      p_jump1[ii,]<-c(update1$pj_alpha1, update1$pj_gamma1, update1$pj_sigma1)
      
      sims2[ii,] <- c (alpha2,  psi, sigma2)
      gamma2_sims[ii,]<-gamma2
      p_jump2[ii,]<-c(update2$pj_alpha2, update2$pj_psi, update2$pj_gamma2, update2$pj_sigma2)
      lik[ii,] <- pd(alpha1, alpha2, gamma1, gamma2, psi)
    }
  }
  
  inits_update <- list(alpha1=alpha1, sigma1=sigma1, gamma1=gamma1, alpha2=alpha2, psi=psi, sigma2=sigma2, gamma2=gamma2)
  
  return(list(sims1=sims1, p_jump1=p_jump1, sims2=sims2, gamma1_sims=gamma1_sims, gamma2_sims=gamma2_sims, p_jump2=p_jump2, lik=lik, inits_update=inits_update))
  
}

#chain_run(inits,  iter=10, this.chain, thin=1)

#RUN IT

zib_ri_run <- function(inits, iter=2000, chains=4, thin=1){   
  #number of chains  
  
  iter.keep <- floor(iter/thin)
  
  
  p_jump1 <- array (NA, c(iter.keep, chains, 3))
  sims1 <- array (NA, c(iter.keep, chains, 8))
  dimnames (sims1) <- list (NULL, NULL,
                            c (paste ("beta1[", 1:7, "]", sep=""),
                               "sigma1"
                            ))
  
  #Adjust this for gamma2 matrix
  p_jump2 <- array (NA, c(iter.keep, chains, 4))
  sims2 <- array (NA, c(iter.keep, chains, 9))
  dimnames (sims2) <- list (NULL, NULL,
                            c (paste ("beta2[", 1:7, "]", sep=""),
                               "psi",
                               "sigma2"
                            ))
  
  
  gamma2_sims<- array(NA, c(iter.keep, chains, 718))
  
  dimnames(gamma2_sims) <- list(NULL, NULL, c(paste("gamma2[", 1:718, "]", sep="")))
  
  
  gamma1_sims<- array(NA, c(iter.keep, chains, 718))
  dimnames(gamma1_sims) <- list(NULL, NULL, c(paste("gamma1[", 1:718, "]", sep="")))
  
  
  
  
  
  lik <- array(NA, c(iter.keep, chains, length(y)))
  
  inits_update <- list(NA)
  
  chain.list <- 1:chains
  result <- mclapply(chain.list, function(x){chain_run(inits=inits,  iter=iter, this.chain=x, thin=thin)}, mc.cores=4)
  
  for(jj in 1:chains){
    sims1[,jj,] <- result[[jj]]$sims1
    sims2[,jj,] <- result[[jj]]$sims2
    p_jump1[,jj,] <- result[[jj]]$p_jump1
    p_jump2[,jj,] <- result[[jj]]$p_jump2
    gamma1_sims[,jj,] <- result[[jj]]$gamma1_sims
    gamma2_sims[,jj,] <- result[[jj]]$gamma2_sims
    
    lik[,jj,] <- result[[jj]]$lik
    inits_update[[jj]] <- result[[jj]]$inits_update
  }
  
  return(list(sims1=sims1, p_jump1=p_jump1, sims2=sims2, gamma1_sims=gamma1_sims, gamma2_sims=gamma2_sims, p_jump2=p_jump2, lik=lik, inits_update=inits_update) )
}










##################################################################################################################
############# LOADING IN DATA AND RUNNING SAMPLER










#Load in Data
# setwd("~/Projects/TVRE ZIB")

sbirt<-read.csv(file="sbirt2baseline.csv", header=T)

sbirt$Group=sbirt$Group-1
sbirt$Site=sbirt$Site-1

heavy_dat <-sbirt[!is.na(sbirt$Heavy_Drinkingdays),]

K=90

heavy_dat <- subset(heavy_dat, select = c("ID1", "VISITNUM", "Group", "Heavy_Drinkingdays"))


wide_dat <- reshape(heavy_dat, idvar ="ID1", timevar = "VISITNUM", direction = "wide")



group <- cbind(1:718,wide_dat$Group.0)
colnames(group) = c("id", "group")

y <- wide_dat[,c(3, 5, 7, 9)]

y <- c(t(y))

ll <- rep(1:718, each =4)

groups <- group[ll,]

zmat <- rbind(c(1,0,0,0), c(0,1,0,0),c(0,0,1,0), c(0,0,0,1))
Z<- matrix(rep(t(zmat), 718), ncol=ncol(zmat), byrow=TRUE)
X <- cbind(Z, Z[,-1]*groups[,2])

P <- ncol(X)
N <- length(unique(ll))
J <- 4 

inits <- readRDS("InitsRI.rds")


########## PRIORS AND TUNING PARAMETERS
  prior.mn1 = 0
  prior.sd1 = 1
  prior.sigma1.mn = 0
  prior.sigma1.sd = .5
  
  prior.mn2 = 0
  prior.sd2 = 1
  prior.psi.mn = 0
  prior.psi.sd = 1
  prior.sigma2.mn = 0
  prior.sigma2.sd = .5

  
  steps.a1 <- 2
  steps.a2 <- 2
  steps.g1 <- 3
  steps.g2 <- 3
  sigma1_jump <- .1
  sigma2_jump <- .2
  psi_jump <- .05





start <- Sys.time()
binom_hurdle_ri  <- zib_ri_run(inits, iter=1000, chains=4, thin=1)
Sys.time()-start

#alpha1, gamma1, sigma1
apply(binom_hurdle_ri$p_jump1,c(2,3), mean)

#alpha2, psi2, gamma2, sigma2
apply(binom_hurdle_ri$p_jump2,c(2,3), mean)



binom_hurdle_ri$sims1[,1,]

saveRDS(binom_hurdle_ri, file="ResultsRI.rds")
saveRDS(binom_hurdle_ri$inits_update, file="InitsRI.rds")






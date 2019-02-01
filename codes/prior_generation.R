### Generate data from the prior
nu = 3
a1 = 3; a2 = 4
as = 1; bs = 0.3

S = 5000
k = 100
p = 10

beta_X = matrix(0,p,S)
Omega_X = array(0,c(p,p,S))

beta2_X = matrix(0,p,S)
Omega2_X = array(0,c(p,p,S))


for (i in 1:S){
   
   delta = c(rgamma(1,a1,1),rgamma(k-1,a2,1))
   tau = cumprod(delta)
   phi = rgamma(p*k,nu/2,nu/2)
   phi = matrix(phi,p,k)
   ps = rgamma(p,as,bs)
   
   Lambda = matrix(0,p,k)
   Lambda2 = matrix(0,p,k)
   for(j in 1:p){
      for(h in 1:k){
         Lambda[j,h] = rnorm(1,0,1/sqrt(tau[h]*phi[j,h]))    
         Lambda2[j,h] = rnorm(1,0,1/sqrt(phi[j,h]))  
      }
   }
   
   omega = rnorm(k,0,1)
   Omega = matrix(0,k,k)
   Omega[lower.tri(Omega)] = rnorm(k*(k-1)/2,0,1)
   diag(Omega) = rnorm(k,0,1)
   Omega = (Omega + t(Omega))/2

   A = solve(t(Lambda)%*%diag(ps)%*%Lambda + diag(k))%*%t(Lambda)%*%diag(ps)
   A2 = solve(t(Lambda2)%*%diag(ps)%*%Lambda + diag(k))%*%t(Lambda2)%*%diag(ps)
   
   beta_X[,i] = t(A)%*%omega
   Omega_X[,,i] = t(A)%*%Omega%*%A  
   
   beta2_X[,i] = t(A2)%*%omega
   Omega2_X[,,i] = t(A2)%*%Omega%*%A2  
   
   if (i%%1000==0){
      print(i)
   }
}


hist(Omega_X[1,1,],freq = F,breaks = seq(-15,15,by =0.5))
lines(density(Omega_X[1,1,]),col = "red")
hist(Omega2_X[1,1,],freq = F)
lines(density(Omega2_X[1,1,]),col = "red")

plot(density(Omega_X[1,1,]),col = "red")
plot(density(Omega_X[1,2,]),col = "red")

par(mfrow = c(2,1))
plot(density(beta_X[1,]),col = "red")
plot(density(beta2_X[1,]),col = "red")

# the coefficients of quadratic effects has heavier tails than the one of interactions
# both have a spike in zero

# super heavy tails if we do not consider the shrinkage parameters tau
# probabibly the distribution for beta_X and Omega_X is improper
# might be a good prior distribution

eps_seq = seq(1,0.01,by = -0.01)
probs = numeric(length(eps_seq))

for (i in 1:length(eps_seq)){
   epsilon = eps_seq[i]
   probs[i] = mean(Omega_X > -epsilon & Omega_X < epsilon)
   print(i/length(eps_seq))
}
   
plot(eps_seq,probs)

sigma1 = rgamma(50000,2,1)
sigma2 = 1/sigma1
par(mfrow = c(1,2))
plot(density(sigma1))
plot(density(sigma2))


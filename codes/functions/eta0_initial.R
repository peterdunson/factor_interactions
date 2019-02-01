eta0_initial = function(y, X, S = 500, burn = 1){
   
   X = cbind(X,y)
   p = ncol(X)
   n = nrow(X)
   
# gibbs sampler
   nrun = S
   burn = burn
   kinit = floor(log(p)*3)
   b0 = 1
   b1 = 0.0005
   epsilon = 1e-3
   prop = 1.00
   
#scaling
   M = apply(X,2,mean)
   VX = apply(X, 2, var)
   X = as.matrix(scale(cbind(X)))
   
   num = 0
   k = kinit

# prior parameters
   as = 1;bs = 0.3;
   df = 3
   ad1 = 2.1;bd1 = 1
   ad2 = 3.1;bd2 = 1
   adf = 1; bdf = 1

# initial values
   ps = rgamma(p,as,bs)
   Sigma = diag(1/ps)
   Lambda = matrix(0,p,k)
   eta = matrix(rnorm(n*k),n,k)
   meta = matrix(0,n,k)
   veta = diag(rep(1,k))

   psijh = matrix(rgamma(p*k,df/2,df/2),p,k)
   delta = c(rgamma(1,ad1,1/bd1),rgamma(k-1,ad2,1/bd2))
   tauh = cumprod(delta)
   Plam = psijh*matrix(rep(tauh,p),p,k,byrow=T)

for(i in 1:nrun){
   
   #update eta
   Lmsg = Lambda*matrix(rep(ps,k),p,k)
   Veta1 = diag(rep(1,k)) + t(Lmsg)%*%Lambda
   TT = chol(Veta1)
   QR = qr(TT)
   R = qr.R(QR)
   S = solve(R)
   Veta = S%*%t(S)
   Meta = X%*%Lmsg%*%Veta
   eta = Meta + matrix(rnorm(n*k,0,1),n,k)%*%t(S)
   
   #update Lambda
   eta2 = t(eta)%*%eta
   for(j in 1:p){
      Qlam = diag(Plam[j,])+ps[j]*eta2
      blam = ps[j]*(t(eta)%*%X[,j])
      Llam = t(chol(Qlam))
      zlam = rnorm(k)
      vlam = solve(Llam,blam)
      mlam = solve(t(Llam),vlam)
      ylam = solve(t(Llam),zlam)
      Lambda[j,] = (ylam + mlam)
   }
   
   #update psijh's
   psijh = rgamma(p*k,df/2 + 0.5,1)
   psijh = psijh*(1/(df/2+(Lambda^2)*matrix(rep(tauh,p),p,k,byrow = T)))
   
   #update delta e tauh
   mat = psijh*(Lambda^2)
   ad = ad1 + 0.5*p*k
   bd = bd1 + 0.5*(1/delta[1])*(sum(tauh*apply(mat,2,sum)))
   delta[1] = rgamma(1,ad,bd)
   tauh = cumprod(delta)
   
   for(h in 2:(k)){
      ad = ad2 + 0.5*p*(k-h+1)
      bd = bd2 + 0.5*(1/delta[h])*
         sum(tauh[h:length(tauh)]*apply(as.matrix(mat[,h:ncol(mat)]),2,sum))
      delta[h] = rgamma(1,ad,bd)
      tauh = cumprod(delta)
   }
   
   #update Sigma
   Ytil = X - eta%*%t(Lambda)
   ps = rgamma(p,as+0.5*n,1)
   ps = (1/(bs + 0.5*apply(Ytil^2,2,sum)))*ps
   Sigma = diag(1/ps)
   
   #update precision parameters
   Plam = psijh*matrix(rep(tauh,p),p,k,byrow=T)
   
   #make adaptations
   prob = 1/(exp(b0+b1*i))
   uu = runif(1)
   lind = apply((abs(Lambda)< epsilon)/p,2,sum)
   vec = (lind >= prop)
   num = sum(vec)
   
   if (uu<prob){
      if ((i>20)&(num==0)&all(lind<0.995)){
         k = k+1
         Lambda = cbind(Lambda,0)
         eta = cbind(eta,rnorm(n))
         psijh = cbind(psijh,rgamma(p,df/2,df/2))
         delta = c(delta,rgamma(1,ad2,bd2))
         tauh = cumprod(delta)
         Plam = psijh*matrix(rep(tauh,p),p,k,byrow=T)
      }else if (num > 0){
         nonred = setdiff(1:k,which(vec !=0))
         k = max(k - num,1)
         Lambda = Lambda[,nonred]
         psijh = psijh[,nonred]
         eta = eta[,nonred]
         delta = delta[nonred]
         tauh = cumprod(delta)
         Plam = psijh*matrix(rep(tauh,p),p,k,byrow=T)
      }
   }
}
output = list(eta = eta, k0 = k)
return(output)
}
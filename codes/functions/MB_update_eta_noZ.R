MB_update_eta_noZ = function(phi,sigmasq_y,Lambda,ps,k,X,y,eta){
   #update of eta_i as the quadratic term in the y part does not exist
   #shrinks the interactions to zero
   n = nrow(X)
   p = ncol(X)
   
   Vn = solve(phi%*%t(phi)/sigmasq_y+t(Lambda)%*%diag(ps)%*%Lambda+diag(rep(1,k)))
   for (i in 1:n){
      mun = Vn%*%(t(Lambda)%*%diag(ps)%*%X[i,] + phi*(y[i])/sigmasq_y)
      eta_star = as.numeric(rmvnorm(1,mun,Vn))
      
      # Modular bayes update
      eta[i,] = eta_star
   }
   return(eta)
}

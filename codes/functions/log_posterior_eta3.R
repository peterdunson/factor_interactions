log_posterior_eta3 = function(X_i,y_i,eta,eta3,phi,phi3,Lambda,ps,k,sigmasq_y,Psi){
   #log posterior of eta including the quadratic terms
   
   -0.5*( t(eta)%*%(phi%*%t(phi)/sigmasq_y + t(Lambda)%*%diag(ps)%*%Lambda + diag(k) - 
                       2*Psi*(y_i-as.vector(t(eta3)%*%phi3))/sigmasq_y )%*%eta - 
             2*eta%*%(t(Lambda)%*%diag(ps)%*%X_i + phi*(y_i-as.vector(t(eta3)%*%phi3))/sigmasq_y) + 
             2*t(eta)%*%phi*(eta%*%Psi%*%eta)/sigmasq_y + 
             (1/sigmasq_y)*(eta%*%Psi%*%eta)^2 + t(eta3)%*%(phi3%*%t(phi3)/sigmasq_y)%*%eta3 -
             2*t(eta3)%*%phi3*y_i/sigmasq_y)
}

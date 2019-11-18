return_G = function(results,norm = T){
   
   MSE_beta = apply(results$MSE_beta,2,mean); (MSE_beta)
   FR = apply(results$FR,2,mean); (FR)
   err_pred = apply(results$err_pred,2,mean); (err_pred)
   err_test = apply(results$err_test,2,mean); (err_test)
   TP_main = apply(results$TP_main,2,mean); (TP_main)
   TN_main = apply(results$TN_main,2,mean); (TN_main)
   TP_int = apply(results$TP_int,2,mean); (TP_int)
   TN_int = apply(results$TN_int,2,mean); (TN_int)
   

   
   G = rbind(err_pred,FR,MSE_beta,
             TP_main,TN_main,TP_int,TN_int)
   
   if (norm == T){
      
      min_pred = min(err_pred[1:5])
      G[1,1:5] = G[1,1:5]/min_pred
      
      min_FR = min(FR[1:5])
      G[2,1:5] = G[2,1:5]/min_FR
      
      min_beta = min(MSE_beta[1:5])
      G[3,1:5] = G[3,1:5]/min_beta
      
   }
   
   return(G)
   
}
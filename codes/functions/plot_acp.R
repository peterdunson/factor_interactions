plot_acp = function(nrun,burn,gibbs){
   acp = gibbs$acp
   hist(acp/(nrun-burn),freq = F)
   lines(density(acp/(nrun-burn)),col="red")
}

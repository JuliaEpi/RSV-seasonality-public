# Converts the output of Turing (read as a DF into R) to a coda object for further plotting
func_tocoda <- function(chainsdf) {
  
  codaobj <- vector("list", length=length(unique(chainsdf$chain)))
  for (i in unique(chainsdf$chain)) {
    codaobj[[which(unique(chainsdf$chain)==i)]] <- as.mcmc(as.matrix(chainsdf[chainsdf$chain==i,c(3:(ncol(chainsdf)-13))]))
  }
  codaobj <- as.mcmc.list(codaobj)
  
  return(codaobj)
  
}

# Function for plotting densities (inspired by coda package)
calcdensity <- function (x, parno) {
  
  xx <- as.matrix(x)
  y <- xx[, parno, drop = TRUE]
  
  # calculate bandwidth
  y2 <- y[!is.na(as.vector(y))]
  bw <- 1.06 * min(sd(y2), IQR(y2)/1.34) * length(y2)^-0.2
  width <- 4 * bw
  scale <- "open"
  
  # is it a proportion?
  if (max(y) <= 1 && 1 - max(y) < 2 * bw) {
    if (min(y) >= 0 && min(y) < 2 * bw) {
      scale <- "proportion"
      y <- c(y, -y, 2 - y)
    }
  }
  # is it positive?
  else if (min(y) >= 0 && min(y) < 2 * bw) {
    scale <- "positive"
    y <- c(y, -y)
  }
  # not proportion and not exclusively positive
  else scale <- "open"
  
  # make density object
  dens <- density(y, width = width)
  
  if (scale == "proportion") {
    dens$y <- 3 * dens$y[dens$x >= 0 & dens$x <= 1]
    dens$x <- dens$x[dens$x >= 0 & dens$x <= 1]
  }
  else if (scale == "positive") {
    dens$y <- 2 * dens$y[dens$x >= 0]
    dens$x <- dens$x[dens$x >= 0]
  }
  dens <- data.frame(x=dens$x, y=dens$y)
  dens$parno <- parno
  return(dens)
  
}


func_compdensity <- function(chains, chains_cos, chains_2019, chains_noar) {
  
  densities <- lapply(1:nvar(chains), function(x) calcdensity(chains, x))
  densities <- plyr::ldply(densities, rbind)
  densities$model <- "Von Mises 2010-2022"
  densities$par <- factor(densities$parno, labels=c("rho0","psi","q","beta0", "eta", "omega","phi","k", "delta1", "delta2", "delta3"))
  
  densities_2019 <- lapply(1:nvar(chains_2019), function(x) calcdensity(chains_2019, x))
  densities_2019 <- plyr::ldply(densities_2019, rbind)
  densities_2019$model <- "Von Mises 2010-2019"
  densities_2019$par <- factor(densities_2019$parno, labels=c("rho0","psi","q","beta0", "eta", "omega","phi","k", "delta1", "delta2", "delta3"))
  
  densities_cos <- lapply(1:nvar(chains_cos), function(x) calcdensity(chains_cos, x))
  densities_cos <- plyr::ldply(densities_cos, rbind)
  densities_cos$model <- "Cosine 2010-2022"
  densities_cos$par <- factor(densities_cos$parno, labels=c("rho0","psi","q","beta0", "eta", "omega","phi", "delta1", "delta2", "delta3"))

  densities_noar <- lapply(1:nvar(chains_noar), function(x) calcdensity(chains_noar, x))
  densities_noar <- plyr::ldply(densities_noar, rbind)
  densities_noar$model <- "Von Mises 2010-2022, no AR"
  densities_noar$par <- factor(densities_noar$parno, labels=c("rho0","psi","q","beta0", "eta", "omega","phi","k", "delta1", "delta2", "delta3"))
  
  
  comp_dens <-  merge(densities, 
                      densities_2019,
                      by=intersect(names(densities), names(densities_2019)), 
                      all=T)
  
  comp_dens <-  merge(comp_dens, 
                      densities_cos,
                      by=intersect(names(comp_dens), names(densities_cos)), 
                      all=T)
  
  comp_dens <-  merge(comp_dens, 
                      densities_noar,
                      by=intersect(names(comp_dens), names(densities_noar)), 
                      all=T)
  return(comp_dens)
  
}

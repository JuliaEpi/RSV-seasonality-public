func_xcorr <- function(codaobj) {
  
  corr <- crosscorr(codaobj)
  corr[lower.tri(corr, diag=TRUE)] <- NA
  corr <- data.frame(corr)
  colnames(corr) <- rownames(corr)
  corr$param1 <- rownames(corr)
  corr$order <- 1:nrow(corr)
  corr <- corr %>% pivot_longer(1:(ncol(corr)-2))
  corr <- corr[!is.na(corr$value),]
  colnames(corr)[3:4] <- c("param2", "r")
  
  return(corr)
}



func_compcorr <- function(chains, chains_2019, chains_cos, chains_noar) {
  
  corr <- func_xcorr(chains)
  corr$model <- "Von Mises 2010-2022"
  corr_2019 <- func_xcorr(chains_2019)
  corr_2019$model <- "Von Mises 2010-2019"
  corr_cos <- func_xcorr(chains_cos)
  corr_cos$model <- "Cosine 2010-2022"
  corr_noar <- func_xcorr(chains_noar)
  corr_noar$model <- "Von Mises 2010-2022, no AR"
  
  corr_all <- merge(corr, corr_2019, 
                    by=intersect(names(corr), names(corr_2019)), 
                    all=T)
  
  corr_all <- merge(corr_all, corr_cos, 
                    by=intersect(names(corr_all), names(corr_cos)), 
                    all=T)
  
  corr_all <- merge(corr_all, corr_noar, 
                    by=intersect(names(corr_all), names(corr_noar)), 
                    all=T)
  
  return(corr_all)
  
}



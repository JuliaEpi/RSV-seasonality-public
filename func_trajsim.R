# Calculates quantiles from all sampled trajectories
trajsim_quantiles <- function(trajsim, var) {
  
  quantiles <- plyr::ddply(.data=trajsim, .variables="time", 
                           function(x) quantile(x[, var], 
                                                prob = c(0.025, 0.5, 0.975), na.rm=T))
  colnames(quantiles)[2:4] <- c("low95CI", "median", "up95CI")
  
  quantiles <- merge(unique(trajsim[,c("time", "date", "week", "month", "doy", "year", "npos")]), 
                     quantiles, by="time", all=T)
  
  return(quantiles)
  
}

# This function loads the trajecories simulated from the posteriors and produces
# all the datasets derived from it. It is done in a function to save memory

func_trajsim <- function(file, data, chainno=NULL) {
  
  trajsim <- read_csv(paste0("output/", file, ".csv"))
  
  #subset to specified chains
  if (!is.null(chainno)) {
    trajsim <- trajsim[trajsim$chain %in% chainno,]
  }
  
  trajsim <- merge(data[,c("time", "date", "week", "month", "year", "doy", "npos")], 
                   trajsim, by="time", all=T)
  colnames(trajsim) <- str_replace_all(colnames(trajsim), "¹", "\\_1")
  colnames(trajsim) <- str_replace_all(colnames(trajsim), "²", "\\_2")

  fit <- trajsim_quantiles(trajsim, "simobserror")
  
  reporting <- trajsim_quantiles(trajsim, "rhopct")
  
  beta <- trajsim %>% 
          filter(year==2011 & !is.na(date)) %>% 
          trajsim_quantiles(., "betaseason") %>% 
          dplyr::mutate(pct_bl_median = median*100/min(median)-100,
                        pct_bl_low95 = low95CI*100/min(low95CI)-100,
                        pct_bl_up95 = up95CI*100/min(up95CI)-100)
  
  traj_year <- trajsim %>% 
                dplyr::filter(year==2011) %>% 
                select(c("replicate", "doy",  matches("^[SEIR][0-3].*?"), matches("inc[0-9]"))) %>% 
                pivot_longer(cols=c(matches("^[SEIR][0-3].*?"), matches("inc[0-9]")), names_to = "state", values_to = "n") %>% 
                group_by(doy, replicate) %>% 
                dplyr::mutate(order=row_number())
              
  out <- list("fit" = fit,
              "reporting" = reporting,
              "beta" = beta,
              "traj_year"= traj_year)
  return(out)
  
}

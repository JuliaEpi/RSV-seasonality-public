func_calclevels <- function(trajsim) {
  
  # Calculate N and proportion of population in each state
  pop_state <- trajsim$traj_year %>% 
    filter(grepl(paste(c("^S", "^E", "^I", "^R"), collapse = "|"), state)) %>% 
    mutate(level = as.numeric(substr(state,2,2))) %>% 
    dplyr::group_by(doy, replicate) %>%
    arrange(order) %>% 
    dplyr::mutate(pct_pop = n*100/popsize,
                  age = cumsum(80*pct_pop/100)) 
  
  # Calculate N and proportion of population in each level
  pop_levels <- pop_state %>% 
    filter(doy==364.0) %>% 
    group_by(replicate, level) %>% 
    dplyr::summarise(state="N",
                     pop=sum(n))
  
  
  # Level-specific incidence and attack rates
  ar_data <- trajsim$traj_year %>% 
                filter(grepl("inc[0-3]", state)) %>% 
                group_by(replicate, state) %>% 
                dplyr::summarise(inc=sum(n)) %>% 
                dplyr::mutate(level = as.numeric(gsub("(inc)([0-3])", "\\2", state))) %>% 
                merge(., pop_levels[,c("replicate", "level", "pop")], 
                      by=c("replicate", "level"), all=T) %>% 
                mutate(ar = inc/pop) 
    
   ar_summary <- ar_data %>%  
                group_by(level, state) %>% 
                dplyr::summarise(quantile = c("low95CI", "median","up95CI"),
                                 inc = quantile(inc, c(0.025, 0.5, 0.975)),
                                 ar = quantile(ar, c(0.025, 0.5, 0.975))) %>% 
                dplyr::mutate(ar_pct = round(ar * 100)) %>% 
                pivot_wider(id_cols=c("state"), 
                            names_from = "quantile", values_from = c("ar_pct"))
              
              
  # Average age at first, second and third infection
  # (assuming 80 years of lifeexpectancy) (summarised across doy)
  age_at_inf <- pop_state %>% 
              filter(state %in% c("I0", "I1", "I2")) %>% 
              group_by(state) %>% 
              dplyr::summarise(quantile = c("low95CI", "median","up95CI"),
                               age = quantile(age, c(0.025, 0.5, 0.975))) %>% 
              pivot_wider(id_cols=c("state"), 
                          names_from = "quantile", values_from = "age")
            
  
  return(list(ar_data = ar_data, ar_summary=ar_summary, age=age_at_inf))
  
}

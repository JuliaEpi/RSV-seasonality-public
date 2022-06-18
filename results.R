# Script for
# "The magnitude of the seasonal forcing of RSV and implications for vaccination strategies"
# Author: Fabienne Krauer et al
# Date last updated: June 18, 2022


# Housekeeping ----------------------------------------------------------------

library(tidyverse)
library(coda)
library(BayesianTools)
library(viridis)
library(grid)
library(gridExtra)
library(readxl)
library(plyr)
library(scales)
library(stringr)

source("func_tocoda.R")
source("func_xcorr.R")
source("func_calcdensity.R")
source("func_trajsim.R")
source("func_calclevels.R")


theme_set(theme_minimal())
fignosize <- 18

popsize <- 8164000
ar_num <- c(33, 52, (73+38+9))

version <- "2022-02-06" #"2021-04-11"
models <- c("mises_2010_2022", "mises_2010_2019", "cosine_2010_2022", "mises_2010_2022_noar")

# Data ----------------------------------------------------------------

# RSV surveillance data
data <- readRDS(paste0("data/data_NSW_RSV_weekly_", version,".rds"))
data$month <- format(data$date, "%B")
data$order <- match(data$month, month.name)
data$monthabb <- month.abb[data$order]

# Bronchiolitis data
bronch <- read_excel("data/bronchiolitis.xlsx")

# Posterior chains
chains <- read_csv("output/trace_final_mises_2010_2022.csv")
chains <- chains[order(chains$iteration),]
chains <- func_tocoda(chains)

chains_2019 <- read_csv("output/trace_final_mises_2010_2019.csv")
chains_2019 <- chains_2019[order(chains_2019$iteration),]
# Select only chains from the set with the higher log posterior density
chains_2019 <- chains_2019 %>% group_by(chain) %>% 
                dplyr::mutate(max_lp=round(max(lp))) %>% 
                dplyr::ungroup() %>% dplyr::filter(lp > min(max_lp)) %>% 
                select(-max_lp)
chainno_2019 <- sort(unique(chains_2019$chain))
chains_2019 <- func_tocoda(chains_2019)

chains_cos <- read_csv("output/trace_final_cosine_2010_2022.csv")
chains_cos <- chains_cos[order(chains_cos$iteration),]
chains_cos <- func_tocoda(chains_cos)

chains_noar <- read_csv("output/trace_final_mises_2010_2022_noar.csv")
chains_noar <- chains_noar[order(chains_noar$iteration),]
# Select only chains from the set with the higher log posterior density
chains_noar <- chains_noar %>% group_by(chain) %>% 
  dplyr::mutate(max_lp=round(max(lp))) %>% 
  dplyr::ungroup() %>% dplyr::filter(lp > min(max_lp)) %>% 
  select(-max_lp)
chains_noar <- sort(unique(chains_noar$chain))
chains_noar <- func_tocoda(chains_noar)


# Posterior predictive for time series data (trajectory)
trajsim <- func_trajsim("traj_sim_mises_2010_2022", data)
trajsim_cos <- func_trajsim("traj_sim_cosine_2010_2022", data)
trajsim_2019 <- func_trajsim("traj_sim_mises_2010_2019", data, chainno_2019)
trajsim_noar <- func_trajsim("traj_sim_mises_2010_2022_noar", data)

# Posterior predictive for AR data
ar <- read_csv("output/ar_sim_mises_2010_2022.csv")
ar <- pivot_longer(ar, cols=1:3, names_to = "level", 
                   values_to = "num")
ar$truth <- ifelse(ar$level=="ar1", ar_num[1], ifelse(ar$level=="ar2", ar_num[2], ar_num[3]))
ar$level <- recode(ar$level, "ar1"="level 0", "ar2" = "level 1/2", "ar3"= "level 3")

# Posterior quantiles
posteriors <- read_csv("output/posterior_quantiles_mises_2010_2022.csv")
#posteriors_2019 <- read_csv("output/posterior_quantiles_mises_2010_2019.csv")
#posteriors_cos <- read_csv("output/posterior_quantiles_cosine_2010_2022.csv")
#posteriors_noar <- read_csv("output/posterior_quantiles_mises_2010_2022_noar.csv")

# Simulated vaccination scenarios
vaccsim <- read_csv("output/vaccsim.csv")
vaccsim <- vaccsim %>% 
  mutate(vaccstart = vaccstart + 8) %>% # correct vaccstart + 7 days such that vaccstart=1 corresponds to doy=1
  group_by(coverage, vaccstart, duration, strategy) %>% 
    dplyr::arrange(time) %>% 
    dplyr::mutate(doy = rep(1:365, 10),
                  year = rep(1:10, each=365),
                  month = format(as.Date(doy, origin="2000-01-01"), "%B"))

# Base case scenario (no vaccination)
base <- read_csv("output/vaccsim_base.csv")
base <- base %>% 
  dplyr::arrange(time) %>% 
  dplyr::mutate(vaccstart = vaccstart + 8,
                doy = rep(1:365, 10),
                year = rep(1:10, each=365),
                month = format(as.Date(doy, origin="2000-01-01"), "%B"))


# Google mobility data
mobility_full <- readRDS(paste0("data/mobility_full_", version, ".rds"))
mobility_fit <- read.csv(paste0("data/mobility_fit_", version ,".csv"))
mobility_fit$date <- as.Date(mobility_fit$date)

# other data
studies <- read_excel("data/studies.xlsx")
stringency <- readRDS("data/stringency.rds")
contacts <- readRDS("data/contacts.rds")


# Results ----------------------------------------------------------------

# Duration of Season
data_my <- data %>% dplyr::group_by(year, month, order, monthabb) %>% 
  dplyr::filter(!is.na(npos)) %>% 
  dplyr::summarise(cases = sum(npos)) %>% 
  dplyr::ungroup() %>% dplyr::group_by(year) %>% 
  dplyr::arrange(order) %>% 
  dplyr::mutate(casespct = cases*100/sum(cases),
                casespct_cumul = cumsum(casespct))

data_my %>% dplyr::filter(order %in% c(4,5,6,7,8,9) & year<2020) %>% 
  dplyr::group_by(year) %>% 
  dplyr::summarise(casespct = sum(casespct))

ggplot(data_my[data_my$year<2020,]) + 
  geom_point(aes(x=reorder(monthabb, order), y=casespct_cumul, group=as.factor(year), colour=as.factor(year))) +
  geom_hline(yintercept=50-(75/2)) +
  geom_hline(yintercept=50+(75/2)) 

# Peak week
data_y <- data %>% dplyr::group_by(year) %>% 
  dplyr::filter(!is.na(npos)) %>% 
  dplyr::summarise(peakmonth = month[npos==max(npos)],
                    peakweek = isoweek[npos==max(npos)],
                   peakdate = date[npos==max(npos)])
summary(data_y$peakweek[data_y$year<2020])
unique(data_y$peakmonth[data_y$year<2020])


# Fitting results  --------------------------------------------------

round(posteriors[4,c(2,4,6)],2) # baseline beta0
round(posteriors[5,c(2,4,6)],2) # eta
round(posteriors[5,c(2,4,6)] + posteriors[4,c(2,4,6)],2) # maximum forcing
max(trajsim$beta$median)/min(trajsim$beta$median) # Ratio of seasonal forcing: max over min
max(trajsim$beta$low95CI)/min(trajsim$beta$low95CI) # Ratio of seasonal forcing: max over min
max(trajsim$beta$up95CI)/min(trajsim$beta$up95CI) # Ratio of seasonal forcing: max over min

trajsim$beta$doy[trajsim$beta$median==max(trajsim$beta$median)] # beta_eff peak time
trajsim$beta$week[trajsim$beta$median==max(trajsim$beta$median)] # beta_eff peak time
trajsim$beta$month[trajsim$beta$median==max(trajsim$beta$median)] # beta_eff peak time

posteriors[6,c(2,4,6)] # omega

round(posteriors[1,c(2,4,6)],4) # rho0
posteriors[3,c(2,4,6)] # q
round(trajsim$reporting[1,8:10],2) # % reported at start and end of time series
round(trajsim$reporting[588,8:10],2)

round(posteriors[2,c(2,4,6)],2) # psi
round(posteriors[8,c(2,4,6)],2) # k
round(posteriors[7,c(2,4,6)] + 7) # phi

# reduction of susceptibility
round(posteriors[9,c(2,4,6)],2) # level 1
round(posteriors[10,c(2,4,6)],2) # delta 2
round(posteriors[11,c(2,4,6)],2) # delta 3

round(posteriors[9,c(2,4,6)] * posteriors[10,c(2,4,6)],2) # level 2
round(posteriors[9,c(2,4,6)] * posteriors[10,c(2,4,6)] * posteriors[11,c(2,4,6)],2) # level 3


# Calculate level-specific attack rates and age-at-infection
by_level <- func_calclevels(trajsim)
by_level$ar_summary
by_level$age

# Total yearly incidence and proportion of population
trajsim$traj_year %>% 
  filter(grepl("inc", state)) %>% 
  group_by(replicate) %>% 
  dplyr::summarise(n=sum(n)) %>% 
  ungroup() %>% 
  dplyr::summarise(quantile = c("low95CI", "median","up95CI"),
                   n = quantile(n, c(0.025, 0.5, 0.975)),
                   pct = n*100/popsize) %>% 
  pivot_wider(names_from = "quantile", values_from = c("n", "pct"))



# Normalised observed pre-pandemic incidence
inc_norm <- trajsim$fit %>% 
  filter(year<2020) %>% 
  dplyr::group_by(year) %>% 
  dplyr::mutate(npos_norm = (npos - min(npos, na.rm = T))/(max(npos, na.rm = T) - min(npos, na.rm = T)),
                median_norm = (median - min(median, na.rm = T))/(max(median, na.rm = T) - min(median, na.rm = T)))


# Fig 1 ====================================================================

# 1a
fig1a <- ggplot(trajsim$fit) + 
  geom_point(aes(x=date, y=npos)) + 
  geom_line(aes(x=date, y=median), color="red") + 
  geom_ribbon(aes(x=date, ymin=low95CI, ymax=up95CI), alpha=0.3, fill="red") +
  geom_vline(xintercept=as.Date("2020-03-16"), color="blue") + 
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme(panel.grid.minor.x = element_blank()) +
  ylab("\nDetected RSV cases") + xlab(NULL)

fig1a

fig1b <- ggplot(inc_norm) + 
  geom_point(aes(x=doy, y=npos_norm), colour="grey") +
  geom_line(aes(x=doy, y=median_norm, group=as.factor(year)), colour="red") +
  ylab("\nnormalised incidence") + xlab("day of the year")
fig1b 
  
fig1c <- ggplot(trajsim$beta) + 
  geom_line(aes(x=week, y=median), color="turquoise4") + 
  geom_ribbon(aes(x=week, ymin=low95CI, ymax=up95CI), alpha=0.3, fill="turquoise4") +
  ylab(~paste(trajsim$beta[eff])) + xlab("week") + 
  coord_cartesian(ylim=c(0.0,4))

fig1c

  
fig1d <- ggplot(trajsim$reporting) + 
  geom_line(aes(x=date, y=median), color="blue") + 
  geom_ribbon(aes(x=date, ymin=low95CI, ymax=up95CI), alpha=0.3, fill="blue") +
  ylab("\n% symptomatic reported") + xlab(NULL)

fig1d

fig1e <- ggplot(ar) +
        geom_histogram(aes(x=num), binwidth=1, alpha=0.5, fill="blue") + 
        geom_vline(aes(xintercept=truth), colour="red") +
        facet_wrap(~level, scale="free", nrow=1) +
        xlab("numerator of attack rate") + ylab("density")
fig1e

# Combine
fig1 <- grid.arrange(arrangeGrob(fig1a, top=textGrob("A", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                     just=c("left", "top"), gp=gpar(col="black", fontsize=fignosize))), 
                     arrangeGrob(fig1b, top=textGrob("B", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                     just=c("left", "top"), gp=gpar(col="black",fontsize=fignosize))),
                     arrangeGrob(fig1c, top=textGrob("C", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                     just=c("left", "top"), gp=gpar(col="black",fontsize=fignosize))),
                     arrangeGrob(fig1d, top=textGrob("D", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                     just=c("left", "top"), gp=gpar(col="black",fontsize=fignosize))),
                     arrangeGrob(fig1e, top=textGrob("E", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                     just=c("left", "top"), gp=gpar(col="black",fontsize=fignosize))),
                     
                     
                     layout_matrix=rbind(c(1,1,1,1,1,1), 
                                         c(1,1,1,1,1,1),
                                         c(1,1,1,1,1,1),
                                         c(2,2,3,3,4,4), 
                                         c(2,2,3,3,4,4),
                                         c(5,5,5,5,5,5),
                                         c(5,5,5,5,5,5)))

ggsave("output/fig1.tiff", plot=fig1, units="cm", width=25, height=20, dpi=300)
ggsave("output/fig1.png", plot=fig1, units="cm", width=25, height=20, dpi=96)


# Vaccination scenarios ------------------------------------------------------

# Reshape to long
base_long <- base %>% 
  filter(time >= max(base$time)-364) %>% 
  select(time, doy, year, month, matches("inc")) %>% 
  pivot_longer(cols=grep("_", colnames(.)), names_to = "state", values_to = "n")  %>% 
  tidyr::separate(state, c("state", "stats"), "_") %>% 
  pivot_wider(names_from = "stats", values_from = "n")
colnames(base_long)[6:8] <- paste0(colnames(base_long)[6:8], "_base")

vaccsim_long <- vaccsim %>% 
  filter(year==10) %>% 
  select(matches(c("strategy", "^time$", "^year$", "^doy$", "^month$", "coverage", "vaccstart", "duration", "vaccswitch", "^inc", "^nvacc"))) %>% 
  pivot_longer(cols = grep("_", colnames(.)), names_to = "state", values_to = "n") %>% 
  separate(state, c("state", "stats"), "_") %>% 
  pivot_wider(names_from = "stats", values_from = "n")
colnames(vaccsim_long)[11:13] <- paste0(colnames(vaccsim_long)[11:13], "_vacc")

# Combine base and vacc scenarios
sim_y10 <- merge(vaccsim_long, base_long, 
                 by=c("time", "state", "doy", "month", "year"), all=T)

# Make helper vars for start and stop of vaccination window
sim_y10 <- sim_y10 %>% group_by(strategy, coverage, vaccstart, duration, state) %>% 
  dplyr::mutate(window_start = min(doy[vaccswitch==1]),
                window_end = max(doy[vaccswitch==1]))


# Calculate the % shift in amplitude and the shift in peak timing and yearly cumulative incidence
sim_agg <- sim_y10 %>% 
  group_by(strategy, coverage, vaccstart, duration, state) %>% 
  dplyr::mutate(tpeak_vacc = min(doy[median_vacc==max(median_vacc)]),
                tpeak_base = min(doy[median_base==max(median_base)]),
                npeak_vacc = max(median_vacc),
                npeak_base = max(median_base)) %>% 
  dplyr::summarise(tpeak_vacc = tpeak_vacc[1],
                   tpeak_base = tpeak_base[1],
                   npeak_vacc = npeak_vacc[1],
                   npeak_base = npeak_base[1],
                   ncumul_vacc = sum(median_vacc),
                   ncumul_base = sum(median_base),
                   cases_averted = ncumul_base - ncumul_vacc,
                   cases_averted_pct = cases_averted*100/ncumul_base,
                   peakdiff_pct = npeak_vacc*100/npeak_base,
                   peakshift_d = tpeak_vacc - tpeak_base) 


# maximum increase by level
sim_agg %>% 
  filter(state!="nvacc") %>% 
  filter((duration==210 & vaccstart==61 & strategy=="seasonal") | (strategy=="continuous")) %>% 
  dplyr::group_by(strategy, state) %>% 
  dplyr::summarise(increase = max(peakdiff_pct, na.rm=T))


# maximum delay
sim_agg %>% 
  filter(state!="nvacc") %>% 
  filter((duration==210 & vaccstart==61 & strategy=="seasonal") | (strategy=="continuous")) %>% 
  dplyr::group_by(strategy, state) %>% 
  dplyr::summarise(max_shift = max(peakshift_d, na.rm=T),
                   min_shift = min(peakshift_d, na.rm=T))


# Fig 2 ===============================================================

# Seasonal dynamics for continuous and seasonal vaccination compared to baseline
labeller <- as_labeller(c("continuous" =  "continuous", "seasonal" = "seasonal", "inc0" = "level 0", "inc1" = "level 1", "inc2" = "level 2", "inc3" = "level 3"))

fig2 <- sim_y10 %>% 
        filter(!(state %in% c("nvacc", "inc"))) %>% 
        filter(coverage %in% c(0.2, 0.4, 0.6, 0.8, 1.0)) %>% 
        filter((strategy=="seasonal" & duration==210 & vaccstart==61) | (strategy=="continuous" & is.na(duration) & vaccstart==1)) %>% 
        ggplot() +
        scale_color_viridis_d(direction=-1) +
        scale_fill_viridis_d(direction=-1) +
        geom_rect(aes(xmin=window_start, xmax=window_end, ymin=0, ymax=Inf), alpha=0.01, fill="grey") +
        geom_line(aes(x=doy, y=median_vacc, colour=as.factor(coverage))) +
        geom_ribbon(aes(x=doy, ymin=low95CI_vacc, ymax=up95CI_vacc, fill=as.factor(coverage)), alpha=0.3, show.legend=F) +
        geom_line(aes(x=doy, y=median_base), colour="black",  linetype="dashed") +
        facet_grid(rows = state ~ strategy, scales="free_y", labeller = labeller) +
        xlab("day of the year") + ylab("daily incidence") + 
        labs(colour="coverage") +
        theme(legend.position = "bottom")
fig2  

ggsave("output/fig2.tiff", plot=fig2, units="cm", width=15, height=15, dpi=300)
ggsave("output/fig2.png", plot=fig2, units="cm", width=15, height=15, dpi=96)



# Fig 3 ===============================================================


# Effect of vaccination coverage on increases in Amplitude and tshift at year 10
fig3a <- sim_agg %>% 
        filter(state!="nvacc") %>% 
        filter((duration==210 & vaccstart==61 & strategy=="seasonal") | (strategy=="continuous")) %>% 
        ggplot(.) +
        facet_wrap(~strategy) + theme_light() +
        geom_point(aes(x=coverage, y=peakdiff_pct, color=as.factor(state))) +
        scale_color_viridis_d() +
        geom_hline(aes(yintercept=100), color="red") +
        theme(legend.position = "right") + labs(color=NULL) +
        ylab("% incidence at peak")
fig3a

fig3b <- sim_agg %>% 
  filter(state!="nvacc") %>% 
  filter((duration==210 & vaccstart==61 & strategy=="seasonal") | (strategy=="continuous")) %>% 
  ggplot(.) +
  facet_wrap(~strategy) + theme_light() +
  geom_point(aes(x=coverage, y=peakshift_d, color=as.factor(state))) +
  geom_hline(aes(yintercept=0), color="red") +
  theme(legend.position = "right") + labs(color=NULL) +
  scale_color_viridis_d() +
  ylab("delay in epidemic peak (days)")
fig3b

fig3 <- grid.arrange(arrangeGrob(fig3a, top=textGrob("A", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                     just=c("left", "top"), gp=gpar(col="black", fontsize=fignosize))), 
                     arrangeGrob(fig3b, top=textGrob("B", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                     just=c("left", "top"), gp=gpar(col="black",fontsize=fignosize))),
                     
                     nrow=2)

ggsave("output/fig3.tiff", plot=fig3, units="cm", width=15, height=18, dpi=300)
ggsave("output/fig3.png", plot=fig3, units="cm", width=15, height=18, dpi=96)

# Fig 4 ===============================================================

fig4 <- sim_agg %>% 
  filter(state!="nvacc") %>% 
  filter((duration==210 & vaccstart==61 & strategy=="seasonal") | (strategy=="continuous")) %>% 
  ggplot() + ylab("% cases averted") + theme_light() +
  geom_bar(aes(x=coverage, y=cases_averted_pct), stat="identity", alpha=0.5, fill="#482677FF", colour="#482677FF") +
  facet_grid(rows = state ~ strategy, labeller = labeller) 
fig4

ggsave("output/fig4.tiff", plot=fig4, units="cm", width=15, height=18, dpi=300)
ggsave("output/fig4.png", plot=fig4, units="cm", width=15, height=18, dpi=96)



# Fig S1 ==================================================================

# RSV cases and influenza samples
figS1a <- ggplot(data) + 
  geom_line(aes(x=date, y=npos)) + 
  geom_line(aes(x=date, y=ntested_influenza/10), color="blue") + 
  xlab(NULL) + ylab("RSV cases") +
  scale_y_continuous(sec.axis = sec_axis(~.*10, 
                                         name="samples tested for influenza")) +
  geom_vline(xintercept=as.Date("2020-03-16"), color="red") +
  theme(axis.text.y.right = element_text(color = "blue"), 
        axis.title.y.right = element_text(color = "blue"))
figS1a    


# Bronchiolitis cases
figS1b <- ggplot(bronch) + 
  geom_line(aes(x=month, y=cases, color=as.factor(year))) +
  ylab("weekly bronchiolitis visits to ED") + 
  scale_x_continuous(breaks=c(1,7,12), labels=month.abb[c(1,7,12)]) +
  labs(color="year") + xlab(NULL)
figS1b


# Percentage distribution by month and year (AAP)
figS1c <- ggplot(data_my[data_my$year<2020,]) + 
  geom_tile(aes(x=reorder(monthabb, order), y=as.character(year), fill=casespct)) +
  xlab(NULL) + ylab(NULL) + labs(fill="% cases") +
  scale_fill_viridis() +
  geom_vline(aes(xintercept=3.5), colour="red")  +
  geom_vline(aes(xintercept=8.5), colour="red")

figS1c


# Pre- vs. pandemic RSV seasonality
figS1d <- ggplot() +
  geom_line(data=data[data$year<2020,], aes(x=doy, y=npos, group=as.factor(year), color="pre-pandemic")) + 
  geom_line(data=data[data$year>=2020,], aes(x=doy, y=npos, group=as.factor(year), color="pandemic")) +
  labs(color=NULL) + ylab("RSV cases") + xlab("day of the year") +
  theme(legend.position = c(0.2,0.8))

figS1d

# combine plots
figS1 <- grid.arrange(arrangeGrob(figS1a, top=textGrob("A", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                       just=c("left", "top"), gp=gpar(col="black", fontsize=fignosize))), 
                      arrangeGrob(figS1b, top=textGrob("B", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                       just=c("left", "top"), gp=gpar(col="black",fontsize=fignosize))),
                      arrangeGrob(figS1c, top=textGrob("C", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                       just=c("left", "top"), gp=gpar(col="black", fontsize=fignosize))), 
                      arrangeGrob(figS1d, top=textGrob("D", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                       just=c("left", "top"), gp=gpar(col="black",fontsize=fignosize))),
                      
                      ncol=2)

ggsave("output/figS1.tiff", plot=figS1, units="cm", width=22, height=18, dpi=300)
ggsave("output/figS1.png", plot=figS1, units="cm", width=22, height=18, dpi=96)


# Fig S2 =====================================================================

comparison_work <- merge(stringency[stringency$category=="c2_workplace_closing",c("date", "index")],
                         mobility_fit[,c("date", "scaling")],
                         by="date", all=T)

comparison_school <- merge(stringency[stringency$category=="c1_school_closing",c("date", "index")],
                           mobility_fit[,c("date", "scaling")],
                           by="date", all=T)

figS2a <- ggplot(mobility_full) +
  geom_point(aes(x=date, y=index), alpha=0.1) +
  geom_line(aes(x=date, y=fit), color="red") +
  scale_x_date(labels = date_format("%m-%Y")) +
  theme(plot.margin = margin(10, 10, 10, 10)) +
  xlab(NULL) + 
  ylab("Google mobility index")

figS2a

figS2b <- ggplot() + 
  geom_line(data=contacts, aes(x=date, y=contacts_scaled, color="contact survey")) +
  geom_line(data=mobility_fit[mobility_fit$date<=max(contacts$date),], aes(x=date, y=scaling, color="mobility data")) +
  labs(color="data source") + 
  scale_x_date(labels = date_format("%m-%Y")) +
  ylab("relative contact reduction") + xlab(NULL) +
  theme(legend.position = c(0.8, 0.2),
        plot.margin = margin(10, 10, 10, 10)) + 
  scale_y_continuous(limits=c(0,1))
figS2b


figS2c <- ggplot(comparison_work[!is.na(comparison_work$index),]) + 
  geom_boxplot(aes(x=as.factor(index), y=scaling)) +
  xlab("Oxford stringency index: work places") + 
  ylab("Google mobility (work places)") +
  theme(plot.margin = margin(10, 10, 10, 10))

figS2c

figS2d <- ggplot(comparison_school[!is.na(comparison_work$index),]) + 
  geom_boxplot(aes(x=as.factor(index), y=scaling)) +
  xlab("Oxford stringency index: schools") + 
  ylab("Google mobility (work places)") +
  theme(plot.margin = margin(10, 10, 10, 10))

figS2d


figS2 <- grid.arrange(arrangeGrob(figS2a, top=textGrob("A", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                       just=c("left", "top"), gp=gpar(col="black", fontsize=fignosize))), 
                      arrangeGrob(figS2b, top=textGrob("B", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                       just=c("left", "top"), gp=gpar(col="black",fontsize=fignosize))),
                      arrangeGrob(figS2c, top=textGrob("C", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                       just=c("left", "top"), gp=gpar(col="black", fontsize=fignosize))), 
                      arrangeGrob(figS2d, top=textGrob("D", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                       just=c("left", "top"), gp=gpar(col="black", fontsize=fignosize))), 
                      
                      layout_matrix=rbind(c(1,2), c(3,4)))


ggsave("output/figS2.tiff", plot=figS2, units="cm", width=19, height=19, dpi=300)
ggsave("output/figS2.png", plot=figS2, units="cm", width=19, height=19, dpi=96)


# Fig S3  ===============================================================

# generated in postprocessing.jl

# Fig S4  ===============================================================

# generated in postprocessing.jl


# Sensitivity  ===============================================================


# Compare cross correlations
corr <- func_compcorr(chains, chains_2019, chains_cos, chains_noar)

# Compare posterior densities
dens <- func_compdensity(chains, chains_cos, chains_2019, chains_noar)

# Compare fit
comp_fit <- merge(trajsim_cos$fit, 
                  trajsim_2019$fit[,c("time", "low95CI", "median", "up95CI")], 
                  by="time", all=T)
comp_fit <- merge(comp_fit, 
                  trajsim_noar$fit[,c("time", "low95CI", "median", "up95CI")], 
                  by="time", all=T)
colnames(comp_fit) <- gsub("\\.x", "_cosine2022", colnames(comp_fit))
colnames(comp_fit) <- gsub("\\.y", "_mises2019", colnames(comp_fit))
colnames(comp_fit)[(ncol(comp_fit)-2):ncol(comp_fit)] <- paste0(colnames(comp_fit)[(ncol(comp_fit)-2):ncol(comp_fit)], "_misesNoar")
comp_fit <- merge(comp_fit, trajsim$fit[,c("time", "low95CI", "median", "up95CI")], 
                  by="time", all=T)

# compare beta_eff
comp_beta <- merge(trajsim$beta, 
                   trajsim_cos$beta[,c("time", "low95CI", "median", "up95CI")],
                   by="time", all=T)
colnames(comp_beta) <- gsub("\\.x", "_mises", colnames(comp_beta))
colnames(comp_beta) <- gsub("\\.y", "_cosine", colnames(comp_beta))

# Calculate level-specific attack rates for model not fitted on AR
by_level_noar <- func_calclevels(trajsim_noar)
by_level_noar$ar_summary

comp_ar <- merge(by_level$ar_data[,c("replicate", "state", "level", "ar")], 
                 by_level_noar$ar_data[,c("replicate", "state", "level", "ar")], 
                 by=c("replicate", "state", "level"), all=T)
colnames(comp_ar) <- gsub("\\.x", "_mises", colnames(comp_ar))
colnames(comp_ar) <- gsub("\\.y", "_misesNoAR", colnames(comp_ar))
comp_ar <- pivot_longer(comp_ar, cols=c(4:5), names_to = "model", values_to = "ar")
comp_ar$model <- recode(comp_ar$model, "ar_mises"="Von Mises", "ar_misesNoAR"= "Von Mises, no AR")
comp_ar$ar <- comp_ar$ar*100


# Fig S5 ===============================================================

figS5A <- ggplot(comp_fit) + 
  geom_point(aes(x=date, y=npos)) + 
  geom_line(aes(x=date, y=median), color="black") + 
  geom_line(aes(x=date, y=low95CI), color="black", linetype="solid") +
  geom_line(aes(x=date, y=up95CI), color="black", linetype="solid") +
  geom_line(aes(x=date, y=median_mises2019), color="salmon") + 
  geom_ribbon(aes(x=date, ymin=low95CI_mises2019, ymax=up95CI_mises2019), alpha=0.5, fill="salmon") +
  ylab("Detected RSV cases") + xlab(NULL)

figS5A


#re-scale the limits to center the palette at 0
limit <- max(abs(corr$r[corr$model %in% c("Von Mises 2010-2022", "Von Mises 2010-2019")])) * c(-1, 1)

figS5B <- ggplot(corr[corr$model %in% c("Von Mises 2010-2022", "Von Mises 2010-2019"),]) + 
  facet_wrap(~model, ncol=1, dir="h") +
  geom_tile(aes(x=reorder(param1, order), y=reorder(param2,order), fill=r)) +
  geom_tile(data=corr[(corr$r>0.7 | corr$r < -0.7) & corr$model %in% c("Von Mises 2010-2022", "Von Mises 2010-2019"),], 
            aes(x=reorder(param1, order), y=reorder(param2,order)), fill=NA, colour="red") +
  #scale_fill_viridis() +
  scale_fill_distiller(type = "div", limit = limit, palette="RdBu") + #PiYG
  #scale_fill_gradient2() +
  ylab(NULL) + xlab(NULL) +
  scale_x_discrete(labels=c("rho0", "psi", "q", "beta0", "eta", "omega", "phi", "k", "delta2", "delta3")) +
  scale_y_discrete(labels=c("psi", "q", "beta0", "eta", "omega", "phi", "k", "delta1", "delta2", "delta3")) +
  labs(fill="correlation") +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom")

figS5B

figS5C <- ggplot(dens[dens$model %in% c("Von Mises 2010-2022", "Von Mises 2010-2019"),]) +
  facet_wrap(~par, scale="free", ncol=3) + 
  geom_area(aes(x=x, y=y, fill=as.factor(model)), alpha=0.6) +
  labs(fill="model") +
  theme(legend.position = "bottom") + ylab("density") + xlab(NULL)
figS5C


# Combine
figS5 <- grid.arrange(arrangeGrob(figS5A, top=textGrob("A", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                       just=c("left", "top"), gp=gpar(col="black", fontsize=fignosize))), 
                      arrangeGrob(figS5B, top=textGrob("B", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                       just=c("left", "top"), gp=gpar(col="black",fontsize=fignosize))),
                      arrangeGrob(figS5C, top=textGrob("C", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                       just=c("left", "top"), gp=gpar(col="black",fontsize=fignosize))),
                      
                      
                      layout_matrix=rbind(c(1,1,1,1,1),
                                          c(1,1,1,1,1),
                                          c(2,2,3,3,3), 
                                          c(2,2,3,3,3), 
                                          c(2,2,3,3,3)))

figS5

ggsave("output/figS5.tiff", plot=figS5, units="cm", width=22, height=24, dpi=300)
ggsave("output/figS5.png", plot=figS5, units="cm", width=22, height=24, dpi=96)


# Fig S6 ===============================================================

figS6a <- ggplot(comp_fit) + 
  geom_point(aes(x=date, y=npos)) + 
  geom_line(aes(x=date, y=median), color="black") + 
  geom_line(aes(x=date, y=low95CI), color="black", linetype="solid") +
  geom_line(aes(x=date, y=up95CI), color="black", linetype="solid") +
  geom_line(aes(x=date, y=median_cosine2022), color="salmon") + 
  geom_ribbon(aes(x=date, ymin=low95CI_cosine2022, ymax=up95CI_cosine2022), alpha=0.5, fill="salmon") +
  ylab("Detected RSV cases") + xlab(NULL)

figS6a


figS6b <- ggplot(comp_beta) + 
  geom_line(aes(x=week, y=median_cosine), color="salmon") + 
  geom_ribbon(aes(x=week, ymin=low95CI_cosine, ymax=up95CI_cosine), alpha=0.5, fill="salmon") +
  geom_line(aes(x=week, y=median_mises), linetype="solid") + 
  geom_line(aes(x=week, y=low95CI_mises), linetype="dashed") + 
  geom_line(aes(x=week, y=up95CI_mises), linetype="dashed") + 
  ylab("seasonal beta") + xlab("week") 

figS6b


figS6c <- ggplot(dens[dens$model %in% c("Von Mises 2010-2022", "Cosine 2010-2022") & !(dens$par %in% c("beta0", "k", "eta")),]) +
  facet_wrap(~par, scale="free", ncol=3) + 
  geom_area(aes(x=x, y=y, fill=as.factor(model)), alpha=0.6) +
  labs(fill="model") +
  theme(legend.position = "bottom") + ylab("density") + xlab(NULL)
figS6c 


# Combine
figS6 <- grid.arrange(arrangeGrob(figS6a, top=textGrob("A", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                       just=c("left", "top"), gp=gpar(col="black", fontsize=fignosize))), 
                      arrangeGrob(figS6b, top=textGrob("B", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                       just=c("left", "top"), gp=gpar(col="black",fontsize=fignosize))),
                      arrangeGrob(figS6c, top=textGrob("C", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                       just=c("left", "top"), gp=gpar(col="black",fontsize=fignosize))),
                      
                      
                      layout_matrix=rbind(c(1,1,1,1,1,1),
                                          c(1,1,1,1,1,1),
                                          c(2,2,3,3,3,3), 
                                          c(2,2,3,3,3,3), 
                                          c(2,2,3,3,3,3)))

ggsave("output/figS6.tiff", plot=figS6, units="cm", width=22, height=20, dpi=300)
ggsave("output/figS6.png", plot=figS6, units="cm", width=22, height=20, dpi=96)



# Fig S7 ===============================================================

figS7A <- ggplot(comp_fit) + 
  geom_point(aes(x=date, y=npos)) + 
  geom_line(aes(x=date, y=median), color="black") + 
  geom_line(aes(x=date, y=low95CI), color="black", linetype="solid") +
  geom_line(aes(x=date, y=up95CI), color="black", linetype="solid") +
  geom_line(aes(x=date, y=median_misesNoar), color="salmon") + 
  geom_ribbon(aes(x=date, ymin=low95CI_misesNoar, ymax=up95CI_misesNoar), alpha=0.5, fill="salmon") +
  ylab("Detected RSV cases") + xlab(NULL)

figS7A

figS7B <- ggplot(comp_ar) +
  geom_boxplot(aes(x=as.factor(level), y=ar, fill=model)) +
  xlab("level") + ylab("attack rate (%)") +
  theme(legend.position = "bottom")
figS7B

#re-scale the limits to center the palette at 0
limit <- max(abs(corr$r[corr$model %in% c("Von Mises 2010-2022", "Von Mises 2010-2022, no AR")])) * c(-1, 1)

figS7C <- ggplot(corr[corr$model %in% c("Von Mises 2010-2022", "Von Mises 2010-2022, no AR"),]) + 
  facet_wrap(~model, ncol=1, dir="h") +
  geom_tile(aes(x=reorder(param1, order), y=reorder(param2,order), fill=r)) +
  geom_tile(data=corr[(corr$r>0.7 | corr$r < -0.7) & corr$model %in% c("Von Mises 2010-2022", "Von Mises 2010-2022, no AR"),], 
            aes(x=reorder(param1, order), y=reorder(param2,order)), fill=NA, colour="red") +
  #scale_fill_viridis() +
  scale_fill_distiller(type = "div", limit = limit, palette="RdBu") + #PiYG
  #scale_fill_gradient2() +
  ylab(NULL) + xlab(NULL) +
  scale_x_discrete(labels=c("rho0", "psi", "q", "beta0", "eta", "omega", "phi", "k", "delta2", "delta3")) +
  scale_y_discrete(labels=c("psi", "q", "beta0", "eta", "omega", "phi", "k", "delta1", "delta2", "delta3")) +
  labs(fill="correlation") +
  theme(axis.text.x = element_text(angle = 90))
        #legend.position = "bottom")

figS7C

figS7D <- ggplot(dens[dens$model %in% c("Von Mises 2010-2022", "Von Mises 2010-2022, no AR"),]) +
  facet_wrap(~par, scale="free", ncol=4) + 
  geom_area(aes(x=x, y=y, fill=as.factor(model)), alpha=0.6) +
  labs(fill="model") +
  theme(legend.position = "bottom") + ylab("density") + xlab(NULL)
figS7D



# Combine
figS7 <- grid.arrange(arrangeGrob(figS7A, top=textGrob("A", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                       just=c("left", "top"), gp=gpar(col="black", fontsize=fignosize))), 
                      arrangeGrob(figS7B, top=textGrob("B", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                       just=c("left", "top"), gp=gpar(col="black",fontsize=fignosize))),
                      arrangeGrob(figS7C, top=textGrob("C", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                       just=c("left", "top"), gp=gpar(col="black",fontsize=fignosize))),
                      arrangeGrob(figS7D, top=textGrob("D", x=unit(0, "npc"), y=unit(0, "npc"), 
                                                       just=c("left", "top"), gp=gpar(col="black",fontsize=fignosize))),
                      
                      
                      layout_matrix=rbind(c(1,1,1,1),
                                          c(1,1,1,1),
                                          c(2,2,3,3), 
                                          c(2,2,3,3), 
                                          c(2,2,3,3), 
                                          c(4,4,4,4),
                                          c(4,4,4,4),
                                          c(4,4,4,4)))

figS7

ggsave("output/figS7.tiff", plot=figS7, units="cm", width=22, height=24, dpi=300)
ggsave("output/figS7.png", plot=figS7, units="cm", width=22, height=24, dpi=96)



# Studies ===============================
studies$label <- paste(studies$ref, studies$country, sep=": ")
studies$avg_transrate_d <- ifelse(studies$unit=="per day", studies$avg_transrate,
                                  ifelse(studies$unit=="per year", studies$avg_transrate/365,
                                         ifelse(studies$unit=="per week", studies$avg_transrate/7, NA)))
studies$beta_min <-  ifelse(studies$forcing_func=="cosine" | studies$forcing_func=="sine",
                            (studies$avg_transrate_d - studies$avg_transrate_d*studies$amplitude), 
                            ifelse(studies$forcing_func=="cosine*", (studies$avg_transrate - studies$amplitude), NA))
studies$beta_max <-  ifelse(studies$forcing_func=="cosine" | studies$forcing_func=="sine",
                            (studies$avg_transrate_d + studies$avg_transrate_d*studies$amplitude), 
                            ifelse(studies$forcing_func=="cosine*", (studies$avg_transrate + studies$amplitude), NA))
studies$max_over_min <- studies$beta_max/studies$beta_min
studies$max_over_min <- ifelse(studies$ref=="Hodgson 2021",studies$amplitude, studies$max_over_min)
summary(studies$max_over_min)

# Summary estimates of duration of immunity and seasonal beta for comparison
max(studies$immunity_up[studies$source=="data"])
min(studies$immunity_low[studies$source=="data"])
studies$ref[studies$source=="data"]
max(studies$immunity_up[studies$source=="model"], na.rm=T)
min(studies$immunity_low[studies$source=="model"], na.rm=T)
max(studies$immunity_d[studies$source=="model"], na.rm=T)
min(studies$immunity_d[studies$source=="model"], na.rm=T)
studies$ref[studies$source=="model" & (!is.na(studies$immunity_d) | !is.na(studies$immunity_low))]

studies$ref[!is.na(studies$immunity_d)]


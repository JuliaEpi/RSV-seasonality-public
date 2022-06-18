# Script for
# "The magnitude of the seasonal forcing of RSV and implications for vaccination strategies"
# Author: Fabienne Krauer et al
# Date last updated: June 18, 2022

# Housekeeping ----------------------------------------------------------------

library(tidyverse)
library(readxl)
library(plyr)
library(grid)
library(gridExtra)
library(zoo)

theme_set(theme_minimal())

fignosize <- 16

version <- "2022-02-06" #"2021-04-11"


# RSV data ---------------------------------------------------------

rsv <- readRDS(paste0("data/data_NSW_RSV_weekly_", version,".rds"))
rsv <- rsv[order(rsv$date),]

# Hospitalization incidence (Homaira 2016)
# 0-1 year olds:
(4579 + 2869 + 289) / (178525 + 171202 + 617877)
# 1-2 year olds:
1643 / 583872
# 2-4 year olds:
385 / 1378996

# # Fit seasonal regression model to estimate coefficient q for exponentially increasing reporting rate
data <- rsv[rsv$date<=as.Date("2019-12-31"), c("date", "npos", "isoweek", "year")]
data <- data[order(data$date),]
data$time <- julian(data$date, orig=data$date[1]) #+1
data <- data[!is.na(data$npos),]

fit <- glm(npos ~ time + sin(2*pi*time/365) + cos(2*pi*time/365), family=poisson,
          data=data)
summary(fit)

data$predict <- exp(predict(fit))
ggplot(data) +
  geom_point(aes(x=time, y=npos)) +
  geom_line(aes(x=time, y=predict), color="red") +
  geom_line(aes(x=time, y=exp(fit$coefficients[1] + time*fit$coefficients[2])), color="blue")

# q parameter and 95% CI:
fit$coefficients[2]
fit$coefficients[2] + 1.96*summary(fit)$coefficients[2,2]
fit$coefficients[2] - 1.96*summary(fit)$coefficients[2,2]
# 0.000857155 [0.0008505084-0.0008638015]

data %>% group_by(year) %>% dplyr::filter(predict==max(predict))


# Google mobility NSW data -------------------------------------------
# data from: https://www.google.com/covid19/mobility/
# This may take 3 mins or so
#mobility <- read_csv("https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv")
mobility <- read_csv(paste0("data/Global_Mobility_Report_", version,".csv"))

colnames(mobility)[10:15] <- c("retail", "grocery", "parks", "transit", "work", "residential")
mobility <- mobility %>% 
            dplyr::filter(country_region_code == "AU" & 
                          sub_region_1=="New South Wales" & 
                          !is.na(sub_region_1) &
                          !is.na(date)) %>% 
            select(date, place_id, work)
colnames(mobility)[3] <- "index"
# limit to extent of RSV data and add time (t=0 at start of data, by definition)
mobility$time <- as.numeric(julian(mobility$date, origin=rsv$date[rsv$time==0])-2)
mobility <- mobility[mobility$date<=(max(rsv$date) + 3),] # + 3  add three days to avoid NA at the end due to padding from rolling average

# Fit 7-days rolling average
mobility_fit <- mobility %>% dplyr::group_by(date, time) %>%
                dplyr::summarise(index_median = mean(index, na.rm=T))
mobility_fit$fit <- rollmedian(mobility_fit$index_median, k=7, fill=NA)
mobility_fit$scaling <- mobility_fit$fit - max(mobility_fit$fit, na.rm=T)
mobility_fit$scaling <- 1+(mobility_fit$scaling/100)

mobility <- merge(mobility[,c("date", "index", "time")],
                  mobility_fit[!is.na(mobility_fit$fit),c("date", "time", "fit")],
                  by=c("date", "time"), all=T)

saveRDS(mobility, paste0("data/mobility_full_", version,".rds"))

mobility_fit <- mobility_fit[!is.na(mobility_fit$scaling),c("date", "time", "scaling")]
mobility_fit$scaling <- ifelse(mobility_fit$date<=as.Date("2020-02-29"),1.0, mobility_fit$scaling)
write.csv(mobility_fit, paste0("data/mobility_fit_", version ,".csv"), row.names = F, na="")
saveRDS(mobility_fit, paste0("data/mobility_fit_", version ,".rds"))


ggplot(mobility_fit) + 
  #geom_point(aes(x=date, y=index)) +
  geom_line(aes(x=date, y=scaling), color="red")


# NSW contact survey -----------------------------------------------------

contacts <- read_excel("data/NSW_contacts.xlsx")
contacts$month <- floor(contacts$time)
contacts$day <- contacts$time - contacts$month
contacts$day <- round(contacts$day*30)
contacts$date <- as.Date(paste0("2020-", contacts$month, "-", contacts$day))
contacts <- contacts[order(contacts$date),]
contacts$contacts_scaled <- contacts$contacts/max(contacts$contacts)
saveRDS(contacts, "data/contacts.rds")


# Estimate scaling parameter between google mobility data and contact survey data:
contact_mobility <- merge(contacts[,c("date", "scaled")], 
                        mobility_fit[,c("date", "scaling")], by="date", all.x=T)

contact_mobility$factor <- contact_mobility$scaled / contact_mobility$scaling
contact_mobility$factor <- ifelse(contact_mobility$factor>1,1,contact_mobility$factor)
ggplot(contact_mobility) + geom_line(aes(x=date, y=factor))


# Oxford stringency index -----------------------------------------------------
#(data source: https://github.com/OxCGRT/covid-policy-tracker/blob/master/data/timeseries/OxCGRT_timeseries_all.xlsx)
path <- "data/OxCGRT_timeseries_all.xlsx"
stringency <- vector("list", length(excel_sheets(path)))
for (i in excel_sheets(path)) {
  foo <- read_excel(path, sheet = i)
  foo <- foo[foo$country_name=="Australia",c(3:ncol(foo))]
  foo <- stack(foo)
  colnames(foo) <- c("index", "date")
  foo$category <- i
  stringency[[which(excel_sheets(path)==i)]] <- foo
}
stringency <- ldply(stringency, rbind)
stringency$date <- as.Date(stringency$date, "%d%B%Y")
sort(unique(stringency$category))
stringency <- stringency[!grepl("_flag$", stringency$category),]
stringency <- stringency[grepl("^c[0-9]_", stringency$category) | stringency$category=="stringency_index",] # keep only containment measures
# limit stringency data to extent of RSV data
stringency <- stringency[stringency$date<=max(rsv$date ,na.rm=T),]
stringency$label <- recode(stringency$category, 
                           "c1_school_closing" = "schools",
                           "c2_workplace_closing" = "workplace",
                           "c3_cancel_public_events" = "public events",
                           "c4_restrictions_on_gatherings" = "gatherings",
                           "c5_close_public_transport" = "public transport",
                           "c6_stay_at_home_requirements" = "stay-at-home",
                           "c7_movementrestrictions" = "internal movements",
                           "c8_internationaltravel" = "cross-border movements",
                           "stringency_index" = "stringency index")

saveRDS(stringency, "data/stringency.rds")
write.csv(stringency, "data/stringency.csv", row.names = F, na="")


# Plots
ggplot(stringency[stringency$category!="stringency_index",]) + 
  geom_tile(aes(x=date, y=label, fill=as.factor(index))) +
  scale_fill_manual(values=rev(c("#d53e4f","#f46d43","#fdae61","#fee08b","#e6f598")),na.value = "grey90") +
  ylab(NULL) + xlab(NULL) + labs(fill="index") + theme_minimal() 


# Compare stringency with Google mobility
comparison <- merge(stringency[stringency$category=="stringency_index",c("date", "index")],
                  mobility_fit[,c("date", "scaling")],
                  by="date", all=T)
comparison$scaling <- ifelse(is.na(comparison$scaling), 1, comparison$scaling)

ggplot(comparison) +
  geom_line(aes(x=date, y=index), color="red") +
  geom_line(aes(x=date, y=scaling*100), color="black")

ggplot(comparison) + 
  geom_point(aes(x=index, y=scaling)) +
  xlab("Oxford stringency index") + ylab("Google mobility (work places)")


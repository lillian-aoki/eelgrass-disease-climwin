# climwin with transect level data
# first try with epifauna data from Emmett and Carmen

# note, there are lots of gaps in epifauna (no AK in 2020, only 1/2 OR in 2020, no BB in 2021 yet)

library(tidyverse)
library(lme4)
library(climwin)
library(lubridate)
library(glmmTMB)
library(ggeffects)
library(sjPlot)
library(readxl)
region_order=c("AK", "BC", "WA", "OR", "BB", "SD")

# read in and summarize epifauna data
epi <- readxl::read_excel("data/EGWD_epifauna_data_v20220309.xlsx", sheet = "data")
epi_summ <- epi %>%
  group_by(SampleId=transect_unique_code) %>%
  summarise(EpiAbun=mean(as.numeric(epifauna_abundance_total, na.rm=TRUE))) %>%
  separate(SampleId, into = c("Region","SiteCode","Transect","Year"), sep = "[.]")
epi_summ$TidalHeight <- str_extract(epi_summ$Transect, "[U,L]")
epi_summ$Transect <- as.numeric(str_extract(epi_summ$Transect, "[1-6]"))
epi_summ$Year <- as.numeric(epi_summ$Year)

# read in environmental data at transect level
env19 <- read_csv("data/all_survey_metrics_transect_2019.csv")
env20 <- read_csv("data/all_survey_metrics_transect_2020.csv")
env21 <- read_csv("data/all_survey_metrics_transect_2021.csv")

env19 <- select(env19,-c("PrevalenceMean","SeverityMean","LesionAreaMean"))
env20 <- select(env20,-c("PrevalenceMean","SeverityMean","LesionAreaMean", "EpiphytePerAreaSe"))
env21 <- select(env21,-c("PrevalenceMean","SeverityMean","LesionAreaMean", "EpiphytePerAreaSe"))
env_c <- rbind(env19,env20, env21)

# read in disease data at blade level
disease19 <- read_csv("data/disease_2019.csv")
disease20 <- read_csv("data/disease_metrics_2020.csv")
disease21 <- read_csv("data/disease_metrics_2021.csv")

disease19 <- select(disease19, c("SampleDate","SampleId","Region","SiteCode","TidalHeight", "Transect", "Blade",
                                 "Year", "Prevalence", "LesionArea", "BladeArea", "Severity"))
disease20 <- select(disease20, c("SampleDate","SampleId","Region","SiteCode","TidalHeight", "Transect", "Blade",
                                 "Year", "Prevalence", "LesionArea", "BladeArea", "Severity"))
disease21 <- select(disease21, c("SampleDate","SampleId","Region","SiteCode","TidalHeight", "Transect", "Blade",
                                 "Year", "Prevalence", "LesionArea", "BladeArea", "Severity"))

dis_c <- rbind(disease19,disease20, disease21)

dis_summ <- dis_c %>%
  group_by(Year,Region, SiteCode,Transect, SampleDate) %>%
  summarise(TransectPrevalence=mean(Prevalence), TransectLesionArea=mean(LesionArea), TransectSeverity=mean(Severity),
            TransectBladeArea=mean(BladeArea))

dis_env <- full_join(dis_summ, env_c, by=c("Year","Region","SiteCode","Transect"))
# add epifauna
dis_env <- full_join(dis_env, epi_summ, by=c("Year","Region","SiteCode","TidalHeight","Transect"))
# limit to observations with no missing data
dis_env <- select(dis_env, c("SampleDate","Region","SiteCode","TidalHeight","Transect","TransectPrevalence",
                             "TransectLesionArea","TransectBladeArea","TransectSeverity","Year",
                             "DensityShootsMean", "CountBlades","EpiAbun"))
dis_env <- na.omit(dis_env)
# add IDs
dis_env$Meadow <- paste(dis_env$Region, dis_env$SiteCode, sep = "_")

# read in MUR data 
# use daily data?
mur <- read_csv("data/MUR_daily.csv")
meadows <- tibble(Meadow=unique(mur$Meadow))
# convert dates for climwin
mur$Day <- factor(format(mur$Date, "%d/%m/%Y"))

# limit disease data to MUR sites
dis_env <- left_join(meadows,dis_env, by="Meadow")
# scale continuous variables based on restricted data set
dis_env$sBladeArea <- scale(dis_env$TransectBladeArea, center=TRUE,scale=TRUE)
dis_env$sDensityShootsMean <- scale(dis_env$DensityShootsMean, center=TRUE, scale=TRUE)
dis_env$sEpiAbun <- scale(dis_env$EpiAbun, center=TRUE, scale=TRUE)
# convert dates for climwin
dis_env$SampleDay <- factor(format(dis_env$SampleDate, "%d/%m/%Y"))
# convert year to factor
dis_env$fYear <- as.factor(dis_env$Year)

# epifauna as predictor of disease ####
## prevalence ####
fit_prev1 <- glmer(TransectPrevalence ~ sBladeArea + sDensityShootsMean + sEpiAbun + TidalHeight + fYear + 
                     (1|Region) + (1|Meadow),
                   data = dis_env,
                   family = "binomial",
                   weights=CountBlades)
summary(fit_prev1)
names(mur)
# median sample date each year is ~July 20 (July 18 in 2019, July 22 in 2020 and 2021)
clim_mur_prev_tr <- slidingwin(xvar = list(TempAnomaly = mur$DiffMean, 
                                           TempAnomalyHeat = mur$DiffMeanHeat, 
                                           Temp = mur$analysed_sst),
                               cdate = mur$Day,
                               bdate = dis_env$SampleDay,
                               baseline = fit_prev1,
                               cohort = as.factor(dis_env$Year),
                               cinterval = "month",
                               range = c(6, 1),
                               type = c("absolute","relative"),
                               refday = c(20, 7),
                               cmissing = "method1",
                               stat = c("mean"),
                               func = "lin",
                               spatial = list(dis_env$Meadow, mur$Meadow))

clim_mur_prev_tr$combos
# best model has TempAnomaly in the month of sampling (immediately prior?)
plotdelta(dataset = clim_mur_prev_tr[[5]]$Dataset)
summary(clim_mur_prev_tr[[5]]$BestModel)
plot_model(clim_mur_prev_tr[[5]]$BestModel, show.values = TRUE, show.p = TRUE)
# epifauna has a positive relationship with disease prevalence at transect level

## severity ####

fit_sev1 <- lmer(TransectSeverity ~ sBladeArea + sDensityShootsMean + sEpiAbun + TidalHeight + fYear + 
                   (1|Region) + (1|Meadow),
                 data = dis_env)
summary(fit_sev1)
clim_mur_sev_tr <- slidingwin(xvar = list(TempAnomaly = mur$DiffMean,
                                           TempAnomalyHeat = mur$DiffMeanHeat,
                                           Temp = mur$analysed_sst),
                         cdate = mur$Day,
                         bdate = dis_env$SampleDay,
                         baseline = fit_sev1,
                         cohort = as.factor(dis_env$Year),
                         cinterval = "month",
                         range = c(6, 0),
                         type = c("absolute","relative"),
                         refday = c(20, 7),
                         cmissing = "method1",
                         stat = c("mean"),
                         func = "lin",
                         spatial = list(dis_env$Meadow, mur$Meadow))

clim_mur_sev_tr$combos
# adding epifauan improves severity model but delta AIC is smaller than for prevalence (<9) 
# and similar for all anomaly + anomalyheat metrics
# best AIC is with temp anomaly over April and May
plotdelta(dataset = clim_mur_sev_tr[[1]]$Dataset)
summary(clim_mur_sev_tr[[1]]$BestModel)
plot_model(clim_mur_sev_tr[[1]]$BestModel,show.values = TRUE, show.p = TRUE)
performance::model_performance(clim_mur_sev_tr[[1]]$BestModel)

# epifauna sensitivity to temp ####
# for purposes of climate window exploration, use a linear model
# for full model fitting, will likely use a gamma model (because epifauna abundance is always greater than 0)
fit_epi1 <- lmer(EpiAbun ~ sBladeArea + sDensityShootsMean + TidalHeight + fYear + 
                     (1|Region) + (1|Meadow),
                   data = dis_env)
summary(fit_epi1)
clim_mur_epi_tr <- slidingwin(xvar = list(TempAnomaly = mur$DiffMean, 
                                           TempAnomalyHeat = mur$DiffMeanHeat, 
                                           Temp = mur$analysed_sst),
                               cdate = mur$Day,
                               bdate = dis_env$SampleDay,
                               baseline = fit_epi1,
                               cohort = as.factor(dis_env$Year),
                               cinterval = "month",
                               range = c(6, 0),
                               type = c("absolute","relative"),
                               refday = c(20, 7),
                               cmissing = "method1",
                               stat = c("mean"),
                               func = "lin",
                               spatial = list(dis_env$Meadow, mur$Meadow))

clim_mur_epi_tr$combos
# best model has the epifauna abundance sensitive to temperature anomalies 3-2 months prior to sampling (or April)
plotdelta(dataset = clim_mur_epi_tr[[5]]$Dataset)
summary(clim_mur_epi_tr[[5]]$BestModel)
plot_model(clim_mur_epi_tr[[5]]$BestModel, show.values = TRUE, show.p = TRUE)

# so epifauna are sensitive to climate, more than many other factors
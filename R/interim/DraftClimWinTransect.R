# climwin with transect level data

library(tidyverse)
library(lme4)
library(climwin)
library(lubridate)
library(sjPlot)
region_order=c("AK", "BC", "WA", "OR", "BB", "SD")

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
# limit to observations with no missing data
dis_env <- select(dis_env, c("SampleDate","Region","SiteCode","TidalHeight","Transect","TransectPrevalence",
                              "TransectLesionArea","TransectBladeArea","TransectSeverity","Year",
                             "DensityShootsMean", "CountBlades"))
dis_env <- na.omit(dis_env)
# add IDs
dis_env$Meadow <- paste(dis_env$Region, dis_env$SiteCode, sep = "_")

# prepare 

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
# convert dates for climwin
dis_env$SampleDay <- factor(format(dis_env$SampleDate, "%d/%m/%Y"))
# convert year to factor
dis_env$fYear <- as.factor(dis_env$Year)

# fit baseline model
fit_prev1 <- glmer(TransectPrevalence ~ sBladeArea + sDensityShootsMean + TidalHeight + fYear + 
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
                            range = c(6, 0),
                            type = c("absolute","relative"),
                            refday = c(20, 7),
                            cmissing = "method1",
                            stat = c("mean"),
                            func = "lin",
                            spatial = list(dis_env$Meadow, mur$Meadow))

clim_mur_prev_tr$combos
# best model has TempAnomaly in the month immediately prior to sampling

plotdelta(dataset = clim_mur_prev_tr[[1]]$Dataset)
summary(clim_mur_prev_tr[[4]]$BestModel)
saveRDS(clim_mur_prev_tr, file = "output/climwin_prev_daily_MUR_transect.rds")
plot_model(clim_mur_prev_tr[[4]]$BestModel)
# now look at monthly mur data for transects
mur_monthly <- read_csv("data/MUR_monthly.csv")
clim_mur_prev_tr_month <- slidingwin(xvar = list(MeanTempAnomaly = mur_monthly$mDiffMean, 
                                           CTempAnomaly = mur_monthly$cDiffMeanHeat, 
                                           CPTempAnomaly = mur_monthly$cDiffMeanHeat),
                               cdate = mur_monthly$month,
                               bdate = dis_env$SampleDay,
                               baseline = fit_prev1,
                               cohort = as.factor(dis_env$Year),
                               cinterval = "month",
                               range = c(6, 0),
                               type = c("absolute","relative"),
                               refday = c(20, 7),
                               cmissing = "method1",
                               stat = c("mean"),
                               func = "lin",
                               spatial = list(dis_env$Meadow, mur_monthly$Meadow))
clim_mur_prev_tr_month$combos
# best model has MeanTempAnomaly in the month immediately prior to sampling
summary(clim_mur_prev_tr_month[[4]]$BestModel)
saveRDS(clim_mur_prev_tr_month, file = "output/climwin_prev_daily_MUR_transect_month.rds")
# plot prevalence vs shoot density
ggplot(dis_env, aes(x=DensityShootsMean, y=TransectPrevalence, color=Region))+geom_point()+
  facet_wrap(~Year)+
  scale_x_log10()


## now look for climate window for shoot density

fit_den <- lmer(DensityShootsMean ~ sBladeArea + TidalHeight + fYear + 
                  (1|Region) + (1|Meadow),
                data = dis_env)
E <- resid(fit_den)
F.den <- fitted(fit_den)
plot(E~F.den)

plot(resid(fit_den)~dis_env$DensityShootsMean)
clim_mur_den_tr <- slidingwin(xvar = list(TempAnomaly = mur$DiffMean, 
                                           TempAnomalyHeat = mur$DiffMeanHeat, 
                                           Temp = mur$analysed_sst),
                               cdate = mur$Day,
                               bdate = dis_env$SampleDay,
                               baseline = fit_den,
                               cohort = as.factor(dis_env$Year),
                               cinterval = "month",
                               range = c(6, 0),
                               type = c("absolute","relative"),
                               refday = c(20, 7),
                               cmissing = "method1",
                               stat = c("mean"),
                               func = "lin",
                               spatial = list(dis_env$Meadow, mur$Meadow))
clim_mur_den_tr$combos
# best models include terms for TempAnomalyHeat for month opening 4 months before July 20 and closing 3 months before and 
# Temp Anomaly Heat for 1 month window opening 5 months prior and closing 4 months prior to site-specific sampling date
summary(clim_mur_den_tr[[5]]$BestModel)
plotdelta(dataset = clim_mur_den_tr[[2]]$Dataset)
plot_model(clim_mur_den_tr[[5]]$BestModel)
# check with randomization
clim_mur_den_tr_randomized <- randwin(repeats=100,
                                      xvar = list(TempAnomaly = mur$DiffMean, 
                                                  TempAnomalyHeat = mur$DiffMeanHeat, 
                                                  Temp = mur$analysed_sst),
                                      cdate = mur$Day,
                                      bdate = dis_env$SampleDay,
                                      baseline = fit_den,
                                      cohort = as.factor(dis_env$Year),
                                      cinterval = "month",
                                      range = c(6, 0),
                                      type = c("absolute","relative"),
                                      refday = c(20, 7),
                                      cmissing = "method1",
                                      stat = c("mean"),
                                      func = "lin",
                                      spatial = list(dis_env$Meadow, mur$Meadow))
saveRDS(clim_mur_den_tr_randomized, file = "output/climwin_density_daily_MUR_transect.rds")
clim_mur_den_tr_randomized[[5]]
plothist(dataset=clim_mur_den_tr[[5]]$Dataset, datasetrand = clim_mur_den_tr_randomized[[5]])
pvalue(dataset=clim_mur_den_tr[[5]]$Dataset, datasetrand = clim_mur_den_tr_randomized[[5]], metric="AIC")
ggsave(filename = "output/hist_AIC_rand_clim_mur_den_transect.jpg")
names(dis_env)
ggplot(mur, aes(x=time, y=analysed_sst, color=Region))+geom_line()+theme_bw()


# add secondary analysis of VIIRS data ####

# read in VIIRS data 
# use daily data?
viirs <- read_csv("data/VIIRS_project_data_daily.csv")
meadows <- tibble(Meadow=unique(viirs$Meadow))
# convert dates for climwin
viirs$Day <- factor(format(viirs$Date, "%d/%m/%Y"))
# remove two sites? WA_E and SD_C have very limited data... enought for SD_C (173 points) maybe but not WA_E (81 points)
meadows <- subset(meadows,Meadow!="WA_E")
# limit disease data to MUR sites
dis_env <- left_join(meadows,dis_env, by="Meadow")
# scale continuous variables based on restricted data set
dis_env$sBladeArea <- scale(dis_env$TransectBladeArea, center=TRUE,scale=TRUE)
dis_env$sDensityShootsMean <- scale(dis_env$DensityShootsMean, center=TRUE, scale=TRUE)
# convert dates for climwin
dis_env$SampleDay <- factor(format(dis_env$SampleDate, "%d/%m/%Y"))
# convert year to factor
dis_env$fYear <- as.factor(dis_env$Year)

## prevalence ####
# fit baseline model
fit_prev1 <- glmer(TransectPrevalence ~ sBladeArea + sDensityShootsMean + TidalHeight + fYear + 
                     (1|Region) + (1|Meadow),
                   data = dis_env,
                   family = "binomial",
                   weights=CountBlades)
summary(fit_prev1)
names(viirs)
# median sample date each year is ~July 20 (July 18 in 2019, July 22 in 2020 and 2021)
clim_viirs_prev_tr <- slidingwin(xvar = list(TempAnomaly = viirs$DiffMean, 
                                           TempAnomalyHeat = viirs$DiffMeanHeat, 
                                           Temp = viirs$sst),
                               cdate = viirs$Date,
                               bdate = dis_env$SampleDay,
                               baseline = fit_prev1,
                               cohort = as.factor(dis_env$Year),
                               cinterval = "month",
                               range = c(6, 0),
                               type = c("absolute","relative"),
                               refday = c(20, 7),
                               cmissing = "method1",
                               stat = c("mean"),
                               func = "lin",
                               spatial = list(dis_env$Meadow, viirs$Meadow))

clim_viirs_prev_tr$combos
# best model has TempAnomalyHeat in 6 months prior to sampling (~Jan-Feb)

plotdelta(dataset = clim_viirs_prev_tr[[5]]$Dataset)
summary(clim_viirs_prev_tr[[5]]$BestModel)
saveRDS(clim_viirs_prev_tr, file = "output/climwin_prev_daily_viirs_transect.rds")
plot_model(clim_viirs_prev_tr[[5]]$BestModel)
# as with mur data, have a negative coefficient - prevalence decreasing with increasing temp anomaly - but now in Jan?

## severity ####
fit_sev1 <- lmer(TransectSeverity ~ sBladeArea + sDensityShootsMean + TidalHeight + fYear + (1|Meadow),
                 data = dis_env)
summary(fit_sev1)
clim_viirs_sev_tr <- slidingwin(xvar = list(TempAnomaly = viirs$DiffMean, 
                                             TempAnomalyHeat = viirs$DiffMeanHeat, 
                                             Temp = viirs$sst),
                                 cdate = viirs$Date,
                                 bdate = dis_env$SampleDay,
                                 baseline = fit_sev1,
                                 cohort = as.factor(dis_env$Year),
                                 cinterval = "month",
                                 range = c(5, 0),
                                 type = c("absolute","relative"),
                                 refday = c(01, 7),
                                 cmissing = "method1",
                                 stat = c("mean"),
                                 func = "lin",
                                 spatial = list(dis_env$Meadow, viirs$Meadow))

clim_viirs_sev_tr$combos
plotdelta(dataset = clim_viirs_sev_tr[[1]]$Dataset)
# best model is TempAnomaly 1 month before sampling/June
summary(clim_viirs_sev_tr[[1]]$BestModel)
plot_model(clim_viirs_sev_tr[[1]]$BestModel,show.p = TRUE, show.values = TRUE)
# weak negative effect size although significant

fit_les <- lmer(TransectLesionArea ~ sBladeArea + sDensityShootsMean + TidalHeight + fYear + (1|Meadow),
                data = dis_env)
summary(fit_les)
clim_viirs_les_tr <- slidingwin(xvar = list(TempAnomaly = viirs$DiffMean, 
                                            TempAnomalyHeat = viirs$DiffMeanHeat, 
                                            Temp = viirs$sst),
                                cdate = viirs$Date,
                                bdate = dis_env$SampleDay,
                                baseline = fit_les,
                                cohort = as.factor(dis_env$Year),
                                cinterval = "month",
                                range = c(5, 0),
                                type = c("absolute","relative"),
                                refday = c(20, 7),
                                cmissing = "method1",
                                stat = c("mean"),
                                func = "lin",
                                spatial = list(dis_env$Meadow, viirs$Meadow))

clim_viirs_les_tr$combos
# similarly improved models for TempAnomaly 2-3 months prior (April-May) and TempAnomalyHeat 5 months (Feb) prior
plotdelta(dataset = clim_viirs_les_tr[[3]]$Dataset)
# best model is TempAnomaly 1 month before sampling/June
summary(clim_viirs_les_tr[[4]]$BestModel)
plot_model(clim_viirs_les_tr[[3]]$BestModel,show.p = TRUE, show.values = TRUE)
# positive effect of climate on lesion area
dis_env$Month <- floor_date(dis_env)
ggplot(dis_env,aes(x=DiffMeanHeat,y=TransectSeverity, color=Meadow))+geom_point()+facet_wrap()

# logger data ####

# maybe logger data will come out in climwin?
hobo <- read_csv("data/NSFWD_HOBO_all_daily_site.csv")
hobo$Meadow <- paste(hobo$Region, hobo$SiteCode, sep="_")
names(hobo)
hobo <- subset(hobo, Year!=2019)
hobo <- hobo[-which(hobo$Year=="2021" & hobo$Region=="OR"),]
hobo <- hobo[-which(hobo$Year=="2020" & hobo$Region=="BB"),]
hobo$MeadowYear <- paste(hobo$Meadow, hobo$Year, sep="_")
hobo <- hobo[-which(hobo$MeadowYear== "AK_C_2021" | hobo$MeadowYear=="BC_C_2021" | hobo$MeadowYear=="WA_C_2021"|
                hobo$MeadowYear=="SD_A_2021" | hobo$MeadowYear=="SD_E_2021" | hobo$Meadow=="BB_A"),]
meadows <- tibble(Meadow=unique(hobo$Meadow))
# convert dates for climwin
hobo$Day <- factor(format(hobo$Day, "%d/%m/%Y"))
# # remove two sites? WA_E and SD_C have very limited data... enought for SD_C (173 points) maybe but not WA_E (81 points)
# meadows <- subset(meadows,Meadow!="WA_E")
# limit disease data to MUR sites
dis_env <- left_join(meadows,dis_env, by="Meadow")
# scale continuous variables based on restricted data set
dis_env$sBladeArea <- scale(dis_env$TransectBladeArea, center=TRUE,scale=TRUE)
dis_env$sDensityShootsMean <- scale(dis_env$DensityShootsMean, center=TRUE, scale=TRUE)
# convert dates for climwin
dis_env$SampleDay <- factor(format(dis_env$SampleDate, "%d/%m/%Y"))
# convert year to factor
dis_env$fYear <- as.factor(dis_env$Year)
names(dis_env)
dis_env$MeadowYear <- paste(dis_env$Meadow, dis_env$Year, sep="_")
unique(dis_env$MeadowYear)
meadowyears <- tibble(MeadowYear=unique(hobo$MeadowYear))
dis_env <- left_join(meadowyears, dis_env, by="MeadowYear")
dis_env <- na.omit(dis_env)
hobo$rangeTempC <- hobo$maxTempC-hobo$minTempC
# prevalence ####
# fit baseline model
fit_prev1 <- glmer(TransectPrevalence ~ sBladeArea + sDensityShootsMean + TidalHeight + fYear + 
                      (1|Meadow),
                   data = dis_env,
                   family = "binomial",
                   weights=CountBlades)
summary(fit_prev1)

# median sample date each year is ~July 20 (July 18 in 2019, July 22 in 2020 and 2021)
clim_hobo_prev_tr <- slidingwin(xvar = list(MeanTemp = hobo$meanTempC, 
                                             MaxTemp = hobo$maxTempC,
                                            MinTemp = hobo$minTempC,
                                            RangeTemp = hobo$rangeTempC),
                                 cdate = hobo$Day,
                                 bdate = dis_env$SampleDay,
                                 baseline = fit_prev1,
                                 cohort = as.factor(dis_env$Year),
                                 cinterval = "month",
                                 range = c(6, 0),
                                 type = c("absolute","relative"),
                                 refday = c(20, 7),
                                 cmissing = "method1",
                                 stat = c("mean"),
                                 func = "lin",
                                 spatial = list(dis_env$Meadow, hobo$Meadow))

clim_hobo_prev_tr$combos
# so no improvement by putting raw temperatures in
# severity ####
fit_sev1 <- lmer(TransectSeverity ~ sBladeArea + sDensityShootsMean  + fYear + (1|Meadow),
                 data = dis_env)
summary(fit_sev1)
clim_viirs_sev_tr <- slidingwin(xvar = list(MeanTemp = hobo$meanTempC, 
                                            MaxTemp = hobo$maxTempC,
                                            MinTemp = hobo$minTempC,
                                            RangeTemp = hobo$rangeTempC),
                                cdate = hobo$Day,
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
                                spatial = list(dis_env$Meadow, hobo$Meadow))

clim_viirs_sev_tr$combos
# absolutely no improvement
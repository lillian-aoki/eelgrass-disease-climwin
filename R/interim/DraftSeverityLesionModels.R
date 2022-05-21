# transect level models

library(glmmTMB)
library(DHARMa)
library(sjPlot)
library(tidyverse)
region_order=c("AK", "BC", "WA", "OR", "BB", "SD")
Region_order1 <- c("Alaska","British Columbia","Washington","Oregon","California - Bodega Bay","California - San Diego")


# import data for transect level analyses ####
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
# order regions
dis_env$Region <- ordered(dis_env$Region,levels=region_order)
# prepare temp data

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
####


# severity ####
# climate metric is temp anomaly over 1 month window opening 2 months prior to sampling date 
# same metric as for lesion area model below
# prep climate windows
samp <- distinct(select(dis_env, c("fYear","Meadow","SampleDate")))
write.csv(samp,"output/sample_date.csv",row.names = F)
samp$OpenWindow <- samp$SampleDate-days(60)
samp$CloseWindow <- samp$OpenWindow+days(30)
samp$MeadowYear <- paste(samp$Meadow, samp$fYear, sep="_")
meadow_years <- unique(samp$MeadowYear)
mur$fYear <- as.factor(format(floor_date(mur$Date, unit="years"), "%Y"))
mur$MeadowYear <- paste(mur$Meadow,mur$fYear, sep="_")
meadows <- unique(samp$Meadow)
years <- unique(as.character(samp$fYear))
out_years <- list()
out_meadow <- list()
for(i in seq_along(meadow_years)){
  #for(i in 19){
  mur_meadow_year <- subset(mur,MeadowYear==meadow_years[i])
  window_year <- subset(samp,MeadowYear==meadow_years[i])
  mur_window <- subset(mur_meadow_year, Date>window_year$OpenWindow & Date<window_year$CloseWindow)
  out_years[[i]] <- mur_window
}

mur_may_window <- bind_rows(out_years)
mur_may_summ <- mur_may_window %>%
  group_by(fYear, Region, Site, Meadow) %>%
  summarise(cDiffMean=sum(DiffMean,na.rm = TRUE), mDiffMean=mean(DiffMean))
head(mur_may_summ)
# now join the calculated anomaly to the diseae data
dis_may <- full_join(dis_env, mur_may_summ, by=c("fYear","Region","SiteCode"="Site","Meadow"))
dis_may$Region <- ordered(dis_may$Region, levels=region_order)
ggplot(dis_may, aes(x=mDiffMean, y=TransectSeverity))+geom_point()

# build model

fit_sev1 <- glmmTMB(TransectSeverity ~ TransectBladeArea + DensityShootsMean + TidalHeight + fYear + mDiffMean +
                      (1|Region) + (1|Meadow),
                    #dispformula = ~Region,
                    ziformula = ~1,
                    family= ziGamma('log'),
                    data = dis_may)

Esevsim <- simulateResiduals(fit_sev1)
plot(Esevsim)
# wow these residuals look beautiful
plot(Esevsim$scaledResiduals~dis_may$TransectBladeArea)
plot(Esevsim$scaledResiduals~dis_may$DensityShootsMean)
plot(Esevsim$scaledResiduals~dis_may$mDiffMean)
plot(Esevsim$scaledResiduals~as.factor(dis_may$TidalHeight))
plot(Esevsim$scaledResiduals~as.factor(dis_may$fYear))
# overall very nice. keep this model

sev_predict <- ggpredict(fit_sev1)
plot(sev_predict, rawdata = T)
ggplot()+
  geom_point(data=dis_may,aes(x=mDiffMean,y=TransectSeverity, color=Region),alpha=0.6, size=4)+
  geom_ribbon(data=sev_predict$mDiffMean, aes(x=x, ymin=conf.low, ymax=conf.high), fill="grey")+
  geom_line(data=sev_predict$mDiffMean, aes(x=x, y=predicted))+
  xlab("Average temperature anomalies in ~May (ºC)")+
  ylab("Wasting disease severity\n (% leaf area damaged)")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),expand = expansion(mult = c(0, .05)), limits=c(0,0.525))+
  scale_color_viridis_d(labels=c("Alaska","British Columbia","Washington","Oregon","California\n San Diego"))+
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12))
ggsave("output/severity_anomaly_plot.jpg", width=8, height=6)
# same plot without raw data
ggplot()+
#  geom_point(data=dis_may,aes(x=mDiffMean,y=TransectSeverity, color=Region),alpha=0.6, size=4)+
  geom_ribbon(data=sev_predict$mDiffMean, aes(x=x, ymin=conf.low, ymax=conf.high), fill="grey")+
  geom_line(data=sev_predict$mDiffMean, aes(x=x, y=predicted))+
  xlab("Average temperature anomalies in ~May (ºC)")+
  ylab("Wasting disease severity\n (% leaf area damaged)")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),expand = expansion(mult = c(0, .05)), limits=c(0,0.1))+
  theme_bw(base_size = 14)
ggsave("output/severity_anomaly_plot_model_only.jpg", width=6, height=4)
# lesion area ####

# climate metric is temp anomaly 2 months prior to sampling for 1 month

samp <- distinct(select(dis_env, c("fYear","Meadow","SampleDate")))
samp$OpenWindow <- samp$SampleDate-days(60)
samp$CloseWindow <- samp$OpenWindow+days(30)
samp$MeadowYear <- paste(samp$Meadow, samp$fYear, sep="_")
meadow_years <- unique(samp$MeadowYear)
mur$fYear <- as.factor(format(floor_date(mur$Date, unit="years"), "%Y"))
mur$MeadowYear <- paste(mur$Meadow,mur$fYear, sep="_")
meadows <- unique(samp$Meadow)
years <- unique(as.character(samp$fYear))
out_years <- list()
out_meadow <- list()
for(i in seq_along(meadow_years)){
  #for(i in 19){
  mur_meadow_year <- subset(mur,MeadowYear==meadow_years[i])
  window_year <- subset(samp,MeadowYear==meadow_years[i])
  mur_window <- subset(mur_meadow_year, Date>window_year$OpenWindow & Date<window_year$CloseWindow)
  out_years[[i]] <- mur_window
}

mur_may_window <- bind_rows(out_years)
mur_may_summ <- mur_may_window %>%
  group_by(fYear, Region, Site, Meadow) %>%
  summarise(cDiffMean=sum(DiffMean,na.rm = TRUE), mDiffMean=mean(DiffMean))
head(mur_may_summ)
# now join the calculated anomaly to the disease data
dis_may <- full_join(dis_env, mur_may_summ, by=c("fYear","Region","SiteCode"="Site","Meadow"))
ggplot(dis_may, aes(x=mDiffMean, y=TransectLesionArea))+geom_point()

# there's a question of which metric is easier to understand - cumulative anomaly or average anomaly? 
# Cumulative is more consistent with other work (e.g. for shoot densities) - but here it is net, not just positive

# first lesion area model with gaussian link
fit_les1 <- lme(TransectLesionArea ~ TransectBladeArea + DensityShootsMean + TidalHeight + fYear + mDiffMean,
                  random = ~ 1|Region/SiteCode,
                         data = dis_may)
E1 <- resid(fit_les1)
hist(E1)
qqnorm(E1)
qqline(E1)
plot(E1~dis_may$TransectBladeArea)
plot(E1~dis_may$DensityShootsMean)
plot(E1~dis_may$mDiffMean)
plot(E1~as.factor(dis_may$TidalHeight))
plot(E1~as.factor(dis_may$fYear))
# add a varaince structure, correlated with either mDiffMean or Density Shoots

fit_les2 <- lme(TransectLesionArea ~ TransectBladeArea + DensityShootsMean + TidalHeight + fYear + mDiffMean,
                random = ~ 1|Region/SiteCode,
                weights = varConstPower(form=~DensityShootsMean),
                data = dis_may)

E2 <- resid(fit_les2)
hist(E2)
qqnorm(E2)
qqline(E2)
plot(E2~dis_may$TransectBladeArea)
plot(E2~dis_may$DensityShootsMean)
plot(E2~dis_may$mDiffMean)

fit_les3 <- lme(TransectLesionArea ~ TransectBladeArea + DensityShootsMean + TidalHeight + fYear + mDiffMean,
                random = ~ 1|Region/SiteCode,
                weights = varConstPower(form=~TransectBladeArea),
                data = dis_may)

E2 <- resid(fit_les2)
hist(E2)
qqnorm(E2)
qqline(E2)
plot(E2~dis_may$TransectBladeArea)
plot(E2~dis_may$DensityShootsMean)
plot(E2~dis_may$mDiffMean)
fit_les4 <- lme(TransectLesionArea ~ TransectBladeArea + DensityShootsMean + TidalHeight + fYear + mDiffMean,
                random = ~ 1|Region/SiteCode,
                weights = varConstPower(form=~mDiffMean),
                data = dis_may)
E4 <- resid(fit_les4)
hist(E4)
qqnorm(E4)
qqline(E4)
plot(E4~dis_may$TransectBladeArea)
plot(E4~dis_may$DensityShootsMean)
plot(E4~dis_may$mDiffMean)
AIC(fit_les1, fit_les2, fit_les3, fit_les4)
# not seeing a huge change in residuals from adding the variance structure, but there is a big improvement by AIC
# Use fit_les4 to fit?
dis_may <- as.data.frame(dis_may)
full_les <- lme(TransectLesionArea ~ TransectBladeArea + DensityShootsMean + TidalHeight + fYear + mDiffMean,
                random = ~ 1|Region/SiteCode,
                weights = varConstPower(form=~mDiffMean),method = "ML",
                data = dis_may)
drop1(full_les)
summary(full_les)
les_pred <- ggeffect(full_les)
plot(les_pred, rawdata = T)

# alternatively, maybe this should be a gamma model? Because lesion area can't be <0?

fit_lesg <- glmmTMB(TransectLesionArea ~ TransectBladeArea + DensityShootsMean + TidalHeight + fYear + mDiffMean +
                      (1|Region) + (1|Meadow),
                    #dispformula = ~Region,
                    ziformula = ~1,
                    family= ziGamma('log'),
                    data = dis_may)

summary(fit_lesg)
Esim <- simulateResiduals(fit_lesg)
plot(Esim$scaledResiduals~dis_may$TransectBladeArea)
plot(Esim$scaledResiduals~dis_may$DensityShootsMean)
plot(Esim$scaledResiduals~dis_may$mDiffMean)
plot(Esim)
# okay this looks better in fact than the variance weight model so use this

lesg_predi <- ggpredict(fit_lesg)
plot(lesg_predi$mDiffMean, rawdata = T)
ggplot()+
  geom_point(data=dis_may,aes(x=mDiffMean,y=TransectLesionArea, color=Region),alpha=0.6, size=4)+
  geom_ribbon(data=lesg_predi$mDiffMean, aes(x=x, ymin=conf.low, ymax=conf.high), fill="grey")+
  geom_line(data=lesg_predi$mDiffMean, aes(x=x, y=predicted))+
  xlab("Average temperature anomalies in ~May (ºC)")+
  ylab(expression(paste("Wasting disease lesion area (cm"^"2"~")")))+
  #scale_y_continuous(labels = scales::percent_format(accuracy = 1),expand = expansion(mult = c(0, .05)), limits=c(0,0.525))+
  scale_color_viridis_d(labels=c("Alaska","British Columbia","Washington","Oregon","California\n San Diego"))+
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12))
ggsave("output/lesion_area_anomaly_plot.jpg", width=8, height=6)
# same plot without raw data looks a lot more significant?
ggplot()+
  #geom_point(data=dis_may,aes(x=mDiffMean,y=TransectLesionArea, color=Region),alpha=0.6, size=4)+
  geom_ribbon(data=lesg_predi$mDiffMean, aes(x=x, ymin=conf.low, ymax=conf.high), fill="grey")+
  geom_line(data=lesg_predi$mDiffMean, aes(x=x, y=predicted))+
  xlab("Average temperature anomalies in ~May (ºC)")+
  ylab(expression(paste("Wasting disease lesion area (cm"^"2"~")")))+
  #scale_y_continuous(labels = scales::percent_format(accuracy = 1),expand = expansion(mult = c(0, .05)), limits=c(0,0.525))+
  theme_bw(base_size = 14)
ggsave("output/lesion_area_anomaly_plot_model_only.jpg", width=6, height=4)
# 
# may <- togo[togo$Month.w=="May",]
# may <- select(may,c("Meadow","Year","Month","CPTA", "SST","TempAnomaly") )
# may_dis <- full_join(dis_env_c, may, by=c("Meadow","Year"))
# may_dis <- na.omit(may_dis)
# names(may_dis)
# may_dis$sSST <- scale(may_dis$SST)
# may_dis$sTempAnomaly <- scale(may_dis$TempAnomaly)
# 
# fit_prev <- glmmTMB(Prevalence ~ sBladeArea + sDensityShootsMean + TidalHeight + fYear + sTempAnomaly +
#                       (1|Region) + (1|Meadow) + (1|TransectID),
#                     #dispformula = ~Region,
#                     ziformula = ~1,
#                     family= ziGamma('log'),
#                     data = may_dis)
# summary(fit_prev)
# summary(may_dis$LesionArea)
# fit_les2 <- glmmTMB(LesionArea ~ sBladeArea + sDensityShootsMean + TidalHeight + fYear + sSST +
#                       (1|Region) + (1|Meadow) + (1|TransectID),
#                     #dispformula = ~Region,
#                     ziformula = ~1,
#                     family= ziGamma('log'),
#                     data = may_dis)
# summary(fit_les2)
# drop1(fit_les2)
# E.sim <- simulateResiduals(fit_les1)
# plot(E.sim)
# plotResiduals(E.sim$scaledResiduals, may_dis$sBladeArea)
# plotResiduals(E.sim$scaledResiduals, may_dis$sDensityShootsMean)
# plotResiduals(E.sim$scaledResiduals, may_dis$TidalHeight)
# plotResiduals(E.sim$scaledResiduals, may_dis$fYear)
# plotResiduals(E.sim$scaledResiduals, may_dis$sSST)
# 
# plot(E.sim$scaledResiduals~may_dis$sDensityShootsMean)
# plot(E.sim$scaledResiduals~sevTdat$sEpiphytePerAreaMean)
# plot(E.sim$scaledResiduals~as.factor(may_dis$TidalHeight))
# plot(E.sim$scaledResiduals~sevTdat$sCDiffMeanHeat)
# plot(BladeArea~LesionArea, may_dis)
# plot_model(fit_les2,show.values = TRUE, show.p = TRUE) + scale_y_log10(limits=c(0.2,5))
# ggsave("lesion_area_MaySST_effect_sizes.jpg", width=4, height=6)
# 
# may_dis_1 <- may_dis[may_dis$Prevalence==1,]
# 
# # fit_les3 <- glmmTMB(LesionArea ~ sBladeArea + sDensityShootsMean + TidalHeight + fYear + SST +
# #                       (1|Region) + (1|Meadow) + (1|TransectID),
# #                     #dispformula = ~Region,
# #                     family = Gamma('log'),
# #                     data = may_dis_1)
# # summary(fit_les3)
# # E3.sim <- simulateResiduals(fit_les3)
# # plot(E3.sim)
# # plot_model(fit_les2,show.values = TRUE, show.p = TRUE)
# 
# 
# # Severity ####
# 
# fit_sev1 <- glmmTMB(Severity ~ sBladeArea + sDensityShootsMean + TidalHeight + fYear + SST +
#                       (1|Region) + (1|Meadow) + (1|TransectID),
#                      #dispformula = ~Region,
#                      ziformula = ~1,
#                      family= ziGamma('log'),
#                      data = may_dis)
# summary(fit_sev1)
# Es1 <- simulateResiduals(fit_sev1)
# plot(Es1)
# 
# plotResiduals(Es1$scaledResiduals, may_dis$sBladeArea)
# plotResiduals(Es1$scaledResiduals, may_dis$sDensityShootsMean)
# plotResiduals(Es1$scaledResiduals, may_dis$TidalHeight)
# plotResiduals(Es1$scaledResiduals, may_dis$fYear)
# plotResiduals(Es1$scaledResiduals, may_dis$sSST)
# plot_model(fit_sev1, show.values = TRUE, show.p = TRUE) + scale_y_log10(limits=c(0.2,5))
# ggsave("severity_MaySST_effect_sizes.jpg", width=4, height=6)

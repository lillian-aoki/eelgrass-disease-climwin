## Refine Density model

## from climwin, use the climate window of the positive temperature anomaly starting five months prior to sampling date
# and ending 4 months prior to sampling

library(glmmTMB)
library(tidyverse)
library(lme4)
library(nlme)
library(climwin)
library(lubridate)
library(sjPlot)
library(RcppRoll)
library(effects)
library(ggeffects)
library(performance)
region_order=c("AK", "BC", "WA", "OR", "BB", "SD")

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
####

# identify the sample dates ####
samp <- distinct(select(dis_env, c("fYear","Meadow","SampleDate")))
samp$OpenWindow <- samp$SampleDate-days(152)
samp$CloseWindow <- samp$OpenWindow+days(60) # for density model including prevalence, 2 month window
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

mur_febmar_window <- bind_rows(out_years)
mur_febmar_summ <- mur_febmar_window %>%
  group_by(fYear, Region, Site, Meadow) %>%
  summarise(cDiffMeanHeat=sum(DiffMeanHeat,na.rm = TRUE), mDiffMeanHeat=mean(DiffMeanHeat))
head(mur_febmar_summ)
# now join the calculated anomaly to the diseae data
dis_febmar <- full_join(dis_env, mur_febmar_summ, by=c("fYear","Region","SiteCode"="Site","Meadow"))
ggplot(dis_febmar, aes(x=cDiffMeanHeat, y=DensityShootsMean))+geom_point()
which(is.na(dis_febmar))
# fit initial density model, using gaussian distribution

fit_den0 <- lmer(DensityShootsMean ~ TransectBladeArea + TidalHeight + fYear + cDiffMeanHeat +
                   (1|Region) + (1|Meadow),
                 data = dis_febmar)
E0 <- resid(fit_den0)
F0 <- fitted(fit_den0)
hist(E0)
# a little skewed
qqnorm(E0)
qqline(E0)
# residuals are skewed at larger values
plot(E0~dis_febmar$TransectBladeArea)
# pattern in residuals at low blade area 
plot(E0~as.factor(dis_febmar$TidalHeight))
plot(E0~dis_env$fYear)
summary(fit_den0)
# fit using nlme::lme() because can include varying variance structure
fit_den1 <- lme(DensityShootsMean ~ TransectBladeArea + TidalHeight + fYear + cDiffMeanHeat,
                random = ~ 1|Region/SiteCode,
              data = dis_febmar)
summary(fit_den1) # same AIC etc as lmer

# Add variance structure, correlated with Blade Area
fit_den2 <- lme(DensityShootsMean ~ TransectBladeArea + TidalHeight + fYear + cDiffMeanHeat,
                random = ~ 1|Region/SiteCode,
                weights = varPower(form=~TransectBladeArea),
                data = dis_febmar)
E2 <- resid(fit_den2, type = "normalized")
hist(E2)
# resids are more normal
qqnorm(E2)
qqline(E2)
# not perfect though
F2 <- fitted(fit_den2)
plot(E2~F2)
# cone shape is mostly gone
plot(E2~dis_febmar$TransectBladeArea)
plot(E2~as.factor(dis_febmar$Region))
plot(E2~dis_febmar$cDiffMeanHeat)
# try variance structure that varies by reigon
fit_den3 <- lme(DensityShootsMean ~ TransectBladeArea + TidalHeight + fYear + cDiffMeanHeat,
                random = ~ 1|Region/SiteCode,
                weights = varIdent(~Region),
                data = dis_febmar)
E3 <- resid(fit_den3, type = "normalized")
hist(E3)
# still skewed
qqnorm(E3)
qqline(E3)
F3 <- fitted(fit_den3)
plot(E3~F3)
plot(E3~dis_febmar$TransectBladeArea)
plot(E3~as.factor(dis_febmar$Region))
# Add variance structure, correlated with Blade Area
fit_den4 <- lme(DensityShootsMean ~ TransectBladeArea + TidalHeight + fYear+ cDiffMeanHeat,
                random = ~ 1|Region/SiteCode,
                weights = varExp(form=~TransectBladeArea),
                data = dis_febmar)
E4 <- resid(fit_den4, type = "normalized")
hist(E4)
# resids are more normal
qqnorm(E4)
qqline(E4)
# not perfect though
F4 <- fitted(fit_den4)
plot(E4~F4)
# more extreme positive residual values 
plot(E4~dis_febmar$TransectBladeArea)
plot(E4~as.factor(dis_febmar$Region))
AIC(fit_den4)
# worse AIC then with varPower
# try combined variance structure?
# Add variance structure, correlated with Blade Area
fit_den5 <- lme(DensityShootsMean ~ TransectBladeArea + TidalHeight + fYear + cDiffMeanHeat,
                random = ~ 1|Region/SiteCode,
                weights = varComb(varIdent (form=~Region), varPower(form=~TransectBladeArea)),
                data = dis_febmar)
E5 <- resid(fit_den5, type = "normalized")
hist(E5)
# resids are more normal
qqnorm(E5)
qqline(E5)
# not perfect though
F5 <- fitted(fit_den5)
plot(E5~F5)
# cone shape is mostly gone
plot(E5~dis_febmar$TransectBladeArea)
plot(E5~as.factor(dis_febmar$Region))

# Add variance structure, correlated with Blade Area
fit_den6 <- lme(DensityShootsMean ~ TransectBladeArea + TidalHeight + fYear+ TransectPrevalence + cDiffMeanHeat,
                random = ~ 1|Region/SiteCode,
                weights = varConstPower(form=~TransectBladeArea),
                data = dis_febmar)
E6 <- resid(fit_den6, type = "normalized")
hist(E6)
# resids are more normal
qqnorm(E6)
qqline(E6)
# not perfect though
F6 <- fitted(fit_den6)
plot(E6~F6)
# cone shape is mostly gone
plot(E6~dis_febmar$TransectBladeArea)
plot(E6~as.factor(dis_febmar$Region))
plot(E6~dis_febmar$cDiffMeanHeat)
plot(E6~dis_febmar$TransectPrevalence)
AIC(fit_den2, fit_den6)
# const power is the best variance structure by AIC
# residuals are reasonably normal
rm(full_model)
dat <- select(dis_febmar,c("DensityShootsMean","TransectBladeArea","TidalHeight","fYear","cDiffMeanHeat","Region","SiteCode"))
dat$TidalHeight <- as.factor(dat$TidalHeight)
dat$Region <- as.factor(dat$Region)
dat$SiteCode <- as.factor(dat$SiteCode)
summary(dat)
dis_febmar <- as.data.frame(dis_febmar)
class(dat)
full_model <- lme(DensityShootsMean ~ TransectBladeArea + TidalHeight + fYear + TransectPrevalence + mDiffMeanHeat,
                  random = ~ 1|Region/SiteCode,
                  weights = varConstPower(form=~TransectBladeArea),
                  data = dis_febmar,
                  method = "ML"
                  )
drop1(full_model)
summary(full_model)
model_performance(full_model)
plot_model(full_model, type = "std")
fitgg_lme <- ggeffect(model=full_model)
ggplot()+
  geom_point(data=dis_febmar,aes(x=mDiffMeanHeat,y=DensityShootsMean, color=fYear),alpha=0.6, size=4)+
  geom_ribbon(data=fitgg_lme$mDiffMeanHeat, aes(x=x, ymin=conf.low, ymax=conf.high), fill="grey")+
  geom_line(data=fitgg_lme$mDiffMeanHeat, aes(x=x, y=predicted))+
  scale_y_log10()+
  xlab("Average positive temperature anomalies \n in ~Feb-Mar (ºC)")+
  ylab(expression(paste("Seagrass shoot density (shoots m"^"-2"~")")))+
  scale_color_discrete(name="Year")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())
ggsave("output/density_anomaly_mean_plot.jpg",width = 6, height = 4)
ggplot()+
  #geom_point(data=dis_febmar,aes(x=mDiffMeanHeat,y=DensityShootsMean, color=fYear),alpha=0.6, size=4)+
  geom_ribbon(data=fitgg_lme$mDiffMeanHeat, aes(x=x, ymin=conf.low, ymax=conf.high), fill="grey")+
  geom_line(data=fitgg_lme$mDiffMeanHeat, aes(x=x, y=predicted))+
  #scale_y_log10()+
  xlab("Average positive temperature anomalies \n in ~Feb-Mar (ºC)")+
  ylab(expression(paste("Seagrass shoot density (shoots m"^"-2"~")")))+
  scale_color_discrete(name="Year")+
  theme_bw(base_size = 14)
ggsave("output/density_anomaly_mean_plot_model_only.jpg",width = 6, height = 4)
plot_model(full_model)

ggplot()+
  geom_point(data=dis_febmar,aes(x=TransectPrevalence,y=DensityShootsMean, color=fYear),alpha=0.6, size=4)+
  geom_ribbon(data=fitgg_lme$TransectPrevalence, aes(x=x, ymin=conf.low, ymax=conf.high), fill="grey")+
  geom_line(data=fitgg_lme$TransectPrevalence, aes(x=x, y=predicted))+
  #scale_y_log10()+
  #xlab("Cumulative positive temperature anomalies \n in ~Feb-Mar (ºC)")+
  #ylab(expression(paste("Seagrass shoot density (shoots m"^"-2"~")")))+
  scale_color_discrete(name="Year")+
  theme_bw(base_size = 14)

plot(fitgg_lme$TransectPrevalence, rawdata = T)#+scale_y_log10()
fitgg <- ggeffect(model=fit_den0) # works for the lmer object
plot(fitgg, rawdata = T)
full_model_lmer


den_gam <- glmmTMB(DensityShootsMean ~ TransectBladeArea + TidalHeight + fYear + TransectPrevalence + cDiffMeanHeat +
                     (1|Region)+ (1|Meadow),
                   ziformula = ~1,
                   family= ziGamma('log'),
                   data = dis_febmar)
summary(den_gam)
Eden <- simulateResiduals(den_gam)
plot(Eden$scaledResiduals~dis_febmar$TransectBladeArea)
plot(Eden$scaledResiduals~dis_febmar$DensityShootsMean)
plot(Eden$scaledResiduals~dis_febmar$cDiffMeanHeat)
plot(Eden)
# okay these resids look very bad so would need to work on this a bunch


### density vs prevalence for sites
names(dis_febmar)
ggplot(dis_febmar, aes(x=DensityShootsMean, y=TransectPrevalence, color=Region))+geom_point()+facet_wrap(~Year)
summ <- dis_febmar %>%
  group_by(fYear, Year, Region, Meadow) %>%
  summarise(MeadowPrevalence=mean(TransectPrevalence), MeadowSeverity=mean(TransectSeverity),MeadowDensity=mean(DensityShootsMean))


ggplot(summ, aes(x=MeadowDensity, y=MeadowPrevalence, color=Region))+geom_point()


# # should this be gaussian because non-zero and also not strictly integers/counts
# # rm(newdata)
# # newdata <- with(dis_febmar, 
# #                 expand.grid(fYear=as.factor(seq(2019, 2021, 1)), 
# #                             #SiteCode=c("A","B"), 
# #                             #Region=c("WA","OR"), 
# #                             TidalHeight=c("U","L"), 
# #                             cDiffMeanHeat=seq(from = min(cDiffMeanHeat), to = max(cDiffMeanHeat), length.out = 100),
# #                             TransectBladeArea=seq(from = min(TransectBladeArea), to = max(TransectBladeArea), length.out = 100)))#,
# #                             #TransectBladeArea=median(TransectBladeArea)))
# 
# 
# fit1<-predict(full_model, type ="response", newdata,  se=TRUE)
# newdata<-cbind(newdata, fit1$fit)
# newdata<-cbind(newdata, fit$se.fit)
# trydata<- tibble(fit1)
# 
# fit <- cbind(
#   response=predict(full_model,newdata,type='response'),
#   variance=predict(full_model,newdata,type='variance'),
#   q_025=predict(full_model,newdata,type='quantile',at=c(0.025)),
#   q_975=predict(full_model,newdata,type='quantile',at=c(0.975)))
# mod_pred <- as.data.frame(fit)
# mod_pred <- cbind(mod_pred, newdata)
# ggplot(mod_pred,aes(x=cDiffMeanHeat, color=fYear))+
#   geom_line(aes(y=response))+
#   geom_line(aes(y=q_025),linetype="dashed")+
#   geom_line(aes(y=q_975),linetype="dashed")+
#   facet_wrap(~TidalHeight)
# ggplot(newdata, aes(x=cDiffMeanHeat, y=fit, color=fYear))+geom_line()
# p<-ggplot(newdata, aes(x=cDiffMeanHeat, y=TransectBladeArea, fill=fit))+
#   geom_tile()+
#   scale_fill_gradient2(low = "blue", mid="yellow", high = "red", midpoint=0.5)
# p+facet_grid(.~fYear)+theme_bw()+ 
#   theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank())+ylab(bquote('Leaf Area (' ~cm^2* ')'))+
#   xlab(expression(paste("Cumulative positive temperature anomaly (", degree,"C) for a 1 month window, beginning 5 months prior to sampling", )))+
#   theme(panel.spacing.x = unit(4, "mm"))+labs(fill='Shoot Density')
# 
# 
# 
# # First simulate change in CPTA
# jvaluesPa <- with(dis_febmar, seq(from = min(cDiffMeanHeat), to = max(cDiffMeanHeat), length.out = 100))
# # create new data and hold blade area at the median value of the dataset (29.1 cm2)
# b.dataPa <- data.frame(cDiffMeanHeat=jvaluesPa,TransectBladeArea=median(dis_febmar$TransectBladeArea),
#                        Region="WA",SiteCode="A", TidalHeight="U",fYear=as.factor("2020"))
# predPa <- cbind(
#   response=predict(full_model,newdata=b.dataPa,type='response'),
#   variance=predict(full_model,newdata=b.dataPa,type='variance'),
#   predict(full_model,newdata=b.dataPa,type='quantile',at=c(0.025,0.975)))
# preva <- as.data.frame(predPa)
# preva <- cbind(preva,b.dataPa)
# a <- ggplot(preva,aes(x=CPTempAnomaly))+
#   geom_line(aes(y=response))+
#   geom_line(aes(y=q_0.025),linetype="dashed")+
#   geom_line(aes(y=q_0.975),linetype="dashed")+
#   geom_point(data=dat,aes(x=CPTempAnomaly,y=PrevalenceMean,color=Region),size=2)+
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ 
#   scale_color_viridis_d()+
#   xlab("Cumulative positive temperature anomaly (ºC)")+
#   ylab("Wasting disease prevalence\n (% individuals infected)")+
#   theme_bw()+
#   theme(panel.grid = element_blank())
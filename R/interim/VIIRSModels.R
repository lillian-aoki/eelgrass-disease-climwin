# transect level models for VIIRS data 
# using climate metrics identified with climwin for MUR data

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

# climate metric is temp anomaly over 1 month window opening 2 months prior to sampling date 
# same metric as for lesion area model below
# prep climate windows
samp <- distinct(select(dis_env, c("fYear","Meadow","SampleDate")))
# write.csv(samp,"output/sample_date.csv",row.names = F)
samp$OpenWindow <- samp$SampleDate-days(60)
samp$CloseWindow <- samp$OpenWindow+days(30)
samp$MeadowYear <- paste(samp$Meadow, samp$fYear, sep="_")
meadow_years <- unique(samp$MeadowYear)
viirs$fYear <- as.factor(format(floor_date(viirs$Date, unit="years"), "%Y"))
viirs$MeadowYear <- paste(viirs$Meadow,viirs$fYear, sep="_")
meadows <- unique(samp$Meadow)
years <- unique(as.character(samp$fYear))
out_years <- list()
out_meadow <- list()
for(i in seq_along(meadow_years)){
  #for(i in 19){
  viirs_meadow_year <- subset(viirs,MeadowYear==meadow_years[i])
  window_year <- subset(samp,MeadowYear==meadow_years[i])
  viirs_window <- subset(viirs_meadow_year, Date>window_year$OpenWindow & Date<window_year$CloseWindow)
  out_years[[i]] <- viirs_window
}

viirs_may_window <- bind_rows(out_years)
viirs_may_summ <- viirs_may_window %>%
  group_by(fYear, Region, Site, Meadow) %>%
  summarise(cDiffMean=sum(DiffMean,na.rm = TRUE), mDiffMean=mean(DiffMean, na.rm=TRUE))
head(viirs_may_summ)
# now join the calculated anomaly to the diseae data
dis_may <- full_join(dis_env, viirs_may_summ, by=c("fYear","Region","SiteCode"="Site","Meadow"))
dis_may$Region <- ordered(dis_may$Region, levels=region_order)
dis_may <- na.omit(dis_may)
ggplot(dis_may, aes(x=cDiffMean, y=TransectSeverity, color=Region, shape=fYear))+geom_point(size=4)+
  xlab("Cumulative May anomaly - VIIRS")+
  ylab("Transect level severity")+
  theme_bw()+
  theme(panel.grid = element_blank())
ggsave("output/severity_VIIRS_anomaly.jpg")

ggplot(dis_may, aes(x=cDiffMean, y=TransectLesionArea, color=Region, shape=fYear))+geom_point(size=4)+
  xlab("Cumulative May anomaly - VIIRS")+
  ylab("Transect level lesion area")+
  theme_bw()+
  theme(panel.grid = element_blank())
ggsave("output/lesion_VIIRS_anomaly.jpg")
# negative relationship = bad news

# mur ####
#join with mur data and see how it looks
## import murs and arrange 
# read in MUR data 
# use daily data?
mur <- read_csv("data/MUR_daily.csv")
meadows <- tibble(Meadow=unique(mur$Meadow))
# convert dates for climwin
mur$Day <- factor(format(mur$Date, "%d/%m/%Y"))

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

mur_may_summ <- rename(mur_may_summ, "MUR_mDiffMean"="mDiffMean", "MUR_cDiffMean" = "cDiffMean")
viirs_may_summ <- rename(viirs_may_summ, "VIIRS_mDiffMean"="mDiffMean", "VIIRS_cDiffMean" = "cDiffMean")
may_combo <- full_join(mur_may_summ, viirs_may_summ, by=c("fYear", "Region", "Site", "Meadow"))
may_combo <- na.omit(may_combo)
ggplot(may_combo, aes(x=MUR_cDiffMean, y=VIIRS_cDiffMean, color=Meadow, shape=fYear))+geom_point(size=4)+
  geom_abline(slope=1, intercept = 0)+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_vline(xintercept = 0, linetype="dashed")+
  labs(title = "Compare anomalies from VIIRS and MUR datasets")+
  xlab("Cumulative May anomaly - MUR")+
  ylab("Cumulative May anomaly - VIIRS")+
  theme_bw()+
  theme(panel.grid = element_blank())
ggsave(filename = "output/VIIRS_MUR_May_anomaly.jpg")

# compare MUR and VIIRS long-term means
names(mur)
mur_Tmean_ma <- select(mur, c(Meadow, Region, Site, Julian, Tmean_ma))
mur_Tmean_ma <- distinct(mur_Tmean_ma)
names(viirs)
viirs_ma <- select(viirs, c(Meadow, Region, Site, Julian, Tmean_ma))
viirs_ma <- distinct(viirs_ma)
ma <- right_join(mur_Tmean_ma, viirs_ma, by=c("Meadow", "Region", "Site", "Julian"))
ma <- rename(ma, "MUR_Tmean_ma"="Tmean_ma.x", "VIIRS_Tmean_ma"="Tmean_ma.y")
## add ghr to plots
longterm <- read_csv("data/MUR_G1SST_9Y_Mean_90th.csv") # this is the long-term means for MUR and GHR
ghr_ma <- select(longterm, c(Meadow, Region, Site, Julian, Tmean_ma))
ghr_ma <- distinct(ghr_ma)
ghr_ma <- rename(ghr_ma, "G1SST_Tmean_ma"="Tmean_ma")
ma <- right_join(ghr_ma, ma, by=c("Meadow", "Region", "Site", "Julian"))

ggplot(ma[ma$Meadow!="SD_C",])+geom_line(aes(x=Julian, y=MUR_Tmean_ma, color="MUR"))+
  geom_line(aes(x=Julian, y=G1SST_Tmean_ma, color="G1SST"))+
  geom_line(aes(x=Julian, y=VIIRS_Tmean_ma, color="VIIRS"), show.legend = TRUE)+
  facet_wrap(~Meadow, scales="free")+
  labs(title = "Compare Long-term means from VIIRS and MUR/G1SST")+
  xlab("Day of year")+
  ylab("Sea surface temperature")+
  theme_bw()
ggsave("output/compare_viirs_mur_longterm_mean.jpg")
# not that different


# so is it the actual temps that are different?
mur_lim <- select(mur, c(Meadow, Region, Site, Julian, fYear, Date, analysed_sst))
viirs_lim <- select(viirs, c(Meadow, Region, Site, Julian, fYear, Date, sst))
viirs_meadows <- unique(viirs_lim$Meadow)
combo <- right_join(viirs_lim, mur_lim, by=c("Meadow", "Region", "Site", "Julian", "fYear", "Date"))
combo <- combo[combo$Meadow %in% viirs_meadows, ]
ggplot(combo)+geom_line(aes(x=Date, y=analysed_sst))+
  geom_line(aes(x=Date, y=sst), color="blue")+
  facet_wrap(~Meadow,scales = "free_y")+
  ylab("Temperature")+
  labs(title = "Compare VIIRS and MUR temperatures",
       subtitle = "Blue line = VIIRS, black line = MUR")+
  theme_bw()
ggsave(filename = "output/compare_VIIRS_MUR_temp_records.jpg")

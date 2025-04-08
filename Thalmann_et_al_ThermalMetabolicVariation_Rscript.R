### All Analyses for Thalmann et al. 2025: Pacific Cod metabolism and swimming performance are similar across temperatures following prolonged thermal acclimation
#Updated on 4/8/2025

#Load Libraries
library(rstudioapi) #set working directory
library(tidyverse) #plotting, data wrangling
library(ggpubr) #plots
library(lme4) # Run mixed effects linear models
library(car) #ANOVA function for models
library(corrr) #Correlations
library(respirometry) #scaled mass
library(vegan) #PCA
library(RColorBrewer) #Color palette
library(investr) #Confidence Intervals for COT
library(MuMIn) #Runs the R.squaredGLMM function
library(effects) #Marginal Means
library(ggeffects) #Visualize Marginal Means

#Set Working Directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

SwimData <- read.csv("Thalmann_et_al_PacificCodThermalMetabolicVariation_Data.csv") %>%
  mutate(Expt = as.factor(Expt)) %>%
  mutate(GSI = (GonadWWTg/WholeFishWWT_g)*100) %>%
  mutate(SalinityRatio = Salinity_DayofTrial/MeanSalinity_3Months) %>%
  mutate(SwimSpeed_mm = SwimSpeed * SL_mm) %>%
  mutate(Ucrit_mm = Ucrit * SL_mm) %>%
  mutate(MO2_J = MO2 * 3.24 * 4.18) %>% #convert MO2 to Joules
  mutate(COT = (((MO2_J / SwimSpeed_mm)/3600) * 1000000)) #### Calculate Cost of Transport (J*kg^-1*km^-1)

#### Mass Scale Performance ####
#Use the quarter-power scaling coefficient b = -0.25
#Experiment 1 Population Mean: 95.0
#Experiment 2 Population Mean: 88.7

ScaledMassExpt1 <- SwimData %>%
  filter(Expt == 1) %>%
  mutate(SMR_scaled = scale_MO2(WholeFishWWT_g, SMR, 95.0, b =-0.25)) %>%
  mutate(RMR_scaled = scale_MO2(WholeFishWWT_g, RMR, 95.0, b =-0.25)) %>%
  mutate(MMR_scaled = scale_MO2(WholeFishWWT_g, MMR, 95.0, b =-0.25)) %>%
  mutate(FAS_scaled = scale_MO2(WholeFishWWT_g, FAS, 95.0, b =-0.25)) %>%
  mutate(AAS_scaled = scale_MO2(WholeFishWWT_g, AAS, 95.0, b =-0.25)) %>%
  mutate(Ucrit_scaled = scale_MO2(WholeFishWWT_g, Ucrit, 95.0, b =-0.25))%>%
  mutate(Ucrit_mm_scaled = scale_MO2(WholeFishWWT_g, Ucrit_mm, 95.0, b =-0.25)) %>%
  mutate(MO2_scaled = scale_MO2(WholeFishWWT_g, MO2, 95.0, b = -0.25)) %>%
  mutate(Temp_Expt = case_when( Temp_Trmt ==  2 ~"2°C Expt 1", Temp_Trmt == 4 ~"4°C Expt 1",  Temp_Trmt == 6 ~ "6°C Expt 1", Temp_Trmt == 8 ~ "8°C Expt 1")) 


ScaledMassExpt2 <- SwimData %>%
  filter(Expt == 2) %>%
  mutate(SMR_scaled = scale_MO2(WholeFishWWT_g, SMR, 88.7, b =-0.25)) %>%
  mutate(RMR_scaled = scale_MO2(WholeFishWWT_g, RMR, 88.7, b =-0.25)) %>%
  mutate(MMR_scaled = scale_MO2(WholeFishWWT_g, MMR, 88.7, b =-0.25)) %>%
  mutate(FAS_scaled = scale_MO2(WholeFishWWT_g, FAS, 88.7, b =-0.25)) %>%
  mutate(AAS_scaled = scale_MO2(WholeFishWWT_g, AAS, 88.7, b =-0.25)) %>%
  mutate(Ucrit_scaled = scale_MO2(WholeFishWWT_g, Ucrit, 88.7, b =-0.25)) %>%
  mutate(Ucrit_mm_scaled = scale_MO2(WholeFishWWT_g, Ucrit_mm, 88.7, b =-0.25)) %>%
  mutate(MO2_scaled = scale_MO2(WholeFishWWT_g, MO2, 88.7, b = -0.25)) %>%
  mutate(Temp_Expt = case_when( Temp_Trmt ==  6 ~"6°C Expt 2", Temp_Trmt == 10 ~"10°C Expt 2",  Temp_Trmt == 14 ~ "14°C Expt 2")) 

SwimData_Scaled <- rbind(ScaledMassExpt1, ScaledMassExpt2) %>%
  mutate(TankID_Expt = factor(TankID_Expt, levels = c("2", "4", "6", "8",  "6A", "6B", "10A", "10B", "14A", "14B"))) %>%
  mutate(FishID = as.factor(FishID)) %>%
  group_by(FishID, SwimSpeed) %>% 
  mutate(TrialNo = row_number()) %>%
  mutate(MO2_scaled_J = MO2_scaled * 3.24 * 4.18) %>% #convert MO2 to Joules
  mutate(COT_scaled = (((MO2_scaled_J / SwimSpeed_mm)/3600) * 1000000)) %>%
  mutate(Expt = as.factor(Expt)) %>%
  mutate(Expt = factor(Expt, levels = c("1", "2"))) %>%
  mutate(Temp_Expt = factor(Temp_Expt, levels = c("2°C Expt 1", "4°C Expt 1", "6°C Expt 1", "6°C Expt 2", "8°C Expt 1",   "10°C Expt 2", "14°C Expt 2"))) %>%
  mutate(Temp= case_when( Temp_Trmt ==  2 ~"2°C", Temp_Trmt == 4 ~"4°C",  Temp_Trmt == 6 ~ "6°C", Temp_Trmt == 8 ~ "8°C", Temp_Trmt == 10 ~ "10°C", Temp_Trmt == 14 ~ "14°C")) %>%
  mutate(Temp = factor(Temp, levels = c("2°C", "4°C", "6°C", "8°C",  "10°C", "14°C"))) 


##### MeanMO2 Per Speed ####
#Calculate one MO2 value for each swimming speed for each individual
MeanMO2 <- SwimData_Scaled %>%
  select(FishID, Expt, Temp_Trmt, Temp, Temp_Expt, TankID_Expt, MeanTemp_3Months, MO2_scaled, COT_scaled, COT, SwimSpeed, SwimSpeed_mm, SMR_AdjRsq, TankID_Expt) %>%
  group_by(FishID, Expt, MeanTemp_3Months, Temp_Trmt, Temp, Temp_Expt, SwimSpeed_mm, SwimSpeed, SMR_AdjRsq, TankID_Expt) %>%
  reframe(MO2_scaled_mean = mean(MO2_scaled), COT_scaled_mean = mean(COT_scaled), COT_mean = mean(COT)) %>%
  mutate(Expt = case_when(Expt == 1 ~ "Expt 1", Expt == 2 ~ "Expt 2"))

#### Swim Data by Individual Fish #####
#### Summarise data from a full swimming trial into one row of data per fish
SwimData_ByFish <- SwimData_Scaled %>%
  select(FishID, Expt, Temp_Trmt, TankID_Expt, SMR_scaled, SMR_AdjRsq, RMR_scaled, MMR_scaled, Ucrit_mm_scaled, FAS_scaled, AAS_scaled, SMR, RMR, MMR, Ucrit_mm, FAS, AAS, WholeFishWWT_g, HSI, LiverLipids, MuscleLipids, Growth, GSI, ParasiteWWT, MeanTemp_3Months, SalinityRatio, DaysSinceFeeding, DaylightHours, COT) %>%
  group_by(FishID, Expt, Temp_Trmt, TankID_Expt) %>%
  reframe(SMR_scaled = mean(SMR_scaled), RMR_scaled = mean(RMR_scaled), MMR_scaled = mean(MMR_scaled), Ucrit_mm_scaled = mean(Ucrit_mm_scaled), FAS_scaled = mean(FAS_scaled), AAS_scaled = mean(AAS_scaled), SMR = mean(SMR), RMR = mean(RMR), MMR = mean(MMR), Ucrit_mm = mean(Ucrit_mm), FAS = mean(FAS), AAS = mean(AAS), COT = mean(COT), WholeFishWWT_g = mean(WholeFishWWT_g), HSI = mean(HSI), LiverLipids = mean(LiverLipids), MuscleLipids = mean(MuscleLipids), Growth = mean(Growth), GSI = mean(GSI), ParasiteWWT = mean(ParasiteWWT), MeanTemp_3Months = mean(MeanTemp_3Months), SalinityRatio = mean(SalinityRatio), DaysSinceFeeding = mean(DaysSinceFeeding), DaylightHours = mean(DaylightHours), SMR_AdjRsq = mean(SMR_AdjRsq) ) %>%
  mutate(Expt = as.factor(Expt))


##### Calculate Tank Means: Table 3 ####
TankMeans <- SwimData_ByFish %>%
  select(-FishID) %>%
  group_by(TankID_Expt, Temp_Trmt, Expt) %>%
  summarise_all(list(mean = ~mean(., na.rm = T), sd = ~sd(., na.rm = T), cv = ~sd(., na.rm = T)/mean(., na.rm =T),  length = ~length(.)),na.rm = TRUE)  %>%
  mutate(Expt = as.factor(Expt))

##### Calculate Experiment Means ####
ExptMeans <- SwimData_ByFish %>%
  select(-FishID, -TankID_Expt, -Temp_Trmt) %>%
  group_by(Expt) %>%
  summarise_all(list(mean = ~mean(., na.rm = T), sd = ~sd(., na.rm = T), cv = ~sd(., na.rm = T)/mean(., na.rm =T),  length = ~length(.)),na.rm = TRUE) 


##### Tank ANOVAS: Supplemental Table 1 ####
#To evaluate similarities between the two tanks at each temperature treatment in Expt. 2

Tank6 <- SwimData_ByFish %>%
  filter(TankID_Expt == "6A" | TankID_Expt == "6B") 

SMR_6_aov <- aov(SMR_scaled ~ TankID_Expt, data = Tank6)
summary(SMR_6_aov)
RMR_6_aov <- aov(RMR_scaled ~ TankID_Expt, data = Tank6)
summary(RMR_6_aov)
MMR_6_aov <- aov(MMR_scaled ~ TankID_Expt, data = Tank6)
summary(MMR_6_aov)
Ucrit_6_aov <- aov(Ucrit_mm_scaled ~ TankID_Expt, data = Tank6)
summary(Ucrit_6_aov)
FAS_6_aov <- aov(FAS_scaled ~ TankID_Expt, data = Tank6)
summary(FAS_6_aov)
AAS_6_aov <- aov(AAS_scaled ~ TankID_Expt, data = Tank6)
summary(AAS_6_aov)


Tank10 <- SwimData_ByFish %>%
  filter(TankID_Expt == "10A" | TankID_Expt == "10B") 

SMR_10_aov <- aov(SMR_scaled ~ TankID_Expt, data = Tank10)
summary(SMR_10_aov)
RMR_10_aov <- aov(RMR_scaled ~ TankID_Expt, data = Tank10)
summary(RMR_10_aov)
MMR_10_aov <- aov(MMR_scaled ~ TankID_Expt, data = Tank10)
summary(MMR_10_aov)
Ucrit_10_aov <- aov(Ucrit_mm_scaled ~ TankID_Expt, data = Tank10)
summary(Ucrit_10_aov)
FAS_10_aov <- aov(FAS_scaled ~ TankID_Expt, data = Tank10)
summary(FAS_10_aov)
AAS_10_aov <- aov(AAS_scaled ~ TankID_Expt, data = Tank10)
summary(AAS_10_aov)


Tank14 <- SwimData_ByFish %>%
  filter(TankID_Expt == "14A" | TankID_Expt == "14B") 

SMR_14_aov <- aov(SMR_scaled ~ TankID_Expt, data = Tank14)
summary(SMR_14_aov)
RMR_14_aov <- aov(RMR_scaled ~ TankID_Expt, data = Tank14)
summary(RMR_14_aov)
MMR_14_aov <- aov(MMR_scaled ~ TankID_Expt, data = Tank14)
summary(MMR_14_aov)
Ucrit_14_aov <- aov(Ucrit_mm_scaled ~ TankID_Expt, data = Tank14)
summary(Ucrit_14_aov)
FAS_14_aov <- aov(FAS_scaled ~ TankID_Expt, data = Tank14)
summary(FAS_14_aov)
AAS_14_aov <- aov(AAS_scaled ~ TankID_Expt, data = Tank14)
summary(AAS_14_aov)


#### SMR by Temperature ####

##### Linear Model: SMR vs. Temp ####

LM_SMR <- lm(SMR_scaled ~ Expt * scale(MeanTemp_3Months), data = SwimData_ByFish)
summary(LM_SMR)
#plot(LM_SMR)

##### Plot SMR ####

marginalmeans_SMR <- ggemmeans(LM_SMR , terms = c("MeanTemp_3Months", "Expt"))
marginalmeans_SMR
plot(marginalmeans_SMR)

marginalmeans_SMR_df <- as.data.frame(marginalmeans_SMR)

Trim_Expt1_SMR <- marginalmeans_SMR_df %>%
  filter(group == 1) %>%
  filter(x <= 8) %>%
  mutate(Expt = "1")
Trim_Expt2_SMR <- marginalmeans_SMR_df %>%
  filter(group == 2) %>%
  filter(x >= 5.5) %>%
  mutate(Expt = "2")

marginalmeans_SMR_Trim <- rbind(Trim_Expt1_SMR, Trim_Expt2_SMR)


SMR_Temp_Expt <- ggplot(data = TankMeans, aes(y=SMR_scaled_mean, x=MeanTemp_3Months_mean)) +
  geom_smooth(data =marginalmeans_SMR_Trim, method = lm,  aes(x = x, y = predicted, group = group, color = Expt), linewidth = 1.5)+
  geom_ribbon(data = marginalmeans_SMR_Trim, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, group = group, fill = Expt), alpha = 0.25, show.legend = FALSE) +
  geom_point(data =  SwimData_ByFish, aes(y=SMR_scaled, x=MeanTemp_3Months, shape = Expt), color = "grey50", size = 3)+
  geom_point(data = TankMeans, aes(y=SMR_scaled_mean, x = MeanTemp_3Months_mean, shape = Expt), color = "black", size = 6, alpha = 1) +
  geom_errorbarh(data = TankMeans, aes(xmin=MeanTemp_3Months_mean - MeanTemp_3Months_sd, xmax=MeanTemp_3Months_mean + MeanTemp_3Months_sd),color = "black", height=.05)+ 
  geom_errorbar(data = TankMeans, aes(ymax = SMR_scaled_mean + SMR_scaled_sd, ymin = SMR_scaled_mean - SMR_scaled_sd), color = "black", width=0.5) +
  theme_classic()+
  scale_color_manual(values=c("#402aaa", "#c7576a")) +
  scale_fill_manual(values=c("#b6abea",  "#e5b2bb")) +
  ylab(expression(bold("Standard Metabolic Rate (mgO"[2]*" kg"^-1*" h"^-1*")"))) +
  xlab(expression(bold("Temperature " ( degree~"C"))))+
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 14, face = "bold"))+
  guides(color = guide_legend(override.aes = list(size = 2))) +
  guides(shape = guide_legend(override.aes = list(size = 3))) 
print(SMR_Temp_Expt)


#### RMR by Temperature ####

##### Linear Model: RMR vs. Temp ####

LM_RMR <- lm(RMR_scaled ~ Expt * scale(MeanTemp_3Months), data = SwimData_ByFish)
summary(LM_RMR)
#plot(LM_RMR)


##### Plot RMR ####

marginalmeans_RMR <- ggemmeans(LM_RMR, terms = c("MeanTemp_3Months", "Expt"))
marginalmeans_RMR
plot(marginalmeans_RMR)

marginalmeans_RMR_df <- as.data.frame(marginalmeans_RMR)

Trim_Expt1_RMR <- marginalmeans_RMR_df %>%
  filter(group == 1) %>%
  filter(x <= 8) %>%
  mutate(Expt = "1")
Trim_Expt2_RMR <- marginalmeans_RMR_df %>%
  filter(group == 2) %>%
  filter(x >= 5.5) %>%
  mutate(Expt = "2")

marginalmeans_RMR_Trim <- rbind(Trim_Expt1_RMR, Trim_Expt2_RMR)


RMR_Temp_Expt <- ggplot(data = TankMeans, aes(y=RMR_scaled_mean, x=MeanTemp_3Months_mean)) +
  geom_smooth(data =marginalmeans_RMR_Trim, method = lm,  aes(x = x, y = predicted, group = group, color = Expt), linewidth = 1.5)+
  geom_ribbon(data = marginalmeans_RMR_Trim, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, group = group, fill = Expt), alpha = 0.25, show.legend = FALSE) +
  geom_point(data =  SwimData_ByFish, aes(y=RMR_scaled, x=MeanTemp_3Months, shape = Expt), color = "grey50", size = 3)+
  geom_point(data = TankMeans, aes(y=RMR_scaled_mean, x = MeanTemp_3Months_mean, shape = Expt), color = "black", size = 6, alpha = 1) +
  geom_errorbarh(data = TankMeans, aes(xmin=MeanTemp_3Months_mean - MeanTemp_3Months_sd, xmax=MeanTemp_3Months_mean + MeanTemp_3Months_sd),color = "black", height=.05)+ 
  geom_errorbar(data = TankMeans, aes(ymax = RMR_scaled_mean + RMR_scaled_sd, ymin = RMR_scaled_mean - RMR_scaled_sd), color = "black", width=0.5) +
  theme_classic()+
  scale_color_manual(values=c("#402aaa", "#c7576a")) +
  scale_fill_manual(values=c("#b6abea",  "#e5b2bb")) +
  ylab(expression(bold("Routine Metabolic Rate (mgO"[2]*" kg"^-1*" h"^-1*")"))) +
  xlab(expression(bold("Temperature " ( degree~"C"))))+
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 14, face = "bold"))+
  guides(color = guide_legend(override.aes = list(size = 2))) +
  guides(shape = guide_legend(override.aes = list(size = 3))) 
print(RMR_Temp_Expt)

#### MMR by Temperature ####

##### Linear Model: MMR vs. Temp ####

LM_MMR <- lm(MMR_scaled ~ Expt * scale(MeanTemp_3Months), data = SwimData_ByFish)
summary(LM_MMR)
#plot(LM_MMR)


##### Plot MMR ####

marginalmeans_MMR <- ggemmeans(LM_MMR, terms = c("MeanTemp_3Months", "Expt"))
marginalmeans_MMR
plot(marginalmeans_MMR)

marginalmeans_MMR_df <- as.data.frame(marginalmeans_MMR)

Trim_Expt1_MMR <- marginalmeans_MMR_df %>%
  filter(group == 1) %>%
  filter(x <= 8)  %>%
  mutate(Expt = "1")
Trim_Expt2_MMR <- marginalmeans_MMR_df %>%
  filter(group == 2) %>%
  filter(x >= 5.5) %>%
  mutate(Expt = "2")

marginalmeans_MMR_Trim <- rbind(Trim_Expt1_MMR, Trim_Expt2_MMR)


MMR_Temp_Expt <- ggplot(data = TankMeans, aes(y=MMR_scaled_mean, x=MeanTemp_3Months_mean)) +
  geom_smooth(data =marginalmeans_MMR_Trim, method = lm,  aes(x = x, y = predicted, group = group, color = Expt), linewidth = 1.5)+
  geom_ribbon(data = marginalmeans_MMR_Trim, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, group = group, fill = Expt), alpha = 0.25, show.legend = FALSE) +
  geom_point(data =  SwimData_ByFish, aes(y=MMR_scaled, x=MeanTemp_3Months, shape = Expt), color = "grey50", size = 3)+
  geom_point(data = TankMeans, aes(y=MMR_scaled_mean, x = MeanTemp_3Months_mean, shape = Expt), color = "black", size = 6, alpha = 1) +
  geom_errorbarh(data = TankMeans, aes(xmin=MeanTemp_3Months_mean - MeanTemp_3Months_sd, xmax=MeanTemp_3Months_mean + MeanTemp_3Months_sd),color = "black", height=.05)+ 
  geom_errorbar(data = TankMeans, aes(ymax = MMR_scaled_mean + MMR_scaled_sd, ymin = MMR_scaled_mean - MMR_scaled_sd), color = "black", width=0.5) +
  theme_classic()+
  scale_color_manual(values=c("#402aaa", "#c7576a")) +
  scale_fill_manual(values=c("#b6abea",  "#e5b2bb")) +
  ylab(expression(bold("Maximum Metabolic Rate (mgO"[2]*" kg"^-1*" h"^-1*")"))) +
  xlab(expression(bold("Temperature " ( degree~"C"))))+
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 14, face = "bold"))+
  guides(color = guide_legend(override.aes = list(size = 2))) +
  guides(shape = guide_legend(override.aes = list(size = 3))) 
print(MMR_Temp_Expt)

#### Ucrit by Temperature ####

##### Linear Model: Ucrit vs. Temp ####

LM_Ucrit <- lm(Ucrit_mm_scaled ~ Expt * scale(MeanTemp_3Months), data = SwimData_ByFish)
summary(LM_Ucrit)

#### Plot Ucrit ####

marginalmeans_Ucrit <- ggemmeans(LM_Ucrit, terms = c("MeanTemp_3Months", "Expt"))
marginalmeans_Ucrit
plot(marginalmeans_Ucrit)

marginalmeans_Ucrit_df <- as.data.frame(marginalmeans_Ucrit)

Trim_Expt1_Ucrit <- marginalmeans_Ucrit_df %>%
  filter(group == 1) %>%
  filter(x <= 8) %>%
  mutate(Expt = "1")
Trim_Expt2_Ucrit <- marginalmeans_Ucrit_df %>%
  filter(group == 2) %>%
  filter(x >= 5.5) %>%
  mutate(Expt = "2")

marginalmeans_Ucrit_Trim <- rbind(Trim_Expt1_Ucrit, Trim_Expt2_Ucrit)


Ucrit_Temp_Expt <- ggplot(data = TankMeans, aes(y=Ucrit_mm_scaled_mean, x=MeanTemp_3Months_mean)) +
  geom_smooth(data =marginalmeans_Ucrit_Trim, method = lm,  aes(x = x, y = predicted, group = group, color = Expt), linewidth = 1.5)+
  geom_ribbon(data = marginalmeans_Ucrit_Trim, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, group = group, fill = Expt), alpha = 0.25, show.legend = FALSE) +
  geom_point(data =  SwimData_ByFish, aes(y=Ucrit_mm_scaled, x=MeanTemp_3Months, shape = Expt), color = "grey50", size = 3)+
  geom_point(data = TankMeans, aes(y=Ucrit_mm_scaled_mean, x = MeanTemp_3Months_mean, shape = Expt), color = "black", size = 6, alpha = 1) +
  geom_errorbarh(data = TankMeans, aes(xmin=MeanTemp_3Months_mean - MeanTemp_3Months_sd, xmax=MeanTemp_3Months_mean + MeanTemp_3Months_sd),color = "black", height=.05)+ 
  geom_errorbar(data = TankMeans, aes(ymax = Ucrit_mm_scaled_mean + Ucrit_mm_scaled_sd, ymin = Ucrit_mm_scaled_mean - Ucrit_mm_scaled_sd), color = "black", width=0.5) +
  theme_classic()+
  scale_color_manual(values=c("#402aaa", "#c7576a")) +
  scale_fill_manual(values=c("#b6abea",  "#e5b2bb")) +
  ylab(expression(bold("Critical Swimming Speed (mm s"^-1*")"))) +
  xlab(expression(bold("Temperature " ( degree~"C"))))+
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 14, face = "bold"))+
  guides(color = guide_legend(override.aes = list(size = 2))) +
  guides(shape = guide_legend(override.aes = list(size = 3))) 
print(Ucrit_Temp_Expt)

#### FAS by Temperature ####

##### Linear Model: FAS vs. Temp ####


LM_FAS <- lm(FAS_scaled ~ Expt * scale(MeanTemp_3Months), data = SwimData_ByFish)
summary(LM_FAS)

##### Plot FAS ####

marginalmeans_FAS <- ggemmeans(LM_FAS, terms = c("MeanTemp_3Months", "Expt"))
marginalmeans_FAS
plot(marginalmeans_FAS)

marginalmeans_FAS_df <- as.data.frame(marginalmeans_FAS)

Trim_Expt1_FAS <- marginalmeans_FAS_df %>%
  filter(group == 1) %>%
  filter(x <= 8) %>%
  mutate(Expt = "1")
Trim_Expt2_FAS <- marginalmeans_FAS_df %>%
  filter(group == 2) %>%
  filter(x >= 5.5) %>%
  mutate(Expt = "2")

marginalmeans_FAS_Trim <- rbind(Trim_Expt1_FAS, Trim_Expt2_FAS)

FAS_Temp_Expt <- ggplot(data = TankMeans, aes(y=FAS_scaled_mean, x=MeanTemp_3Months_mean)) +
  geom_smooth(data =marginalmeans_FAS_Trim, method = lm,  aes(x = x, y = predicted, group = group, color = Expt), linewidth = 1.5)+
  geom_ribbon(data = marginalmeans_FAS_Trim, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, group = group, fill = Expt), alpha = 0.25, show.legend = FALSE) +
  geom_point(data =  SwimData_ByFish, aes(y=FAS_scaled, x=MeanTemp_3Months, shape = Expt), color = "grey50", size = 3)+
  geom_point(data = TankMeans, aes(y=FAS_scaled_mean, x = MeanTemp_3Months_mean, shape = Expt), color = "black", size = 6, alpha = 1) +
  geom_errorbarh(data = TankMeans, aes(xmin=MeanTemp_3Months_mean - MeanTemp_3Months_sd, xmax=MeanTemp_3Months_mean + MeanTemp_3Months_sd),color = "black", height=.05)+ 
  geom_errorbar(data = TankMeans, aes(ymax = FAS_scaled_mean + FAS_scaled_sd, ymin = FAS_scaled_mean - FAS_scaled_sd), color = "black", width=0.5) +
  theme_classic()+
  scale_color_manual(values=c("#402aaa", "#c7576a")) +
  scale_fill_manual(values=c("#b6abea",  "#e5b2bb")) +
  ylab(expression(bold("Factorial Aerobic Scope (mgO"[2]*" kg"^-1*" h"^-1*")"))) +
  xlab(expression(bold("Temperature " ( degree~"C"))))+
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 14, face = "bold"))+
  guides(color = guide_legend(override.aes = list(size = 2))) +
  guides(shape = guide_legend(override.aes = list(size = 3))) 
print(FAS_Temp_Expt)

#### AAS by Temperature ####

##### Linear Model: AAS vs. Temp ####

LM_AAS <- lm(AAS_scaled ~ Expt * scale(MeanTemp_3Months), data = SwimData_ByFish)
summary(LM_AAS)


#### Plot AAS ####

marginalmeans_AAS <- ggemmeans(LM_AAS, terms = c("MeanTemp_3Months", "Expt"))
marginalmeans_AAS
plot(marginalmeans_AAS)

marginalmeans_AAS_df <- as.data.frame(marginalmeans_AAS)

Trim_Expt1_AAS <- marginalmeans_AAS_df %>%
  filter(group == 1) %>%
  filter(x <= 8) %>%
  mutate(Expt = "1")
Trim_Expt2_AAS <- marginalmeans_AAS_df %>%
  filter(group == 2) %>%
  filter(x >= 5.5) %>%
  mutate(Expt = "2")

marginalmeans_AAS_Trim <- rbind(Trim_Expt1_AAS, Trim_Expt2_AAS)


AAS_Temp_Expt <- ggplot(data = TankMeans, aes(y=AAS_scaled_mean, x=MeanTemp_3Months_mean)) +
  geom_smooth(data =marginalmeans_AAS_Trim, method = lm,  aes(x = x, y = predicted, group = group, color = Expt), linewidth = 1.5)+
  geom_ribbon(data = marginalmeans_AAS_Trim, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, group = group, fill = Expt), alpha = 0.25, show.legend = FALSE) +
  geom_point(data =  SwimData_ByFish, aes(y=AAS_scaled, x=MeanTemp_3Months, shape = Expt), color = "grey50", size = 3)+
  geom_point(data = TankMeans, aes(y=AAS_scaled_mean, x = MeanTemp_3Months_mean, shape = Expt), color = "black", size = 6, alpha = 1) +
  geom_errorbarh(data = TankMeans, aes(xmin=MeanTemp_3Months_mean - MeanTemp_3Months_sd, xmax=MeanTemp_3Months_mean + MeanTemp_3Months_sd),color = "black", height=.05)+ 
  geom_errorbar(data = TankMeans, aes(ymax = AAS_scaled_mean + AAS_scaled_sd, ymin = AAS_scaled_mean - AAS_scaled_sd), color = "black", width=0.5) +
  theme_classic()+
  scale_color_manual(values=c("#402aaa", "#c7576a")) +
  scale_fill_manual(values=c("#b6abea",  "#e5b2bb")) +
  ylab(expression(bold("Absolute Aerobic Scope (mgO"[2]*" kg"^-1*" h"^-1*")"))) +
  xlab(expression(bold("Temperature " ( degree~"C"))))+
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 14, face = "bold"))+
  guides(color = guide_legend(override.aes = list(size = 2))) +
  guides(shape = guide_legend(override.aes = list(size = 3))) 
print(AAS_Temp_Expt)


#### Figure 1 ####
Performance_AllFish <- ggarrange(SMR_Temp_Expt, RMR_Temp_Expt, MMR_Temp_Expt, Ucrit_Temp_Expt, FAS_Temp_Expt, AAS_Temp_Expt, common.legend = T, ncol = 2, nrow = 3, legend = "right", labels = "auto")
#dev.copy(jpeg,'AllMetabolicMetrics_Temp_ALLFISH.jpg', width=11, height = 13, units='in', res=300)
#dev.off()


#### Q10: Table 4 ####
#Use the equation: Q10 = (R2/R1)^(10/(T2 - T1)), where R is a swim metric and T is the temp
#Using values from the tank means calculations

#Expt 1 SMR
Expt1_Q10_SMR<- (47.95/40.51)^(10/(8-2))
#Q10 = 1.32

#Expt 2  SMR: 6-14
Expt2_Q10_SMR <- (((90.73+ 79.22)/2) / ((74.84 + 75.36)/2))^(10/(14-6))
#Q10 = 1.17


#Expt 1 RMR:
Expt1_Q10_RMR <- (62.88/61.94 )^(10/(8-2))
#Q10 = 1.03

#Expt 2 RMR: 
Expt2_Q10_RMR <-  (((102.02+ 91.56)/2) / ((95.25+ 90.54)/2))^(10/(14-6))
#Q10 = 1.05

#Expt 1 MMR: 
Expt1_Q10_MMR <- (143.29/126.60 )^(10/(8-2))
#Q10 = 1.23

#Expt 2 MMR: 
Expt2_Q10_MMR <- (((177.77 + 162.69)/2) / ((161.01 + 169.63)/2))^(10/(14-6))
#Q10 = 1.04


#Expt 1 Ucrit:
Expt1_Q10_Ucrit <- (347.55/293.24)^(10/(8-2))
#Q10 = 1.33

#Expt 2 Ucrit:
Expt2_Q10_Ucrit <- (((406.26+ 379.78)/2) / ((343.21 + 346.97)/2))^(10/(14-6))
#Q10 = 1.18

#Expt 1 FAS:
Expt1_Q10_FAS <- (2.33/2.05  )^(10/(8-2))
#Q10 = 1.24

#Expt 2 FAS: 
Expt2_Q10_FAS <- (((1.94+ 1.94)/2) / ((1.90 + 1.75)/2))^(10/(14-6))
#Q10 = 1.08

#Expt 1 AAS: 
Expt1_Q10_AAS <- (80.40/64.66  )^(10/(8-2))
#Q10 = 1.44

#Expt 2 AAS: 
Expt2_Q10_AAS <-  (((71.13+ 75.75)/2) / ((79.09 + 65.76)/2))^(10/(14-6))
#Q10 = 1.02

#Expt 1 COT: 
Expt1_Q10_COT <- (1701.63/1713.17 )^(10/(8-2))
#Q10 = 0.99

#Expt 2 COT: 
Expt2_Q10_COT <-  (((1945.54 + 2413.64)/2) / ((2212.87 + 2135.46)/2))^(10/(14-6))
#Q10 = 1.00

#### Plot Body Mass: Supplemental Figure 1  ####

SMR_Mass <-  ggplot(data = TankMeans , aes(y=SMR_mean, x=WholeFishWWT_g_mean, color = TankID_Expt)) +
  geom_point(data = TankMeans, aes(y=SMR_mean, x= WholeFishWWT_g_mean, color =  TankID_Expt, shape = Expt), size = 6)+
  geom_point(data =  SwimData_ByFish, aes(y=SMR, x=WholeFishWWT_g, color = TankID_Expt, shape = Expt), size = 2, alpha = 0.5)+
  geom_errorbar(data = TankMeans, aes(ymin=SMR_mean - SMR_sd, ymax=SMR_mean + SMR_sd), width=.5)+ 
  geom_errorbarh(data =  TankMeans, aes(xmax = WholeFishWWT_g_mean + WholeFishWWT_g_sd, xmin = WholeFishWWT_g_mean- WholeFishWWT_g_sd), height=3) +
  scale_color_brewer(name="Tank ID", palette = "Paired") + 
  theme_classic()+
  ylab(expression(bold("Standard Metabolic Rate (mgO"[2]*"kg"^-1*"h"^-1*")"))) +
  xlab(expression(bold("Body Mass (g)")))+
  ggtitle("Standard Metabolic Rate", subtitle = "Fish: n = 71 \n Tanks: n = 10") +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  ylim(25, 150) +
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 14, face = "bold"))
print(SMR_Mass)

RMR_Mass <-  ggplot(data = TankMeans , aes(y=RMR_mean, x=WholeFishWWT_g_mean, color = TankID_Expt)) +
  geom_point(data = TankMeans, aes(y=RMR_mean, x= WholeFishWWT_g_mean, color =  TankID_Expt, shape = Expt), size = 6)+
  geom_point(data =  SwimData_ByFish, aes(y=RMR, x=WholeFishWWT_g, color = TankID_Expt, shape = Expt), size = 2, alpha = 0.5)+
  geom_errorbar(data = TankMeans, aes(ymin=RMR_mean - RMR_sd, ymax=RMR_mean + RMR_sd), width=.5)+ 
  geom_errorbarh(data =  TankMeans, aes(xmax = WholeFishWWT_g_mean + WholeFishWWT_g_sd, xmin = WholeFishWWT_g_mean- WholeFishWWT_g_sd), height=3) +
  scale_color_brewer(name="Tank ID", palette = "Paired") + 
  theme_classic()+
  ylab(expression(bold("Routine Metabolic Rate (mgO"[2]*"kg"^-1*"h"^-1*")"))) +
  xlab(expression(bold("Body Mass (g)")))+
  ggtitle("Routine Metabolic Rate", subtitle = "Fish: n = 100 \n Tanks: n = 10") +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  ylim(25, 150) +
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 14, face = "bold"))
print(RMR_Mass)

MMR_Mass <-  ggplot(data = TankMeans , aes(y=MMR_mean, x=WholeFishWWT_g_mean, color = TankID_Expt)) +
  geom_point(data = TankMeans, aes(y=MMR_mean, x= WholeFishWWT_g_mean, color =  TankID_Expt, shape = Expt), size = 6)+
  geom_point(data =  SwimData_ByFish, aes(y=MMR, x=WholeFishWWT_g, color = TankID_Expt, shape = Expt), size = 2, alpha = 0.5)+
  geom_errorbar(data = TankMeans, aes(ymin=MMR_mean - MMR_sd, ymax=MMR_mean + MMR_sd), width=.5)+ 
  geom_errorbarh(data =  TankMeans, aes(xmax = WholeFishWWT_g_mean + WholeFishWWT_g_sd, xmin = WholeFishWWT_g_mean- WholeFishWWT_g_sd), height=3) +
  scale_color_brewer(name="Tank ID", palette = "Paired") + 
  theme_classic()+
  ylab(expression(bold("Maximum Metabolic Rate (mgO"[2]*"kg"^-1*"h"^-1*")"))) +
  xlab(expression(bold("Body Mass (g)")))+
  ggtitle("Maximum Metabolic Rate", subtitle = "Fish: n = 100 \n Tanks: n = 10") +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 14, face = "bold"))
print(MMR_Mass)

Ucrit_Mass <-  ggplot(data = TankMeans , aes(y=Ucrit_mm_mean, x=WholeFishWWT_g_mean, color = TankID_Expt)) +
  geom_point(data = TankMeans, aes(y=Ucrit_mm_mean, x= WholeFishWWT_g_mean, color =  TankID_Expt, shape = Expt), size = 6)+
  geom_point(data =  SwimData_ByFish, aes(y=Ucrit_mm, x=WholeFishWWT_g, color = TankID_Expt, shape = Expt), size = 2, alpha = 0.5)+
  geom_errorbar(data = TankMeans, aes(ymin=Ucrit_mm_mean - Ucrit_mm_sd, ymax=Ucrit_mm_mean + Ucrit_mm_sd), width=.5)+ 
  geom_errorbarh(data =  TankMeans, aes(xmax = WholeFishWWT_g_mean + WholeFishWWT_g_sd, xmin = WholeFishWWT_g_mean- WholeFishWWT_g_sd), height=3) +
  scale_color_brewer(name="Tank ID", palette = "Paired") + 
  theme_classic()+
  ylab(expression(bold("Critical Swimming Speed (mm/sec)"))) +
  xlab(expression(bold("Body Mass (g)")))+
  ggtitle("Critical Swimming Speed", subtitle = "Fish: n = 100 \n Tanks: n = 10") +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 14, face = "bold"))
print(Ucrit_Mass)

FAS_Mass <-  ggplot(data = TankMeans , aes(y=FAS_mean, x=WholeFishWWT_g_mean, color = TankID_Expt)) +
  geom_point(data = TankMeans, aes(y=FAS_mean, x= WholeFishWWT_g_mean, color =  TankID_Expt, shape = Expt), size = 6)+
  geom_point(data =  SwimData_ByFish, aes(y=FAS, x=WholeFishWWT_g, color = TankID_Expt, shape = Expt), size = 2, alpha = 0.5)+
  geom_errorbar(data = TankMeans, aes(ymin=FAS_mean - FAS_sd, ymax=FAS_mean + FAS_sd), width=.5)+ 
  geom_errorbarh(data =  TankMeans, aes(xmax = WholeFishWWT_g_mean + WholeFishWWT_g_sd, xmin = WholeFishWWT_g_mean- WholeFishWWT_g_sd), height=0.05) +
  scale_color_brewer(name="Tank ID", palette = "Paired") + 
  theme_classic()+
  ylab(expression(bold("Factorial Aerobic Scope (mgO"[2]*"kg"^-1*"h"^-1*")"))) +
  xlab(expression(bold("Body Mass (g)")))+
  ggtitle("Factorial Aerobic Scope", subtitle = "Fish: n = 100 \n Tanks: n = 10") +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 14, face = "bold"))
print(FAS_Mass)


AAS_Mass <-  ggplot(data = TankMeans , aes(y=AAS_mean, x=WholeFishWWT_g_mean, color = TankID_Expt)) +
  geom_point(data = TankMeans, aes(y=AAS_mean, x= WholeFishWWT_g_mean, color =  TankID_Expt, shape = Expt), size = 6)+
  geom_point(data =  SwimData_ByFish, aes(y=AAS, x=WholeFishWWT_g, color = TankID_Expt, shape = Expt), size = 2, alpha = 0.5)+
  geom_errorbar(data = TankMeans, aes(ymin=AAS_mean - AAS_sd, ymax=AAS_mean + AAS_sd), width=.5)+ 
  geom_errorbarh(data =  TankMeans, aes(xmax = WholeFishWWT_g_mean + WholeFishWWT_g_sd, xmin = WholeFishWWT_g_mean- WholeFishWWT_g_sd), height=0.3) +
  scale_color_brewer(name="Tank ID", palette = "Paired") + 
  theme_classic()+
  ylab(expression(bold("Absolute Aerobic Scope (mgO"[2]*"kg"^-1*"h"^-1*")"))) +
  xlab(expression(bold("Body Mass (g)")))+
  ggtitle("Absolute Aerobic Scope", subtitle = "Fish: n = 100 \n Tanks: n = 10") +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 14, face = "bold"))
print(AAS_Mass)

Performance_Mass <- ggarrange(SMR_Mass, RMR_Mass, MMR_Mass, Ucrit_Mass, FAS_Mass, AAS_Mass, common.legend = T, ncol = 2, nrow = 3, legend = "right", labels = "auto")
#dev.copy(jpeg,'AllMetabolicMetrics_Mass.jpg', width=12, height = 12, units='in', res=300)
#dev.off()

#### Plot Growth: Supplemental Fig. 2 #### 

Growth_Temp_Expt <- ggplot(data = TankMeans, aes(y=Growth_mean, x=MeanTemp_3Months_mean)) +
  geom_point(data = TankMeans, aes(y=Growth_mean, x = MeanTemp_3Months_mean, color = TankID_Expt, shape = Expt), size = 6, alpha = 1) +
  geom_point(data =  SwimData_ByFish, aes(y=Growth, x=MeanTemp_3Months, color = TankID_Expt, shape = Expt), size = 2, alpha = 0.5)+
  geom_errorbarh(data = TankMeans, aes(xmin=MeanTemp_3Months_mean - MeanTemp_3Months_sd, xmax=MeanTemp_3Months_mean + MeanTemp_3Months_sd, color = TankID_Expt), height=.05)+ 
  geom_errorbar(data = TankMeans, aes(ymax = Growth_mean + Growth_sd, ymin = Growth_mean -Growth_sd, color = TankID_Expt), width=0.5) +
  theme_classic()+
  scale_color_brewer(name="Tank", palette = "Paired") + 
  ylab(expression(bold("Growth (mm/day)")))+
  xlab(expression(bold("Temperature " ( degree~"C"))))+
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 14, face = "bold"))+
  guides(color = guide_legend(override.aes = list(size = 2))) +
  guides(shape = guide_legend(override.aes = list(size = 3))) 
print(Growth_Temp_Expt )


#### Oxygen Consumption, Swimming Speed #####

#LMM between MO2 and swimming speed, with a mean MO2 for each swim speed
lmm_MO2 <- lmer(log(MO2_scaled_mean) ~ scale(SwimSpeed_mm) * scale(MeanTemp_3Months) * Expt + (1|FishID), MeanMO2, na.action = na.exclude, REML = TRUE, control = lmerControl(optimizer = "Nelder_Mead") )
model_summary_MO2 <- summary(lmm_MO2)
Anova(lmm_MO2, type = 3)
r.squaredGLMM(lmm_MO2)

#### Run again with Temp as a factor
lmm_MO2_asfactor <- lmer(log(MO2_scaled_mean) ~ scale(SwimSpeed_mm) * Temp_Expt * Expt + (1|FishID), MeanMO2, na.action = na.exclude, REML = TRUE, control = lmerControl(optimizer = "Nelder_Mead"))

#Calculate Marginal Means of MO2 ~ SwimSpeed * Expt 

mydf_MO2_swim <- ggemmeans(lmm_MO2_asfactor, terms = c("SwimSpeed_mm", "Temp_Expt"))
mydf_MO2_swim 

plot(mydf_MO2_swim )+
  theme_classic(base_size = 17, base_family = "")

## Trim lines that extend beyond data
#Min and Max Sizes
MinMax <- MeanMO2%>% 
  group_by(Temp_Expt) %>%
  reframe(min_speed = min(SwimSpeed_mm), max_speed = max(SwimSpeed_mm))

Expt1_2Degrees <- mydf_MO2_swim %>%
  filter(group == '2°C Expt 1') %>%
  filter(x<= 380) 
Expt1_4Degrees <- mydf_MO2_swim %>%
  filter(group == '4°C Expt 1') %>%
  filter(x<= 410)  
Expt1_6Degrees <- mydf_MO2_swim %>%
  filter(group == '6°C Expt 1') %>%
  filter(x<= 435)  
Expt1_8Degrees <- mydf_MO2_swim %>%
  filter(group == '8°C Expt 1') %>%
  filter(x<= 430) 
Expt2_6Degrees <- mydf_MO2_swim %>%
  filter(group == '6°C Expt 2') %>%
  filter(x<= 430) 
Expt2_10Degrees <- mydf_MO2_swim %>%
  filter(group == '10°C Expt 2') %>%
  filter(x<= 540) 
Expt2_14Degrees <- mydf_MO2_swim %>%
  filter(group == '14°C Expt 2') %>%
  filter(x<= 540) 

mydf_TRIM <- rbind(Expt1_2Degrees, Expt1_4Degrees, Expt1_6Degrees, Expt1_8Degrees, Expt2_6Degrees, Expt2_10Degrees, Expt2_14Degrees)
mydf_TRIM <- mydf_TRIM %>%
  mutate(Temp_Expt = as.factor(group)) %>%
  mutate(MO2_scaled_mean = predicted) %>%
  mutate(SwimSpeed_mm = x) %>%
  mutate(Temp= case_when( Temp_Expt ==  '2°C Expt 1' ~"2°C", Temp_Expt == '4°C Expt 1' ~"4°C", Temp_Expt == '6°C Expt 1' | Temp_Expt == '6°C Expt 2'  ~ "6°C", Temp_Expt== '8°C Expt 1' ~ "8°C", Temp_Expt ==  '10°C Expt 2'~ "10°C", Temp_Expt== '14°C Expt 2' ~ "14°C")) %>%
  mutate(Temp = factor(Temp, levels = c("2°C", "4°C", "6°C", "8°C",  "10°C", "14°C"))) %>%
  mutate(Expt = case_when( Temp_Expt ==  '2°C Expt 1' | Temp_Expt ==  '4°C Expt 1' | Temp_Expt ==  '6°C Expt 1' | Temp_Expt ==  '8°C Expt 1' ~"Expt 1", Temp_Expt == '6°C Expt 2' | Temp_Expt == '10°C Expt 2' | Temp_Expt == '14°C Expt 2'~ "Expt 2"))




##### Figure 2 ####
MO2_MarginalMeans_byTemp <- ggplot() +
  facet_wrap(~Temp, ncol = 6) +
  geom_point(data = MeanMO2, aes(x =SwimSpeed_mm, y =MO2_scaled_mean, color = Expt), alpha = 0.7) +
  geom_ribbon(data = mydf_TRIM, aes(x = SwimSpeed_mm, y = MO2_scaled_mean, ymin = conf.low, ymax = conf.high, fill = Expt, group = Temp_Expt), alpha = 0.5, show.legend = FALSE) +
  geom_smooth(data =mydf_TRIM , aes(x = SwimSpeed_mm, y = MO2_scaled_mean, color = Expt), linewidth = 1.25)+
  facet_wrap(~Temp, ncol = 6) +
  scale_color_manual(values=c("#402aaa", "#c7576a")) +
  scale_fill_manual(values=c("#b6abea",  "#e5b2bb")) +
  theme_bw() + 
  theme(
    panel.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()) +  
  theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 0.5, hjust=.5),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 14, face = "bold"), 
        strip.background = element_rect(colour="black", fill="white"),
        strip.text = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))+
  xlab(expression(bold("Swim Speed (mm s"^-1*")"))) +
  ylab(expression(bold(" Oxygen Consumption (mgO"[2]*" kg"^-1*" h"^-1*")"))) 

print(MO2_MarginalMeans_byTemp )


#### Cost of Transport ####
#Fit Non-linear models for each temp. 
#Using COT = a/SwimSpeed + (b*SwimSpeed)^(c-1)
#From Claireaux et al. 2006

##### Expt. 1: 2 degrees only ####
COT_Expt1_Data_2degrees <- SwimData %>%
  filter(Expt == 1) %>%
  filter(Temp_Trmt == 2) %>%
  filter(SMR_AdjRsq > 0.3) #Remove poor swimmers to facilitate model fit

StartingValues_a <- list(a = 3000000, b = 0.01, c = 3)
Equation <- COT ~ (a/SwimSpeed_mm) + (b * (SwimSpeed_mm^(c-1)))
COT_Expt1_2degrees <- nls(Equation,data=COT_Expt1_Data_2degrees,start=StartingValues_a)
summary(COT_Expt1_2degrees)

Uopt_Expt1_2 <- (236700/((3.047 - 1)*0.005264))^(1/3.047 )

#COT @ Uopt
(236700 /256.75) + (0.005264     * (256.75^(3.047   -1)))

#Calculate an R^2
(RSS.p1 <- sum(residuals(COT_Expt1_2degrees)^2))  # Residual sum of squares
(TSS.1 <- sum((COT_Expt1_Data_2degrees$COT- mean(COT_Expt1_Data_2degrees$COT))^2))  # Total sum of squares
1 - (RSS.p1/TSS.1)  # R-squared measure

#Get fit and confidence intervals
test <- predFit(COT_Expt1_2degrees, interval = "confidence")
interval_1 <- as.data.frame(predFit(COT_Expt1_2degrees , interval = "confidence", level= 0.95))
ModelFit_Expt1_2degrees <- cbind(COT_Expt1_Data_2degrees, interval_1) %>%
  mutate(Temp_Expt = "2°C Expt 1")


##### Expt. 1: 4 degrees only ####
COT_Expt1_Data_4degrees <- SwimData %>%
  filter(Expt == 1) %>%
  filter(Temp_Trmt == 4) %>%
  filter(SMR_AdjRsq > 0.3) #Remove poor swimmers to facilitate model fit

StartingValues_a <- list(a = 3000000, b = 0.01, c = 3)
Equation <- COT ~ (a/SwimSpeed_mm) + (b * (SwimSpeed_mm^(c-1)))
COT_Expt1_4degrees <- nls(Equation,data=COT_Expt1_Data_4degrees,start=StartingValues_a)
summary(COT_Expt1_4degrees)

Uopt_Expt1_4 <- (270000/((2.746- 1)*0.02111))^(1/2.746)

#COT @ Uopt 
(270000 /316.18) + (0.02111     * (316.18^(2.746   -1)))

#Get fit and confidence intervals
test <- predFit(COT_Expt1_4degrees, interval = "confidence")
interval_1 <- as.data.frame(predFit(COT_Expt1_4degrees , interval = "confidence", level= 0.95))
ModelFit_Expt1_4degrees <- cbind(COT_Expt1_Data_4degrees, interval_1) %>%
  mutate(Temp_Expt = "4°C Expt 1")

#Calculate an R^2
(RSS.p1 <- sum(residuals(COT_Expt1_4degrees)^2))  # Residual sum of squares
(TSS.1 <- sum((COT_Expt1_Data_4degrees$COT- mean(COT_Expt1_Data_4degrees$COT))^2))  # Total sum of squares
1 - (RSS.p1/TSS.1)  # R-squared measure

##### Expt. 1: 6 degrees only ####
COT_Expt1_Data_6degrees <- SwimData %>%
  filter(Expt == 1) %>%
  filter(Temp_Trmt == 6) %>%
  filter(SMR_AdjRsq > 0.3) #Remove poor swimmers to facilitate model fit

StartingValues_a <- list(a = 3000000, b = 0.01, c = 3)
Equation <- COT ~ (a/SwimSpeed_mm) + (b * (SwimSpeed_mm^(c-1)))
COT_Expt1_6degrees <- nls(Equation,data=COT_Expt1_Data_6degrees,start=StartingValues_a)
summary(COT_Expt1_6degrees)

Uopt_Expt1_6 <- (244700/((3.271  - 1)*0.001349 ))^(1/3.271 )

#COT @ Uopt
(244700/260.55) + (0.001349    * (260.55^(3.271  -1)))

#Get fit and confidence intervals
#Cannot get confidence intervals for this one
#Error in solve.default(crossprod(R1)) : 
# system is computationally singular: reciprocal condition number = 4.81006e-17

test <- predFit(COT_Expt1_6degrees)
interval_1 <- as.data.frame(predFit(COT_Expt1_6degrees)) %>%
  mutate(fit = predFit(COT_Expt1_6degrees)) %>%
  select(fit) %>%
  mutate(lwr = 0) %>%
  mutate(upr = 0) %>%
  mutate(Temp_Expt = "6°C Expt 1")

ModelFit_Expt1_6degrees <- cbind(COT_Expt1_Data_6degrees, interval_1)

#Calculate an R^2
(RSS.p1 <- sum(residuals(COT_Expt1_6degrees)^2))  # Residual sum of squares
(TSS.1 <- sum((COT_Expt1_Data_6degrees$COT- mean(COT_Expt1_Data_6degrees$COT))^2))  # Total sum of squares
1 - (RSS.p1/TSS.1)  # R-squared measure

##### Expt. 1: 8 degrees only ####
COT_Expt1_Data_8degrees <- SwimData %>%
  filter(Expt == 1) %>%
  filter(Temp_Trmt == 8) %>%
  filter(SMR_AdjRsq > 0.3) #Remove poor swimmers to facilitate model fit

StartingValues_a <- list(a = 3000000, b = 0.01, c = 3)
Equation <- COT ~ (a/SwimSpeed_mm) + (b * (SwimSpeed_mm^(c-1)))
COT_Expt1_8degrees <- nls(Equation,data=COT_Expt1_Data_8degrees,start=StartingValues_a)
summary(COT_Expt1_8degrees)

Uopt_Expt1_8 <- (270100/((2.893- 1)*0.008995))^(1/2.893)

#COT @ Uopt
(270100/308.24) + (0.008995   * (308.24^(2.893 -1)))

#Get fit and confidence intervals
test <- predFit(COT_Expt1_8degrees, interval = "confidence")
interval_1 <- as.data.frame(predFit(COT_Expt1_8degrees , interval = "confidence", level= 0.95))
ModelFit_Expt1_8degrees <- cbind(COT_Expt1_Data_8degrees, interval_1) %>%
  mutate(Temp_Expt = "8°C Expt 1")

#Calculate an R^2
(RSS.p1 <- sum(residuals(COT_Expt1_8degrees)^2))  # Residual sum of squares
(TSS.1 <- sum((COT_Expt1_Data_8degrees$COT- mean(COT_Expt1_Data_8degrees$COT))^2))  # Total sum of squares
1 - (RSS.p1/TSS.1)  # R-squared measure


### Expt. 2
COT_Expt2_Data <- SwimData %>%
  filter(Expt == 2) %>%
  filter(SMR_AdjRsq > 0.3) #Remove poor swimmers to facilitate model fit

##### Expt. 2: 6 degrees only ####
COT_Expt2_Data_6degrees <- SwimData %>%
  filter(Expt == 2) %>%
  filter(Temp_Trmt == 6) %>%
  filter(SMR_AdjRsq > 0.3) #Remove poor swimmers to facilitate model fit

StartingValues_a <- list(a = 3000000, b = 0.01, c = 3)
Equation <- COT ~ (a/SwimSpeed_mm) + (b * (SwimSpeed_mm^(c-1)))
COT_Expt2_6degrees <- nls(Equation,data=COT_Expt2_Data_6degrees,start=StartingValues_a)
summary(COT_Expt2_6degrees)

Uopt_Expt2_6 <- (359000/((2.645- 1)*0.04063))^(1/2.645)

#COT @ Uopt
(359000/350.31) + (0.04063  * (350.31^(2.645 -1)))

#Get fit and confidence intervals

test <- predFit(COT_Expt2_6degrees, interval = "confidence")
interval_1 <- as.data.frame(predFit(COT_Expt2_6degrees , interval = "confidence", level= 0.95))
ModelFit_Expt2_6degrees <- cbind(COT_Expt2_Data_6degrees, interval_1) %>%
  mutate(Temp_Expt = "6°C Expt 2")

#Calculate an R^2
(RSS.p1 <- sum(residuals(COT_Expt2_6degrees)^2))  # Residual sum of squares
(TSS.1 <- sum((COT_Expt2_Data_6degrees$COT- mean(COT_Expt2_Data_6degrees$COT))^2))  # Total sum of squares
1 - (RSS.p1/TSS.1)  # R-squared measure


##### Expt. 2: 10 degrees only ####
COT_Expt2_Data_10degrees <- SwimData %>%
  filter(Expt == 2) %>%
  filter(Temp_Trmt == 10) %>%
  filter(SMR_AdjRsq > 0.3) #Remove poor swimmers to facilitate model fit

StartingValues_a <- list(a = 3000000, b = 0.01, c = 2)
Equation <- COT ~ (a/SwimSpeed_mm) + (b * (SwimSpeed_mm^(c-1)))
COT_Expt2_10degrees <- nls(Equation,data=COT_Expt2_Data_10degrees,start=StartingValues_a)
#Cannot Converge

NoFit_Expt2_10degrees <- COT_Expt2_Data_10degrees %>%
  mutate(fit = 0) %>%
  mutate(lwr = 0) %>%
  mutate(upr = 0) %>%
  mutate(Temp_Expt = "10°C Expt 2")




##### Expt. 2: 14 degrees only ####
COT_Expt2_Data_14degrees <- SwimData %>%
  filter(Expt == 2) %>%
  filter(Temp_Trmt == 14) %>%
  filter(SMR_AdjRsq > 0.3) #Remove poor swimmers to facilitate model fit


StartingValues_a <- list(a = 3000000, b = 0.01, c = 2)
Equation <- COT ~ (a/SwimSpeed_mm) + (b * (SwimSpeed_mm^(c-1)))
COT_Expt2_14degrees <- nls(Equation,data=COT_Expt2_Data_14degrees,start=StartingValues_a)

#Cannot fit a c value here
#Error in nls(Equation, data = COT_Expt2_Data_14degrees, start = StartingValues_a) : singular gradient

NoFit_Expt2_14degrees <- COT_Expt2_Data_14degrees %>%
  mutate(fit = 0) %>%
  mutate(lwr = 0) %>%
  mutate(upr = 0) %>%
  mutate(Temp_Expt = "14°C Expt 2")


##### Figure 3 ####

AllCOT_byTemp <- rbind(ModelFit_Expt1_2degrees, ModelFit_Expt1_4degrees, ModelFit_Expt1_6degrees, ModelFit_Expt1_8degrees, ModelFit_Expt2_6degrees,NoFit_Expt2_10degrees ,NoFit_Expt2_14degrees ) %>%
  select(-Expt) %>%
  mutate(Temp_Expt = factor(Temp_Expt, levels = c("2°C Expt 1", "4°C Expt 1", "6°C Expt 1", "6°C Expt 2", "8°C Expt 1",   "10°C Expt 2", "14°C Expt 2"))) %>%
  mutate(Temp= case_when( Temp_Trmt ==  2 ~"2°C", Temp_Trmt == 4 ~"4°C",  Temp_Trmt == 6 ~ "6°C", Temp_Trmt == 8 ~ "8°C", Temp_Trmt == 10 ~ "10°C", Temp_Trmt == 14 ~ "14°C")) %>%
  mutate(Temp = factor(Temp, levels = c("2°C", "4°C", "6°C", "8°C",  "10°C", "14°C")))  %>%
  mutate(Expt = case_when( Temp_Expt ==  '2°C Expt 1' | Temp_Expt ==  '4°C Expt 1' | Temp_Expt ==  '6°C Expt 1' | Temp_Expt ==  '8°C Expt 1' ~"Expt 1", Temp_Expt == '6°C Expt 2' | Temp_Expt == '10°C Expt 2' | Temp_Expt == '14°C Expt 2'~ "Expt 2")) 


COT_FittedModel_byTEMP <- ggplot(data = AllCOT_byTemp, aes(x = SwimSpeed_mm, y = fit)) +
  geom_point(data = AllCOT_byTemp, aes(x = SwimSpeed_mm, y = COT, color = Expt),  alpha = 0.5) +
  geom_ribbon(data=AllCOT_byTemp, aes(x=SwimSpeed_mm, ymin= lwr, ymax=upr, group = Expt, fill = Expt), alpha=0.2, na.rm = T, show_guide = FALSE ) +
  geom_line(data = AllCOT_byTemp, aes(x = SwimSpeed_mm, y = fit, color = Expt), linewidth = 1.5) +
  facet_wrap(~Temp, ncol = 6) +
  scale_color_manual(values=c("#402aaa", "#c7576a")) +
  scale_fill_manual(values=c("#b6abea",  "#e5b2bb")) +
  theme_bw() + 
  theme(
    panel.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()) +  
  ylab(expression(bold(" Cost of Transport (J kg"^-1*" km"^-1*")"))) +
  theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 0.5, hjust=.5),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 14, face = "bold"), 
        strip.background = element_rect(colour="black", fill="white"),
        strip.text = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))+
  xlab(expression(bold("Swim Speed (mm s"^-1*")"))) +
  coord_cartesian(ylim = c(500, 7000)) 
COT_FittedModel_byTEMP 



#### Tank Correlation of Performance with State Dependent Variables: Supplemental Table 4 ####
# Only showing code for correlating one performance metric (AAS) with Expt. 1 tank 2, other correlations obtained by manipulating this code by performance metric and tank
Correlations <- SwimData_ByFish %>%
   filter(Expt == "1") %>%
  filter(TankID_Expt == "2") %>%
  select(AAS_scaled, HSI, LiverLipids, MuscleLipids, Growth,  GSI, ParasiteWWT, SalinityRatio, DaysSinceFeeding, DaylightHours ) %>%
  correlate() %>% 
  focus(AAS_scaled) 
Correlations


#### Principal Components Analysis ####
##### PCA Expt 1 #####

Expt1_PCA_Data <- SwimData_ByFish %>%
  filter(Expt == 1) %>% 
  select( -ParasiteWWT, -MeanTemp_3Months, -WholeFishWWT_g) %>%
  drop_na(LiverLipids)


Expt1_PCA <- rda(Expt1_PCA_Data[,c(18:25)], scale = TRUE)
summary(Expt1_PCA)

Loadings_Expt1 <- as.data.frame(Expt1_PCA$CA$u[,1:3])
Loadings_Vars_Expt1 <- as.data.frame(Expt1_PCA$CA$v[,1:3])
#Table S5

Expt1_Scores <- cbind(Expt1_PCA_Data, Loadings_Expt1) %>%
  mutate(PC1_Reverse = PC1 * -1)

## Expt. 1 Bi-plot
StateVars1 <- Expt1_Scores %>%
  mutate("Liver Lipids" = LiverLipids) %>%
  mutate("Muscle Lipids" = MuscleLipids) %>%
  mutate("GSI" = GSI) %>%
  mutate("Days Since Feeding" = DaysSinceFeeding) %>%
  mutate("Salinity" = SalinityRatio) %>%
  mutate("Photoperiod" = DaylightHours) %>%
  select( "Liver Lipids", HSI, "Muscle Lipids", Growth, "Salinity", "Days Since Feeding", "GSI", "Photoperiod") 
en_StateVars1 <- envfit(Expt1_PCA, StateVars1, permutations = 999, na.rm = T)
en_StateVars_coords1 = as.data.frame(scores(en_StateVars1, "vectors")) * ordiArrowMul(en_StateVars1) * 0.28

en_StateVars_coords1 <- en_StateVars_coords1 %>%
  mutate(PC1_Reverse = PC1 * -1)

Expt1 <-  ggplot(data=Expt1_Scores, aes(PC1, PC2))+
  geom_point(data=Expt1_Scores, aes(PC1, PC2, color= TankID_Expt), show.legend=T, size = 4)+
  theme_classic()+
  scale_linetype_manual(values = "solid") +
  scale_color_manual("Temp.", values=c("#a6cee3","#1f78b4", "#6a3d9a", "#33a02c")) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), data = en_StateVars_coords1, linewidth = 1, alpha = 0.8, colour = "black", arrow = arrow()) +
  geom_label(data =en_StateVars_coords1, aes(x=PC1, y = PC2), colour = "black", fontface = "bold", label = row.names(en_StateVars_coords1), position=position_jitter(0.1),  alpha = 0.6)+
  xlab("PC1 (34.3% of variance explained)") +
  ylab("PC2 (19.9% of variance explained)") +
  theme(axis.text = element_text( size = 12), axis.title = element_text(size = 14, face = "bold"),     plot.title = element_text(hjust = 0.5, face = "bold", size = 14), strip.text = element_text(face = "bold", size = 10), plot.subtitle= element_text(hjust = 0.5, size = 14), legend.text = element_text(size = 12), legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("Experiment 1", subtitle = "All Fish: n = 43")
print(Expt1)

##### PCA Expt. 2 ####
Expt2_PCA_Data <- SwimData_ByFish %>%
  filter(Expt == "2") %>% 
  select( -LiverLipids, -MuscleLipids, -MeanTemp_3Months, -WholeFishWWT_g) 


Expt2_PCA <- rda(Expt2_PCA_Data[,c(16:22)], scale = TRUE)
summary(Expt2_PCA)

Loadings_Expt2 <- as.data.frame(Expt2_PCA$CA$u[,1:3])
Loadings_Vars_Expt2 <- as.data.frame(Expt2_PCA$CA$v[,1:3])
#Table S5

Expt2_Scores <- cbind(Expt2_PCA_Data, Loadings_Expt2)

## Expt. 2 Bi-plot
StateVars2 <- Expt2_Scores %>%
  mutate("Parasites" = ParasiteWWT) %>%
  mutate("GSI" = GSI) %>%
  mutate("Days Since Feeding" = DaysSinceFeeding) %>%
  mutate("Salinity" = SalinityRatio) %>%
  mutate("Photoperiod" = DaylightHours) %>%
  select( "Parasites", HSI,  Growth, "Salinity", "Days Since Feeding", "GSI",  "Photoperiod") 
en_StateVars2 <- envfit(Expt2_PCA, StateVars2, permutations = 999, na.rm = T)
en_StateVars_coords2 = as.data.frame(scores(en_StateVars2, "vectors")) * ordiArrowMul(en_StateVars2) * .2

Expt2 <-  ggplot(data=Expt2_Scores, aes(PC1, PC2))+
  geom_point(data=Expt2_Scores, aes(PC1, PC2, color=as.factor(Temp_Trmt)), show.legend=T, size = 4)+
  theme_classic()+
  scale_linetype_manual(values = "solid") +
  scale_color_manual("Temp.", values=c("#6a3d9a", "#ff7f00", "#e31a1c")) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), data = en_StateVars_coords2, linewidth = 1, alpha = 0.8, colour = "black", arrow = arrow()) +
  geom_label(data =en_StateVars_coords2, aes(x=PC1, y = PC2), colour = "black", fontface = "bold", label = row.names(en_StateVars_coords2), position=position_jitter(0.06),  alpha = 0.6)+
  xlab("PC1 (29.8% of variance explained)") +
  ylab("PC2 (20.2% of variance explained)") +
  theme(axis.text = element_text( size = 12), axis.title = element_text(size = 14, face = "bold"),     plot.title = element_text(hjust = 0.5, face = "bold", size = 14), strip.text = element_text(face = "bold", size = 10), plot.subtitle= element_text(hjust = 0.5, size = 14), legend.text = element_text(size = 12), legend.title = element_text(size = 12, face = "bold")) +
  ggtitle("Experiment 2", subtitle = "All Fish: n = 57")
print(Expt2)

##### Figure 4 ####
PCA <- ggarrange(Expt1, Expt2, nrow = 1, ncol = 2, legend = "right", align = 'h', labels = "auto" )
#dev.copy(jpeg,'PCAs__ColoredbyTemp.jpg', width=14, height = 7, units='in', res=300)
#dev.off()

#### Tank Correlations with the PCA ####
# Supplemental Table S6, manipulate the following code for each tank 
Expt1_PCA2 <- Expt1_Scores %>% 
  filter(TankID_Expt == 8) %>%
  select(PC2, SMR_scaled, RMR_scaled, MMR_scaled, Ucrit_mm_scaled, FAS_scaled, AAS_scaled) %>%
  correlate() %>% 
  focus(PC2) 
Expt1_PCA2

Expt2_PCA2 <- Expt2_Scores %>% 
  filter(TankID_Expt == '14B') %>%
  select(PC2, SMR_scaled, RMR_scaled, MMR_scaled, Ucrit_mm_scaled, FAS_scaled, AAS_scaled) %>%
  correlate() %>% 
  focus(PC2) 
Expt2_PCA2

#### Correlations Plots ####

Expt1_SMR_PCA1 <- ggplot(Expt1_Scores, aes(x= PC1, y=SMR_scaled, color=TankID_Expt)) + 
  geom_point(size=3) + 
  theme_classic() +
  annotate(geom="text", x=0.2, y=26, size = 5, label="r = -0.14", fontface = "bold") +
  annotate(geom="text", x=0.2, y=20, size = 4, label="n = 31") +
  ylab(expression(bold("SMR (mgO"[2]*" kg"^-1*" h"^-1*")"))) +
  xlab("PC1") +
  scale_color_manual("Expt. 1 Temp.", values=c("#a6cee3","#1f78b4", "#6a3d9a","#33a02c")) +
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 14, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold")) 
Expt1_SMR_PCA1

Expt1_RMR_PCA1 <- ggplot(Expt1_Scores, aes(x= PC1, y=RMR_scaled, color=TankID_Expt)) + 
  geom_point(size=3) + 
  theme_classic() +
  annotate(geom="text", x=0.2, y=48, size = 5, label="r = -0.21", fontface = "bold") +
  annotate(geom="text", x=0.2, y=40, size = 4, label="n = 43") +
  ylab(expression(bold("RMR (mgO"[2]*" kg"^-1*" h"^-1*")"))) +
  xlab("PC1") +
  scale_color_manual("Expt. 1 Temp.", values=c("#a6cee3","#1f78b4", "#6a3d9a","#33a02c")) +
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 14, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold")) +
  theme(legend.position="none")
Expt1_RMR_PCA1

Expt1_MMR_PCA1 <- ggplot(Expt1_Scores, aes(x= PC1, y=MMR_scaled, color=TankID_Expt)) + 
  geom_point(size=3) + 
  theme_classic() +
  annotate(geom="text", x=0.2, y=110, size = 5, label="r = 0.49", fontface = "bold") +
  annotate(geom="text", x=0.2, y=98, size = 4, label="n = 43") +
  ylab(expression(bold("MMR (mgO"[2]*" kg"^-1*" h"^-1*")"))) +
  xlab("PC1") +
  scale_color_manual("Expt. 1 Temp.", values=c("#a6cee3","#1f78b4", "#6a3d9a","#33a02c")) +
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 14, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold")) +
  theme(legend.position="none")
Expt1_MMR_PCA1

Expt1_Ucrit_PCA1 <- ggplot(Expt1_Scores, aes(x= PC1, y=Ucrit_mm_scaled, color= TankID_Expt)) + 
  geom_point(size=3) + 
  theme_classic() +
  annotate(geom="text", x=0.2, y=215, size = 5, label="r = 0.55", fontface = "bold") +
  annotate(geom="text", x=0.2, y=185, size = 4, label="n = 43") +
  ylab(expression(bold("U"["crit"] * " (mm s"^-1*")"))) +
  xlab("PC1") +
  scale_color_manual("Expt. 1 Temp.", values=c("#a6cee3","#1f78b4", "#6a3d9a","#33a02c")) +
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 14, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold")) +
  theme(legend.position="none")
Expt1_Ucrit_PCA1

Expt1_FAS_PCA1 <- ggplot(Expt1_Scores, aes(x= PC1, y=FAS_scaled, color=TankID_Expt)) + 
  geom_point(size=3) + 
  theme_classic() +
  annotate(geom="text", x=0.2, y=1.7, size = 5, label="r = 0.56", fontface = "bold") +
  annotate(geom="text", x=0.2, y=1.4, size = 4, label="n = 43") +
  ylab(expression(bold("FAS (mgO"[2]*" kg"^-1*" h"^-1*")"))) +
  xlab("PC1") +
  scale_color_manual("Expt. 1 Temp.", values=c("#a6cee3","#1f78b4", "#6a3d9a","#33a02c")) +
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 14, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold")) +
  theme(legend.position="none")
Expt1_FAS_PCA1

Expt1_AAS_PCA1 <- ggplot(Expt1_Scores, aes(x= PC1, y=AAS_scaled, color=TankID_Expt)) + 
  geom_point(size=3) + 
  theme_classic() +
  annotate(geom="text", x=0.2, y=40, size = 5, label="r = 0.56", fontface = "bold") +
  annotate(geom="text", x=0.2, y=30, size = 4, label="n = 43") +
  ylab(expression(bold("AAS (mgO"[2]*" kg"^-1*" h"^-1*")"))) +
  xlab("PC1") +
  scale_color_manual("Expt. 1 Temp.", values=c("#a6cee3","#1f78b4", "#6a3d9a","#33a02c")) +
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 14, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold")) +
  theme(legend.position="none")
Expt1_AAS_PCA1

##### Expt. 2 ####

Expt2_SMR_PCA1 <- ggplot(Expt2_Scores, aes(x= PC1, y=SMR_scaled, color=as.factor(Temp_Trmt))) + 
  geom_point(size=3) + 
  theme_classic() +
  annotate(geom="text", x=0.2, y=45, size = 5, label="r = 0.26", fontface = "bold") +
  annotate(geom="text", x=0.2, y=38, size = 4, label="n = 40") +
  ylab(expression(bold("SMR (mgO"[2]*" kg"^-1*" h"^-1*")"))) +
  xlab("PC1") +
  scale_color_manual("Expt. 2 Temp.", values=c("#6a3d9a", "#ff7f00", "#e31a1c")) +
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 14, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold"))
Expt2_SMR_PCA1

Expt2_RMR_PCA1 <- ggplot(Expt2_Scores, aes(x= PC1, y=RMR_scaled, color=as.factor(Temp_Trmt))) + 
  geom_point(size=3) + 
  theme_classic() +
  annotate(geom="text", x=0.2, y=55, size = 5, label="r = 0.30", fontface = "bold") +
  annotate(geom="text", x=0.2, y=40, size = 4, label="n = 57") +
  ylab(expression(bold("RMR (mgO"[2]*" kg"^-1*" h"^-1*")"))) +
  xlab("PC1") +
  scale_color_manual("Expt. 2 Temp.", values=c("#6a3d9a", "#ff7f00", "#e31a1c")) +
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 14, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold")) +
  theme(legend.position="none")
Expt2_RMR_PCA1

Expt2_MMR_PCA1 <- ggplot(Expt2_Scores, aes(x= PC1, y=MMR_scaled, color=as.factor(Temp_Trmt)))+ 
  geom_point(size=3) + 
  theme_classic() +
  annotate(geom="text", x=0.2, y=115, size = 5, label="r = 0.32", fontface = "bold") +
  annotate(geom="text", x=0.2, y=95, size = 4, label="n = 57") +
  ylab(expression(bold("MMR (mgO"[2]*" kg"^-1*" h"^-1*")"))) +
  xlab("PC1") +
  scale_color_manual("Expt. 2 Temp.", values=c("#6a3d9a", "#ff7f00", "#e31a1c")) +
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 14, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold")) +
  theme(legend.position="none")
Expt2_MMR_PCA1

Expt2_Ucrit_PCA1 <- ggplot(Expt2_Scores, aes(x= PC1, y=Ucrit_mm_scaled, color=as.factor(Temp_Trmt))) + 
  geom_point(size=3) + 
  theme_classic() +
  annotate(geom="text", x=0.2, y=285, size = 5, label="r = 0.21", fontface = "bold") +
  annotate(geom="text", x=0.2, y=255, size = 4, label="n = 57") +
  ylab(expression(bold("U"["crit"] * " (mm s"^-1*")"))) +
  xlab("PC1") +
  scale_color_manual("Expt. 2 Temp.", values=c("#6a3d9a", "#ff7f00", "#e31a1c")) +
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 14, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold")) +
  theme(legend.position="none")
Expt2_Ucrit_PCA1

Expt2_FAS_PCA1 <- ggplot(Expt2_Scores, aes(x= PC1, y=FAS_scaled, color=as.factor(Temp_Trmt)))+ 
  geom_point(size=3) + 
  theme_classic() +
  annotate(geom="text", x=0.2, y=1, size = 5, label="r = -0.16", fontface = "bold") +
  annotate(geom="text", x=0.2, y=0.65, size = 4, label="n = 57") +
  ylab(expression(bold("FAS (mgO"[2]*" kg"^-1*" h"^-1*")"))) +
  xlab("PC1") +
  scale_color_manual("Expt. 2 Temp.", values=c("#6a3d9a", "#ff7f00", "#e31a1c")) +
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 14, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold")) +
  theme(legend.position="none")
Expt2_FAS_PCA1


Expt2_AAS_PCA1 <- ggplot(Expt2_Scores, aes(x= PC1, y=AAS_scaled, color=as.factor(Temp_Trmt))) + 
  geom_point(size=3) + 
  theme_classic() +
  annotate(geom="text", x=0.2, y=36, size = 5, label="r = 0.07", fontface = "bold") +
  annotate(geom="text", x=0.2, y=20, size = 4, label="n = 57") +
  ylab(expression(bold("AAS (mgO"[2]*" kg"^-1*" h"^-1*")"))) +
  xlab("PC1") +
  scale_color_manual("Expt. 2 Temp.", values=c("#6a3d9a", "#ff7f00", "#e31a1c")) +
  theme(axis.text.x = element_text( vjust = 0.5, size = 10),  axis.text.y = element_text( vjust = 0.5, size = 10), axis.title = element_text(size = 14, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle= element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold")) +
  theme(legend.position="none")
Expt2_AAS_PCA1

PCA_Correlations_Expt1 <- ggarrange(Expt1_SMR_PCA1, Expt1_RMR_PCA1, Expt1_MMR_PCA1, Expt1_Ucrit_PCA1, Expt1_FAS_PCA1, Expt1_AAS_PCA1, ncol = 3, nrow = 2,  common.legend = T,legend = "right", align = "v")

PCA_Correlations_Expt2 <- ggarrange(Expt2_SMR_PCA1, Expt2_RMR_PCA1, Expt2_MMR_PCA1, Expt2_Ucrit_PCA1, Expt2_FAS_PCA1, Expt2_AAS_PCA1, ncol = 3, nrow = 2,  common.legend = T,legend = "right", align = "v")

##### Figure 5 ####
PCA_Corr <- ggarrange(PCA_Correlations_Expt1, PCA_Correlations_Expt2, nrow = 2, ncol = 1, align = 'v', legend = "right")
#dev.copy(jpeg,'PCA_Correlations.jpg', width=18, height = 12, units='in', res=300)
#dev.off()



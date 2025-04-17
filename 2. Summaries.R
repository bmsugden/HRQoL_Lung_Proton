rm(list = ls())

# Load data from CSV file (1. Data Loading and Cleaning)
MergedData <- read.csv("MergedData1.csv", stringsAsFactors = TRUE, check.names = FALSE, na.strings = c("", " ","NA", "<NA>", "#NUM!"))
str(MergedData)
library(summarytools)
library(colorRamps)
library(ggpubr)
library(ggplot2)

#Treatment30 determines proton patients if at least 30% of total fractions were delivered with protons. Otherwise, patient classed as photons. 
# in case of missing in number of fractions, take corresponding "Label" value (i.e., assigned treatment)
MergedData$Treatment30 <- ifelse(
  !is.na(MergedData$Fractions_proton) & !is.na(MergedData$Fractions_total), 
  ifelse(MergedData$Fractions_proton / MergedData$Fractions_total >= 0.3, "Protons", "Photons"),
  ifelse(MergedData$Label == "Protons (PV)", "Protons",
         ifelse(MergedData$Label == "Photons (PV)", "Photons", NA))
)
#Treatment50 determines proton patients if at least 50% of total fractions were delivered with protons. Otherwise, patient classed as photons. 
# in case of missing in number of fractions, take corresponding "Label" value (i.e., assigned treatment)
MergedData$Treatment50 <- ifelse(
  !is.na(MergedData$Fractions_proton) & !is.na(MergedData$Fractions_total), 
  ifelse(MergedData$Fractions_proton / MergedData$Fractions_total >= 0.5, "Protons", "Photons"),
  ifelse(MergedData$Label == "Protons (PV)", "Protons",
         ifelse(MergedData$Label == "Photons (PV)", "Photons", NA))
)

#Treatment80 determines proton patients if at least 80% of total fractions were delivered with protons. Otherwise, patient classed as photons. 
# in case of missing in number of fractions, take corresponding "Label" value (i.e., assigned treatment)
MergedData$Treatment80 <- ifelse(
  !is.na(MergedData$Fractions_proton) & !is.na(MergedData$Fractions_total), 
  ifelse(MergedData$Fractions_proton / MergedData$Fractions_total >= 0.8, "Protons", "Photons"),
  ifelse(MergedData$Label == "Protons (PV)", "Protons",
         ifelse(MergedData$Label == "Photons (PV)", "Photons", NA))
)

#Treatment100 determines proton patients if at least 100% of total fractions were delivered with protons. Otherwise, patient classed as photons. 
# in case of missing in number of fractions, take corresponding "Label" value (i.e., assigned treatment)
MergedData$Treatment100 <- ifelse(
  !is.na(MergedData$Fractions_proton) & !is.na(MergedData$Fractions_total), 
  ifelse(MergedData$Fractions_proton / MergedData$Fractions_total >= 1.0, "Protons", "Photons"),
  ifelse(MergedData$Label == "Protons (PV)", "Protons",
         ifelse(MergedData$Label == "Photons (PV)", "Photons", NA))
)

# Create variable that splits patients into: protons (only), photons (only), and combination (of fractions as protons and photons)
MergedData$TreatmentSplit <- ifelse(
  !is.na(MergedData$Fractions_proton) & !is.na(MergedData$Fractions_total),
  ifelse(MergedData$Fractions_proton == MergedData$Fractions_total, "Protons",
         ifelse(MergedData$Fractions_proton == 0, "Photons", "Combination")),
  ifelse(MergedData$Label == "Protons (PV)", "Protons",
         ifelse(MergedData$Label == "Photons (PV)", "Photons", NA))
)



PatSummary <- dfSummary(MergedData)
sink(file = "Pat_summary.txt") #Save summary
cat("\n")
print(PatSummary)
cat("\n")
sink()

# PROMs Data summary - expected vs actual

ExpVsActual <- MergedData[, c("PatientID", "Date_StartRT", "Date_lastRT", "Periode", "Treatment30", "Survival_Status", "Survival_status.date")]
ExpVsActual$Periode <- as.factor(ExpVsActual$Periode)

ExpVsActual$Date_StartRT <- as.Date(ExpVsActual$Date_StartRT)
target_date <- as.Date("2022-10-12")
# Calculate the number of months between Date_StartRT and target_date
ExpVsActual$Months_fromStartRT <- with(ExpVsActual, 
                                   (as.numeric(format(target_date, "%Y")) - as.numeric(format(Date_StartRT, "%Y"))) * 12 +
                                     (as.numeric(format(target_date, "%m")) - as.numeric(format(Date_StartRT, "%m")))
)


#sink(file = "MissingPROMs_summary.txt") #Save Missing PROMs summary
#cat("\n")
#print(MissingPROMs)
#cat("\n")
#sink()

#sink(file = "Missing_PROMs_Proton_summary.txt") #Save Missing Proton PROMs summary
#cat("\n")
#print(MissingPROMs_Proton)
#cat("\n")
#sink()

#sink(file = "Missing_PROMs_Photon_summary.txt") #Save Missing Photon PROMs summary
#cat("\n")
#print(MissingPROMs_Photon)
#cat("\n")
#sink()

# Baseline utilities
summary(MergedData[MergedData$Periode == 0,]$`EQ5DVAS`) #Summary EQVAS
summary(MergedData[MergedData$Periode == 0,]$`GlobalHealthStatus`) #Summary GHS
summary(MergedData[MergedData$Periode == 0,]$`EQindex`) #Summary EQ5D index
sd(MergedData[MergedData$Periode == 0, ]$`EQ5DVAS`, na.rm = TRUE)
sd(MergedData[MergedData$Periode == 0, ]$`GlobalHealthStatus`, na.rm = TRUE)
sd(MergedData[MergedData$Periode == 0, ]$`EQindex`, na.rm = TRUE)

#Create new columns for baseline utility
baseline_df <- subset(MergedData, Periode == 0, select = c(PatientID, `EQ5DVAS`, GlobalHealthStatus, EQindex))
names(baseline_df) <- c("PatientID", "BaselineEQVAS", "BaselineGHS", "BaselineEQ5D")
MergedData <- merge(MergedData, baseline_df, by = "PatientID", all.x = TRUE)


#Baseline utility summaries
BaselineUtilitySummary <- dfSummary(baseline_df)
sink(file = "Baseline_utilities_summary.txt") #Save Missing Photon PROMs summary
cat("\n")
print(BaselineUtilitySummary)
cat("\n")
sink()


#Exploring considered fixed effects variables
FE_Variables <- MergedData[, c("PatientID", "Label", "Treatment30", "Treatment50", "Treatment80", 
                               "Treatment100", "TreatmentSplit", "Date_StartRT",
                               "Fractions_proton", "Fractions_photon", "Fractions_total",
                               "Stage", "Fraction_Dailyfrequency", 
                               "Pulmonary_comorbidity", "Durvalumab", "Chemo_sequence", 
                               "Smoking_status", "Surgery", "Adapted", "TumorLocation", 
                               "TumorCellType", "Radiation_DoseReceived", "Gender", 
                               "Age_startRT", "WHO_PS", "Days_StartStop", "GTV_volume",
                               "Dysphagia_BL_PHYS_CTCAE", "Dyspnea_PHYS_BL_CTCAE",
                               "BaselineEQVAS", "BaselineEQ5D", "BaselineGHS")]


str(FE_Variables)
FE_Overall_Summary <- dfSummary(FE_Variables)
sink(file = "FE_Overall_Summary.txt") #Save Missing Photon PROMs summary
cat("\n")
print(FE_Overall_Summary)
cat("\n")
sink()

UniquePat_FE_Variables <- FE_Variables[!duplicated(FE_Variables$PatientID), ]
FE_UniquePat_Summary <- dfSummary(UniquePat_FE_Variables)
sink(file = "FE_UniquePat_Summary.txt") #Save Missing Photon PROMs summary
cat("\n")
print(FE_UniquePat_Summary)
cat("\n")
sink()

ProtonPatSummary <- UniquePat_FE_Variables[UniquePat_FE_Variables$Treatment30 == "Protons", ]
PhotonPatSummary <- UniquePat_FE_Variables[UniquePat_FE_Variables$Treatment30 == "Photons", ]
view(dfSummary(ProtonPatSummary)) #Proton Patient Summary
view(dfSummary(PhotonPatSummary)) #Photon Patient Summary
view(FE_UniquePat_Summary) #Total Patient Summary

summary(ProtonPatSummary) #Summary GHS
summary(PhotonPatSummary$BaselineEQ5D) #Summary EQ5D index

# Plot combination treatments vs proton and photon only against date of start of treatment
UniquePat_FE_Variables$Date_StartRT <- as.Date(UniquePat_FE_Variables$Date_StartRT)

ggplot(UniquePat_FE_Variables, aes(x = Date_StartRT, y = TreatmentSplit)) +
  geom_jitter(width = 0, height = 0.2, alpha = 0.7) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  labs(x = "Start RT Date", y = "Treatment Split", 
       title = "Treatment Split over Time") +
  theme_minimal()



# Number of patients with a combination of protons and photons (out of 240 that data is available - 2 missing both proton patients)
sum(UniquePat_FE_Variables$Fractions_proton > 0 & UniquePat_FE_Variables$Fractions_photon > 0, 
    na.rm = TRUE) #total
sum(ProtonPatSummary$Fractions_proton > 0 & ProtonPatSummary$Fractions_photon > 0,
    na.rm = TRUE) #proton patients
sum(PhotonPatSummary$Fractions_proton > 0 & PhotonPatSummary$Fractions_photon > 0,
    na.rm = TRUE) #photon patients

missing_proton_patients <- UniquePat_FE_Variables[is.na(UniquePat_FE_Variables$Fractions_proton), ]
print(missing_proton_patients)
missing_photon_patients <- UniquePat_FE_Variables[is.na(UniquePat_FE_Variables$Fractions_photon), ]
print(missing_photon_patients)

#Baseline data (only observations at timepoint/Periode = 0)
BaselineData <- MergedData[MergedData$Periode == 0, ]
sum(duplicated(BaselineData$PatientID) | duplicated(BaselineData$PatientID, fromLast = TRUE)) # check for patients with multiple baseline values
str(MergedData)
MergedData$WHO_PS <- as.factor(MergedData$WHO_PS) 
dfSummary(BaselineData)

########################################################
################ Plot boxplots QoL data ################
give.n <- function(x){    
  return(c(y = mean(x), label = length(x)))
}
palette(rainbow(3, s = 0.5)) #color ramp 
#palette(gray.colors(3)) #color ramp 
MergedData$Periode <- as.factor(MergedData$Periode) 
dfSummary(MergedData$EQindex)

#EQ5D5L
bxp_EQindex <- ggplot(MergedData, aes(x = Periode, y = EQindex, fill = Treatment30)) +
  scale_fill_manual(values=c("1", "2", "3")) +
  geom_boxplot(alpha=0.7) + 
  stat_summary(fun.data = give.n, geom = "text", position = position_dodge(width = 0.75), size =2, aes(y=-0.5)) +
  scale_y_continuous(name = "EQ5D5L Utility", limits = c(-0.5, 1)) +
  scale_x_discrete(name = "Time (months)") +
  labs(fill = "") +
  stat_summary(fun = mean, colour="black", geom="point", shape=21, size=2, show.legend = FALSE, position = position_dodge(width = 0.75)) 

#GlobalHealthStatus
bxp_GHS <- ggplot(MergedData, aes(x = Periode, y = GlobalHealthStatus, fill = Treatment30)) +
  scale_fill_manual(values=c("1", "2", "3")) +
  geom_boxplot(alpha=0.7) + 
  stat_summary(fun.data = give.n, geom = "text", position = position_dodge(width = 0.75), size =2, aes(y=-0.5)) +
  scale_y_continuous(name = "Global Health Status score", limits = c(-0.5, 1)) +
  scale_x_discrete(name = "Time (months)") +
  labs(fill = "") +
  stat_summary(fun = mean, colour="black", geom="point", shape=21, size=2, show.legend = FALSE, position = position_dodge(width = 0.75)) 
#EQ VAS
bxp_EQVAS <- ggplot(MergedData, aes(x = Periode, y = EQ5DVAS, fill = Treatment30)) +
  scale_fill_manual(values=c("1", "2", "3")) +
  geom_boxplot(alpha=0.7) + 
  stat_summary(fun.data = give.n, geom = "text", position = position_dodge(width = 0.75), size =2, aes(y=-0.5)) +
  scale_y_continuous(name = "EQ VAS Utility", limits = c(-0.5, 1)) +
  scale_x_discrete(name = "Time (months)") +
  labs(fill = "") +
  stat_summary(fun = mean, colour="black", geom="point", shape=21, size=2, show.legend = FALSE, position = position_dodge(width = 0.75)) 

figure12 <- ggarrange(bxp_EQindex, bxp_EQVAS, bxp_GHS,
                      labels = c("", "", ""),
                      ncol = 1, nrow = 3)
figure12 #note warnings due to remove missing data

Abstract_bxp <- ggarrange(bxp_EQindex, bxp_GHS,
                          labels = c("", ""),
                          ncol = 1, nrow = 2)
Abstract_bxp #note warnings due to remove missing data

# HRQoL scores at each time point: 
## Overall
aggregate(EQindex ~ Periode, data = MergedData, FUN = mean)
aggregate(EQ5DVAS ~ Periode, data = MergedData, FUN = mean)
aggregate(GlobalHealthStatus ~ Periode, data = MergedData, FUN = mean)
## Protons:
aggregate(EQindex ~ Periode, data = subset(MergedData, Treatment30 == "Protons"), FUN = mean)
aggregate(EQ5DVAS ~ Periode, data = subset(MergedData, Treatment30 == "Protons"), FUN = mean)
aggregate(GlobalHealthStatus ~ Periode, data = subset(MergedData, Treatment30 == "Protons"), FUN = mean)
# photons:
aggregate(EQindex ~ Periode, data = subset(MergedData, Treatment30 == "Photons"), FUN = mean)
aggregate(EQ5DVAS ~ Periode, data = subset(MergedData, Treatment30 == "Photons"), FUN = mean)
aggregate(GlobalHealthStatus ~ Periode, data = subset(MergedData, Treatment30 == "Photons"), FUN = mean)

######################################

#Data required for analyses
MergedData <- MergedData[, c("PatientID", "Periode", "Label", "Treatment30", "Treatment50", "Treatment80", "Treatment100",
                             "TreatmentSplit", "Fractions_proton", "Fractions_photon", "Fractions_total",
                             "Stage", "Fraction_Dailyfrequency", "Pulmonary_comorbidity", "Durvalumab", "Chemo_sequence", 
                              "Smoking_status", "Surgery", "Adapted", "TumorLocation", 
                              "TumorCellType", "Radiation_DoseReceived", "Gender", "Age_startRT", "WHO_PS", "WHO_PS_matching", "Days_StartStop", 
                             "GTV_volume", "Dysphagia_BL_PHYS_CTCAE", 
                             "Dyspnea_PHYS_BL_CTCAE", "BaselineEQVAS", "BaselineEQ5D", "BaselineGHS", "EQindex", "GlobalHealthStatus", "EQ5DVAS")]
dfSummary(MergedData)

# Save data to CSV file
write.csv(MergedData, "MergedData2.csv", row.names = FALSE)

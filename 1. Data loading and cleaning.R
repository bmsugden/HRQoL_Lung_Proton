rm(list = ls())

library(eq5d) # See https://cran.r-project.org/web/packages/eq5d/vignettes/eq5d.html
library(summarytools)

# Load Patient and PROMs Data

PatData <- read.csv("./Database._total_09062023.csv", na.strings = c("", " ","NA"))
PROMsData <- read.csv("./Lung Cancer PROMS Maastro.csv", na.strings = c("", " ","NA"))

str(PROMsData) # inspect
str(PatData) # inspect

###### Patient Data ###### 
PatData2 <- PatData[, c(-1, -4, -5, -6, -7, -8, -9, -10, -17, -21, -22, 
                        -38, -39, -40, -41, -42, -43, -44, -45, -46, 
                        -48, -49, -50, -51, -52, -53, -54, -55, -56, 
                        -57, -58, -59, -60, -61, -62, -66, -67, -68, 
                        -69, -70, -71, -72, -76, -77, -78, -79, -80, 
                        -81, -82, -84, -85, -86, -87, -88, 
                        -90, -91, -115, -116, -117, -118, -119, -120, -121)] # select variables needed

str(PatData2) # inspect again

# Clean Patient Data (PatData)
PatData2$WHO_PS[PatData2$WHO_PS %in% c("2", "3")] <- "2"

PatData2$Stage[PatData2$Stage %in% c("IVA", "IVB")] <- "IV"
PatData2$Stage[PatData2$Stage %in% c("IIIA", "IIIB", "IIIC")] <- "III"
PatData2$Stage[PatData2$Stage %in% c("IA", "IB", "IIA", "IIB", "IIC")] <- "I/II" #COMBINE STAGE I and II

PatData2$Pulmonary_comorbidity[PatData2$Pulmonary_comorbidity %in% c("No", "no")] <- "no"
PatData2$Pulmonary_comorbidity[PatData2$Pulmonary_comorbidity %in% c("Yes", "yes")] <- "yes"
PatData2$Durvalumab[PatData2$Durvalumab %in% c("No", "No ")] <- "no"
PatData2$Durvalumab[PatData2$Durvalumab %in% c("Yes", "yes")] <- "yes"
PatData2$Chemo_sequence[PatData2$Chemo_sequence %in% c("None", "none")] <- "none"
PatData2$Chemo_sequence[PatData2$Chemo_sequence %in% c("Concurrent", "concurrent", "Induction and concurrent")] <- "concurrent"

PatData2$Surgery[PatData2$Surgery == "no"] <- "No"
PatData2$Surgery[PatData2$Surgery == "yes"] <- "Yes"
PatData2$Surgery[PatData2$Surgery == "no"] <- "No"
PatData2$Adapted[PatData2$Adapted == "Not Adapted"] <- "Not adapted"
PatData2$TumorLocation[PatData2$TumorLocation %in% c("trachea", "distale trachea", "ROK", "RMK", "RMK  ", "RMK  ", "RBK", "R", 
                                                     "voorste mediast","thymus", "retrosternale laesie")] <- "Other tumours"
PatData2$TumorLocation[PatData2$TumorLocation %in% c("LOK", "LBK", "LBK ", "LBK  ", "L", "ROK/RMK", "ROK/RBK", "ROK/LBK", "RMK/ROK", "RMK/RBK", "RBK/ROK", "RBK/RMK", 
                                                         "RBK,LBK", "LOK, RBK", "LOK, ROK", "LOK, ROK", "LBK, ROK", "LBK, LOK", "ROK, RMK", "ROK, RBK", "RBK, LBK", "LBK, RBK", "ROK, LBK")] <- "Left/Cross-lung"
PatData2$TumorLocation[PatData2$TumorLocation == "n/a"] <- NA
PatData2$TumorCellType[PatData2$TumorCellType == "n/a"] <- NA
PatData2$TumorCellType[PatData2$TumorCellType == "No PA"] <- NA
PatData2$TumorCellType[PatData2$TumorCellType %in% c("Squamous cell", "NSCLC NOS", "Adenocarcinoma")] <- "NSCLC"
PatData2$Smoking_status[PatData2$Smoking_status %in% c("Former smoker <3 months", "Former smoker >3 months", "Gestopt")] <- "Former smoker"

PatData2$Dysphagia_BL_PHYS_CTCAE[PatData2$Dysphagia_BL_PHYS_CTCAE == "0"] <- "No"
PatData2$Dysphagia_BL_PHYS_CTCAE[PatData2$Dysphagia_BL_PHYS_CTCAE == "1"] <- "Yes"
PatData2$Dysphagia_BL_PHYS_CTCAE[PatData2$Dysphagia_BL_PHYS_CTCAE == "2"] <- "Yes"

PatData2$Dyspnea_PHYS_BL_CTCAE[PatData2$Dyspnea_PHYS_BL_CTCAE == "0"] <- "No"
PatData2$Dyspnea_PHYS_BL_CTCAE[PatData2$Dyspnea_PHYS_BL_CTCAE == "1"] <- "Yes"
PatData2$Dyspnea_PHYS_BL_CTCAE[PatData2$Dyspnea_PHYS_BL_CTCAE == "2"] <- "Yes"

str(PatData2)


# Restructure patient data
PatData2$Label <- as.factor(PatData2$Label)
PatData2$Date_StartRT <- as.Date(PatData2$Date_StartRT, format = "%d/%m/%Y")
PatData2$Date_lastRT <- as.Date(PatData2$Date_lastRT, format = "%d/%m/%Y")
PatData2$Gender <- as.factor(PatData2$Gender)
PatData2$Smoking_status <- as.factor(PatData2$Smoking_status)
PatData2$Pulmonary_comorbidity <- as.factor(PatData2$Pulmonary_comorbidity)
PatData2$TumorCellType <- as.factor(PatData2$TumorCellType)
PatData2$T.stage <- as.factor(PatData2$T.stage)
PatData2$N.stage <- as.factor(PatData2$N.stage)
PatData2$M.stage <- as.factor(PatData2$M.stage)
PatData2$Surgery <- as.factor(PatData2$Surgery)
PatData2$Chemo_sequence <- as.factor(PatData2$Chemo_sequence)
PatData2$Specific_chemotherapeutic <- as.factor(PatData2$Specific_chemotherapeutic)
PatData2$PostRT_IMT <- as.factor(PatData2$PostRT_IMT)
PatData2$Durvalumab <- as.factor(PatData2$Durvalumab)
PatData2$Adapted <- as.factor(PatData2$Adapted)
PatData2$Stage <- as.factor(PatData2$Stage)
PatData2$proton.adapted <- as.factor(PatData2$proton.adapted)
PatData2$photon.adapted <- as.factor(PatData2$photon.adapted)
PatData2$Fraction_Dailyfrequency <- as.factor(PatData2$Fraction_Dailyfrequency)
PatData2$DateLastFU <- as.Date(PatData2$DateLastFU)
PatData2$X1yr.survival.date <- as.Date(PatData2$X1yr.survival.date, format = "%d/%m/%Y")
PatData2$X2yr.survival.date <- as.Date(PatData2$X2yr.survival.date, format = "%d/%m/%Y")
PatData2$Observed.1yr.mortality <- as.factor(PatData2$Observed.1yr.mortality)
PatData2$Observed.2yr.mortality <- as.factor(PatData2$Observed.2yr.mortality)
PatData2$Date_death <- as.Date(PatData2$Date_death)
PatData2$Survival_Status <- as.factor(PatData2$Survival_Status)
PatData2$Survival_status.date <- as.Date(PatData2$Survival_status.date, format = "%d/%m/%Y")
PatData2$WHO_PS <- as.factor(PatData2$WHO_PS)
PatData2$TumorLocation <- as.factor(PatData2$TumorLocation)

str(PatData2)
###### PROMs Data ###### 

# Cleaning: EQ5D responses

PROMsData$Vraagstelling[PROMsData$Vraagstelling %in% c("EQ5D 1. MOBILITEIT", "Mobiliteit", "EQ5D Mobiliteit")] <- "EQ-5D Mobility"

PROMsData$Vraagstelling[PROMsData$Vraagstelling %in% c("EQ5D 2. ZELFZORG", "Zelfzorg", "EQ5D Zelfzorg")] <- "EQ-5D Self care"

PROMsData$Vraagstelling[PROMsData$Vraagstelling %in% c("EQ5D 3. DAGELIJKSE ACTIVITEITEN (bijv. werk, studie, huishouden, gezins- en ", 
                                                         "Dagelijkse activiteiten", "EQ5D Dagelijkse activiteiten")] <- "EQ-5D Daily activities"

PROMsData$Vraagstelling[PROMsData$Vraagstelling %in% c("EQ5D 4. PIJN/ONGEMAK", 
                                                         "Pijn - ongemak", "EQ5D Pijn")] <- "EQ-5D Pain"

PROMsData$Vraagstelling[PROMsData$Vraagstelling %in% c("EQ5D 5. ANGST/SOMBERHEID", 
                                                         "Stemming", "EQ5D Stemming")] <- "EQ-5D Anxiety/depression"

PROMsData$Vraagstelling[PROMsData$Vraagstelling %in% c("EQ5D VAS", "Gezondheidsthermometer", 
                                                       "EQ5D VAS: We willen weten hoe goed of slecht uw gezondheid VANDAAG is (0-100).")] <- "EQ-5D VAS"

str(PROMsData)


# Filter the dataset to include EQ5D5L, VAS, and GHS PROMs only
PROMsData <- subset(PROMsData, Vraagstelling %in% c("EQ-5D Mobility", 
                                                     "EQ-5D Self care", 
                                                     "EQ-5D Daily activities", 
                                                     "EQ-5D Pain", 
                                                     "EQ-5D Anxiety/depression", 
                                                     "EQ-5D VAS", 
                                                     "EQ-5D Profiel", 
                                                     "GlobalHealthStatus"))

PROMsData$Periode[PROMsData$Periode %in% c("tijdens RT")] <- "1"
PROMsData$Periode[PROMsData$Periode %in% c("Baseline")] <- "0"
PROMsData$Periode[PROMsData$Periode %in% c("0 - 3 md")] <- "3"
PROMsData$Periode[PROMsData$Periode %in% c("4 - 6 md")] <- "6"
PROMsData$Periode[PROMsData$Periode %in% c("1 jr")] <- "12"
PROMsData$Periode[PROMsData$Periode %in% c("2 jr")] <- "24"
PROMsData$Periode[PROMsData$Periode %in% c("3 jr")] <- "36"
PROMsData$Periode[PROMsData$Periode %in% c("4 jr")] <- "48"
PROMsData$Periode[PROMsData$Periode %in% c("5 jr")] <- "60"
PROMsData$Periode[PROMsData$Periode %in% c(">5 jr")] <- ">60"

# Formatting data
PROMsData$PBD.aanmaakdatum <- as.Date(PROMsData$PBD.aanmaakdatum)
PROMsData$Datum.laatste.bestraling <- as.Date(PROMsData$Datum.laatste.bestraling, format = "%Y-%m-%d %H:%M:%S")
PROMsData$Datum.eerste.bestraling <- as.Date(PROMsData$Datum.eerste.bestraling, format = "%Y-%m-%d %H:%M:%S")
PROMsData$Vragenlijstcode <- as.factor(PROMsData$Vragenlijstcode)
PROMsData$Vraagtype <- as.factor(PROMsData$Vraagtype)
PROMsData$Zorglijn <- as.factor(PROMsData$Zorglijn)
PROMsData$Zorgplan <- as.factor(PROMsData$Zorgplan)
PROMsData$Baseline.Follow.up <- as.factor(PROMsData$Baseline.Follow.up)
PROMsData$Vraagcode <- as.factor(PROMsData$Vraagcode)
PROMsData$Vragenlijst <- as.factor(PROMsData$Vragenlijst)
PROMsData$Vraagstelling <- as.factor(PROMsData$Vraagstelling)
PROMsData$PBD <- as.factor(PROMsData$PBD)
PROMsData$Periode <- as.factor(PROMsData$Periode)
PROMsData$Invoerdatum <- as.Date(PROMsData$Invoerdatum, format = "%Y-%m-%d %H:%M:%S")

PROMsData2 <- PROMsData[, c(-2, -3, -4, -8, -9, -10, -12, -15, -16, -18)]
anyDuplicated(PROMsData2)

PROMsData2 <- subset(PROMsData2, Vragenlijst %in% c("Longen", "Longen BL", "Longen FU")) #PROMs relating to treatment for lung cancer only

# View the first few rows of the filtered dataset
head(PROMsData2)
str(PROMsData2)

dfSummary(PROMsData2)

#Reshape PROMs

reshaped_PROMs <- reshape(PROMsData2, idvar = c("Patientnr", "Periode", "PBD", "PBD.aanmaakdatum", "Datum.eerste.bestraling", 
                                                "Dagen.tussen.Laatste.RT.en.CTC", "Vragenlijst"), 
                          timevar = "Vraagstelling", 
                          direction = "wide",
                          v.names = "Antwoord") #Added additional identifiers to prevent duplicates

reshaped_PROMs <- reshaped_PROMs[, c(-2, -3, -5, -6)] #remove planning variables
reshaped_PROMs <- reshaped_PROMs[!duplicated(reshaped_PROMs), ] # Skip duplicates entered multiple times (i.e., PROMs were recorded for each treatment plan)
BaselinePROMs <- reshaped_PROMs[reshaped_PROMs$Periode == 0, ]
sum(duplicated(BaselinePROMs$PatientID) | duplicated(BaselinePROMs$PatientID, fromLast = TRUE)) # check for patients with multiple baseline values

names(reshaped_PROMs)[names(reshaped_PROMs) == "Antwoord.EQ-5D Mobility"] <- "MO" # Required for utility package
names(reshaped_PROMs)[names(reshaped_PROMs) == "Antwoord.EQ-5D Self care"] <- "SC" # Required for utility package
names(reshaped_PROMs)[names(reshaped_PROMs) == "Antwoord.EQ-5D Daily activities"] <- "UA" # Required for utility package
names(reshaped_PROMs)[names(reshaped_PROMs) == "Antwoord.EQ-5D Pain"] <- "PD" # Required for utility package
names(reshaped_PROMs)[names(reshaped_PROMs) == "Antwoord.EQ-5D Anxiety/depression"] <- "AD" # Required for utility package

reshaped_PROMs$MO[reshaped_PROMs$MO == 9] <- NA
reshaped_PROMs$SC[reshaped_PROMs$SC == 9] <- NA
reshaped_PROMs$UA[reshaped_PROMs$UA == 9] <- NA
reshaped_PROMs$PD[reshaped_PROMs$PD == 9] <- NA
reshaped_PROMs$AD[reshaped_PROMs$AD == 9] <- NA

reshaped_PROMs$`MO` <- as.numeric(reshaped_PROMs$`MO`)
reshaped_PROMs$`SC` <- as.numeric(reshaped_PROMs$`SC`)
reshaped_PROMs$`UA` <- as.numeric(reshaped_PROMs$`UA`)
reshaped_PROMs$`PD` <- as.numeric(reshaped_PROMs$`PD`)
reshaped_PROMs$`AD` <- as.numeric(reshaped_PROMs$`AD`)

levels(reshaped_PROMs$PBD)

names(reshaped_PROMs)[names(reshaped_PROMs) == "Antwoord.EQ-5D VAS"] <- "EQ5DVAS"
names(reshaped_PROMs)[names(reshaped_PROMs) == "Antwoord.GlobalHealthStatus"] <- "GlobalHealthStatus"
names(reshaped_PROMs)[names(reshaped_PROMs) == "Datum.eerste.bestraling"] <- "Date_StartRT"

colnames(reshaped_PROMs)[colnames(reshaped_PROMs) == "Patientnr"] <- "PatientID" # To match PatData2 for Merging

anyDuplicated(reshaped_PROMs) #Check for duplicates



### Merge datasets ###
MergedData <- merge(reshaped_PROMs, PatData2, by = c("PatientID", "Date_StartRT")) #to handle patients treated are separate time points and to handle PROMs being inputted multiple times for each treatment plan
MergedData <- MergedData[MergedData$Label %in% c("Protons (PV)", "Photons (PV)"), ] # Remove photon patients without planning comparisons
length(unique(MergedData$PatientID))
 
### Utility calculations  ###
# EQ VAS and GHS utilities
MergedData$EQ5DVAS <- as.numeric(as.character(MergedData$EQ5DVAS))
MergedData$EQ5DVAS[MergedData$EQ5DVAS == 999] <- NA
MergedData$EQ5DVAS <- MergedData$EQ5DVAS / 100 #Derive EQ VAS Utility Score
MergedData$GlobalHealthStatus <- as.numeric(as.character(MergedData$GlobalHealthStatus))
MergedData$GlobalHealthStatus <- MergedData$GlobalHealthStatus / 100 #Derive GHS Utility

# EQ5D Utilities (Dutch Tariff) - EuroQol inputs

MergedData$EQindex <- eq5d(MergedData, country="Netherlands", version="5L", type="VT", ignore.invalid=TRUE)
str(MergedData)

# Preparation for gen matching
dfSummary(MergedData)
MergedData$WHO_PS_matching <- MergedData$WHO_PS #Create new variable for matching
MergedData$WHO_PS_matching <- as.character(MergedData$WHO_PS_matching)
MergedData$WHO_PS_matching[MergedData$WHO_PS_matching == "0"] <- "0/1" #WHO PS of 0 now "0/1"
MergedData$WHO_PS_matching[MergedData$WHO_PS_matching == "1"] <- "0/1" #WHO PS of 1 now "0/1"
MergedData$WHO_PS_matching[MergedData$WHO_PS_matching == "2"] <- "2+" #WHO PS of 2+ now "2+"
MergedData$WHO_PS_matching <- as.factor(MergedData$WHO_PS_matching)


# Save data to CSV file
write.csv(MergedData, "MergedData1.csv", row.names = FALSE)



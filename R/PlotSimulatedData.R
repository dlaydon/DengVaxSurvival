### This script makes summary plots of simulated data.

require(here)


MDValue = 9999
ProjectDirectory 			= here()
R_ScriptDirectory			= file.path(ProjectDirectory, "R")
CppRootDirectory 			= here()
CppOutputDirectory 			= file.path(CppRootDirectory, "ParamFiles", "Output") 
ProcessedDataDirectory 		= file.path(CppRootDirectory, "ParamFiles", "Data") 
ProcessedOutputDirectory	= file.path(ProjectDirectory, "ProcessedOutput")

source(file.path(R_ScriptDirectory, "SubsetData.R"))
source(file.path(R_ScriptDirectory, "SurvivalCurves.R"))

Data = read.table(file = file.path(ProcessedDataDirectory, "SimData.txt"			), header = T, sep = "\t") 

AddProportions = function(SDFrame, NAMES)
{
	SDFrame[, paste0("Prop", NAMES)] = SDFrame[, NAMES] / SDFrame[, c("StrataSizes")] 
	return(SDFrame)
}

CYD_14_indices		= which(Data$Country <= 4)	### CYD14
CYD_15_indices		= which(Data$Country >= 5)	### CYD15
Data_CYD_14			= 	Data[CYD_14_indices, ]
Data_CYD_15			= 	Data[CYD_15_indices, ]
str(Data_CYD_14)
str(Data_CYD_15)


length(CYD_14_indices)
length(CYD_15_indices)
Trial_Num			= rep(NA, dim(Data)[1])
Trial_Arm			= rep(NA, dim(Data)[1])
Trial_Arm[as.logical(Data$Arm)]		= "Vaccine"
Trial_Arm[!as.logical(Data$Arm)]	= "Control"
Trial_Num[CYD_14_indices] = "CYD14"
Trial_Num[CYD_15_indices] = "CYD15"



### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** 
### *** ### *** ### *** 	CASES BY ARM, TRIAL, AGE GROUP and PHASE
### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** 

### Num Cases (any serotype/phase/severity) by arm and age group for each trial.
ReducedDFrame 	= cbind(StrataSizes = rep(1, dim(Data)[1]), Cases = Data$IsCase, ActiveCases = Data$IsCaseActive, PassiveCases = Data$IsCasePassive)
colnames(ReducedDFrame)
FactorNames 	= c("Trial", "Arm")
Summary_NoAge 	= aggregate(x = ReducedDFrame, by = list(Trial = Trial_Num, Arm = Data$Arm), FUN = sum)
Summary_NoAge	= cbind(Summary_NoAge[, FactorNames], AgeGroup = rep(0, dim(Summary_NoAge)[1]), Summary_NoAge[, colnames(ReducedDFrame)])
Summary_wAge 	= aggregate(x = ReducedDFrame, by = list(Trial = Trial_Num, Arm = Data$Arm, AgeGroup = Data$AgeGroup1), FUN = sum)
SummaryTable 	= rbind(Summary_NoAge, Summary_wAge)
SummaryTable 	= AddProportions(SummaryTable, NAMES = c("Cases", "ActiveCases", "PassiveCases"))

Mat_Active_CYD14 	= 100 * matrix(SummaryTable[SummaryTable$Trial == "CYD14", "PropActiveCases"], ncol = 4, byrow = FALSE, dimnames = list(c("Control", "Vaccine"), c("All ages", "2-5 years", "6-11 years", "12-14 years")))
Mat_Active_CYD15 	= 100 * matrix(SummaryTable[SummaryTable$Trial == "CYD15", "PropActiveCases"], ncol = 3, byrow = FALSE, dimnames = list(c("Control", "Vaccine"), c("All ages", "9-11 years", "12-16 years")))
Mat_Passive_CYD14 	= 100 * matrix(SummaryTable[SummaryTable$Trial == "CYD14", "PropPassiveCases"], ncol = 4, byrow = FALSE, dimnames = list(c("Control", "Vaccine"), c("All ages", "2-5 years", "6-11 years", "12-14 years")))
Mat_Passive_CYD15 	= 100 * matrix(SummaryTable[SummaryTable$Trial == "CYD15", "PropPassiveCases"], ncol = 3, byrow = FALSE, dimnames = list(c("Control", "Vaccine"), c("All ages", "9-11 years", "12-16 years")))

png(file = file.path(ProjectDirectory, "SimDataFigures", "CasesByPhaseAgeGroup.png"), res = 300, units = "in", width = 10, height = 8) ### individual ones are , width = 7, height= 5 
YLAB = "Cases (%)"; XLAB = "Age Group"; COLS = c("lightgreen", "darkblue"); CEXLAB = 1.35
par (mfrow = c(2,2))
barplot(Mat_Active_CYD14 , main = "Active phase cases (CYD14)"	, ylab = YLAB, ylim = c(0, 12)	, border = NA, xlab= XLAB, col=COLS, beside=TRUE, legend.text = TRUE, args.legend = list(x = "topright", cex = 1.1), cex.lab = CEXLAB, space = c(0, 0.5)) 
barplot(Mat_Active_CYD15 , main = "Active phase cases (CYD15)"	, ylab = YLAB, ylim = c(0, 7)	, border = NA, xlab= XLAB, col=COLS, beside=TRUE, legend.text = TRUE, args.legend = list(x = "topright", cex = 1.1), cex.lab = CEXLAB, space = c(0, 0.5)) 
barplot(Mat_Passive_CYD14, main = "Passive phase cases (CYD14)", ylab = YLAB, ylim = c(0, 1)	, border = NA, xlab= XLAB, col=COLS, beside=TRUE, legend.text = TRUE, args.legend = list(x = "topright", cex = 1.1), cex.lab = CEXLAB, space = c(0, 0.5)) 
barplot(Mat_Passive_CYD15, main = "Passive phase cases (CYD15)", ylab = YLAB, ylim = c(0, 0.4)	, border = NA, xlab= XLAB, col=COLS, beside=TRUE, legend.text = TRUE, args.legend = list(x = "topright", cex = 1.1), cex.lab = CEXLAB, space = c(0, 0.5)) 
dev.off()




### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** 
### *** ### *** ### *** 	ImSub Figures
### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** 

Data_CYD_14_ImSub	= 	Data_CYD_14[Data_CYD_14$SeroStatus != MDValue, ]
Data_CYD_15_ImSub	= 	Data_CYD_15[Data_CYD_15$SeroStatus != MDValue, ]
YLAB = "percent (%)"; XLAB = "Baseline serostatus"; ImSubCols = c("skyblue2", "skyblue3", "skyblue4"); CEXLAB = 1.35; CEXMAIN = 1.5
png(file = file.path(ProjectDirectory, "SimDataFigures", "BaselineSerosatusByTrial.png"), res = 300, units = "in", width = 10, height = 5) ### individual ones are , width = 7, height= 5 
par (mfrow = c(1,2))
barplot(table(Data_CYD_14$SeroStatus)/dim(Data_CYD_14)[1] * 100, main = "CYD14", border = NA, ylab = YLAB, xlab= XLAB, names = c("seronegative","seropositive","unknownn"), col = ImSubCols, ylim = c(0,100), cex.lab = CEXLAB, cex.main = CEXMAIN)
barplot(table(Data_CYD_15$SeroStatus)/dim(Data_CYD_15)[1] * 100, main = "CYD15", border = NA, ylab = YLAB, xlab= XLAB, names = c("seronegative","seropositive","unknownn"), col = ImSubCols, ylim = c(0,100), cex.lab = CEXLAB, cex.main = CEXMAIN)
dev.off()
png(file = file.path(ProjectDirectory, "SimDataFigures", "BaselineSerosatusByTrial_ImSubOnly.png"), res = 300, units = "in", width = 10, height = 5) ### individual ones are , width = 7, height= 5 
CEXMAIN = 1.3   
par (mfrow = c(1,2))
barplot(table(Data_CYD_14_ImSub$SeroStatus)/dim(Data_CYD_14_ImSub)[1] * 100, main = "CYD14: Immunogenicity subset", border = NA, ylab = YLAB, xlab= XLAB, names = c("seronegative", "seropositive"), col = ImSubCols, ylim = c(0,100), cex.lab = CEXLAB, cex.main = CEXMAIN)
barplot(table(Data_CYD_15_ImSub$SeroStatus)/dim(Data_CYD_15_ImSub)[1] * 100, main = "CYD15: Immunogenicity subset", border = NA, ylab = YLAB, xlab= XLAB, names = c("seronegative", "seropositive"), col = ImSubCols, ylim = c(0,100), cex.lab = CEXLAB, cex.main = CEXMAIN)
dev.off()



### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** 
### *** ### *** ### *** 	CASES BY TRIAL AND SEROTYPE
### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** 

Data$CaseSeroType[which(Data$CaseSeroType == 9999)] = NA #### don't tabulate non-cases. 
CasesBySerotype_CYD14 			= table(Data_CYD_14$CaseSeroType)[1:4]  #### don't tabulate non-cases.
CasesBySerotype_CYD15 			= table(Data_CYD_15$CaseSeroType)[1:4]  #### don't tabulate non-cases.
names(CasesBySerotype_CYD14) 	= paste0("serotype\n", 1:4) 
names(CasesBySerotype_CYD15) 	= paste0("serotype\n", 1:4) 

YLAB = ""; XLAB = ""; SerotypeCols = c("coral", "coral1", "coral2", "coral3"); CEXLAB = 1.35; CEXMAIN = 1.5
png(file = file.path(ProjectDirectory, "SimDataFigures", "CasesByTrialAndSerotype.png"), res = 300, units = "in", width = 12, height = 5)  
par (mfrow = c(1,2))
barplot(CasesBySerotype_CYD14, main = "Cases by serotype (CYD14)", border = NA, ylab = YLAB, xlab= XLAB, col = SerotypeCols, cex.lab = CEXLAB, cex.main = CEXMAIN, ylim = c(0,250))
barplot(CasesBySerotype_CYD15, main = "Cases by serotype (CYD15)", border = NA, ylab = YLAB, xlab= XLAB, col = SerotypeCols, cex.lab = CEXLAB, cex.main = CEXMAIN, ylim = c(0,250))
dev.off()




### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** 
### *** ### *** ### *** 	CASES BY ARM, TRIAL AND Disease Severity
### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** 

Data_MildSevere 			= read.table(file = file.path(ProcessedDataDirectory, "SimData_Severe.txt")	, header = T, sep = "\t")
Data_MildSevere_Hosp 		= read.table(file = file.path(ProcessedDataDirectory, "SimData_Hosp.txt")	, header = T, sep = "\t")
Data_MildSevere_CYD14 		= Data_MildSevere[CYD_14_indices, ]
Data_MildSevere_CYD15 		= Data_MildSevere[CYD_15_indices, ]
Data_MildSevere_Hosp_CYD14 	= Data_MildSevere_Hosp[CYD_14_indices, ]
Data_MildSevere_Hosp_CYD15 	= Data_MildSevere_Hosp[CYD_15_indices, ]

ReducedDFrame 	= cbind(StrataSizes = rep(1, dim(Data)[1]), Cases = Data$IsCase, 
		NonHospCases 	= Data_MildSevere_Hosp	$IsCaseMild, HospitalisedCases 	= Data_MildSevere_Hosp	$IsCaseSevere,
		NonSevereCases 	= Data_MildSevere		$IsCaseMild, SevereCases 		= Data_MildSevere		$IsCaseSevere)

FactorNames 	= c("Trial", "Arm")
SummaryTable_Disease 	= aggregate(x = ReducedDFrame, by = list(Trial = Trial_Num, Arm = Data$Arm), FUN = sum)
SummaryTable_Disease 	= AddProportions(SummaryTable_Disease, NAMES = colnames(ReducedDFrame))

NamesForMatSelection 	= paste0("Prop", c("Cases", "NonHospCases", "HospitalisedCases", "NonSevereCases", "SevereCases"))
NamesForPlot 			= c("any", "non-\nhospitalised", "hospitalised", "non-\nsevere", "severe")
Mat_Disease_CYD14 		= as.matrix(SummaryTable_Disease[SummaryTable_Disease$Trial == "CYD14", NamesForMatSelection])
Mat_Disease_CYD15 		= as.matrix(SummaryTable_Disease[SummaryTable_Disease$Trial == "CYD15", NamesForMatSelection])
rownames(Mat_Disease_CYD14) = c("Control", "Vaccine")
colnames(Mat_Disease_CYD14) = NamesForPlot
rownames(Mat_Disease_CYD15) = c("Control", "Vaccine")
colnames(Mat_Disease_CYD15) = NamesForPlot
YLAB = "Cases (%)"; XLAB = "Disease severity"; DiseaseCols = c("indianred1", "red4"); CEXLAB = 1.35
png(file = file.path(ProjectDirectory, "SimDataFigures", "CasesByDiseaseSeverity.png"), res = 300, units = "in", width = 14, height = 5) ### individual ones are , width = 7, height= 5 
par (mfrow = c(1,2))
barplot(Mat_Disease_CYD14, main = "Cases by disease severity (CYD14)", border = NA, cex.main = 1.5, ylab = YLAB, xlab = XLAB, col=DiseaseCols, beside=TRUE, legend.text = TRUE, args.legend = list(x = "topright", cex = 1.5), cex.lab = CEXLAB, space = c(0, 0.5), ylim = c(0, 0.1))
barplot(Mat_Disease_CYD15, main = "Cases by disease severity (CYD15)", border = NA, cex.main = 1.5, ylab = YLAB, xlab = XLAB, col=DiseaseCols, beside=TRUE, legend.text = TRUE, args.legend = list(x = "topright", cex = 1.5), cex.lab = CEXLAB, space = c(0, 0.5), ylim = c(0, 0.1))
dev.off()

png(file = file.path(ProjectDirectory, "SimDataFigures", "CasesByPhaseAgeGroup.png"), res = 300, units = "in", width = 10, height = 8) ### individual ones are , width = 7, height= 5 
par (mfrow = c(2,2))
barplot(Mat_Active_CYD14 , main = "Active phase cases (CYD14)"	, border = NA, ylab = YLAB, ylim = c(0, 12)	, xlab= XLAB, col=COLS, beside=TRUE, legend.text = TRUE, args.legend = list(x = "topright", cex = 1.1), cex.lab = CEXLAB, space = c(0, 0.5)) 
barplot(Mat_Active_CYD15 , main = "Active phase cases (CYD15)"	, border = NA, ylab = YLAB, ylim = c(0, 7)	, xlab= XLAB, col=COLS, beside=TRUE, legend.text = TRUE, args.legend = list(x = "topright", cex = 1.1), cex.lab = CEXLAB, space = c(0, 0.5)) 
barplot(Mat_Passive_CYD14, main = "Passive phase cases (CYD14)", border = NA, ylab = YLAB, ylim = c(0, 1)	, xlab= XLAB, col=COLS, beside=TRUE, legend.text = TRUE, args.legend = list(x = "topright", cex = 1.1), cex.lab = CEXLAB, space = c(0, 0.5)) 
barplot(Mat_Passive_CYD15, main = "Passive phase cases (CYD15)", border = NA, ylab = YLAB, ylim = c(0, 0.4)	, xlab= XLAB, col=COLS, beside=TRUE, legend.text = TRUE, args.legend = list(x = "topright", cex = 1.1), cex.lab = CEXLAB, space = c(0, 0.5)) 
dev.off()













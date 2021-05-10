## This scripts generates various plots and summary tables from Cpp output. 

# Load dependencies
require(shape)
require(Hmisc)
require(ggplot2)
require(sp)
require(png)
library(vioplot)
require(here)

options(width = 172L)
ParDefaults 					= par()

############################################################################
OrigMAI = par ("mai")
OrigMAR = par ("mar")

## In R indices, hazard groups are: 1 = SeroPos Control, 2 = SeroNeg Control, 3 = SeroPos Vaccine, 4 = SeroNeg Vaccine
## 5 = Control,  6 = Vaccine,  7 = SeroPositive, 8 = SeroNegative, 
HazStrings 					= c("SPosControl"	, "SNegControl"	, "SPosVaccine"	, "SNegVaccine"	, "Control"	, "Vaccine"	, "SeroPositive", "SeroNegative", "" )
HazStrings_short 			= c("SPosCon"		, "SNegCon"		, "SPosVac"		, "SNegVac"		, "Con"		, "Vac"		, "SPos"		, "SNeg"		, "" )
HazardGroupNumbers 			= 1:9
names(HazardGroupNumbers) 	= HazStrings ; names(HazardGroupNumbers)[9] = "All"
CountryNames_Long			= c("Indonesia", "Malaysia", "Philippines", "Thailand", "Vietnam", "Brazil", "Colombia", "Honduras", "Mexico", "Puerto Rico", "CYD-14", "CYD-15", "")
NoCountries 				= 10

SummStat = function(PostSamples, IncludeMode = FALSE) ### returns summary statistics of vector
{
	Summary = c(mean(PostSamples), median(PostSamples), quantile(PostSamples, c(0.025,0.975), na.rm = T))
	NAMES 	= c("Mean", "Median", "LowerCrI", "UpperCrI")
	if (IncludeMode) 
	{
		Summary = c(Summary	, mlv(PostSamples, method = "mfv")$M)
		NAMES 	= c(NAMES	, "Mode")
	}
	names(Summary) = NAMES
	return (Summary)
}
GetTransparentColours 		= function(ColourString, Alpha = 0.5)
{
	Cols 			= col2rgb(ColourString, alpha = FALSE) ### colours to be set.
	TransparentCols = rgb(red = Cols[1,]/255, green = Cols[2,]/255, blue = Cols[3,]/255, alpha = Alpha)
	
	return (TransparentCols)
}
ChooseCountryBitOfPNGTitle 	= function(country)
{
	if (country <  10)	CountryBitOfPNGTitle = paste0("_c", country) 	else 
	if (country == 10)	CountryBitOfPNGTitle = "_CYD14" 				else 
	if (country == 11)	CountryBitOfPNGTitle = "_CYD15" 				else 
	if (country == 12)	CountryBitOfPNGTitle = "_CYD1415" 
	
	return(CountryBitOfPNGTitle)
}
Alpha_ParamGuessCols = alpha("pink", 0.02)
CloseOpenPlotDevices = function()	if (!is.null(dev.list()))	for (openplot in 1:length(dev.list())) dev.off()
PlotPolygon = function(XValues, LowerLine, UpperLine, MeanLine = NULL, PolyCOL = "pink", MiddleCol = "red", ALPHA = 0.5, LWD = 4, type = "l")
{
	polygon	(c(XValues,rev(XValues)),c(LowerLine,rev(UpperLine)), col= PolyCOL, border = NA)
	if (!is.null(MeanLine)) lines(x = XValues, y = MeanLine, col = MiddleCol, lwd = LWD, type = "b")
}

PosteriorCol 	= "green"
PNG_res 		= 300
MDValue 		= 9999

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##  
## ## ## ## ## ## 		Decide project, dataset and model variant for which to produce plots 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##  

RemoveNansFromAttackTables 			= function (Attack_Table)
{
	#### Save rownames
	ROWNAMES = rownames(Attack_Table)
	
	#### find indices with dodgy parts to them. Reset to zero. Note the arr.ind part. 
	Attack_Table[which(Attack_Table == "-nan(ind)", arr.ind = TRUE)] = 0
	
	#### make matrix into numeric matrix rather than characters.
	Attack_Table = as.data.frame(lapply(Attack_Table, as.numeric)) 
	
	#### add rownames back in
	rownames(Attack_Table) = ROWNAMES
	
	return(Attack_Table)
}
Choose_WBIC_String 					= function(WBICs)
{
	if (class(WBICs) != "logical") stop ("Choose_WBIC_String error: class(WBICs) != logical")
	if (WBICs) WBIC_String = "WBIC_" else WBIC_String = ""
	return(WBIC_String)
}
GetImSubString 						= function(ImSub)
{
	if (!ImSub) ImSubString = "" else if (ImSub) ImSubString = "_ImSub" else stop("GetImSubString error: ImSub argument not logical/bool")
	return (ImSubString)
}
ImportSurvivalTables 				= function(modelrun) 
{
	### function imports all survival tables: i) create correct filename; ii) assign filename to correct object in global environment, if file exists. 
	
	PASSIVE_SURVIVAL_TABLES_EXIST	 			= FALSE
	IMMUNE_SUBSET_SURVIVAL_TABLES_EXIST 		= FALSE
	IMMUNE_SUBSET_PASSIVE_SURVIVAL_TABLES_EXIST = FALSE
	SUMMARY_SURVIVAL_TABLES_EXIST 				= FALSE
	WBIC_SUMMARY_SURVIVAL_TABLES_EXIST 			= FALSE
	
	if (DoAttackRates || DoSurvCurves)
		for (WBIC_String in c("", "WBIC_"))
			for (ImSubString in c("", "_ImSub"))
				for (whichdiseaseplot in WhichDiseasePlots)
					for (TrialPhaseString in c("", "Passive"))
						for (StatString in c("Mean_", "LowerCI_", "UpperCI_", "Modal_", "MaxLike_"))
						{
							### Stitch together FILENAME
							if (StatString %in% c("Mean_", "LowerCI_", "UpperCI_")) #### naming inconsistent between different tables. 
								FILENAME = file.path(CppOutputDirectory, paste0(WBIC_String, StatString, TrialPhaseString, "SurvivalTable", ImSubString, OutputString, whichdiseaseplot, ".txt"))
							else 
								FILENAME = file.path(CppOutputDirectory, paste0(WBIC_String, StatString, TrialPhaseString, "SurvivalTable", whichdiseaseplot, ImSubString, OutputString, ".txt"))
							### Stitch together VariableName - e.g. VariableName = Mean_PassiveSurvivalTable_Severe
							VariableName = paste0(WBIC_String, StatString, TrialPhaseString, "SurvivalTable", whichdiseaseplot, ImSubString)
							
							### import filename and assign it to variable name (if filename exists)
							if (file.exists(FILENAME))
							{
								assign(x = VariableName, value = read.table(file = FILENAME, header = TRUE, sep = "\t", row.names = 1), envir = .GlobalEnv)
								
								if (StatString %in% c("Modal_", "MaxLike_")) 							SUMMARY_SURVIVAL_TABLES_EXIST 		= TRUE
								if (WBIC_String == "WBIC_" & StatString %in% c("Modal_", "MaxLike_")) 	WBIC_SUMMARY_SURVIVAL_TABLES_EXIST 	= TRUE
								
								if (TrialPhaseString == "Passive") PASSIVE_SURVIVAL_TABLES_EXIST = TRUE
								if (ImSubString == "_ImSub") 
								{
									IMMUNE_SUBSET_SURVIVAL_TABLES_EXIST = TRUE
									if (TrialPhaseString == "Passive") IMMUNE_SUBSET_PASSIVE_SURVIVAL_TABLES_EXIST = TRUE
								}
							}
						}
	PASSIVE_SURVIVAL_TABLES_EXIST 				<<- PASSIVE_SURVIVAL_TABLES_EXIST
	IMMUNE_SUBSET_SURVIVAL_TABLES_EXIST 		<<- IMMUNE_SUBSET_SURVIVAL_TABLES_EXIST 		
	IMMUNE_SUBSET_PASSIVE_SURVIVAL_TABLES_EXIST <<- IMMUNE_SUBSET_PASSIVE_SURVIVAL_TABLES_EXIST
	SUMMARY_SURVIVAL_TABLES_EXIST 				<<- SUMMARY_SURVIVAL_TABLES_EXIST
	WBIC_SUMMARY_SURVIVAL_TABLES_EXIST 			<<- WBIC_SUMMARY_SURVIVAL_TABLES_EXIST
}
ImportSeroPrevTables 				= function()
{
	if (DoSeroPrevs)
	{
		MEANSEROPREVFILE 		= file.path(CppOutputDirectory, paste0("Mean_SeroPrevOutput"	, OutputString, ".txt"))
		MEDIANSEROPREVFILE 		= file.path(CppOutputDirectory, paste0("Median_SeroPrevOutput"	, OutputString, ".txt"))
		LOWERCISEROPREVFILE 	= file.path(CppOutputDirectory, paste0("LowerCI_SeroPrevOutput"	, OutputString, ".txt"))
		UPPERCISEROPREVFILE 	= file.path(CppOutputDirectory, paste0("UpperCI_SeroPrevOutput"	, OutputString, ".txt"))
		
		if (!file.exists(MEANSEROPREVFILE)) return()
		
		Mean_SeroPrevOutput	 	= read.table(file = MEANSEROPREVFILE 	, header = T, sep = "\t", row.names = 1)
		Median_SeroPrevOutput 	= read.table(file = MEDIANSEROPREVFILE 	, header = T, sep = "\t", row.names = 1)
		LowerCI_SeroPrevOutput 	= read.table(file = LOWERCISEROPREVFILE , header = T, sep = "\t", row.names = 1)
		UpperCI_SeroPrevOutput 	= read.table(file = UPPERCISEROPREVFILE , header = T, sep = "\t", row.names = 1)
		colnames(Mean_SeroPrevOutput)
		ColsToDrop = dim(Mean_SeroPrevOutput)[2]
		
		Mean_SeroPrevOutput	 	= Mean_SeroPrevOutput	 [, -ColsToDrop]
		Median_SeroPrevOutput 	= Median_SeroPrevOutput	 [, -ColsToDrop]
		LowerCI_SeroPrevOutput 	= LowerCI_SeroPrevOutput [, -ColsToDrop]
		UpperCI_SeroPrevOutput 	= UpperCI_SeroPrevOutput [, -ColsToDrop]
		
		Mean_SeroPrevOutput	 	<<- Mean_SeroPrevOutput	 	
		Median_SeroPrevOutput 	<<- Median_SeroPrevOutput 	
		LowerCI_SeroPrevOutput 	<<- LowerCI_SeroPrevOutput 	
		UpperCI_SeroPrevOutput 	<<- UpperCI_SeroPrevOutput 	
	}
}
ImportAttackRatePosteriorSamples 	= function(modelrun)
{
	if (DoAttackRates)
	{
		for (WBIC_String in c("", "WBIC_"))
			for (ImSubString in c("", "_ImSub"))
				for (MeanModeString in c("", "_MeanMode"))
					for (whichdiseaseplot in WhichDiseasePlots)
						for (whichtrialphase in 1:length(WhichTrialPhases))
						{
							### Create RootAttackFileName, to which you add. 
							RootAttackFileName = file.path(CppOutputDirectory, paste0(WBIC_String, "AttackRates", MeanModeString, ImSubString, OutputString))
							
							#### create logical values as to whether that specific table exists. Assign to false first, amend later if filename exists
							assign(x = paste0(sub("_", "", WhichTrialPhases[whichtrialphase]), "AttackFilesExist"), value = FALSE, envir = .GlobalEnv) #### e.g. AttackFilesExist 	= FALSE - will be reset if file exists. 
							
							#### create FileName from RootAttackFileName
							FileName = paste0(RootAttackFileName, sub("Only", "PhaseOnly", WhichTrialPhases[whichtrialphase]), whichdiseaseplot, ".txt")
							
							### if file exists import file and remove NaNs.  
							if (file.exists(FileName)) #### e.g. AttackRateTable_PassiveOnly_Mild	= read.table(file = paste0(CppOutputDirectory, "AttackRates", OutputString, "_PassivePhaseOnly_Mild.txt")	, header = T, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
								assign(x = paste0(WBIC_String, "AttackRateTable", WhichTrialPhases[whichtrialphase], whichdiseaseplot, MeanModeString, ImSubString), 
										value = RemoveNansFromAttackTables(read.table(file = FileName, header = T, sep = "\t", row.names = 1, stringsAsFactors = FALSE))	, 
										envir = .GlobalEnv)
							
							if (file.exists(FileName))
								assign(x = paste0(sub("_", "", WhichTrialPhases[whichtrialphase]), MeanModeString, ImSubString, "AttackFilesExist"), value = TRUE, envir = .GlobalEnv) #### e.g. PassiveOnlyAttackFilesExist 	= TRUE
						}
		if (modelrun$phase == "_ACTIVE") ParamSurvivePostSampleRatio = 10 else	### horribly hacky, written in late on sunday night so you could go home. Try to delete. 
		if (AttackFilesExist)
			ParamSurvivePostSampleRatio <<- round(dim(CHAINS)[1] / dim(get(paste0("AttackRateTable", WhichDiseasePlots[1])))[2] ) ## e.g. ParamSurvivePostSampleRatio = round(dim(CHAINS)[1] / dim(AttackRateTable_Mild)[2])
	}
}
GetHazRatioTableName 				= function(StatString = "Mean", HR_serotype = 1, DiseaseSeverity = "", ImSubString = "", WBIC_String = "")
{
	return(paste0(WBIC_String, StatString, "_HazRatios", HRs_SeroNames[HR_serotype], DiseaseSeverity, ImSubString))
}
ImportHazRatPosts					= function()
{
	for (WBIC_String in c(""))#, "WBIC_")) ### add in the WBIC later if you need it. 
		for (ImSubString in c("", "_ImSub"))
			for (DiseaseSeverity in WhichDiseasePlots)
				for (StatString in c("Mean", "LowerCI", "UpperCI"))
					for (HR_serotype in 1:HRs_NumSTypes)
					{
						FILENAME 		= file.path(CppOutputDirectory, paste0(WBIC_String, StatString, "_HazRatios", HRs_SeroNames[HR_serotype], ImSubString, OutputString, DiseaseSeverity, ".txt"))
						VariableName	= GetHazRatioTableName(StatString = StatString, HR_serotype = HR_serotype, DiseaseSeverity = DiseaseSeverity, ImSubString = ImSubString, WBIC_String = WBIC_String)
						if (file.exists(FILENAME))
							assign(x = VariableName, value = read.table(file = FILENAME, header = TRUE, sep = "\t", row.names = 1), envir = .GlobalEnv) else 
							warning(paste0("ImportHazRatPosts error: ", FILENAME, " does not exist"))
					}
}
ShouldCountryBeSkipped 				= function(modelrun)
{
	TorF = FALSE ## default is that country shouldn't be skipped. 
	
	#### skipping coniditons
	
	## i.e. if countries are pooled, only plot full trials. 
	if (modelrun$PooledCountries == 1 	& country < 10	) TorF = TRUE 					else 
	if (modelrun$PooledTrials == 1 		& country != 12	) TorF = TRUE 					else #### may be something dodgy about this as I can't remember if PooledTrials outputs to "country" 10, or "country" 12
	
	if (country <  10 && (!any(country == CountriesFitted))) 	TorF = TRUE 	else  
	#### for CYD14 and CYD15 (need their own conditions as "countries" 10, 11, and 12 are not countries but trials/pooled trials, and so will never be in CountriesFitted. 
	if (country == 10 & !any(1:4 %in% CountriesFitted))  		TorF = TRUE		else  # i.e. if plotting for CYD14 but no CYD14 countries were fitted in this model run, skip!  
	if (country == 11 & !any(5:9 %in% CountriesFitted))  		TorF = TRUE		else  # i.e. if plotting for CYD15 but no CYD15 countries were fitted in this model run, skip!  
	
	if (modelrun$PooledTrials == 1 		& country != 12	) cat ("Check your output country index for PooledTrials - think either 10 or 12. ")
	
	return (TorF)
}

ProjectDirectory 			= here()
R_ScriptDirectory			= file.path(ProjectDirectory, "R")
CppRootDirectory 			= here()
CppOutputDirectory 			= file.path(CppRootDirectory, "ParamFiles", "Output") 
ProcessedDataDirectory 		= file.path(CppRootDirectory, "ParamFiles", "Data") 
ProcessedOutputDirectory	= file.path(ProjectDirectory, "ProcessedOutput")

source(file.path(R_ScriptDirectory, "DirectoriesEtc.R"))
source(file.path(R_ScriptDirectory, "Layout.R"))
source(file.path(R_ScriptDirectory, "AttackRates.R"))
source(file.path(R_ScriptDirectory, "SubsetData.R"))
source(file.path(R_ScriptDirectory, "Splines.R"))
source(file.path(R_ScriptDirectory, "SurvivalCurves.R"))

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##  
## ## ## ## ## ## 		What to plot?
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##  

DoPosteriors 				= TRUE				 
DoChains 					= TRUE
DoAttackRates				= TRUE
DoSurvCurves				= TRUE			
DoParamCorrelations			= TRUE	
DoSeroPrevs 				= TRUE
DoHazRatios					= TRUE			
DoAgeEffectPlots			= TRUE			
DoRunSummaries				= TRUE			

ShortenOutputString = TRUE

ModelRuns = DefineModelRuns( 
		
#		ModelVariants 									= c("K_SEROPOS")		, 
#		MildSeveres 									= c("_MILDSEVERE")		,
#		Hosp_And_NonHosp 								= c("", "_hosp")		,
		SS_VEs_AND_Not									= 1						,
		ASVE_Options 									= c("CATEGORICAL")		, #   "INDEPENDENT", "SPLINE", "CATEGORICAL", "CUBIC", "HILL" 
		AS_Haz_Options 									= c("CATEGORICAL")		, #   "INDEPENDENT", "SPLINE", "CATEGORICAL", "CUBIC", "HILL" 
)
ModelRuns$OutputFolderNames
length(ModelRuns$OutputFolderNames)

MR_index = 1
ModelRun = ModelRuns[MR_index,]
ModelRun$OutputFolderNames

FilenamesCompleted = c()

for (MR_index in 1:dim(ModelRuns)[1])
{
	# remove output from previous model runs
	if (exists("DIC_table"					)) rm(DIC_table					)
	if (exists("WBIC_CHAINS"				)) rm(WBIC_CHAINS				)
	if (exists("TotalLogLikeChain"			)) rm(TotalLogLikeChain			)
	if (exists("WBIC_TotalLogLikeChain"		)) rm(WBIC_TotalLogLikeChain	)
	if (exists("SPrev_LL_Chains"			)) rm(SPrev_LL_Chains			)
	if (exists("SPrev_LL_Chains_ImSub"		)) rm(SPrev_LL_Chains_ImSub		)
	if (exists("WBIC_SPrev_LL_Chains"		)) rm(WBIC_SPrev_LL_Chains		)
	if (exists("WBIC_SPrev_LL_Chains_ImSub"	)) rm(WBIC_SPrev_LL_Chains_ImSub)
	
	CloseOpenPlotDevices()	 ## close any open plot devices
	ModelRun = ModelRuns[MR_index,]
	DefineVariousQuantitiesForModelRun(ModelRun)
	DefineFollowUpQuantities		(ModelRun)
	
	OutputString 	= ChooseOutputString(ModelRun, Folder = FALSE	, PrintToConsole = FALSE)
	OutputSubDir 	= ChooseOutputString(ModelRun, Folder = TRUE	, PrintToConsole = FALSE)
	cat(paste0("\n", MR_index, "/", dim(ModelRuns)[1], "\tOutputSubDir ", OutputSubDir, "\n"))
	
	#### get rid of duplicates 
	if (any(OutputString == FilenamesCompleted)) cat(paste0("Skipping ", OutputString, " done already\n"))
	if (any(OutputString == FilenamesCompleted)) next
	
	
	##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
	##### ##### ##### 		Import Chains
	
	ChainFileName 						= file.path(CppOutputDirectory, paste0("ParameterChainOutput" 		, OutputString, ".txt"))
	WBIC_ChainFileName 					= file.path(CppOutputDirectory, paste0("WBIC_ParameterChainOutput" 	, OutputString, ".txt"))
	SPrev_LL_Chain_FileName 			= file.path(CppOutputDirectory, paste0("SPrev_LL_Chain" 			, OutputString, ".txt"))
	SPrev_LL_Chain_ImSub_FileName 		= file.path(CppOutputDirectory, paste0("SPrev_LL_Chain_ImSub" 		, OutputString, ".txt"))
	WBIC_SPrev_LL_Chain_FileName 		= file.path(CppOutputDirectory, paste0("WBIC_SPrev_LL_Chain" 		, OutputString, ".txt"))
	WBIC_SPrev_LL_Chain_ImSub_FileName 	= file.path(CppOutputDirectory, paste0("WBIC_SPrev_LL_Chain_ImSub" 	, OutputString, ".txt"))
	TotalLogLikeChain_FileName			= file.path(CppOutputDirectory, paste0("TotalLogLikeChain" 			, OutputString, ".txt"))
	WBIC_TotalLogLikeChain_FileName		= file.path(CppOutputDirectory, paste0("WBIC_TotalLogLikeChain" 	, OutputString, ".txt"))
	
	if (file.exists(SPrev_LL_Chain_FileName 			)) SPrev_LL_Chains 				= read.table(file = SPrev_LL_Chain_FileName 			, header = T, sep = "\t")
	if (file.exists(SPrev_LL_Chain_ImSub_FileName 		)) SPrev_LL_Chains_ImSub 		= read.table(file = SPrev_LL_Chain_ImSub_FileName 		, header = T, sep = "\t")
	if (file.exists(WBIC_SPrev_LL_Chain_FileName 		)) WBIC_SPrev_LL_Chains 		= read.table(file = WBIC_SPrev_LL_Chain_FileName 		, header = T, sep = "\t")
	if (file.exists(WBIC_SPrev_LL_Chain_ImSub_FileName 	)) WBIC_SPrev_LL_Chains_ImSub 	= read.table(file = WBIC_SPrev_LL_Chain_ImSub_FileName 	, header = T, sep = "\t")
	if (file.exists(TotalLogLikeChain_FileName			)) TotalLogLikeChain			= read.table(file = TotalLogLikeChain_FileName			, header = T, sep = "\t")[,1] #### imported as a data fram hence the [,1] as you'll use it here as a vector. 
	if (file.exists(WBIC_TotalLogLikeChain_FileName		)) WBIC_TotalLogLikeChain		= read.table(file = WBIC_TotalLogLikeChain_FileName		, header = T, sep = "\t")[,1]
	
	### if results not done, skip fit 
	if(!file.exists(ChainFileName)) cat("\n", OutputSubDir, " does not exist ", "\n")
	if(!file.exists(ChainFileName)) next
	
	CHAINS 					= read.table(file = ChainFileName, header = T, sep = "\t")	; ChainNames = colnames(CHAINS) 
	if (file.exists(WBIC_ChainFileName)) 
		WBIC_CHAINS 		= read.table(file = WBIC_ChainFileName			, header = T, sep = "\t")	; ChainNames = colnames(CHAINS) 
	dim(CHAINS)
	colnames(CHAINS)
	NoParams = dim(CHAINS)[2] - 1 ##### minus one because of log likelihood. 
	Acceptance_Filename 		= file.path(CppOutputDirectory, paste0("AcceptanceArray" 		, OutputString, ".txt"))
	AcceptanceProbs				= try(read.table(file = Acceptance_Filename, header = T, sep = "\t"), silent = TRUE)	
	WBIC_Acceptance_Filename 	= file.path(CppOutputDirectory, paste0("WBIC_AcceptanceArray" 	, OutputString, ".txt"))
	if (file.exists(WBIC_Acceptance_Filename)) 
		WBIC_AcceptanceProbs 	= try(read.table(file = WBIC_Acceptance_Filename, header = T, sep = "\t"), silent = TRUE)	
	if (class(AcceptanceProbs) == "try-error") cat (paste0("AcceptanceProb: ", OutputString, "dodgy, skipping to next run"))
	if (class(AcceptanceProbs) == "try-error") next
	
	
	#### Find 95%CrIs of likelihood (used for SeroPrevs, Predicted hazards and ASVEs)
	#### for ASVEs etc, will need to change these to WBIC specific ones. Best to put them in functions even though small waste. 
	LikeCI 					= quantile(CHAINS$logLikelihood, probs = c(0.025, 0.975), na.rm = FALSE)
	WhichParamsToPlot 		= which(CHAINS$logLikelihood >= LikeCI[1] & CHAINS$logLikelihood <= LikeCI[2])
	MaxLikeParam 			= which(CHAINS$logLikelihood == max(CHAINS$logLikelihood))
	HowManyParamSetsToPlot 	= length(WhichParamsToPlot)
	MeanParamValues 		= colMeans(CHAINS[,1:NoParams]) ## Get the means of all parameters in the chain. 
	
	if (ModelRun$SS_VEs & ModelRun$SType_Equiv)
	{
		### find all rho params and ensure they are always set to 0.25.
		RhoNames = grep("rho", colnames(CHAINS))
		if (all(CHAINS[, RhoNames] == 0.25)) cat (paste0("Rhos cool ", OutputString, "\n")) else for (i in 1:10) cat (paste0("Rhos AWFUL ", OutputString, "\n"))
	}
	
	##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
	##### ##### ##### Import Survival, SeroPrev and attack rate tables
	
	ImportSurvivalTables(ModelRun)
	ImportAttackRatePosteriorSamples(ModelRun)
	ImportSeroPrevTables()

	##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
	##### ##### ##### Import simulated data. 
	
	if (ModelRun$disease == "")
	{
		if (ModelRun$SS_VEs || ModelRun$SS_Ks)
			DATA = read.table(file = file.path(ProcessedDataDirectory, "SimData.txt"			), header = T, sep = "\t") else 
			DATA = read.table(file = file.path(ProcessedDataDirectory, "SimData_NoSerotype.txt"), header = T, sep = "\t")  
		
	} else 
	{
		if (ModelRun$hosp == "_hosp")
			DATA = read.table(file = file.path(ProcessedDataDirectory, "SimData_Hosp.txt"		), header = T, sep = "\t") else 
			DATA = read.table(file = file.path(ProcessedDataDirectory, "SimData_Severe.txt"		), header = T, sep = "\t")  
		
		if (DoAttackRates)
		{
			#### Also must import non-mild-severe data for attack rate calculations. For active phase (say), need to know which cases were in active phase in addition to which cases were severe etc.
			FilenameNonMildSevere 	= file.path(ProcessedDataDirectory, "SimData.txt")
			DATA_nonMildSevere 		= read.table(file = FilenameNonMildSevere, header = T, sep = "\t")
		}
	}
	
	## Import DIC / LL
	DIC_Filename = file.path(CppOutputDirectory, paste0("DIC" , OutputString, ".txt"))
	if (file.exists(DIC_Filename)) 
		DIC_table = read.table(file = DIC_Filename, header = TRUE, sep = "\t", row.names = 1)
	
	
	##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
	##### ##### ##### 		Make Plots
	##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
	
	### Create directory for this model run
	ModelPlotSubDir = file.path(ProcessedOutputDirectory, OutputSubDir)
	dir.create(ModelPlotSubDir, showWarnings = FALSE)
	
	if (DoRunSummaries)
	{
		cat ("RunSummaries ")
		
		CalcBIC = function() #### leave as function of global variables CHAINS and DATA for now. 
		{
			### Get number of data points. Use subsampling 
			N_DataPoints 	= length(ChooseRowIndices_Country(country = 12, CountriesFitted)) 
			### Get Maximum Likelihood 
			MaxLogLikelihood 	= max(CHAINS$logLikelihood)
			#### Get Number of parameters -
			N_Params		= dim(CHAINS)[2] - length(grep("rho_", colnames(CHAINS))) - 1 ### you add rho values to CHAINS, even though qvals are the actual parameters, and - 1 so you don't include log likelihood. 
			
			BIC  = round((log(N_DataPoints) * N_Params) - 2 * MaxLogLikelihood, 2)
			return (BIC)
		}
		CalcWBIC = function()
		{
			if (file.exists(WBIC_ChainFileName)) return (round(mean(WBIC_CHAINS[, "logLikelihood"]), 2) ) else return (NA) 
		}
		
		ChooseAgeEffParamRootName 			= function(modelrun = ModelRun)
		{
			if (modelrun$ASVE == "CATEGORICAL") RootName = "Eff_AgeGroup_" else RootName = "Eff_AgeKnot_"	
			return(RootName)
		}
		Convert_Multipliers_ToPosteriors 	= function(BaselineChain, MultiplierChain) ## for MULTIPLICATIVE	SSASVEs
		{
			return(BaselineChain * MultiplierChain)
		}
		Convert_Intercepts_ToPosteriors 	= function(BaselineChain, MultiplierChain) ## for ADDITIVE			SSASVEs
		{
			return(BaselineChain + MultiplierChain)
		}
		SSASVE_Chain 						= function(BaselineSeroStatus = "SeroNeg", serotype = 1, AgeGroupOrAgeKnot = 1, modelrun = ModelRun)
		{
			BSS_Char 			= Choose_BSS_Char_Splines("Age_Efficacy", modelrun, BaselineSeroStatus)
			AgeEffChain 		= CHAINS[, paste0(BSS_Char, ChooseAgeEffParamRootName(modelrun), AgeGroupOrAgeKnot)]
			SerotypeChain 		= CHAINS[, paste0(BSS_Char, "Eff_", serotype)]
			if (ModelRun$SSASVE_Additive)	
				SerotypeEffChain 	= Convert_Intercepts_ToPosteriors	(AgeEffChain, SerotypeChain) else 
				SerotypeEffChain 	= Convert_Multipliers_ToPosteriors	(AgeEffChain, SerotypeChain)
			return(SerotypeEffChain)
		}
		Make95CrI_Char						= function(Mean, Lower, Upper) as.factor(paste0(Mean, " (", Lower, ", ", Upper, ")"))
		Make95CrI_CharVecVerison			= function(VectorOf_Mean_Lower_Upper) as.factor(paste0(VectorOf_Mean_Lower_Upper["Mean"], " (", VectorOf_Mean_Lower_Upper["LowerCRI"], ", ", VectorOf_Mean_Lower_Upper["UpperCrI"], ")"))
		GetChangeIn_RR 						= function(PhaseSeverity = "AM", PrevInf_Con = 0, Mult = FALSE, IncludeMode = FALSE, serotype = 1, SSKs = FALSE)
		{
			if (SSKs) SerotypeString = paste0("_sero", serotype) else SerotypeString = ""
	
			Con_Name = paste0("K", PhaseSeverity, "_", PrevInf_Con		, SerotypeString)
			Vac_Name = paste0("K", PhaseSeverity, "_", PrevInf_Con + 1	, SerotypeString)
			K_Con_Chain = CHAINS[, Con_Name]
			K_Vac_Chain = CHAINS[, Vac_Name]
			
			if (Mult & PhaseSeverity == "PS")
			{
				Con_Name_AM 	= paste0("KAM_", PrevInf_Con		, SerotypeString)
				Vac_Name_AM 	= paste0("KAM_", PrevInf_Con + 1	, SerotypeString)
				K_Con_Chain_AM 	= CHAINS[, Con_Name_AM]
				K_Vac_Chain_AM 	= CHAINS[, Vac_Name_AM]
				
				Ratio_Chain = (K_Vac_Chain * K_Vac_Chain_AM) / (K_Con_Chain * K_Con_Chain_AM)  ## if PhaseSeverity == "PS" then K_Vac_Chain and K_Con_Chain will refer to passive severe chains for vaccine and control groups respectively.
				
			} else	Ratio_Chain = K_Vac_Chain / K_Con_Chain
			
			return(SummStat(Ratio_Chain, IncludeMode = IncludeMode))
		}
		Get_RelRisk_Hosp_Summary 			= function(PrevInf = 0, IncludeMode = FALSE, serotype = 1, SSKs = FALSE) ### PS_Ks_Multiply_AM_Ks then passive Ks are no longer relative risks - they are proportions of symptomatic disease that is hospitalised. So to get relative risk you need to multiply active K by passive K. 
		{
			if (SSKs) SerotypeString = paste0("_sero", serotype) else SerotypeString = ""
			AM_Name = paste0("KAM_", PrevInf, SerotypeString)
			PS_Name = paste0("KPS_", PrevInf, SerotypeString)
			KAM_i_Chain = CHAINS[, AM_Name]
			KPS_i_Chain = CHAINS[, PS_Name]
			return(SummStat(KPS_i_Chain * KAM_i_Chain, IncludeMode = IncludeMode))
		}
		SummaryTable_filename = file.path(ModelPlotSubDir, "RunSummary.txt")
		
		Means			= rep(NA, dim(CHAINS)[2]) ### will include Loglike chain. 
		Medians     	= rep(NA, dim(CHAINS)[2])
		Modes			= rep(NA, dim(CHAINS)[2])
		Lower_CrIs  	= rep(NA, dim(CHAINS)[2])
		Upper_CrIs  	= rep(NA, dim(CHAINS)[2])
		
		SummaryTableDummy			= matrix(nrow = dim(CHAINS)[2], ncol = 5)
		#rownames(SummaryTableDummy) = colnames(CHAINS)
		colnames(SummaryTableDummy) = c("Mean", "Mode", "Median", "LowerCrI", "UpperCrI") 
		ROWNAMES					= colnames(CHAINS)
		
		param_no = 1
		for (param_no in 1:dim(CHAINS)[2])
		{
			ParamChain 			= CHAINS[,param_no]
			SummStat_Dummy		= SummStat(ParamChain, IncludeMode = FALSE)
			
			if (ChainNames[param_no] == "logLikelihood") 
				SummStat_Dummy = round(SummStat_Dummy) 		else 
				SummStat_Dummy = signif(SummStat_Dummy, 2)
			SummaryTableDummy[param_no,] = SummStat_Dummy[colnames(SummaryTableDummy)] ## need colnames so quantities go in right place. 
		}
		
		### get Serotype/age effiacies (or knot params), just add rows to SummaryTableDummy
		if (ModelRun$ASVE != "INDEPENDENT" & ModelRun$SS_VEs == 1)
		{
			if (ModelRun$ASVE == "CATEGORICAL") AgeParamNos = 1:3 else AgeParamNos = 0:3
			for (BS in c("SeroNeg", "SeroPos"))
				for (serotype in 1:N_STypes_VEs) #### make this start from 2. 1 is there just as a check. 
					for (AgeParam in AgeParamNos)
					{
						ParamChain 			= SSASVE_Chain(BaselineSeroStatus = BS, serotype = serotype, AgeGroupOrAgeKnot = AgeParam, modelrun = ModelRun)
						SummStat_Dummy		= SummStat(ParamChain, IncludeMode = FALSE)
						SummStat_Dummy 		= signif(SummStat_Dummy, 2)
						SummaryTableDummy	= rbind (SummaryTableDummy, SummStat_Dummy[colnames(SummaryTableDummy)])
						ROWNAMES			= c(ROWNAMES	, paste0(Choose_BSS_Char_Splines("Age_Efficacy", ModelRun, BS), ChooseAgeEffParamRootName(ModelRun), AgeParam, "_sero", serotype) )
					}
		}
		
		if (exists("WBIC_CHAINS"))
		{
			Cooled_LL_Chain = WBIC_CHAINS[,"logLikelihood"]
			
			NumExtraQuantitiesToInclude = 1
			if (exists("TotalLogLikeChain"		)) NumExtraQuantitiesToInclude = NumExtraQuantitiesToInclude + 1
			if (exists("WBIC_TotalLogLikeChain"	)) NumExtraQuantitiesToInclude = NumExtraQuantitiesToInclude + 1
			for (QQQ in NumExtraQuantitiesToInclude:1) #### doing in reverse order so that I add to top of table. Note this is different to the above. 
			{
				if (QQQ == 1) Chain = WBIC_CHAINS[,"logLikelihood"] else if (QQQ == 2) Chain = TotalLogLikeChain 	else if (QQQ == 3) Chain = WBIC_TotalLogLikeChain
				if (QQQ == 1) RName = "Cooled_logLikelihood" 					else if (QQQ == 2) RName = "NonAug_logLikelihood" 			else if (QQQ == 3) RName =  "Cooled_NonAug_logLikelihood"
				
				SummStat_Dummy		= SummStat(Chain, IncludeMode = FALSE)
				SummStat_Dummy 		= round(SummStat_Dummy) 
				SummaryTableDummy	= rbind (SummStat_Dummy[colnames(SummaryTableDummy)], SummaryTableDummy) #### reverse order of above - summary of likelihood measures go on top of table. 
				ROWNAMES			= c(RName		, ROWNAMES		)
			}
			rm(QQQ, Chain, RName)
		}
		length(ROWNAMES)
		dim(SummaryTableDummy)
		
		##### add in Proportion Hospitalised (do this before Change in Relative risk from Vaccination stuff as the latter should be on the top)
		if (ModelRun$PS_Ks_Multiply_AM_Ks)
		{
			PrevInfs = 2:0
			for (serotype in N_STypes_Ks:1) ### note (reverse) order
				for (PrevInf in PrevInfs)  ### note (reverse) order
					SummaryTableDummy = rbind (		
							signif(Get_RelRisk_Hosp_Summary(PrevInf = PrevInf, serotype = serotype, SSKs = ModelRun$SS_Ks)[colnames(SummaryTableDummy)], 2), 
							SummaryTableDummy)
			if (!ModelRun$SS_Ks) 
				RelRisk_Hosp_RowNames = 				paste0("RelRisk_Hosp_", 0:2, "_PriorInfs") else
				RelRisk_Hosp_RowNames = paste0( rep(	paste0("RelRisk_Hosp_", 0:2, "_PriorInfs") , N_STypes_Ks), "_sero", rep( 1:N_STypes_Ks, each = length(PrevInfs))) ## each
			ROWNAMES		= c(RelRisk_Hosp_RowNames		 , ROWNAMES)
			rm(PrevInf, PrevInfs)
		}
		
		##### add in Change in Relative risk from Vaccination stuff
		PhSevs 		= c("AM", "PS")
		PrevInfs 	= 0:1
		for (serotype in N_STypes_Ks:1) ### note reverse order
			for (PhSev in rev(PhSevs)) ### note reverse order
				for (PrevInf in rev(PrevInfs))
					SummaryTableDummy = rbind (		signif(GetChangeIn_RR(PhaseSeverity = PhSev, PrevInf_Con = PrevInf, Mult = ModelRun$PS_Ks_Multiply_AM_Ks, serotype = serotype, SSKs = ModelRun$SS_Ks)[colnames(SummaryTableDummy)], 2), SummaryTableDummy)
		
		if (!ModelRun$SS_Ks) 
			ChangeInRR_RowNames = 				paste0(	rep(paste0("ChangeInRR_", PhSevs ), each = 2), "_PrevInf_", PrevInfs)	else 
			ChangeInRR_RowNames = paste0( rep(	paste0(	rep(paste0("ChangeInRR_", PhSevs ), each = 2), "_PrevInf_", PrevInfs)		, N_STypes_Ks), "_sero", rep( 1:N_STypes_Ks, each = length(PhSevs) * length(PrevInfs))) 
		
		ROWNAMES		= c(ChangeInRR_RowNames			 , ROWNAMES)
		rm(PhSevs, PrevInfs)
		
		### Make final summary table
		CharacterVec = Make95CrI_Char(SummaryTableDummy[, "Mean"], SummaryTableDummy[, "LowerCrI"], SummaryTableDummy[, "UpperCrI"])  #### this is short hand for tables you'll fill in manually. 
		SummaryTable = cbind(data.frame(SummaryTableDummy), Character = CharacterVec)
		rownames(SummaryTable) = ROWNAMES
		
		### add in DICs
		if (file.exists(DIC_Filename))
		{
			DIC_Value 		= as.numeric(DIC_table["DIC_defn3_MaxLike_LLChain", ])
			BIC_Value 		= CalcBIC()
			WBIC_Value 		= CalcWBIC()
			Values 			= c(DIC_Value, BIC_Value, WBIC_Value)
			LittleDFrame 	= data.frame(Values, NA, NA, NA, NA, as.factor(Values))
			rownames(LittleDFrame) = as.factor(c("DIC", "BIC", "WBIC"))
			colnames(LittleDFrame) = colnames(SummaryTable)
			SummaryTable = rbind(LittleDFrame, SummaryTable)
			rm(LittleDFrame, DIC_Value, BIC_Value, WBIC_Value, Values)
		}
		
		### rearrange so logLikelihood on the top.
		LL_Rownames 	= grep("logLikelihood", rownames(SummaryTable) )
		SummaryTable 	= SummaryTable[c(LL_Rownames, (1:dim(SummaryTable)[1])[-LL_Rownames]), ]
		LL_Rownames 	= which("logLikelihood" == rownames(SummaryTable) )
		SummaryTable 	= SummaryTable[c(LL_Rownames, (1:dim(SummaryTable)[1])[-LL_Rownames]), ]
		
		write.table(SummaryTable, file = SummaryTable_filename,	row.names = T, col.names = NA, quote = F, sep = "\t") ### col.names = NA ensures that colnames aren't moved one to the left (i.e. and are then on top of rownames)
		rm(ParamChain, ROWNAMES, SummaryTable)
	}
		
	##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
	##### ##### ##### Posteriors
	
	PlotPosterior = function(param_no = 2, modelrun = ModelRun, PlotTitle = NULL, FileName = NULL, VecToHist = NULL, Chains = CHAINS, AcceptArray = AcceptanceProbs, PlotPrefix = "", 
			MEAN = TRUE, XLAB = NULL, SavePlot = TRUE, Overwrite = OVERWRITE_FILES, OutputStringInPlotTitle = FALSE, ShortenOutputString = TRUE, RES = PNG_res, XLIM = NULL, CexAxis = 1.2, CexLab = 1.4) ### param_no is 1 out of sync with rest of numbering. Default set to 2 as first column (labelled zero) is log likelihood so that if doing a combined posterior (say K0/K1) then we still round posteriors to significant figures. 
	{
		if (is.null(FileName)) FileName = paste0(PlotPrefix, "P", param_no - 1, "_Post") ## minus 1 due to Cpp notation - Param_0 is log likelihood. 
		FileName = paste0(FileName, ".png")
		
		PNGFILENAME = file.path(PosteriorsPlotDirectory, FileName)
		
		if (is.null(VecToHist)) VecToHist = na.omit(as.numeric(Chains[,param_no]))
		if (SavePlot) png(file = PNGFILENAME, res = RES, units = "in", width = 6, height = 5)
		
		Param_95CrI = quantile(VecToHist,c(0.025,0.975), na.rm = T)
		Lower_CrI 	= Param_95CrI[1]
		Upper_CrI 	= Param_95CrI[2]
		
		if (is.null(PlotTitle)) ParamString = ChainNames[param_no] else ParamString = PlotTitle
		if (MEAN) 				
		{
			MeanOrMedianString 	= "Mean"
			MeanOrMedian		= mean(VecToHist)
		}
		else 
		{
			MeanOrMedianString = "Median"
			MeanOrMedian		= median(VecToHist)
		}
		
		if (ChainNames[param_no] != "logLikelihood")
		{
			Lower_CrI 		= signif(Lower_CrI		, 3)
			Upper_CrI 		= signif(Upper_CrI		, 3)
			MeanOrMedian	= signif(MeanOrMedian	, 3)
		
		} else 
		{
			Lower_CrI 		= round(Lower_CrI)
			Upper_CrI 		= round(Upper_CrI)
			MeanOrMedian 	= round(MeanOrMedian)
		}
		
		if (ShortenOutputString) ModelRun_string = AbbreviateOutputString(OutputSubDir, modelrun = modelrun) else ModelRun_string = OutputSubDir
		
		if (is.null(PlotTitle)) 
		{
			PlotTitle = paste0(PlotPrefix, ParamString, " Posterior\n") 
			if (OutputStringInPlotTitle) PlotTitle = paste0(PlotTitle, ModelRun_string, "\n", )
			PlotTitle = paste0(PlotTitle, MeanOrMedianString, " ", MeanOrMedian, ", 95%CrI ", Lower_CrI," - ", Upper_CrI)
			if (param_no != 1 & (param_no - 1) < length(AcceptArray))  PlotTitle = paste0(PlotTitle, ",   Acpt Pr = ", signif(AcceptArray[param_no - 1], 3))
		}
		if (is.null(XLAB)) XLAB = ChainNames[param_no]
		if (is.null(XLIM)) 	
			hist(	VecToHist,	xlab = XLAB, cex.lab = CexLab, cex.axis = CexAxis, freq = FALSE, col = PosteriorCol,	main = PlotTitle) else				### NULL values of xlim argument in hist not accepted, and default value not easy to change (see hist documentation if you doubt this). Hence just repeat code. 
			hist(	VecToHist,	xlab = XLAB, cex.lab = CexLab, cex.axis = CexAxis, freq = FALSE, col = PosteriorCol,	main = PlotTitle, xlim = XLIM, breaks = 50)
		
		if (SavePlot) dev.off()
	}
	
	if (DoPosteriors)	
	{
		cat ("Posteriors ")
		PosteriorsPlotDirectory = file.path(ModelPlotSubDir, "Posteriors")
		dir.create(PosteriorsPlotDirectory, showWarnings = FALSE)
		
		param_no = 2
		for (param_no in 1:(NoParams + 1)) 
		{
			#### plot different xlims for SNegEffs - so you can see the posterior. 
			### find all SNegEff posteriors, excluding the AgeKnots 
			IsParamSNegEfficacyParam 	= length(grep("SNegEff", colnames(CHAINS)[param_no])) != 0
			IsParamNotAnAgeKnotParam 	= length(grep("AgeKnot", colnames(CHAINS)[param_no])) == 0
			ParamCanBeLessThanZero		= min(CHAINS[, param_no]) < 0
			if (IsParamSNegEfficacyParam & IsParamNotAnAgeKnotParam & ParamCanBeLessThanZero) XLIM = c(-1, 1) else XLIM = NULL
			try(PlotPosterior(param_no = param_no, RES = 200, XLIM = XLIM), silent = TRUE)
			rm (XLIM)
		}
		
		if (exists("TotalLogLikeChain"		) & exists("WBIC_TotalLogLikeChain"	))
		{
			CombinedNonAugLikeFileName = file.path(PosteriorsPlotDirectory, "CombNonAugLike.png")
			png(file = CombinedNonAugLikeFileName, res = 200, units = "in", width = 12, height = 5)
			par(mfrow = c(1,2))
			PlotPosterior(param_no = 1, PlotTitle = "NonAugLike"		, FileName = paste0("NonAugLike_Post")		, VecToHist = TotalLogLikeChain		, SavePlot = FALSE)
			PlotPosterior(param_no = 1, PlotTitle = "WBIC_NonAugLike"	, FileName = paste0("NonAugLike_WBIC_Post")	, VecToHist = WBIC_TotalLogLikeChain, SavePlot = FALSE)
			dev.off()
		}
		CloseOpenPlotDevices()	

		cat ("\t")
	}
	## end DoPosteriors if statement
	
	## close any open plot devices
	CloseOpenPlotDevices()	
		

	##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
	##### ##### ##### HazRatios
	
	if (DoHazRatios) 
	{
		cat ("HazRatios ")
		
		HazRatioPlotDirectory = ModelPlotSubDir
		dir.create(HazRatioPlotDirectory, showWarnings = FALSE)
		
		### Definitions
		ActiveMild 			= 1	; 	PassiveSevere 		= 2	; BothPhasesTogether= 3		## In R PhaseSeverity = 1 refers to ActiveMild, PhaseSeverity = 2 is PassiveSevere
		SeroNeg 			= 1	; 	SeroPos 			= 2	; Both_Sero			= 3		## In R SeroNeg = 1, SeroPos = 2; 
		MonoTypic			= 4	; 	MultiTypic 			= 5;  #### even though you don't augment to non-binary serostatus 
		ControlArm			= 0	; 	VaccineArm			= 1	; EitherArm			= 2
		BaselineSeroStatuses 	= c(Both_Sero, SeroNeg, SeroPos)
		BS_Names 				= c("SeroNeg", "SeroPos", "Any Serostatus", "1 prior infection", "> 1 prior infection")
		
		ChoosePhaseSeverityChar 	= function(modelrun = ModelRun, PhaseSeverity = ActiveMild)
		{
			if (modelrun$disease == "")  
				PhaseSeverity_char = "" else 
				PhaseSeverity_char = Convert_DiseaseString(WhichDiseasePlots[PhaseSeverity], modelrun$disease, modelrun$hosp)

			return(PhaseSeverity_char)
		}
		ChooseSerotypeString		= function(modelrun = ModelRun, serotype = 1)
		{
			if (!modelrun$SS_VEs & !modelrun$SS_Ks) SerotypeString = ""
			else 
			{
				if (serotype == 1) 	
					SerotypeString = "" else 
					SerotypeString = paste0("_sero", serotype - 1) 		### For Hazard Ratios - serotype 1 (0 in Cpp) coded as All serotypes, so Actual serotype = serotype + 1. 
			}	
			return(SerotypeString)
		}
		ChooseHazRatioFileName 		= function(country, modelrun = ModelRun, PhaseSeverity = ActiveMild, 
				Directory = HazRatioPlotDirectory, serotype = 1, withDir = TRUE, Whiskers = TRUE)
		{
			### Build Haz Ratio FileName
			CountryBitOfPNGTitle 				= ChooseCountryBitOfPNGTitle(country)
			PhaseSeverity_char 					= ChoosePhaseSeverityChar	(modelrun = modelrun, PhaseSeverity = PhaseSeverity)
			SerotypeString 						= ChooseSerotypeString		(modelrun = modelrun, serotype = serotype)
			WhichRatiosString 					= paste0(CountryBitOfPNGTitle, PhaseSeverity_char, SerotypeString)
			
			RootName = "HRPs"
			HRPlot_FileName 					= paste0(RootName, WhichRatiosString)
			if (withDir)	HRPlot_FileName		= file.path(Directory, HRPlot_FileName)	
			if (!Whiskers) HRPlot_FileName		= paste0(HRPlot_FileName, "_poly")
			HRPlot_FileName 					= paste0(HRPlot_FileName, ".png") ### do this here otherwise file.exists etc won't work
			return(HRPlot_FileName)
		}

		CloseOpenPlotDevices()	## end DoParamCorrelations if statement
		
		ImportHazRatPosts()
		ls()[grep("HazRat", ls())]
		
		GetHazRatioValue = function(StatString = "Mean", AgeGroups = 0, Countries = 12, SeroStatuses = Both_Sero, Days = 0, ImSub = FALSE, HR_serotype = 1, DiseaseSeverity = "")
		{
			ImSubString 		= GetImSubString(ImSub)
			TableVariableName 	= GetHazRatioTableName(StatString = StatString, HR_serotype = HR_serotype, DiseaseSeverity = DiseaseSeverity, ImSubString)
			Values 				= get(TableVariableName)[paste0("A_", AgeGroups, "_C_", Countries, "_BS_", SeroStatuses - 1), paste0("D_", Days) ]
			return(Values)
		}	
		
		PlotHazRatioPosts = function(Country = 10, AgeGroups = 0:3, SeroStatus = Both_Sero, HR_serotype = 1, DiseaseSeverity = "", Days = c(0, 365, 760, 1095), 
				LowerYLim = 0, UpperYlim = NULL, MaxYLim = 12, PLOTTITLE = "", Whiskers = FALSE, AddTimeArrows = FALSE, YCoord_TimeArrows = -0.05,
				CEXAXIS = 3, LWD = 3, AB_LWD = 7)
		{
			### note: providing value for UpperYlim will override MaxYLim
			
			if (Whiskers) AgeGroupOffset 	= 0.3 else 
			{
				#### if each age group plot centered around AgeGroupIndex (i.e. 1,2,3...) then want this value to be less than 0.5 (otherwise haz rat progressions will overlap)
				if (!AddTimeArrows)		AgeGroupOffset = 0.45 else AgeGroupOffset = 0.3
			}
			
			MeanValues 		= GetHazRatioValue(StatString = "Mean"		, AgeGroups = AgeGroups, Countries = Country, SeroStatuses = SeroStatus, Days = Days, HR_serotype = HR_serotype, DiseaseSeverity = DiseaseSeverity)
			LowerCIValues 	= GetHazRatioValue(StatString = "LowerCI"	, AgeGroups = AgeGroups, Countries = Country, SeroStatuses = SeroStatus, Days = Days, HR_serotype = HR_serotype, DiseaseSeverity = DiseaseSeverity)
			UpperCIValues 	= GetHazRatioValue(StatString = "UpperCI"	, AgeGroups = AgeGroups, Countries = Country, SeroStatuses = SeroStatus, Days = Days, HR_serotype = HR_serotype, DiseaseSeverity = DiseaseSeverity)
			
			if (is.null(UpperYlim)) UpperYlim = min(MaxYLim, max(c(1, max(rbind(MeanValues, LowerCIValues, UpperCIValues))))) ### i.e. no higher than MaxYLim, but at least as high as 1. 
			plot(NA, ylim = c(LowerYLim, UpperYlim), 
					xlim = c(1 - AgeGroupOffset, length(AgeGroups) + AgeGroupOffset), xaxt = "n", xlab = "", 
					ylab = "", main = PLOTTITLE, cex.axis = CEXAXIS)
			abline(h = 1, col = "grey", lwd = AB_LWD)
			AgeGroupIndex = 1
			abline(h = 0, col = "black", lwd = 1)
			
			for (AgeGroupIndex in 1:length(AgeGroups))
			{
				HR_xvalues 		= seq(AgeGroupIndex - AgeGroupOffset, AgeGroupIndex + AgeGroupOffset,  length.out = length(Days))
				
				if (Whiskers)
				{
					Arrows(HR_xvalues, as.numeric(LowerCIValues[AgeGroupIndex,]), HR_xvalues, as.numeric(UpperCIValues[AgeGroupIndex,]), code = 3 , arr.type = "T" , col = "red", lwd = LWD)
					points(HR_xvalues, MeanValues[AgeGroupIndex,], col = "red", lwd = LWD, pch = 4, cex = 2)
					
				} else PlotPolygon(XValues = HR_xvalues, LowerLine = LowerCIValues[AgeGroupIndex,], UpperLine = UpperCIValues[AgeGroupIndex,], MeanLine = MeanValues[AgeGroupIndex,])
				
				if (AddTimeArrows)
				{
					text(HR_xvalues,  rep(YCoord_TimeArrows, length(HR_xvalues)), 
							labels = paste0("Y", 1:length(Days) - 1), pos = 1, cex = 1.6)
				}
				
			}
			axis(side = 1, labels= rep("", length(AgeGroups)), at = 1:length(AgeGroups))
			mtext(side = 1, text = AgeGroupNames[AgeGroups + 1], at = 1:length(AgeGroups), cex = CEXAXIS, line = 2)
		}
		
		PlotMultiHazRatios = function(HR_serotype = 1, DiseaseSeverity = "", WholePlotTitle = NULL, AddWholePlotTitle = TRUE, 
				Whiskers = FALSE, SavePlot = TRUE, Overwrite = OVERWRITE_FILES, ShortenOutputString = TRUE, OutputStringInTitle = TRUE,  
				modelrun = ModelRun, CEXAXIS = 3, AddTimeArrows = FALSE, AddMetaLeftAxes = TRUE,
				RowNames = c("Either Serostatus", "seronegative", "seropositive"), ### can be any one or combo of these, but nothing else. (or can have \n before each)
				ColNames = c("CYD14", "CYD15"), 
				SwitchRowAndColNames = FALSE, PNG_Height = NULL, PNG_Width = NULL, LowerYLim = 0, YCoord_TimeArrows = -0.05,
				Ratio_Plots_To_Title = 4, Ratio_Plots_To_ColNames = 8, Ratio_Plots_To_RowNames = 4, 
				ColName_Cex = 3, RowName_Cex = 3, MultiPlot_Title_Cex = 3)
		{
			SerotypeString = ChooseSerotypeString(serotype = HR_serotype)
			
			#### choose PhaseSeverity (only needed for ChooseHazRatioFileName function and WholePlotTitle (really hacky and probably quite circular but otherwise requires a rewrite)
			if (DiseaseSeverity == ""			) PhaseSeverity = ActiveMild 			else 	## arbitrary. Value ignored for disease == "" 
			if (DiseaseSeverity == "_Mild"		) PhaseSeverity = ActiveMild 			else
			if (DiseaseSeverity == "_Severe"	) PhaseSeverity = PassiveSevere 		else
			if (DiseaseSeverity == "_Either"	) PhaseSeverity = BothPhasesTogether 			### remember WhichDiseasePlots = c("_Mild", "_Severe", "_Either")
			
			### ChooseVacBenefitCountries
			### Apart from Both_Sero (statuses), hazard ratios will vary between countries if either: phase is passive (or plotting both together) & you are not fixing relative risks and you are modelling different hosp K values. 
			### Only variation will be in Both_Sero, which reflects only small differences in K+ country posteriors and is not worht it's own plot. 
			### In that one scenario, plot for CYD14 and CYD15 separately (i.e. countries 10 and 11). For everything else, plot combined CYD14 and CYD15 (i.e. country = 12)
			
			
			HRPost_FileName	= ChooseHazRatioFileName(country = 12, modelrun = ModelRun, PhaseSeverity = PhaseSeverity, Directory = HazRatioPlotDirectory, serotype = HR_serotype, Whiskers = Whiskers)
			cat (paste0 ("\n", "HRPosts ST ", HR_serotype, " PhaseSeverity ", PhaseSeverity, " "))
			
			if (is.null(PNG_Height)) 	
			{
				PNG_Height = (length(RowNames) + (1 / Ratio_Plots_To_ColNames)) * 4
				if (AddWholePlotTitle) PNG_Height = PNG_Height + 2
			}
			if (is.null(PNG_Width)) 	PNG_Width = 14
			
			if (SwitchRowAndColNames)
			{
				## store
				HDummy = PNG_Height
				WDummy = PNG_Width
				CNamesDummy = ColNames
				RNamesDummy = RowNames
				
				# reset and switch
				PNG_Height 	= WDummy
				PNG_Width	= HDummy
				ColNames =  RNamesDummy
				RowNames =  CNamesDummy
				
				BY_ROW = FALSE
				
			}	else BY_ROW = TRUE
			
			
			png(file = HRPost_FileName, res = PNG_res, units = "in"	, width = PNG_Width, height = PNG_Height) 
			#par(mfrow = c(length(BaselineSeroStatuses), length(Days)), omi = c(0.4, 0.5, 1.3, 0)) #(1=bottom, 2=left, 3=top, 4=right).
			
			par(mar = OrigMAR)
			if (is.null(WholePlotTitle) & AddWholePlotTitle)
			{
				WholePlotTitle = paste0(ChoosePhaseSeverityChar(PhaseSeverity = PhaseSeverity), SerotypeString)
				if (OutputStringInTitle) WholePlotTitle = paste0(WholePlotTitle,  "\n", AbbreviateOutputString(OutputSubDir))
				while (substr(WholePlotTitle, 1, 1) == "_") WholePlotTitle = sub("_", "", WholePlotTitle)
			}
			SetUpMultiPlot (#NoRows = length(BaselineSeroStatuses), NoCols = length(Days), 						
					ColNames = ColNames, 
					RowNames = RowNames, 
					MultiPlot_Title = WholePlotTitle, AddTitle = AddWholePlotTitle, BY_ROW = BY_ROW, 
					Ratio_Plots_To_Title = Ratio_Plots_To_Title, Ratio_Plots_To_ColNames = Ratio_Plots_To_ColNames, 
					Ratio_Plots_To_RowNames = Ratio_Plots_To_RowNames, 
					ColName_Cex = ColName_Cex, RowName_Cex = RowName_Cex, MultiPlot_Title_Cex = MultiPlot_Title_Cex, OG_Mar = OrigMAR, modelrun = ModelRun
			)
			par (mar = c(5.1, 3.1, 1.1, 1.1)) #(1=bottom, 2=left, 3=top, 4=right).
			
			if (SwitchRowAndColNames) NamesToCheck = ColNames else NamesToCheck = RowNames 
			
			if (length(grep("Either Serostatus", NamesToCheck)) != 0)
			{
				PlotHazRatioPosts(Country = 10, AgeGroups = 0:3		, SeroStatus = Both_Sero, HR_serotype = HR_serotype, Whiskers = Whiskers, DiseaseSeverity = DiseaseSeverity, UpperYlim = 1.5	, LowerYLim = LowerYLim, CEXAXIS = CEXAXIS, AddTimeArrows = AddTimeArrows, YCoord_TimeArrows = YCoord_TimeArrows)
				PlotHazRatioPosts(Country = 11, AgeGroups = c(0,4,5), SeroStatus = Both_Sero, HR_serotype = HR_serotype, Whiskers = Whiskers, DiseaseSeverity = DiseaseSeverity, UpperYlim = 1.5	, LowerYLim = LowerYLim, CEXAXIS = CEXAXIS, AddTimeArrows = AddTimeArrows, YCoord_TimeArrows = YCoord_TimeArrows)
			}
			if (length(grep("seronegative", NamesToCheck)) != 0)
			{
				PlotHazRatioPosts(Country = 10, AgeGroups = 0:3		, SeroStatus = SeroNeg	, HR_serotype = HR_serotype, Whiskers = Whiskers, DiseaseSeverity = DiseaseSeverity, UpperYlim = 8.5	, LowerYLim = LowerYLim, CEXAXIS = CEXAXIS, AddTimeArrows = AddTimeArrows, YCoord_TimeArrows = YCoord_TimeArrows)
				PlotHazRatioPosts(Country = 11, AgeGroups = c(0,4,5), SeroStatus = SeroNeg	, HR_serotype = HR_serotype, Whiskers = Whiskers, DiseaseSeverity = DiseaseSeverity, UpperYlim = 8.5	, LowerYLim = LowerYLim, CEXAXIS = CEXAXIS, AddTimeArrows = AddTimeArrows, YCoord_TimeArrows = YCoord_TimeArrows)
			}
			if (length(grep("seropositive", NamesToCheck)) != 0)
			{
				PlotHazRatioPosts(Country = 10, AgeGroups = 0:3		, SeroStatus = SeroPos	, HR_serotype = HR_serotype, Whiskers = Whiskers, DiseaseSeverity = DiseaseSeverity, UpperYlim = 1.15	, LowerYLim = LowerYLim, CEXAXIS = CEXAXIS, AddTimeArrows = AddTimeArrows, YCoord_TimeArrows = YCoord_TimeArrows)
				PlotHazRatioPosts(Country = 11, AgeGroups = c(0,4,5), SeroStatus = SeroPos	, HR_serotype = HR_serotype, Whiskers = Whiskers, DiseaseSeverity = DiseaseSeverity, UpperYlim = 1.15	, LowerYLim = LowerYLim, CEXAXIS = CEXAXIS, AddTimeArrows = AddTimeArrows, YCoord_TimeArrows = YCoord_TimeArrows)
				
			}
			par (mar = OrigMAR) #(1=bottom, 2=left, 3=top, 4=right).
			par (mfrow = c(1, 1), omi = rep(0, 4)) #(1=bottom, 2=left, 3=top, 4=right).
			if (AddMetaLeftAxes)
				mtext( "Hazard Ratio: Vaccine / Control"	, side = 2, line = 2.2,	cex = 2.4, adj = 0.28) 	#(1=bottom, 2=left, 3=top, 4=right).
			dev.off()
		}
		
		for (HR_serotype in 1:HRs_NumSTypes)
			for (DiseaseSeverity in WhichDiseasePlots)
				PlotMultiHazRatios(HR_serotype = HR_serotype, DiseaseSeverity = DiseaseSeverity, Whiskers = F,
						Ratio_Plots_To_RowNames = 4, ColName_Cex = 4, RowName_Cex = 2.1, CEXAXIS = 2.5, AddWholePlotTitle = FALSE, LowerYLim = -0.3, YCoord_TimeArrows = 0, AddTimeArrows = T) ### polygones (rather than whiskers) are easier to read - but the lines imply continunity when values actually discrete. 
	} 	## end DoHazRatios if statement
		
	CloseOpenPlotDevices()	
		
	##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
	##### ##### ##### Chains
	
	if (DoChains)	
	{
		cat ("Chains ")
		
		ChainsPlotDirectory = file.path(ModelPlotSubDir, "Chains")
		dir.create(ChainsPlotDirectory, showWarnings = FALSE)
		
		PlotChains = function(Chains = NULL, param_no = 2, modelrun = ModelRun, PlotTitle = NULL, FileName = NULL, ChainToPlot = NULL, Dir = ChainsPlotDirectory, ChainPlotTitle = NULL, 
				EveryN_iter = 100, AddMeanLine = TRUE, MeanLineCol = "red", MeanLine_LWD = 2.5, YLIM = NULL, YLAB = NULL, 
				MEAN = TRUE, XLAB = NULL, SavePlot = TRUE, OutputStringInPlotTitle = FALSE, ShortenOutputString = TRUE, RES = 160, XLIM = NULL, CexAxis = 1.2, CexLab = 1.4, WBICs = FALSE, AcptProb = NULL)
		{
			WBIC_String = Choose_WBIC_String(WBICs)
			
			### Choose defaults. 
			if (is.null(Chains))
			{
				if (WBICs) Chains = WBIC_CHAINS else Chains = CHAINS
			}
			if (is.null(FileName)) 		FileName = file.path(Dir, paste0(WBIC_String, "P", param_no - 1, "_Chain.png"))
			if (is.null(ChainToPlot)) 	ChainToPlot = na.omit(as.numeric(Chains[,param_no]))
			if (SavePlot) png(file = FileName, res = RES, units = "in", width = 6, height = 5)
			if (is.null(AcptProb))
			{
				if (WBICs) AcptProb = WBIC_AcceptanceProbs else  AcptProb = AcceptanceProbs
			}
			if (is.null(YLAB)) YLAB = paste0(WBIC_String, ChainNames[param_no])
			
			if (is.null(ChainPlotTitle))
			{
				ChainPlotTitle = paste0(WBIC_String, ChainNames[param_no]," chain")
				if (param_no != 1 && param_no != (NoParams + 1) && length(grep("rho", ChainNames[param_no])) == 0)   ChainPlotTitle = paste0(ChainPlotTitle, ",   Acpt Pr = ", signif(AcptProb[param_no - 1], 3)) ### second condition is because you don't have acceptance probabities for serotype rho parameters, but you do output their values in the chains.
				if (ShortenOutputString) OutputStringDummy  = AbbreviateOutputString(OutputSubDir, modelrun)  else OutputStringDummy = OutputSubDir
				if (OutputStringInPlotTitle) ChainPlotTitle = paste0(ChainPlotTitle, "\n", OutputStringDummy)
			}
			plot(1:length(ChainToPlot) * EveryN_iter, ChainToPlot, type = "l", ylim = YLIM, 
					ylab = YLAB, xlab = "Iteration", 
					cex.lab = CexLab, cex.axis = CexAxis, main = ChainPlotTitle)
			if (AddMeanLine) abline(h = mean(ChainToPlot), col = MeanLineCol, lwd = MeanLine_LWD)
			if (SavePlot) dev.off()
		}

		if (DoChains) for (param_no in 1:(NoParams + 1)) PlotChains(param_no = param_no)
		cat ("\t")
		
		if (ModelRun$SS_VEs || ModelRun$SS_Ks)
		{
			COMBINED_RHO_PLOT = FALSE
			#COMBINED_RHO_PLOT = TRUE
			ParamNames 		= ChainNames[-1]
			Param_PCHs	 	= rep(1, NoParams)
			ParamColours 	= rep("black", NoParams)
			
			SerotypeCols 	= c("red", "grey", "blue", "orange")
			
			country = 0
			Xaxis = seq(0, dim(CHAINS)[1]*100, length.out = 11)
			if (COMBINED_RHO_PLOT)
			{
				png(file = file.path(ChainsPlotDirectory, "CombRhoStacks.png"), units = "in", width = 12, height = 30, res = PNG_res)
				par(mfcol = c(5,2))
			}
			for (country in 0:9)
			{
				StackedRho_ChainFileName = file.path(ChainsPlotDirectory, paste0("c", country, "_RhoStack.png"))
				if (!COMBINED_RHO_PLOT) png(file = StackedRho_ChainFileName, units = "in", width = 6, height = 6, res = PNG_res)
				barplot(t(CHAINS[, paste("rho", country, 1:4, sep = "_" )]), col = SerotypeCols, border = NA, space = 0, cex.axis = 1.2, cex.lab = 1.4, 
						ylab = "", xlab = "", main = CountryNames_Long[country + 1], cex.main = 2.5)
				mtext(side = 2, text = "Proportion (stacked)"	, cex = 1.4, line = 2.7)
				mtext(side = 1, text = "Iteration"				, cex = 1.4, line = 1)
				StackedRho = rep(0, dim(CHAINS)[1])
				serotype = 1
				for (serotype in 1:4)
				{
					PreviousRho = StackedRho
					RhoChain 	= CHAINS[, paste("rho", country, serotype, sep = "_" )]
					StackedRho 	= StackedRho + RhoChain
					text(x = 0.5 * dim(CHAINS)[1], y = mean((PreviousRho + StackedRho)/2), adj = 0.5,
							labels = paste0("serotype ", serotype), col = "white", cex = 1.3)
				}
				if (!COMBINED_RHO_PLOT) dev.off()
			}
			if (COMBINED_RHO_PLOT) dev.off()
		}
	} 
	 	
	CloseOpenPlotDevices()	## end DoChains if statement
	
	##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
	##### ##### ##### SurvivalCurves
	
	if (DoSurvCurves)	
	{
		SCurveDirectory = file.path(ModelPlotSubDir, "SurvCurvs")
		dir.create(SCurveDirectory, showWarnings = FALSE)
		
		if (DoSurvCurves)
		{
			cat ("Survival Curves: ")
			
			### want to plot: 
			# a) For given age group, all countries in that age group
			# 		i) 		All ages 		AgeGroup 0 ==> Countries = countriesfitted - plus CYD14 combined, CYD15 combined and both trials combined
			# 		ii) 	CYD14 2-5: 		AgeGroup 1 ==> Countries = c(0:4, 10)
			# 		iii)	CYD14 6-11: 	AgeGroup 2 ==> Countries = c(0:4, 10)
			# 		iv)		CYD14 12-14: 	AgeGroup 3 ==> Countries = c(0:4, 10)
			# 		v)		CYD15 9-11: 	AgeGroup 4 ==> Countries = c(5:9, 11)
			# 		vi)		CYD15 12-16: 	AgeGroup 5 ==> Countries = c(5:9, 11)
			# b) For given country, plot all age groups that apply to that country
			# 		i) 		if country <=4 or country == 10 (CYD15) AgeGroups = 0:3
			# 		ii) 	if country >=5 or country == 11 (CYD15) AgeGroups = c(0,4,5)
			
			PlotMultipleSurvivalCurves = function(Countries, AgeGroups = 0, Phases = "", DiseaseSeverities = "", ImSubs = FALSE, CurvesToPlot = "ByTrialArm", 
					DATA = DATA, countriesfitted = CountriesFitted, ObsMethod = "Default", Passive = FALSE, SavePlot = TRUE, Overwrite = OVERWRITE_FILES, 
					ShortenOutputString = TRUE, AgeMin = NULL, AgeMax = NULL, modelrun = ModelRun, 
					Either_FullDuration = TRUE, ExpectedCasesInPlotTitles = FALSE, 
					Directory = SCurveDirectory, SubDir = NULL, FileName = NULL, Plot_ModalPost = FALSE, Plot_MaxLike = FALSE, WBIC_String = "",
					AddMultiPlotTitle = FALSE, MultiPlot_Title = NULL, MultiPlot_Title_Cex = 2, ColName_Cex = 1.4, RowName_Cex = ColName_Cex, 
					RowNames = NULL, ColNames = NULL, 
					PlotToTitleHeightRatio = 7, PlotToColNameHeightRatio = PlotToTitleHeightRatio/2, PlotToRowNameWidthRatio = PlotToColNameHeightRatio,
					N_Rows = NULL, N_Cols = NULL, IncABlines = FALSE, Ind_PlotHeadings = NULL, Inc_Ind_AxesLabels = TRUE, IncDataLegend = TRUE,
					CEXLAB = 1.4, CEXAXIS = 1.2, LOWERYLIM = NULL, IndividualTitleCex = 2, Inc_IndividualPlotTitles = TRUE, LEGCEX = 1.3)
			{
				CloseOpenPlotDevices()
				
				GroupLengths = c(length(Countries), length(AgeGroups))	
				NoPlotsTotal = prod(GroupLengths) 						 
				
				### Build FILENAME
				AGs_String = paste0("_AGs", paste0(AgeGroups, collapse = ""))
				if (all(Countries == 0:9)		) Countries_String 	= "_All" 										else 
				if (all(Countries == 10)		) Countries_String 	= "_CYD14" 										else 
				if (all(Countries == c(10, 0:4))) Countries_String 	= "_CYD14countries" 							else 
				if (all(Countries == 11)		) Countries_String 	= "_CYD15" 										else  
				if (all(Countries == c(11, 5:9)	) | all(Countries == 11)) Countries_String 	= "_CYD15countries" 	else  
					Countries_String = paste0("_Cs", paste0(Countries, collapse = ""))
				
				DiseaseSeverities_String 	= paste0(DiseaseSeverities	, collapse = "")
				Phases_String 				= paste0(Phases				, collapse = "")
				
				if (!is.null(SubDir))  Directory = file.path(Directory, SubDir)
				dir.create(Directory, showWarnings = FALSE)
				
				if (Plot_ModalPost) ModPostString = "_MP" else ModPostString = ""
				if (Plot_MaxLike) 	MaxLikeString = "_ML" else MaxLikeString = ""
				FileName_NoDir 	= paste0(WBIC_String, "SCs", AGs_String, Countries_String, DiseaseSeverities_String, Phases_String, ModPostString, MaxLikeString, ".png")
				FileName 		= file.path(Directory, FileName_NoDir)
				
				if (is.null(N_Rows)) N_Rows = length(AgeGroups) ## put in defaults if need be.
				if (is.null(N_Cols)) N_Cols = length(Countries)
				
				OrigSurvivalCurveWidth = 5; OrigSurvivalCurveHeight = 5; 
				if (SavePlot) png(file = FileName  , res = PNG_res, units = "in", width = OrigSurvivalCurveWidth * N_Cols, height = OrigSurvivalCurveHeight * N_Rows + (OrigSurvivalCurveHeight/PlotToTitleHeightRatio))
				
				if (is.null(RowNames))  AddRowNames = FALSE else AddRowNames = TRUE 
				if (is.null(ColNames))  AddColNames = FALSE else AddColNames = TRUE 
				
				if (is.null(MultiPlot_Title) && AddMultiPlotTitle)
				{
					if (ShortenOutputString) MultiPlot_Title = AbbreviateOutputString(OutputSubDir, modelrun) else MultiPlot_Title = OutputSubDir
					MultiPlot_Title = paste0(MultiPlot_Title, "\n", AGs_String, Countries_String, DiseaseSeverities_String, Phases_String)
				}
				
				SetUpMultiPlot(NoCols = N_Cols, NoRows = N_Rows, AddRowNames = AddRowNames, AddColNames = AddColNames, AddTitle = AddMultiPlotTitle,
						ColNames = ColNames, RowNames = RowNames, MultiPlot_Title = MultiPlot_Title, BY_ROW = TRUE, RemoveCorner = TRUE, 
						Ratio_Plots_To_Title = PlotToTitleHeightRatio, Ratio_Plots_To_ColNames = PlotToColNameHeightRatio, Ratio_Plots_To_RowNames = PlotToRowNameWidthRatio, 
						ColName_Cex = ColName_Cex, RowName_Cex = RowName_Cex, MultiPlot_Title_Cex = MultiPlot_Title_Cex)
				
				par(mar = OrigMAR) ### reset default margins
				if (AddMultiPlotTitle) OutputStringInTitle = FALSE else OutputStringInTitle = TRUE
				plotcounter = 1
				for (AgeGroup in AgeGroups)
					for (country in Countries)
						for (DiseaseSeverity in DiseaseSeverities)
						{
							if (is.null(Ind_PlotHeadings)) PLOTTITLE = NULL else PLOTTITLE = Ind_PlotHeadings[plotcounter]
							
							PlotSurvivalCurves(country = country, AgeGroup = AgeGroup, SavePlot = FALSE, 
									Plot_ModalPost = Plot_ModalPost, Plot_MaxLike = Plot_MaxLike, WBIC_String = WBIC_String, PLOTTITLE = PLOTTITLE,
									WhichDiseasePlot = DiseaseSeverity, OutputStringInTitle = OutputStringInTitle, Either_FullDuration = Either_FullDuration, 
									ExpectedCasesInPlotTitle = ExpectedCasesInPlotTitles, TitleCex = IndividualTitleCex, IncludePlotTitle = Inc_IndividualPlotTitles, IncDataLegend = IncDataLegend,
									IncABlines = IncABlines, CEXLAB = CEXLAB, CEXAXIS = CEXAXIS, LOWERYLIM = LOWERYLIM, 
									IncAxesLabels = Inc_Ind_AxesLabels, LEGCEX = LEGCEX)
							plotcounter = plotcounter + 1
						}
				if (SavePlot) dev.off()
			}

			if (WBIC_SUMMARY_SURVIVAL_TABLES_EXIST) WBICString = "WBIC_" else WBICString = ""
			
			### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === 
			### CountryComps: Countries plotted side by side. Disease and AgeGroup constant
			cat  ("Countries... ")
			Countries = c(10, 0:4)
			for (AGs in 0:3)	
			{
				if (ModelRun$disease == "") DiseasePlots = "" else if (AGs == 0) DiseasePlots = WhichDiseasePlots else DiseasePlots = "_Either"  ## i.e. don't plot the country breakdowns for say disease == "_Severe" for 12-16 year olds. Too much detail. 
				for (Disease in DiseasePlots)
					PlotMultipleSurvivalCurves(Countries = Countries, AgeGroups = AGs, N_Rows = 2, N_Cols = 3,
							Ind_PlotHeadings = paste0(CountryNames_Long[Countries + 1], ": ", AgeGroupNamesVerbose[AGs + 1]), 
							DiseaseSeverities = Disease, SubDir = "CountryComps") ## for each relevant age group, do all countries in CYD14
			}
			Countries = c(11, 5:9)
			for (AGs in c(0,4,5))	
			{
				if (ModelRun$disease == "") DiseasePlots = "" else if (AGs == 0) DiseasePlots = WhichDiseasePlots else DiseasePlots = "_Either"  ## i.e. don't plot the country breakdowns for say disease == "_Severe" for 12-16 year olds. Too much detail. 
				for (Disease in DiseasePlots)
					PlotMultipleSurvivalCurves(Countries = Countries, AgeGroups = AGs, N_Rows = 2, N_Cols = 3, 
							Ind_PlotHeadings = paste0(CountryNames_Long[Countries + 1], ": ", AgeGroupNamesVerbose[AGs + 1]), 
							DiseaseSeverities = Disease, SubDir = "CountryComps") ## for each relevant age group, do all countries in CYD15
			}
			rm(Countries)
			
			#### Plot all countries, for all "disease"
			if (ModelRun$disease == "") DiseasePlots = "" else DiseasePlots = WhichDiseasePlots
			for (Disease in DiseasePlots)
				PlotMultipleSurvivalCurves(Countries = 0:9, AgeGroups = 0, N_Rows = 2, N_Cols = 5, DiseaseSeverities = Disease, 
						SubDir = "CountryComps",
						IncABlines = FALSE, Ind_PlotHeadings = CountryNames_Long[1:10], 
						CEXAXIS = 2.3, LOWERYLIM = 0.85, IndividualTitleCex = 4, Inc_IndividualPlotTitles = TRUE, LEGCEX = 2.5) ## for each relevant age group, do all countries in CYD15
			rm(DiseasePlots) 
			
			### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === 
			### AgeComps: AgeGroups plotted side by side. Country and Disease constant
			cat  ("Age groups... ")
			
			for (Country in c(10)) 
			{
				if (ModelRun$disease == "") DiseasePlots = "" else if (Country == 10) DiseasePlots = WhichDiseasePlots else DiseasePlots = "_Either"  ## i.e. don't plot the age breakdowns for say disease == "_Severe" for country 0. Too much detail. 
				for (Disease in DiseasePlots)
					PlotMultipleSurvivalCurves(
							Countries = Country, AgeGroups = 0:3, N_Rows = 2, N_Cols = 2, 
							IncABlines = FALSE, Ind_PlotHeadings = paste0("CYD-14: ", c("All Ages", "2-5 years", "6-11 years", "12-14 years")), #Inc_Ind_AxesLabels = FALSE,
							AddMultiPlotTitle = FALSE, CEXAXIS = 1.6, LOWERYLIM = 0.85, 
							IndividualTitleCex = 2.7, Inc_IndividualPlotTitles = TRUE, LEGCEX = 1.8, IncDataLegend = TRUE,
							DiseaseSeverities = Disease, SubDir = "AgeComps")  ## for each CYD14 country, do all relevant age groups
			}
			for (Country in c(11))
			{
				if (ModelRun$disease == "") DiseasePlots = "" else if (Country == 11) DiseasePlots = WhichDiseasePlots else DiseasePlots = "_Either"  ## i.e. don't plot the age breakdowns for say disease == "_Severe" for country 5. Too much detail. 
				for (Disease in DiseasePlots)
					PlotMultipleSurvivalCurves(Countries = Country, AgeGroups = c(0,4,5), N_Rows = 2, N_Cols = 2, DiseaseSeverities = Disease, 
							IncABlines = FALSE, Ind_PlotHeadings = paste0("CYD-15: ", c("All Ages", "9-11 years", "12-16 years")), Inc_Ind_AxesLabels = FALSE,
							AddMultiPlotTitle = FALSE, CEXAXIS = 1.6, LOWERYLIM = 0.928, 
							IndividualTitleCex = 2.7, Inc_IndividualPlotTitles = TRUE, LEGCEX = 1.8, IncDataLegend = TRUE,
							SubDir = "AgeComps")  		 ## for each CYD15 country, do all relevant age groups
			}
			
			
			### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === 
			### DiseaseComps: Diseases plotted side by side. Country and AgeGroup constant
			cat  ("Disease severities ")
			if (ModelRun$disease == "_MILDSEVERE")
			{
				if (ModelRun$hosp == "_hosp") 
					Headings = c("any disease", "non-hospitalised", "hospitalised") else
					Headings = c("any disease", "non-severe", "severe")
				for (Cs in c(10, 11)) ### only do this for whole trials, all ages. 
					PlotMultipleSurvivalCurves(Countries = Cs, AgeGroups = 0, N_Rows = 1, N_Cols = 3, DiseaseSeverities = c("_Either", "_Mild", "_Severe"), 
							Ind_PlotHeadings = Headings, IncDataLegend = TRUE, AddMultiPlotTitle = FALSE, SubDir = "DiseaseComps")
			}
			
			
			
			### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === 
			### Full country and age comparisons
			for (DiseaseSeverities in WhichDiseasePlots)
			{
				#### CYD-14
				AgeGroups 	= 0:3			;	Row_Names 	= paste0(paste("AgeGroup", AgeGroups), "\n", AgeGroupNames[AgeGroups + 1]); 
				Countries 	= c(10, 0:4)	;	Col_Names 	= CountryNames_Long[Countries + 1]
				PlotMultipleSurvivalCurves(Countries = Countries, AgeGroups = AgeGroups, Inc_IndividualPlotTitles = FALSE, DiseaseSeverities = DiseaseSeverities,
						PlotToTitleHeightRatio = 7, PlotToColNameHeightRatio = 3.5, MultiPlot_Title_Cex = 4, ColName_Cex = 3,
						RowNames = Row_Names, ColNames = Col_Names)
				rm(AgeGroups, Countries, Row_Names, Col_Names)
				
				#### CYD-15
				AgeGroups 	= c(0,4,5)		;	Row_Names	= paste0(paste("AgeGroup", AgeGroups), "\n", AgeGroupNames[AgeGroups + 1]); 
				Countries 	= c(11, 5:9)	;	Col_Names 	= CountryNames_Long[Countries + 1]
				PlotMultipleSurvivalCurves(Countries = Countries, AgeGroups = AgeGroups, Inc_IndividualPlotTitles = FALSE, DiseaseSeverities = DiseaseSeverities, 
						PlotToTitleHeightRatio = 7, PlotToColNameHeightRatio = 3.5, MultiPlot_Title_Cex = 4, ColName_Cex = 3,	RowNames = Row_Names, ColNames = Col_Names)
				rm(AgeGroups, Countries, Row_Names, Col_Names)
			}
			rm(DiseaseSeverities, WBICString)
			#cat (paste0("CombinedSurvivalCurves not coded for MILDSEVERE\nCombinedSurvivalCurves not coded for MILDSEVERE\n"))
		}
		
		CloseOpenPlotDevices()	## end DoSurvivalCurves if statement
		
		cat (" DONE\n")
	}
	
	CloseOpenPlotDevices()	
	
	##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
	##### ##### ##### Seroprevalence curves

	if (DoSeroPrevs)	
	{
		CloseOpenPlotDevices()	
		
		cat ("SeroPrevs ")
		
		SeroPrevsPlotDirectory = ModelPlotSubDir
		dir.create(SeroPrevsPlotDirectory, showWarnings = FALSE)
		
		PlotSeroPrevs = function (country, Plot_All_Patients = TRUE, Plot_Augmented_Patients = TRUE, Plot_HH_Posteriors = FALSE, SavePlot = TRUE, Overwrite = OVERWRITE_FILES, 
				ShortenOutputString = TRUE, modelrun = ModelRun, OutputStringInTitle = TRUE, TitleCex = 1.2, LegCex = 1.1, IncLegend = TRUE, 
				PlotPrefix = "Predicted Seroprevalence ", FileNamePrefix = "", 
				Mean_Output = Mean_SeroPrevOutput, UpperCrI_Output = UpperCI_SeroPrevOutput, LowerCrI_Output = LowerCI_SeroPrevOutput)
		{
			if (country < 10	)	PNGFILENAME = file.path(SeroPrevsPlotDirectory, paste0(FileNamePrefix, "SeroPrevs_c"		, country	))
			if (country == 10	)	PNGFILENAME = file.path(SeroPrevsPlotDirectory, paste0(FileNamePrefix, "SeroPrevs_CYD14"				))		## CYD 14
			if (country == 11	)	PNGFILENAME = file.path(SeroPrevsPlotDirectory, paste0(FileNamePrefix, "SeroPrevs_CYD15"				))		## CYD 15
			if (country == 12	)	PNGFILENAME = file.path(SeroPrevsPlotDirectory, paste0(FileNamePrefix, "SeroPrevs_CYD1415"				))		## CYD 14 and 15
			
			## i.e. if plotting proper posterior samples of seroprevalence (whether just augmented or augmented plus immunosubset) AND historical hazard samples
			if ((Plot_All_Patients | Plot_Augmented_Patients) & Plot_HH_Posteriors) 	PNGFILENAME = paste0(PNGFILENAME, "_Combo") 	else
			if (!(Plot_All_Patients | Plot_Augmented_Patients) & Plot_HH_Posteriors)	PNGFILENAME = paste0(PNGFILENAME, "_HHplot")
			PNGFILENAME = paste0(PNGFILENAME, ".png")
			
			if (SavePlot) png(file = PNGFILENAME, res = 160, units = "in", width = 5, height= 5)
			
			if (ShortenOutputString) ModelRun_string = AbbreviateOutputString(OutputSubDir, modelrun) else ModelRun_string = OutputSubDir
			PLOTTITLE = paste0(PlotPrefix, CountryNames_Long[country+1])
			if (OutputStringInTitle) PLOTTITLE = paste0(PLOTTITLE, "\n", ModelRun_string)
			
			### Find correct rows
			NonAugmented_Index 		= (country * 3) + 1
			Augmented_Index 		= (country * 3) + 2	 
			All_Index 				= (country * 3) + 3	 
			LegendLabels 			= c("All patients", "Augmented", "Data (Immune Subset)")
			LegendColours			= c("forestgreen", "deeppink", "blue")
			LegendLabels 			= c("Data (Immune Subset)")
			LegendColours			= c("blue")
			
			if ((country <= 4) | (country == 10)) 					ColsToDrop = c(0:1, 15:17	) + 1 ## plus 1 becuase of age zero.
			if ((country >= 5 & country <= 9) | (country == 11)) 	ColsToDrop = c(0:8, 17		) + 1
			if ((country == 12))									ColsToDrop = c(0:1, 17		) + 1
			
			MeanAug 		= na.omit(as.numeric(Mean_Output		[Augmented_Index,-ColsToDrop] ))
			UpperCIAug 		= na.omit(as.numeric(UpperCrI_Output	[Augmented_Index,-ColsToDrop] ))
			LowerCIAug 		= na.omit(as.numeric(LowerCrI_Output	[Augmented_Index,-ColsToDrop] ))
			
			Mean_All 		= na.omit(as.numeric(Mean_Output		[All_Index,-ColsToDrop] ))
			UpperCI_All 	= na.omit(as.numeric(UpperCrI_Output	[All_Index,-ColsToDrop] ))
			LowerCI_All 	= na.omit(as.numeric(LowerCrI_Output	[All_Index,-ColsToDrop] ))
			
			MeanNonAug		= na.omit(as.numeric(Mean_Output		[NonAugmented_Index,-ColsToDrop] ))
			UpperCINonAug 	= na.omit(as.numeric(UpperCrI_Output	[NonAugmented_Index,-ColsToDrop] ))
			LowerCINonAug 	= na.omit(as.numeric(LowerCrI_Output	[NonAugmented_Index,-ColsToDrop] ))
			
			ages 		= (0:17)[-ColsToDrop]
			
			# data (immune subset) - plot  again so more easily visible
			plot( ages, MeanNonAug, col = "blue", lwd = 4, 
					xlab = "Age", ylab = "proportion seropositive at baseline", 
					cex.lab = 1.8, cex.axis = 2, main= PLOTTITLE, type = "l", cex.main = TitleCex,
					ylim = c(0,1),)
			
			# augmented (non-immune subset)
			if (Plot_Augmented_Patients)
			{
				Aug_Col 			= "deeppink"
				LegendLabels 		= c(LegendLabels, "Augmented")
				LegendColours 		= c(LegendColours, Aug_Col)
				
				polygon	(c(ages,rev(ages)),c(LowerCIAug,rev(UpperCIAug)),col =  GetTransparentColours("pink") , border = NA)
				lines(ages, MeanAug, col = Aug_Col, lwd = 4)
			}
			
			# all patients (augmented and immuno subset)
			if (Plot_All_Patients)
			{
				All_patients_Col 	= "forestgreen"
				LegendLabels 		= c(LegendLabels, "All patients")
				LegendColours 		= c(LegendColours, All_patients_Col)
				
				polygon	(c(ages,rev(ages)),c(LowerCI_All,rev(UpperCI_All)),col= GetTransparentColours("darkseagreen1"), border = NA)
				lines(ages, Mean_All, col = All_patients_Col, lwd = 4)
			}
			
			if (Plot_HH_Posteriors & country < 9)
			{
				HH_col 				= "orange"
				LegendLabels 		= c(LegendLabels, "HHposteriors")
				LegendColours 		= c(LegendColours, HH_col)
				
				#### plot historical hazards by age. 
				Hist_Haz_chain_95 = CHAINS[WhichParamsToPlot, paste0("h_", country)]
				
				Min_HH 	= min (Hist_Haz_chain_95)
				Max_HH 	= max (Hist_Haz_chain_95)
				Mean_HH = mean(CHAINS[,paste0("h_", country)])
				
				polygon	(c(ages,rev(ages)),c(1 - exp( - Min_HH * ages),rev(1 - exp( - Max_HH * ages))),col = GetTransparentColours("peachpuff"), border = NA)
				lines(ages, 1 - exp( - Mean_HH * ages), col = "orange", lwd = 4)
			}
			
			# data (immune subset) - plot  again so more easily visible
			lines(ages, MeanNonAug, col = "blue", lwd = 4)
			if (IncLegend)
				legend("bottomright", LegendLabels,
						pch = rep(NA,length(LegendLabels)) ,
						col = LegendColours, lty = rep(1,length(LegendLabels)), lwd = rep(2,length(LegendLabels)) , cex = LegCex	)
			
			if (SavePlot) dev.off()
		}

		PlotMultipleSeroPrevs = function(Countries = NULL, Plot_All_Patients = TRUE, Plot_Augmented_Patients = TRUE, Plot_HH_Posteriors = FALSE, SavePlot = TRUE, Overwrite = OVERWRITE_FILES, 
				ShortenOutputString = TRUE, modelrun = ModelRun, N_Rows = 3, N_Cols = 4, AddMultiPlotTitle = TRUE, PlotToTitleHeightRatio = 3, IndividualTitleCex = 3, MultiPlot_Title = NULL, FileNamePrefix = "", 
				Mean_Output = Mean_SeroPrevOutput, UpperCrI_Output = UpperCI_SeroPrevOutput, LowerCrI_Output = LowerCI_SeroPrevOutput)
		{
			if (is.null(Countries)) 
			{
				Countries 		= 0:11 
				CountriesString = "_AllCountries"

			}	else	
			{
				CountriesString = paste0("_Cs", paste0(Countries, collapse = ""))
			}
			
			PNGFILENAME = file.path(SeroPrevsPlotDirectory, paste0(FileNamePrefix, "SeroPrevs", CountriesString))
			## i.e. if plotting proper posterior samples of seroprevalence (whether just augmented or augmented plus immunosubset) AND historical hazard samples
			if ((Plot_All_Patients | Plot_Augmented_Patients) & Plot_HH_Posteriors) PNGFILENAME = paste0(PNGFILENAME, "_Combo") else
			if (!(Plot_All_Patients | Plot_Augmented_Patients) & Plot_HH_Posteriors)PNGFILENAME = paste0(PNGFILENAME, "_HHplot")
			PNGFILENAME = paste0(PNGFILENAME, ".png")
			
			IndividualPlotWidth = IndividualPlotHeight = 5
			
			if (SavePlot) png(file = PNGFILENAME, res = PNG_res, units = "in", width = IndividualPlotWidth * N_Cols, height = IndividualPlotHeight * N_Rows + (IndividualPlotHeight/PlotToTitleHeightRatio))
			
			LayoutMatrix = MakeLayoutMatrix(NoCols = N_Cols, NoRows = N_Rows, AddTitle = AddMultiPlotTitle, AddColNames = FALSE, AddRowNames = FALSE, BY_ROW = TRUE)
			if (AddMultiPlotTitle) HEIGHTS = c(1,rep(PlotToTitleHeightRatio, N_Rows)) else HEIGHTS = NULL
			layout(LayoutMatrix, heights = HEIGHTS)
			
			if (AddMultiPlotTitle)
			{	
				par(mar = rep(0, 4)) ### without this you often get the annoying "Error in plot.new() : figure margins too large" error. Set margins to zero for the title, then reset them to OrigMAR (global variable defined at start of script). 
				plot.new()
				if (is.null(MultiPlot_Title))	
				{	
					MultiPlot_Title = ""
					if (ShortenOutputString) MultiPlot_Title = AbbreviateOutputString(OutputSubDir, modelrun) else MultiPlot_Title = OutputSubDir	
					MultiPlot_Title = paste0( MultiPlot_Title, "\nPredicted Seroprevalences")
				}
				text(0.5, 0.5, MultiPlot_Title , cex = 5, font = 2)
			}
			
			### Include legend only in final plot
			IncLegends = rep(FALSE, length(Countries))
			IncLegends[length(Countries)] = TRUE
			
			par(mar = OrigMAR) ### reset default margins or S_Curves will look horrid. 
			if (AddMultiPlotTitle) OutputStringInTitle = FALSE else OutputStringInTitle = TRUE
			for (country in Countries)
				PlotSeroPrevs(country, Plot_All_Patients = Plot_All_Patients, Plot_Augmented_Patients = Plot_Augmented_Patients, Plot_HH_Posteriors = Plot_HH_Posteriors, 
						modelrun = modelrun, Overwrite = OVERWRITE_FILES, OutputStringInTitle = OutputStringInTitle, SavePlot = FALSE, TitleCex = IndividualTitleCex, PlotPrefix = "", 
						FileNamePrefix = FileNamePrefix, LegCex = 2.2, IncLegend = IncLegends[country + 1],
						Mean_Output = Mean_Output, UpperCrI_Output = UpperCrI_Output, LowerCrI_Output = LowerCrI_Output)
			
			if (SavePlot) dev.off()
		}
		
		PlotMultipleSeroPrevs(MultiPlot_Title = "Age-specific seroprevalence by country & trial")
		cat ("\t")
	} 	## end do seroprevalences if statement
	
	CloseOpenPlotDevices()	
	
	##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
	##### ##### ##### Attack Rates
	
	if (DoAttackRates)	
	{		
		AttackRatePlotDirectory = file.path(ModelPlotSubDir, "AttackRates")
		dir.create(AttackRatePlotDirectory, showWarnings = FALSE)
		
		cat ("Attack Rates: ")

		#### Attack rate functions in script AttackRates.R
		whichdiseaseplot = WhichDiseasePlots[1]
		for (whichdiseaseplot in WhichDiseasePlots)
		{
			country = 0
			for (country in 0:11)
			{
				if (country == 12 & ModelRun$CountriesFittedString == "_Cs01234") next
				if (country == 12 & ModelRun$CountriesFittedString == "_Cs56789") next
				
				SkipCountry = ShouldCountryBeSkipped(ModelRun)
				if (SkipCountry) next
				
				PhaseChar = "Both"
				for (PhaseChar in c("Both", "Active", "Passive"))
				{
					try(PlotAttackRates(country, Phase = PhaseChar, WhichDisease = whichdiseaseplot), silent = TRUE)
					if (country %in% 10:11)
					{
						for (BS in c("SeroNeg", "SeroPos"))
							try(PlotAttackRates(country, Phase = PhaseChar, BaselineSeroStatus = BS, ImSub = TRUE, WhichDisease = whichdiseaseplot), silent = TRUE) ## error-trapped as some combinations of trial phase, disease severity, country and serostatus won't work (e.g. mild cases in passive phase). 
						rm(BS)	
					}
				}
				rm(PhaseChar)	
			}
		}
		
		cat ("\t")
	}  ## end do attack rates if statement
	
	CloseOpenPlotDevices()	
	
	##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
	##### ##### ##### AgeEfficacy Plots
		
	if (DoAgeEffectPlots & (ModelRun$ASVE  != "INDEPENDENT"  || ModelRun$AS_Haz  != "INDEPENDENT" || ModelRun$AS_Waning != "INDEPENDENT")) 
	{
		AgeEffectsDirectory = ModelPlotSubDir
		dir.create(AgeEffectsDirectory, showWarnings = FALSE)		
		
		## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
		## Age specific Efficacies
		if (ModelRun$ASVE 	!= "INDEPENDENT" )
		{
			if (ModelRun$ASVE == "HILL") 
			{
				cat ("AgeEfficacyPlots: HILL\n")
				BaselineSeroStatus 		= "SeroNeg"
				serotype 				= 1
				### Make "dual" plot 
				if (ModelRun$SS_VEs & ModelRun$SType_Equiv) N_STypes_VEs_dummy = 1 else N_STypes_VEs_dummy = N_STypes_VEs ## this save doing unneccessary costly runs when serotypes all the same. 
				
				for (serotype in 1:N_STypes_VEs_dummy)
				{
					if (ModelRun$SS_VEs & !ModelRun$SType_Equiv) SerotypeString = paste0("_serotype_", serotype) else SerotypeString = ""
					ASVEHill_PNGFILENAME = file.path(AgeEffectsDirectory, paste0("AgeEffHill_DualPlot", SerotypeString, ".png"))
					png (file = ASVEHill_PNGFILENAME, res = PNG_res, units = "in", width = 10, height = 5)
					par(mar=rep(0,4)) # no margins order: bottom, left, top, and right. The default is c(5.1, 4.1, 4.1, 2.1).
					LayoutMatrix = matrix (c(1,2,1,3), ncol = 2)
					layout(LayoutMatrix, heights = c(1,7))
					plot.new()
					text(0.5, 0.5, AbbreviateOutputString(OutputSubDir, ModelRun) , cex = 1.5, font = 2)
					par(mar = c(5.1, 4.1, 4.1, 2.1))
					for (BaselineSeroStatus in c("SeroNeg", "SeroPos"))
						PlotHill(PlotType = "Age_Efficacy", BaselineSeroStatus = BaselineSeroStatus, serotype = serotype, modelrun = ModelRun, SavePlot = FALSE, OutputSubDirInPlotTite = FALSE)	
					dev.off()
				}
				
			} else 
			{
				cat (paste0("AgeEfficacyPlots: SPLINE\n"))
				serotype = 1
				for (serotype in 1:N_STypes_VEs)
				{
					if (ModelRun$SS_VEs == 1) SerotypeString = paste0("_serotype_", serotype) else SerotypeString = ""
					ASVE_SplinePNGFILENAME = file.path(AgeEffectsDirectory, paste0("ASVEs", SerotypeString, ".png"))						
					png (file  = ASVE_SplinePNGFILENAME, res = PNG_res, units = "in", width = 10, height= 5 )
					par(mar=rep(0,4)) # no margins. order: bottom, left, top, and right. default is c(5.1, 4.1, 4.1, 2.1).
					LayoutMatrix = matrix (c(1,2,1,3), ncol = 2)
					layout(LayoutMatrix, heights = c(1,7))
					plot.new()
					if (ModelRun$SS_VEs == 1) text(0.5, 0.5, paste("serotype", serotype) , cex = 3, font = 2)
					par(mar = c(5.1, 4.1, 4.1, 2.1))
					for (BaselineSeroStatus in c("SeroNeg","SeroPos"))
					{
						if (ModelRun$ParamRangeFileName == "prs3_2") LOWERYLIM = 0 else LOWERYLIM = -1
						if (BaselineSeroStatus == "SeroNeg") BSS_string = "seronegative" 	else BSS_string = "seropositive" 
						if (BaselineSeroStatus == "SeroPos") IncLegend 	= TRUE 				else IncLegend 	= FALSE
						PlotSplines(BaselineSeroStatus = BaselineSeroStatus, modelrun = ModelRun, PlotType = "Age_Efficacy", serotype = serotype, 
								SavePlot = FALSE, ## need SavePlot = FALSE here for dual plot. 
								PLOTTITLE = BSS_string, LOWERYLIM = LOWERYLIM, IncLegend = IncLegend, LegPosChar = "bottomright", 
								#NumParamSetsToPlot = 10, 
								Directory = AgeEffectsDirectory, OutputSubDirInPlotTite = FALSE)
					}
					dev.off()
				}
			}
		}
		
		## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
		## Age specific hazards
		if (ModelRun$AS_Haz != "INDEPENDENT")
		{
			cat (paste0("AS_Haz: SPLINE or HILL\n"))
			
			if (ModelRun$AS_Haz != "HILL")		
				PlotSplines(modelrun = ModelRun, PlotType = "Age_Haz_Mult", Directory = AgeEffectsDirectory, OutputSubDirInPlotTite = FALSE) else
				PlotHill(PlotType = "Age_Haz_Mult", modelrun = ModelRun)
		}
		CloseOpenPlotDevices()		
		
		## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
		## Age specific vaccine durations / waning
		if (ModelRun$AS_Waning 	!= "INDEPENDENT")
		{
			if (ModelRun$AS_Waning == "HILL") ### otherwise just then you just use posteriors. 
			{
				cat ("AgeWaningPlots: HILL\n")
				
				BaselineSeroStatus 		= "SeroNeg"
				### Make "dual" plot 
				ASWaningHill_PNGFILENAME = file.path(AgeEffectsDirectory, paste0("AgeDurationHill_DualPlot.png"))
				
				png (file  = ASWaningHill_PNGFILENAME, res = PNG_res, units = "in", width = 10, height= 5 )
				par(mar=rep(0,4)) # no margins order: bottom, left, top, and right. The default is c(5.1, 4.1, 4.1, 2.1).
				LayoutMatrix = matrix (c(1,2,1,3), ncol = 2)
				layout(LayoutMatrix, heights = c(1,7))
				plot.new()
				text(0.5, 0.5, AbbreviateOutputString(OutputSubDir, ModelRun) , cex = 1.5, font = 2)
				par(mar= c(5.1, 4.1, 4.1, 2.1))
				for (BaselineSeroStatus in c("SeroNeg","SeroPos"))
					PlotHill(PlotType = "Age_VacDuration", BaselineSeroStatus = BaselineSeroStatus, modelrun = ModelRun, SavePlot = FALSE, OutputSubDirInPlotTite = FALSE)	
				dev.off()
				
			} else
			{
				cat ("AgeWaningPlots: SPLINE \n")
				AS_WaningPNGFILENAME = file.path(AgeEffectsDirectory, paste0("AS_Waning_Splines", ".png"))
				
				png (file  = AS_WaningPNGFILENAME, res = PNG_res, units = "in", width = 10, height= 5 )
				par(mar=rep(0,4)) # no margins order: bottom, left, top, and right. The default is c(5.1, 4.1, 4.1, 2.1).
				LayoutMatrix = matrix (c(1,2,1,3), ncol = 2)
				layout(LayoutMatrix, heights = c(1,7))
				plot.new()
				par(mar= c(5.1, 4.1, 4.1, 2.1))
				for (BaselineSeroStatus in c("SeroNeg","SeroPos"))
				{
					if (BaselineSeroStatus == "SeroNeg") BSS_string = "seronegative" 	else BSS_string = "seropositive" 
					if (BaselineSeroStatus == "SeroPos") IncLegend 	= TRUE 				else IncLegend 	= FALSE
					PlotSplines(BaselineSeroStatus = BaselineSeroStatus, modelrun = ModelRun, # NumParamSetsToPlot = 1000,
							PlotType = "Age_VacDuration", SavePlot = FALSE, 
							PLOTTITLE = BSS_string, 
							OutputSubDirInPlotTite = FALSE, Directory = AgeEffectsDirectory)
				}
				dev.off()
			}
		}
		
		CloseOpenPlotDevices()		
	}	
	
	CloseOpenPlotDevices()	
	cat ("DONE \n")
	
	FilenamesCompleted = c(FilenamesCompleted, OutputString)
}

print ("All Done")
#FilenamesCompleted	

### which Filenames didn't run/work?
setdiff(paste0("_", as.character(ModelRuns$OutputFolderNames)) , FilenamesCompleted)




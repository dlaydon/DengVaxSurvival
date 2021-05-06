# TODO: Add comment
# 
# Author: dlaydon
###############################################################################


GetOutputtedPredicted_SProbs 	= function(Statistic = "Mean", Passive = FALSE, WhichDiseasePlot = "", ImSub = FALSE, HazGroup_Indices, WBIC_String = "")
{
	if (!(WBIC_String %in% c("", "WBIC_"))) 										stop("GetOutputtedPredicted_SProbs error: WBIC_String not recognized")
	if (!(Statistic %in% c("Mean", "LowerCI", "UpperCI", "Modal", "MaxLike"))) 	stop("GetOutputtedPredicted_SProbs error: Statistic not recognized")
	ImSubString = GetImSubString(ImSub)
	if (Passive) 		PassiveOrNot_Char	= "Passive" else PassiveOrNot_Char = ""  
	
	return(na.omit(as.numeric(get(paste0(WBIC_String, Statistic, "_", PassiveOrNot_Char, "SurvivalTable"	, WhichDiseasePlot, ImSubString))	[HazGroup_Indices,] ))) #### note that even though you named the FILES inconsistently between Modal, MaxLike and canonical Mean/LowerCI/UpperCI, you named the VARIABLES in this R session consistently. Therefore this line doesn't need to accont for differences between potential values of Statistic argument string
}
Adjust_Either_PassiveSProbs 	= function(Either_SProbVec, Statistic, Passive = FALSE, ImSub = FALSE, HazGroup_Indices, modelrun = ModelRun, WBIC_String = "")
{
	if (!(Statistic %in% c("Mean", "LowerCI", "UpperCI"))) stop("Adjust_Either_PassiveSProbs error: Statistic not recognized")
	if (modelrun$disease == "") stop("Adjust_Either_PassiveSProbs error: modelrun$disease != Mild severe. Function shouldn't be called ")
	
	## Amend all of them post NoDaysActive and Passive. 
	### Get final probability of surviving mild disease. Multiply all Severe disease by this probabilities
	
	StartIndex 		= NoDaysActive + 2	## + 2 as you want the first day of passive onwards (i.e. one day after last active day), and also the indexes go from day 0.
	EndIndex 		= length(Either_SProbVec) 
	IndicesToChange = StartIndex:EndIndex
	
	### Get final probability of surviving mild disease.
	Mild_SProbVec 	= GetOutputtedPredicted_SProbs(Statistic = Statistic, Passive = Passive, WhichDiseasePlot = "_Mild"		, ImSub = ImSub, HazGroup_Indices = HazGroup_Indices, WBIC_String = WBIC_String)
	FinalMildProb 	= Mild_SProbVec[NoDaysActive + 1]
	
	### Get Prob Severe disease vec. 
	Severe_SProbVec = GetOutputtedPredicted_SProbs(Statistic = Statistic, Passive = Passive, WhichDiseasePlot = "_Severe"	, ImSub = ImSub, HazGroup_Indices = HazGroup_Indices, WBIC_String = WBIC_String)
	
	## Amend All indices post active phase
	Either_SProbVec[IndicesToChange] =  Severe_SProbVec[IndicesToChange] * FinalMildProb
	
	return(Either_SProbVec)
}


ChooseMaxTi_ForSurvalCurve 			= function(Whichdiseaseplot = "", modelrun = ModelRun, MaxTi = NULL, Either_FullDuration = FALSE)
{
	if (is.null(MaxTi))
	{
		if 	(modelrun$disease == "_MILDSEVERE")	
		{
			if 	(modelrun$phase == "_PASSIVE"		)
			{
				if (Whichdiseaseplot == "_Either"	) 
				{
					if (!Either_FullDuration) MaxTi = NoDaysActive 	else  	MaxTi = NoDays_Active_And_Passive
				}
				if (Whichdiseaseplot == "_Mild"		) MaxTi = NoDaysActive 	#### NoDays_Active_And_Passive 	takes account of SFU
				if (Whichdiseaseplot == "_Severe"	) MaxTi = NoDays_Active_And_Passive 	#### NoDays_Active_And_Passive 	takes account of SFU
				
			} else if (modelrun$phase == "_ACTIVE"		)	MaxTi = NoDaysActive else if 	(modelrun$phase == "_PASSIVE_ONLY"	)	MaxTi = NoDaysPassive
			
		} else	
		{
			if 	(modelrun$phase == "_PASSIVE"		)	MaxTi = NoDays_Active_And_Passive	else
			if 	(modelrun$phase == "_ACTIVE"		)	MaxTi = NoDaysActive				else	#### NoDaysActive 				takes account of SFU
			if 	(modelrun$phase == "_PASSIVE_ONLY"	)	MaxTi = NoDaysPassive						#### NoDaysPassive 				takes account of SFU
		}
	}
	return( MaxTi )
}
ProbSurviveOverTime_stratum 		= function(country, Arm = "Control", AgeGroup = 0, Whichdiseaseplot = "", modelrun = ModelRun, MaxTi = ChooseMaxTi_ForSurvalCurve(Whichdiseaseplot, modelrun), 
		SeroStatus = "Both_Sero", countriesfitted = CountriesFitted, Passive = FALSE, AgeMin = NULL, AgeMax = NULL)
{
	
	### subset data: 
	PatientsToSelect 	= ChooseRowIndices	(country = country, CountriesFitted = countriesfitted, AgeGroup = AgeGroup, Arm = Arm, SeroStatus = SeroStatus, AgeMin = AgeMin, AgeMax = AgeMax)
	if (length(PatientsToSelect) == 0 ) { warning("ProbSurviveOverTime_stratum: length(PatientsToSelect) == 0") ; return(NA) }
	DummyData 			= DATA[PatientsToSelect, ]
	
	if (!Passive)	
	{
		CaseColName 				= ChooseCaseColName (Whichdiseaseplot, "Both", modelrun$disease, modelrun$hosp)
		NumParticipants_stratum 	= length(PatientsToSelect)
		DaysPostStart_Phase			= DummyData$DaysPostFirstDose #### DaysPostFirstDose is relative to StartFollowUpDate_InDays, i.e. start of trial / active phase. Need another vector that gives DaysPostStartPassive
		Start_Phase					= DummyData$StartFollowUpDate_InDays
		
	}	else
	{
		#### reset MaxTi
		MaxTi 						= NoDays_Active_And_Passive - NoDaysActive
		
		CaseColName 				= ChooseCaseColName (Whichdiseaseplot, "Passive", modelrun$disease, modelrun$hosp)
		NumParticipants_stratum 	= length(which(DummyData$IsCase_Active == 0 & DummyData$LTFU == 0)) ### i.e. for passive phase survival curves exclude LTFU and active phase cases. 
		
		#### DaysPostFirstDose is relative to StartFollowUpDate_InDays, i.e. start of trial / active phase. Need another vector that gives DaysPostStartPassive
		DaysPostStart_Phase 		= DummyData$EndPassivePhase_InDays - DummyData$StartPassivePhase_InDays
		Start_Phase					= DummyData$StartPassivePhase_InDays
	}
	
	if (any(DaysPostStart_Phase < 0)) stop("ProbSurviveOverTime_stratum error: any(DaysPostStart_Phase < 0) TRUE")
	CaseTimes_stratum			= DaysPostStart_Phase[which(DummyData[, CaseColName] == 1)] 	
	
	NumSurvivorsLeft_stratum 	= rep(NA, MaxTi+1)
	for (day in 0:MaxTi)
		NumSurvivorsLeft_stratum [day+1] = NumParticipants_stratum - length(which(CaseTimes_stratum <= day) )
	
	return(NumSurvivorsLeft_stratum / NumParticipants_stratum )
}


#### At some point, make all these irritating control arguments into a list. Too many
PlotSurvivalCurves = function(country, AgeGroup = 0, WhichDiseasePlot = "", modelrun = ModelRun, MaxTi = NULL, ImSub = FALSE,
		CurveToPlot = "ByTrialArm", countriesfitted = CountriesFitted, Passive = FALSE, 
		SavePlot = TRUE, 
		IncludePlotTitle = TRUE, PLOTTITLE = NULL, 
		OutputStringInTitle = TRUE,  #### NOTE: I changed code so that OutputStringInTitle = FALSE is the default to avoid typing it here 10 billion times. It's necessary for model comparison though, so change it back. 
		#OutputStringInTitle = FALSE,  #### NOTE: I changed code so that OutputStringInTitle = FALSE is the default to avoid typing it here 10 billion times. It's necessary for model comparison though, so change it back. 
		ExpectedCasesInPlotTitle = TRUE, ShortenOutputString = TRUE,
		CEXLAB = 1.4, CEXAXIS = 1.2, LOWERYLIM = NULL, TitleCex = 1, LEGCEX = 1.3,
		AgeMin = NULL, AgeMax = NULL, Either_FullDuration = TRUE, 
		Directory = SurvivalCurvePlotDirectory, FileName = NULL, IncAxesLabels = TRUE, Plot_Data = TRUE, Plot_Model = TRUE, 
		CIsObs = FALSE, Plot_ModalPost = FALSE, Plot_MaxLike = FALSE, WBIC_String = "", IncABlines = FALSE, IncDataLegend = FALSE, LWD = 4)
{
	# SavePlot 	= FALSE will make the plot but not save it.  
	# Either_FullDuration argument only used when modelrun$disease == "_MILDSEVERE" and WhichDiseasePlot = "_Either". Changes plot title and filename as well as those detailed below
	# 	i) TRUE uses outputted predicted CPP survival probabilites, with no amendment. Also sets MaxTi to NoDaysActive using ChooseMaxTi_ForSurvalCurve function
	# 	ii) FALSE (new default) uses outputted predicted CPP survival probabilites but amends after NoDaysActive so that P(SurviveEither) = P(Survived Active Phase Mild disease) * P(Survive Severe). Also sets MaxTi to NoDays_Active_And_Passive using ChooseMaxTi_ForSurvalCurve function
	
	if (!any( c("ByTrialArm", "BySeroStatus", "_SPos", "_SNeg") == CurveToPlot)) 	stop("PlotSurvivalCurves error: CurveToPlot argument not recognized. ")
	if (!any( c(TRUE, FALSE) == Passive)) 											stop("PlotSurvivalCurves error: Passive argument not recognized. ")
	
	if (is.null(MaxTi)) MaxTi = ChooseMaxTi_ForSurvalCurve(Whichdiseaseplot = WhichDiseasePlot, modelrun = modelrun, Either_FullDuration = Either_FullDuration)
	
	### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
	### Build PNGFILENAME
	### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
	
	if (Passive		 ) 	PassiveOrNot_Char 		= "_Pass" 					else PassiveOrNot_Char = ""  
	CountryBitOfPNGTitle = ChooseCountryBitOfPNGTitle(country)
	ImSubString = GetImSubString(ImSub)
	if (CIsObs) 		CIsObsString 			= "_ObsCIs" 	else CIsObsString = ""
	if (CurveToPlot == "ByTrialArm") WhichCurve_String = "" else WhichCurve_String = CurveToPlot
	
	if (WhichDiseasePlot == "_Either" & Either_FullDuration)	DiseaseString = paste0(WhichDiseasePlot, "FullTrial") else DiseaseString = WhichDiseasePlot
	
	RootFileName	= paste0(WBIC_String, "SC_AG", AgeGroup, CountryBitOfPNGTitle, DiseaseString, WhichCurve_String, PassiveOrNot_Char, ImSubString, CIsObsString)
	if (Plot_MaxLike) 	RootFileName = paste0(RootFileName, "_ML")
	if (Plot_ModalPost) RootFileName = paste0(RootFileName, "_MP")
	if (SavePlot) FileName 		= file.path(Directory, paste0(RootFileName, ".png"))
	
	
	### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
	### ### ### ### ### 	Choose Data and Predictions to plot.  
	### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
	
	CountryAgeGroupIndex = (AgeGroup * 13 * 9) + (country * 9) #  * 13 countries * 9 hazard groups
	
	## In R indices, hazard groups are: 1 = SeroPos Control, 2 = SeroNeg Control, 3 = SeroPos Vaccine, 4 = SeroNeg Vaccine, 5 = Control,  6 = Vaccine,  7 = SeroPositive, 8 = SeroNegative, 
	
	if (CurveToPlot == "ByTrialArm")
	{
		HazGroup_1_Indices = CountryAgeGroupIndex + 5 ## control 
		HazGroup_2_Indices = CountryAgeGroupIndex + 6 ## vaccine
		
		ArmArg_1 =  "Control"
		ArmArg_2 =  "Vaccine"
		
		if (!ImSub) SeroStatusArg_1 = "Both_Sero" else SeroStatusArg_1 = "_ImSub" 	#### For ImSub, predicted survival curves are based only on immune subset - even "vaccine" and "control" arms which are nominally independent of serostatus. "Both_Sero" will select everybody in the stratum regardless of serostatus, including patients with unknown serostatus. "_ImSub" will select patients in immune subset only. 
		if (!ImSub) SeroStatusArg_2 = "Both_Sero" else SeroStatusArg_2 = "_ImSub"	#### For ImSub, predicted survival curves are based only on immune subset - even "vaccine" and "control" arms which are nominally independent of serostatus. "Both_Sero" will select everybody in the stratum regardless of serostatus, including patients with unknown serostatus. "_ImSub" will select patients in immune subset only. 
		
		LegLabels = c( "Vaccine", "Control")
		
	} else if (CurveToPlot == "BySeroStatus")
	{
		HazGroup_1_Indices = CountryAgeGroupIndex + 7 ## SeroPositive 
		HazGroup_2_Indices = CountryAgeGroupIndex + 8 ## SeroNegative
		
		ArmArg_1 =  "Either"
		ArmArg_2 =  "Either"
		
		SeroStatusArg_1 = "SeroNeg"
		SeroStatusArg_2 = "SeroPos"
		
		LegLabels = c("SeroNegative", "SeroPositive")
		
	} else if (CurveToPlot == "_SPos")
	{
		HazGroup_1_Indices = CountryAgeGroupIndex + 1 ## SeroPositive control
		HazGroup_2_Indices = CountryAgeGroupIndex + 3 ## SeroPositive vaccine
		
		ArmArg_1 =  "Control"
		ArmArg_2 =  "Vaccine"
		
		SeroStatusArg_1 = "SeroPos"
		SeroStatusArg_2 = "SeroPos"
		
		LegLabels = c( "SeroPos Vaccine", "SeroPos Control")
		
	} else if (CurveToPlot == "_SNeg")
	{
		HazGroup_1_Indices = CountryAgeGroupIndex + 2 ## SeroNegative control 
		HazGroup_2_Indices = CountryAgeGroupIndex + 4 ## SeroNegative vaccine
		
		ArmArg_1 =  "Control"
		ArmArg_2 =  "Vaccine"
		
		SeroStatusArg_1 = "SeroNeg"
		SeroStatusArg_2 = "SeroNeg"
		
		LegLabels = c( "SeroNeg Vaccine", "SeroNeg Control")
	}		
	
	#### calculate observed
	if (Plot_Data)
	{
		Data_Haz_Group_1 = ProbSurviveOverTime_stratum(country = country, Arm = ArmArg_1, AgeGroup = AgeGroup, SeroStatus = SeroStatusArg_1, Whichdiseaseplot = WhichDiseasePlot, MaxTi = MaxTi, countriesfitted = countriesfitted, Passive = Passive, AgeMin = AgeMin, AgeMax = AgeMax, modelrun = modelrun)        
		Data_Haz_Group_2 = ProbSurviveOverTime_stratum(country = country, Arm = ArmArg_2, AgeGroup = AgeGroup, SeroStatus = SeroStatusArg_2, Whichdiseaseplot = WhichDiseasePlot, MaxTi = MaxTi, countriesfitted = countriesfitted, Passive = Passive, AgeMin = AgeMin, AgeMax = AgeMax, modelrun = modelrun)         
		NumParticipants_ControlArm_country 	= length(ChooseRowIndices	(country = country, CountriesFitted = countriesfitted, AgeGroup = AgeGroup, Arm = ArmArg_1, SeroStatus = SeroStatusArg_1, AgeMin = AgeMin, AgeMax = AgeMax))
		NumParticipants_VaccineArm_country 	= length(ChooseRowIndices	(country = country, CountriesFitted = countriesfitted, AgeGroup = AgeGroup, Arm = ArmArg_2, SeroStatus = SeroStatusArg_2, AgeMin = AgeMin, AgeMax = AgeMax))
		if (CIsObs)
		{
			CIs_Data_Haz_Group_1 = binconf(Data_Haz_Group_1 * NumParticipants_ControlArm_country, NumParticipants_ControlArm_country, method = "exact")
			CIs_Data_Haz_Group_2 = binconf(Data_Haz_Group_2 * NumParticipants_VaccineArm_country, NumParticipants_VaccineArm_country, method = "exact")
		}
	}
	
	#### choose predicted
	if (Plot_Model)
	{
		Mean_HazGroup_1			= GetOutputtedPredicted_SProbs(Statistic = "Mean"	, Passive = Passive, WhichDiseasePlot = WhichDiseasePlot, ImSub = ImSub, HazGroup_Indices = HazGroup_1_Indices) ### for now at least, don't add in WBIC_String here - you aren't outputting Mean and 95%CrI's for the WBIC chains.
		LowerCI_HazGroup_1 		= GetOutputtedPredicted_SProbs(Statistic = "LowerCI", Passive = Passive, WhichDiseasePlot = WhichDiseasePlot, ImSub = ImSub, HazGroup_Indices = HazGroup_1_Indices)
		UpperCI_HazGroup_1 		= GetOutputtedPredicted_SProbs(Statistic = "UpperCI", Passive = Passive, WhichDiseasePlot = WhichDiseasePlot, ImSub = ImSub, HazGroup_Indices = HazGroup_1_Indices)
		Mean_HazGroup_2 		= GetOutputtedPredicted_SProbs(Statistic = "Mean"	, Passive = Passive, WhichDiseasePlot = WhichDiseasePlot, ImSub = ImSub, HazGroup_Indices = HazGroup_2_Indices)
		LowerCI_HazGroup_2 		= GetOutputtedPredicted_SProbs(Statistic = "LowerCI", Passive = Passive, WhichDiseasePlot = WhichDiseasePlot, ImSub = ImSub, HazGroup_Indices = HazGroup_2_Indices)
		UpperCI_HazGroup_2 		= GetOutputtedPredicted_SProbs(Statistic = "UpperCI", Passive = Passive, WhichDiseasePlot = WhichDiseasePlot, ImSub = ImSub, HazGroup_Indices = HazGroup_2_Indices)
		
		if (Plot_ModalPost)
		{
			ModalPost_HazGroup_1	= GetOutputtedPredicted_SProbs(Statistic = "Modal", Passive = Passive, WhichDiseasePlot = WhichDiseasePlot, ImSub = ImSub, HazGroup_Indices = HazGroup_1_Indices, WBIC_String = WBIC_String)
			ModalPost_HazGroup_2	= GetOutputtedPredicted_SProbs(Statistic = "Modal", Passive = Passive, WhichDiseasePlot = WhichDiseasePlot, ImSub = ImSub, HazGroup_Indices = HazGroup_2_Indices, WBIC_String = WBIC_String)
		}
		if (Plot_MaxLike)
		{
			MaxLike_HazGroup_1	= GetOutputtedPredicted_SProbs(Statistic = "MaxLike", Passive = Passive, WhichDiseasePlot = WhichDiseasePlot, ImSub = ImSub, HazGroup_Indices = HazGroup_1_Indices, WBIC_String = WBIC_String)
			MaxLike_HazGroup_2	= GetOutputtedPredicted_SProbs(Statistic = "MaxLike", Passive = Passive, WhichDiseasePlot = WhichDiseasePlot, ImSub = ImSub, HazGroup_Indices = HazGroup_2_Indices, WBIC_String = WBIC_String)
		}
		
		if (WhichDiseasePlot == "_Either" & Either_FullDuration)
		{
			### Adjust survial probabilities from "_Either" survival table
			Mean_HazGroup_1			= Adjust_Either_PassiveSProbs	(Either_SProbVec = Mean_HazGroup_1		, Statistic = "Mean"	, Passive = Passive, ImSub = ImSub, HazGroup_Indices = HazGroup_1_Indices)
			LowerCI_HazGroup_1 		= Adjust_Either_PassiveSProbs	(Either_SProbVec = LowerCI_HazGroup_1 	, Statistic = "LowerCI"	, Passive = Passive, ImSub = ImSub, HazGroup_Indices = HazGroup_1_Indices)
			UpperCI_HazGroup_1 		= Adjust_Either_PassiveSProbs	(Either_SProbVec = UpperCI_HazGroup_1 	, Statistic = "UpperCI"	, Passive = Passive, ImSub = ImSub, HazGroup_Indices = HazGroup_1_Indices)
			Mean_HazGroup_2 		= Adjust_Either_PassiveSProbs	(Either_SProbVec = Mean_HazGroup_2 		, Statistic = "Mean"	, Passive = Passive, ImSub = ImSub, HazGroup_Indices = HazGroup_2_Indices)
			LowerCI_HazGroup_2 		= Adjust_Either_PassiveSProbs	(Either_SProbVec = LowerCI_HazGroup_2 	, Statistic = "LowerCI"	, Passive = Passive, ImSub = ImSub, HazGroup_Indices = HazGroup_2_Indices)
			UpperCI_HazGroup_2 		= Adjust_Either_PassiveSProbs	(Either_SProbVec = UpperCI_HazGroup_2 	, Statistic = "UpperCI"	, Passive = Passive, ImSub = ImSub, HazGroup_Indices = HazGroup_2_Indices)
			
			if (Plot_ModalPost)
			{
				ModalPost_HazGroup_1	= Adjust_Either_PassiveSProbs(Either_SProbVec = ModalPost_HazGroup_1, Statistic = "Modal", Passive = Passive, WhichDiseasePlot = WhichDiseasePlot, ImSub = ImSub, HazGroup_Indices = HazGroup_1_Indices, WBIC_String = WBIC_String)
				ModalPost_HazGroup_2	= Adjust_Either_PassiveSProbs(Either_SProbVec = ModalPost_HazGroup_2, Statistic = "Modal", Passive = Passive, WhichDiseasePlot = WhichDiseasePlot, ImSub = ImSub, HazGroup_Indices = HazGroup_2_Indices, WBIC_String = WBIC_String)
			}
			if (Plot_MaxLike)
			{
				MaxLike_HazGroup_1	= Adjust_Either_PassiveSProbs(Either_SProbVec = MaxLike_HazGroup_1, Statistic = "MaxLike", Passive = Passive, WhichDiseasePlot = WhichDiseasePlot, ImSub = ImSub, HazGroup_Indices = HazGroup_1_Indices, WBIC_String = WBIC_String)
				MaxLike_HazGroup_2	= Adjust_Either_PassiveSProbs(Either_SProbVec = MaxLike_HazGroup_2, Statistic = "MaxLike", Passive = Passive, WhichDiseasePlot = WhichDiseasePlot, ImSub = ImSub, HazGroup_Indices = HazGroup_2_Indices, WBIC_String = WBIC_String)
			}
		} 
		### Expected number of cases
		ExpNoCasesAtEnd_Control 			= NumParticipants_ControlArm_country * (1 - tail(Mean_HazGroup_1, 1))
		ExpNoCasesAtEnd_Vaccine 			= NumParticipants_VaccineArm_country * (1 - tail(Mean_HazGroup_2, 1))
	}
	
	### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
	### Build PlotTitle
	### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
	
	DiseaseString = Convert_DiseaseString(WhichDiseasePlot, modelrun$disease, modelrun$hosp)
	if (WhichDiseasePlot == "_Either" & Either_FullDuration)	DiseaseString = paste0(DiseaseString, "FullTrial")
	
	AgeDiseasePhaseBit = paste(AgeGroupNames[AgeGroup + 1], DiseaseString, PassiveOrNot_Char)
	if (country < 10)	CountryBitOfPlotTitle = paste("Country", country, CountryNames_Long[country+1], " ", AgeDiseasePhaseBit)    else 
	if (country == 10)	CountryBitOfPlotTitle = paste("CYD-14 "		, AgeDiseasePhaseBit)                                      else 
	if (country == 11)	CountryBitOfPlotTitle = paste("CYD-15 "		, AgeDiseasePhaseBit)                                      else 
	if (country == 12)	CountryBitOfPlotTitle = paste("CYD-14/15 "	, AgeDiseasePhaseBit)
	
	if (!OutputStringInTitle) ModelRun_string = "" else
	{
		if (ShortenOutputString) ModelRun_string = AbbreviateOutputString(OutputSubDir, modelrun = modelrun) else ModelRun_string = OutputSubDir
		ModelRun_string = paste0("\n", ModelRun_string)
	}
	if(ExpectedCasesInPlotTitle & Plot_Model) CasesString = paste0("\nE[Vac. Case] = ", round(ExpNoCasesAtEnd_Vaccine), "   E[Con Case] = ", round(ExpNoCasesAtEnd_Control)) 	else  CasesString = "" 
	if (!IncludePlotTitle)  PLOTTITLE = "" else		{	if (is.null(PLOTTITLE)) PLOTTITLE = paste0(CountryBitOfPlotTitle, WhichCurve_String, ImSubString, ModelRun_string, CasesString)		}
	
	### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
	### Plot everything. 
	### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
	
	if (SavePlot) png(file = FileName, res = PNG_res, units = "in", width = 5, height= 5)
	
	## new version with limits
	if (is.null(LOWERYLIM))
	{
		### Grow Vector to take minimum value of depending on what you're plotting. 
		VecToTakeMin = c()
		if (Plot_Data		)	VecToTakeMin = c(VecToTakeMin, Data_Haz_Group_1		, Data_Haz_Group_2		)
		if (Plot_Model		)	VecToTakeMin = c(VecToTakeMin, LowerCI_HazGroup_1	, LowerCI_HazGroup_2	) ### should only need LowerCI, and not mean or upper CI, unless something quite wrong. 
		if (Plot_MaxLike	) 	VecToTakeMin = c(VecToTakeMin, MaxLike_HazGroup_1	, MaxLike_HazGroup_2	)
		if (Plot_ModalPost	) 	VecToTakeMin = c(VecToTakeMin, ModalPost_HazGroup_1	, ModalPost_HazGroup_2	)
		
		LOWERYLIM = min(VecToTakeMin, na.rm = TRUE)
	}
	XLIM 		= c(0,MaxTi + 1)
	
	if (Plot_Model		) PlaceHolder = LowerCI_HazGroup_2 		else
	if (Plot_Data		) PlaceHolder = Data_Haz_Group_1 		else
	if (Plot_MaxLike	) PlaceHolder = MaxLike_HazGroup_1 		else
	if (Plot_ModalPost	) PlaceHolder = ModalPost_HazGroup_1	else 	stop("PlaceHolder if/else tree out of options")
	
	days 		= 1:length(PlaceHolder) - 1
	if (IncAxesLabels)	{	XLAB = "days post first dose" ; YLAB = "Proportion disease free"	} 	else 	{	XLAB = "" ; YLAB = ""	}
	
	plot(	days, PlaceHolder, ### PlaceHolder necessary as you may plot only data, or only mean values etc. 
			col = "white", ### col white so it doesn't affect transparency below
			xlab = XLAB, ylab = YLAB, 
			cex.lab = CEXLAB, cex.axis = CEXAXIS, main = PLOTTITLE, type = "l", cex.main = TitleCex,
			ylim = c(LOWERYLIM,1), xlim = XLIM,
			cex = 0.5)	
	
	if (Plot_Model) 
	{
		### Vaccine Mean and 95% CrI
		VacCrICols = col2rgb("pink", alpha = FALSE) ### colours to be set (pretty convoluted: Pick name of colour ("pink" here), get rgb values of that colour. Input rgb values into "col" argument, as then you can use a value of alpha for transparency.  
		polygon	(c(days,rev(days)),c(LowerCI_HazGroup_2,rev(UpperCI_HazGroup_2)), col=rgb(red = VacCrICols[1,]/255, green = VacCrICols[2,]/255, blue = VacCrICols[3,]/255, alpha = 0.5) ,	border = NA)
		points(days, Mean_HazGroup_2, col = "deeppink", type = "l", lty = 4, lwd = 8)
		
		ModalPost_Col_Vac = "deeppink4"
		if (Plot_MaxLike && Plot_ModalPost) MaxLike_Col_Vac = "orange" else MaxLike_Col_Vac = ModalPost_Col_Vac
		if (Plot_ModalPost) points(days, ModalPost_HazGroup_2	, col = ModalPost_Col_Vac)
		if (Plot_MaxLike) 	points(days, MaxLike_HazGroup_2		, col = MaxLike_Col_Vac)
		
		### Control Mean and 95% CrI
		ConCrICols = col2rgb("cadetblue1", alpha = FALSE) ### colours to be set.
		polygon	(c(days,rev(days)),c(LowerCI_HazGroup_1,rev(UpperCI_HazGroup_1)), col=rgb(red = ConCrICols[1,]/255, green = ConCrICols[2,]/255, blue = ConCrICols[3,]/255, alpha = 0.5) ,	border = NA)
		points(days, Mean_HazGroup_1, col = "skyblue", type = "l", lty = 4, lwd = 8)
		
		ModalPost_Col_Con = "skyblue3"
		if (Plot_MaxLike && Plot_ModalPost) MaxLike_Col_Con = "green" else MaxLike_Col_Con = ModalPost_Col_Con
		if (Plot_ModalPost) points(days, ModalPost_HazGroup_1	, col = ModalPost_Col_Con)
		if (Plot_MaxLike) 	points(days, MaxLike_HazGroup_1		, col = MaxLike_Col_Con)
	}
	
	if (Plot_Data)
	{
		### Vaccine Observed
		points(0:(length(Data_Haz_Group_2) - 1), Data_Haz_Group_2, col = "red", type = "l", lwd = LWD)
		### Control Observed
		points(0:(length(Data_Haz_Group_1) - 1), Data_Haz_Group_1, col = "blue", type = "l", lwd = LWD)
		
		if (CIsObs) ## plot polygons of observed
		{
			OBS_ConCrICols = col2rgb("blue", alpha = FALSE) ### colours to be set.
			OBS_VacCrICols = col2rgb("red", alpha = FALSE) ### colours to be set (pretty convoluted: Pick name of colour ("pink" here), get rgb values of that colour. Input rgb values into "col" argument, as then you can use a value of alpha for transparency.  
			
			days_Data = 0:(length(Data_Haz_Group_1) - 1)
			polygon	(c(days_Data,rev(days_Data)),c(CIs_Data_Haz_Group_1[,"Lower"],rev(CIs_Data_Haz_Group_1[,"Upper"])), col=rgb(red = OBS_ConCrICols[1,]/255, green = OBS_ConCrICols[2,]/255, blue = OBS_ConCrICols[3,]/255, alpha = 0.25) ,	border = NA)
			polygon	(c(days_Data,rev(days_Data)),c(CIs_Data_Haz_Group_2[,"Lower"],rev(CIs_Data_Haz_Group_2[,"Upper"])), col=rgb(red = OBS_VacCrICols[1,]/255, green = OBS_VacCrICols[2,]/255, blue = OBS_VacCrICols[3,]/255, alpha = 0.25) ,	border = NA)
		}
	}
	
	if (IncDataLegend)
		legend("bottomleft", c(LegLabels, "Data", "Model"), pch = rep(NA,4), col = c("red", "blue", "black", "black"), lty = c(rep(1,2), 1,4), lwd = c(rep(2,2), LWD, 8) , cex = LEGCEX	) else
		legend("bottomleft", LegLabels, pch = rep(NA,2), col = c("red", "blue"), lty = rep(1,2), lwd = rep(2,2) , cex = LEGCEX	)
	if (IncABlines)
	{
		DayLinesToCutXaxis = c(NoDaysActive + 1, NoDays_Active_And_Passive)
		abline (v = DayLinesToCutXaxis, col = "grey")
		if (Plot_Model)		abline (h = c(Mean_HazGroup_2[DayLinesToCutXaxis], Mean_HazGroup_1[DayLinesToCutXaxis]), col = rep(c("deeppink", "skyblue"), each = 2) )
	}
	
	if (SavePlot) dev.off()
	
}















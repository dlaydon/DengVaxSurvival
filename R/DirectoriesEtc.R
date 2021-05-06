
DefaultParamRangeFileName	<<- "ParamRanges_Default"
ChooseCountriesForAgeGroup 			= function(AgeGroup)
{
	if (AgeGroup == 0) Countries = 0:12
	if (AgeGroup == 1) Countries = c(0:4, 10)
	if (AgeGroup == 2) Countries = c(0:4, 10)
	if (AgeGroup == 3) Countries = c(0:4, 10)
	if (AgeGroup == 4) Countries = c(5:9, 11)
	if (AgeGroup == 5) Countries = c(5:9, 11)
	if (AgeGroup == 6) Countries = c(0:4, 10)
	if (AgeGroup == 7) Countries = c(0:4, 10)
	if (AgeGroup == 8) Countries = 0:12
	
	return(Countries)
}

DefineVariousQuantitiesForModelRun 	= function(modelrun = ModelRun)
{
	### Function makes global variables associated with particular model run
	
	if (modelrun$disease == "")					WhichDiseasePlots = c("")	### default, i.e. doesn't distinguish between mild and severe disease
	if (modelrun$disease == "_MILDSEVERE") 		WhichDiseasePlots = c("_Mild", "_Severe", "_Either")
	
	if (modelrun$CountriesFittedString == "") CountriesFitted = 0:9 else 
	{
		CountriesFitted = c()
		if (length(grep("0", modelrun$CountriesFittedString) == 1)) CountriesFitted = c(CountriesFitted, 0)
		if (length(grep("1", modelrun$CountriesFittedString) == 1)) CountriesFitted = c(CountriesFitted, 1)
		if (length(grep("2", modelrun$CountriesFittedString) == 1)) CountriesFitted = c(CountriesFitted, 2)
		if (length(grep("3", modelrun$CountriesFittedString) == 1)) CountriesFitted = c(CountriesFitted, 3)
		if (length(grep("4", modelrun$CountriesFittedString) == 1)) CountriesFitted = c(CountriesFitted, 4)
		if (length(grep("5", modelrun$CountriesFittedString) == 1)) CountriesFitted = c(CountriesFitted, 5)
		if (length(grep("6", modelrun$CountriesFittedString) == 1)) CountriesFitted = c(CountriesFitted, 6)
		if (length(grep("7", modelrun$CountriesFittedString) == 1)) CountriesFitted = c(CountriesFitted, 7)
		if (length(grep("8", modelrun$CountriesFittedString) == 1)) CountriesFitted = c(CountriesFitted, 8)
		if (length(grep("9", modelrun$CountriesFittedString) == 1)) CountriesFitted = c(CountriesFitted, 9)
	}
	
	if (modelrun$SS_Ks || modelrun$SS_VEs) N_STypes = 4 else N_STypes = 1
	if (modelrun$SS_VEs) N_STypes_VEs = 4 else N_STypes_VEs = 1
	if (modelrun$SS_Ks ) N_STypes_Ks  = 4 else N_STypes_Ks  = 1
	
	
	if (modelrun$SS_Ks || modelrun$SS_VEs) HRs_SeroNames = c("", paste0("_sero", 1:N_STypes)) 	else HRs_SeroNames = ""
	if (modelrun$SS_Ks || modelrun$SS_VEs) HRs_NumSTypes = 1 + N_STypes 						else HRs_NumSTypes = 1
	
	if (modelrun$phase == "_ACTIVE" 		& modelrun$disease == "") 		PhaseSeverityNames = c("Active")			 	else
	if (modelrun$phase == "_PASSIVE_ONLY" 	& modelrun$disease == "") 		PhaseSeverityNames = c("Passive")			 	else
	if (modelrun$disease == "" ) 											PhaseSeverityNames = c("Active" , "Passive") 	else 
	if (modelrun$disease == "_MILDSEVERE" 	& modelrun$hosp == "") 			PhaseSeverityNames = c("Mild" 	, "Severe") 	else 
	if (modelrun$disease == "_MILDSEVERE" 	& modelrun$hosp == "_hosp") 	PhaseSeverityNames = c("non-hosp", "hosp") 	
	
	BaselineSeroStatusNames = c("SeroNeg", "SeroPos", "Either_Serostatus")
	AgeGroupNames 			= c("All Ages", "2-5", "6-11", "12-14", "9-11", "12-16", "<9yrs", ">9yrs", ">9yrs")
	AgeGroupNamesVerbose 	= c("All Ages", "2-5 years", "6-11 years", "12-14 years", "9-11 years", "12-16 years", "Under 9 years", "Over 9 Years")
	
	
	if (modelrun$PSVEs == 1) EffPrefixStrings = c("AM_", "PS_") else EffPrefixStrings = ""
	
	
	#/*
	# 0		= any age group
	# 1,2,3	= CYD14 trial 2 to 5 years, 6 to 11 years, 12 to 14 years
	# 4,5	= CYD15 trial 9 to 11 years, 12 to 16 years.
	# 6,7	= CYD14 trial under 9 years, over 9 years
	# 8		= CYD14 and CYD15 trial over 9 years (no need for CYD15 under 9 years as all CYD15 participants are under 9 years, so this is just the same as CYD15 any age group).
	#*/
	
	if (modelrun$SFU)
	{
		if (modelrun$phase == "_ACTIVE"	)	KnotsPerCountry = 8;
		if (modelrun$phase == "_PASSIVE"	)	KnotsPerCountry = 11; 
		
	}	else
	{
		if (modelrun$phase == "_ACTIVE"		)	KnotsPerCountry = 6;
		if (modelrun$phase == "_PASSIVE"	)	KnotsPerCountry = 9; 
	}
	if (modelrun$phase == "_PASSIVE_ONLY"	)	KnotsPerCountry = 3 
	
	if (modelrun$phase == "_ACTIVE"			) WhichTrialPhases 			= "_ActiveOnly" 						else
	if (modelrun$phase == "_PASSIVE_ONLY"	) WhichTrialPhases 			= "" 									else 
	if (modelrun$phase == "_PASSIVE"		) WhichTrialPhases 			= c("", "_ActiveOnly", "_PassiveOnly") 
	
	SeroStringForData = modelrun$SerotypeString
	if (modelrun$SS_Ks || modelrun$SS_VEs) SeroStringForData = "_SeroSpec"
	
	
	#### redefine globally. 
	N_STypes 					<<- N_STypes
	N_STypes_VEs 				<<- N_STypes_VEs
	N_STypes_Ks 				<<- N_STypes_Ks
	HRs_SeroNames				<<- HRs_SeroNames
	HRs_NumSTypes				<<- HRs_NumSTypes
	CountriesFitted 			<<- CountriesFitted
	WhichDiseasePlots 			<<- WhichDiseasePlots
	WhichTrialPhases			<<- WhichTrialPhases
	PhaseSeverityNames 			<<- PhaseSeverityNames
	BaselineSeroStatusNames 	<<- BaselineSeroStatusNames
	AgeGroupNames			 	<<- AgeGroupNames
	AgeGroupNamesVerbose	 	<<- AgeGroupNamesVerbose
	KnotsPerCountry			 	<<- KnotsPerCountry
	SeroStringForData			<<- SeroStringForData
	EffPrefixStrings			<<- EffPrefixStrings
}
DefineFollowUpQuantities 			= function(modelrun = ModelRun)
{
	##### Below definitions not totally true. 
	##### If doing proper Sanofi_FollowUp_Definitions then there is no fixed number of active or passive surveillance days - it varies by patient. 
	##### However, your survival probabilities are aggregated across patients, so all you need here are reasonable averages from which to plot survival curves and draw boundaries. 
	##### Laurent said active phase is 25 months so input that here if doing Sanofi_FollowUp_Definitions
	
	if (modelrun$SFU	) NoDaysActive 		= round((25/12)*365) 	else 	NoDaysActive 		= 2*365	
	NoDays_Active_And_Passive = 3*365; 
	
	if (modelrun$Include_Late_Cases 	& 	modelrun$FakeExtObs 	& 	modelrun$SFU)
	{
		NoDaysActive 				= round((25/12)*365)
		NoDays_Active_And_Passive 	= 1372; 
	}
	
	if (modelrun$phase == "_ACTIVE") NoDays_Active_And_Passive = NoDaysActive ### i.e. redefine if you're not doing passive phase. 
	
	NoDaysPassive				= NoDays_Active_And_Passive - NoDaysActive
	PassivePhaseDuration_Years 	= NoDaysPassive	/ 365
	
	NoDaysActive 				<<- NoDaysActive
	NoDaysPassive 				<<- NoDaysPassive
	NoDays_Active_And_Passive 	<<- NoDays_Active_And_Passive
	PassivePhaseDuration_Years 	<<- PassivePhaseDuration_Years
}
ChooseOutputString 					= function(ModelRun, Folder = FALSE, PrintToConsole = TRUE, AddExtraOutputToDirectories = TRUE, StringPrefix = "")
{
	OutputSubDir = StringPrefix
	
	OutputSubDir = paste0(OutputSubDir, ModelRun$ModelVariant)
	OutputSubDir = paste0(OutputSubDir, ModelRun$phase)
	if (ModelRun$EffNegWane == "FROM_ZERO"	)			OutputSubDir = paste0(OutputSubDir, "_AWFZ"	) else  #### Alternative waning. efficaies that start out negative don't "wane" to higher values of less absolute magnitude. They "wane from zero" to their negative value. 
	if (ModelRun$EffNegWane == "NO_WANE" 	)			OutputSubDir = paste0(OutputSubDir, "_nENW"	) 		#### stands for no eff neg waning. See notes in EffNegWane_Option definition. 
	
	OutputSubDir = paste0(OutputSubDir, ModelRun$disease)
	
	### For disease == "", don't want CYD14 or CYD15 runs having hosp on them (as they make no difference) but the name complicates other scripts. Get rid of here. 
	OutputSubDir = paste0(OutputSubDir, ModelRun$hosp) ### 1st way of doing it
	
	if ( (ModelRun$phase == "_PASSIVE" || ModelRun$disease == "_MILDSEVERE") && ModelRun$PS_Ks_Multiply_AM_Ks) OutputSubDir = paste0(OutputSubDir, "X")
	
	
	if (!ModelRun$SS_VEs & !ModelRun$SS_Ks)	OutputSubDir = paste0(OutputSubDir, ModelRun$SerotypeString)
	
	
	if (ModelRun$SS_VEs && ModelRun$SS_Ks)	OutputSubDir = paste0(OutputSubDir, "_SSVEs" ) else ## see note below as to why you've lopped off the final "s". 
	if (ModelRun$SS_VEs)					OutputSubDir = paste0(OutputSubDir, "_SS_VEs")
	if (ModelRun$SS_Ks)
	{
		if (!ModelRun$SS_VEs) OutputSubDir = paste0(OutputSubDir, "_SS_")
		#### A little convulted. Default for SSKs (SSVEsKs) is to have output string appended by either SS_Ks (SSVEsKs). Howevever if only some Ks are serospecific (e.g. KAM_0 and KPS_1), you want it to read _SSKAM_0PS_1. 
		#### However you still want the "s" at the end. So that happens after this if statement (regardless of whether SSVEs or not). Check this as will matter for output storage. 
		AnySSKsFalse = FALSE; 
		if (!ModelRun$SS_KAM_0)  AnySSKsFalse = TRUE
		if (!ModelRun$SS_KAM_2)  AnySSKsFalse = TRUE
		if (!ModelRun$SS_KPS_0)  AnySSKsFalse = TRUE
		if (!ModelRun$SS_KPS_1)  AnySSKsFalse = TRUE
		if (!ModelRun$SS_KPS_2)  AnySSKsFalse = TRUE
		
		if (AnySSKsFalse)
		{
			if (ModelRun$SS_KAM_0)  OutputSubDir = paste0(OutputSubDir, "KAM0s")
			if (ModelRun$SS_KAM_2)  OutputSubDir = paste0(OutputSubDir, "KAM2s")
			if (ModelRun$SS_KPS_0)  OutputSubDir = paste0(OutputSubDir, "KPS0s")
			if (ModelRun$SS_KPS_1)  OutputSubDir = paste0(OutputSubDir, "KPS1s")
			if (ModelRun$SS_KPS_2)  OutputSubDir = paste0(OutputSubDir, "KPS2s")
		}
		else				OutputSubDir = paste0(OutputSubDir, "Ks")
	}
	
	if (ModelRun$SType_Equiv & (ModelRun$SS_VEs | ModelRun$SS_Ks))
											OutputSubDir = paste0(OutputSubDir, "Equiv"	) 
	
	if (ModelRun$dose != "") 				OutputSubDir = paste0(OutputSubDir, "_SINGLE_DOSE")
	
	OutputSubDir = paste0(OutputSubDir, ModelRun$CountriesFittedString)
	
	if (ModelRun$LinKnts)							OutputSubDir = paste0(OutputSubDir, "_LinKnts")
	if (ModelRun$SimulatedAnnealing)				OutputSubDir = paste0(OutputSubDir, "_SIM_ANNEALING")
	if (ModelRun$FitWaningRate)						OutputSubDir = paste0(OutputSubDir, "_FLIP_WANING")
	if (ModelRun$MH_Aug_OR_GIBBS_AUG == "_MH_AUG")	OutputSubDir = paste0(OutputSubDir, "_MH_AUG")
	if (!ModelRun$AreWeAugmenting)					OutputSubDir = paste0(OutputSubDir, "_NotAug") 
	if (ModelRun$SingleEff) 						OutputSubDir = paste0(OutputSubDir, "_SingleEff")
	if (ModelRun$HillWaning)						OutputSubDir = paste0(OutputSubDir, "_HillWaning")
	if (ModelRun$ResidEffs)  						OutputSubDir = paste0(OutputSubDir, "_ResidEffs")
	
	if (ModelRun$Weighting_Pass_Sev)
	{
		if (as.numeric(ModelRun$PS_Weight) > 1) 	OutputSubDir = paste0(OutputSubDir, "_PSweight_" , as.character(ModelRun$PS_Weight)) else 
													OutputSubDir = paste0(OutputSubDir, "_PSweight_" , as.character(paste0(as.character(as.integer(1/ModelRun$PS_Weight)), "inv")))
	} 
	if (ModelRun$Fixed_Severe_RelRisks)	
	{
		if (ModelRun$FSKs_Ratio_SetNum == 0) OutputSubDir = paste0(OutputSubDir, "_Fixed_Ks") else
		if (ModelRun$FSKs_Ratio_SetNum == 1) OutputSubDir = paste0(OutputSubDir, "_FSKsV2"	) else
		if (ModelRun$FSKs_Ratio_SetNum == 2) OutputSubDir = paste0(OutputSubDir, "_FSKs3"	) else
		stop("ChooseOutputString Error: ModelRun$FSKs_Ratio_SetNum value not recognized")
	}			
	
	if (ModelRun$ASVE != "INDEPENDENT")
	{
		if (ModelRun$ASVE == "HILL"			) OutputSubDir = paste0(OutputSubDir,  "_ASVE"		)	else	####  "age-specific vaccine efficacy" 
		if (ModelRun$ASVE == "CATEGORICAL"	) OutputSubDir = paste0(OutputSubDir,  "_AGSVE"		)	else	####  "age-group-specific vaccine efficacy"
		if (ModelRun$ASVE == "SPLINE"		) OutputSubDir = paste0(OutputSubDir,  "_ASVESpln"	)	else	####  "age-specific vaccine efficacy with Spline"
		if (ModelRun$ASVE == "SPLINE_LINE"	) OutputSubDir = paste0(OutputSubDir,  "_ASVELine"	)	else	####  "age-specific vaccine efficacy with Spline"
		if (ModelRun$ASVE == "SPLINE_STEP"	) OutputSubDir = paste0(OutputSubDir,  "_ASVEStep"	)	else	####  "age-specific vaccine efficacy with Spline"
		if (ModelRun$ASVE == "CUBIC"		) OutputSubDir = paste0(OutputSubDir,  "_ASVECubic"	)	else	####  "age-specific vaccine efficacy with Spline"
			stop("ChooseOutputString Error: ModelRun$ASVE not recognized")
		
		if (ModelRun$ASVE == "SPLINE" || ModelRun$ASVE == "SPLINE_LINE" || ModelRun$ASVE == "SPLINE_STEP" || ModelRun$ASVE == "CUBIC")
			if (ModelRun$ASVE_AdditionalKnots	) OutputSubDir = paste0(OutputSubDir,  "AddKnots")
		
		if (!ModelRun$AS_VE_Homogeneous) 									OutputSubDir = paste0(OutputSubDir, "hetero") 
		if (ModelRun$SSASVE_Additive & ModelRun$SS_VEs)						OutputSubDir = paste0(OutputSubDir, "Add") 
		if ( ModelRun$ASVE_OnlyOneSeroStatus == 1 && ModelRun$ASVE_BS == 0) OutputSubDir = paste0(OutputSubDir, "SNeg")
		if ( ModelRun$ASVE_OnlyOneSeroStatus == 1 && ModelRun$ASVE_BS == 1) OutputSubDir = paste0(OutputSubDir, "SPos")
	}

	if (ModelRun$AS_Haz != "INDEPENDENT")
	{
		if (ModelRun$AS_Haz == "HILL"			) 	OutputSubDir = paste0(OutputSubDir, "_AS_Haz"		)	else 
		if (ModelRun$AS_Haz == "CATEGORICAL"	)	OutputSubDir = paste0(OutputSubDir, "_AS_Hazmult"	) 	else
		if (ModelRun$AS_Haz == "SPLINE"			)	OutputSubDir = paste0(OutputSubDir, "_AS_HazSpln"	) 	else
		if (ModelRun$AS_Haz == "SPLINE_LINE"	)	OutputSubDir = paste0(OutputSubDir, "_AS_HazLine"	) 	else
		if (ModelRun$AS_Haz == "SPLINE_STEP"	)	OutputSubDir = paste0(OutputSubDir, "_AS_HazStep"	) 	else
		if (ModelRun$AS_Haz == "CUBIC"			)	OutputSubDir = paste0(OutputSubDir, "_AS_HazCubic"	) 	else
			stop("ChooseOutputString Error: ModelRun$AS_Haz not recognized")
	}
	
	if (ModelRun$AS_Haz == "SPLINE" || ModelRun$AS_Haz == "SPLINE_LINE" || ModelRun$AS_Haz == "SPLINE_STEP" || ModelRun$AS_Haz == "CUBIC")
		if (ModelRun$AS_Haz_AdditionalKnots	) OutputSubDir = paste0(OutputSubDir,  "AddKnots")
	
	if (ModelRun$AS_Waning != "INDEPENDENT")
	{
		OutputSubDir = paste0(OutputSubDir, "_AS_Wane")
		if (ModelRun$AS_Waning_OnlyOneSeroStatus)
		{
			if (ModelRun$AS_Waning_OneSeroBS == 0) OutputSubDir = paste0(OutputSubDir, "SNeg") else 
			if (ModelRun$AS_Waning_OneSeroBS == 1) OutputSubDir = paste0(OutputSubDir, "SPos")
		}
	}
	if (ModelRun$AS_Waning == "CATEGORICAL"	) OutputSubDir = paste0(OutputSubDir, "Cat"	) else
	if (ModelRun$AS_Waning == "SPLINE" ||ModelRun$AS_Waning == "SPLINE_STEP" ||ModelRun$AS_Waning == "SPLINE_LINE" ||ModelRun$AS_Waning == "CUBIC") 
	{
		if (ModelRun$AS_Waning == "SPLINE"		)	OutputSubDir = paste0(OutputSubDir, "Spln")
		if (ModelRun$AS_Waning == "SPLINE_LINE"	) 	OutputSubDir = paste0(OutputSubDir, "Line")
		if (ModelRun$AS_Waning == "SPLINE_STEP"	) 	OutputSubDir = paste0(OutputSubDir, "Step")
		if (ModelRun$AS_Waning == "CUBIC"		) 	OutputSubDir = paste0(OutputSubDir, "Cubic")
		if (ModelRun$AS_Waning_KnotSet!= 0) 	
		{
			OutputSubDir = paste0(OutputSubDir, "AddKnots")
			if (ModelRun$AS_Waning_KnotSet != 1) OutputSubDir = paste0(OutputSubDir, ModelRun$AS_Waning_KnotSet)
		}
	}
	if (ModelRun$AS_Waning != "INDEPENDENT")
		if (ModelRun$AS_Waning_Homogeneous)
			OutputSubDir = paste0(OutputSubDir, "Homog")
	
	
	if (ModelRun$PSVEs)									OutputSubDir = paste0(OutputSubDir, "_PSVEs")
	if (ModelRun$AdjHaz)
	{
		if (ModelRun$AdjHazFitted)		OutputSubDir = paste0(OutputSubDir, "_fAdjHaz") else OutputSubDir = paste0(OutputSubDir, "_AdjHaz")
	}
	if (ModelRun$PooledCountries)						OutputSubDir = paste0(OutputSubDir 	, "_PooledCountries"	)	else
	if (ModelRun$PooledTrials)							OutputSubDir = paste0(OutputSubDir 	, "_PooledTrials"		)
	if (ModelRun$Empirical_SeroPrevs) 					OutputSubDir = paste0(OutputSubDir	, "_Empirical_SeroPrevs")
	
	if (ModelRun$AllDosesRequired_SNeg & ModelRun$AllDosesRequired_SPos) OutputSubDir = paste0(OutputSubDir, "_BothSeroNeed3Dose") else
	{
		if (ModelRun$AllDosesRequired_SNeg)							OutputSubDir = paste0(OutputSubDir, "_SNegNeed3Dose") else
		{
			AnySNegNeed3 	= FALSE
			AG_SNeg_String	= ""
			if (ModelRun$AllDosesRequired_SNeg_AgeGroup1)		{ 	AnySNegNeed3 = TRUE		; AG_SNeg_String = paste0(AG_SNeg_String, "1");	} 
			if (ModelRun$AllDosesRequired_SNeg_AgeGroup2)		{ 	AnySNegNeed3 = TRUE		; AG_SNeg_String = paste0(AG_SNeg_String, "2");	} 
			if (ModelRun$AllDosesRequired_SNeg_AgeGroup3)		{ 	AnySNegNeed3 = TRUE		; AG_SNeg_String = paste0(AG_SNeg_String, "3");	} 
			if (ModelRun$AllDosesRequired_SNeg_AgeGroup4)		{ 	AnySNegNeed3 = TRUE		; AG_SNeg_String = paste0(AG_SNeg_String, "4");	} 
			if (ModelRun$AllDosesRequired_SNeg_AgeGroup5)		{ 	AnySNegNeed3 = TRUE		; AG_SNeg_String = paste0(AG_SNeg_String, "5");	} 
			
			if (AnySNegNeed3) OutputSubDir = paste0(OutputSubDir, "_SNegAG", AG_SNeg_String, "_Need3Dose") 
		}
		if (ModelRun$AllDosesRequired_SPos)							OutputSubDir = paste0(OutputSubDir, "_SPosNeed3Dose") else
		{
			AnySPosNeed3 	= FALSE
			AG_SPos_String	= ""
			if (ModelRun$AllDosesRequired_SPos_AgeGroup1)		{ 	AnySPosNeed3 = TRUE		; AG_SPos_String = paste0(AG_SPos_String, "1");	} 
			if (ModelRun$AllDosesRequired_SPos_AgeGroup2)		{ 	AnySPosNeed3 = TRUE		; AG_SPos_String = paste0(AG_SPos_String, "2");	} 
			if (ModelRun$AllDosesRequired_SPos_AgeGroup3)		{ 	AnySPosNeed3 = TRUE		; AG_SPos_String = paste0(AG_SPos_String, "3");	} 
			if (ModelRun$AllDosesRequired_SPos_AgeGroup4)		{ 	AnySPosNeed3 = TRUE		; AG_SPos_String = paste0(AG_SPos_String, "4");	} 
			if (ModelRun$AllDosesRequired_SPos_AgeGroup5)		{ 	AnySPosNeed3 = TRUE		; AG_SPos_String = paste0(AG_SPos_String, "5");	} 
			
			if (AnySPosNeed3) OutputSubDir = paste0(OutputSubDir, "_SPosAG", AG_SPos_String, "_Need3Dose") 
		}
	}
	if (ModelRun$ExtImSub == "AS_DATA") OutputSubDir = paste0(OutputSubDir, "_ExtImData") else
	if (ModelRun$ExtImSub == "AS_PROB") OutputSubDir = paste0(OutputSubDir, "_ExtImProb") 
		
	
	if (ModelRun$ParamRangeFileName != DefaultParamRangeFileName) OutputSubDir = paste0(OutputSubDir, "_", ModelRun$ParamRangeFileName)
	
	
	if (ModelRun$SFU)  OutputSubDir	= paste0(OutputSubDir, "_SFU")  #### this deals with old and new versions where you had SFU defined by strings, and new version where SFU defined by 0 or 1. 
	
	
	if (ModelRun$Include_Late_Cases) 
	{
		OutputSubDir = paste0(OutputSubDir, "_IncLate")
		
		if (ModelRun$FakeExtObs) OutputSubDir = paste0(OutputSubDir, "_FakeExtObs")
	}
	
	if (AddExtraOutputToDirectories) #### here you 	DO 		use the name of the output as the directory (e.g. for NoWaning where you want to separate this run from runs where Waning rate fitted) 
	{
		OutputSubDir 	= paste0(OutputSubDir, ModelRun$OutputStringExtra)
		if (!Folder)	OutputSubDir 				= paste0("_", OutputSubDir)
		
	} else							#### here you 	DON'T 	use the name of the output as the directory (e.g. for _LongChain where you want a particular run to be "canonical") 
	{
		if (!Folder) OutputSubDir 				= paste0("_", OutputSubDir, ModelRun$OutputStringExtra)
	}
	
	if (PrintToConsole)
	{
		if (Folder) 	cat(paste("OutputSubDir", OutputSubDir, "\n"))
		if (!Folder) 	cat(paste("OutputString", OutputSubDir, "\n"))
		
	}
	
	
	return(OutputSubDir)
	#OutputString				<<- OutputString
	#OutputSubDir	<<- OutputSubDir
}
AbbreviateOutputString 				= function(OutputSubDir_Or_OutputString, modelrun = ModelRun)
{
	ShortenedString = OutputSubDir_Or_OutputString
	
	#### remove all the things that are (effectively) defaults 
	if (modelrun$ModelVariant == "VAC_SILENT"	) ShortenedString = sub(modelrun$ModelVariant	, "", ShortenedString)
	if (modelrun$phase == "_PASSIVE"			) ShortenedString = sub(modelrun$phase			, "", ShortenedString)
	if (modelrun$SFU							) ShortenedString = sub("_SFU"					, "", ShortenedString)
	
	#### Abbreviate 
	ShortenedString = sub("_AGSVEhetero_AS_Hazmult_AS_WaneCat"				, "_AgeCat"		, ShortenedString)
	ShortenedString = sub("_ASVECubichetero_AS_HazCubic_AS_WaneCubic"		, "_AgeCubic"	, ShortenedString)
	ShortenedString = sub("_ASVESplnhetero_AS_HazSpln_AS_WaneSpln"			, "_AgeSpline"	, ShortenedString)
	ShortenedString = sub("_AGSVEheteroAdd_AS_Hazmult_AS_WaneCat"			, "_AgeCatAdd"		, ShortenedString)
	ShortenedString = sub("_ASVECubicheteroAdd_AS_HazCubic_AS_WaneCubic"	, "_AgeCubicAdd"	, ShortenedString)
	ShortenedString = sub("_ASVESplnheteroAdd_AS_HazSpln_AS_WaneSpln"		, "_AgeSplineAdd"	, ShortenedString)
	ShortenedString = sub("PooledCountries", "PC", ShortenedString)
	ShortenedString = sub("SS_VEs", "SSVEs", ShortenedString)
	ShortenedString = sub("NO_AUGMENTATION", "NoAug", ShortenedString)
	ShortenedString = sub("Fixed_Ks", "FSKs", ShortenedString)
	ShortenedString = sub("hetero", "", ShortenedString)
	ShortenedString = sub("_Cs01234", "_CYD14", ShortenedString)
	ShortenedString = sub("_Cs56789", "_CYD15", ShortenedString)
	ShortenedString = sub("_AS_Hazmult", "_ASHazCat", ShortenedString)
	ShortenedString = sub("AGSVEhetero_AS_Hazmult_AS_WaneCat", "_AgeCat", ShortenedString)
	ShortenedString = sub("_AS_Haz", "_ASHaz", ShortenedString)
	while (substr(ShortenedString, 1, 1) == "_") ShortenedString = sub("_", "", ShortenedString)
	
	return(ShortenedString)
}


DefineModelRuns = function( 
		
		RemoveDuplicates 			= TRUE											, 
		AgeEffectsSame				= FALSE											, 
		AgeEffectsSame_FOI_VEs		= TRUE											, #### only relevant if AgeEffectsSame = FALSE. For debugging runs where you want Age specific hazards (AS_Haz) and Age specific efficacy / transient immunity, but not Age specific waning. 
		PrintToConsole				= FALSE											,
		FromDateAndTime				= NULL											,
		
		ModelVariants 								= c("VAC_SILENT")				,
		DataFilenames								= c("")							,
		SerotypeStrings 							= c("")          				,
		CountriesFittedStrings 						= c("")          				,
		Doses 										= c("")          				,
		Phases 										= c("_PASSIVE")  				,
		EffNegWane_Options							= "NO_WANE"						, 
		MildSeveres 								= c("")          				,
		PS_Ks_Multiply_AM_Ks_AND_NOT				= 1								, 
		Hosp_And_NonHosp 							= c("")          				,
		LogAndNormalScales 							= c("")          				,
		LinKntsVec									= 0								, 
		Normal_AND_SimulatedAnnealing 				= 0              				,
		AreWeAugmenting_And_Not 					= 1              				,
		SFU_and_not_SFU 							= 1              				,
		SS_VEs_AND_Not 								= 1              				,
		SS_Ks_AND_Not 								= 0              				,
		
		SS_KAM_0_AND_NOT   							= 1              				, ### should be set to 1 as that is default, and is only relevant if SS_Ks
		SS_KAM_2_AND_NOT   							= 1              				, ### should be set to 1 as that is default, and is only relevant if SS_Ks
		SS_KPS_0_AND_NOT   							= 1              				, ### should be set to 1 as that is default, and is only relevant if SS_Ks
		SS_KPS_1_AND_NOT   							= 1              				, ### should be set to 1 as that is default, and is only relevant if SS_Ks
		SS_KPS_2_AND_NOT   							= 1              				, ### should be set to 1 as that is default, and is only relevant if SS_Ks
		
		SType_Equiv_And_Not							= 0								,
		LikeWeightingAndNot 						= 0              				,
		PS_Weights 									= c(0.1)       					,
		MH_Aug_AND_GIBBS_AUG 						= c("")          				,
		Fixed_Severe_RelRisks_AND_not 				= 1								, ### default used to be 0 but you rarely use that now 
		FSKs_Ratio_SetNumbers 						= 2								, 
		ASVE_Options 								= "CATEGORICAL"					,
		AS_VE_Homogeneous_And_Not					= 0								,
		SSASVE_Additive_And_Not						= 1								, 
		ASVE_OnlyOneSeroStatus_And_Not				= 0								,
		ASVE_BSs									= 0								,
		ASVE_AdditionalKnots_And_Not				= 0								, # 0 means FALSE here. not that the default number of knots is zero. 
		AS_Haz_Options 								= "CATEGORICAL"					,
		AS_Haz_AdditionalKnots_And_Not				= 0								, # 0 means FALSE here. not that the default number of knots is zero. 
		AS_Waning_Options							= "INDEPENDENT"					,
		AS_Waning_OnlyOneSeroStatus_And_Not			= 0								,
		AS_Waning_OneSeroBSs						= 0								,
		AS_Waning_KnotSets							= 0								, # 0 means FALSE here. not that the default number of knots is zero. 
		AS_Waning_Homogeneous_And_Not				= 0								, 
		FitWaningRate_AND_Not 						= 0								,
		HillWaning_And_Not 							= 0								,
		Empirical_SeroPrevs_And_Not 				= 0								,
		PSVEs_And_Not 								= 0								,
		AdjHazAndNot		 						= 1								,
		AdjHazFitted_And_Not						= 1								, 
		Include_Late_Cases	 						= 0								,
		Fake_Extended_Observation_And_Not 			= 0								,
		ExtImSub_Options							= "IGNORED"						,
		PooledCountries_And_Not 					= 0								,
		PooledTrials_And_Not 						= 0								,
		PassiveOnly_And_Not 						= 0								,
		ResidEffs_And_Not 							= 0								,
		SingleEff_And_Not 							= 0								,
		AllDosesRequired_SNeg_And_Not 				= 0								,
		AllDosesRequired_SNeg_AgeGroup1_And_Not		= 0								, 
		AllDosesRequired_SNeg_AgeGroup2_And_Not		= 0								, 
		AllDosesRequired_SNeg_AgeGroup3_And_Not		= 0								, 
		AllDosesRequired_SNeg_AgeGroup4_And_Not		= 0								, 
		AllDosesRequired_SNeg_AgeGroup5_And_Not		= 0								, 
		AllDosesRequired_SPos_And_Not 				= 0								,
		AllDosesRequired_SPos_AgeGroup1_And_Not		= 0								, 
		AllDosesRequired_SPos_AgeGroup2_And_Not		= 0								, 
		AllDosesRequired_SPos_AgeGroup3_And_Not		= 0								, 
		AllDosesRequired_SPos_AgeGroup4_And_Not		= 0								, 
		AllDosesRequired_SPos_AgeGroup5_And_Not		= 0								, 
		ParamRangeFileNames							= "prs1_2"						, 
		OutputStringExtras							= ""							, 
		OutputStringExtras_PrevChains           	= ""   
)
{
	
	#
	
	### Make all possible versions of ModelRuns - note that expand.grid will not account for duplicates or for model runs that have not actually been run. 
	ModelRuns = expand.grid		(
			
			ModelVariant 						= ModelVariants								,
			DataFilename						= DataFilenames								,
			
			disease 							= MildSeveres								,
			PS_Ks_Multiply_AM_Ks				= PS_Ks_Multiply_AM_Ks_AND_NOT				,
			hosp 								= Hosp_And_NonHosp							, 
			phase 								= Phases									, 
			EffNegWane							= EffNegWane_Options							,
			LinKnts		 						= LinKntsVec								,
			SerotypeString 						= SerotypeStrings							,
			SFU 								= SFU_and_not_SFU							,
			SimulatedAnnealing 					= Normal_AND_SimulatedAnnealing				,
			CountriesFittedString 				= CountriesFittedStrings					,
			dose 								= Doses										,
			FitWaningRate 						= FitWaningRate_AND_Not						,
			MH_Aug_OR_GIBBS_AUG 				= MH_Aug_AND_GIBBS_AUG						,
			SS_VEs 								= SS_VEs_AND_Not							,
			SS_Ks 								= SS_Ks_AND_Not								,
			
			SS_KAM_0   							= SS_KAM_0_AND_NOT                          ,
			SS_KAM_2   							= SS_KAM_2_AND_NOT                          ,
			SS_KPS_0   							= SS_KPS_0_AND_NOT                          ,
			SS_KPS_1   							= SS_KPS_1_AND_NOT                          ,
			SS_KPS_2   							= SS_KPS_2_AND_NOT                          ,
			
			SType_Equiv							= SType_Equiv_And_Not						, 
			HillWaning 							= HillWaning_And_Not						,
			Weighting_Pass_Sev 					= LikeWeightingAndNot						,	
			PS_Weight 							= PS_Weights								,
			ASVE 								= ASVE_Options								,
			SSASVE_Additive						= SSASVE_Additive_And_Not					,
			AS_Haz 								= AS_Haz_Options							,
			AS_Haz_AdditionalKnots				= AS_Haz_AdditionalKnots_And_Not			, # 0 means FALSE here. not that the default number of knots is zero. 
			AS_VE_Homogeneous		 			= AS_VE_Homogeneous_And_Not 				,  
			ASVE_OnlyOneSeroStatus 				= ASVE_OnlyOneSeroStatus_And_Not 			, 
			ASVE_BS 							= ASVE_BSs									,
			ASVE_AdditionalKnots				= ASVE_AdditionalKnots_And_Not				, # 0 means FALSE here. not that the default number of knots is zero. 
			AS_Waning							= AS_Waning_Options							,
			AS_Waning_OnlyOneSeroStatus			= AS_Waning_OnlyOneSeroStatus_And_Not		,
			AS_Waning_OneSeroBS					= AS_Waning_OneSeroBSs						,
			AS_Waning_KnotSet		 			= AS_Waning_KnotSets						,
			AS_Waning_Homogeneous				= AS_Waning_Homogeneous_And_Not				,
			PSVEs 								= PSVEs_And_Not								,
			AdjHaz 								= AdjHazAndNot								, 
			AdjHazFitted						= AdjHazFitted_And_Not						,
			Fixed_Severe_RelRisks				= Fixed_Severe_RelRisks_AND_not				,
			FSKs_Ratio_SetNum					= FSKs_Ratio_SetNumbers						, 
			Include_Late_Cases					= Include_Late_Cases						,
			Empirical_SeroPrevs					= Empirical_SeroPrevs_And_Not				,
			FakeExtObs							= Fake_Extended_Observation_And_Not			,
			PooledTrials						= PooledTrials_And_Not						,
			ResidEffs							= ResidEffs_And_Not							,
			AllDosesRequired_SNeg 				= AllDosesRequired_SNeg_And_Not				, 
			AllDosesRequired_SNeg_AgeGroup1		= AllDosesRequired_SNeg_AgeGroup1_And_Not	,
			AllDosesRequired_SNeg_AgeGroup2		= AllDosesRequired_SNeg_AgeGroup2_And_Not	,
			AllDosesRequired_SNeg_AgeGroup3		= AllDosesRequired_SNeg_AgeGroup3_And_Not	,
			AllDosesRequired_SNeg_AgeGroup4		= AllDosesRequired_SNeg_AgeGroup4_And_Not	,
			AllDosesRequired_SNeg_AgeGroup5		= AllDosesRequired_SNeg_AgeGroup5_And_Not	,
			AllDosesRequired_SPos				= AllDosesRequired_SPos_And_Not				,
			AllDosesRequired_SPos_AgeGroup1		= AllDosesRequired_SPos_AgeGroup1_And_Not	,
			AllDosesRequired_SPos_AgeGroup2		= AllDosesRequired_SPos_AgeGroup2_And_Not	,
			AllDosesRequired_SPos_AgeGroup3		= AllDosesRequired_SPos_AgeGroup3_And_Not	,
			AllDosesRequired_SPos_AgeGroup4		= AllDosesRequired_SPos_AgeGroup4_And_Not	,
			AllDosesRequired_SPos_AgeGroup5		= AllDosesRequired_SPos_AgeGroup5_And_Not	,
			SingleEff							= SingleEff_And_Not							,
			AreWeAugmenting						= AreWeAugmenting_And_Not					,
			ExtImSub 							= ExtImSub_Options							,
			ParamRangeFileName					= ParamRangeFileNames						,
			PooledCountries 					= PooledCountries_And_Not					, 
			OutputStringExtra					= OutputStringExtras						,
			OutputStringExtra_PrevChain			= OutputStringExtras_PrevChains				, 
			
			stringsAsFactors = FALSE
	)
	
	#### Create OutputStrings for all possible model runs...
	OutputFolderNames 	= rep(NA, dim(ModelRuns)[1])
	for (run in 1:dim(ModelRuns)[1])	OutputFolderNames[run] = ChooseOutputString(ModelRuns[run,], PrintToConsole = PrintToConsole, Folder = TRUE)
	#### .... bind them to data.frame
	ModelRuns = cbind(ModelRuns, OutputFolderNames)
	
	#### Remove duplicates
	if (RemoveDuplicates)			ModelRuns = ModelRuns[!duplicated(ModelRuns$OutputFolderNames), ]
	#### Remove runs where age effects modelled differently (e.g. categorical AS_Haz, Hill ASVE and spline AS_Wane). 
	if (AgeEffectsSame) 			
		ModelRuns = ModelRuns[which(ModelRuns$ASVE == ModelRuns$AS_Haz & ModelRuns$ASVE == ModelRuns$AS_Waning), ] else ### and therefore AS_Haz == AS_Waning
	if (AgeEffectsSame_FOI_VEs)
		ModelRuns = ModelRuns[which(ModelRuns$ASVE == ModelRuns$AS_Haz & 												#### But not necessarily ASVE == AS_Waning or AS_Haz == AS_Waning.....
								(ModelRuns$ASVE == ModelRuns$AS_Waning | ModelRuns$AS_Waning == "INDEPENDENT")), ] 	#### However, if ASVE != AS_Waning, can only be because AS_Waning == "INDEPENDENT" (i.e. don't want mixing of splines with categoricals with cubics etc.)					 
	
	if (dim(ModelRuns)[1] == 0) stop("DefineModelRuns Error: ModelRuns data.frame empty. Is NewRunsOnly on unintentionally?")
	
	return(ModelRuns)
}




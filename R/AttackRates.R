Choose_AR_TableChar = function(WhichDisease = "", Phase = "Both", ImSub = FALSE, RootAttackRateTable_Char = "AttackRateTable", WBIC = FALSE)
{
	ImSubString 			= GetImSubString(ImSub)
	WBIC_String				= Choose_WBIC_String(WBIC)
	AttackRateTable_Char 	= paste0(WBIC_String, RootAttackRateTable_Char)
	if (any( c("Active", "Passive") == Phase)) AttackRateTable_Char = paste0(AttackRateTable_Char, "_", Phase, "Only")
	AttackRateTable_Char = paste0(AttackRateTable_Char, WhichDisease, ImSubString)
	return(AttackRateTable_Char)
}
Get_AR_Samples = function(Country = country, AgeGroup = 0, SeroStatus = "Both_Sero", Arm = "Either", WhichDisease = "", Phase = "Both", ImSub = FALSE, stratum = NULL)
{
	if (is.null(stratum)) 
		stratum 			= ChooseStratum(AgeGroup = AgeGroup, Country = country, Arm = Arm, SeroStatus = SeroStatus)
	AttackRateTable_Char 	= Choose_AR_TableChar(WhichDisease = WhichDisease, Phase = Phase, ImSub = ImSub, RootAttackRateTable_Char = "AttackRateTable")
	Rates 					= get(AttackRateTable_Char)[stratum, ]
	return(Rates)
}

ObservedAttackRate 	= function(country = country, AgeGroup = 0, SeroStatus = "Both_Sero", Arm = "Either", 
		WhichDisease = whichdiseaseplot, Phase = "Both", modelrun = ModelRun, Units = "years", percent = TRUE)
{
	if (modelrun$disease == "_MILDSEVERE")
	{
		if (Phase == "Passive" && WhichDisease == "_Mild"	) stop("ObservedAttackRate error: Phase == Passive && WhichDisease == _Mild"	)
		if (Phase == "Passive" && WhichDisease == "_Either"	) stop("ObservedAttackRate error: Phase == Passive && WhichDisease == _Either"	)
	}
	
	stratum 			= ChooseStratum(AgeGroup, country, Arm, SeroStatus)
	
	### Choose PopSizes 
	## Using ObservedPopulationSizes won't work either for passive phase (as you need to exclude Active Phase cases or people LTFU), or for a particular serostatus (i.e. not Both_Sero), as in the latter case you don't/can't record them. If however you want to calculate attack rates in immunosubset then you need serostatus specific attack rates. 
	
	CaseColName 		= ChooseCaseColName (DiseaseSeverity = WhichDisease, Phase = Phase, ModelRunDisease = modelrun$disease, ModelRunHosp = modelrun$hosp)
	PatientsToSelect 	= ChooseRowIndices	(country = country, CountriesFitted = CountriesFitted, AgeGroup = AgeGroup, SeroStatus = SeroStatus, Arm = Arm)
	DummyData 			= DATA[PatientsToSelect, ]
	if (WhichDisease != "")	DummyData_nonMS = DATA_nonMildSevere[PatientsToSelect, ] #### for mild severe runs, subset the non-mildsevere data so you can determine cases in active or passive phase easily (e.g. IsCase_Severe has no info on phase in which case occured, and you need this info). 
	
	### Choose Number of participants (basically the number of patients to select, but if doing passive phase you must subtract those who were either lost to follow up or who were cases in the active phase. 
	if (Phase == "Passive")
	{
		if (WhichDisease == "") 
			NumParticipants_stratum = length(which(DummyData$IsCaseActive 		== 0 	& DummyData$LTFU == 0)) else ## in MildSevere, there is no column called IsCase_Active
			NumParticipants_stratum = length(which(DummyData_nonMS$IsCaseActive == 0 	& DummyData$LTFU == 0)) 
		
	} else NumParticipants_stratum 	= dim(DummyData)[1]
	
	PopSizes = NumParticipants_stratum
	if (PopSizes == 0) { warning("ObservedAttackRate PopSizes = 0"); 	return(NA) 		}
	
	#### Subset data
	PatientsToSelect = ChooseRowIndices(country = country, CountriesFitted = CountriesFitted, AgeGroup = AgeGroup, SeroStatus = SeroStatus, Arm = Arm)
	
	if (length(PatientsToSelect) == 0 ) stop("specified data subset of dim 0")
	DummyData = DATA[PatientsToSelect, ]
	
	if (WhichDisease == "") ## i.e. Either, from non-mild-severe runs
	{
		if (Phase == "Both"		)	Numerator = sum(DummyData$IsCase)			else 	## WholeTrial
		if (Phase == "Active"	)	Numerator = sum(DummyData$IsCaseActive)	else 	## Active Phase only
		if (Phase == "Passive"	)	Numerator = sum(DummyData$IsCasePassive)			## Passive Phase only
		
	} else if (WhichDisease == "_Mild")
	{
		if (Phase == "Both"		)	Numerator 	= sum(DummyData$IsCaseMild)			else	## WholeTrial
		if (Phase == "Active"	)	Numerator 	= sum(DummyData$IsCaseMild)			else	## Active Phase only
		if (Phase == "Passive"	)	Numerator 	= length(which(DummyData$IsCaseMild == 1 & DummyData_nonMS$IsCasePassive == 1))	## Passive Phase only	
		
	} else if (WhichDisease == "_Severe")
	{
		if (Phase == "Both"		)	Numerator 	= sum(DummyData$IsCaseSevere)			                               			 	else    ## WholeTrial         
		if (Phase == "Active"	)	Numerator 	= length(which(DummyData$IsCaseSevere == 1 & DummyData_nonMS$IsCaseActive == 1))	else    ## Active Phase only  
		if (Phase == "Passive"	)	Numerator 	= length(which(DummyData$IsCaseSevere == 1 & DummyData_nonMS$IsCasePassive == 1))	        ## Passive Phase only 
		
	} else if (WhichDisease == "_Either")
	{
		if (Phase == "Both"		)	Numerator 	= sum(DummyData$IsCase)				                                                else 	## WholeTrial               
		if (Phase == "Active"	)	Numerator 	= length(which(DummyData$IsCase == 1 & DummyData_nonMS$IsCaseActive == 1))			else    ## Active Phase only        
		if (Phase == "Passive"	)	Numerator 	= length(which(DummyData$IsCase == 1 & DummyData_nonMS$IsCasePassive == 1))		        ## Passive Phase only       
	}
	
	if (Phase == "Both"		) Denominator = sum(DummyData$EndPassive 	- DummyData$Start)		else  ## WholeTrial         
	if (Phase == "Active"	) Denominator = sum(DummyData$EndActive 	- DummyData$Start)		else  ## Active Phase only 
	if (Phase == "Passive"	) Denominator = sum(DummyData$EndPassive 	- DummyData$EndActive)        ## Passive Phase only
	
	if (Units == "years") Denominator 	= Denominator	/ 365 ### convert rate from days to years
	Obs_Attack = Numerator / Denominator
	
	#### Get confidence intervals
	BBB 				= binconf(Obs_Attack * PopSizes, PopSizes, method = "exact")
	Obs_Attack_LowerCI 	= BBB[, "Lower"] 
	Obs_Attack_UpperCI 	= BBB[, "Upper"] 
	
	if (percent)
	{
		Obs_Attack_LowerCI 	= Obs_Attack_LowerCI 	* 100
		Obs_Attack_UpperCI 	= Obs_Attack_UpperCI 	* 100
		Obs_Attack		 	= Obs_Attack 			* 100
	}
	
	MeanAndCIs 			= c(Obs_Attack, Obs_Attack_LowerCI, Obs_Attack_UpperCI)
	names(MeanAndCIs) 	= c("Mean", "LowerCI", "UpperCI")
	
	return (MeanAndCIs)
}
PredictedAttackRate	= function(country = country, AgeGroup = 0, SeroStatus = "Both_Sero", Arm = "Either", WhichDisease = "", Phase = "Both", ImSub = FALSE,
		Units = "years", percent = TRUE, IncMeanModeEtc = FALSE, WBIC = FALSE)
{
	if (!any( c("Active", "Passive", "Both") == Phase)) 					stop("PredictedAttackRate error: Phase argument not recognized.")
	if (!any( c("SeroNeg", "SeroPos", "Both_Sero") == SeroStatus)) 			stop("PredictedAttackRate error: SeroStatus argument not recognized.")
	if (!any( c("Control", "Vaccine", "Either") == Arm)) 					stop("PredictedAttackRate error: Arm argument not recognized.")
	if (!any( c("", "_Mild", "_Severe", "_Either") == WhichDisease)) 		stop("PredictedAttackRate error: WhichDisease argument not recognized.")
	
	stratum 	= ChooseStratum(AgeGroup, country, Arm, SeroStatus)
	ImSubString = GetImSubString(ImSub)
	Rates 		= Get_AR_Samples(Country = country, AgeGroup = AgeGroup, SeroStatus = SeroStatus, Arm = Arm, WhichDisease = WhichDisease, Phase = Phase, ImSub = ImSub, stratum = stratum)
	
	Attack_Mean		= rowMeans(Rates)
	Quantiles 		= quantile(Rates, c(0.025, 0.975))
	Attack_LowerCI 	= Quantiles[,1]
	Attack_UpperCI 	= Quantiles[,2]
	
	MeanAndCIs 			= c(Attack_Mean, Attack_LowerCI, Attack_UpperCI)
	names(MeanAndCIs) 	= c("Mean", "LowerCI", "UpperCI")
	
	if (Units != "years") 	MeanAndCIs 	= MeanAndCIs / 365
	### Cpp output in percentages, so must convert from percentages to probabilities. 
	if (!percent) 			MeanAndCIs 	= MeanAndCIs / 100
	
	if (IncMeanModeEtc) 
	{
		AttackRateTable_MeanMode_Char = paste0(WBIC_String, "AttackRateTable")
		if (any( c("Active", "Passive") == Phase)) AttackRateTable_MeanMode_Char = paste0(AttackRateTable_MeanMode_Char, "_", Phase, "Only")
		AttackRateTable_MeanMode_Char = paste0(AttackRateTable_MeanMode_Char, WhichDisease, "_MeanMode", ImSubString)
		SummStats			= as.numeric(get(AttackRateTable_MeanMode_Char)[stratum, ])
		names(SummStats) 	= c("Mean_Post", "Modal_Post", "Max_Like")
		### Using survival curve, we have probabilities, so must convert from probabilities to percentages. 
		if (Units != "years") 	SummStats 	= SummStats / 365
		### Cpp output in percentages, so must convert from percentages to probabilities. 
		if (!percent) 			SummStats 	= SummStats / 100
	}
	
	if (IncMeanModeEtc) VectorToReturn = c(MeanAndCIs, SummStats) else VectorToReturn = MeanAndCIs
	return (VectorToReturn)
}

PlotAttackRates 	= function (country, BaselineSeroStatus = "Both_Sero", Phase = "Both", WhichDisease = "", modelrun = ModelRun, ImSub = FALSE,
		SavePlot = TRUE, ShortenOutputString = TRUE,
		ConHazGroup_Obs = NULL, VacHazGroup_Obs = NULL, ConHazGroup_Pred = ConHazGroup_Obs, VacHazGroup_Pred = VacHazGroup_Obs, ### will usually be governed by BaselineSeroStatus argument (NULL/default), but can plot different ones if you want.  
		Units = "years", percent = TRUE, 
		TitleAddOn = "", Directory = AttackRatePlotDirectory, PNGFILENAME = NULL, PlotPredicted = TRUE, PLOTTITLE = "",
		ViolinPlot = FALSE, Inc_All_MeanModeEtc = FALSE, Inc_Mean_Post = FALSE, Inc_Modal_Post = FALSE, Inc_Max_Like = FALSE, WBIC = FALSE) 							### do this first with global variables to make it easy - tidy up later. 
{
	if (Inc_All_MeanModeEtc) {	Inc_Modal_Post = TRUE; Inc_Max_Like = TRUE 	} 
	if (any(c(Inc_Mean_Post, Inc_Modal_Post, Inc_Max_Like))) Inc_Any_MeanModeEtc = TRUE else Inc_Any_MeanModeEtc = FALSE
	
	if (!any( c("Both", "Active", "Passive") 		== Phase)) 					stop("PlotAttackRates error: Phase argument not recognized.")
	if (!any( c("Both_Sero", "SeroNeg", "SeroPos") 	== BaselineSeroStatus)) 	stop("PlotAttackRates error: BaselineSeroStatus argument not recognized")
	
	#### Build Filenames and plottitle. 
	if (Phase == "Both" 	) Addtional_TitleString = ""
	if (Phase == "Active" 	) Addtional_TitleString = "_Act"
	if (Phase == "Passive" 	) Addtional_TitleString = "_Pass"
	
	if (BaselineSeroStatus == "SeroPos"	) 	Addtional_TitleString = paste0(Addtional_TitleString, "_SPos")
	if (BaselineSeroStatus == "SeroNeg"	) 	Addtional_TitleString = paste0(Addtional_TitleString, "_SNeg")	
	if (ViolinPlot)							Addtional_TitleString = paste0(Addtional_TitleString, "_Violin") else
	if (Inc_All_MeanModeEtc)				Addtional_TitleString = paste0(Addtional_TitleString, "_ML_MP") 	
	if (Inc_All_MeanModeEtc & WBIC)			Addtional_TitleString = paste0(Addtional_TitleString, "_WBIC") 	
	
	CountryBitOfPNGTitle 	= ChooseCountryBitOfPNGTitle(country)
	WBIC_String				= Choose_WBIC_String(WBIC)
	ImSubString 			= GetImSubString(ImSub)
	if (is.null(PNGFILENAME)) ### if summarising data
		PNGFILENAME = paste0("AR", CountryBitOfPNGTitle, WhichDisease, ImSubString, Addtional_TitleString)
	PNGFILENAME 	= file.path(Directory, PNGFILENAME)
	PNGFILENAME 	= paste0(PNGFILENAME, TitleAddOn)
	PNGFILENAME 	= paste0(PNGFILENAME, ".png")
	
	#### Choose AgeGroups
	if (country <= 4 | country == 10) WhichAgeGroups = 0:3  		else	#### CYD 14 grouping 1
	if (country >= 5 | country == 11) WhichAgeGroups = c(0, 4,5)	 		#### CYD 15 grouping 1
	
	if (is.null(PLOTTITLE))
	{
		PLOTTITLE = paste0("Country ", country, " ", CountryNames[country+1], WhichDisease, ImSubString, Addtional_TitleString)
		PLOTTITLE = paste0(PLOTTITLE, TitleAddOn)
		
		if (ShortenOutputString) ModelRun_string = AbbreviateOutputString(OutputSubDir, modelrun) else ModelRun_string = OutputSubDir
		PLOTTITLE = paste0(PLOTTITLE, "\n", ModelRun_string)
	}
	
	### Observed & Predicted Attack rates (though wasteful, keep Predicted_Attack_rates even if doing ViolinPlot)
	ColNames_Predicted_Attack_rates = c("Mean", "LowerCI", "UpperCI")
	if (Inc_Any_MeanModeEtc) 	ColNames_Predicted_Attack_rates = c(ColNames_Predicted_Attack_rates, "Mean_Post", "Modal_Post", "Max_Like")
	Predicted_Attack_rates 	= matrix(nrow = length(WhichAgeGroups) * 2, ncol = length(ColNames_Predicted_Attack_rates))
	Observed_Attack_rates	= matrix(nrow = length(WhichAgeGroups) * 2, ncol = 3) 
	colnames(Predicted_Attack_rates) = ColNames_Predicted_Attack_rates
	colnames(Observed_Attack_rates) = c("Mean", "LowerCI", "UpperCI")
	
	if (ImSub & BaselineSeroStatus == "Both_Sero") SeroStatusArg = "_ImSub" else SeroStatusArg = BaselineSeroStatus
	for (age in 1:length(WhichAgeGroups))
	{
		Predicted_Attack_rates[(age-1) * 2 + 1, ] 	= PredictedAttackRate(country = country, AgeGroup = WhichAgeGroups[age], SeroStatus = BaselineSeroStatus, Arm = "Control", WhichDisease = WhichDisease, Phase = Phase, ImSub = ImSub, Units = Units, percent = percent, IncMeanModeEtc = Inc_Any_MeanModeEtc, WBIC = WBIC)
		Predicted_Attack_rates[(age-1) * 2 + 2, ] 	= PredictedAttackRate(country = country, AgeGroup = WhichAgeGroups[age], SeroStatus = BaselineSeroStatus, Arm = "Vaccine", WhichDisease = WhichDisease, Phase = Phase, ImSub = ImSub, Units = Units, percent = percent, IncMeanModeEtc = Inc_Any_MeanModeEtc, WBIC = WBIC)
		
		Observed_Attack_rates[(age-1) * 2 + 1, ] 	= ObservedAttackRate(country = country, AgeGroup = WhichAgeGroups[age], SeroStatus = SeroStatusArg, Arm = "Control", WhichDisease = WhichDisease, Phase = Phase, Units = "years", percent = percent, modelrun = modelrun)
		Observed_Attack_rates[(age-1) * 2 + 2, ] 	= ObservedAttackRate(country = country, AgeGroup = WhichAgeGroups[age], SeroStatus = SeroStatusArg, Arm = "Vaccine", WhichDisease = WhichDisease, Phase = Phase, Units = "years", percent = percent, modelrun = modelrun)
	}
	
	Obs_Attack			= Observed_Attack_rates	[,"Mean"]
	Obs_Attack_LowerCI 	= Observed_Attack_rates	[,"LowerCI"]
	Obs_Attack_UpperCI 	= Observed_Attack_rates	[,"UpperCI"]
	
	Attack_Mean			= Predicted_Attack_rates[,"Mean"]   
	Attack_LowerCI 		= Predicted_Attack_rates[,"LowerCI"]
	Attack_UpperCI 		= Predicted_Attack_rates[,"UpperCI"]
	if (Inc_Any_MeanModeEtc) 
	{
		Attack_Mean_Post	= Predicted_Attack_rates[,"Mean_Post"]   
		Attack_Modal_Post 	= Predicted_Attack_rates[,"Modal_Post"]
		Attack_Max_Like 	= Predicted_Attack_rates[,"Max_Like"]
	}
	
	LWD 	= 4	 	
	CEXAXIS = 1
	CEX.LAB = 1.4
	YLIM 	= c(0, min(c(10, max(c(Attack_Mean, Attack_LowerCI, Attack_UpperCI, Obs_Attack, Obs_Attack_UpperCI, Obs_Attack_LowerCI))) * 1.05) ) ### i.e. no higher than 10, but otherwise big enough to show all the upper and lower and mean attack rates of either data or model
	COLS 	= rep(c("cadetblue1", "pink"), length(Attack_Mean))
	LABELS 	= paste(rep(AgeGroupNames[WhichAgeGroups + 1], each = 2), "\n",  c("Control", "Vaccine"), sep = "")
	
	if (PlotPredicted) Offset = 0.15 else Offset = 0  ### if summarising data don't plot predicted so observed can be exactly centred. 
	if (SavePlot)	png(file = PNGFILENAME  , res = PNG_res, units = "in", width = 7, height= 5)
	plot(	NA, 
			ylab = "Attack Rate (%)", cex.lab = CEX.LAB, cex.main = 1.4, xaxt = "n", xlab = "", ylim = YLIM, 
			col = COLS, lwd = LWD, cex = 2, xlim = c(0.5, length(Attack_Mean) + 0.5),	
			main = PLOTTITLE )
	
	axis(side = 1, labels= rep("", length(Attack_Mean)), at = 1:length(Attack_Mean), cex.axis = CEXAXIS)
	mtext(side = 1, text= LABELS, at = 1:length(Attack_Mean), cex.axis = CEXAXIS, line = 2)
	mtext(side = 1, text= "Age Group", at = mean(1:length(Attack_Mean)), cex =  1.9, line = 4)
	points(1:length(Obs_Attack) - Offset, Obs_Attack, col = c("blue", "red"), pch = 4, cex = 2, lwd = LWD)
	Arrows(1:length(Obs_Attack) - Offset, Obs_Attack_LowerCI 	, 	1:length(Attack_Mean) - Offset,    Obs_Attack_UpperCI, code = 3 , arr.type = "T" , col = c("blue", "red"), lwd = LWD)
	
	if (PlotPredicted)
	{
		if (!ViolinPlot)
		{
			Arrows( 1:length(Attack_Mean) + Offset,  Attack_LowerCI 	, 	1:length(Attack_Mean) + Offset,    Attack_UpperCI, code = 3 , arr.type = "T" , col = COLS, lwd = LWD)
			points(1:length(Attack_Mean) + Offset, Attack_Mean, col = COLS, pch = 1, cex = 2, lwd = LWD)
			if (Inc_Any_MeanModeEtc) 	{	PCH_MeanPost = 7; PCH_Modal = 6; PCH_MaxLike = 5;	}
			if (Inc_Mean_Post)	points(1:length(Attack_Mean_Post) 	+ (2*Offset), Attack_Mean_Post	, col = c("darkseagreen1", "darksalmon"), pch = PCH_MeanPost	, cex = 2, lwd = LWD)
			if (Inc_Modal_Post)	points(1:length(Attack_Modal_Post) 	+ (2*Offset), Attack_Modal_Post	, col = c("green", "orange")			, pch = PCH_Modal		, cex = 2, lwd = LWD)
			if (Inc_Max_Like)	points(1:length(Attack_Max_Like) 	+ (2*Offset), Attack_Max_Like	, col = c("darkgreen", "darkorange4")	, pch = PCH_MaxLike		, cex = 2, lwd = LWD)
			
		}	else	
		{
			counter = 1
			for (AgeGroup in WhichAgeGroups)
				for (TrialArm in c("Control", "Vaccine"))
				{
					vec = as.numeric(Get_AR_Samples(Country = country, AgeGroup = AgeGroup, SeroStatus = BaselineSeroStatus, Arm = TrialArm, WhichDisease = WhichDisease, Phase = Phase, ImSub = ImSub))
					vioplot(vec, at = counter + Offset, wex = 0.5, col=COLS[counter], add = TRUE)
					counter = counter + 1
				}
		}
		legend("topleft", 
				c(	 "Trial", "Model"), pt.cex = 1.8 ,
				pch = c(4, 1) ,
				lty = c(1, 1), lwd = rep(LWD, 2) , cex = 1.1	)
		
		if (Inc_Any_MeanModeEtc & !ViolinPlot) 
			legend("topright", 
					c(	 "ModalPosterior", "MaxLike"), pt.cex = 1.8 ,
					pch = c(PCH_Modal, PCH_MaxLike) ,
					lty = c(NA, NA), lwd = rep(LWD, 2) , cex = 1.1	)
		
	} else	legend("topleft", c( "Control", "Vaccine"), pt.cex = 1.8 , col = c("blue", "red"), pch = c(4,4), lty = c(1, 1), lwd = rep(LWD, 2) , cex = 1.1	)
	
	if (SavePlot)	dev.off()
}

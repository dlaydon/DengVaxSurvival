

Convert_DiseaseString 		= function(WhichDiseasePlot, ModelRunDisease, ModelRunHosp) 
{
	DiseaseString = WhichDiseasePlot
	if (ModelRunDisease == "_MILDSEVERE" & ModelRunHosp == "_hosp")
	{
		if (WhichDiseasePlot == "_Mild"		) DiseaseString = "_NonHosp"
		if (WhichDiseasePlot == "_Severe"	) DiseaseString = "_Hosp"
	}
	return(DiseaseString)
}
ChooseCaseColName 			= function(DiseaseSeverity, Phase = "Both", ModelRunDisease, ModelRunHosp)
{
	if (!any( c("", "_Mild", "_Severe", "_Either") == DiseaseSeverity)) 	stop("ChooseCaseColName error: DiseaseSeverity argument not recognized.")
	if (!any( c("Active", "Passive", "Both") == Phase)) 					stop("ChooseCaseColName error: Phase argument not recognized.")
	
	CaseColName = "IsCase"
	
	if (ModelRunDisease == "")
	{
		if (Phase == "Active"	) CaseColName = paste0(CaseColName, "Active"	)
		if (Phase == "Passive"	) CaseColName = paste0(CaseColName, "Passive"	)
		
	} else if (ModelRunDisease == "_MILDSEVERE")
	{
		if (DiseaseSeverity == "_Mild"	) CaseColName = paste0(CaseColName	, "Mild")
		if (DiseaseSeverity == "_Severe") CaseColName = paste0(CaseColName	, "Severe"		)
	}
	return(CaseColName)
}
ChooseRowIndices_Country 	= function(country, CountriesFitted)
{
	if (country <  10) RowIndices 	= which(DATA$Country == country) 								else
	if (country == 10) RowIndices 	= which(DATA$Country <= 4 & DATA$Country %in% CountriesFitted)	else
	if (country == 11) RowIndices 	= which(DATA$Country >= 5 & DATA$Country %in% CountriesFitted)	else 
	if (country == 12) RowIndices 	= which(DATA$Country %in% CountriesFitted)						else
		cat("ChooseRowIndices_Country error: country argument not chosen correctly")
	return(RowIndices)
}
ChooseRowIndices_AgeGroup 	= function(AgeGroup = 0)
{
	if (AgeGroup == 0 || is.null(AgeGroup) )	
		RowIndices = 1:dim(DATA)[1]	else	 
		RowIndices = which(	DATA$AgeGroup1 == AgeGroup | DATA$AgeGroup2 == AgeGroup)
	
	return(RowIndices)
}
ChooseRowIndices_AgeRange 	= function(AgeMin = NULL, AgeMax = NULL)
{
	if ( is.null(AgeMin) & is.null(AgeMax)) 	RowIndices = 1:dim(DATA)[1] 									else
	if (!is.null(AgeMin)& !is.null(AgeMax))		RowIndices = which (DATA$Age >= AgeMin & DATA$ais <= AgeMax) 	else
	if (!is.null(AgeMin)) 						RowIndices = which (DATA$Age >= AgeMin) 						else
	if (!is.null(AgeMax)) 						RowIndices = which (DATA$Age <= AgeMax)
	
	return(RowIndices)
}
ChooseRowIndices_SeroStatus	= function(SeroStatus = "Both_Sero")
{
	if (SeroStatus == "Both_Sero"	) RowIndices = 1:dim(DATA)[1] 						else 	
	if (SeroStatus == "SeroNeg"		) RowIndices = which(DATA$SeroStatus == 0)			else	
	if (SeroStatus == "SeroPos"		) RowIndices = which(DATA$SeroStatus == 1) 			else
	if (SeroStatus == "_ImSub"		) RowIndices = which(DATA$SeroStatus != MDValue) 	else 	### in case you want vaccine and control arms of only the immune subset, without reference to their actual serostatus. 
		cat("ChooseRowIndices_SeroStatus error: SeroStatus argument not chosen correctly")
	
	return(RowIndices)
}
ChooseRowIndices_Arm		= function(Arm = "Either")
{
	if (Arm == "Either"	) RowIndices = 1:dim(DATA)[1] 		else 	
	if (Arm == "Control") RowIndices = which(DATA$Arm == 0) else 			
	if (Arm == "Vaccine") RowIndices = which(DATA$Arm == 1) else 
		cat("ChooseRowIndices_Arm error: Arm argument not chosen correctly")
	
	return(RowIndices)
}
ChooseRowIndices 			= function(country, CountriesFitted, AgeGroup = 0, SeroStatus = "Both_Sero", Arm = "Either", AgeMin = NULL, AgeMax = NULL)
{
	if (AgeGroup != 0 & (!is.null(AgeMin) || !is.null(AgeMax))) warning("ChooseRowIndices: Age range selected for AgeGroup != 0")
	
	PatientsToSelect 		= 1:dim(DATA)[1]
	RowIndicesCountry 		= ChooseRowIndices_Country		(country, CountriesFitted)
	if (is.null(AgeMin) && is.null(AgeMax)) 
		RowIndicesAgeGroup 		= ChooseRowIndices_AgeGroup	(AgeGroup) else RowIndicesAgeGroup = 1:dim(DATA)[1]
	RowIndicesSeroStatus 	= ChooseRowIndices_SeroStatus	(SeroStatus)
	RowIndicesArm 			= ChooseRowIndices_Arm			(Arm)
	RowIndices_AgeRange		= ChooseRowIndices_AgeRange		(AgeMin, AgeMax)
	
	PatientsToSelect 		= intersect (PatientsToSelect, RowIndicesCountry 		)	
	PatientsToSelect 		= intersect (PatientsToSelect, RowIndicesAgeGroup 		)	
	PatientsToSelect 		= intersect (PatientsToSelect, RowIndicesSeroStatus 	)	
	PatientsToSelect 		= intersect (PatientsToSelect, RowIndicesArm 			)	
	PatientsToSelect 		= intersect (PatientsToSelect, RowIndices_AgeRange 		)	
	length(PatientsToSelect)
	return(PatientsToSelect)
}
ChooseHazardNumber 			= function(Arm, SeroStatus)
{
	if (!any( c("SeroNeg", "SeroPos", "Both_Sero") == SeroStatus)) 			stop("ChooseHazardNumber error: SeroStatus argument not recognized.")
	if (!any( c("Control", "Vaccine", "Either") == Arm)) 					stop("ChooseHazardNumber error: Arm argument not recognized.")
	
	if (Arm == "Control")
	{
		if (SeroStatus == "SeroPos") HazardNumber = 1 else if (SeroStatus == "SeroNeg") HazardNumber = 2 else HazardNumber = 5
		
	} else if (Arm == "Vaccine")
	{
		if (SeroStatus == "SeroPos") HazardNumber = 3 else if (SeroStatus == "SeroNeg") HazardNumber = 4 else HazardNumber = 6
		
	} else if (Arm == "Either") 
	{
		if (SeroStatus == "SeroPos") HazardNumber = 7 else if (SeroStatus == "SeroNeg") HazardNumber = 8 else HazardNumber = 9
	} 
	return(HazardNumber)
}
ChooseStratum 				= function(AgeGroup, Country, Arm, SeroStatus)
{
	if (!any( c("SeroNeg", "SeroPos", "Both_Sero") == SeroStatus)) 			stop("ChooseStratum error: SeroStatus argument not recognized.")
	if (!any( c("Control", "Vaccine", "Either") == Arm)) 					stop("ChooseStratum error: Arm argument not recognized.")
	
	CountryAgeGroupIndex = (AgeGroup * 13 * 9) + (Country * 9) #  * 13 countries * 9 hazard groups
	
	### loads of possibilities for arm and serostatus - given haphazard way hazard groups are numbered, best to do this with an exhaustive list. 
	## In R indices, hazard groups are: 1 = SeroPos Control, 2 = SeroNeg Control, 3 = SeroPos Vaccine, 4 = SeroNeg Vaccine, 5 = Control,  6 = Vaccine,  7 = SeroPositive, 8 = SeroNegative, 9 = everyone in that country. 
	
	HazardNumber 		= ChooseHazardNumber(Arm, SeroStatus)
	
	HazardGroupIndex 	= CountryAgeGroupIndex + HazardNumber
	
	return(HazardGroupIndex)
}
ChooseCatAgeGroup			= function(Age)
{
	#### this function used to get age group (CYD-14) from Age. Saves writing the code below all the time. Works for CYD-15, 
	if (Age < 0) stop("ChooseCatAgeGroup error: Age < 0") else 
	if (Age < 6) 	AG = 1 else
	if (Age < 12) 	AG = 2 else 	
		AG = 3
	return(AG)
}


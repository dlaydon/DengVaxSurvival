Quadratic					= function(x, C_0, C_1, C_2)
{
	return(C_0 + C_1*x + C_2*x*x)
}
Cubic						= function(x, C_0, C_1, C_2, C_3)
{
	return(C_0 + C_1*x + C_2*x*x + C_3*x*x*x)
}
FindCoefficients_Cubic			= function(x0, x1, x2, x3, y0, y1, y2, y3)
{
	#### Definitions: Makes rest of code MUCH easier. 
	L_20 = (x1 - x0) * (y2				- y0)				- (x2 - x0) * (y1			- y0);				## was labelled K20 in notes
	L_30 = (x1 - x0) * (y3				- y0)				- (x3 - x0) * (y1			- y0);				## was labelled K30 in notes
	
	Q_20 = (x1 - x0) * (x2 * x2			- x0 * x0)			- (x2 - x0) * (x1 * x1		- x0 * x0);			## was labelled D20 in notes
	Q_30 = (x1 - x0) * (x3 * x3			- x0 * x0)			- (x3 - x0) * (x1 * x1		- x0 * x0);			## was labelled D30 in notes
	
	C_20 = (x1 - x0) * (x2 * x2 * x2	- x0 * x0* x0)		- (x2 - x0) * (x1 * x1* x1	- x0 * x0 * x0);	## was labelled Q20 in notes
	C_30 = (x1 - x0) * (x3 * x3 * x3	- x0 * x0* x0)		- (x3 - x0) * (x1 * x1* x1	- x0 * x0 * x0);	## was labelled Q30 in notes
	
	a3 = (L_30*Q_20 - L_20*Q_30) / (C_30*Q_20 - C_20*Q_30); 
	a2 = (L_20 - a3 * C_20) / Q_20;	## these should be equivalent.
	#a2 = (L_30 - a3 * C_30) / Q_30;	## these should be equivalent. 
	a1 = (y1-y0 - a2*(x1*x1 - x0*x0) -a3*(x1*x1*x1 - x0*x0*x0)) / (x1-x0);
	a0 = y0 - a1*x0 - a2*x0*x0 - a3* x0*x0*x0; 
	
	return (c(a0, a1, a2, a3))
}
FindCoefficients			= function(x1, x2, x3, y1, y2, y3)
{
	a2 = (y3-y2 - ((x3-x2) * (y1-y2)/(x1-x2) )  ) / (x3*x3 - x2*x2 - (x1+x2)*(x3-x2)) 
	a1 = ((y1-y2)/(x1-x2)) - (x1+x2)* a2
	a0 = y1-(a1*x1)-(a2*x1*x1)
	
	return (c(a0, a1, a2))
}
Choose_xKnots 				= function(Type, modelrun)
{
	if (Type == "Hazards") 												xKnots = (1/3 * 1:KnotsPerCountry)				else 
	if (Type == "Age_VacDuration" && modelrun$AS_Waning_KnotSet == 0) 	xKnots = c(2, 6, 12, 16) 						else
	if (Type == "Age_VacDuration" && modelrun$AS_Waning_KnotSet == 1) 	xKnots = c(2.0, 5.9, 6.1, 11.9, 12.1, 16.0) 	else
	if (Type == "Age_VacDuration" && modelrun$AS_Waning_KnotSet == 2) 	xKnots = c(2.0, 4.8, 7.6, 10.4, 13.2, 16.0) 	else
	if (Type == "Age_VacDuration" && modelrun$AS_Waning_KnotSet == 3) 	xKnots = c(2.0, 5.5, 9.0, 12.5, 16.0) 			else
	if (any(c("Age_Efficacy", "Age_Haz_Mult") == Type)) 				xKnots = c(2, 6, 12, 16) 					else
		stop("Choose_xKnots error: Type argument not recognized - xKnots not chosen")
}
Choose_BSS_Char_Splines 	= function(Type, modelrun = ModelRun, BaselineSeroStatus = NULL)
{
	if (Type == "Age_Efficacy" && BaselineSeroStatus == "SeroNeg") 			BSS_Char = "SNeg" else ## note the double &, works as a circuit breaker if BaselineSeroStatus = NULL
	if (Type == "Age_Efficacy" && BaselineSeroStatus == "SeroPos") 			BSS_Char = "SPos" else ## note the double &, works as a circuit breaker if BaselineSeroStatus = NULL
	if (Type == "Age_VacDuration" & modelrun$AS_Waning_Homogeneous == 0) 	BSS_Char = paste(BaselineSeroStatus) else	BSS_Char = ""
	
	return(BSS_Char)
}
ChooseKnotIndices 			= function(Type, modelrun = ModelRun, BaselineSeroStatus = NULL, country = NULL)
{
	BSS_Char = Choose_BSS_Char_Splines(Type, modelrun, BaselineSeroStatus)
	
	## Choose Knot indices from CHAINS
	if (Type == "Hazards"			) KnotIndices 	= grep(paste0("Knot_", country, "_")				, colnames(CHAINS)) else
	if (Type == "Age_Efficacy"		) 
	{
		if (modelrun$ASVE == "CATEGORICAL")
			KnotIndices 	= which(colnames(CHAINS) %in% paste0(BSS_Char, "Eff_AgeGroup_"	, 1:3)) else
			KnotIndices 	= which(colnames(CHAINS) %in% paste0(BSS_Char, "Eff_AgeKnot_"	, 0:3))
		
	}	else if (Type == "Age_Haz_Mult") 
	{
		if (modelrun$AS_Haz == "CATEGORICAL") 
			KnotIndices 	= which(colnames(CHAINS) %in% paste0("AS_Haz_AgeGroupMult_", 1:3)) else	### weirdly, for categorical variables, you have two age group parameters (for age group 0 and 1) fixed to 1. You'll want 1(baseline), 2 and 3 for plots, not 0 (was because you were originally going to refer to participants age group in Cpp likelihood, but that didn't pan out - you're now stuck with this. 
			KnotIndices 	= which(colnames(CHAINS) %in% paste0("AS_Haz_Knot_", 0:3)) 				
		
	}	else if (Type == "Age_VacDuration") 
	{
		if (modelrun$AS_Waning == "CATEGORICAL") ParamNameType = "_AG_" else ParamNameType = "_AgeKnot_" 
		KnotIndices 	= grep(paste0(BSS_Char, "Waning", ParamNameType), colnames(CHAINS))
		
	}	else stop("ChooseKnotIndices error: KnotIndices if else tree wrong. ")
	
	return(KnotIndices)
}
FindWhichPolynomial			= function(xValue, xKnots)
{
	## not the same as finding previous knot - wouldn't work for final polynomial as this connects 3 knots, not 2. 
	WhichPolynomial = length(xKnots) - 2;  
	while (WhichPolynomial > 0 && xKnots[WhichPolynomial] > xValue) WhichPolynomial = WhichPolynomial - 1
	return (WhichPolynomial)
}
Quadratic_wrapper			= function(xValue, xKnots, yKnots, Coeffs)
{
	## calculate value of spline, given coefficients, x and y knots, and knots per spline. Essentially wrapper of Quadratic that extends past the knots at the extremes. 
	#### If xValue is less than / greater than extremes, spline takes constant value of the corresponding yKnot. 
	if (xValue <= xKnots[1]			)	yValue = yKnots[1]				else 
	if (xValue >= tail(xKnots,1)	)	yValue = tail(yKnots,1)			else 	
		yValue = Quadratic (xValue, Coeffs[1], Coeffs[2], Coeffs[3])
	return (yValue)
}
Choose_yKnots				= function(Type, ParamSetNumber, BaselineSeroStatus = NULL, modelrun = ModelRun, country = NULL)
{
	KnotIndices = ChooseKnotIndices(Type, modelrun, BaselineSeroStatus, country)
	yKnots 		= CHAINS[ParamSetNumber, KnotIndices]
	return(as.numeric(yKnots))
}
Calculate_Mean_yKnots		= function(Type, BaselineSeroStatus = NULL, modelrun = ModelRun, country = NULL)
{
	KnotIndices 		= ChooseKnotIndices(Type, modelrun, BaselineSeroStatus, country)
	Mean_yKnots 		= colMeans(CHAINS[, KnotIndices])
	return(as.numeric(Mean_yKnots))
}
Spline_x_sample 			= function(xValue, Type, modelrun = ModelRun, ParamSetNumber = 1, BaselineSeroStatus = NULL, xKnots = NULL, yKnots = NULL, yMin = NULL, yMax = NULL, Multiplier = 1, country = NULL)
{
	# Multiplier is used when doing Age_Efficacy splines with SS_VEs, the "serotype efficacy" terms are multipliers of the age efficacy profile. 
	#### For a given plot type (e.g. baseline hazard, seropositive cross-immunity / duration, splines), and given posterior sample, calculate value of spline. 
	# # # # # # # # # # 		Error traps
	# # # # # Wrong arguments
	if (!any( c("Hazards", "Age_Efficacy", "Age_Haz_Mult", "Age_VacDuration") == Type)) 		stop("Spline_x_sample error: Type argument not recognized.")
	if (!any(c("SeroNeg", "SeroPos") == BaselineSeroStatus) & !is.null(BaselineSeroStatus)) 	stop("Spline_x_sample error: BaselineSeroStatus argument not recognized")
	
	if (Type == "Age_Efficacy" 										& is.null(yMax)) 		yMax = 1 								
	if (Type %in% c("Age_VacDuration", "Hazards", "Age_Haz_Mult") 	& is.null(yMin)) 		yMin = 0							
	
	## Find knots (if not provided)
	if (is.null(xKnots))	xKnots	= Choose_xKnots(Type, modelrun)
	if (is.null(yKnots))	yKnots	= Choose_yKnots(Type, ParamSetNumber, BaselineSeroStatus, modelrun, country)
	
	#if (length(xKnots) != length(yKnots)) stop("Spline error: length(xKnots) != length(yKnots)") ## comment out as if doing categoricals you ignore final xKnot.
	KnotsPerSpline 	= length(xKnots)
	
	## Find WhichPolynomial
	WhichPolynomial = FindWhichPolynomial(xValue, xKnots)
	
	## Calculate relevant coefficients
	Coeffs = FindCoefficients(	xKnots[WhichPolynomial], 	xKnots[WhichPolynomial + 1], 	xKnots[WhichPolynomial + 2], ## x knots
								yKnots[WhichPolynomial], 	yKnots[WhichPolynomial + 1], 	yKnots[WhichPolynomial + 2]) ## y knots
	
	## Calculate yValue
	yValue = Quadratic_wrapper(xValue, xKnots, yKnots, Coeffs) * Multiplier
	
	if (!is.null(yMax)) yValue = min(c(yValue, yMax))
	if (!is.null(yMin)) yValue = max(c(yValue, yMin))
	
	return(yValue)
}
Spline_XVec_sample 			= function(xValues, Type, modelrun = ModelRun, ParamSetNumber = 1, BaselineSeroStatus = NULL, xKnots = NULL, yKnots = NULL, yMin = NULL, yMax = NULL, country = NULL)
{
	#### Wrapper of Spline_x_sample that allows you to plot for multiple ages / dates / x values at once, without having to re-write mapply. 
	Vec = mapply (	Spline_x_sample, xValue = xValues, MoreArgs = list(Type = Type, BaselineSeroStatus = BaselineSeroStatus, modelrun = ModelRun, ParamSetNumber = ParamSetNumber))	
	return(Vec)
}
Spline_XVec_Mean 			= function(xValues, Type, modelrun = ModelRun, BaselineSeroStatus = NULL, xKnots = NULL, Mean_yKnots = NULL, yMin = NULL, yMax = NULL)
{
	#### Wrapper of Spline_x_sample that allows you to plot for multiple ages / dates / x values at once, without having to re-write mapply. Also uses the mean knots
	if (is.null(Mean_yKnots)) Mean_yKnots = Calculate_Mean_yKnots(Type = Type, BaselineSeroStatus = BaselineSeroStatus, modelrun = modelrun)
	
	Vec = mapply (	Spline_x_sample, xValue = xValues, MoreArgs = list(Type = Type, BaselineSeroStatus = BaselineSeroStatus, modelrun = ModelRun, yKnots = Mean_yKnots))	
	return(Vec)
}
ChooseDegreeFrom_AgeOption 	= function(AgeOptionString)
{
	if (AgeOptionString == "CUBIC"		) 										Degree = 3 else 
	if (AgeOptionString == "SPLINE"		) 										Degree = 2 else 
	if (AgeOptionString == "SPLINE_LINE"	) 									Degree = 1 else 
	if (AgeOptionString == "CATEGORICAL" || AgeOptionString == "SPLINE_STEP") 	Degree = 0 else 
		stop("ChooseDegreeFrom_AgeOption error, AgeOptionString value not recognized")
	
	return (Degree)
}

PlotSplines = function(country = NULL, BaselineSeroStatus = NULL, modelrun = ModelRun, COLOUR = Alpha_ParamGuessCols, NumParamSetsToPlot = HowManyParamSetsToPlot, OutputSubDirInPlotTite = TRUE,
		SavePlot = TRUE, ShortenOutputString = TRUE, LOWERYLIM = NULL, UPPERYLIM = NULL, Directory = HazardPlotDirectory, DoubleSpline = FALSE, NoPointsPerKnotGap = 5,
		TitleAddOn = "", PlotType = "Hazards", NumDaysAroundKnots = 30, Degree = NULL, serotype = 1, PLOTTITLE = NULL, IncLegend = TRUE, LegPosChar = NULL)
{
	### Degree refers to degree of polynomial. Default is quadratic spline, but in ASWaneStep, Degree = 0, and in ASWaneLine, Degree = 1. 
	
	# # # # # # # # # # 		Error traps
	# # # # # Wrong arguments
	
	if (!any( c("Hazards", "Age_Efficacy", "Age_Haz_Mult", "Age_VacDuration") == PlotType)) 	stop("PlotSplines error: PlotType argument not recognized.")
	if (!any(c("SeroNeg", "SeroPos") == BaselineSeroStatus) & !is.null(BaselineSeroStatus)) 	stop("PlotSplines error: BaselineSeroStatus argument not recognized")
	if (any(c("red", "blue") == COLOUR)) 	stop("PlotSplines error: COLOUR of posterior samples can't be red or blue")
	
	if (is.null(Degree))
	{
		if (PlotType == "Hazards"			)	Degree = 2													else
		if (PlotType == "Age_Efficacy"		)	Degree = ChooseDegreeFrom_AgeOption(modelrun$ASVE)			else
		if (PlotType == "Age_Haz_Mult"		)	Degree = ChooseDegreeFrom_AgeOption(modelrun$AS_Haz)		else
		if (PlotType == "Age_VacDuration" 	)	Degree = ChooseDegreeFrom_AgeOption(modelrun$AS_Waning)    
	}
	if (Degree > 3) stop("PlotSplines error: Degree > 3")
	if (Degree == 0 && PlotType == "Age_Efficacy"		&& (modelrun$ASVE 		!= "CATEGORICAL" && modelrun$ASVE 		!= "SPLINE_STEP")) stop("ASVE 		Splines: Step function for non-step modelrun")
	if (Degree == 0 && PlotType == "Age_Haz_Mult"		&& (modelrun$AS_Haz 	!= "CATEGORICAL" && modelrun$AS_Haz 	!= "SPLINE_STEP")) stop("AS_Haz 	Splines: Step function for non-step modelrun")
	if (Degree == 0 && PlotType == "Age_VacDuration" 	&& (modelrun$AS_Waning 	!= "CATEGORICAL" && modelrun$AS_Waning 	!= "SPLINE_STEP")) stop("AS_Waning 	Splines: Step function for non-step modelrun")
	if (Degree == 1 && PlotType == "Age_Efficacy"		&& (modelrun$ASVE 		!= "SPLINE_LINE")) stop("ASVE 		Splines: Line function for non-step modelrun")
	if (Degree == 1 && PlotType == "Age_Haz_Mult"		&& (modelrun$AS_Haz 	!= "SPLINE_LINE")) stop("AS_Haz 	Splines: Line function for non-step modelrun")
	if (Degree == 1 && PlotType == "Age_VacDuration" 	&& (modelrun$AS_Waning 	!= "SPLINE_LINE")) stop("AS_Waning 	Splines: Line function for non-step modelrun")
	
	if (NumParamSetsToPlot != HowManyParamSetsToPlot ) warning("PlotSplines: NumParamSetsToPlot != HowManyParamSetsToPlot")
	
	# # # # # Wrong combinations of arguments
	if (any(c("Hazards", "Age_Haz_Mult") == PlotType) & !is.null(BaselineSeroStatus)) 				stop("PlotSplines error: BaselineSeroStatus != NULL for hazards or hazard multiplier spline")
	if (any(c("Age_Efficacy", "Age_VacDuration") == PlotType) & is.null(BaselineSeroStatus)) 		stop("PlotSplines error: BaselineSeroStatus == NULL for Efficacies or duration spline")
	if ("Hazards" == PlotType & is.null(country)) 													stop("PlotSplines error: country == NULL for hazards or hazard multiplier spline")
	if (any(c("Age_Efficacy", "Age_Haz_Mult", "Age_VacDuration") == PlotType) & !is.null(country)) 	stop("PlotSplines error: country != NULL for Efficacies, duration, or hazard multiplier spline")
	
	## Choose Filename
	BSS_Char = Choose_BSS_Char_Splines(PlotType, modelrun, BaselineSeroStatus)
	if (PlotType == "Age_Efficacy" & modelrun$SS_VEs == 1) SerotypeString = paste0("_sero", serotype) else SerotypeString = ""
	
	## Create file name	(use filname and date modified to determine whether it is worthwhile making new plot)
	if (PlotType == "Hazards"			) FILENAME = file.path(Directory, paste0("PredictedHazard_country_", country, 		".png")) else 
	if (PlotType == "Age_Efficacy"		) FILENAME = file.path(Directory, paste0("AgeEffSpline_", BSS_Char, SerotypeString,	".png")) else 
	if (PlotType == "Age_Haz_Mult"		) FILENAME = file.path(Directory, paste0("HazMultSpline", 							".png")) else 
	if (PlotType == "Age_VacDuration"	) FILENAME = file.path(Directory, paste0("DurationSpline_", BSS_Char, 	 			".png")) 
	
	if (file.exists(FILENAME) && (file.mtime(FILENAME) < file.mtime(ChainFileName))) 
		print(paste0("Chain new (mod ", file.mtime(ChainFileName), ") - Overwriting plot (mod ", file.mtime(FILENAME), ")" ))
	
	if (SavePlot) png(file = FILENAME, res = PNG_res, units = "in", width = 6, height= 6)
	
	## Choose Knot indices from CHAINS
	KnotIndices 	= ChooseKnotIndices(PlotType, modelrun, BaselineSeroStatus)
	## Get mean value of yKnots (parameters)
	MeanKnots 		= colMeans(CHAINS[, KnotIndices])
	if (PlotType == "Age_Efficacy" & modelrun$SS_VEs == 1) 
	{
		MeanMultOrAdd 	= colMeans(CHAINS[paste0(BSS_Char, "Eff_", serotype)])
		if (modelrun$SSASVE_Additive) MeanKnots		= MeanKnots + MeanMultOrAdd else MeanKnots		= MeanKnots * MeanMultOrAdd
	}
	
	## Define xKnots
	## easy for non-hazards, but because of dates being converted to weird integers, you plot y-axis values of splines calculated from regular xKnots and yknots, but x-axis values are dates). 
	## Further, to save time you put them in a list so you don't calculated them for every parameter set. 
	xKnots	= Choose_xKnots(PlotType, modelrun)
	COLS 	= rep(c("red"), 2*(length(xKnots)))
	
	#### remove final knot if doing categorical (or step?). Don't do this. Still need final knot so that you can seq to it for plotting. 
	if (PlotType == "Hazards") 
	{
		## create xKnots with units in calendar days. 
		Day_0CYD14 			= "2011-06-03"
		Day_0 				= as.Date(Day_0CYD14)
		xKnotsDummy 		= xKnots * 365
		xKnotsDummy 		= as.Date(Day_0 + xKnotsDummy)
		xKnotsDates 		= as.Date(xKnotsDummy, origin = Day_0)
	}
	KnotSpline_Xaxis_List = list()
	if (PlotType == "Hazards")	KnotSpline_Xaxis_AsDates_List 	= list()
	
	if (Degree == 3) LastKnot = 2 else if (Degree == 2) LastKnot = length(xKnots) - 1 else LastKnot = length(xKnots) #### For Degree == 3, you are assuming that there are only 4 knots. Not true for other Degree values. 
	
	for (WhichKnot in 2:LastKnot)
	{
		x1 = xKnots[WhichKnot-1]
		x2 = xKnots[WhichKnot  ]
		
		if (Degree == 3)
		{
			x3 = xKnots[WhichKnot+1]
			x4 = xKnots[WhichKnot+2]
			KnotSpline_Xaxis_List[[WhichKnot]]  = seq(x1, x4, length.out = 	3 *	NoPointsPerKnotGap)
			
		} else if (Degree == 2)
		{
			x3 = xKnots[WhichKnot+1]
			
			if (DoubleSpline)						KnotSpline_Xaxis_List[[WhichKnot]] 	= seq(x1, x3, length.out = 		NoPointsPerKnotGap) else
			if (WhichKnot == length(xKnots) - 1)	KnotSpline_Xaxis_List[[WhichKnot]] 	= seq(x1, x3, length.out = 2 * 	NoPointsPerKnotGap) else
				KnotSpline_Xaxis_List[[WhichKnot]]  = seq(x1, x2, length.out = 		NoPointsPerKnotGap) 
			
		} else 										KnotSpline_Xaxis_List[[WhichKnot]]  = seq(x1, x2, length.out = 		NoPointsPerKnotGap) 
		
		if (PlotType == "Hazards")				KnotSpline_Xaxis_AsDates_List[[WhichKnot]] = KnotSpline_Xaxis_List[[WhichKnot]] * 365 + as.numeric(Day_0)
	}
	
	## Get x and y axis limits (SeroNegative efficacies should be a) allowed to go negative and b) but not too negative so that plot is unreadable, say -100% 
	if (is.null(LOWERYLIM))
	{
		if (PlotType == "Age_Efficacy" && BaselineSeroStatus == "SeroNeg") 	 
			LOWERYLIM = max(c(-1, min(CHAINS[WhichParamsToPlot, KnotIndices]))) else LOWERYLIM = 0 ### need the double &&, as a circuit breaker in case is.null(BaselineSeroStatus)
	}
	if (is.null(UPPERYLIM))		UPPERYLIM = max(CHAINS[WhichParamsToPlot, KnotIndices])	else 	UPPERYLIM = min(c(UPPERYLIM, max(CHAINS[WhichParamsToPlot, KnotIndices])))
	
	YLIM = c(LOWERYLIM,	UPPERYLIM)
	if (is.null(LOWERYLIM) & is.null(UPPERYLIM))		if (PlotType == "Age_Efficacy") YLIM = c(-1, 1) #### redefine for efficacies. 
	
	if (PlotType == "Hazards")	XLIM = c(   min(as.numeric(xKnotsDates)) - NumDaysAroundKnots   , max(as.numeric(xKnotsDates)) + NumDaysAroundKnots ) else XLIM = c(2, 16)
	
	## Choose PLOTTITLE
	
	if (ShortenOutputString) ModelRun_string = AbbreviateOutputString(OutputSubDir, modelrun) else ModelRun_string = OutputSubDir
	if (is.null(PLOTTITLE))
	{
		if (PlotType == "Hazards"			) 	PLOTTITLE = paste0("Predicted Hazard ", CountryNames[country+1])	else
		if (PlotType == "Age_Efficacy"		) 	PLOTTITLE = paste0(BaselineSeroStatus, SerotypeString) 				else
		if (PlotType == "Age_Haz_Mult"		) 	PLOTTITLE = paste0("Pred. Age Hazard Multiplier")					else
		if (PlotType == "Age_VacDuration"	) 	PLOTTITLE = paste0("Pred. Age Vac Duration Spline ", BSS_Char)
		if (OutputSubDirInPlotTite) 			PLOTTITLE = paste0(PLOTTITLE, "\n", ModelRun_string)
	}
	
	## Choose axis labels
	if (PlotType == "Hazards"			) XLAB = "Date" else XLAB = "Age (Years)"
	if (PlotType == "Hazards"			) YLAB = "Hazard" 										else
	if (PlotType == "Age_Efficacy"		) YLAB = paste0("Initial transient immunity") 			else
	if (PlotType == "Age_Haz_Mult"		) YLAB = "Hazard Multiplier" 							else
	if (PlotType == "Age_VacDuration"	) YLAB = paste(BaselineSeroStatus, "Vaccine Duration")
	
	## Plots Mean values of Knots
	if (PlotType == "Hazards") xKnotsOrKnotDates = xKnotsDates else xKnotsOrKnotDates = xKnots
	if (PlotType == "Hazards") KnotCol = "blue" else KnotCol = "white" 
		plot(	xKnotsOrKnotDates[1:length(MeanKnots)], MeanKnots, cex = 2, main = PLOTTITLE, ylim = YLIM, xlim = XLIM, cex.main = 1.5, xaxt = "n", xlab = XLAB, ylab = "", col = KnotCol, cex.lab = 2, cex.axis = 1.5, lwd = 4, pch = 4)
	#Plot "end values"
	if (PlotType == "Hazards") 	
	{
		XaxisBeginning 	= seq(		min(as.numeric(xKnotsDates)) - NumDaysAroundKnots	,  	min(as.numeric(xKnotsDates))		, length.out = NoPointsPerKnotGap)
		XaxisEnd		= seq(		max(as.numeric(xKnotsDates))		, 	max(as.numeric(xKnotsDates)) + NumDaysAroundKnots	, length.out = NoPointsPerKnotGap)
		points(XaxisBeginning	, rep(MeanKnots[1]		, length(XaxisBeginning))		, col = "blue")
		points(XaxisEnd			, rep(tail(MeanKnots, 1), length(XaxisBeginning))		, col = "blue")
	}
	mtext(YLAB, side = 2, line = 3, cex = 1.2)
	abline (h = 0)
	if (PlotType == "Hazards")	
		axis(1, xKnotsOrKnotDates, format (xKnotsDates, "%b %y"), cex.axis = 1) else
		axis(1, xKnots, cex.axis = 1)
	
	## Plot each param guess
	param_set = 1
	for (param_set in 1:NumParamSetsToPlot)
	{
		if (param_set %% 1000 == 0) cat (paste(PlotType, BSS_Char, SerotypeString, "Splines : param_set" , param_set, Sys.time(), '\n'))
		yKnots	= Choose_yKnots(PlotType, WhichParamsToPlot[param_set], BaselineSeroStatus, modelrun)
		
		if (PlotType == "Hazards") 	
		{
			points(XaxisBeginning	, rep(yKnots[1]				, NoPointsPerKnotGap)	, col = COLOUR, type = "l")
			points(XaxisEnd			, rep(yKnots[length(xKnots)], NoPointsPerKnotGap)	, col = COLOUR, type = "l")
		}
		
		WhichKnot = 2
		for (WhichKnot in 2:LastKnot)
		{
			x1 = xKnots[WhichKnot-1]						### need x1, y1 for step, line and quad and cubic
			y1 = as.numeric(yKnots[WhichKnot-1])
			
			if (Degree == 3)
			{
				x2 = xKnots[WhichKnot  ]					### need x2, y2 for line and quad and cubic
				y2 = as.numeric(yKnots[WhichKnot  ])
				
				x3 = xKnots[WhichKnot+1]					### need x3, y3 for quad and cubic
				y3 = as.numeric(yKnots[WhichKnot+1])
				
				x4 = xKnots[WhichKnot+2]					### need x4, y4 for cubic
				y4 = as.numeric(yKnots[WhichKnot+2])
				
				Coeffs = FindCoefficients_Cubic(x1,x2,x3,x4, y1,y2,y3,y4)
				a3 = Coeffs[4]; a2 = Coeffs[3]; a1 = Coeffs[2]; a0 = Coeffs[1]
				
				rm(Coeffs)
				
				KnotSpline_Yaxis = a0 + (a1 * KnotSpline_Xaxis_List[[WhichKnot]]) + (a2*KnotSpline_Xaxis_List[[WhichKnot]]^2) + (a3*KnotSpline_Xaxis_List[[WhichKnot]]^3)
				
			} else if (Degree == 2)
			{
				x2 = xKnots[WhichKnot  ]					### need x2, y2 for line and quad and cubic
				y2 = as.numeric(yKnots[WhichKnot  ])
				
				x3 = xKnots[WhichKnot+1]					### need x3, y3 for quad and cubic
				y3 = as.numeric(yKnots[WhichKnot+1])
				
				# find coefficients a0, a1, and a2. Solution for this in notes
				a2 = (y3-y2 - ((x3-x2) * (y1-y2)/(x1-x2) )  ) / (x3*x3 - x2*x2 - (x1+x2)*(x3-x2)) 
				a1 = ((y1-y2)/(x1-x2)) - (x1+x2)* a2
				a0 = y1-(a1*x1)-(a2*x1*x1)
				
				KnotSpline_Yaxis = a0 + (a1 * KnotSpline_Xaxis_List[[WhichKnot]]) + (a2*KnotSpline_Xaxis_List[[WhichKnot]]^2)
				
			} else if (Degree == 1)
			{
				x2 = xKnots[WhichKnot  ]					### need x2, y2 for line and quad
				y2 = as.numeric(yKnots[WhichKnot  ])
				
				a1 = (y2-y1)/(x2-x1)
				a0 = y2 - a1 * x2 
				
				KnotSpline_Yaxis = a0 + (a1 * KnotSpline_Xaxis_List[[WhichKnot]]) 
				
			} else KnotSpline_Yaxis = rep(y1, length(KnotSpline_Xaxis_List[[WhichKnot]])) ### i.e. if Degree == 0, have a step function. 
			
			if (PlotType == "Age_Efficacy" & modelrun$SS_VEs == 1)
			{
				if (modelrun$SSASVE_Additive)
					KnotSpline_Yaxis = KnotSpline_Yaxis + CHAINS[WhichParamsToPlot[param_set], paste0(BSS_Char, "Eff_", serotype)] else ### i.e. same except + rather than *
					KnotSpline_Yaxis = KnotSpline_Yaxis * CHAINS[WhichParamsToPlot[param_set], paste0(BSS_Char, "Eff_", serotype)]
			}
			
			# Guard against impossible values. 
			if (PlotType == "Age_Efficacy")		
				KnotSpline_Yaxis[KnotSpline_Yaxis > 1] = 1 	else 
				KnotSpline_Yaxis[KnotSpline_Yaxis < 0] = 0 	### efficacy cannot be greater than 1. Hazards (multipliers) and durations cannot be less than zero. 
			
			if (PlotType == "Hazards") 	points(KnotSpline_Xaxis_AsDates_List[[WhichKnot]], KnotSpline_Yaxis, col = COLOUR, type = "l") else
				points(KnotSpline_Xaxis_List		[[WhichKnot]], KnotSpline_Yaxis, col = COLOUR, type = "l") 
		}
	}
	
	## plot mean line again (in red)
	yKnots = MeanKnots #### redefine yKnots so code below works - turn this into a function. 
	points(xKnotsOrKnotDates[1:length(yKnots)]	, yKnots, cex = 2, col = COLOUR, lwd = 2)
	
	if (PlotType == "Hazards") 	
	{
		points(XaxisBeginning	, rep(yKnots[1]				, NoPointsPerKnotGap)	, col = "red", type = "l")
		points(XaxisEnd			, rep(yKnots[length(xKnots)], NoPointsPerKnotGap)	, col = "red", type = "l")
	}
	for (WhichKnot in 2:LastKnot)
	{
		x1 = xKnots[WhichKnot-1]					### need x1, y1 for step, line and quad and cubic
		y1 = as.numeric(yKnots[WhichKnot-1])
		
		if (Degree == 3)
		{
			x2 = xKnots[WhichKnot  ]					### need x2, y2 for line and quad and cubic
			y2 = as.numeric(yKnots[WhichKnot  ])
			
			x3 = xKnots[WhichKnot+1]					### need x3, y3 for quad and cubic
			y3 = as.numeric(yKnots[WhichKnot+1])
			
			x4 = xKnots[WhichKnot+2]					### need x4, y4 for cubic
			y4 = as.numeric(yKnots[WhichKnot+2])
			
			Coeffs = FindCoefficients_Cubic(x1,x2,x3,x4, y1,y2,y3,y4)
			a3 = Coeffs[4]; a2 = Coeffs[3]; a1 = Coeffs[2]; a0 = Coeffs[1]
			
			rm(Coeffs)
			
			KnotSpline_Yaxis = a0 + (a1 * KnotSpline_Xaxis_List[[WhichKnot]]) + (a2*KnotSpline_Xaxis_List[[WhichKnot]]^2) + (a3*KnotSpline_Xaxis_List[[WhichKnot]]^3)
			
			
		} else if (Degree == 2)
		{
			x2 = xKnots[WhichKnot  ]					### need x2, y2 for line and quad and cubic
			y2 = as.numeric(yKnots[WhichKnot  ])
			
			x3 = xKnots[WhichKnot+1]					### need x3, y3 for quad and cubic
			y3 = as.numeric(yKnots[WhichKnot+1])
			
			# find coefficients a0, a1, and a2. Solution for this in notes
			a2 = (y3-y2 - ((x3-x2) * (y1-y2)/(x1-x2) )  ) / (x3*x3 - x2*x2 - (x1+x2)*(x3-x2)) 
			a1 = ((y1-y2)/(x1-x2)) - (x1+x2)* a2
			a0 = y1-(a1*x1)-(a2*x1*x1)
			
			KnotSpline_Yaxis = a0 + (a1 * KnotSpline_Xaxis_List[[WhichKnot]]) + (a2*KnotSpline_Xaxis_List[[WhichKnot]]^2)
			
			
		} else if (Degree == 1)
		{
			x2 = xKnots[WhichKnot  ]				### need x2, y2 for line and quad and cubic
			y2 = as.numeric(yKnots[WhichKnot  ])
			
			a1 = (y2-y1)/(x2-x1)
			a0 = y2 - a1 * x2 
			
			KnotSpline_Yaxis = a0 + (a1 * KnotSpline_Xaxis_List[[WhichKnot]]) 
			
		} else KnotSpline_Yaxis = rep(y1, length(KnotSpline_Xaxis_List[[WhichKnot]])) ### i.e. if Degree == 0, have a step function.
		
		# Guard against impossible values. 
		if (PlotType == "Age_Efficacy")		
			KnotSpline_Yaxis[KnotSpline_Yaxis > 1] = 1 	else 
			KnotSpline_Yaxis[KnotSpline_Yaxis < 0] = 0 	### efficacy cannot be greater than 1 - hazards (multipliers) and durations cannot be less than zero. 
		
		if (PlotType == "Hazards") 	points(KnotSpline_Xaxis_AsDates_List[[WhichKnot]], KnotSpline_Yaxis, col = "red", type = "l", cex = 2, lwd = 10) else
			points(KnotSpline_Xaxis_List		[[WhichKnot]], KnotSpline_Yaxis, col = "red", type = "l", cex = 2, lwd = 10) 
	}
	## plot mean knots on top
	if (PlotType == "Hazards") points(xKnotsOrKnotDates[1:length(MeanKnots)]	, MeanKnots, cex = 2, pch = 4, col = "blue", cex.lab = 2, cex.axis = 1.5, lwd = 4)
	
	if (is.null(LegPosChar))
		if ((PlotType == "Age_Efficacy" && BaselineSeroStatus == "SeroNeg") | PlotType == "Age_Haz_Mult") 
			LegPosChar = "bottomright" else 
			LegPosChar = "topleft"
	
	if (IncLegend)	legend(LegPosChar, c("Mean", "95% CrI"), pch = rep(NA,2), col = c("red", "pink"), lty = rep(1,2), lwd = c(7,2), cex = 1.5)
	if (SavePlot) dev.off()

}



PlotHill = function(PlotType = "Age_Efficacy", BaselineSeroStatus = NULL, serotype = 1, modelrun, 
		COLOUR = Alpha_ParamGuessCols, NumParamSetsToPlot = HowManyParamSetsToPlot, OutputSubDirInPlotTite = TRUE, Ages = 2:18,
		SavePlot = TRUE, ShortenOutputString = TRUE, Directory = AgeEffectsDirectory)
{
	if (!any( c("Age_Efficacy", "Age_Haz_Mult", "Age_VacDuration") == PlotType)) 				stop("PlotHill error: PlotType argument not recognized.")
	if (!any(c("SeroNeg", "SeroPos") == BaselineSeroStatus) & !is.null(BaselineSeroStatus)) 	stop("PlotHill error: BaselineSeroStatus argument not recognized")
	if (any(c("red", "blue") == COLOUR)) 														stop("PlotHill error: COLOUR of posterior samples can't be red or blue")
	
	# # # # # Wrong combinations of arguments
	if ("Age_Haz_Mult" == PlotType & !is.null(BaselineSeroStatus)) 									stop("PlotHill error: BaselineSeroStatus != NULL for hazard multiplier")
	if (any(c("Age_Efficacy", "Age_VacDuration") == PlotType) & is.null(BaselineSeroStatus)) 		stop("PlotHill error: BaselineSeroStatus == NULL for Efficacies or duration spline")
	
	## Choose Filename
	if (PlotType == "Age_Efficacy" && BaselineSeroStatus == "SeroNeg") 	BSS_Char = "SNeg" else ## note the double &, works as a circuit breaker if BaselineSeroStatus = NULL
	if (PlotType == "Age_Efficacy" && BaselineSeroStatus == "SeroPos") 	BSS_Char = "SPos" else ## note the double &, works as a circuit breaker if BaselineSeroStatus = NULL
	if (PlotType == "Age_VacDuration") 									BSS_Char = paste(BaselineSeroStatus) else
		BSS_Char = ""
	
	## Create file name	(use filname and date modified to determine whether it is worthwhile making new plot)
	if (ModelRun$SS_VEs					) SerotypeString = paste0("_serotype_", serotype) else 	SerotypeString = ""
	
	if (PlotType == "Age_Efficacy"		) Hill_PNGFILENAME = file.path(Directory, paste0("AgeEffHill_"	, BSS_Char, SerotypeString, ".png")) else 
	if (PlotType == "Age_Haz_Mult"		) Hill_PNGFILENAME = file.path(Directory, paste0("HazMultHill"	, 							".png")) else 
	if (PlotType == "Age_VacDuration"	) Hill_PNGFILENAME = file.path(Directory, paste0("DurationHill_", BSS_Char, 				".png")) 
	
	if (SavePlot)	png(file = Hill_PNGFILENAME  , res = PNG_res, units = "in", width = 7, height= 5)
	
	if (PlotType == "Age_Efficacy"		) PLOTTITLE = paste0("Efficacy Multiplier ", BaselineSeroStatus, SerotypeString) 	else
	if (PlotType == "Age_Haz_Mult"		) PLOTTITLE = paste0("Hazard Multiplier") 											else
	if (PlotType == "Age_VacDuration"	) PLOTTITLE = paste0(BaselineSeroStatus, " Vaccine Duration")
	
	if (ShortenOutputString) ModelRun_string = AbbreviateOutputString(OutputSubDir, modelrun) else ModelRun_string = OutputSubDir
	if (OutputSubDirInPlotTite) PLOTTITLE = paste0(PLOTTITLE, "\n", ModelRun_string)
	
	### Get Hill variable names
	if (PlotType == "Age_Efficacy")
	{
		Power_Char      = "ASVE_Power"		
		Halflife_Char   = "ASVE_HalfLife"	### note cases different between AS_Haz "Halflife"'s and ASVE "HalfLifes"	
		Prop_Char       = "ASVE_Prop"		
		Eff_Char		= paste0(BSS_Char, "Eff_", serotype)
		
		if (!modelrun$AS_VE_Homogeneous)
		{
			Power_Char      = paste0(BSS_Char, "_", Power_Char     )	
			Halflife_Char   = paste0(BSS_Char, "_", Halflife_Char  )	
			Prop_Char       = paste0(BSS_Char, "_", Prop_Char      )	
		}
		
	} else if (PlotType == "Age_Haz_Mult") 
	{
		Power_Char      = "AS_Haz_Power"		
		Halflife_Char   = "AS_Haz_Halflife"	### note cases different between AS_Haz "Halflife"'s and ASVE "HalfLifes"
		Prop_Char       = "AS_Haz_Prop"		
		
	} else if (PlotType == "Age_VacDuration") 
	{
		Power_Char      	= paste0(BaselineSeroStatus, "WaningPower")		
		Halflife_Char   	= paste0(BaselineSeroStatus, "WaningHalflife"	)
		MaxDuration_Char	= paste0(BaselineSeroStatus, "Waning")		
	}
	
	# Choose (mean) variable values
	Power		= MeanParamValues[Power_Char      ]
	HLife 		= MeanParamValues[Halflife_Char   ]
	### Choose proportion of whatever quantity (efficacy, duration, hazard) that is affected by hill function (default is 100%)
	if (PlotType == "Age_Efficacy" || PlotType == "Age_Haz_Mult") 	Proportion = MeanParamValues[Prop_Char] 	else Proportion = 1
	### Choose Multiplier of hill function (for Hill efficay, this is (serotype and) serostatus specific efficacy. For duration, it's the duration, otherwise it's 1. 
	if (PlotType == "Age_Efficacy" )	Multiplier = MeanParamValues[Eff_Char] 			else 
	if (PlotType == "Age_VacDuration")	Multiplier = MeanParamValues[MaxDuration_Char]	else 	Multiplier = 1
	
	### Choose y-axis limits
	if (PlotType == "Age_Efficacy"		) YLIM = c(min(CHAINS[WhichParamsToPlot, Eff_Char]), 1)				else
	if (PlotType == "Age_Haz_Mult"		) YLIM = c(0,1) 													else
	if (PlotType == "Age_VacDuration"	) YLIM = c(0, max(CHAINS[WhichParamsToPlot, MaxDuration_Char]))
	
	### Choose y-axis labels
	if (PlotType == "Age_Efficacy"		) YLAB = paste0(BaselineSeroStatus, "Efficacy", SerotypeString)		else
	if (PlotType == "Age_Haz_Mult"		) YLAB = "Hazard Multiplier" 										else
	if (PlotType == "Age_VacDuration"	) YLAB = paste0(BaselineSeroStatus, " Vaccine Duration")
	
	### plot mean values first
	plot(Ages, Multiplier * ((1 - Proportion) + Proportion * Hill(Ages, Power, HLife)), col = "red", lwd = 3, type = "l", xlab = "Age", ylab = YLAB, 
			ylim = YLIM, main = PLOTTITLE)
	
	for (param in WhichParamsToPlot[1:NumParamSetsToPlot])
	{
		# Choose variable values
		Power		= CHAINS[param, Power_Char      ]
		HLife 		= CHAINS[param, Halflife_Char   ]
		### Choose proportion of whatever quantity (efficacy, duration, hazard) that is affected by hill function (default is 100%)
		if (PlotType == "Age_Efficacy" || PlotType == "Age_Haz_Mult") 	Proportion = CHAINS[param, Prop_Char] 		else Proportion = 1
		### Choose Multiplier of hill function (for Hill efficay, this is (serotype and) serostatus specific efficacy. For duration, it's the duration, otherwise it's 1. 
		if (PlotType == "Age_Efficacy" )	Multiplier = CHAINS[param, Eff_Char] 			else 
		if (PlotType == "Age_VacDuration")	Multiplier = CHAINS[param, MaxDuration_Char]	else 	Multiplier = 1
		
		points(Ages, Multiplier * ((1 - Proportion) + Proportion * Hill(Ages, Power, HLife)), col = COLOUR, type = "l")
	}
	
	### plot mean values again
	Power		= MeanParamValues[Power_Char      ]
	HLife 		= MeanParamValues[Halflife_Char   ]
	### Choose proportion of whatever quantity (efficacy, duration, hazard) that is affected by hill function (default is 100%)
	if (PlotType == "Age_Efficacy" || PlotType == "Age_Haz_Mult") 	Proportion = MeanParamValues[Prop_Char] 	else Proportion = 1
	### Choose Multiplier of hill function (for Hill efficay, this is (serotype and) serostatus specific efficacy. For duration, it's the duration, otherwise it's 1. 
	if (PlotType == "Age_Efficacy" )	Multiplier = MeanParamValues[Eff_Char] 			else 
	if (PlotType == "Age_VacDuration")	Multiplier = MeanParamValues[MaxDuration_Char]	else 	Multiplier = 1
	points(Ages, Multiplier * ((1 - Proportion) + Proportion * Hill(Ages, Power, HLife)), col = "red", lwd = 3, type = "l")
	
	if (SavePlot)	dev.off()
	
}









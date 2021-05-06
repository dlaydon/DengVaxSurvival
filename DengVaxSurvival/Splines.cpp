
#include "HeaderAndSwitches.h"

DType Hill						(DType TimePostDose, DType HillCoeff, DType Halflife)
{
	DType ToverH_RaisedToPower = pow(TimePostDose/Halflife, HillCoeff);

	return(ToverH_RaisedToPower / (1 + ToverH_RaisedToPower));
}
int FindLastPolyBefore_t		(DType t, DType *xKnots, const int &PolynomialsPerCountry)
{
	int WhichPolynomial = PolynomialsPerCountry - 1;  
	while (xKnots[WhichPolynomial] > t) WhichPolynomial--;

	return WhichPolynomial;
}
void CoefficientsFromKnots		(DType* xKnots, DType* yKnots, int FirstKnot, DType *Coeffs)
{
	//// Note that in PARAMS structure, xKnots/yKnots are 2D arrays, with the first dimension being country. To make things more general just use the knots and input PARAMS.xKnots[country]
	DType x0, x1, x2, y0, y1, y2;

	x0 = xKnots[FirstKnot];
	x1 = xKnots[FirstKnot + 1];
	x2 = xKnots[FirstKnot + 2];

	y0 = yKnots[FirstKnot];
	y1 = yKnots[FirstKnot + 1];
	y2 = yKnots[FirstKnot + 2];

	Coeffs[2] = (y2 - y1 - ((x2 - x1) * (y0 - y1) / (x0 - x1))) / ((x2*x2) - (x1*x1) - (x0 + x1)*(x2 - x1));
	Coeffs[1] = ((y0 - y1) / (x0 - x1)) - (x0 + x1)* Coeffs[2];
	Coeffs[0] = y0 - (Coeffs[1] * x0) - (Coeffs[2] * (x0*x0)); 
}
void CoefficientsFromKnots_Line	(DType* xKnots, DType* yKnots, int FirstKnot, DType *Coeffs)
{
	//// Note that in PARAMS structure, xKnots/yKnots are 2D arrays, with the first dimension being country. To make things more general just use the knots and input PARAMS.xKnots[country]
	DType x0, x1, y0, y1;

	x0 = xKnots[FirstKnot		];
	x1 = xKnots[FirstKnot + 1	];

	y0 = yKnots[FirstKnot		];
	y1 = yKnots[FirstKnot + 1	];

	Coeffs[1] = (y1 - y0) / (x1 - x0);
	Coeffs[0] = y1 - Coeffs[1] * x1; 
}
void CoefficientsFromKnots_Cubic(DType* xKnots, DType* yKnots, int FirstKnot, DType *Coeffs)
{
	//// Note that in PARAMS structure, xKnots/yKnots are 2D arrays, with the first dimension being country. To make things more general just use the knots and input PARAMS.xKnots[country]
	DType x0, x1, x2, x3, y0, y1, y2, y3;

	x0 = xKnots[FirstKnot];
	x1 = xKnots[FirstKnot + 1];
	x2 = xKnots[FirstKnot + 2];
	x3 = xKnots[FirstKnot + 3];

	y0 = yKnots[FirstKnot];
	y1 = yKnots[FirstKnot + 1];
	y2 = yKnots[FirstKnot + 2];
	y3 = yKnots[FirstKnot + 3];

	//// Definitions: Makes rest of code MUCH easier. 
	DType L_20, L_30, Q_20, Q_30, C_20, C_30; //// L = linear, Q = quadratic and C = cubic terms from notes. Annoyingly you've changed the letters but the ones here make more sense and the pattern is easier to see. 

	L_20 = (x1 - x0) * (y2				- y0)				- (x2 - x0) * (y1			- y0);				//// was labelled K20 in notes
	L_30 = (x1 - x0) * (y3				- y0)				- (x3 - x0) * (y1			- y0);				//// was labelled K30 in notes
																												  
	Q_20 = (x1 - x0) * (x2 * x2			- x0 * x0)			- (x2 - x0) * (x1 * x1		- x0 * x0);			//// was labelled D20 in notes
	Q_30 = (x1 - x0) * (x3 * x3			- x0 * x0)			- (x3 - x0) * (x1 * x1		- x0 * x0);			//// was labelled D30 in notes
																												  
	C_20 = (x1 - x0) * (x2 * x2 * x2	- x0 * x0* x0)		- (x2 - x0) * (x1 * x1* x1	- x0 * x0 * x0);	//// was labelled Q20 in notes
	C_30 = (x1 - x0) * (x3 * x3 * x3	- x0 * x0* x0)		- (x3 - x0) * (x1 * x1* x1	- x0 * x0 * x0);	//// was labelled Q30 in notes


	DType a3 = (L_30*Q_20 - L_20*Q_30) / (C_30*Q_20 - C_20*Q_30); 
	DType a2 = (L_20 - a3 * C_20) / Q_20;	//// these should be equivalent.
	//DType a2 = (L_30 - a3 * C_30) / Q_30;	//// these should be equivalent. 
	DType a1 = (y1-y0 - a2*(x1*x1 - x0*x0) -a3*(x1*x1*x1 - x0*x0*x0))/(x1-x0) ;
	DType a0 = y0 - a1*x0 - a2*x0*x0 - a3* x0*x0*x0; 

	Coeffs[0] = a0;
	Coeffs[1] = a1;
	Coeffs[2] = a2;
	Coeffs[3] = a3;
}

///// Because there are different limits on spline (e.g. efficacies can be negative, param ranges allowing, but multipliers of the hazard cannot, and efficacies cannot exceed 1, but durations (say) can), you have essentially one spline, but put the limits for the others on wrapper functions, below.
inline DType Spline				(DType xValue, DType* xKnots, DType* yKnots, DType **Coeffs, int KnotsPerSpline, int PolynomialsPerSpline)
{
	DType yValue;

	//// If xValue is less than or equal to minimum xKnot, spline takes constant value of the corresponding yKnot. 
		 if (xValue <= xKnots[0]					)	yValue = yKnots[0];
	else if (xValue >= xKnots[KnotsPerSpline - 1]	)	yValue = yKnots[KnotsPerSpline - 1];
	else
	{
		int WhichPolynomial = FindLastPolyBefore_t(xValue, xKnots, PolynomialsPerSpline);

		yValue = Coeffs[WhichPolynomial][0] + Coeffs[WhichPolynomial][1] * xValue + Coeffs[WhichPolynomial][2] * xValue*xValue;
	}
	//return ((yValue < 0) ? 0 : yValue); //// efficacy splines may need to be negative - instead put a guard in the waning and baseline hazard splines
	return yValue;
}
inline DType Line				(DType xValue, DType* xKnots, DType* yKnots, DType **Coeffs, int KnotsPerSpline, int PolynomialsPerSpline)
{
	DType yValue;

	//// If xValue is less than or equal to minimum xKnot, spline takes constant value of the corresponding yKnot. 
		 if (xValue <= xKnots[0]					)	yValue = yKnots[0];
	else if (xValue >= xKnots[KnotsPerSpline - 1]	)	yValue = yKnots[KnotsPerSpline - 1];
	else
	{
		int WhichPolynomial = FindLastPolyBefore_t(xValue, xKnots, PolynomialsPerSpline);

		yValue = Coeffs[WhichPolynomial][0] + Coeffs[WhichPolynomial][1] * xValue;
	}
	//return ((yValue < 0) ? 0 : yValue); //// efficacy splines may need to be negative - instead put a guard in the waning and baseline hazard splines
	return yValue;
}
inline DType Step				(DType xValue, DType* xKnots, DType* yKnots, DType **Coeffs, int KnotsPerSpline, int PolynomialsPerSpline)
{
	DType yValue;

	//// If xValue is less than or equal to minimum xKnot, spline takes constant value of the corresponding yKnot. 
		 if (xValue <= xKnots[0]					)	yValue = yKnots[0];
	else if (xValue >= xKnots[KnotsPerSpline - 1]	)	yValue = yKnots[KnotsPerSpline - 1];
	else
	{
		int WhichPolynomial = FindLastPolyBefore_t(xValue, xKnots, PolynomialsPerSpline);

		yValue = yKnots[WhichPolynomial];
	}
	//return ((yValue < 0) ? 0 : yValue); //// efficacy splines may need to be negative - instead put a guard in the waning and baseline hazard splines
	return yValue;
}
inline DType Cubic				(DType xValue, DType* xKnots, DType* yKnots, DType **Coeffs, int KnotsPerSpline, int PolynomialsPerSpline)
{
	DType yValue;

	//// If xValue is less than or equal to minimum xKnot, spline takes constant value of the corresponding yKnot. 
		 if (xValue <= xKnots[0]					)	yValue = yKnots[0];
	else if (xValue >= xKnots[KnotsPerSpline - 1]	)	yValue = yKnots[KnotsPerSpline - 1];
	else
	{
		int WhichPolynomial = FindLastPolyBefore_t(xValue, xKnots, PolynomialsPerSpline);

		yValue = Coeffs[WhichPolynomial][0] + Coeffs[WhichPolynomial][1]*xValue + Coeffs[WhichPolynomial][2]*xValue*xValue + Coeffs[WhichPolynomial][3]*xValue*xValue*xValue;
	}
	return yValue;
}

DType BaselineHazard			(DType xValue, const int &country, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	DType yValue = Spline(xValue, PARAMS.xKnots[country], PARAMS.yKnots[country], PARAMS.SplineCoeffs[country], HOUSE.KnotsPerCountry, HOUSE.PolynomialsPerCountry);

	return ((yValue < 0) ? 0 : yValue);
}
DType WaningDurationSpline		(int AgeInYears, int BaselineSeroStatus, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	DType yValue; 
	if (HOUSE.AS_Waning==Age_Option::SPLINE)
		yValue = Spline((DType)AgeInYears,	PARAMS.Age_DurationRate_xKnots		[BaselineSeroStatus],
											PARAMS.WaningParams					[BaselineSeroStatus],
											PARAMS.Age_SplineCoeffs_DurationRate[BaselineSeroStatus],
											HOUSE.KnotsPerSpline_WaningDuration,
											HOUSE.PolynomialsPerSpline_WaningDuration				);
	else if (HOUSE.AS_Waning == Age_Option::SPLINE_LINE)
		yValue = Line((DType)AgeInYears,	PARAMS.Age_DurationRate_xKnots		[BaselineSeroStatus],
											PARAMS.WaningParams					[BaselineSeroStatus],
											PARAMS.Age_SplineCoeffs_DurationRate[BaselineSeroStatus],
											HOUSE.KnotsPerSpline_WaningDuration,
											HOUSE.PolynomialsPerSpline_WaningDuration				);
	else if (HOUSE.AS_Waning == Age_Option::SPLINE_STEP)
		yValue = Step((DType)AgeInYears, PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus],
			PARAMS.WaningParams[BaselineSeroStatus],
			PARAMS.Age_SplineCoeffs_DurationRate[BaselineSeroStatus],
			HOUSE.KnotsPerSpline_WaningDuration,
			HOUSE.PolynomialsPerSpline_WaningDuration);
	else if (HOUSE.AS_Waning == Age_Option::CUBIC)
		yValue = Cubic((DType)AgeInYears,	PARAMS.Age_DurationRate_xKnots		[BaselineSeroStatus],
											PARAMS.WaningParams					[BaselineSeroStatus],
											PARAMS.Age_SplineCoeffs_DurationRate[BaselineSeroStatus],
											HOUSE.KnotsPerSpline_WaningDuration,
											HOUSE.PolynomialsPerSpline_WaningDuration				);
	else std::cerr << "WaningDurationSpline ERROR: Age_Option not recognized" << endl;
	return ((yValue < 0) ? 0 : yValue);
}
DType AS_HazMultSpline			(int AgeInYears, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	DType yValue;
	if (HOUSE.AS_Haz == Age_Option::SPLINE)
		yValue = Spline((DType)AgeInYears,	PARAMS.Age_HazMult_xKnots					,
											PARAMS.ASHaz_Params							,
											PARAMS.Age_SplineCoeffs_HazMult				, 
											HOUSE.KnotsPerSpline_HazMultiplier			, 
											HOUSE.PolynomialsPerSpline_HazMultiplier	);
	else if (HOUSE.AS_Haz == Age_Option::SPLINE_LINE)
		yValue = Line((DType)AgeInYears,	PARAMS.Age_HazMult_xKnots					,
											PARAMS.ASHaz_Params							,
											PARAMS.Age_SplineCoeffs_HazMult				, 
											HOUSE.KnotsPerSpline_HazMultiplier			, 
											HOUSE.PolynomialsPerSpline_HazMultiplier	);
	else if (HOUSE.AS_Haz == Age_Option::SPLINE_STEP)
		yValue = Step((DType)AgeInYears,	PARAMS.Age_HazMult_xKnots					,
											PARAMS.ASHaz_Params							,
											PARAMS.Age_SplineCoeffs_HazMult				, 
											HOUSE.KnotsPerSpline_HazMultiplier			, 
											HOUSE.PolynomialsPerSpline_HazMultiplier	);
	else if (HOUSE.AS_Haz == Age_Option::CUBIC)
		yValue = Cubic((DType)AgeInYears,	PARAMS.Age_HazMult_xKnots					,
											PARAMS.ASHaz_Params							,
											PARAMS.Age_SplineCoeffs_HazMult				, 
											HOUSE.KnotsPerSpline_HazMultiplier			, 
											HOUSE.PolynomialsPerSpline_HazMultiplier	);
	else std::cerr << "AS_HazMultSpline ERROR: Age_Option not recognized" << endl;  
	return ((yValue < 0) ? 0 : yValue);
}
DType EfficacyMultiplierSpline	(DType age, int BaselineSeroStatus, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	DType Efficacy;
	//// Depending on value of age option, call different function
	if (HOUSE.ASVE == Age_Option::SPLINE)
		Efficacy = Spline(age, PARAMS.Age_Effs_xKnots[BaselineSeroStatus],
			PARAMS.ASVE_Params[BaselineSeroStatus],
			PARAMS.Age_SplineCoeffs_Effs[BaselineSeroStatus],
			HOUSE.KnotsPerSpline_EffMultiplier,
			HOUSE.PolynomialsPerSpline_EffMultiplier);
	else if (HOUSE.ASVE == Age_Option::SPLINE_LINE)
		Efficacy = Line(age, PARAMS.Age_Effs_xKnots[BaselineSeroStatus],
			PARAMS.ASVE_Params[BaselineSeroStatus],
			PARAMS.Age_SplineCoeffs_Effs[BaselineSeroStatus],
			HOUSE.KnotsPerSpline_EffMultiplier,
			HOUSE.PolynomialsPerSpline_EffMultiplier);
	else if (HOUSE.ASVE == Age_Option::SPLINE_STEP)
		Efficacy = Step(age, PARAMS.Age_Effs_xKnots[BaselineSeroStatus],
			PARAMS.ASVE_Params[BaselineSeroStatus],
			PARAMS.Age_SplineCoeffs_Effs[BaselineSeroStatus],
			HOUSE.KnotsPerSpline_EffMultiplier,
			HOUSE.PolynomialsPerSpline_EffMultiplier);
	else if (HOUSE.ASVE == Age_Option::CUBIC)
		Efficacy = Cubic(age, PARAMS.Age_Effs_xKnots[BaselineSeroStatus],
			PARAMS.ASVE_Params[BaselineSeroStatus],
			PARAMS.Age_SplineCoeffs_Effs[BaselineSeroStatus],
			HOUSE.KnotsPerSpline_EffMultiplier,
			HOUSE.PolynomialsPerSpline_EffMultiplier);
	else std::cerr << "AgeEfficacySpline ERROR: Age_Option not recognized" << endl;

	//// Make sure not lower than lower bound (and not higher than 1)
		 if (BaselineSeroStatus == SeroNeg && Efficacy < HOUSE.IntialParRanges.SNegEff_1[LowerBound]) return HOUSE.IntialParRanges.SNegEff_1[LowerBound];
	else if (BaselineSeroStatus == SeroPos && Efficacy < HOUSE.IntialParRanges.SPosEff_1[LowerBound]) return HOUSE.IntialParRanges.SPosEff_1[LowerBound];
	else if (Efficacy > 1) return (DType)1; /// and not higher than 100% (either sero)
	else return Efficacy; 
}

DType IntBaseHaz				(int Lower_t, int Upper_t, const int &country, const Params_Struct &PARAMS)	
{
	DType IntBaseHaz = PARAMS.IntBaseHazLookUp[country][Upper_t] - PARAMS.IntBaseHazLookUp[country][Lower_t];
	return IntBaseHaz;
}
DType VaccineHazard				(int CalendarDay, int Dose_1_Day, int Dose_2_Day, int Dose_3_Day, int AgeInYears, char BaselineSeroStatus, const int &country, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	DType VaccineHazard = 0; 
	if (HOUSE.SingleOrMultiDose == MULTI_DOSE)
			 if	(CalendarDay < Dose_2_Day)	VaccineHazard = PARAMS.BaselineHazardValues[country][CalendarDay] * PARAMS.WaningMults[AgeInYears][BaselineSeroStatus][CalendarDay - Dose_1_Day];
		else if (CalendarDay < Dose_3_Day)	VaccineHazard = PARAMS.BaselineHazardValues[country][CalendarDay] * PARAMS.WaningMults[AgeInYears][BaselineSeroStatus][CalendarDay - Dose_2_Day];
		else								VaccineHazard = PARAMS.BaselineHazardValues[country][CalendarDay] * PARAMS.WaningMults[AgeInYears][BaselineSeroStatus][CalendarDay - Dose_3_Day];
	else									VaccineHazard = PARAMS.BaselineHazardValues[country][CalendarDay] * PARAMS.WaningMults[AgeInYears][BaselineSeroStatus][CalendarDay - Dose_1_Day];

	return VaccineHazard; 
}
DType IntVacHaz_SingleDay		(int CalendarDay, int Dose_1_Day, int Dose_2_Day, int Dose_3_Day, int AgeInYears, char BaselineSeroStatus, const int &country, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE) //// Used for survival curve where you cumulate the hazard by day, saves loads of loops.
{
	DType IntVacHaz = HOUSE.TimeInterval * VaccineHazard(CalendarDay, Dose_1_Day,  Dose_2_Day,  Dose_3_Day, AgeInYears, BaselineSeroStatus, country, PARAMS, HOUSE);
	return IntVacHaz;
}
DType IntVacHaz					(int Lower_t, int Upper_t, int Dose_1_Day, int Dose_2_Day, int Dose_3_Day, int AgeInYears, char BaselineSeroStatus, const int &country, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	DType IntVacHaz = 0;
	if (Upper_t >= Lower_t) 
		if (HOUSE.SingleOrMultiDose == SINGLE_DOSE) 		
			for (int CalendarDay = Lower_t; CalendarDay < Upper_t; CalendarDay++)
				IntVacHaz += PARAMS.BaselineHazardValues[country][CalendarDay] * PARAMS.WaningMults[AgeInYears][BaselineSeroStatus][CalendarDay - Dose_1_Day];
		else 
			if (Lower_t < Dose_2_Day)  
				if (Upper_t < Dose_2_Day)
					for (int CalendarDay = Lower_t; CalendarDay < Upper_t; CalendarDay++)
						IntVacHaz += PARAMS.BaselineHazardValues[country][CalendarDay] * PARAMS.WaningMults[AgeInYears][BaselineSeroStatus][CalendarDay - Dose_1_Day];
				else		// i.e. Upper_t >= Dose_2_Date
					if (Upper_t < Dose_3_Day)
					{
						for (int CalendarDay = Lower_t; CalendarDay < Dose_2_Day; CalendarDay++)	// add in instantaneous hazard days between 1st dose and 2nd dose
							IntVacHaz += PARAMS.BaselineHazardValues[country][CalendarDay] * PARAMS.WaningMults[AgeInYears][BaselineSeroStatus][CalendarDay - Dose_1_Day];
						for (int CalendarDay = Dose_2_Day; CalendarDay < Upper_t; CalendarDay++)	// add in instantaneous hazard days between 2nd dose and Upper_t
							IntVacHaz += PARAMS.BaselineHazardValues[country][CalendarDay] * PARAMS.WaningMults[AgeInYears][BaselineSeroStatus][CalendarDay - Dose_2_Day];
					}
					else 	// i.e. Upper_t >= Dose_3_Date
					{
						for (int CalendarDay = Lower_t; CalendarDay < Dose_2_Day; CalendarDay++)	// add in instantaneous hazard days between 1st dose and 2nd dose
							IntVacHaz += PARAMS.BaselineHazardValues[country][CalendarDay] * PARAMS.WaningMults[AgeInYears][BaselineSeroStatus][CalendarDay - Dose_1_Day];
						for (int CalendarDay = Dose_2_Day; CalendarDay < Dose_3_Day; CalendarDay++)	// add in instantaneous hazard days between 2nd dose and 3rd dose
							IntVacHaz += PARAMS.BaselineHazardValues[country][CalendarDay] * PARAMS.WaningMults[AgeInYears][BaselineSeroStatus][CalendarDay - Dose_2_Day];
						for (int CalendarDay = Dose_3_Day; CalendarDay < Upper_t; CalendarDay++)	// add in instantaneous hazard days between 3rd dose and Upper_t
							IntVacHaz += PARAMS.BaselineHazardValues[country][CalendarDay] * PARAMS.WaningMults[AgeInYears][BaselineSeroStatus][CalendarDay - Dose_3_Day];
					}
			else			// i.e. Lower_t >= Dose_2_Date (therefore Upper_t >= Dose_2_date).
				if (Lower_t < Dose_3_Day)	// lower_t between dose 2 and 3
					if (Upper_t < Dose_3_Day)	// upper_t between dose 2 and 3
						for (int CalendarDay = Lower_t; CalendarDay < Upper_t; CalendarDay++)
							IntVacHaz += PARAMS.BaselineHazardValues[country][CalendarDay] * PARAMS.WaningMults[AgeInYears][BaselineSeroStatus][CalendarDay - Dose_2_Day];
					else					// lower_t between dose 2 and 3, upper_t on or after dose 3 
					{	
						for (int CalendarDay = Lower_t; CalendarDay < Dose_3_Day; CalendarDay++)
							IntVacHaz += PARAMS.BaselineHazardValues[country][CalendarDay] * PARAMS.WaningMults[AgeInYears][BaselineSeroStatus][CalendarDay - Dose_2_Day];
						for (int CalendarDay = Dose_3_Day; CalendarDay < Upper_t; CalendarDay++)	// add in instantaneous hazard days between 3rd dose and Upper_t
							IntVacHaz += PARAMS.BaselineHazardValues[country][CalendarDay] * PARAMS.WaningMults[AgeInYears][BaselineSeroStatus][CalendarDay - Dose_3_Day];
					}
				else						// i.e. Lower_t >= Dose_3_date (therefore Upper_t >= Dose_3_Date)
					for (int CalendarDay = Lower_t; CalendarDay < Upper_t; CalendarDay++)
						IntVacHaz += PARAMS.BaselineHazardValues[country][CalendarDay] * PARAMS.WaningMults[AgeInYears][BaselineSeroStatus][CalendarDay - Dose_3_Day];

	/// ultiply answer by TimeInterval delta t (if not doing analytical solution).
	IntVacHaz *= HOUSE.TimeInterval;
	return IntVacHaz;
}



#include "HeaderAndSwitches.h"
#include "Probability.h"
#include "Splines.h"
#include "Initialize.h"
#include "WriteOutput.h"
#include "SurvivalCurves.h"
#include "ConsolePrinting.h"
#include "Augmentation.h"
#include "StructureDefs.h"

void DoOrDontUpdateLowerTail				(DType &NewElement, DType *LowerTail, DType &MaxValueOfLower, int &MaxValueOfLower_index, int &NoElementsInTail)
{
	if (NewElement < MaxValueOfLower)	 // i.e. element should be in Lower Tail
	{
		//// replace maximum element of Lower Tail...
		LowerTail[MaxValueOfLower_index] = NewElement;

		//// ... and find new MaxValueOfLower and MaxValueOfLower_index
		MaxValueOfLower = -DBL_MAX; MaxValueOfLower_index = 0;
		for (int i = 0; i < NoElementsInTail; i++) if (LowerTail[i] >= MaxValueOfLower) { MaxValueOfLower = LowerTail[i]; MaxValueOfLower_index = i; }
	}
}
void DoOrDontUpdateUpperTail				(DType &NewElement, DType *UpperTail, DType &MinValueOfUpper, int &MinValueOfUpper_index, int &NoElementsInTail)
{
	if (NewElement > MinValueOfUpper)	 // i.e. element should be in UpperTail
	{
		//// replace minimum element of Upper Tail...
		UpperTail[MinValueOfUpper_index] = NewElement;

		//// ... and find new MinValueOfUpper and MinValueOfUpper_index
		MinValueOfUpper = DBL_MAX; MinValueOfUpper_index = 0;
		for (int i = 0; i < NoElementsInTail; i++) if (UpperTail[i] <= MinValueOfUpper) { MinValueOfUpper = UpperTail[i]; MinValueOfUpper_index = i; }
	}
}

void DivideSummedProbabilies				(DType **PosteriorSample_SurvivalTable, DType **Strata_Sizes, int NoDaysOfFollowUp, int NumStrata, Survival_Struct &SURVIVE, const Housekeeping_Struct &HOUSE)
{
	/*
		Function is for any survival table: i) mild; ii) severe or iii) either severity. For full survival table or passive.
		Definitions:
		PosteriorSample_SurvivalTable	is the survival table for this parameter set.
		NoDaysOfFollowUp				Normally this will be SURVIVE.NoDaysOfFollowUp + 1. However if doing passive phase survival curve then will be SURVIVE.NoPassiveDaysFollowUp
		Strata_Sizes					Number of people in each strata for each daypostdose. 
	*/

	//// Divide PosteriorSample_SurvivalTableb by approrpriate denominator. 
#pragma omp parallel for schedule(static,1)
	for (int th = 0; th < HOUSE.max_threads; th++)
		for (int category = th; category < NumStrata; category += HOUSE.max_threads)
			for (int daypostdose = 0; daypostdose < NoDaysOfFollowUp; daypostdose++)
				if (!(Strata_Sizes[category][daypostdose] == 0))	PosteriorSample_SurvivalTable[category][daypostdose] /= Strata_Sizes[category][daypostdose];
}
void DivideSummedProbabilies				(SCurves &SetOfSCurves, int NumDays, int NumStrata, Survival_Struct &SURVIVE, const Housekeeping_Struct &HOUSE)
{
	for (int DiseaseSeverity = 0; DiseaseSeverity < SURVIVE.HowManyDiseaseSeverities; DiseaseSeverity++)
		DivideSummedProbabilies(SetOfSCurves.PostSample_SurvivalTables[DiseaseSeverity], SetOfSCurves.Strata_Sizes, NumDays, NumStrata, SURVIVE, HOUSE);
}
void Process_LTFU_SCurves					(DType **PosteriorSample_SurvivalTable, int NoDaysOfFollowUp, int NumStrata, Survival_Struct &SURVIVE)
{
	// to avoid zeros associated with given strata having nobody in follow up post a certain follow up day, redefine so that it takes the previous day's survival probability. 
	// Can't parallelize this though (must do day t, then t+1 and can't be out of order), and  needs to be done before adding to mean and tails. 
	for (int daypostdose = 1; daypostdose < NoDaysOfFollowUp; daypostdose++)
		for (int category = 0; category < NumStrata; category++)
			if (PosteriorSample_SurvivalTable[category][daypostdose] == 0)
				PosteriorSample_SurvivalTable[category][daypostdose] = PosteriorSample_SurvivalTable[category][daypostdose - 1];
}
void AdjustTailsAndAddToFinalMeanPosterior	(DType **PosteriorSample_SurvivalTable, DType **Mean_PosteriorSurvivalCurves,
	DType ****CurveTails, DType *** MaxMinValsForLowerUpperSurvivaltails, int *** MaxMin_Indices_ForLowerUpperSurvivalTails, int NoDaysOfFollowUp, int NumStrata,
	Survival_Struct &SURVIVE, const Housekeeping_Struct &HOUSE)
{
	/*

		This function is for any survival table: i) mild; ii) severe or iii) either severity. For full survival table or passive.

		Definitions

		PosteriorSample_SurvivalTable					is the survival table for this parameter set.
		Mean_PosteriorSurvivalCurves					is the outputted mean posterior survival table, to which you add at each iteration. (can be i) mild; ii) severe or iii) either severity).
		SurvivalCurveTails								store the tail values for upper and lower CrIs (again, either i) mild; ii) severe or iii) either severity)
		MaxMinValsForLowerUpperSurvivaltails			stores the minimum and maximum values of the lower and upper tails (again, either i) mild; ii) severe or iii) either severity)
		MaxMin_Indices_ForLowerUpperSurvivalTails		stores the indices of the minimum and maximum values of the lower and upper tails (again, either i) mild; ii) severe or iii) either severity)
		NoDaysOfFollowUp								Normally this will be SURVIVE.NoDaysOfFollowUp + 1. However if doing passive phase survival curve then will be SURVIVE.NoPassiveDaysFollowUp

	*/
	///// add to mean and tails. 
#pragma omp parallel for schedule(static,1)
	for (int th = 0; th < HOUSE.max_threads; th++)
		for (int category = th; category < NumStrata; category += HOUSE.max_threads)
			for (int daypostdose = 0; daypostdose < NoDaysOfFollowUp; daypostdose++)
			{
				//// add to mean
				Mean_PosteriorSurvivalCurves[category][daypostdose] += PosteriorSample_SurvivalTable[category][daypostdose];

				//// update lower tail
				DoOrDontUpdateLowerTail(PosteriorSample_SurvivalTable[category][daypostdose], CurveTails[LowerTail_Index][category][daypostdose],
					MaxMinValsForLowerUpperSurvivaltails[LowerTail_Index][category][daypostdose], MaxMin_Indices_ForLowerUpperSurvivalTails[LowerTail_Index][category][daypostdose], SURVIVE.NumElementsOutside_CrI_Tails);

				//// update upper tail. 
				DoOrDontUpdateUpperTail(PosteriorSample_SurvivalTable[category][daypostdose], CurveTails[UpperTail_Index][category][daypostdose],
					MaxMinValsForLowerUpperSurvivaltails[UpperTail_Index][category][daypostdose], MaxMin_Indices_ForLowerUpperSurvivalTails[UpperTail_Index][category][daypostdose], SURVIVE.NumElementsOutside_CrI_Tails);
			}
}
void AdjustTailsAndAddToFinalMeanPosterior	(SCurves &SetOfSCurves, int NumDays, int NumStrata,	Survival_Struct &SURVIVE, const Housekeeping_Struct &HOUSE)
{
	//// Wrapper of AdjustTailsAndAddToFinalMeanPosterior. Applies function to all "DiseaseSeverities". Used for HRs (hazard ratios), as SC_WT and SC_PP done in different overload. 
	for (int DiseaseSeverity = 0; DiseaseSeverity < SURVIVE.HowManyDiseaseSeverities; DiseaseSeverity++)
		AdjustTailsAndAddToFinalMeanPosterior(SetOfSCurves.PostSample_SurvivalTables[DiseaseSeverity], SetOfSCurves.FinalPosteriorSurvivalCurves[DiseaseSeverity][MEAN_POST], 
			SetOfSCurves.CurveTails[DiseaseSeverity], SetOfSCurves.MaxMin_Vals_ForTails[DiseaseSeverity], SetOfSCurves.MaxMin_Indices_ForTails[DiseaseSeverity], NumDays, NumStrata, SURVIVE, HOUSE);
}
void AddToSurvivalTailsAndMean				(DType **PosteriorSample_SurvivalTable, DType **Mean_PosteriorSurvivalCurves, DType **Strata_Sizes, 
	DType ****CurveTails					 ,  DType *** MaxMinValsForLowerUpperSurvivaltails, int *** MaxMin_Indices_ForLowerUpperSurvivalTails, int NoDaysOfFollowUp, int NumStrata,
	Survival_Struct &SURVIVE, const Housekeeping_Struct &HOUSE)
{
	/*		
	
	This function is for any survival table: i) mild; ii) severe or iii) either severity. For full survival table or passive. It uses DivideSummedProbabilies and AdjustTailsAndAddToFinalMeanPosterior
	Definitions:
	PosteriorSample_SurvivalTable					is the survival table for this parameter set. 
	Mean_PosteriorSurvivalCurves					is the outputted mean posterior survival table, to which you add at each iteration. (can be i) mild; ii) severe or iii) either severity). 
	SurvivalCurveTails								store the tail values for upper and lower CrIs (again, either i) mild; ii) severe or iii) either severity)
	MaxMinValsForLowerUpperSurvivaltails			stores the minimum and maximum values of the lower and upper tails (again, either i) mild; ii) severe or iii) either severity)
	MaxMin_Indices_ForLowerUpperSurvivalTails		stores the indices of the minimum and maximum values of the lower and upper tails (again, either i) mild; ii) severe or iii) either severity)
	NoDaysOfFollowUp								Normally this will be SURVIVE.NoDaysOfFollowUp + 1. However if doing passive phase survival curve then will be SURVIVE.NoPassiveDaysFollowUp

	*/
	//// Divide PosteriorSample_SurvivalTableb by approrpriate denominator. 
	DivideSummedProbabilies(PosteriorSample_SurvivalTable, Strata_Sizes, NoDaysOfFollowUp, NumStrata, SURVIVE, HOUSE);

	//// to avoid zeros associated with given strata having nobody in follow up post a certain follow up day, redefine so that it takes the previous days survival probability. 
	///// Can't parallelize this though (must do day t, then t+1 and can't be out of order), and it needs to be done before adding to mean and tails.
	if (HOUSE.LTFU_SurvivalCurves)
		Process_LTFU_SCurves(PosteriorSample_SurvivalTable, NoDaysOfFollowUp, NumStrata, SURVIVE);

	///// add to Final Mean Posterior and tails. 
	AdjustTailsAndAddToFinalMeanPosterior(PosteriorSample_SurvivalTable, Mean_PosteriorSurvivalCurves, CurveTails, 
		MaxMinValsForLowerUpperSurvivaltails, MaxMin_Indices_ForLowerUpperSurvivalTails, 
		NoDaysOfFollowUp, NumStrata, SURVIVE, HOUSE);
}
void AddToSurvivalTailsAndMean				(SCurves &SetOfSCurves, int NumDays, int NumStrata, Survival_Struct &SURVIVE, const Housekeeping_Struct &HOUSE)
{
	///// Wrapper of earlier AddToSurvivalTailsAndMean function. Applies function to all "DiseaseSeverities"
	for (int DiseaseSeverity = 0; DiseaseSeverity < SURVIVE.HowManyDiseaseSeverities; DiseaseSeverity++)
		AddToSurvivalTailsAndMean(SetOfSCurves.PostSample_SurvivalTables[DiseaseSeverity], SetOfSCurves.FinalPosteriorSurvivalCurves[DiseaseSeverity][MEAN_POST], SetOfSCurves.Strata_Sizes,
			SetOfSCurves.CurveTails[DiseaseSeverity], SetOfSCurves.MaxMin_Vals_ForTails[DiseaseSeverity], SetOfSCurves.MaxMin_Indices_ForTails[DiseaseSeverity], NumDays, NumStrata, SURVIVE, HOUSE);
}
void AddToSurvivalTailsAndMean				(Survival_Struct &SURVIVE, const Housekeeping_Struct &HOUSE)
{
	///// Function is a wrapper of (a wrapper of) earlier AddToSurvivalTailsAndMean function. Applies function to multiple sets of SCurves (i.e. whole trial and passive phase). 
	AddToSurvivalTailsAndMean(SURVIVE.SC_WT, SURVIVE.NoDaysOfFollowUp + 1, SURVIVE.NoSubjectCategories, SURVIVE, HOUSE);

	if (SURVIVE.CalculateSeparatePassivePhaseCurves)
		AddToSurvivalTailsAndMean(SURVIVE.SC_PP, SURVIVE.NoPassiveDaysFollowUp + 1, SURVIVE.NoSubjectCategories, SURVIVE, HOUSE);
}
void Update_Hazards_Etc_For_SurvivalCurves	(const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// 
	////// Update the coefficients hazard values etc. that are not used in likelihood calculation but that are used in survival curve calculation (i.e. the ones after maximum calendar time for DATA, but not for prediction, e.g. a patient who enrolls late and still requires 3 years of survival curve prediction). 
	////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// 

	bool HazConditionToSkip = 0;
	for (int country = 0; country < HOUSE.TotalCountries; country++)
	{
		HazConditionToSkip = !(std::any_of(HOUSE.WhichCountries.begin(), HOUSE.WhichCountries.end(), [&](int i) {return i == country; }));
		if (HazConditionToSkip) continue;

		//// Amend SplineCoeffs
		for (int poly = 0; poly < HOUSE.PolynomialsPerCountry; poly++)
			CoefficientsFromKnots(PARAMS.xKnots[country], PARAMS.yKnots[country], poly, PARAMS.SplineCoeffs[country][poly]); 

		//// Change Baseline Hazard Values. 
#pragma omp parallel for schedule(static,1) 
		for (int th = DATA.NumCalendarDaysFollowUp[country]; th < HOUSE.max_threads; th++)
			for (int CalendarDay = th; CalendarDay <= (DATA.MaxStartFollowUp_CalendarTime + DATA.NoDaysOfFollowUp + 1); CalendarDay += HOUSE.max_threads)
				PARAMS.BaselineHazardValues[country][CalendarDay] = BaselineHazard((CalendarDay * HOUSE.TimeInterval), country, PARAMS, HOUSE);

		//// Change Integrated Hazard Values
		DType IntBaseHazCumulative = PARAMS.IntBaseHazLookUp[country][DATA.NumCalendarDaysFollowUp[country]];
		for (int CalendarDay = DATA.NumCalendarDaysFollowUp[country] + 1; CalendarDay <= (DATA.MaxStartFollowUp_CalendarTime + DATA.NoDaysOfFollowUp); CalendarDay++)
		{
			// add to cumulative / approximate integrated hazard. 
			PARAMS.IntBaseHazLookUp[country][CalendarDay] = IntBaseHazCumulative;
			IntBaseHazCumulative += PARAMS.BaselineHazardValues[country][CalendarDay] * HOUSE.TimeInterval;	// tempting to multiply by TimeInterval once at end as common factor, but as TimeInterval interval required for each day then don't do this.
		}
	}
}
void ResetSurvivalTable						(DType **Table, int NumStrata, int NumDays) //// use this function instead of "Populate" functions because of first column. 
{
	for (int stratum = 0; stratum < NumStrata; stratum++)
	{
		Table[stratum][0] = 1;
		for (int daypostdose = 1; daypostdose < NumDays; daypostdose++)		Table[stratum][daypostdose] = 0;
	}
}
void ResetStuff								(Survival_Struct &SURVIVE, const DATA_struct &DATA)
{
	//// NOTE:  used for Survival tables only: not for MeanHazVals or HRs (hazard ratios). 

	//// reset Strata_Sizes
	Populate_2D_Array(SURVIVE.SC_WT.Strata_Sizes, (DType)0, SURVIVE.NoSubjectCategories, DATA.NoDaysOfFollowUp + 1);

	//////// clear non-thread-specific survival table. 
	for (int DiseaseSeverity = 0; DiseaseSeverity < SURVIVE.HowManyDiseaseSeverities; DiseaseSeverity++)
		ResetSurvivalTable(SURVIVE.SC_WT.PostSample_SurvivalTables[DiseaseSeverity], SURVIVE.NoSubjectCategories, SURVIVE.NoDaysOfFollowUp + 1);

	//// //// //// //// //// //// //// //// Passive Survival Curve (basically same as above). 
	if (SURVIVE.CalculateSeparatePassivePhaseCurves)
	{
		//// reset Strata_Sizes
		Populate_2D_Array(SURVIVE.SC_PP.Strata_Sizes, (DType)0, SURVIVE.NoSubjectCategories, SURVIVE.NoPassiveDaysFollowUp + 1);

		for (int DiseaseSeverity = 0; DiseaseSeverity < SURVIVE.HowManyDiseaseSeverities; DiseaseSeverity++)
			ResetSurvivalTable(SURVIVE.SC_PP.PostSample_SurvivalTables[DiseaseSeverity], SURVIVE.NoSubjectCategories, SURVIVE.NoPassiveDaysFollowUp + 1);
	}
}
void ResetThreadedStuff						(int th, Survival_Struct &SURVIVE, const DATA_struct &DATA)
{
	//// NOTE: used this for Survival tables and attack rates only: not for MeanHazVals or HRs (hazard ratios). 
	
	//// reset Strata_Sizes
	Populate_2D_Array(SURVIVE.SC_WT.Threaded_Strata_Sizes[th], (DType)0, SURVIVE.NoSubjectCategories, DATA.NoDaysOfFollowUp + 1);

	//// initialize survival table to zero (except for first column because you won't amend this in rest of loop).
	for (int DiseaseSeverity = 0; DiseaseSeverity < SURVIVE.HowManyDiseaseSeverities; DiseaseSeverity++)
		ResetSurvivalTable(SURVIVE.SC_WT.Threaded_PostSample_SurvivalTables[th][DiseaseSeverity], SURVIVE.NoSubjectCategories, SURVIVE.NoDaysOfFollowUp + 1);

	//// reset Sum_FU_Durations (denominators of attack rate). 
	Populate_2D_Array(SURVIVE.Sum_FU_Durations[th], (DType)0, SURVIVE.HowManyTrialPhases, SURVIVE.NoSubjectCategories); 

	//// reset threaded Attack rates
	Populate_3D_Array(SURVIVE.Threaded_AttackRates[th], (DType)0, SURVIVE.HowManyTrialPhases, SURVIVE.HowManyDiseaseSeverities, SURVIVE.NoSubjectCategories);

	//// //// //// //// //// //// //// //// Passive Survival Curve (basically same as above). 

	if (SURVIVE.CalculateSeparatePassivePhaseCurves)
	{
		//// reset Strata_Sizes
		Populate_2D_Array(SURVIVE.SC_PP.Threaded_Strata_Sizes[th], (DType)0, SURVIVE.NoSubjectCategories, SURVIVE.NoPassiveDaysFollowUp + 1);

		for (int DiseaseSeverity = 0; DiseaseSeverity < SURVIVE.HowManyDiseaseSeverities; DiseaseSeverity++)
			ResetSurvivalTable(SURVIVE.SC_PP.Threaded_PostSample_SurvivalTables[th][DiseaseSeverity], SURVIVE.NoSubjectCategories, SURVIVE.NoPassiveDaysFollowUp + 1);
	}
}
int  DetermineNumAgeGroups					(int country)
{
	int NoAgeGroupsToAddTo; 
	if (country <= 4)	/*CYD-14*/	NoAgeGroupsToAddTo = 4;		else	/*CYD-15*/	NoAgeGroupsToAddTo = 3;
	return NoAgeGroupsToAddTo; 
}
void DetermineStrata						(int patient, Survival_Struct &SURVIVE, const DATA_struct &DATA, int *ThisPatientsStrata)
{
	int country = DATA.ci_s[patient];

	int *CountryGrouping	= new int[3]();	// always add to three "country groups": 0 = patient's country			  ;	1 = patient's trial;   2 = combined CYD14 and 15 trial.  
	int *HazardGrouping		= new int[4](); // always add to four "hazard groups"	: 0 = patient's arm and serostatus;	1 = patient's arm;	   2 = patient's serostatus,	and 3 = either arm, either status. 
	int *AgeGroupGrouping	= new int[4]();	// add a constant amount of memory (at least as big as ever needed), then only access the bits needed for a particular patient. 
	int *HazRatioGrouping	= new int[2](); // always add to two "serostatuses": 0 = patient's serostatus			  ;	1 = either status

	HazRatioGrouping[0] = DATA.Ii_s[patient];	// patient's serostatus 
	HazRatioGrouping[1] = EitherSeroStatus;		// either arm, either serostatus 


	/*		
		AGE-GROUP NUMBERS / DEFINITIONS

		0		= any age group
		1,2,3	= CYD14 trial 2 to 5 years, 6 to 11 years, 12 to 14 years
		4,5		= CYD15 trial 9 to 11 years, 12 to 16 years.
		6,7		= CYD14 trial under 9 years, over 9 years
		8		= CYD14 and CYD15 trial over 9 years (no need for CYD15 under 9 years as all CYD15 participants are under 9 years, so this is just the same as CYD15 any age group).
	*/

	/*	
		HAZARD-GROUP NUMBERS / DEFINITIONS
	
		0		= SeroPos Control					(1 in R notation)
		1		= SeroNeg Control					(2 in R notation)
		2		= SeroPos Vaccine					(3 in R notation)
		3		= SeroNeg Vaccine					(4 in R notation)
		4		= Control							(5 in R notation)
		5		= Vaccine							(6 in R notation)
		6		= SeroPos 							(7 in R notation)
		7		= SeroNeg 							(8 in R notation)
		8		= Either Arm, Either Serostatus		(9 in R notation)
	*/

	if (country <= 4)		// CYD-14
	{
		CountryGrouping [1] = 10;				// CYD-14

		AgeGroupGrouping[1] = DATA.AgeGroup1[patient];
		AgeGroupGrouping[2] = DATA.AgeGroup2[patient];
		AgeGroupGrouping[3] = 8; 
	}
	else					// CYD-15
	{
		CountryGrouping [1] = 11;				// CYD-15

		AgeGroupGrouping[1] = DATA.AgeGroup1[patient];
		AgeGroupGrouping[2] = 8;
	}
	CountryGrouping[0] = country; CountryGrouping[2] = 12;	// CYD 14 and 15 pooled


	HazardGrouping[3] = 8; //// either arm, either serostatus 
	if (DATA.Vi_s[patient] == 0) 
	{
		HazardGrouping[1] = 4; // control arm.
		if (DATA.Ii_s[patient] == 1)	{	HazardGrouping[0] = 0;	/* SeroPos at baseline	- control arm. */	HazardGrouping[2] = 6; }
		else							{	HazardGrouping[0] = 1;	/* SeroNeg at baseline  - control arm. */	HazardGrouping[2] = 7; }
	}
	else 
	{
		HazardGrouping[1] = 5; // vaccine arm.
		if (DATA.Ii_s[patient] == 1)	{	HazardGrouping[0] = 2;	/* SeroPos at baseline	- vaccine arm. */	HazardGrouping[2] = 6; }
		else							{	HazardGrouping[0] = 3;	/* SeroNeg at baseline  - vaccine arm. */	HazardGrouping[2] = 7; }
	}

	int NoAgeGroupsToAddTo = DetermineNumAgeGroups(country); 

	int strata_counter = 0;
	for (int ageDummy = 0; ageDummy < NoAgeGroupsToAddTo; ageDummy++)
		for (int countryDummy = 0; countryDummy < 3; countryDummy++)
			for (int hazardDummy = 0; hazardDummy < 4; hazardDummy++)
			{
				ThisPatientsStrata[strata_counter] = (SURVIVE.NoSubjectCategoriesPerAgeGroup * AgeGroupGrouping[ageDummy]) + (SURVIVE.TotalGroupsPerCountry * CountryGrouping[countryDummy]) + HazardGrouping[hazardDummy];
				strata_counter++;
			}

	delete[] AgeGroupGrouping;
	delete[] CountryGrouping;
	delete[] HazardGrouping;
}

void AddToStrataSizes_fromThreads			(SCurves &SetOfSCurves, int NumDays, const DATA_struct &DATA, const Housekeeping_Struct &HOUSE, Survival_Struct &SURVIVE)
{
	///// Add to strata sizes from different threads
	for (int category = 0; category < SURVIVE.NoSubjectCategories; category++)
		for (int th = 0; th < HOUSE.max_threads; th++)
			for (int daypostdose = 0; daypostdose < NumDays; daypostdose++)
				SetOfSCurves.Strata_Sizes[category][daypostdose] += SetOfSCurves.Threaded_Strata_Sizes[th][category][daypostdose];
}
void AddToStrataSizes_fromThreads			(const DATA_struct &DATA, const Housekeeping_Struct &HOUSE, Survival_Struct &SURVIVE)
{
	//// Wrapper of overloaded function above.
	///// Add to strata sizes from different threads
	AddToStrataSizes_fromThreads(SURVIVE.SC_WT, DATA.NoDaysOfFollowUp + 1, DATA, HOUSE, SURVIVE);

	if (SURVIVE.CalculateSeparatePassivePhaseCurves)///// Add to Passive Survival tables for this posterior sample from different threads
		AddToStrataSizes_fromThreads(SURVIVE.SC_PP, SURVIVE.NoPassiveDaysFollowUp + 1, DATA, HOUSE, SURVIVE);
}
void AddToSurvivalTables_fromThreads		(SCurves &SetOfSCurves, int NoDays, const Housekeeping_Struct &HOUSE, Survival_Struct &SURVIVE, int StartDay)
{
	//// Start Day - for survival tables you start from day 1, not zero. Survival probability day 0 is 1

	///// Add to Survival tables for this posterior sample from different threads
	for (int DiseaseSeverity = 0; DiseaseSeverity < SURVIVE.HowManyDiseaseSeverities; DiseaseSeverity++)
		for (int category = 0; category < SURVIVE.NoSubjectCategories; category++)
			for (int daypostdose = StartDay; daypostdose < NoDays; daypostdose++)
				for (int th = 0; th < HOUSE.max_threads; th++)
					SetOfSCurves.PostSample_SurvivalTables[DiseaseSeverity][category][daypostdose] +=
					SetOfSCurves.Threaded_PostSample_SurvivalTables[th][DiseaseSeverity][category][daypostdose];
}
void AddToSurvivalTables_fromThreads		(const Housekeeping_Struct &HOUSE, Survival_Struct &SURVIVE)
{
	// Wrapper of overloaded function AddToSurvivalTables_fromThreads above.
	// Add to Survival tables for this posterior sample from different threads
	AddToSurvivalTables_fromThreads(SURVIVE.SC_WT, SURVIVE.NoDaysOfFollowUp + 1, HOUSE, SURVIVE, 1);

	if (SURVIVE.CalculateSeparatePassivePhaseCurves)///// Add to Passive Survival tables for this posterior sample from different threads
		AddToSurvivalTables_fromThreads(SURVIVE.SC_PP, SURVIVE.NoPassiveDaysFollowUp + 1, HOUSE, SURVIVE, 1);
}
void AddToAttackRates_fromThreads			(int AR_PostSampleNo, const Housekeeping_Struct &HOUSE, Survival_Struct &SURVIVE, DType **** AttackRateContainer)
{
	//// for mean and modal posterior attack rate estimates, dd to MeanModeAttackRates. Otherwise (default), add to SURVIVE.MetaAttackRates. 
	for (int TrialPhase = 0; TrialPhase < SURVIVE.HowManyTrialPhases; TrialPhase++) //// note SURVIVE.HowManyTrialPhases not HOUSE.HowManyTrialPhases
		for (int DiseaseSeverity = 0; DiseaseSeverity < SURVIVE.HowManyDiseaseSeverities; DiseaseSeverity++)
			for (int category = 0; category < SURVIVE.NoSubjectCategories; category++)
				for (int th = 0; th < HOUSE.max_threads; th++)
					AttackRateContainer[TrialPhase][DiseaseSeverity][category][AR_PostSampleNo] += SURVIVE.Threaded_AttackRates[th][TrialPhase][DiseaseSeverity][category];
}


void AddToSum_FU_Durations_fromThreads_AndCalculateAttackRates(int AR_PostSampleNo, const Housekeeping_Struct &HOUSE, Survival_Struct &SURVIVE, DType **** AttackRateContainer)
{
	int Sum_FU_duration_phase_stratum = 0; //// i.e. the denominator of attack rates for a given Trial Phase and stratum / category
	for (int TrialPhase = 0; TrialPhase < SURVIVE.HowManyTrialPhases; TrialPhase++)
		for (int stratum = 0; stratum < SURVIVE.NoSubjectCategories; stratum++)
		{
			//// (reset and) sum the FU_duration for this trial phase and stratum. 
			Sum_FU_duration_phase_stratum = 0;
			for (int th = 0; th < HOUSE.max_threads; th++)
				Sum_FU_duration_phase_stratum += SURVIVE.Sum_FU_Durations[th][TrialPhase][stratum];

			//// for each DiseaseSeverity, divide MetaAttackRates (which is currently only the numerator) by Sum_FU_duration_phase_stratum (i.e. the denominator). 
			for (int DiseaseSeverity = 0; DiseaseSeverity < SURVIVE.HowManyDiseaseSeverities; DiseaseSeverity++)
				if (Sum_FU_duration_phase_stratum != 0)
					AttackRateContainer[TrialPhase][DiseaseSeverity][stratum][AR_PostSampleNo] /= Sum_FU_duration_phase_stratum;

			//// Convert to percentages
			for (int DiseaseSeverity = 0; DiseaseSeverity < SURVIVE.HowManyDiseaseSeverities; DiseaseSeverity++)
				AttackRateContainer[TrialPhase][DiseaseSeverity][stratum][AR_PostSampleNo] *= 100;
		}
}
int ChooseOldNotationStratumNumber			(int AgeGroup, int Country, int TrialArm, int BaselineSeroStatus, Survival_Struct &SURVIVE)
{
	int HazardNumberOldNotation = SURVIVE.Strata_SCurves[TrialArm][BaselineSeroStatus];
	int CountryAgeGroupIndex	= (AgeGroup * SURVIVE.TotalCountries * SURVIVE.TotalGroupsPerCountry) + (Country * SURVIVE.TotalGroupsPerCountry);

	int Stratum = CountryAgeGroupIndex + HazardNumberOldNotation; 

	return Stratum;
}
void CalculateRatioOfMeanHazards			(const Housekeeping_Struct &HOUSE, Survival_Struct &SURVIVE)
{
	int VacStratum, ConStratum, HR_Stratum; //// if you parallelize, remember to put these within threading loop. 

	for (int AgeGroup = 0; AgeGroup < SURVIVE.HRs_NumAgeGroups; AgeGroup++)
		for (int Country = 0; Country < SURVIVE.TotalCountries; Country++)
			for (int BaselineSeroStatus = 0; BaselineSeroStatus <= HOUSE.HowManySeroStatuses; BaselineSeroStatus++) /// note the <= sign. Want EitherSeroStatus as well. 
			{
				//// Choose strata indices for vaccine and control group (given AgeGroup, Country, and BaselineSeroStatus) for mean hazard ratios. 
				VacStratum = ChooseOldNotationStratumNumber(AgeGroup, Country, VaccineGroup, BaselineSeroStatus, SURVIVE);
				ConStratum = ChooseOldNotationStratumNumber(AgeGroup, Country, ControlGroup, BaselineSeroStatus, SURVIVE);

				//// Choose stratum index for Hazard Ratios
				HR_Stratum = (AgeGroup * SURVIVE.HRs_NumCountries * SURVIVE.HRs_NumHRsPerAgeGroupAndCountry) + (Country * SURVIVE.HRs_NumHRsPerAgeGroupAndCountry) + BaselineSeroStatus;

				for (int HR_serotype = 0; HR_serotype < SURVIVE.HRs_NumSTypes; HR_serotype++)
					for (int DiseaseSeverity = 0; DiseaseSeverity < SURVIVE.HowManyDiseaseSeverities; DiseaseSeverity++)
						for (int DayPostDoseIndex = 0; DayPostDoseIndex < SURVIVE.HRs_DaysPostDose.size(); DayPostDoseIndex++)
							if (SURVIVE.MeanHazVals_BySerotype	[HR_serotype].PostSample_SurvivalTables[DiseaseSeverity][ConStratum][DayPostDoseIndex] != 0)
								SURVIVE.HRs_BySerotype			[HR_serotype].PostSample_SurvivalTables[DiseaseSeverity][HR_Stratum][DayPostDoseIndex] =
								SURVIVE.MeanHazVals_BySerotype	[HR_serotype].PostSample_SurvivalTables[DiseaseSeverity][VacStratum][DayPostDoseIndex] /
								SURVIVE.MeanHazVals_BySerotype	[HR_serotype].PostSample_SurvivalTables[DiseaseSeverity][ConStratum][DayPostDoseIndex] ;
			}
}
void Process_MeanHaz_And_HRs				(const Housekeeping_Struct &HOUSE, Survival_Struct &SURVIVE)
{
	//// Call this function after daypostdose and thread loops in GenerateSurvivalCurves function. 
	//// Mean Hazard Values (for Hazard Ratios)
	for (int HR_serotype = 0; HR_serotype < SURVIVE.HRs_NumSTypes; HR_serotype++)
		AddToSurvivalTables_fromThreads(SURVIVE.MeanHazVals_BySerotype[HR_serotype], SURVIVE.HRs_NumDaysPostDoseToCalculate, HOUSE, SURVIVE, 0);

	//// Divide the summed hazard values (done in GenerateSurvivalCurves function) by strata sizes (for all DiseaseSeverities you're doing). 
	for (int HR_serotype = 0; HR_serotype < SURVIVE.HRs_NumSTypes; HR_serotype++)
		DivideSummedProbabilies(SURVIVE.MeanHazVals_BySerotype[HR_serotype], SURVIVE.HRs_NumDaysPostDoseToCalculate, SURVIVE.NoSubjectCategories, SURVIVE, HOUSE); //// note here you do use SURVIVE.NoSubjectCategories. 

	//// now have mean hazards or every strata. Calculate Ratios: this changes values in HRs (PostSample_SurvivalTables) from MeanHazVals. 
	CalculateRatioOfMeanHazards(HOUSE, SURVIVE); 

	//// Now have ratio of mean hazards HRs. Add to posterior samples (mean and tails). 
	for (int HR_serotype = 0; HR_serotype < SURVIVE.HRs_NumSTypes; HR_serotype++)
		AdjustTailsAndAddToFinalMeanPosterior(SURVIVE.HRs_BySerotype[HR_serotype], SURVIVE.HRs_NumDaysPostDoseToCalculate, SURVIVE.HRs_NumStrata, SURVIVE, HOUSE);

	//// reset (i.e. threaded and non-threaded MeanHazVals). 
	for (int HR_serotype = 0; HR_serotype < SURVIVE.HRs_NumSTypes; HR_serotype++)
		for (int DiseaseSeverity = 0; DiseaseSeverity < SURVIVE.HowManyDiseaseSeverities; DiseaseSeverity++)
		{
			//// MeanHazVals_BySerotype[HR_serotype].PostSample_SurvivalTables (non-threaded). Populate to 0. 
			Populate_2D_Array(SURVIVE.MeanHazVals_BySerotype[HR_serotype].PostSample_SurvivalTables[DiseaseSeverity], (DType)0, SURVIVE.NoSubjectCategories, SURVIVE.HRs_NumDaysPostDoseToCalculate);

			//// MeanHazVals_BySerotype[HR_serotype].Threaded_PostSample_SurvivalTables (threaded). Populate to 0. 
			for (int thread = 0; thread < HOUSE.max_threads; thread++)
				Populate_2D_Array(SURVIVE.MeanHazVals_BySerotype[HR_serotype].Threaded_PostSample_SurvivalTables[thread][DiseaseSeverity], (DType)0, SURVIVE.NoSubjectCategories, SURVIVE.HRs_NumDaysPostDoseToCalculate);
		}
}

void GenerateSurvivalCurves(const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE,
	Survival_Struct &SURVIVE, int iter, const Chains_Struct &CHAINS, Augment_Struct &AUG, bool AllPatients)
{
	//// overload of GenerateSurvivalCurves function that essentially applies defaults. 

	string ImSubString							= AllPatients ? "" : "_ImSub"; 
	string ParamSet_SCurve_RootFilename			= "ParamSet_" + std::to_string((iter - CHAINS.BurnIn) / SURVIVE.AddToSurvivalCurvesEveryHowManyIterations) + "_SurvivalTable"			+ ImSubString;
	string ParamSet_PassiveSCurve_RootFilename	= "ParamSet_" + std::to_string((iter - CHAINS.BurnIn) / SURVIVE.AddToSurvivalCurvesEveryHowManyIterations) + "_PassiveSurvivalTable"	+ ImSubString;

	///// Value of iter in this function overload is value of Iteration_OR_AR_Post_Stat_Index in general function (i.e. will refer to chain update number as before). 
	///// When doing Mean and modal post, don't need chain update number, but  do need an index. Therefore use same argument to mean different things as required. 
	GenerateSurvivalCurves(DATA, PARAMS, HOUSE, SURVIVE, iter, CHAINS, AUG, AllPatients, CHAINS.OutputIndividualSurvivalTables, CHAINS.OutputIndividualSurvivalTables, SURVIVE.MetaAttackRates,
		ParamSet_SCurve_RootFilename, ParamSet_PassiveSCurve_RootFilename, true, false);
}

void GenerateSurvivalCurves					(const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, 
	Survival_Struct &SURVIVE, int Iteration_OR_AR_Post_Stat_Index, const Chains_Struct &CHAINS, Augment_Struct &AUG, bool AllPatients,
	//// these arguments were put in for Mean and Modal Posterior SCurves and attack rates. In  main code, these aren't set - that is, overloaded function which gives the defaults is called. 
	bool OutputIndividualSurvivalTables, bool OutputIndividual_Passive_SurvivalTables, DType **** AttackRateContainer,
	string ParamSet_SCurve_RootFilename, string ParamSet_PassiveSCurve_RootFilename, bool Use_Default_AR_Post_Sample_No, bool CleanFirstDay_IndividualCurves)
{
	/*
		Function calculates survival probabilities and various attack rates
		Function is threaded and uses threaded and un-threaded quantities. 
		Structure is: 

			i)		Clear everything from previous iterations
			ii)		Begin threading
			ii)		Calculate threaded survival curves and attack rates. Loop over patients, loop over dayspostdose, 
			iii)	End threading: sum / gather everything
			iv)		Add to posterior samples (mean and maybe change lower/upper CrIs)
	*/

	clock_t timer = clock();

	int AR_PostSampleNo; /// attack rate posterior sample number/index for tables. By default, add to SURVIVE.MetaAttackRates. However for mean and max like post samples, add to specific indices. 
	if (Use_Default_AR_Post_Sample_No)
		AR_PostSampleNo = (Iteration_OR_AR_Post_Stat_Index - CHAINS.BurnIn) / (SURVIVE.AddToSurvivalCurvesEveryHowManyIterations);
	else
		AR_PostSampleNo = Iteration_OR_AR_Post_Stat_Index;

	//// reset Strata_Sizes and clear non-thread-specific survival table. 
	ResetStuff(SURVIVE, DATA); 

	////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// ////// 
	////// Update the coefficients hazard values etc. that are not used in likelihood calculation but that are used in survival curve calculation (i.e. the ones after maximum calendar time for DATA, but not for prediction, e.g. a patient who enrolls really late and still requires 3 years of survival curve prediction). 

	Update_Hazards_Etc_For_SurvivalCurves(DATA, PARAMS, HOUSE); 

	// Parallelize over patients - not days. This way you can increment integral from t-1 to t 
	// but since survival probabilities must be summed and averaged over each patient, multiple patients contribute to a given survival probability, hence multiple threads contribute to given survival probability. 

#pragma omp parallel for schedule(static,1)
	for (int th = 0; th < HOUSE.max_threads; th++)
	{
		//// reset Strata_Sizes, Sum_FU_Durations, threaded Attack rates and initialize survival table to zero (except for first column).
		ResetThreadedStuff(th, SURVIVE, DATA);

		//// //// //// //// //// //// //// ////	//// //// //// //// //// //// //// //// 
		//// Allocate memory and quantities for loop below. 

		int *ThisPatientsStrata				= new int[4*4*3](); //// However many age groups (max 4) * 3 "countries" * 4 "hazard groups". 
		int *ThisPatients_HazRatioStrata	= new int[4*3*2](); //// However many age groups (max 4) * 3 "countries" * 2 "serostatuses". 

		
		int country, ThisDay, PreviousDay, NumStrata_thisPatient;
		DType SurvivalProb, SurvivalProb_Severe, SurvivalProb_Either, AM_IntBaseHaz_ThisDay, IntVacHaz_AM_RunningTotal = 0, PartitionHazardMultiplier = 1; ///// Don't want a million if statements, so set this equal to one if not doing partition method, but change below (by country) otherwise
		bool AllCool = 1, ConditionToSkip, ConditionToCalcHazRatio; 
		//// Hazard ratio quantities. 
		DType BaseHazValue, VacHazValue; 
		DType HazValue_stratum_Either	, BaseHazMult_Either	, VacHazMult_Either	; //// you use these when HOUSE.MildSevere == TREATED_EQUALLY
		DType HazValue_stratum_Mild		, BaseHazMult_Mild		, VacHazMult_Mild	;
		DType HazValue_stratum_Severe	, BaseHazMult_Severe	, VacHazMult_Severe	;
		int WhichDayIndex; 

		int NoActiveDays_ThisPatient = NULL, NoActivePlusPassiveDays_ThisPatient = NULL, PassiveDayPostDose;

		int patient = NULL; 
		int MaxPatientIndex = AllPatients ? NPat : DATA.NonAugmentedIndices.size(); 

		for (int patientindex = th; patientindex < MaxPatientIndex; patientindex += HOUSE.max_threads)
		{
			if (AllPatients) patient = DATA.PatientIndexNumbers[patientindex]; else patient = DATA.NonAugmentedIndices[patientindex]; 

			country			= DATA.ci_s[patient];
			ConditionToSkip = !(std::any_of(HOUSE.WhichCountries.begin(), HOUSE.WhichCountries.end(), [&](int i){return i == country; }));
			if (ConditionToSkip) continue; 

			if (HOUSE.SeroSpecific_K_values || HOUSE.SeroSpecificEfficacies)	if (HOUSE.BaselinePartition) PartitionHazardMultiplier = PARAMS.SumRhos[country]; 

			///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
			///// Choose K values - need one K value for ActiveMild and another for PassiveSevere
			///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
			
			Choose_iHazMults				(th, patient, DATA.Ii_s[patient], DATA, PARAMS, HOUSE, AUG);	//// used for survival curves and "combined serotypes" of hazard ratios
			if (HOUSE.SeroSpecific_K_values || HOUSE.SeroSpecificEfficacies)
				Choose_HazMults_SingleSero	(th, patient, DATA, PARAMS, HOUSE, AUG);	//// used for single serotypes of hazard ratios

			DType K_rho_AM		= AUG.iBaseHaz_Mult	[th][ActiveMild][DATA.Ii_s[patient]];
			DType K_Eff_rho_AM	= AUG.iVacHaz_Mult	[th][ActiveMild][DATA.Ii_s[patient]];
			DType K_rho_PS		= NULL;
			DType K_Eff_rho_PS	= NULL; 

			if (HOUSE.HowManyCaseCategories == 2)
			{
				K_rho_PS		= AUG.iBaseHaz_Mult	[th][PassiveSevere	][DATA.Ii_s[patient]];
				K_Eff_rho_PS	= AUG.iVacHaz_Mult	[th][PassiveSevere	][DATA.Ii_s[patient]];
			}

			//// Which Strata does this patient belong to. 
			NumStrata_thisPatient = DetermineNumAgeGroups(country) * 3 * 4; //// 3 "countries" including trial and CYD-14/15 together; 4 hazard groups. See DetermineStrata function for more details. 
			DetermineStrata(patient, SURVIVE, DATA, ThisPatientsStrata);

			/// add to apprppriate group population size. 
			for (int stratum = 0; stratum < NumStrata_thisPatient; stratum++)
				SURVIVE.SC_WT.Threaded_Strata_Sizes[th][ThisPatientsStrata[stratum]][0]++;

			if (HOUSE.SFU && HOUSE.ActiveOrPassivePhase == DO_ACTIVE_AND_PASSIVE) //// If HOUSE.SFU, trial phase durations differ by patient. However, if doing the active phase only then  do the longest active phase duration, i.e. SURVIVE.NoActiveDaysFollowUp. If !HOUSE.SFU, then patients' trial phase durations do not vary between them. 
			{
				NoActiveDays_ThisPatient			= DATA.FU_Duration_days[ActivePhase	][patient];
				NoActivePlusPassiveDays_ThisPatient = DATA.FU_Duration_days[WholeTrial	][patient];
			}
			else
			{
				NoActiveDays_ThisPatient			= SURVIVE.NoActiveDaysFollowUp;
				NoActivePlusPassiveDays_ThisPatient = SURVIVE.NoDaysOfFollowUp;
			}

			///// Add to attack rate denominators for attack rate
			for (int TrialPhase = 0; TrialPhase < SURVIVE.HowManyTrialPhases; TrialPhase++)
				for (int stratum = 0; stratum < NumStrata_thisPatient; stratum++)
					SURVIVE.Sum_FU_Durations[th][TrialPhase][ThisPatientsStrata[stratum]] += DATA.FU_Duration_years[TrialPhase][patient];

			if (HOUSE.MildAndSevere == TREATED_SEPARATELY)
			{
				DType IntBaseHaz_ThisDay = 0, IntVacHaz_RunningTotal = 0, Haz_3_IntRunningTotal = 0,  Haz_4_IntRunningTotal = 0;
				DType SurviveActive, SurviveActive_Either, SurviveActive_Severe;	/// For passive phase attack rate calculation 
				int NoDaysToLoopOver;

				if (HOUSE.LTFU_SurvivalCurves) //// if accounting for loss to follow up, SURVIVE.NoDaysOfFollowUp won't account for individual variation in follow up times. 
					if (HOUSE.ActiveOrPassivePhase == DO_ACTIVE_AND_PASSIVE)	NoDaysToLoopOver = NoActivePlusPassiveDays_ThisPatient;
					else														NoDaysToLoopOver = NoActiveDays_ThisPatient		;
				else  NoDaysToLoopOver = SURVIVE.NoDaysOfFollowUp; //// but if not accoutning for LTFU, then SURVIVE.NoDaysOfFollowUp accounts for active or passive. 

				for (int daypostdose = 1; daypostdose <= NoDaysToLoopOver; daypostdose++)
				{
					ConditionToCalcHazRatio = (std::any_of(SURVIVE.HRs_DaysPostDose.begin(), SURVIVE.HRs_DaysPostDose.end(), [&](int i) {return i == (daypostdose - 1); })); //// if day is equal to one of HRs_DaysPostDose that we want to calculate. 
					if (CHAINS.HazRatiosOnly & !ConditionToCalcHazRatio) continue;

					/// add to apprppriate group population size. 
					for (int stratum = 0; stratum < NumStrata_thisPatient; stratum++)
						SURVIVE.SC_WT.Threaded_Strata_Sizes[th][ThisPatientsStrata[stratum]][daypostdose]++;

					ThisDay				= DATA.FollowUp[START][ActiveMild][patient] + daypostdose;
					PreviousDay			= DATA.FollowUp[START][ActiveMild][patient] + daypostdose - 1;
					IntBaseHaz_ThisDay	= IntBaseHaz(DATA.FollowUp[START][ActiveMild][patient], ThisDay, country, PARAMS) * PartitionHazardMultiplier;
				
					if (DATA.Vi_s[patient] == 0) 			// country hazard 1 or 2 (0 or 1 c++ indexing)
					{
						SurvivalProb		= exp(-IntBaseHaz_ThisDay * K_rho_AM); 
						SurvivalProb_Severe = exp(-IntBaseHaz_ThisDay *	K_rho_PS);
					}
					else 									// country hazard 3 or 4 (2 or 3 c++ indexing)
					{
						if (!HOUSE.AllDosesRequired_AG_BS[DATA.AgeGroup1[patient]]	[DATA.Ii_s[patient]]	|| PreviousDay >= DATA.ThirdDose[patient]) //// i.e. only add to the vaccine hazard if past the third dose, or if not all three doses are required for this patient's serostatus (latter is default)
							IntVacHaz_RunningTotal += IntVacHaz_SingleDay(PreviousDay, DATA.FollowUp[START][ActiveMild][patient], DATA.SecondDose[patient], DATA.ThirdDose[patient], DATA.ai_s[patient], DATA.Ii_s[patient], country, PARAMS, HOUSE);

						SurvivalProb		= exp(	(K_Eff_rho_AM * IntVacHaz_RunningTotal - (K_rho_AM * IntBaseHaz_ThisDay))	);
						SurvivalProb_Severe = exp(	(K_Eff_rho_PS * IntVacHaz_RunningTotal - (K_rho_PS * IntBaseHaz_ThisDay))	);
					}

					////  Calcualate Prob(Mild Or Severe)
					SurvivalProb_Either = SurvivalProb * SurvivalProb_Severe; //// i.e. the probability survived Both mild and severe disease

					/// add to apprpriate groups
					for (int stratum = 0; stratum < NumStrata_thisPatient; stratum++)
					{
						SURVIVE.SC_WT.Threaded_PostSample_SurvivalTables[th][EITHER			][ThisPatientsStrata[stratum]][daypostdose] += SurvivalProb_Either;
						SURVIVE.SC_WT.Threaded_PostSample_SurvivalTables[th][MILD_NON_HOSP	][ThisPatientsStrata[stratum]][daypostdose] += SurvivalProb;
						SURVIVE.SC_WT.Threaded_PostSample_SurvivalTables[th][SEVERE_HOSP	][ThisPatientsStrata[stratum]][daypostdose] += SurvivalProb_Severe;
					}

					if (daypostdose == NoActiveDays_ThisPatient) ///// add to active phase attack rate for appropriate groups. Needs to be thread-specific
					{
						for (int stratum = 0; stratum < NumStrata_thisPatient; stratum++)
						{
							///// Active Phase Only. 
							SURVIVE.Threaded_AttackRates[th][ActivePhase][EITHER		][ThisPatientsStrata[stratum]] += (DType)1 - SurvivalProb_Either;
							SURVIVE.Threaded_AttackRates[th][ActivePhase][MILD_NON_HOSP	][ThisPatientsStrata[stratum]] += (DType)1 - SurvivalProb;
							SURVIVE.Threaded_AttackRates[th][ActivePhase][SEVERE_HOSP	][ThisPatientsStrata[stratum]] += (DType)1 - SurvivalProb_Severe;
						}
						
						SurviveActive			= SurvivalProb;
						SurviveActive_Severe	= SurvivalProb_Severe;
						SurviveActive_Either	= SurvivalProb_Either;
					}
					else if (daypostdose == NoActivePlusPassiveDays_ThisPatient)
						for (int stratum = 0; stratum < NumStrata_thisPatient; stratum++)
						{
							///// Passive Phase Only. 
							SURVIVE.Threaded_AttackRates[th][PassivePhase][EITHER		][ThisPatientsStrata[stratum]] += SurviveActive_Either	- SurvivalProb_Either	;
							SURVIVE.Threaded_AttackRates[th][PassivePhase][MILD_NON_HOSP][ThisPatientsStrata[stratum]] += SurviveActive			- SurvivalProb			;
							SURVIVE.Threaded_AttackRates[th][PassivePhase][SEVERE_HOSP	][ThisPatientsStrata[stratum]] += SurviveActive_Severe	- SurvivalProb_Severe	;

							///// Full Trial  
							SURVIVE.Threaded_AttackRates[th][WholeTrial][EITHER			][ThisPatientsStrata[stratum]] += (DType)1 - SurvivalProb_Either;
							SURVIVE.Threaded_AttackRates[th][WholeTrial][MILD_NON_HOSP	][ThisPatientsStrata[stratum]] += (DType)1 - SurvivalProb;
							SURVIVE.Threaded_AttackRates[th][WholeTrial][SEVERE_HOSP	][ThisPatientsStrata[stratum]] += (DType)1 - SurvivalProb_Severe;
						}

					if (SURVIVE.CalculateSeparatePassivePhaseCurves && HOUSE.ActiveOrPassivePhase == DO_ACTIVE_AND_PASSIVE && daypostdose > NoActiveDays_ThisPatient) ///// need last two conditions here (and not in TREATED_EQUALLY loop below) as not automatically true here, but they are in loop below. 
						if (!DATA.LTFU[patient]) // i.e. only include patients who were not cases in the active phsae and who have not been lost to follow up. 
							if (!DATA.IsCase_AMandPS[ActiveMild][patient] && !(DATA.FollowUp[END][ActiveMild][patient] < DATA.EndActiveSurveillance[patient])) //// need second condition here (and not in TREATED_EQUALLY loop below)
							{
								//// this would make calculations start on 0th day, which is actually last day of active. 
								//// Must start instead on 1st day of passive, as this is what SurvivalProb will refer to in loop (i.e. will be one day out otherwise). 
								//// PassiveDayPostDose = daypostdose - (NoActiveDays_ThisPatient + 1);
								PassiveDayPostDose = daypostdose - NoActiveDays_ThisPatient;

								if (PassiveDayPostDose <= SURVIVE.NoPassiveDaysFollowUp) //// i.e. don't exceed memory, and keep constant 
								{
									/// add to apprppriate group population size. 
									for (int stratum = 0; stratum < NumStrata_thisPatient; stratum++)	SURVIVE.SC_PP.Threaded_Strata_Sizes[th][ThisPatientsStrata[stratum]][PassiveDayPostDose]++;

									/// add to apprpriate groups
									for (int stratum = 0; stratum < NumStrata_thisPatient; stratum++)
									{
										SURVIVE.SC_PP.Threaded_PostSample_SurvivalTables[th][EITHER			][ThisPatientsStrata[stratum]][PassiveDayPostDose] += (SurvivalProb_Either	/ SurviveActive_Either	);
										SURVIVE.SC_PP.Threaded_PostSample_SurvivalTables[th][MILD_NON_HOSP	][ThisPatientsStrata[stratum]][PassiveDayPostDose] += (SurvivalProb			/ SurviveActive			);
										SURVIVE.SC_PP.Threaded_PostSample_SurvivalTables[th][SEVERE_HOSP	][ThisPatientsStrata[stratum]][PassiveDayPostDose] += (SurvivalProb_Severe	/ SurviveActive_Severe	);
									}
								}
						}

					//// add to appropriate sum of mean hazards - later divide by i) denominator to get mean and ii) other mean hazard ratios to get ratio
					if (ConditionToCalcHazRatio)
					{
						BaseHazValue	= PARAMS.BaselineHazardValues[country][PreviousDay];
						if (DATA.Vi_s[patient] == 1)	
							VacHazValue		= VaccineHazard(PreviousDay, DATA.FollowUp[START][ActiveMild][patient], DATA.SecondDose[patient], DATA.ThirdDose[patient], DATA.ai_s[patient], DATA.Ii_s[patient], country, PARAMS, HOUSE);
						WhichDayIndex = -6; for (int dummyindex = 0; dummyindex < SURVIVE.HRs_DaysPostDose.size(); dummyindex++)
							if (SURVIVE.HRs_DaysPostDose[dummyindex] == (daypostdose - 1)) WhichDayIndex = dummyindex;

						for (int HR_serotype = 0; HR_serotype < SURVIVE.HRs_NumSTypes; HR_serotype++)
						{
							if (HR_serotype == 0)
							{
								BaseHazMult_Mild	= K_rho_AM		;	//// i.e. AUG.iBaseHaz_Mult	[th][ActiveMild		][DATA.Ii_s[patient]];
								VacHazMult_Mild		= K_Eff_rho_AM	;	//// i.e. AUG.iVacHaz_Mult	[th][ActiveMild		][DATA.Ii_s[patient]];
								BaseHazMult_Severe	= K_rho_PS		;	//// i.e. AUG.iBaseHaz_Mult	[th][PassiveSevere	][DATA.Ii_s[patient]];
								VacHazMult_Severe	= K_Eff_rho_PS	;	//// i.e. AUG.iVacHaz_Mult	[th][PassiveSevere	][DATA.Ii_s[patient]];
							}
							else
							{
								BaseHazMult_Mild	= AUG.HR_BaseHaz_seroMults	[th][ActiveMild		][DATA.Ii_s[patient]][HR_serotype - 1]; //// note the - 1. This is due do HR_serotype == 0 refering to combined serotypes. 
								VacHazMult_Mild		= AUG.HR_VacHaz_seroMults	[th][ActiveMild		][DATA.Ii_s[patient]][HR_serotype - 1]; //// note the - 1. This is due do HR_serotype == 0 refering to combined serotypes. ;
								BaseHazMult_Severe	= AUG.HR_BaseHaz_seroMults	[th][PassiveSevere	][DATA.Ii_s[patient]][HR_serotype - 1]; //// note the - 1. This is due do HR_serotype == 0 refering to combined serotypes. 
								VacHazMult_Severe	= AUG.HR_VacHaz_seroMults	[th][PassiveSevere	][DATA.Ii_s[patient]][HR_serotype - 1]; //// note the - 1. This is due do HR_serotype == 0 refering to combined serotypes. ;
							}
							///// add multipliers for "either"
							BaseHazMult_Either	= BaseHazMult_Mild	+ BaseHazMult_Severe;	
							VacHazMult_Either	= VacHazMult_Mild	+ VacHazMult_Severe;	

							if (DATA.Vi_s[patient] == 0)
							{
								HazValue_stratum_Mild	= BaseHazMult_Mild		* BaseHazValue;
								HazValue_stratum_Severe	= BaseHazMult_Severe	* BaseHazValue;
								HazValue_stratum_Either = BaseHazMult_Either	* BaseHazValue;
							}
							else
							{
								HazValue_stratum_Mild	= (BaseHazMult_Mild		* BaseHazValue) - (VacHazMult_Mild		* VacHazValue);
								HazValue_stratum_Severe	= (BaseHazMult_Severe	* BaseHazValue) - (VacHazMult_Severe	* VacHazValue);
								HazValue_stratum_Either = (BaseHazMult_Either	* BaseHazValue) - (VacHazMult_Either	* VacHazValue);
							}

							//// add to apprpriate groups
							for (int stratum = 0; stratum < NumStrata_thisPatient; stratum++)
							{
								SURVIVE.MeanHazVals_BySerotype[HR_serotype].Threaded_PostSample_SurvivalTables[th][MILD_NON_HOSP][ThisPatientsStrata[stratum]][WhichDayIndex] += HazValue_stratum_Mild;
								SURVIVE.MeanHazVals_BySerotype[HR_serotype].Threaded_PostSample_SurvivalTables[th][SEVERE_HOSP	][ThisPatientsStrata[stratum]][WhichDayIndex] += HazValue_stratum_Severe;
								SURVIVE.MeanHazVals_BySerotype[HR_serotype].Threaded_PostSample_SurvivalTables[th][EITHER		][ThisPatientsStrata[stratum]][WhichDayIndex] += HazValue_stratum_Either;
							}
						}


					}
				}		//// end daypostdose loop 
			}
			else
			{
				/////////     /////////     /////////     /////////     /////////     /////////     /////////     /////////     /////////     
				/////////     ACTIVE PHASE LOOP
				IntVacHaz_AM_RunningTotal = 0;

				for (int daypostdose = 1; daypostdose <= NoActiveDays_ThisPatient; daypostdose++) 
				{
					ConditionToCalcHazRatio = (std::any_of(SURVIVE.HRs_DaysPostDose.begin(), SURVIVE.HRs_DaysPostDose.end(), [&](int i) {return i == (daypostdose - 1); })); //// if day is equal to one of HRs_DaysPostDose that we want to calculate. 
					if (CHAINS.HazRatiosOnly & !ConditionToCalcHazRatio) continue;

					/// add to apprppriate group population size. 
					for (int stratum = 0; stratum < NumStrata_thisPatient; stratum++) 
						SURVIVE.SC_WT.Threaded_Strata_Sizes[th][ThisPatientsStrata[stratum]][daypostdose]++;

					ThisDay					= DATA.FollowUp[START][ActiveMild][patient] + daypostdose;
					PreviousDay				= DATA.FollowUp[START][ActiveMild][patient] + daypostdose - 1;
					AM_IntBaseHaz_ThisDay	= IntBaseHaz(DATA.FollowUp[START][ActiveMild][patient], ThisDay, country, PARAMS) * PartitionHazardMultiplier;

					if (DATA.Vi_s[patient] == ControlGroup) 	SurvivalProb = exp(-AM_IntBaseHaz_ThisDay * K_rho_AM);	
					else 
					{
						if (!HOUSE.AllDosesRequired_AG_BS[DATA.AgeGroup1[patient]]	[DATA.Ii_s[patient]] || PreviousDay >= DATA.ThirdDose[patient]) //// i.e. only add to the vaccine hazard if past the third dose, or if not all three doses are required for this patient's serostatus (latter is default)
							IntVacHaz_AM_RunningTotal += IntVacHaz_SingleDay(PreviousDay, DATA.FollowUp[START][ActiveMild][patient], DATA.SecondDose[patient], DATA.ThirdDose[patient], DATA.ai_s[patient], DATA.Ii_s[patient], country, PARAMS, HOUSE);
						SurvivalProb				= exp((K_Eff_rho_AM * IntVacHaz_AM_RunningTotal) - (K_rho_AM * AM_IntBaseHaz_ThisDay));
					}

					/// add to apprpriate groups
					for (int stratum = 0; stratum < NumStrata_thisPatient; stratum++)
						SURVIVE.SC_WT.Threaded_PostSample_SurvivalTables[th][EITHER][ThisPatientsStrata[stratum]][daypostdose] += SurvivalProb;

					//// add to appropriate sum of mean hazards - later divide by i) denominator to get mean and ii) other mean hazard ratios to get ratio
					if (ConditionToCalcHazRatio)
					{
						BaseHazValue	= PARAMS.BaselineHazardValues[country][PreviousDay];
						if (DATA.Vi_s[patient] == VaccineGroup)	
							VacHazValue		= VaccineHazard(PreviousDay, DATA.FollowUp[START][ActiveMild][patient], DATA.SecondDose[patient], DATA.ThirdDose[patient], DATA.ai_s[patient], DATA.Ii_s[patient], country, PARAMS, HOUSE);
						WhichDayIndex = -6; for (int dummyindex = 0; dummyindex < SURVIVE.HRs_DaysPostDose.size(); dummyindex++)
							if (SURVIVE.HRs_DaysPostDose[dummyindex] == (daypostdose - 1)) WhichDayIndex = dummyindex;

						for (int HR_serotype = 0; HR_serotype < SURVIVE.HRs_NumSTypes; HR_serotype++)
						{
							if (HR_serotype == 0)
							{
								BaseHazMult_Either	= K_rho_AM		;	//// i.e. AUG.iBaseHaz_Mult	[th][ActiveMild][DATA.Ii_s[patient]];
								VacHazMult_Either	= K_Eff_rho_AM	;	//// i.e. AUG.iVacHaz_Mult	[th][ActiveMild][DATA.Ii_s[patient]];
							}
							else
							{
								BaseHazMult_Either	= AUG.HR_BaseHaz_seroMults	[th][ActiveMild][DATA.Ii_s[patient]][HR_serotype - 1]; 
								VacHazMult_Either	= AUG.HR_VacHaz_seroMults	[th][ActiveMild][DATA.Ii_s[patient]][HR_serotype - 1]; 
							}
							if (DATA.Vi_s[patient] == 0)	HazValue_stratum_Either =  BaseHazMult_Either * BaseHazValue;
							else							HazValue_stratum_Either = (BaseHazMult_Either * BaseHazValue) - (VacHazMult_Either * VacHazValue);

							//// add to apprpriate groups
							for (int stratum = 0; stratum < NumStrata_thisPatient; stratum++)
								SURVIVE.MeanHazVals_BySerotype[HR_serotype].Threaded_PostSample_SurvivalTables[th][EITHER][ThisPatientsStrata[stratum]][WhichDayIndex] += HazValue_stratum_Either;

						}
					}
				}		//// end ACTIVE daypostdose loop 

				///// add to active phase attack rate for appropriate groups. Needs to be thread-specific
				for (int stratum = 0; stratum < NumStrata_thisPatient; stratum++)
					SURVIVE.Threaded_AttackRates[th][ActivePhase][EITHER][ThisPatientsStrata[stratum]] += (DType)1 - SurvivalProb;

				DType SurviveActive = SurvivalProb; /// For passive phase attack rate calculation 

				/////////     /////////     /////////     /////////     /////////     /////////     /////////     /////////     /////////     
				/////////     PASSIVE PHASE LOOP
				if (HOUSE.ActiveOrPassivePhase == DO_ACTIVE_AND_PASSIVE) 
				{
					//// must add integrated hazards from hospital phase to integrated hazards from active phase. Format is always exp(- (Act_Haz_X + HospitalHazardUpToThisDayPostDose)) for all hazard groups. 
					DType IntVacHaz_PS_RunningTotal = 0, PS_IntBaseHaz_ThisDay;

					///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
					///// pre-compute  integrated hazards for  active phase for each patient
					///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
			
					DType AM_IntHaz = NULL;  //// this will refer to vaccine group or control group depending on patient. And already has K multipliers etc.
							if (DATA.Vi_s[patient] == ControlGroup) AM_IntHaz = K_rho_AM *  AM_IntBaseHaz_ThisDay;
					else	if (DATA.Vi_s[patient] == VaccineGroup) AM_IntHaz = (K_rho_AM * AM_IntBaseHaz_ThisDay) - (K_Eff_rho_AM * IntVacHaz_AM_RunningTotal);

					////// check that SurviveActive is the same as exp(-AM_IntHaz) as it should be!
					if ((SurviveActive != exp(-AM_IntHaz)) & (NoActiveDays_ThisPatient != 0))
						std::cout << "p" << patient << " S_Act " << SurviveActive << " e(-hz) " << exp(-AM_IntHaz) << " diff " << SurviveActive - exp(-AM_IntHaz) << endl;

					int LastDay_This_Patient = (HOUSE.LTFU_SurvivalCurves) ? NoActivePlusPassiveDays_ThisPatient : SURVIVE.NoDaysOfFollowUp; 
					 
					for (int daypostdose = NoActiveDays_ThisPatient + 1; daypostdose <= LastDay_This_Patient; daypostdose++)
					{
						ConditionToCalcHazRatio = (std::any_of(SURVIVE.HRs_DaysPostDose.begin(), SURVIVE.HRs_DaysPostDose.end(), [&](int i) {return i == (daypostdose - 1); })); //// if day is equal to one of HRs_DaysPostDose 
						if (CHAINS.HazRatiosOnly & !ConditionToCalcHazRatio) continue;

						/// add to apprppriate group population size.  
						for (int stratum = 0; stratum < NumStrata_thisPatient; stratum++)	SURVIVE.SC_WT.Threaded_Strata_Sizes[th][ThisPatientsStrata[stratum]][daypostdose]++;

						ThisDay					= DATA.FollowUp[START][ActiveMild][patient] + daypostdose;
						PreviousDay				= DATA.FollowUp[START][ActiveMild][patient] + daypostdose - 1;
						PS_IntBaseHaz_ThisDay	= IntBaseHaz(DATA.FollowUp[START][PassiveSevere][patient], ThisDay, country, PARAMS) * PartitionHazardMultiplier;
						if (DATA.Vi_s[patient] == 0) 	SurvivalProb = exp(-(AM_IntHaz + (PS_IntBaseHaz_ThisDay * K_rho_PS)));	
						else 																									
						{
							IntVacHaz_PS_RunningTotal	+= IntVacHaz_SingleDay(PreviousDay, DATA.FollowUp[START][ActiveMild][patient], DATA.SecondDose[patient], DATA.ThirdDose[patient], DATA.ai_s[patient], DATA.Ii_s[patient], country, PARAMS, HOUSE); 
							SurvivalProb				=  exp(-(AM_IntHaz + (K_rho_PS * PS_IntBaseHaz_ThisDay) - (K_Eff_rho_PS * IntVacHaz_PS_RunningTotal)    ));
						}

						/// add to apprpriate groups
						for (int stratum = 0; stratum < NumStrata_thisPatient; stratum++)
							SURVIVE.SC_WT.Threaded_PostSample_SurvivalTables[th][EITHER][ThisPatientsStrata[stratum]][daypostdose] += SurvivalProb;

						///// add to active phase attack rate for appropriate groups. Needs to be thread-specific
						if (daypostdose == NoActivePlusPassiveDays_ThisPatient)
							for (int stratum = 0; stratum < NumStrata_thisPatient; stratum++)
							{
								SURVIVE.Threaded_AttackRates[th][PassivePhase	][EITHER][ThisPatientsStrata[stratum]] += SurviveActive - SurvivalProb; 
								SURVIVE.Threaded_AttackRates[th][WholeTrial		][EITHER][ThisPatientsStrata[stratum]] += (DType)1 - SurvivalProb; 
							}

						if (SURVIVE.CalculateSeparatePassivePhaseCurves) 
							if (!DATA.LTFU[patient] && !DATA.IsCase_AMandPS[ActiveMild][patient]) //// i.e. we only include patients who were not cases in the active phsae and who have not been lost to follow up.
							{
								PassiveDayPostDose = daypostdose - NoActiveDays_ThisPatient ;

								if (PassiveDayPostDose <= SURVIVE.NoPassiveDaysFollowUp) //// i.e. don't exceed memory, and keep constant 
								{
									/// add to apprppriate group population size. 
									for (int stratum = 0; stratum < NumStrata_thisPatient; stratum++)	SURVIVE.SC_PP.Threaded_Strata_Sizes[th][ThisPatientsStrata[stratum]][PassiveDayPostDose]++;

									/// add to apprpriate groups
									for (int stratum = 0; stratum < NumStrata_thisPatient; stratum++)
										SURVIVE.SC_PP.Threaded_PostSample_SurvivalTables[th][EITHER][ThisPatientsStrata[stratum]][PassiveDayPostDose] += (SurvivalProb / SurviveActive);
								}
							}

						//// add to appropriate sum of mean hazards - later divide by i) denominator to get mean and ii) other mean hazard ratios to get ratio
						if (ConditionToCalcHazRatio)
						{
							BaseHazValue	= PARAMS.BaselineHazardValues[country][PreviousDay];
							if (DATA.Vi_s[patient] == 1)	
								VacHazValue	= VaccineHazard(PreviousDay, DATA.FollowUp[START][ActiveMild][patient], DATA.SecondDose[patient], DATA.ThirdDose[patient], DATA.ai_s[patient], DATA.Ii_s[patient], country, PARAMS, HOUSE);
							WhichDayIndex = -6; for (int dummyindex = 0; dummyindex < SURVIVE.HRs_DaysPostDose.size(); dummyindex++)
								if (SURVIVE.HRs_DaysPostDose[dummyindex] == (daypostdose - 1)) WhichDayIndex = dummyindex;
							for (int HR_serotype = 0; HR_serotype < SURVIVE.HRs_NumSTypes; HR_serotype++)
							{
								if (HR_serotype == 0)
								{
									BaseHazMult_Either	= K_rho_PS		;	//// i.e. AUG.iBaseHaz_Mult	[th][PassiveSevere][DATA.Ii_s[patient]];
									VacHazMult_Either	= K_Eff_rho_PS	;	//// i.e. AUG.iVacHaz_Mult	[th][PassiveSevere][DATA.Ii_s[patient]];
								}
								else
								{
									BaseHazMult_Either	= AUG.HR_BaseHaz_seroMults	[th][PassiveSevere][DATA.Ii_s[patient]][HR_serotype - 1]; 
									VacHazMult_Either	= AUG.HR_VacHaz_seroMults	[th][PassiveSevere][DATA.Ii_s[patient]][HR_serotype - 1]; 
								}
								if (DATA.Vi_s[patient] == 0)	HazValue_stratum_Either =  BaseHazMult_Either * BaseHazValue;
								else							HazValue_stratum_Either = (BaseHazMult_Either * BaseHazValue) - (VacHazMult_Either * VacHazValue);

								//// add to apprpriate groups
								for (int stratum = 0; stratum < NumStrata_thisPatient; stratum++)
									SURVIVE.MeanHazVals_BySerotype[HR_serotype].Threaded_PostSample_SurvivalTables[th][EITHER][ThisPatientsStrata[stratum]][WhichDayIndex] += HazValue_stratum_Either;
							}
						}
					}		//// end PASSIVE daypostdose loop 
				}			//// end if DO_ACITVE_AND_HOSPITAL if statement

			}	//// End if MildAndSevere == TREATED_SEPARATELY else statement (i.e. MildAndSevere == TREATED_EQUALLY)
			if (AllCool == 0) break; 

		}			//// end patient loop 
		delete[] ThisPatientsStrata	; 
		std::fflush(stderr);
	}		//// end pragma loop 
	std::fflush(stderr);

	///// Add to strata sizes from different threads
	AddToStrataSizes_fromThreads(DATA, HOUSE, SURVIVE);

	//// Add to Survival tables for this posterior sample from different threads
	AddToSurvivalTables_fromThreads(HOUSE, SURVIVE);
	
	//// Divide survival probabilities by appropriate denominators (if non-zero), and add to add to ParamSet_SurvivalTable and MetaParamSet_SurvivalTable
	AddToSurvivalTailsAndMean(SURVIVE, HOUSE);

	//// Calculate HRs: sum from threads; divide by strata sizes; calculate ratios; add to mean posterior sample and tails; reset threads. 
	Process_MeanHaz_And_HRs(HOUSE, SURVIVE); 

	/// Add to attack rates from the different threads. 
	AddToAttackRates_fromThreads(AR_PostSampleNo, HOUSE, SURVIVE, AttackRateContainer);

	//// (reset and) sum the FU_duration for this trial phase and stratum. 
	//// Divide MetaAttackRates (which is currently only the numerator) by Sum_FU_duration_phase_stratum(i.e.the denominator). Then Convert To Percentages
	AddToSum_FU_Durations_fromThreads_AndCalculateAttackRates(AR_PostSampleNo, HOUSE, SURVIVE, AttackRateContainer);

	if (OutputIndividualSurvivalTables)
	{
		if (CleanFirstDay_IndividualCurves) /// don't appear to need to do this with passive phase. 
			for (int DiseaseSeverity = 0; DiseaseSeverity < SURVIVE.HowManyDiseaseSeverities; DiseaseSeverity++)
				for (int stratum = 0; stratum < SURVIVE.NoSubjectCategories; stratum++)
					SURVIVE.SC_WT.PostSample_SurvivalTables[DiseaseSeverity][stratum][0] = 1; 
		for (int DiseaseSeverity = 0; DiseaseSeverity < SURVIVE.HowManyDiseaseSeverities; DiseaseSeverity++)
			WriteSurvivalTable(CHAINS.OutputFolder + ParamSet_SCurve_RootFilename + HOUSE.OutputString + SURVIVE.DiseaseNames[DiseaseSeverity] + ".txt",
				SURVIVE.SC_WT.PostSample_SurvivalTables[DiseaseSeverity], SURVIVE.NoSubjectCategories, SURVIVE.NoDaysOfFollowUp + 1, SURVIVE.SurvivalTableRowNames, SURVIVE.SurvivalTableColNames);
	}
	
	if (OutputIndividual_Passive_SurvivalTables)
		if (SURVIVE.CalculateSeparatePassivePhaseCurves)
			for (int DiseaseSeverity = 0; DiseaseSeverity < SURVIVE.HowManyDiseaseSeverities; DiseaseSeverity++)
				WriteSurvivalTable(CHAINS.OutputFolder + ParamSet_PassiveSCurve_RootFilename + HOUSE.OutputString + SURVIVE.DiseaseNames[DiseaseSeverity] + ".txt",
					SURVIVE.SC_PP.PostSample_SurvivalTables[DiseaseSeverity], SURVIVE.NoSubjectCategories, SURVIVE.NoPassiveDaysFollowUp + 1, SURVIVE.SurvivalTableRowNames, SURVIVE.PassivePhaseSurvivalTableColNames);

	//// Add to the counter
	SURVIVE.PosteriorSampleCounter++;
}

void CalculateAttackRates					(DType **SurvivalTable, DType **AttackRates, const Housekeeping_Struct &HOUSE, Survival_Struct &SURVIVE, int &iter, const Chains_Struct &CHAINS, int FinalDay, DType Duration)
{
	//// go through all "final day" survival probabilites in SurvivalTable and calculate attack rate
#pragma omp parallel for schedule(static,1)
	for (int th = 0; th < HOUSE.max_threads; th++)
		for (int stratum = th; stratum < SURVIVE.NoSubjectCategories; stratum += HOUSE.max_threads)
			AttackRates[stratum][(iter - CHAINS.BurnIn) / (SURVIVE.AddToSurvivalCurvesEveryHowManyIterations)] =
				((((DType)1) - SurvivalTable[stratum][FinalDay]) / Duration) * ((DType)100);
}
void CalculateAttackRates_PassiveOnly		(DType **SurvivalTable, DType **AttackRates, const Housekeeping_Struct &HOUSE, Survival_Struct &SURVIVE, int &iter, const Chains_Struct &CHAINS, int FirstDay, int FinalDay, DType Duration)
{
	//// go through all "final day" survival probabilites in SurvivalTable and calculate attack rate
#pragma omp parallel for schedule(static,1)
	for (int th = 0; th < HOUSE.max_threads; th++)
		for (int stratum = th; stratum < SURVIVE.NoSubjectCategories; stratum += HOUSE.max_threads)
			AttackRates[stratum][(iter - CHAINS.BurnIn) / (SURVIVE.AddToSurvivalCurvesEveryHowManyIterations)] =
				((SurvivalTable[stratum][FirstDay] - SurvivalTable[stratum][FinalDay]) / Duration) * ((DType)100);
}
void CalculateSurvivalMeansAndCrIs			(Survival_Struct &SURVIVE, DType **** FinalPosteriorSurvivalCurves, DType ***** CurveTails, int NoDaysOfFollowUp, int NumStrata, bool Passive)
{
	if (SURVIVE.NoSurvivePostSamples > 0)
	{

		//// idea of function is: Go through every "tail vector". That is, for both upper and lower tails, for all days post dose, and for all strata, find the min or max of the tail to get the credible interval. Add it to appropriate table

		int LowerIndex = max(0, int(SURVIVE.NumElementsOutside_CrI_Tails - 1)); /// pick last one of lower  first one of upper. 
		int UpperIndex = 0; 

		DType *ptrStart_LowerTail = NULL; //// address of tail array, will change in loop depending on category and daypostdose and DiseaseSeverity. 
		DType *ptrStart_UpperTail = NULL; //// address of tail array, will change in loop depending on category and daypostdose and DiseaseSeverity. 

		int DiseaseSeverity = EITHER; 

		int StartingDay = Passive ? 0 : 1;

		for (int DiseaseSeverity = 0; DiseaseSeverity < SURVIVE.HowManyDiseaseSeverities; DiseaseSeverity++)
			for (int category = 0; category < NumStrata; category++)
			{
				if (!Passive) 
				{
					FinalPosteriorSurvivalCurves[DiseaseSeverity][MEAN_POST		][category][0] = 1;
					FinalPosteriorSurvivalCurves[DiseaseSeverity][LOWER_CrI_POST][category][0] = 1;
					FinalPosteriorSurvivalCurves[DiseaseSeverity][UPPER_CrI_POST][category][0] = 1;
				}

				for (int daypostdose = StartingDay; daypostdose < NoDaysOfFollowUp; daypostdose++)
				{
					//// divide the Mean
					FinalPosteriorSurvivalCurves[DiseaseSeverity][MEAN_POST][category][daypostdose] /= (SURVIVE.NoSurvivePostSamples);

					//// sort the tails
					ptrStart_LowerTail = CurveTails[DiseaseSeverity][LowerTail_Index][category][daypostdose]; // Makes sorting easier to read. 
					ptrStart_UpperTail = CurveTails[DiseaseSeverity][UpperTail_Index][category][daypostdose]; // Makes sorting easier to read. 

					sort(		ptrStart_LowerTail,		ptrStart_LowerTail + SURVIVE.NumElementsOutside_CrI_Tails		);
					sort(		ptrStart_UpperTail,		ptrStart_UpperTail + SURVIVE.NumElementsOutside_CrI_Tails		);

					//// add appropriate member of sorted vector to Lower Posterior Survival Curve. 
					FinalPosteriorSurvivalCurves[DiseaseSeverity][LOWER_CrI_POST][category][daypostdose] = ptrStart_LowerTail[LowerIndex];
					FinalPosteriorSurvivalCurves[DiseaseSeverity][UPPER_CrI_POST][category][daypostdose] = ptrStart_UpperTail[UpperIndex];
				}
			}
	}
}
void CalculateSurvivalCurveOutput			(Survival_Struct &SURVIVE)
{
#ifdef PRINT_PROGRAM_PROGRESS
	std::cerr << "CalculateSurvivalCurveOutput: "; 
#endif
	//// Whole Trial survival curves
	CalculateSurvivalMeansAndCrIs(		SURVIVE, SURVIVE.SC_WT.FinalPosteriorSurvivalCurves	, SURVIVE.SC_WT.CurveTails, SURVIVE.NoDaysOfFollowUp + 1		, SURVIVE.NoSubjectCategories, false	);
	
	//// Passive Phase survival curves
	if (SURVIVE.CalculateSeparatePassivePhaseCurves)
		CalculateSurvivalMeansAndCrIs(	SURVIVE, SURVIVE.SC_PP.FinalPosteriorSurvivalCurves	, SURVIVE.SC_PP.CurveTails, SURVIVE.NoPassiveDaysFollowUp + 1	, SURVIVE.NoSubjectCategories, false	); 

	//// Hazard ratios
	for (int HR_serotype = 0; HR_serotype < SURVIVE.HRs_NumSTypes; HR_serotype++)
		CalculateSurvivalMeansAndCrIs(SURVIVE, SURVIVE.HRs_BySerotype[HR_serotype].FinalPosteriorSurvivalCurves, SURVIVE.HRs_BySerotype[HR_serotype].CurveTails, SURVIVE.HRs_NumDaysPostDoseToCalculate, SURVIVE.HRs_NumStrata, true);
#ifdef PRINT_PROGRAM_PROGRESS
	std::cerr << "DONE" << endl; ;
#endif
}

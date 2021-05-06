
#ifndef STRUCT_DEFNS_HEADER_INCLUDED
#define STRUCT_DEFNS_HEADER_INCLUDED

#include "stdafx.h"		
#include "Macros.h"

typedef double DType; // for flexibility only in case ever need long doubles everywhere.

template <typename TYPE> void Allocate_2D_Array		(TYPE **    &OBJECT, int Dim1, int Dim2)
{
	OBJECT = new TYPE *[Dim1]();
	for (int row = 0; row < Dim1; row++)	OBJECT[row] = new TYPE[Dim2]();
}
template <typename TYPE> void Allocate_3D_Array		(TYPE ***   &OBJECT, int Dim1, int Dim2, int Dim3)
{
	OBJECT = new TYPE **[Dim1]();
	for (int row = 0; row < Dim1; row++)	Allocate_2D_Array(OBJECT[row], Dim2, Dim3);
}
template <typename TYPE> void Allocate_4D_Array		(TYPE ****  &OBJECT, int Dim1, int Dim2, int Dim3, int Dim4)
{
	OBJECT = new TYPE ***[Dim1]();
	for (int row = 0; row < Dim1; row++)	Allocate_3D_Array(OBJECT[row], Dim2, Dim3, Dim4);
}
template <typename TYPE> void Allocate_5D_Array		(TYPE ***** &OBJECT, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5)
{
	OBJECT = new TYPE ****[Dim1]();
	for (int row = 0; row < Dim1; row++)	Allocate_4D_Array(OBJECT[row], Dim2, Dim3, Dim4, Dim5);
}

template <typename T> bool IsEqual_1D_Arrays	(T*    const &lhs, T*    const &rhs, int dim)
{
	bool All_Equal = 1;
	if (lhs != 0 && rhs != 0)	All_Equal = (memcmp(lhs, rhs, sizeof(T) * dim) == 0);		//// if both arrays non-empty compare them
	else if (lhs == 0 && rhs == 0)	All_Equal = 1;											//// if both arrays empty return true
	else							All_Equal = 0; 											//// if one array empty and other not then return false. 
	return All_Equal;
}
template <typename T> bool IsEqual_2D_Arrays	(T**   const &lhs, T**   const &rhs, int dim1, int dim2)
{
	bool All_Equal = 1; 
	if (lhs != 0 && rhs != 0)														//// if both arrays non-empty compare them
	{
		for (int row = 0; row < dim1; row++)
			if (IsEqual_1D_Arrays(lhs[row], rhs[row], dim2) == false) { All_Equal = 0; break; }
	}
	else if (lhs == 0 && rhs == 0)
	{
		All_Equal = 1; 																//// if both arrays empty return true
	}
	else All_Equal = 0;																//// if one array empty and other not then return false. 

	return All_Equal; 
}
template <typename T> bool IsEqual_3D_Arrays	(T***  const &lhs, T***  const &rhs, int dim1, int dim2, int dim3)
{
	bool All_Equal = 1;
	if (lhs != 0 && rhs != 0)
	{
		for (int row = 0; row < dim1; row++)
			if (IsEqual_2D_Arrays(lhs[row], rhs[row], dim2, dim3) == false) { All_Equal = 0; break; }
	}
	else if (lhs == 0 && rhs == 0)
	{
		All_Equal = 1; 																//// if both arrays empty return true
	}
	else All_Equal = 0;																//// if one array empty and other not then return false.
	return All_Equal;
}
template <typename T> bool IsEqual_4D_Arrays	(T**** const &lhs, T**** const &rhs, int dim1, int dim2, int dim3, int dim4)
{
	bool All_Equal = 1;
	if (lhs != 0 && rhs != 0)
	{
		for (int row = 0; row < dim1; row++)
			if (IsEqual_3D_Arrays(lhs[row], rhs[row], dim2, dim3, dim4) == false) { All_Equal = 0; break; }
	}
	else if (lhs == 0 && rhs == 0)
	{
		All_Equal = 1; 																//// if both arrays empty return true
	}
	else All_Equal = 0;																//// if one array empty and other not then return false.
	return All_Equal;
}

template <typename T> void SetEqual_1D_Arrays	(T*    &ArrayToChange, T*    const &ArrayToCopy, int dim) //// assumes ArrayToChange and ArrayToCopy have same dimensions. 
{
	if (ArrayToChange != 0) if (ArrayToCopy != 0) for (int element = 0; element < dim; element++) ArrayToChange[element] = ArrayToCopy[element];
}
template <typename T> void SetEqual_2D_Arrays	(T**   &ArrayToChange, T**   const &ArrayToCopy, int dim1, int dim2) //// assumes ArrayToChange and ArrayToCopy have same dimensions. 
{
	if (ArrayToChange != 0) if (ArrayToCopy != 0) for (int element = 0; element < dim1; element++) SetEqual_1D_Arrays(ArrayToChange[element], ArrayToCopy[element], dim2);
}
template <typename T> void SetEqual_3D_Arrays	(T***  &ArrayToChange, T***  const &ArrayToCopy, int dim1, int dim2, int dim3)  //// assumes ArrayToChange and ArrayToCopy have same dimensions. 
{
	if (ArrayToChange != 0) if (ArrayToCopy != 0) for (int element = 0; element < dim1; element++) SetEqual_2D_Arrays(ArrayToChange[element], ArrayToCopy[element], dim2, dim3);
}
template <typename T> void SetEqual_4D_Arrays	(T**** &ArrayToChange, T**** const &ArrayToCopy, int dim1, int dim2, int dim3, int dim4)  //// assumes ArrayToChange and ArrayToCopy have same dimensions. 
{
	if (ArrayToChange != 0) if (ArrayToCopy != 0) for (int element = 0; element < dim1; element++) SetEqual_3D_Arrays(ArrayToChange[element], ArrayToCopy[element], dim2, dim3, dim4);
}

template <typename TYPE> void Populate_1D_Array(TYPE *		&OBJECT, TYPE Value, int Dim1)
{
	for (int row = 0; row < Dim1; row++)	OBJECT[row] = Value;
}
template <typename TYPE> void Populate_2D_Array(TYPE **		&OBJECT, TYPE Value, int Dim1, int Dim2)
{
	for (int row = 0; row < Dim1; row++) Populate_1D_Array(OBJECT[row], Value, Dim2);
}
template <typename TYPE> void Populate_3D_Array(TYPE ***	&OBJECT, TYPE Value, int Dim1, int Dim2, int Dim3)
{
	for (int row = 0; row < Dim1; row++) Populate_2D_Array(OBJECT[row], Value, Dim2, Dim3);
}

///////////////////////// ///////////////////////// ///////////////////////// ///////////////////////// ///////////////////////// ///////////////////////// ///////////////////////// ///////////////////////// 
///////////////////////// ///////////////////////// ///////////////////////// ///////////////////////// ///////////////////////// ///////////////////////// ///////////////////////// ///////////////////////// 
///////////////////////// Define Data Structures
///////////////////////// ///////////////////////// ///////////////////////// ///////////////////////// ///////////////////////// ///////////////////////// ///////////////////////// ///////////////////////// 


enum class Age_Option		{	INDEPENDENT	/*i.e. does not vary by age*/, HILL, CATEGORICAL /*varies by age group*/, SPLINE, SPLINE_LINE, SPLINE_STEP, CUBIC }; // if changing, change Convert_AS_String function too. 
enum class ExtImSub_Option	{	IGNORED		/*i.e. Imputed values ignored and model as before*/, 
								AS_DATA		/*i.e. Imputed values maximum likelihood estimates of serostatus*/, 
								AS_PROB		/*i.e. Imputed values simulations used as Prob(SeroPos) where available for Gibbs augmentation (but not in MH MCMC*/}; // if changing, change Convert_ExtImSub_String function too. 
enum class EffNegWane_Option {	DEFAULT		/*i.e. Even negative efficacies increase by "waning" to zero)*/,
								FROM_ZERO	/*i.e. Negative efficacies start from zero and decline to that negative efficacy. */,
								NO_WANE		/*i.e. Negative efficacies don't wane at all: they stay at their initial value. These are canonical results.*/ };


struct DATA_struct {

	///// Follow Up dates / indices	 
	int ***FollowUp					= NULL;						//// indexed by i) START or END (0 or 1 hash defined), ii) PhaseSeverity, iii) patient number. Units are in days. 
	int *EndActiveSurveillance		= new int [NPat]();			//// For predicted survival curves when using Sanofi follow up definitions. Usual follow up variables can be altered by patient being a case, but for survival curves need to know if patient was under active or passive surveillance. 
	int *EndPassiveSurveillance		= new int [NPat]();			//// For predicted survival curves when using Sanofi follow up definitions. Usual follow up variables can be altered by patient being a case, but for survival curves need to know if patient was under active or passive surveillance. 
	bool *LTFU						= new bool[NPat]();			//// For predicted survival curves when using Sanofi follow up definitions. Usual follow up variables can be altered by patient being a case, but for survival curves need to know if patient was under active or passive surveillance. 

	int **FU_Duration_days			= NULL;						//// indexed by i) Trial Phase (active or passive or active + passive); ii) patient number. Units are in days. 
	DType **FU_Duration_years		= NULL;						//// indexed by i) Trial Phase (active or passive or active + passive); ii) patient number. Units are in years (hence DType not int). 

	int *SecondDose					= new int	[NPat]();
	int *ThirdDose					= new int	[NPat]();
	int *TimePost_1st_Dose			= new int	[NPat]();
	int *TimePost_Final_Dose		= new int	[NPat]();

	//// Demographic characteristics
	int *ai_s						= new int	[NPat]();		//// age in years
	int	*Ii_s						= new int	[NPat]();		//// Baseline immunity - keep as int because of missing data. MDValue = 9999 not acceptable as boolean or char
	int	*Imputed_Ii_s				= new int	[NPat]();		//// Baseline immunity - Or more accurately imputed estimates of baseline immunity combined with the immune-subset. 
	DType * Imputed_ProbSPos		= new DType	[NPat]();		//// Baseline immunity - Based on 1000 imputations
	bool *Vi_s						= new bool	[NPat]();		//// Trial Arm - 0 is control group, 1 is vaccine group
	int *ci_s						= new int	[NPat]();		//// country of person i
	
	//// IsCase
	bool *IsCase					= new bool	[NPat]();		//// was patient a case? 
	bool *IsCase_ActiveMild			= new bool	[NPat]();		//// If distinguishing between mild and severe, this will refer to mild		disease
	bool *IsCase_PassiveSevere		= new bool	[NPat]();		//// If distinguishing between mild and severe, this will refer to severe	disease
	bool **IsCase_AMandPS			= NULL;						//// indexed by i) PhaseSeverity, iii) Patient number
	int *Case_PhaseSeverity			= new int	[NPat]();		//// indexed by i) patient. 

	//// Are patients in Immune subset? Simple true/false used in Survival curve generation
	//// We loop over AugmentedIndices and NonAugmentedIndices, but we don't have an easy label for each patient. Initialized to zero (i.e. not in immune subset). Changed when read in data. If data != MDValue, then patient is in Immune Sub set and so value InImmSub[patient] = 1. 
	bool * InImmSub					= new bool	[NPat]();		//// indexed by i) patient. Must be initialized to zero for everybody. Correct when reading in data. 
	

	int *AgeGroup1					= new int	[NPat]();
	int *AgeGroup2					= new int	[NPat]();
	int *CaseSerotype				= NULL;

	int NoDataVariables				= 22;

	DType **CountryMinMaxCalendarTime	= new DType*[NumC](); 
	int *NumCalendarDaysFollowUp		= new int	[NumC]();		//// Number of days of follow-up for each country. Affected by last day of passive/severe followup (and first day of active/mild followup if not starting from Day 0 for each country). 

	int NoDaysOfFollowUp;		//// Refers to maximum follow up duration over all patients. Will be different depending on whether you model Active only or Active+Passive, whether you use Sanofi's follow up definitions, and whether you include late cases. 
	int NoActiveDaysFollowUp;
	int NoAugmented = 0, NoInfected = 0; //// add to these as you read in data. They account for skipped countries. Note that NoInfected is required to calculate probabilites of cases, as opposed to densities. For every case (mild and severe, active and hospital phase), you calculate the density, multiplied by the time interval (1 day, 1/365) to get a probability. Therefore, you do not need to distinguish between how many of the infected are active/mild or hospital/severe, except when you count how many cases there are (e.g. don't want to include hospital phase cases if only modelling active). 

	std::vector<int> ** Set_Array		= new std::vector<int> *[NumC]();	//// indexed by i) country, ii) Hazard Group (includes pooled see hash defines for definitions), iii) person (note that number/order specific to country, hazard group - e.g. if country = 5 and Hazard Group = 1 then person 0 does not refers to 1st person in that group, not first person in trial who would be in country = 0 ;
	std::vector<int> *** Cases_Array	= new std::vector<int> **[NumC]();	//// indexed by i) country, ii) Hazard Group (inc. pooled), iii) PhaseSeverity; iv) person (note that number/order specific to country, hazard group - e.g. if country = 5 and Hazard Group = 1 then person 0 does not refers to 1st person in that group, not first person in trial who would be in country = 0 ;
	int MaxStartFollowUp_CalendarTime	= 0; //// need this for survival curve calculation. You need (integrated) baseline hazard values for all patients for the duration of the trial (whether or not they were a case who therefore dropped out of follow up is irrelevant for predicted survival curves). Therefore you need hazard values from 0 to (Latest Start Time + Trial Duration)
	int MaxEndFollowUp_CalendarTime		= 0;

	std::vector<int> PatientIndexNumbers; //// just 1:NPat. Useful for Survival curves. When doing full trial, will loop over this quantity. When doing immune subset only, will loop over CHAINS.NonAugmentedIndices (which is also a vector and hence switching between two should be easy). Will make more sense in thread loop / in patient loop within GenerateSurvivalCurves function.  

	//// Used for waning values when fitting PASSIVE_PHASE_ONLY. 
	//// If fitting only passive phase, day e.g. 1000 calendar time is redefined to be 1000 - Offset_From_Day0, where Offset_From_Day0 is the earliest date that a subject enters passive surveillance (if active surveillance, Offset_From_Day0 = 0).
	//// However if day 1000 calendar time corresponds to say 900 days of follow up for patient i, you still need the waning value to reflect 900 days, not 900 minus some offset. 
	//// Solution: Subtract Offset_From_Day0 from dose dates as well as follow up dates. So that difference will be still be 900. This leaves only problem of waning value indicies. 
	//// Previously, 3rd waning index referred to 3 days of waning after most recent dose. For passive phase, there are only 335 days of follow up, but you will need waning values for say 900th day post dose. 
	//// Solution: Keep 900th waning day as index 900, even though this will mean calculating waning values that you won't need. Could offset this... 
	//// ... but then this would complicate the IntVacHaz function which would presumably slow down the rest of your code for model variants you actually care about. 
	//// Not true for integrated hazard as you condition on having
	//// Therefore you need to add to exp(-t/Tau) (or Hill/waning function equivalent) to be exp(-t_adjusted + Offset_From_Day0 / Tau). If (!PASSIVE_PHASE_ONLY), leave defined as zero. 

	int N_WaningDays; /// Will usually be same as DATA.NoDaysFollowUp except for PASSIVE_PHASE_ONLY, where for reasons in comment above this will be separate. You use DATA.NoDaysFollowUp to allocate memory for survival tables and (integrated) hazard tables. Therefore good to separate the two. 
	int Offset_From_Day0 = 0;

	std::vector<int> AugmentedIndices, NonAugmentedIndices;
};
struct LikeIndices_Struct {

	////// Following categories of index: Augmentation, IntBaseHaz, IntVacHaz, VacIntBaseHaz, K's, WaningEffs, rhos. Anything associated with survival functions needs an Active/Mild and Passive/Severe Phase/Severity version. 

	int HugeLikeIndex = 41000;		////  pick something that will always be bigger than the length of the LikePartsArray. Don't have the non Pointers as NULL as NULL is evaluated to zero. 

	int *BSs					= NULL;				//// ///// Prob (SNeg) and Prob (SPos) LL components: indexed by i) baseline serostatus.

	int *AgeHazMult				= NULL;				//// indexed by i) PhaseSeverity; 

	int ** rhos					= NULL;				//// indexed by i) serotype ii) PhaseSeverity (NOT country) 	
	//int *** InfEffs			= NULL;				//// indexed by i) serotype ii) PhaseSeverity. iii) BaselineSeroStatus. Used for HOUSE.ResidEffs. (i.e. log(1 - VE_ inf)) part of likelihood for cases	). 
	int *** IntVacHaz			= NULL;				//// indexed by i) PhaseSeverity, ii) Baseline serostatus and iii) serotype
	int *** WaningEffs			= NULL;				//// indexed by i) serotype ii) baseline serostatus iii) PhaseSeverity

	int *** Ks					= NULL;				//// indexed by i) serotype, ii) PhaseSeverity, iii) No. prior infections (K2 taken to mean 2 or 3 prior infections). 
	int ** Kplus				= NULL;				//// indexed by i) serotype, ii) PhaseSeverity
	int ** KPrimes				= NULL;				//// For HOUSE.ModelVariant == AS_PRIME. indexed by i) serotype, ii) PhaseSeverity 
	int ** KplusPrime			= NULL;				//// For HOUSE.ModelVariant == AS_PRIME. indexed by i) serotype, ii) PhaseSeverity

	int *BaseHaz				= NULL;				//// indexed by i) PhaseSeverity; 
	int *** IntBaseHaz			= NULL;				//// indexed by i) PhaseSeverity; ii) Baseline serostatus and iii) trial arm. 

	int **l_BS_BaseHaz_Mults	= NULL;				//// indexed by i) PhaseSeverity; ii) Baseline serostatus 
};
struct ParamRangeStruct {

	DType KAM_0						[2] = { 0,0 };
	DType KAM_1						[2] = { 0,0 };
	DType KAM_2						[2] = { 0,0 };
	DType KPS_0						[2] = { 0,0 };
	DType KPS_1						[2] = { 0,0 };
	DType KPS_2						[2] = { 0,0 };
	DType Hosp_K_multiplier			[2] = { 0,0 };
	DType SNegEff_1					[2] = { 0,0 };
	DType SPosEff_1					[2] = { 0,0 };
	DType SeroNegWaning				[2] = { 0,0 };
	DType SeroPosWaning				[2] = { 0,0 };
	DType HillHalfLife_SNeg			[2] = { 0,0 };
	DType HillHalfLife_SPos			[2] = { 0,0 };
	DType knots						[2] = { 0,0 };
	DType knots_logged				[2] = { 0,0 };
	DType Late_knots				[2] = { 0,0 };
	DType Late_knots_logged			[2] = { 0,0 };
	DType SeroNegWaning_Flipped		[2] = { 0,0 };
	DType SeroPosWaning_Flipped		[2] = { 0,0 };
	DType HillPower_SNeg			[2] = { 0,0 };
	DType HillPower_SPos			[2] = { 0,0 };
	DType HistHaz					[2] = { 0,0 };
	DType ASVEHalflife				[2] = { 0,0 };
	DType ASVEPower					[2] = { 0,0 };
	DType ASVEProp					[2] = { 0,0 };
	DType AS_Haz_Halflife			[2] = { 0,0 };
	DType AS_Haz_Power				[2] = { 0,0 };
	DType AS_Haz_Prop				[2] = { 0,0 };
	DType ASWaning_Halflife			[2] = { 0,0 };
	DType ASWaning_Power			[2] = { 0,0 };
	DType qval						[2] = { 0,0 };
	DType rho						[2] = { 0,0 };
	DType BS_BaseHazMult			[2] = { 0,0 };


};
struct Housekeeping_Struct {

	//// set defaults
	char ModelVariant			= VAC_SILENT				;		//// choose between SIMPLE_ANALYTICAL, SIMPLE_NUMERICAL, K_SEROPOS, DROP_K, and VAC_SILENT. 
	char SingleOrMultiDose		= MULTI_DOSE				;
	char ActiveOrPassivePhase	= DO_ACTIVE_AND_PASSIVE		;
	char MildAndSevere			= TREATED_EQUALLY			;
	bool PS_Ks_Multiply_AM_Ks	= false						;		//// if true, then PassiveSevere Ks are given by KAM_i x KPS_i. if false (default), PassiveSevere Ks are given by KPS_i.
	char LinKnts				= 0							;		
	char HowManyCaseCategories	= 1							;		
	char Aug_MH_or_Gibbs		= GIBBS_AUG					;
	bool FitWaningRate			= 0							;
	bool SeroSpecificEfficacies	= 0							;
	bool SeroSpecific_K_values	= 0							;

	bool SS_KAM_0				= true						;		//// relevant only if SeroSpecific_K_values == true. Should KAM_0s be serotype specific?
	const bool SS_KAM_1			= true						;		//// relevant only if SeroSpecific_K_values == true. Should KAM_1s be serotype specific? NEVER change this as it is the baseline. 
	bool SS_KAM_2				= true						;		//// relevant only if SeroSpecific_K_values == true. Should KAM_2s be serotype specific?
	bool SS_KPS_0				= true						;		//// relevant only if SeroSpecific_K_values == true. Should KPS_0s be serotype specific?
	bool SS_KPS_1				= true						;		//// relevant only if SeroSpecific_K_values == true. Should KPS_1s be serotype specific?
	bool SS_KPS_2				= true						;		//// relevant only if SeroSpecific_K_values == true. Should KPS_2s be serotype specific?

	bool **SSKs_FitMatrix		= NULL						;		//// indexed by i) PhaseSeverity; ii) PrevInf. Do you want K_PhaseSeverity_PrevInf to be serotype specific? Only relevant if SeroSpecific_K_values == true; 

	bool BaselinePartition		= 0							;		//// i.e. 1 + p2 + p3 + p4 where p1 = 1 is the baseline and the others are relative multipliers of baseline hazard. Could also have p1 + 1 + p3 + p4 if serotype 2 is the baseline. 
	int  Baseline_Serotype		= 0							;		//// i.e. serotype 1 for Cpp notation. For serotypes 3 and 4 (2 and 3 Cpp), there aren't always cases in every country, so more dangerous to set as baseline (i.e. zero hazard is legitimate choice and so no multiplier would work). 
	bool UseSyntheticData		= 0							; 
	bool HillWaning				= 0							;		//// i.e. Waning values described by Hill Function, not exponential. 
	bool ModellingHospitalized	= 0							;		//// If true, then have differnt Passive/Severe K's for CYD-15. Can easily extend to other countries. 
	bool ModelHosp_Indie_Ks		= 1							;		//// ModellingHospitalized variable overides this. 
	bool Fixed_Severe_RelRisks  = 0							; 		//// InitializeHousekeeping sets HOUSE.Skip_KH_0 = 1 and HOUSE.Skip_KH_2 = 1 if true. Default is false. 
	int  FSKs_Ratio_SetNum		= 1							;		//// Previously set K0/K1 = 1/4 and K2/K1 = 1/16 (0 = default). Now assigned to 1, where K0/K1 = 1/8 and K2/K1 = 1/16 
	bool Weighting_Pass_Sev		= 0							;		//// Will you weight data from passive phase - Don't be tempted to name this Weighted_PassiveSevere because of PassiveSevere #define macro.
	DType PS_Weight				= 2							;		
	bool Run_MCMC				= 1							;		//// This also encompasses SCurves etc. Only thing that happens if this is set to false is total likelihood evaluated for old chains (provided  set various switches properly). 
	bool WBIC_RunOnly			= 0							; 
	bool Do_WBIC_Runs			= 1							;
	bool SSASVE_Additive		= false						;		//// do serotypes interact multiplicitavely or additively with age-efficacy profile? Default (false) is multiplicative. 
	bool AdditiveSSASVEs		= false						;		//// This looks identical to SSASVE_Additive but it isn't. SSASVE_Additive can be set to true (say when reading in params), but this will only be effective if SS_VEs = TRUE and ASVE != Age_Option::INDEPENDENT. This switch ensures that the others are also flicked and therefore cuts down on irritating code repeated. 

	EffNegWane_Option EffNegWane = EffNegWane_Option::DEFAULT;		///// Because (by default) you apply multiplier of exp(-t/Tau_b) to initial value of efficacy/transient immunity, efficay diminishes over time if initial value positive but INCREASES over time if intial value negative. With this switch, you change waning s.t. if "efficacy" negative, then it starts at zero and declines to whatever the efficacy value is. 

	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** 
	///// Age-specific (AS) effects

	Age_Option ASVE			= Age_Option::INDEPENDENT	;	//// Default is no age dependence
	Age_Option AS_Haz		= Age_Option::INDEPENDENT	;	//// Default is no age dependence
	Age_Option AS_Waning	= Age_Option::INDEPENDENT	; 	//// Default is no age dependence
	bool ASVE_AdditionalKnots		= false; 
	bool AS_Haz_AdditionalKnots		= false; 
	int AS_Waning_KnotSet			= 0;	//// zero is default - corresponding to 2, 6, and 12 years (and 16 years if categorical/step). 1 (used to be used when old variable AS_Waning_AdditionalKnots = true, now set to non-binary). 1 correspondes to 2, 5.9, 6.1, 11.9, 12.1 and 16 (i.e. extra knots put around middle ones). 2 corresponds to 6 knots equidistantly spaced between 2 and 16 years old, i.e. 2, 4.8, 7.6, 10.4, 13.2, 16.0. 

	int NumASVE_Function_Params		= 3									;	/// using a Hill function where parameters are "halflife" and power, and now Proportion of Vaccine Efficacy that is dependent upon age. 
	int NumAS_Haz_Function_Params	= 3									;	/// using a Hill function where parameters are "halflife" and Power, and Proportion of Vaccine Efficacy that is dependent upon age. See Initialize housekeeping for further details. 
	int NumWaningParamsPer_BS		= 1									;	//// Number of waning parameters per baseline serostatus; Default is 1, i.e. where you model either the duration (default), or the waning rate. 
	
	int Num_AS_Priming_ParamsPer_BS = 2;	//// i) Proportion ii) Rate


	//// Spline values
	int	MaxSplineDegree_EffMultiplier	= 2, KnotsPerSpline_EffMultiplier	= 4, PolynomialsPerSpline_EffMultiplier		= 2;
	int	MaxSplineDegree_HazMultiplier	= 2, KnotsPerSpline_HazMultiplier	= 4, PolynomialsPerSpline_HazMultiplier		= 2;
	int	MaxSplineDegree_WaningDuration	= 2, KnotsPerSpline_WaningDuration	= 4, PolynomialsPerSpline_WaningDuration	= 2;
	 
	bool AgeEffectsSame_Waning	= 0;		//// fix age effects? identical to AS_Waning == Age_Option == Age_Option::INDEPENDENT, 
	bool AgeEffectsSame_VE		= 0;		//// fix age effects? identical to ASVE ==  Age_Option == Age_Option::INDEPENDENT, 
	bool ASVE_FitAllSero_VEs	= 0;		//// If doing age efficacy categoricals and splines, then you fix SNegEff_1 and SPosEff_1 equal to 1. To obtain SS_VEs model from SSASVEs model, can fix all ASVE params (and multipliers) to 1, but fit all the SSVEs (i.e. not fix SNegEff_1 and SPosEff_1 to 1). 

	bool AS_VE_Homogeneous		= false;		//// Will you allow the relationship between efficacy and age to differ by baseline serostatus - default is yes. 
	int NumSeroStatuses_ASVEs	= NULL						; 
	int ASVE_BS					= SeroNeg					;		//// used if HOUSE.ASVE_OnlyOneSeroStatus == true
	bool ASVE_OnlyOneSeroStatus = 0							;		//// Do you only want to apply the Age Specific Vaccine Efficacy function/multiplier to only one serostatus. If so, ASVE_Mults applied only to serostatus given by HOUSE.ASVE_BS (other serostatus multipliers are 1. )

	bool AS_Waning_OnlyOneSeroStatus = false				;
	int AS_Waning_OneSeroBS		= SeroNeg					;		//// if AS_Waning_OnlyOneSeroStatus, which serostatus? Default is SeroNeg.
	int NumAges_Waning			= 1							;		//// Default is 1, unless doing age-specific waning. Want to distinguish between number of ages in study (HOUSE.HowManyAges) and number of ages as far as waning concerned.  
	bool AS_Waning_Homogeneous	= false;		//// Will you allow the relationship between efficacy and age to differ by baseline serostatus - default is yes. 
	int NumSeroStatuses_AS_Waning = NULL;


	bool AdjHaz					= 0							;		//// Include multiplier of baseline hazard for particular baselineserostatus (given by HOUSE.Which_BS_BaseHazMult)? Baseline hazard is sum of all hazards for all serotypes. But seropositive individuals have been exposed to at least one serotype already, and therefore have reduced hazard. This boolean decides if you should inclde a multiplier of the baseline hazard for seropositives.  
	int Which_BS_BaseHazMult	= SeroPos					;



	bool Empirical_SeroPrevs	= false	;
	bool AreWeAugmenting		= 1		;			///// default is true. If false, then CHAINS.AreWeAugmenting and WBIC_CHAINS.AreWeAugmenting set to false. 


	bool PSVEs					= false						;
	int NumEffsPer_BS_And_SType = 1							;		//// Need this for param number functions. i.e. if you model different efficacy values for different trial phases (if PSVEs == true). If so, then value is HowManyCaseCategories (and so must be reset if HowManyCaseCategories is)
	int * MinMaxPhaseSeveritiesToLoopOver = new int[2]()	;		//// For phase severity dependant vaccine efficacies. 
	bool PASSIVE_PHASE_ONLY		= false						;

	bool ResidEffs				= false						;		//// i.e. Efficacy declines not to zero but to some other VE_min


	bool AllDosesRequired_SNeg		= false					;		//// Are all three doses reqired before any vaccine efficacy? These values are only used for data input (stored in AllDosesRequired_BS array), and output string. These values trump those that populate AllDosesRequired_AG_BS array. 
	bool AllDosesRequired_SPos		= false					;		//// Are all three doses reqired before any vaccine efficacy? These values are only used for data input (stored in AllDosesRequired_BS array), and output string. These values trump those that populate AllDosesRequired_AG_BS array. 

	bool AllDosesRequired_SNeg_AgeGroup1 = false			;
	bool AllDosesRequired_SNeg_AgeGroup2 = false			;
	bool AllDosesRequired_SNeg_AgeGroup3 = false			;
	bool AllDosesRequired_SNeg_AgeGroup4 = false			;
	bool AllDosesRequired_SNeg_AgeGroup5 = false			;
	bool AllDosesRequired_SPos_AgeGroup1 = false			;
	bool AllDosesRequired_SPos_AgeGroup2 = false			;
	bool AllDosesRequired_SPos_AgeGroup3 = false			;
	bool AllDosesRequired_SPos_AgeGroup4 = false			;
	bool AllDosesRequired_SPos_AgeGroup5 = false			;

	bool **AllDosesRequired_AG_BS = NULL					;		//// indexed by i) AgeGroup (technically DATA.AgeGroup1) and ii) BaselineSeroStatus. 

	bool SingleEff						= false				;		//// Will you fit a single efficacy parameter, or at least keep VE- and VE+ the same? 

	ExtImSub_Option ExtImSub			= ExtImSub_Option::IGNORED;

	string SerotypeString		= "";
	string KnotsInputFilename	= "MeanKnots_Initial.txt"			;
	string HH_InputFilename		= "HistHazardValues.txt"			;
	string DefaultParamRangeFileName	= "ParamRanges_Default"		; 
	string ParamRangeFileName			= DefaultParamRangeFileName	; //// useful for ParamRangeFileName to be separate from DefaultParamRangeFileName in CreateOutputString function. 

	ParamRangeStruct IntialParRanges; 

	string PhaseSeverity_strings[2] = { "AM", "PS" };

	bool Single_SNeg_Eff			= 0;		//// Will use this for checking. Looks same as SeroSpecificEfficacies == false, but isn't as want to keep all	 Seronegative	 efficacies same but keep rest of MCMC machinery for serotype specific efficacies.
	bool Single_SPos_Eff			= 0;		//// Will use this for checking. Looks same as SeroSpecificEfficacies == false, but isn't as want to keep all	 Seropositive	 efficacies same but keep rest of MCMC machinery for serotype specific efficacies.
	int  Which_SNeg_SeroFitsAll		= 0;		//// i.e. do you want to set seronegative efficacies for serotype 2,3,4 equal to efficacy for serotype 1, or e.g. 1,3,4 equal to stype 2 etc. 
	int  Which_SPos_SeroFitsAll		= 0;		//// i.e. do you want to set seropositive efficacies for serotype 2,3,4 equal to efficacy for serotype 1, or e.g. 1,3,4 equal to stype 2 etc. 
	bool RelRisksSameBtwSerotypes	= 0;		//// used to see if SS_Ks improves likelihood. Want to keep machinery for SS_Ks (e.g. rhos param update functions etc.) without actually allowing the model additional flexibility. In this way we do not compare apples and oranges. 
	int  Which_Sero_FitsAll			= 0; 
	
	bool SFU						= true;		//// SFU = Sanofi_FollowUp_Definitions
	bool Include_Late_Cases			= false; 
	bool FakeExtObs					= false;
	bool MakeEveryoneAControl		= false;	//// obviously for checking only. 
	bool MakeEveryoneAVaccinee		= false;	//// obviously for checking only. 
	
	
	bool PooledCountries		= false; //// so all CYD-14 countries aggregated together and all CYD-15 aggregated together
	bool PooledTrials			= false; 

	int Combined_CYD_14_country			= 0;  //// are we combining all CYD14	countries as a single country (i.e. if PooledCountries == true	). If so they are all combined as country Combined_CYD_14_country	;	
	int Combined_CYD_15_country			= 5;  //// are we combining all CYD15	countries as a single country (i.e. if PooledCountries == true	). If so they are all combined as country Combined_CYD_15_country	;	
	int Combined_CYD_14_15_country		= 0;  //// are we combining all			countries as a single country (i.e. if PooledTrials == true). If so they are all combined as country Combined_CYD_14_15_country;		

	bool StartFromPreviousChain				= false	;
	bool StartFromPreviousChain_Aug			= false	;
	bool StartFrom_WBIC_Chain				= false	; //// start from regular chains until you fix pathological maximum likelihood associated with the augmentation. 
	bool AddOutputStringExtraToOldChains	= false	; //// Set to false if e.g. you want use old chains as the starting point for something new (this new thing would require an OutputStringExtra to distinguish it from old input). Set to true if you want to resume chains where you previously used OutputStringExtra
	string OutputStringExtra				= ""	; //// useful if you want to try out things (e.g. when debugging). 
	string OutputStringExtra_PrevChain		= ""	; //// Anytime where you want to use a previous chain, you need to take what (if anything) was different about that chain, and then make a duplicate of that variable in this structure, with the _PrevChain suffix. In your code, when you create new HOUSE structure from which you derive OutputString of Previous chain, you change value of original (e.g. OutputStringExtra) to value of _PrevChian (e.g. OutputStringExtra_PrevChain). In your R code you will set OutputStringExtra_PrevChain to equal OutputStringExtra by default, so that you will only overwrite if you actually want to. 
	string OutputStringForOldChainInput		= ""	; //// 
	string OutputStringForOldChainInput_Aug	= ""	; //// For checking. To fully compare WBICs etc., want to keep augmented data constant between and within model runs. Need to import Parameter values from previous chains (model-specific) but use only one set of augmented data - maximum likelihood
	string OldChainFileNamePrefix_Aug				= "FinalState_AugData"; //// choose from i) "FinalState_AugData" (default); ii) "ModalPost_AugData"; iii) "MaxLike_AugData". Can add in WBICs etc if need be. 

	bool FitAllGlobalParams = 1;
	bool Skip_PosWaning = 0, Skip_NegWaning = 0; //// refers to duration or rate (depending on HOUSE.FitWaningRate), and does not refer to age specific "duration" halflife or "duration" power.
	bool Skip_PosWaning_HLife = 0, Skip_NegWaning_HLife = 0; //// refers to AS_Waning parameters, i.e. these hill functions are for functions of vaccine efficacay by age, and not strictly by day post dose. 
	bool Skip_PosWaning_Power = 0, Skip_NegWaning_Power = 0; //// refers to AS_Waning parameters, i.e. these hill functions are for functions of vaccine efficacay by age, and not strictly by day post dose. 

	bool Skip_KA_0 = 0, Skip_KA_2 = 0, Skip_KH_0 = 0, Skip_KH_1 = 0, Skip_KH_2 = 0, Skip_Hosp_Mult = 0; //// all parameters fitted by default. 
	bool Skip_PosEfficacies = 0, Skip_NegEfficacies = 0, Skip_qvals = 0, SkipRelativeRisks = 0, SkipKnots = 0, SkipHistHazards = 0, Skip_Rhos = 0;
	bool Skip_All_ASVEs = 0, Skip_ASVEsPowers = 0, Skip_ASVE_Halflives = 0, Skip_ASVE_Props = 0, Skip_SNeg_ASVEs = 0, Skip_SPos_ASVEs = 0;
	bool Skip_All_AS_Haz = 0, Skip_AS_Haz_Powers = 0, Skip_AS_Haz_Halflives = 0, Skip_AS_Haz_Props = 0;
	bool Skip_BS_BaseHazMult = 1; ///// default should be 1 (i.e. to skip), as you don't want to fit this. This is just here in case you decide that you do want to fit the multiplier of the baseline hazard for a particular serostatus. 
	bool Skip_AS_Prime_All = 0, Skip_AS_Prime_SNegRate = 0, Skip_AS_Prime_SPosRate = 0, Skip_AS_Prime_SNegProp = 0, Skip_AS_Prime_SPosProp = 0; 
	bool FitAllCountries = 1, Fit_c0 = 1, Fit_c1 = 1, Fit_c2 = 1, Fit_c3 = 1, Fit_c4 = 1, Fit_c5 = 1, Fit_c6 = 1, Fit_c7 = 1, Fit_c8 = 1, Fit_c9 = 1; //// all countries fitted by default. 

	int NoCountriesToFit; //// NoCountriesToFit is how many countries you'll actually fit. 
	const int TotalCountries = 10; //// TotalCountries is total countries in the two trials, (i.e. 10). Does NOT account for skipped countries. 
	int TotalPatients = 0; //// add to this as you read in data (accounting for skipped countries).  
	std::vector<int> WhichCountries, CYD_14_countries, CYD_15_countries, Diff_HospK_s_Countries; 

	int	MaxSplineDegree = 2, KnotsPerCountry, PolynomialsPerCountry;
	int	max_threads;
	int	HowManyAges = 18, HowManyAgeGroups = 6; ///// 6 because there are 6 age groups (for DATA.AgeGroup1, won't work otherwise). 0 (all ages) won't be used and MUST be set to false. 
	int HowManySeroStatuses = 2, HowManyTrialArms = 2, HowManyTrialPhases = 1;
	int NumTrialArms_IBH = 2;	//// number of trial arms for l_IntBaseHaz. For SIMPLE_NUMERICAL and K_SEROPOS, K multipliers are the same between trial arms. Not for VAC_SILENT so you need another HowManySeroStatuses = 2 indices to make e.g. K2 param updates easier. 
	int **Strata = NULL;		//// indexed by i) Trial Arm; ii) BaselineSerostatus. Matrix of various strata. NOT used for survival curves (annoyingly). 

	bool LTFU_SurvivalCurves = 0; 
	bool Use_WithoutCases_FU_Defns = 0; 

	///// NumKs refers to number of PARAMETERS. Num_K_Likes refers to indices. Good to keep these separate as for example there is no log(K2) term in K_SEROPOS model variant, but there are still three K params per trial phase and serotype (baseline notwithstanding). 
	int Num_K_Params	= 0	; /// number of relative risks per serotype and TrialPhase/Disease Severity. DOES NOT INCLUDE K+ values, which are composites of these Ks and historical hazards and patients' ages.  Set to zero as default as this will then cause a crash if not set correctly.  
	int Num_K_Likes		= 0	; /// DOES NOT INCLUDE Kplus values. Number of K mulpiliers per serotype and TrialPhase/Disease Severity In VAC_SILENT, have l(K0), l(K1) (for passive) and l(K2), and l. For K_SEROPOS, have l(K0) only. For Simple numerical, have  l(K0) and l(K1)(for passive). Set to zero as default as this will then cause a crash if not set correctly.  
	int Num_Kplus_Likes = 0	;

	int	LCPerC, No_Parameters; //// LCPerC stands for "LikelihoodComponentsPerCountry". 
	
	LikeIndices_Struct L_Indices;

	int N_STypes					= 1;		//// Num serotypes (used for rho / serotype proportion parameters)
	int N_STypes_Ks					= 1;		//// Num serotypes for	relative risks	(helpful to separate this from N_STypes as if modelling SS_VEs, then 4 serotypes but not as far as	K			parameters concerned)
	int N_STypes_VEs				= 1;		//// Num serotypes for	efficacies		(helpful to separate this from N_STypes as if modelling SS_Ks , then 4 serotypes but not as far as	Efficacy	parameters concerned)
	int KS_SSVEsKs_BaselineSerotype = 0;		//// For K_SEROPOS SS_VEs and SS_Ks. Only need one IntVacHaz component and the rest are zero. This decides which serotype is the non-zero one. Better notes are in l_Int_Vac_Haz function and l_Int_Vac_Haz_sero_i function. Think choice or 0,1,2 or 3 is completely arbitrary: but worth checking. 
	int SSASVEs_BaselineSerotype	= 0;		//// For SSASVEs where SSASVE_Additive == true. Then only need one IntVacHaz component and the rest are zero. This decides which serotype is the non-zero one. Better notes are in l_Int_Vac_Haz function and l_Int_Vac_Haz_sero_i function. Think choice or 0,1,2 or 3 is completely arbitrary: but worth checking. 

	DType	TimeInterval = (DType(1) / DType(365));	
	DType	TimeIntervalReciprocal = DType(365);	

	string DataFilename = ""; 
	string OutputString = "";

	bool OutputSeroPrev_LLChains	= true; 
	bool RandomCases				= 0;
	int RandCase_Seed				= 448712963;
	long seed1 = 547838717ul, seed2 = 943517526ul;

	void init_DosesRequired			()
	{
		Allocate_2D_Array(AllDosesRequired_AG_BS, HowManyAgeGroups, HowManySeroStatuses); ///// 6 because there are 6 age groups (for DATA.AgeGroup1, won't work otherwise). 0 (all ages) won't be used and MUST be set to false. 

		if (AllDosesRequired_SNeg)
			for (int AgeGroup = 1; AgeGroup < HowManyAgeGroups; AgeGroup++) 
				AllDosesRequired_AG_BS[AgeGroup][SeroNeg] = true;
		else 
		{
			if (AllDosesRequired_SNeg_AgeGroup1)		AllDosesRequired_AG_BS[1][SeroNeg] = true;
			if (AllDosesRequired_SNeg_AgeGroup2)		AllDosesRequired_AG_BS[2][SeroNeg] = true;
			if (AllDosesRequired_SNeg_AgeGroup3)		AllDosesRequired_AG_BS[3][SeroNeg] = true;
			if (AllDosesRequired_SNeg_AgeGroup4)		AllDosesRequired_AG_BS[4][SeroNeg] = true;
			if (AllDosesRequired_SNeg_AgeGroup5)		AllDosesRequired_AG_BS[5][SeroNeg] = true;
		}

		if (AllDosesRequired_SPos)
			for (int AgeGroup = 1; AgeGroup < HowManyAgeGroups; AgeGroup++) 
				AllDosesRequired_AG_BS[AgeGroup][SeroPos] = true;
		else 
		{
			//// !!!!!! **** !!!!! NOTE THAT IT STARTS FROM 1, not zero
			if (AllDosesRequired_SPos_AgeGroup1)		AllDosesRequired_AG_BS[1][SeroPos] = true;
			if (AllDosesRequired_SPos_AgeGroup2)		AllDosesRequired_AG_BS[2][SeroPos] = true;
			if (AllDosesRequired_SPos_AgeGroup3)		AllDosesRequired_AG_BS[3][SeroPos] = true;
			if (AllDosesRequired_SPos_AgeGroup4)		AllDosesRequired_AG_BS[4][SeroPos] = true;
			if (AllDosesRequired_SPos_AgeGroup5)		AllDosesRequired_AG_BS[5][SeroPos] = true;
		}
	}
	void Init_Skips					()
	{
		//// decide on skipping conditions and initializations for ASVE parameters
		if (ASVE != Age_Option::INDEPENDENT & !SeroSpecificEfficacies)
		{
			//// In these scenarios,  PosEfficacy and NegEfficacy variables (and the PARAMS.Efficacies array they go on to populate) won't be fitted, as ASVE params will deal with efficacy. 
			//// But to save a million if statements, it's easier to include them in the likelihood calculations etc. but simply fix them at 1 so they alter nothing. 
			if (ASVE_OnlyOneSeroStatus)
			{
					 if (ASVE_BS == SeroPos)	Skip_PosEfficacies = 1;
				else if (ASVE_BS == SeroNeg)	Skip_NegEfficacies = 1;
			}
			else
			{
				Skip_PosEfficacies = 1;
				Skip_NegEfficacies = 1;
			}
		}
		if (Fixed_Severe_RelRisks)	{ Skip_KH_0 = 1;  Skip_KH_2 = 1;}
		if (SingleEff)				{ Skip_NegEfficacies = 1; Skip_NegWaning = 1;  }

		if (AS_Waning != Age_Option::INDEPENDENT && AS_Waning_OnlyOneSeroStatus)
		{
			//// If only one serostatus is subject to age specific waning, then skip the HLife's and Power's associated with AS_Waning. Do not skip the duration
			if (AS_Waning_OneSeroBS == SeroNeg)
			{
				Skip_PosWaning_HLife = 1;				///// notice that	if (HOUSE.AS_Waning_OneSeroBS == SeroNeg) you skip the SEROPOSITIVE parameters, and ....
				Skip_PosWaning_Power = 1;
			}
			else if (AS_Waning_OneSeroBS == SeroPos)
			{
				Skip_NegWaning_HLife = 1;				/////		...		if (HOUSE.AS_Waning_OneSeroBS == SeroPos) you skip the SERONEGATIVE parameters
				Skip_NegWaning_Power = 1;
			}
			else std::cerr << "InitializeWaningParamsAndValues ERROR: AS_Waning_OneSeroBS = " << AS_Waning_OneSeroBS << " not recognized" << endl;
		}
		if (ModelVariant == AS_PRIME)
		{
			if (Skip_AS_Prime_All)
			{
				Skip_AS_Prime_SNegProp = 1;
				Skip_AS_Prime_SNegRate = 1;
				Skip_AS_Prime_SPosProp = 1;
				Skip_AS_Prime_SPosRate = 1;
			}
		}
		if (Single_SNeg_Eff & Single_SPos_Eff & SeroSpecificEfficacies)	Skip_qvals = 1;
	}
	void InitializeCountriesFitted	()
	{
		///// Initialize Countries fitted etc. 
		if (FitAllCountries == 1)  //// if fitting all countries don't have a non-empty CountriesFittedString (i.e. fitting all countries is the default. 
		{
			NoCountriesToFit = 10;
			for (int country = 0; country < NoCountriesToFit; country++)
			{
				WhichCountries.push_back(country);
				if (country < 5)	CYD_14_countries.push_back(country);			//// CYD-14
				else				CYD_15_countries.push_back(country);			//// CYD-15
			}
		}
		else 
		{
			NoCountriesToFit = 0;
			if (Fit_c0) {		WhichCountries.push_back(0);		NoCountriesToFit++;		CYD_14_countries.push_back(0);	}	//// CYD-14
			if (Fit_c1) {		WhichCountries.push_back(1);		NoCountriesToFit++;		CYD_14_countries.push_back(1);	}
			if (Fit_c2) {		WhichCountries.push_back(2);		NoCountriesToFit++;		CYD_14_countries.push_back(2);	}
			if (Fit_c3) {		WhichCountries.push_back(3);		NoCountriesToFit++;		CYD_14_countries.push_back(3);	}
			if (Fit_c4) {		WhichCountries.push_back(4);		NoCountriesToFit++;		CYD_14_countries.push_back(4);	}
			if (Fit_c5) {		WhichCountries.push_back(5);		NoCountriesToFit++;		CYD_15_countries.push_back(5);	}	//// CYD-15
			if (Fit_c6) {		WhichCountries.push_back(6);		NoCountriesToFit++;		CYD_15_countries.push_back(6);	}
			if (Fit_c7) {		WhichCountries.push_back(7);		NoCountriesToFit++;		CYD_15_countries.push_back(7);	}
			if (Fit_c8) {		WhichCountries.push_back(8);		NoCountriesToFit++;		CYD_15_countries.push_back(8);	}
			if (Fit_c9) {		WhichCountries.push_back(9);		NoCountriesToFit++;		CYD_15_countries.push_back(9);	}
		}

		if (PooledCountries)
		{
			WhichCountries.clear(); //// to get rid of previous 

			NoCountriesToFit	= 0;

			if (Fit_c0 || Fit_c1 || Fit_c2 || Fit_c3 || Fit_c4)
			{ 
				WhichCountries.push_back(Combined_CYD_14_country); 
				NoCountriesToFit++;		
				CYD_14_countries.push_back(Combined_CYD_14_country);
			}
			if (Fit_c5 || Fit_c6 || Fit_c7 || Fit_c8 || Fit_c9)
			{
				WhichCountries.push_back(Combined_CYD_15_country);
				NoCountriesToFit++;
				CYD_15_countries.push_back(Combined_CYD_15_country);
			}
		}
		else if (PooledTrials)
		{
			WhichCountries.clear(); //// to get rid of previous 

			NoCountriesToFit	= 1;
			WhichCountries.push_back(Combined_CYD_14_15_country);
			CYD_14_countries.push_back(Combined_CYD_14_15_country); //// think you want to use CYD_14_countries as this is unaffected by hospitalisation K's. 
		}

		if (PooledCountries)		//// note that mean knots only for SFU. 
		{
			KnotsInputFilename	= "MeanKnots_Initial_PooledCountries.txt";
			HH_InputFilename		= "HistHazardValues_PooledCountries.txt";
		}
		std::cerr << "Fitting countries: ";
		for (int country = 0; country < TotalCountries; country++)
		{
			bool condition;
			condition = std::any_of(WhichCountries.begin(), WhichCountries.end(), [&](int i) {	return i == country;	});
			if (condition) std::cerr << country << " ";
		}
#ifdef PRINT_PROGRAM_PROGRESS
		std::cerr << "NoCountriesToFit " << NoCountriesToFit << "/" << TotalCountries << ": FitAllCountries " << FitAllCountries << endl;
#endif
		std::cerr << "Fit Cs:"; for (int i = 0; i < WhichCountries.size(); i++) std::cerr << WhichCountries[i]; std::cerr << '\n';
	}
	void init_Strata				()
	{
		//// Initialize Strata
		Allocate_2D_Array(Strata, HowManyTrialArms + 1, HowManySeroStatuses + 1);
		int strata_counter = 0;
		for (int TrialArm = 0; TrialArm <= HowManyTrialArms; TrialArm++)	//// note the <= as want "either arm". 
			for (int BaselineSeroStatus = 0; BaselineSeroStatus <= HowManySeroStatuses; BaselineSeroStatus++)	//// note the <= as want "either serostatus". 
				Strata[TrialArm][BaselineSeroStatus] = strata_counter++;
	}
	void init_PhaseSevereties		()
	{
		if (PASSIVE_PHASE_ONLY)	ActiveOrPassivePhase = ACTIVE_PHASE_ONLY; ///// This line looks weird. But it is right. Essentially, if fitting only the passive phase, you want to use all the code's machinery that is dependent on their being only a single trial phase to fit (default would obviously be the active). Therefore set everything to use ACTIVE_PHASE_ONLY. This means that throughout the rest of the code, if there is anything that refers specifically to the passive phase (calendar start times, knots etc.) then you'll need to have if(PASSIVE_PHASE_ONLY) before any if (ActiveOrPassivePhase == ACTIVE_PHASE_ONLY) (e.g. where you set KnotsPerCountry below). 

			 if (ActiveOrPassivePhase == DO_ACTIVE_AND_PASSIVE) HowManyTrialPhases = 2; 
		else if (ActiveOrPassivePhase == ACTIVE_PHASE_ONLY	) HowManyTrialPhases = 1;

		if ((ActiveOrPassivePhase == DO_ACTIVE_AND_PASSIVE) || (MildAndSevere == TREATED_SEPARATELY))	HowManyCaseCategories = 2; 

		if (HowManyCaseCategories == 1)		Fixed_Severe_RelRisks = false;

		if (PSVEs) 
				NumEffsPer_BS_And_SType = HowManyCaseCategories;
		else	NumEffsPer_BS_And_SType = 1;

		MinMaxPhaseSeveritiesToLoopOver[0] = 0; 
		MinMaxPhaseSeveritiesToLoopOver[1] = HowManyCaseCategories; //// NOTE: Not PassiveSevere
		if (PASSIVE_PHASE_ONLY) PhaseSeverity_strings[0] = "PS"; //// i.e. because you're using ActiveOrPassivePhase == "ACTIVE_PHASE_ONLY"

	}
	void init_KnotsPolysPerCountry	()
	{
		if (PASSIVE_PHASE_ONLY)	KnotsPerCountry = 3; 
		else
		{
			if (SFU)
			{
						if (ActiveOrPassivePhase == ACTIVE_PHASE_ONLY		) KnotsPerCountry = 8;
				else	if (ActiveOrPassivePhase == DO_ACTIVE_AND_PASSIVE	) KnotsPerCountry = 11;
				else	std::cerr << "InitializeHousekeeping ERROR: ActiveOrPassivePhase argument not recognized" << endl;
			}
			else
			{
						if (ActiveOrPassivePhase == ACTIVE_PHASE_ONLY		) KnotsPerCountry = 6;
				else	if (ActiveOrPassivePhase == DO_ACTIVE_AND_PASSIVE	) KnotsPerCountry = 9;
				else	std::cerr << "InitializeHousekeeping ERROR: ActiveOrPassivePhase argument not recognized" << endl;
			}
		}
		PolynomialsPerCountry = KnotsPerCountry - 2;
	}
	void init_N_STypes				()
	{
		if (SeroSpecificEfficacies || SeroSpecific_K_values)
		{
			N_STypes = 4;
			if (SeroSpecificEfficacies	)  N_STypes_VEs	= 4;
			if (SeroSpecific_K_values	)  N_STypes_Ks	= 4;
		}
		else	N_STypes = 1;
	}
	void init_Num_K_Likes_Params	()
	{
		//// KParams
		if (ModelVariant == VAC_SILENT || ModelVariant == AS_PRIME || ModelVariant == K_SEROPOS)	Num_K_Params = 3; /// even though no log(K2) term in K_SEROPOS variant, still have K2 in K+ values. Also, KA1 counts
		else																						Num_K_Params = 2; /// in SIMPLE model, although there is only apparently K0, in Passive phase there is KH1. Also, KA1 counts 

		//// K-likes (different as in K_SEROPOS model, while you have three K params (per trial phase/ serotype/ trial), you do not ever have a log(K1) or log(K2) term. 
			 if (ModelVariant == VAC_SILENT			)	Num_K_Likes = 3; //// K0,K1,K2
		else if (ModelVariant == K_SEROPOS			)	Num_K_Likes = 1; //// K0 (no explicit K1 or K2 likelihood parts, just a Kplus part). 
		else if (ModelVariant == SIMPLE_NUMERICAL	)	Num_K_Likes = 2; //// K0,K1		
		else if (ModelVariant == AS_PRIME			)	Num_K_Likes = 1; //// K0 (no explicit K1 or K2 likelihood parts). Instead a 
		else std::cerr << "ModelVariant not recognized " << endl; 

		if (ModelVariant == VAC_SILENT || ModelVariant == K_SEROPOS || ModelVariant == AS_PRIME)	Num_Kplus_Likes = 1;	else Num_Kplus_Likes = 0;
	}
	void init_L_Indices				()
	{
#ifdef PRINT_PROGRAM_PROGRESS
		std::cerr << "Initialize_LikeStruture ";
#endif
		////// ////// ////// ////// ////// ////// ////// ////// ////// ////// 
		////// BUILD Likelihood indices etc. according to what you want to fit (e.g. Model Variant, trial phase, serospecific blah). 
		////// Don't want loads of string comparisons for say AmendProposedParams, but also can't use #defines due to parameter numbers varying by ModelVariant, trial phase and whether or not including sero-specific efficacies. 
		////// ////// ////// ////// ////// ////// ////// ////// ////// ////// 

		int LikeIndexCounter = 0;			//// index will be recorded with this, then incremented. 

		L_Indices.BaseHaz = new int[HowManyCaseCategories]();
		for (int PhaseSeverity = 0; PhaseSeverity < HowManyCaseCategories; PhaseSeverity++) L_Indices.BaseHaz[PhaseSeverity] = LikeIndexCounter++;

		///// Prob (SNeg) and Prob (SPos) LL components
		L_Indices.BSs = new int[HowManySeroStatuses]();
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HowManySeroStatuses; BaselineSeroStatus++)
			L_Indices.BSs[BaselineSeroStatus] = LikeIndexCounter++;
	
		///// IntBaseHaz Like indices. //// NumTrialArms_IBH refers to number of trial arms for l_IntBaseHaz. For SIMPLE_NUMERICAL and K_SERPOS, K multipliers are the same between trial arms. Not for VAC_SILENT so you need another HowManySeroStatuses = 2 indices to make e.g. K2 param updates easier. 
		Allocate_3D_Array(L_Indices.IntBaseHaz, HowManyCaseCategories, HowManySeroStatuses, NumTrialArms_IBH); 
		for (int PhaseSeverity = 0; PhaseSeverity < HowManyCaseCategories; PhaseSeverity++)
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HowManySeroStatuses; BaselineSeroStatus++)
				for (int TrialArm = 0; TrialArm < NumTrialArms_IBH; TrialArm++)
					L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus][TrialArm] = LikeIndexCounter++;

		///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
		///// Case (log) K values, including Kplus. 

		Allocate_3D_Array(L_Indices.Ks, N_STypes_Ks, HowManyCaseCategories, Num_K_Likes);
		//// K_likes
		for (int serotype = 0; serotype < N_STypes_Ks; serotype++)
			for (int PhaseSeverity = 0; PhaseSeverity < HowManyCaseCategories; PhaseSeverity++)
				for (int K_like = 0; K_like < Num_K_Likes; K_like++)
					L_Indices.Ks[serotype][PhaseSeverity][K_like] = LikeIndexCounter++;	

		//// K_plus_likes
		if (ModelVariant == VAC_SILENT || ModelVariant == K_SEROPOS || ModelVariant == AS_PRIME)
		{
			Allocate_2D_Array(L_Indices.Kplus, N_STypes_Ks, HowManyCaseCategories);
			for (int serotype = 0; serotype < N_STypes_Ks; serotype++)
				for (int PhaseSeverity = 0; PhaseSeverity < HowManyCaseCategories; PhaseSeverity++)
					L_Indices.Kplus[serotype][PhaseSeverity] = LikeIndexCounter++; 
		}

		if (ModelVariant == AS_PRIME)
		{
			Allocate_2D_Array(L_Indices.KPrimes		, N_STypes_Ks, HowManyCaseCategories);
			Allocate_2D_Array(L_Indices.KplusPrime	, N_STypes_Ks, HowManyCaseCategories);

			for (int serotype = 0; serotype < N_STypes_Ks; serotype++)
				for (int PhaseSeverity = 0; PhaseSeverity < HowManyCaseCategories; PhaseSeverity++)
				{
					//// KPrime_likes
					L_Indices.KPrimes		[serotype][PhaseSeverity] = LikeIndexCounter++;
					//// K_plus_Prime_likes
					L_Indices.KplusPrime	[serotype][PhaseSeverity] = LikeIndexCounter++;
				}
		}

		///// IntVacHaz indices
		Allocate_3D_Array(L_Indices.IntVacHaz, HowManyCaseCategories, HowManySeroStatuses, N_STypes_VEs);
		for (int PhaseSeverity = 0; PhaseSeverity < HowManyCaseCategories; PhaseSeverity++)
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HowManySeroStatuses; BaselineSeroStatus++)
				for (int serotype = 0; serotype < N_STypes_VEs; serotype++)
					L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype] = LikeIndexCounter++;

		///// waning indices
		Allocate_3D_Array(L_Indices.WaningEffs, N_STypes_VEs, HowManySeroStatuses, HowManyCaseCategories);
		for (int serotype = 0; serotype < N_STypes_VEs; serotype++)
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HowManySeroStatuses; BaselineSeroStatus++)
				for (int PhaseSeverity = 0; PhaseSeverity < HowManyCaseCategories; PhaseSeverity++)
					L_Indices.WaningEffs[serotype][BaselineSeroStatus][PhaseSeverity] = LikeIndexCounter++;

		///// rhos indices
		if ((SeroSpecificEfficacies) || (SeroSpecific_K_values))
		{
			Allocate_2D_Array(L_Indices.rhos, N_STypes, HowManyCaseCategories);
			for (int serotype = 0; serotype < N_STypes; serotype++)	
				for (int PhaseSeverity = 0; PhaseSeverity < HowManyCaseCategories; PhaseSeverity++)
					L_Indices.rhos[serotype][PhaseSeverity] = LikeIndexCounter++;
		}

		///// BS_BaseHaz_Mult index
		if (AdjHaz)
		{
			Allocate_2D_Array(L_Indices.l_BS_BaseHaz_Mults, HowManyCaseCategories, HowManySeroStatuses);
			for (int PhaseSeverity = 0; PhaseSeverity < HowManyCaseCategories; PhaseSeverity++)
				for (int BaselineSeroStatus = 0; BaselineSeroStatus < HowManySeroStatuses; BaselineSeroStatus++)
					L_Indices.l_BS_BaseHaz_Mults[PhaseSeverity][BaselineSeroStatus] = LikeIndexCounter++; //// SeroPos or SeroNeg
		}

		///// AgeHazMult index
		if (AS_Haz != Age_Option::INDEPENDENT)
		{
			L_Indices.AgeHazMult = new int[HowManyCaseCategories]();
			for (int PhaseSeverity = 0; PhaseSeverity < HowManyCaseCategories; PhaseSeverity++)
				L_Indices.AgeHazMult[PhaseSeverity] = LikeIndexCounter++; 
		}
		///////////////// Set number of likelihood components per country
		LCPerC = LikeIndexCounter; 

#ifdef PRINT_PROGRAM_PROGRESS
		std::cerr << "DONE " << endl;
#endif
	}
	string CreateDataFileName()
	{
		string DataFileName = ""; 

		if ((SeroSpecificEfficacies) || (SeroSpecific_K_values)) SerotypeString = "_SeroSpec"; //// by definition, if doing case serotypes don't need single serotype data. 
		if (MildAndSevere == TREATED_SEPARATELY)
		{
			DataFileName = SerotypeString + "_MildSevere";
			if (ActiveOrPassivePhase == ACTIVE_PHASE_ONLY)  DataFileName = DataFileName + "_Active";
			if (ModellingHospitalized)	DataFileName = DataFileName + "_UsingHosp";	
		}
		else DataFileName = SerotypeString;

		string FolderString = "Data"; 
		if (!SFU)					FolderString = FolderString + "Backup"; 
		if (Include_Late_Cases)		FolderString = FolderString + "_LateCasesIncluded";
		if (FakeExtObs)				FolderString = FolderString + "_FakeExtendedObservation";
		FolderString = FolderString + "\\";
	
		DataFileName = FolderString + "DataCYD14_15" + DataFileName + ".txt";

		return DataFileName; 
	}
	string CreateOutputString(bool AddOutputStringExtra) // pass CountriesFittedString as an argument because otherwise have to pass Fit_c0, Fit_c1 etc. 
	{
		string OutputString = ""; 

		//// ModelVariant
			 if (ModelVariant == SIMPLE_ANALYTICAL					)	OutputString			= OutputString + "_SIMPLE_ANALYTICAL"	;
		else if (ModelVariant == SIMPLE_NUMERICAL					)	OutputString			= OutputString + "_SIMPLE_NUMERICAL"	;
		else if (ModelVariant == K_SEROPOS							)	OutputString			= OutputString + "_K_SEROPOS"			;
		else if (ModelVariant == DROP_K								)	OutputString			= OutputString + "_DROP_K"				;
		else if (ModelVariant == VAC_SILENT							)	OutputString			= OutputString + "_VAC_SILENT"			;
		else if (ModelVariant == AS_PRIME							)	OutputString			= OutputString + "_AS_PRIME"			;
		else std::cerr << "CreateOutputString: ModelVariant variable incorrectly specified" << endl;

		//// Phase(s)
			 if (PASSIVE_PHASE_ONLY)									OutputString			= OutputString + "_PASSIVE_ONLY"		; 
		else if (ActiveOrPassivePhase == ACTIVE_PHASE_ONLY			)	OutputString			= OutputString + "_ACTIVE"				;
		else if (ActiveOrPassivePhase == DO_ACTIVE_AND_PASSIVE		)	OutputString			= OutputString + "_PASSIVE"				;
		else std::cerr << "CreateOutputString: ActiveOrPassivePhase OR PASSIVE_PHASE_ONLY variable incorrectly specified" << endl; 

			 if (EffNegWane == EffNegWane_Option::FROM_ZERO	)			OutputString = OutputString + "_AWFZ"; //// Alternative waning. efficaies that start out negative don't "wane" to higher values of less absolute magnitude. They "wane from zero" to their negative value. 
		else if (EffNegWane == EffNegWane_Option::NO_WANE	)			OutputString = OutputString + "_nENW"; //// stands for no eff neg waning. See notes in EffNegWane_Option definition. 

		//// disease severity
		if (MildAndSevere == TREATED_SEPARATELY					)	OutputString			= OutputString + "_MILDSEVERE"			;
		if (ModellingHospitalized								)	OutputString			= OutputString + "_hosp"				;		
		if (ModellingHospitalized && !ModelHosp_Indie_Ks		)	OutputString			= OutputString + "SingleMult"			;	

		if (HowManyCaseCategories == 2 && PS_Ks_Multiply_AM_Ks) OutputString = OutputString + "X";

		//// Serotype
			 if (std::regex_match(DataFilename, std::regex("(.*)(_Sero1)(.*)"		)))			OutputString	= OutputString + "_Sero1";
		else if (std::regex_match(DataFilename, std::regex("(.*)(_Sero2)(.*)"		)))			OutputString	= OutputString + "_Sero2";
		else if (std::regex_match(DataFilename, std::regex("(.*)(_Sero3)(.*)"		)))			OutputString	= OutputString + "_Sero3";
		else if (std::regex_match(DataFilename, std::regex("(.*)(_Sero4)(.*)"		)))			OutputString	= OutputString + "_Sero4";
		else if (std::regex_match(DataFilename, std::regex("(.*)(_UnConfSero)(.*)"	)))			OutputString	= OutputString + "_UnConfSero";

			 if (SeroSpecificEfficacies && SeroSpecific_K_values)	OutputString	= OutputString + "_SSVEs"; 
		else if (SeroSpecificEfficacies)							OutputString	= OutputString + "_SS_VEs";
		if (SeroSpecific_K_values)
		{
			if (!SeroSpecificEfficacies) OutputString = OutputString + "_SS_";
			bool AnySSKsFalse = false; 
			if (!SS_KAM_0)  AnySSKsFalse = true;
			if (!SS_KAM_2)  AnySSKsFalse = true;
			if (!SS_KPS_0)  AnySSKsFalse = true;
			if (!SS_KPS_1)  AnySSKsFalse = true;
			if (!SS_KPS_2)  AnySSKsFalse = true;

			if (AnySSKsFalse)
			{
				if (SS_KAM_0)  OutputString = OutputString + "KAM0s";
				if (SS_KAM_2)  OutputString = OutputString + "KAM2s";
				if (SS_KPS_0)  OutputString = OutputString + "KPS0s";
				if (SS_KPS_1)  OutputString = OutputString + "KPS1s";
				if (SS_KPS_2)  OutputString = OutputString + "KPS2s";
			}
			else				OutputString = OutputString + "Ks";
		}
		if (SeroSpecificEfficacies)	
			if (Single_SNeg_Eff && Single_SPos_Eff && (Which_SNeg_SeroFitsAll == Which_SPos_SeroFitsAll))
					OutputString = OutputString + "Equiv";
		if ((SeroSpecificEfficacies || SeroSpecific_K_values) && BaselinePartition)	OutputString	= OutputString + "_BaselinePartition";

		//// dose
		if (SingleOrMultiDose == SINGLE_DOSE)	OutputString	= OutputString + "_SINGLE_DOSE";	// i.e. don't append string if multiple doses (this should be the default)

		//// Countries
		if (FitAllCountries == 0)  //// if fitting all countries don't have a non-empty CountriesFittedString (i.e. fitting all countries is the default. 
		{
			OutputString = OutputString + "_Cs";
			if (Fit_c0)	OutputString	= OutputString + "0"	;
			if (Fit_c1)	OutputString	= OutputString + "1"	;
			if (Fit_c2)	OutputString	= OutputString + "2"	;
			if (Fit_c3)	OutputString	= OutputString + "3"	;
			if (Fit_c4)	OutputString	= OutputString + "4"	;
			if (Fit_c5)	OutputString	= OutputString + "5"	;
			if (Fit_c6)	OutputString	= OutputString + "6"	;
			if (Fit_c7)	OutputString	= OutputString + "7"	;
			if (Fit_c8)	OutputString	= OutputString + "8"	;
			if (Fit_c9)	OutputString	= OutputString + "9"	;
		}

		if (LinKnts)									OutputString	= OutputString + "_LinKnts";
		if (FitWaningRate)							OutputString	= OutputString + "_FLIP_WANING";
		if (Aug_MH_or_Gibbs == MH_AUG)				OutputString	= OutputString + "_MH_AUG";

		if (!AreWeAugmenting)															OutputString	= OutputString + "_NotAug"; 
		if (UseSyntheticData)		if (RandomCases) 									OutputString	= OutputString + "_RandomCases_Seed_" + std::to_string(RandCase_Seed);
		if (SeroSpecific_K_values)	if (RelRisksSameBtwSerotypes)						OutputString	= OutputString + "_SameKsBtwSerotypes_" + std::to_string(Which_Sero_FitsAll);
		if (SingleEff)				OutputString = OutputString + "_SingleEff";

		if (HillWaning) OutputString = OutputString + "_HillWaning";
		if (ResidEffs)  OutputString = OutputString + "_ResidEffs"; //// i.e. residual efficacy at infinity. 

		string WeightString = (PS_Weight > 1) ? std::to_string((int)PS_Weight) : std::to_string((int)((DType)1 / PS_Weight)) + "inv"; ///// this is to stop PS_Weight = 5 or 0.2 (say) being given as 5.00000 or 0.200000. As well as being irritating the decimal point is bound to cause file naming problems. 
		if (Weighting_Pass_Sev) OutputString = OutputString + "_PSweight_" + WeightString;

		if (Fixed_Severe_RelRisks)
		{
				 if (FSKs_Ratio_SetNum == 0)		OutputString = OutputString + "_Fixed_Ks";
			else if (FSKs_Ratio_SetNum == 1)		OutputString = OutputString + "_FSKsV2";
			else if (FSKs_Ratio_SetNum == 2)		OutputString = OutputString + "_FSKs3";
			else std::cerr << "ChooseOutputStringError: FSKs_Ratio_SetNum value not recognized" << endl;
		}

		if (ASVE != Age_Option::INDEPENDENT)
		{
			string ASVE_String = ""; 
				 if (ASVE == Age_Option::HILL			) ASVE_String = "_ASVE"		;	////  "age-specific vaccine efficacy"
			else if (ASVE == Age_Option::CATEGORICAL	) ASVE_String = "_AGSVE"	;	////  "age-group-specific vaccine efficacy"
			else if (ASVE == Age_Option::SPLINE			) ASVE_String = "_ASVESpln"	;	////  "age-specific vaccine efficacy with Spline"
			else if (ASVE == Age_Option::SPLINE_LINE	) ASVE_String = "_ASVELine";	////  "age-specific vaccine efficacy with line"
			else if (ASVE == Age_Option::SPLINE_STEP	) ASVE_String = "_ASVEStep"	;	////  "age-specific vaccine efficacy with Step"
			else if (ASVE == Age_Option::CUBIC			) ASVE_String = "_ASVECubic";	////  "age-specific vaccine efficacy with Step"
			else std::cerr << "WriteOutput Error: ASVE not recognized" << endl; 
		
			if (ASVE == Age_Option::SPLINE || ASVE == Age_Option::SPLINE_LINE || ASVE == Age_Option::SPLINE_STEP || ASVE == Age_Option::CUBIC)
				if (ASVE_AdditionalKnots)	OutputString = OutputString + "AddKnots";

			OutputString = OutputString + ASVE_String;
			if (!AS_VE_Homogeneous && !ASVE_OnlyOneSeroStatus)	OutputString = OutputString + "hetero";
			if (SeroSpecificEfficacies & SSASVE_Additive)		OutputString = OutputString + "Add";
			if (ASVE_OnlyOneSeroStatus && ASVE_BS == SeroNeg)	OutputString = OutputString + "SNeg";
			if (ASVE_OnlyOneSeroStatus && ASVE_BS == SeroPos)	OutputString = OutputString + "SPos" ;
		}
		if (AS_Haz != Age_Option::INDEPENDENT)
		{
				 if (AS_Haz == Age_Option::HILL			)		OutputString = OutputString + "_AS_Haz";
			else if (AS_Haz == Age_Option::CATEGORICAL	)		OutputString = OutputString + "_AS_Hazmult";
			else if (AS_Haz == Age_Option::SPLINE		)		OutputString = OutputString + "_AS_HazSpln";
			else if (AS_Haz == Age_Option::SPLINE_LINE	)		OutputString = OutputString + "_AS_HazLine";
			else if (AS_Haz == Age_Option::SPLINE_STEP	)		OutputString = OutputString + "_AS_HazStep";
			else if (AS_Haz == Age_Option::CUBIC		)		OutputString = OutputString + "_AS_HazCubic";
			else std::cerr << "WriteOutput Error: AS_Haz not recognized" << endl;
		}

		if (AS_Waning != Age_Option::INDEPENDENT)
		{
			OutputString = OutputString + "_AS_Wane";
			if (AS_Waning_OnlyOneSeroStatus)
			{
					 if (AS_Waning_OneSeroBS == SeroNeg) OutputString = OutputString + "SNeg";
				else if (AS_Waning_OneSeroBS == SeroPos) OutputString = OutputString + "SPos";
			}
		}
		if (AS_Waning == Age_Option::CATEGORICAL		)	OutputString = OutputString + "Cat"	;
		if (AS_Waning == Age_Option::SPLINE || AS_Waning == Age_Option::SPLINE_STEP || AS_Waning == Age_Option::SPLINE_LINE || AS_Waning == Age_Option::CUBIC)
		{
			if (AS_Waning == Age_Option::SPLINE			)	OutputString = OutputString + "Spln";
			if (AS_Waning == Age_Option::SPLINE_LINE	)	OutputString = OutputString + "Line";
			if (AS_Waning == Age_Option::SPLINE_STEP	)	OutputString = OutputString + "Step";
			if (AS_Waning == Age_Option::CUBIC			)	OutputString = OutputString + "Cubic";
			if (AS_Waning_KnotSet != 0)
			{
				OutputString = OutputString + "AddKnots"; 
				if (AS_Waning_KnotSet != 1) OutputString = OutputString + std::to_string(AS_Waning_KnotSet);
			}
		}
		if (AS_Waning != Age_Option::INDEPENDENT)
			if (AS_Waning_Homogeneous)
				OutputString = OutputString + "Homog";


		if (PSVEs)			OutputString = OutputString + "_PSVEs";

		if (AdjHaz)
		{
			OutputString = OutputString + "_";
			if (!Skip_BS_BaseHazMult) OutputString = OutputString + "f"; // i.e. fit AdjHaz as opposed to assuming it.

				 if (Which_BS_BaseHazMult == SeroPos) OutputString = OutputString + "AdjHaz"		;
			else if (Which_BS_BaseHazMult == SeroNeg) OutputString = OutputString + "AdjHazSNeg"	;
		}

			 if (PooledCountries	)	OutputString = OutputString + "_PooledCountries"	;
		else if (PooledTrials		)	OutputString = OutputString + "_PooledTrials"		;

		if (Empirical_SeroPrevs	)		OutputString = OutputString + "_Empirical_SeroPrevs";

		if (AllDosesRequired_SNeg && AllDosesRequired_SPos) OutputString = OutputString + "_BothSeroNeed3Dose";
		else
		{
			if (AllDosesRequired_SNeg)						OutputString = OutputString + "_SNegNeed3Dose";
			else
			{
				bool AnySNegNeed3 = false; 
				string AG_SNeg_String = ""; 
				//// if needed, add to AG_SPos_String
				for (int AgeGroup = 1; AgeGroup < HowManyAgeGroups; AgeGroup++) //// !!!!!! **** !!!!! NOTE THAT IT STARTS FROM 1, not zero
					if (AllDosesRequired_AG_BS[AgeGroup][SeroNeg]) { AnySNegNeed3 = true; AG_SNeg_String = AG_SNeg_String + std::to_string(AgeGroup); }

				if (AnySNegNeed3) OutputString = OutputString + "_SNegAG" + AG_SNeg_String + "_Need3Dose"; 
			}
			if (AllDosesRequired_SPos)							OutputString = OutputString + "_SPosNeed3Dose";
			else
			{
				bool AnySPosNeed3 = false;
				string AG_SPos_String = "";
				//// if needed, add to AG_SPos_String
				for (int AgeGroup = 1; AgeGroup < HowManyAgeGroups; AgeGroup++) //// !!!!!! **** !!!!! NOTE THAT IT STARTS FROM 1, not zero
					if (AllDosesRequired_AG_BS[AgeGroup][SeroPos]) { AnySPosNeed3 = true; AG_SPos_String = AG_SPos_String + std::to_string(AgeGroup); }

				if (AnySPosNeed3) OutputString = OutputString + "_SPosAG" + AG_SPos_String + "_Need3Dose";

			}
		}
		if (LTFU_SurvivalCurves	) OutputString = OutputString + "_LTFU_SCurves"; 
	
			 if (ExtImSub == ExtImSub_Option::AS_DATA) OutputString = OutputString + "_ExtImData"; 
		else if (ExtImSub == ExtImSub_Option::AS_PROB) OutputString = OutputString + "_ExtImProb";

		if (ParamRangeFileName != DefaultParamRangeFileName) OutputString = OutputString + "_" + ParamRangeFileName;

		if (SFU)					OutputString = OutputString + "_SFU";
		if (Include_Late_Cases)		OutputString = OutputString + "_IncLate";
		if (FakeExtObs)				OutputString = OutputString + "_FakeExtObs";

		if (AddOutputStringExtra)	OutputString = OutputString + OutputStringExtra;

		//// Guards against horrific misuse of cluster or overwriting of files.
#ifndef USE_CLUSTER
		OutputString = OutputString + "_Dummy";
#endif
#ifndef USE_RANDOM_NUMBERS
		OutputString = OutputString + "_NOT_USING_RANDOM_NUMBERS";
#endif
		return OutputString; 
	}
	string CreateOutputString()
	{
		CreateOutputString(/*AddOutputStringExtra = */ true);
	}

	void init						()
	{
		///// Init will initialize HOUSE structure using all of the above functions. Note that it is called in main AFTER command line parameters are read in - if you do it before you will initialize using defaults only. 
#ifdef PRINT_PROGRAM_PROGRESS
		std::cerr << "InitializeHousekeeping ";
#endif
		Init_Skips					(); //// decide on skipping conditions and initializations for parameters
		init_Strata					(); 
		init_PhaseSevereties		(); 
		InitializeCountriesFitted	(); //// Initialize Countries fitted etc. 
		init_KnotsPolysPerCountry	();	//// For baseline hazard splines
		init_N_STypes				(); //// if SeroSpecificEfficacies || SeroSpecific_K_values
		init_Num_K_Likes_Params		(); 

		///// ///// ///// 
		///// AGESPECIFIC  - later rewrite this to use a class named Spline, which contains coefficients, knots, max spline degree etc., complete with constructors. Right now is messy and irritating. 
		if (ASVE != Age_Option::INDEPENDENT) 
		{
			if (ASVE_OnlyOneSeroStatus)
			{
				AS_VE_Homogeneous		= 0;  //// reset this value from default if VE affected by age for only one serostatus. 
				NumSeroStatuses_ASVEs = 1;
			}
			else
			{
				if (AS_VE_Homogeneous)	NumSeroStatuses_ASVEs = 1;
				else							NumSeroStatuses_ASVEs = HowManySeroStatuses;
			}

			//// Set Number of Knots and MaxSplineDegree
			if (ASVE_AdditionalKnots) KnotsPerSpline_EffMultiplier = 6; else KnotsPerSpline_EffMultiplier = 4;
			if (ASVE == Age_Option::SPLINE_STEP)
			{
				KnotsPerSpline_EffMultiplier--; //// i.e.  don't need the last knot if doing a step function. One less (unneccesary parameter to fit) per serostatus. 
				MaxSplineDegree_EffMultiplier = 0;
			}
			else if (ASVE == Age_Option::SPLINE_LINE	) 	MaxSplineDegree_EffMultiplier = 1;
			else if (ASVE == Age_Option::SPLINE		) 	MaxSplineDegree_EffMultiplier = 2;
			else if (ASVE == Age_Option::CUBIC		) 	MaxSplineDegree_EffMultiplier = 3;

			//// Set Number of Function parameters
				 if (ASVE == Age_Option::HILL)			NumASVE_Function_Params = 3;
			else if (ASVE == Age_Option::CATEGORICAL	)	NumASVE_Function_Params = 4;  //// there are 3 CYD-14 age groups, but  want to use the labels given by AgeGroup1 variable, which starts from 1. Hence one more than 3. Don't fit the first parameter / age group which we'll set to baseline / 1. 
			else 												NumASVE_Function_Params = KnotsPerSpline_EffMultiplier;  ////		//// don't need to list options below - if statement assumes != INDEPENDENT.

			//// Set Number of PolynomialsPerSpline
				 if (ASVE == Age_Option::SPLINE		) PolynomialsPerSpline_EffMultiplier = KnotsPerSpline_EffMultiplier - 2;
			else if (ASVE == Age_Option::SPLINE_LINE	) PolynomialsPerSpline_EffMultiplier = KnotsPerSpline_EffMultiplier - 1;
			else if (ASVE == Age_Option::SPLINE_STEP	) PolynomialsPerSpline_EffMultiplier = KnotsPerSpline_EffMultiplier - 1;
			else if (ASVE == Age_Option::CUBIC		) 
			{
				if (KnotsPerSpline_EffMultiplier == 4) PolynomialsPerSpline_EffMultiplier = 1; 
				else for (int line = 0; line < 100; line++) std::cerr << "ASVE CUBIC and NoKnots don't match - fix" << endl; 
			}
		}

		if (AS_Haz != Age_Option::INDEPENDENT)
		{
			//// Set Number of Knots and MaxSplineDegree
			if (AS_Haz_AdditionalKnots) KnotsPerSpline_HazMultiplier = 6; else KnotsPerSpline_HazMultiplier = 4;
			if (AS_Haz == Age_Option::SPLINE_STEP)
			{
				KnotsPerSpline_HazMultiplier--; //// i.e. you don't need the last knot if doing a step function. One less (unneccesary parameter to fit) per serostatus. 
				MaxSplineDegree_HazMultiplier = 0;
			}
			else if (AS_Haz == Age_Option::SPLINE_LINE) 	MaxSplineDegree_HazMultiplier = 1;
			else if (AS_Haz == Age_Option::SPLINE		)	MaxSplineDegree_HazMultiplier = 2;
			else if (AS_Haz == Age_Option::CUBIC		)	MaxSplineDegree_HazMultiplier = 3;

			//// Set Number of Function parameters
				 if (AS_Haz == Age_Option::HILL		)	NumAS_Haz_Function_Params = 3;
			else if (AS_Haz == Age_Option::CATEGORICAL)	NumAS_Haz_Function_Params = 4;  //// there are 3 CYD-14 age groups, but  want to use the labels given by AgeGroup1 variable, which starts from 1. Hence one more than 3. Don't fit the first parameter / age group which we'll set to baseline / 1. 
			else if (AS_Haz == Age_Option::SPLINE		)	NumAS_Haz_Function_Params = KnotsPerSpline_HazMultiplier;
			else if (AS_Haz == Age_Option::SPLINE_LINE)	NumAS_Haz_Function_Params = KnotsPerSpline_HazMultiplier;
			else if (AS_Haz == Age_Option::SPLINE_STEP)	NumAS_Haz_Function_Params = KnotsPerSpline_HazMultiplier;
			else if (AS_Haz == Age_Option::CUBIC		)	NumAS_Haz_Function_Params = KnotsPerSpline_HazMultiplier;
		
			//// Set Number of PolynomialsPerSpline
				 if (AS_Haz == Age_Option::SPLINE			) PolynomialsPerSpline_HazMultiplier = KnotsPerSpline_HazMultiplier - 2;
			else if (AS_Haz == Age_Option::SPLINE_LINE	) PolynomialsPerSpline_HazMultiplier = KnotsPerSpline_HazMultiplier - 1;
			else if (AS_Haz == Age_Option::SPLINE_STEP	) PolynomialsPerSpline_HazMultiplier = KnotsPerSpline_HazMultiplier - 1;
			else if (AS_Haz == Age_Option::CUBIC)
			{
				if (KnotsPerSpline_HazMultiplier == 4) PolynomialsPerSpline_HazMultiplier = 1;
				else for (int line = 0; line < 100; line++) std::cerr << "AS_Haz CUBIC and NoKnots don't match - fix" << endl;
			}
		}

		//// Set Number of Knots and MaxSplineDegree
	 		 if (AS_Waning_KnotSet == 0)	KnotsPerSpline_WaningDuration = 4;
		else if (AS_Waning_KnotSet == 1)	KnotsPerSpline_WaningDuration = 6;
		else if (AS_Waning_KnotSet == 2)	KnotsPerSpline_WaningDuration = 6;
		else if (AS_Waning_KnotSet == 3)	KnotsPerSpline_WaningDuration = 5;

		if (AS_Waning == Age_Option::SPLINE_STEP)
		{
			KnotsPerSpline_WaningDuration--; //// i.e. you don't need the last knot if doing a step function. One less (unneccesary parameter to fit) per serostatus. 
			MaxSplineDegree_WaningDuration = 0; 
		}
		else if (AS_Waning == Age_Option::SPLINE_LINE	) 	MaxSplineDegree_WaningDuration = 1;
		else if (AS_Waning == Age_Option::SPLINE		) 	MaxSplineDegree_WaningDuration = 2;
		else if (AS_Waning == Age_Option::CUBIC			) 	MaxSplineDegree_WaningDuration = 3;

			 //// Set Number of Function parameters
			 if (AS_Waning == Age_Option::INDEPENDENT)											NumWaningParamsPer_BS = 1;
		else if (AS_Waning == Age_Option::HILL || AS_Waning == Age_Option::CATEGORICAL)		NumWaningParamsPer_BS = 3;  //// Three CYD14 age groups. Unlike AGSVE and AS_Hazmult, you don't intend to use AgeGroup, and therefore don't need them to start from 1. Hence 3 parameters are required, and not 3 + 1 = 4, 
		else /*if (AS_Waning == Age_Option::SPLINE) - // this way all spline versions covered*/	NumWaningParamsPer_BS = KnotsPerSpline_WaningDuration;

		//// Set Number of PolynomialsPerSpline
			 if (AS_Waning == Age_Option::SPLINE		) PolynomialsPerSpline_WaningDuration = KnotsPerSpline_WaningDuration - 2;
		else if (AS_Waning == Age_Option::SPLINE_LINE	) PolynomialsPerSpline_WaningDuration = KnotsPerSpline_WaningDuration - 1; 
		else if (AS_Waning == Age_Option::SPLINE_STEP	) PolynomialsPerSpline_WaningDuration = KnotsPerSpline_WaningDuration - 1; 
		else if (AS_Waning == Age_Option::CUBIC)
		{
			if (KnotsPerSpline_WaningDuration == 4) PolynomialsPerSpline_WaningDuration = 1;
			else for (int line = 0; line < 100; line++) std::cerr << "AS_Waning CUBIC and NoKnots don't match - fix" << endl;
		}

		if (AS_Waning == Age_Option::INDEPENDENT) NumAges_Waning = 1; else NumAges_Waning = HowManyAges;

		if (AS_Waning_Homogeneous)	NumSeroStatuses_AS_Waning = 1;
		else						NumSeroStatuses_AS_Waning = HowManySeroStatuses;

		init_DosesRequired();

		if (DataFilename == "")  //// here if inputting new data set or test data set that is not accounted for in CreateDataFileName function. 
			DataFilename				= CreateDataFileName();		std::cerr << "DataFilename " << DataFilename << endl;
		OutputString					= CreateOutputString(/*AddOutputStringExtra =*/ true);		std::cerr << "OutputString " << OutputString << endl;
		if (OutputStringForOldChainInput		== "")		OutputStringForOldChainInput		= CreateOutputString(/*AddOutputStringExtra =*/ AddOutputStringExtraToOldChains);
		if (OutputStringForOldChainInput_Aug	== "")		OutputStringForOldChainInput_Aug	= CreateOutputString(/*AddOutputStringExtra =*/ AddOutputStringExtraToOldChains);
		if (OutputStringForOldChainInput != OutputString)			
			std::cerr << "Warning: OutputStringForOldChainInput != OutputString: OutputStringForOldChainInput = " << OutputStringForOldChainInput << endl;

		init_L_Indices(); //// this has to come after the other quantities initialized. 

		if (Weighting_Pass_Sev)
			PS_Weight = (DType)1 / PS_Weight; //// Multiplying log likelihood terms equivalent to raising weighted components to a power. Probabilities raised to power > 1 are lower. Therefore to maintain the weighting's intuitive meaning of say "passive components count for say twice as much", I'd like to use a weight of 2, not 1/2. Hence I take reciprocal to make it more intuitive. 

		AdditiveSSASVEs = SeroSpecificEfficacies && ASVE != Age_Option::INDEPENDENT && SSASVE_Additive; ///// this comes up so often - put it in Initialize Housekeeping function and then have only a single bool within functions. 

		Allocate_2D_Array(SSKs_FitMatrix, HowManyCaseCategories, Num_K_Params);
		Populate_2D_Array(SSKs_FitMatrix, true, HowManyCaseCategories, Num_K_Params); //// by default, set all of them to true (Allocate function zeros everything out, i.e. would make false. ). 
		if (SeroSpecific_K_values) ///// i.e. only amend at all if SeroSpecific_K_values
		{
			if (PS_Ks_Multiply_AM_Ks)
			{
				if (!SS_KAM_0) SS_KPS_0 = false;
				if (!SS_KAM_2) SS_KPS_2 = false;
			}

			if (!SS_KAM_0)		SSKs_FitMatrix[ActiveMild][0] = false;
			if (!SS_KAM_1)		SSKs_FitMatrix[ActiveMild][1] = false;
			if (!SS_KAM_2)		SSKs_FitMatrix[ActiveMild][2] = false;
			if (HowManyCaseCategories == 2)
			{
				if (!SS_KPS_0)		SSKs_FitMatrix[PassiveSevere][0] = false;
				if (!SS_KPS_1)		SSKs_FitMatrix[PassiveSevere][1] = false;
				if (!SS_KPS_2)		SSKs_FitMatrix[PassiveSevere][2] = false;
			}
		}

#ifdef PRINT_PROGRAM_PROGRESS
		std::cerr << "InitializeHousekeeping DONE " << endl;
#endif
	}

};
struct ParamNumbers_Struct {

	int HugeParamNo = 40000; //// NULL evaluates to zero for int's. Set default value as higher than No_Parameters will ever be. 

	int MinWaningParam		= HugeParamNo,		MaxWaningParam		= HugeParamNo;		// min and max parameter numbers of waning parameters (durations/rates), plus half lives, powers etc.
	int Min_Knot			= HugeParamNo,		Max_Knot			= HugeParamNo;		// min and max parameter numbers of knots
	int Min_HHaz			= HugeParamNo,		Max_HHaz			= HugeParamNo;		// min and max parameter numbers of historical hazards
	int Min_VacE			= HugeParamNo,		Max_VacE			= HugeParamNo;		// min and max parameter numbers of Vaccine Efficacies
	int Min_VacE_atInf		= HugeParamNo,		Max_VacE_atInf		= HugeParamNo;		// min and max parameter numbers of Vaccine Efficacies at infinity. 
	int Min_qval			= HugeParamNo,		Max_qval			= HugeParamNo;		// min and max parameter numbers of country specific serotype proportions. 
	int Min_K				= HugeParamNo,		Max_K				= HugeParamNo;		// min and max parameter numbers of serotype specific relative risks.  
	int Min_Rho				= HugeParamNo,		Max_Rho				= HugeParamNo;		// min and max parameter numbers of country specific serotype proportions.

	int Min_Hosp_K			= HugeParamNo,		Max_Hosp_K			= HugeParamNo;		//// Hospital K's for CYD15 -  
	int Hosp_K_multiplier	= HugeParamNo;

	int Min_HillHalfLife	= HugeParamNo,		Max_HillHalfLife	= HugeParamNo; 
	int Min_HillPower		= HugeParamNo,		Max_HillPower		= HugeParamNo; 

	int Min_ASVE_Param		= HugeParamNo,		Max_ASVE_Param		= HugeParamNo;
	int Min_ASHaz_Param		= HugeParamNo,		Max_ASHaz_Param		= HugeParamNo;

	int BS_BaseHazMult		= HugeParamNo; //// used if HOUSE.AdjHaz == true

	int Min_AS_PrimingParam = HugeParamNo,		Max_AS_PrimingParam = HugeParamNo; 

	bool IsEqualToAnother(const ParamNumbers_Struct &AnotherStruct)
	{
		return
			HugeParamNo		==		AnotherStruct.HugeParamNo		&&
			Min_Knot		==		AnotherStruct.Min_Knot 			&&
			Min_HHaz		==		AnotherStruct.Min_HHaz 			&&
			Min_VacE		==		AnotherStruct.Min_VacE 			&&
			Min_qval		==		AnotherStruct.Min_qval			&&
			Min_K			==		AnotherStruct.Min_K				&&
			Max_Knot		==		AnotherStruct.Max_Knot			&&
			Max_HHaz		==		AnotherStruct.Max_HHaz			&&
			Max_VacE		==		AnotherStruct.Max_VacE			&&
			Max_qval		==		AnotherStruct.Max_qval			&&
			Max_K			==		AnotherStruct.Max_K				&&
			Min_Rho			==		AnotherStruct.Min_Rho			&&
			Max_Rho			==		AnotherStruct.Max_Rho			;
	}
	void SetEqualToAnother(const ParamNumbers_Struct &AnotherStruct)
	{
		HugeParamNo		=		AnotherStruct.HugeParamNo		;
		Min_Knot		=		AnotherStruct.Min_Knot 			;
		Min_HHaz		=		AnotherStruct.Min_HHaz 			;
		Min_VacE		=		AnotherStruct.Min_VacE 			;
		Min_qval		=		AnotherStruct.Min_qval			;
		Min_K			=		AnotherStruct.Min_K				;
		Max_Knot		=		AnotherStruct.Max_Knot			;
		Max_HHaz		=		AnotherStruct.Max_HHaz			;
		Max_VacE		=		AnotherStruct.Max_VacE			;
		Max_qval		=		AnotherStruct.Max_qval			;
		Max_K			=		AnotherStruct.Max_K				;
		Min_Rho			=		AnotherStruct.Min_Rho			;
		Max_Rho			=		AnotherStruct.Max_Rho			;
	}
};
struct Params_Struct {

	DType KA_0 = 0.8, KA_1 = 1, KA_2 = 0.2, KH_0 = 0.06, KH_1 = 0.17, KH_2 = 0.04;  //// used for data input and initial parameter values only. 
	DType ** Initial_Ks = NULL;				//// indexed by i) PhaseSeverity and ii) No. Previous Infections. Used for parameter input eventually. 
	DType Hosp_K_mult = 0.5;				//// multiply all hospital K's for CYD 15 countries by this number. 
	DType PosEfficacy = 0.5, NegEfficacy = 0.27, PosWaning = 12, NegWaning = 0.5;
	//// these are just here so you can easily read in parameters. 
	DType Initial_Halflife_Hill = 2, Initial_Power_Hill = 2;		
	DType Initial_Halflife_ASVE = 2, Initial_Power_ASVE = 2;		
	DType Initial_Prop_ASVE = 0;									//// Proportion of Efficacy that is dependent on age. Zero value reduces to default model where efficacy unaffected by age. Make this the default. 
	DType Initial_Halflife_AS_Waning = 8, Initial_Power_AS_Waning = 10;		
	DType Initial_Halflife_AS_Haz = 2, Initial_Power_AS_Haz = 2;	
	DType Initial_Prop_AS_Haz = 0;									//// Proportion of hazard that is dependent on age. Zero value reduces to default model where hazard unaffected by age. Make this the default. 
	DType Initial_BS_BaseHazMult = 0.7;
	DType Initial_AS_Prime_PropIndptAge_SNeg	= 0, Initial_AS_Prime_PropIndptAge_SPos = 0; //// i.e. completely dependent on age. 
	DType Initial_AS_Prime_SNeg_rate			= 5, Initial_AS_Prime_SPos_rate = 5;

	DType * Hill_Halflives = NULL, *Hill_Powers = NULL;
	DType * Fixed_SevereK_ratios = NULL;				//// values set according to Science paper - Table S3 Q0 and Q2=Q3. 

	std::vector<DType>	ParamVec; ///// vector of parameter values. 
	std::vector<string> NamesOfParameters;
	std::vector<int> ParamNosYouWillFit;

	ParamNumbers_Struct ParamNumbers;		//// Depending on ModelVariant, ActivePhase and loads of other things the index of say KHS_1 will change, but not during runtime. This data structure stores the number of the parameter name

	//// **** //// Splines etc. //// **** //// 
	//// To model baseline hazard / force of infection
	DType ** xKnots = NULL, ** yKnots = NULL;	//// indexed by i) country ii) knot
	DType *** SplineCoeffs = NULL;				//// indexed by i) country ii) polynomial iii) coefficient
	//// For various age effects (multiplier of hazard / vaccine efficacy / vaccine duration or rate), as alternative to categorical variables and hill function). 
	DType ** Age_Effs_xKnots			= NULL;	//// indexed by i) BaselineSeroStatus ii) knot 
	DType * Age_HazMult_xKnots			= NULL;	//// indexed by i) knot 
	DType ** Age_DurationRate_xKnots	= NULL;	//// indexed by i) BaselineSeroStatus ii) knot  
	DType *** Age_SplineCoeffs_Effs			= NULL;					//// indexed by i) BaselineSeroStatus ii) polynomial iii) coefficient
	DType ** Age_SplineCoeffs_HazMult		= NULL;					//// indexed by							i) polynomial ii) coefficient
	DType *** Age_SplineCoeffs_DurationRate	= NULL;					//// indexed by i) BaselineSeroStatus ii) polynomial iii) coefficient


	int LikeSize;
	DType ** LikeParts = NULL; //// Indexed by i) country ii) component
	DType LikeFull = NULL, Total_logLike = NULL;
	DType LogPosterior	= NULL	; 
	DType * ParamArray_logPriors = NULL, LogPriorFull = 0; ; ////

	DType ** w_LikeParts	= NULL	, w_LikeFull	= NULL; //// Indexed by i) country ii) component. w_LikeParts is same as LikeParts but passive components are weighted. w_LikeFull is full weighted likelihood. 
	DType ** ptr_LikeParts	= NULL	, *ptr_L_Full	= NULL; //// Indexed by i) country ii) component. ptr_LikeParts & ptr_L_Full are pointers to the likelihood array and full likelihood. . Used when weighting particular parts of data: You may wish to point to weighted likelihood. Value set in Initialize_LikeStruture function. 

	DType **** K_s				= NULL;		//// indexed by i) country, ii) Phase/Severity iii) No. Previous Infections, and iv) serotype
	DType *** KplusValues		= NULL;		//// indexed by i) country, ii) age, iii) PhaseSeverity ///// ////// IMPORTANT IMPORTANT IMPORTANT IMPORTANT IMPORTANT IMPORTANT: If using serospecific K's, Plus values will refer to aggregate over all serotypes, and hence used only for survivors, not cases. 
	DType **** Meta_KplusValues = NULL;		//// indexed by i) Trial Arm ii) country, iii) age, iv) PhaseSeverity ///// Used for SSKs and SSVEs together. For K_SEROPOS variant, vaccine arm KplusValues differ from control arm KplusValues as efficacies must be included.     ////// IMPORTANT IMPORTANT IMPORTANT IMPORTANT IMPORTANT IMPORTANT: If using serospecific K's, Plus values will refer to aggregate over all serotypes, and hence used only for survivors, not cases. 

	DType ***** Ks_Prime		= NULL;		//// For HOUSE.ModelVariant == AS_PRIME indexed by i) country, ii) PhaseSeverity iii) BaselineSeroStatus iv) serotype and v) age. Need country to deal with hospital K's. 
	DType *** KPlusPrimeValues	= NULL;		//// For HOUSE.ModelVariant == AS_PRIME indexed by i) country, ii) age, iii) PhaseSeverity
	DType *** SumRhoK0_SNeg_Primes	= NULL;		//// (for serospec Ks && AS_PRIME) indexed by i) country, ii) PhaseSeverity, iii) age. Don't need BaselineSeroStatus as you only ever do this for seronegatives. 

	DType *** Efficacies	= NULL;		//// indexed by i) PhaseSeverity; ii) serotype,	iii) BaselineSeroStatus;  If doing HOUSE.ResidEffs, where efficacies wane to a non zero value at infinity (i.e. not default), these refer to intitial efficacies, not Inf_Effs, which refer to efficacies as t -> inf. 
	DType *** Inf_Effs		= NULL;		//// indexed by i) PhaseSeverity; ii) serotype,	iii) BaselineSeroStatus;  Inf_Effs, which refer to efficacies as t -> inf. 

	DType *** VacHazLikes	= NULL;		//// indexed by i) country ii) BaselineSeroStatus iii) PhaseSeverity. For given trial phase , this value this value same but needed for all serotypes (if doing Serospecific efficacies), makes sense to calculate once then record then multiply accordingly. 
	
	DType **IVH_vals		= NULL;		//// indexed by i) person; ii) PhaseSeverity;Person-specific stored array of integrated vaccine hazards. Affected by change in knots, waning parameters or augmentation. 
	
	DType ** rhos				= NULL;		//// (for serospec efficacies & Ks) indexed by i) country,	ii) serotype. proportions of each serotype in each country.
	DType *** SumRhoEffs		= NULL;		//// indexed by i) PhaseSeverity ii) country	iii) BaselineSeroStatus  This is for the augmentation. Whenever you calculate survival probability, you need to multiply integrated vaccine hazard by sum_d^HOUSE.N_STypes(rho_c_d * Efficacy_d) for each serostatus. This quantity will be the same for all individuals of the same country and serostatus so seems pointless to calculate repeatedly. 
	DType *** SumRhoKs			= NULL;		//// (for serospec Ks) indexed by i) country, ii) PhaseSeverity, iii) No. prior infections (K2 taken to mean 2 or 3 prior infections). Do this just for serospecific K's  
	DType * SumRhos				= NULL;		//// (for serospec efficacies & Ks, Baseline_Partition version) indexed by country. 
	DType **** SumRhoEffNegs	= NULL;		//// Used when HOUSE.AltWaneEffNeg. indexed by i) PhaseSeverity; ii) country; iii) Age; iv) BaselineSerostatus. Value of SumRhoEffNegs[PS][c][A][BS] is given by sum_d (\rho_{cd}^*(A)), where \rho_{cd}*(A) = rho_{cd} (if Eff >= 0) OR rho_cd*(1 - Eff) if Eff < 0; 

	DType *** SumRhoEffKs	= NULL;		//// (for SSKs & SSVEs together) indexed by i) country, ii) PhaseSeverity, iii) BaselineSeroStatus

	DType ** AgeEff_Mult	= NULL;		//// (For age-dependent vaccine efficacies). Indexed by i) BaselineSerostatus; ii) age (although if HOUSE.AS_VE_Homogeneous then first index will point to same memory regardless of BaselinSeroStatus - i.e. will effictively by 1-dimensional)
	DType ** ASVE_Params	= NULL;		//// (For age-dependent vaccine efficacies). Indexed by i) BaselineSerostatus; ii) Parameter; (Again, if HOUSE.AS_VE_Homogeneous then first index will point to same memory regardless of BaselinSeroStatus - i.e. will effictively by 1-dimensional)

	DType * AgeHaz_Mult		= NULL;		//// (for AS_Hazards / age specific hazards). Indexed by i) age. 
	DType * ASHaz_Params	= NULL;		//// (For age-dependent vaccine efficacies). Indexed by i) Parameter; (Again, if HOUSE.AS_VE_Homogeneous then first index will point to same memory regardless of BaselinSeroStatus - i.e. will effictively by 1-dimensional)
	DType * BS_BaseHazMults	= NULL;		//// multipliers of baseline hazard by baseline serostatus. 
	DType **** SeroPrevs	= NULL;		//// indexed by i) log-or-not-log; ii) country; iii) serostatus; iv) age.  Will usually be governed by historical hazards, otherwise by empirical seroprevalence, but will always be required. 

	DType ** AS_Priming_Params	= NULL;	//// For HOUSE.ModelVariant == AS_PRIME. Indexed by i) BaselineSerostatus; ii) Parameter. 

	DType ** ParamRanges		= NULL; //// indexed by i) LowerBound or UpperBound	; ii) param number
	DType *  ProposalStandDevs	= NULL; //// indexed by i) param number

	//// numerical 
	DType **WaningParams			= NULL;		//// indexed by i) BaselineSeroStatus	; ii) ParamType (if (!HOUSE.AS_Waning) then this is simply the duration or the rate (depending on flip waning or not). if (HOUSE.AS_Waning), then this is a) Maximum rate/duration, b) Halflife, c) Power
	DType ***WaningMults			= NULL;		//// indexed by i) Age in years			; ii) BaselineSeroStatus iii) daypostdose. Indexed by age for where we want to have the duration (or rate) of vaccine efficacy depend on age (i.e. if HOUSE.AS_Waning). If !HOUSE.AS_Waning, then set all "age" pointers to point to the same thing. 

	DType **IntBaseHazLookUp		= NULL;		//// indexed by i) country ii) Calendar Time (in days)
	DType **BaselineHazardValues	= NULL;		//// indexed by i) country ii) Calendar Time (in days)

	//// analytical
	std::vector<DType> *** RootsArray = NULL;
	DType **AreasUnderFullPolynomials = NULL;

	bool IsEqualToAnother (const Params_Struct &AnotherStruct, const Housekeeping_Struct &HOUSE, const DATA_struct &DATA)
	{
		//DType **IntBaseHazLookUp = NULL, **BaselineHazardValues = NULL;

		bool BaseHazardsEqual		= 1; 
		bool Int_BaseHazardsEqual	= 1;

		if (BaselineHazardValues != 0 && IntBaseHazLookUp != 0)
			for (int country = 0; country < HOUSE.TotalCountries; country++)
				for (int CalendarDay = 0; CalendarDay <= DATA.NumCalendarDaysFollowUp[country]; CalendarDay++)
				{
					if (BaselineHazardValues[country][CalendarDay] != AnotherStruct.BaselineHazardValues[country][CalendarDay]) BaseHazardsEqual		= 0; 
					if (IntBaseHazLookUp	[country][CalendarDay] != AnotherStruct.IntBaseHazLookUp	[country][CalendarDay]) Int_BaseHazardsEqual	= 0;
				}

		bool ParamNosYouWillFit_Equal = true; 
		if (ParamNosYouWillFit.size() != AnotherStruct.ParamNosYouWillFit.size()) ParamNosYouWillFit_Equal = false; //// check sizes equal
		else 
			for (int i = 0; i < ParamNosYouWillFit.size(); i++) //// check elements equal. 
				if (ParamNosYouWillFit[i] != AnotherStruct.ParamNosYouWillFit[i])  ParamNosYouWillFit_Equal = false;

		bool All_Equal = BaseHazardsEqual && Int_BaseHazardsEqual && ParamNosYouWillFit_Equal &&
			KA_0				== AnotherStruct.KA_0																						&&
			KA_1				== AnotherStruct.KA_1																						&&
			KA_2				== AnotherStruct.KA_2																						&&
			KH_0				== AnotherStruct.KH_0																						&&
			KH_1				== AnotherStruct.KH_1																						&&
			KH_2				== AnotherStruct.KH_2																						&&
			PosEfficacy			== AnotherStruct.PosEfficacy																				&&
			NegEfficacy			== AnotherStruct.NegEfficacy																				&&
			PosWaning			== AnotherStruct.PosWaning																					&&
			NegWaning			== AnotherStruct.NegWaning																					&&
			ParamVec			== AnotherStruct.ParamVec																					&&
			NamesOfParameters	== AnotherStruct.NamesOfParameters		&&

			LikeSize			== AnotherStruct.LikeSize				&&
			LikeFull			== AnotherStruct.LikeFull				&&
			w_LikeFull			== AnotherStruct.w_LikeFull				&&
			ptr_L_Full			== AnotherStruct.ptr_L_Full				&&

			ParamNumbers.IsEqualToAnother(AnotherStruct.ParamNumbers) &&
			IsEqual_2D_Arrays(xKnots, AnotherStruct.xKnots, HOUSE.TotalCountries, HOUSE.KnotsPerCountry) &&
			IsEqual_2D_Arrays(yKnots, AnotherStruct.yKnots, HOUSE.TotalCountries, HOUSE.KnotsPerCountry) &&
			
			IsEqual_3D_Arrays(SplineCoeffs, AnotherStruct.SplineCoeffs, HOUSE.TotalCountries, HOUSE.PolynomialsPerCountry, HOUSE.MaxSplineDegree + 1) &&

			IsEqual_2D_Arrays(LikeParts, AnotherStruct.LikeParts, HOUSE.TotalCountries, HOUSE.LCPerC) &&
			IsEqual_2D_Arrays(w_LikeParts, AnotherStruct.w_LikeParts, HOUSE.TotalCountries, HOUSE.LCPerC) &&
			IsEqual_2D_Arrays(ptr_LikeParts, AnotherStruct.ptr_LikeParts, HOUSE.TotalCountries, HOUSE.LCPerC) &&

			IsEqual_1D_Arrays(SumRhos  , AnotherStruct.SumRhos  , HOUSE.TotalCountries) &&

			IsEqual_3D_Arrays(KplusValues, AnotherStruct.KplusValues, HOUSE.TotalCountries, HOUSE.HowManyAges, HOUSE.HowManyCaseCategories) &&

			IsEqual_3D_Arrays(VacHazLikes		, AnotherStruct.VacHazLikes		, HOUSE.TotalCountries	, HOUSE.HowManySeroStatuses		, HOUSE.HowManyCaseCategories) &&
			IsEqual_2D_Arrays(rhos				, AnotherStruct.rhos			, HOUSE.TotalCountries	, HOUSE.N_STypes				) &&

			IsEqual_3D_Arrays(Efficacies		, AnotherStruct.Efficacies		, HOUSE.HowManyCaseCategories	, HOUSE.N_STypes_VEs			, HOUSE.HowManySeroStatuses	) &&
			IsEqual_3D_Arrays(SumRhoEffs		, AnotherStruct.SumRhoEffs		, HOUSE.HowManyCaseCategories	, HOUSE.TotalCountries			, HOUSE.HowManySeroStatuses	) &&
			IsEqual_3D_Arrays(SumRhoKs			, AnotherStruct.SumRhoKs		, HOUSE.TotalCountries			, HOUSE.HowManyCaseCategories	, HOUSE.Num_K_Params		) &&
			IsEqual_3D_Arrays(SumRhoEffKs		, AnotherStruct.SumRhoEffKs		, HOUSE.TotalCountries			, HOUSE.HowManyCaseCategories	, HOUSE.HowManySeroStatuses	) &&

			IsEqual_4D_Arrays(K_s				, AnotherStruct.K_s				, HOUSE.TotalCountries	, HOUSE.HowManyCaseCategories	, HOUSE.Num_K_Params, HOUSE.N_STypes_Ks) &&
			//IsEqual_2D_Arrays(WaningValues		, AnotherStruct.WaningValues, 2, DATA.NoDaysOfFollowUp + 1) &&
			IsEqual_3D_Arrays(WaningMults		, AnotherStruct.WaningMults, HOUSE.NumAges_Waning, HOUSE.HowManySeroStatuses, DATA.NoDaysOfFollowUp + 1) &&
			IsEqual_4D_Arrays(Meta_KplusValues	, AnotherStruct.Meta_KplusValues, HOUSE.HowManyTrialArms, HOUSE.TotalCountries, HOUSE.HowManyAges, HOUSE.HowManyCaseCategories) &&
			IsEqual_2D_Arrays(ASVE_Params		, AnotherStruct.ASVE_Params, HOUSE.HowManySeroStatuses, HOUSE.NumASVE_Function_Params) &&
			IsEqual_2D_Arrays(AgeEff_Mult		, AnotherStruct.AgeEff_Mult, HOUSE.HowManySeroStatuses, HOUSE.HowManyAges) ;

		return All_Equal; 

	} 
	void SetEqualToAnother(const Params_Struct &AnotherStruct, const Housekeeping_Struct &HOUSE, const DATA_struct &DATA)
	{
		for (int country = 0; country < HOUSE.TotalCountries; country++)
			for (int CalendarDay = 0; CalendarDay <= DATA.NumCalendarDaysFollowUp[country]; CalendarDay++)
			{
				BaselineHazardValues[country][CalendarDay] = AnotherStruct.BaselineHazardValues	[country][CalendarDay]	;
				IntBaseHazLookUp	[country][CalendarDay] = AnotherStruct.IntBaseHazLookUp		[country][CalendarDay]	;
			}
		ParamNosYouWillFit	= AnotherStruct.ParamNosYouWillFit; 

		KA_0				= AnotherStruct.KA_0					;
		KA_1				= AnotherStruct.KA_1					;
		KA_2				= AnotherStruct.KA_2					;
		KH_0				= AnotherStruct.KH_0					;
		KH_1				= AnotherStruct.KH_1					;
		KH_2				= AnotherStruct.KH_2					;
		PosEfficacy			= AnotherStruct.PosEfficacy				;
		NegEfficacy			= AnotherStruct.NegEfficacy				;
		PosWaning			= AnotherStruct.PosWaning				;
		NegWaning			= AnotherStruct.NegWaning				;
		ParamVec			= AnotherStruct.ParamVec				;
		NamesOfParameters	= AnotherStruct.NamesOfParameters		;
																	;
		LikeSize			= AnotherStruct.LikeSize				;
		LikeFull			= AnotherStruct.LikeFull				;

		ParamNumbers.SetEqualToAnother(AnotherStruct.ParamNumbers) 										;
		SetEqual_2D_Arrays(xKnots, AnotherStruct.xKnots, HOUSE.TotalCountries, HOUSE.KnotsPerCountry)	; 
		SetEqual_2D_Arrays(yKnots, AnotherStruct.yKnots, HOUSE.TotalCountries, HOUSE.KnotsPerCountry)	; 


		SetEqual_3D_Arrays(SplineCoeffs, AnotherStruct.SplineCoeffs, HOUSE.TotalCountries, HOUSE.PolynomialsPerCountry, HOUSE.MaxSplineDegree + 1);

		SetEqual_2D_Arrays(LikeParts, AnotherStruct.LikeParts, HOUSE.TotalCountries, HOUSE.LCPerC);

		SetEqual_1D_Arrays(SumRhos, AnotherStruct.SumRhos, HOUSE.TotalCountries);

		SetEqual_4D_Arrays(K_s		, AnotherStruct.K_s		, HOUSE.TotalCountries, HOUSE.HowManyCaseCategories, HOUSE.Num_K_Params, HOUSE.N_STypes_Ks);
		SetEqual_3D_Arrays(KplusValues, AnotherStruct.KplusValues, HOUSE.TotalCountries, HOUSE.HowManyAges, HOUSE.HowManyCaseCategories);
		
		//SetEqual_2D_Arrays(S_typeVEs		, AnotherStruct.S_typeVEs		, HOUSE.N_STypes_VEs		, HOUSE.HowManySeroStatuses		) ;
		SetEqual_3D_Arrays(VacHazLikes		, AnotherStruct.VacHazLikes		, HOUSE.TotalCountries	, HOUSE.HowManySeroStatuses		, HOUSE.HowManyCaseCategories) ;
		SetEqual_2D_Arrays(rhos				, AnotherStruct.rhos			, HOUSE.TotalCountries	, HOUSE.N_STypes				) ;
		//SetEqual_2D_Arrays(SumRhoEffProds	, AnotherStruct.SumRhoEffProds	, HOUSE.TotalCountries	, HOUSE.HowManySeroStatuses		) ;
		
		SetEqual_3D_Arrays(Efficacies		, AnotherStruct.Efficacies		, HOUSE.HowManyCaseCategories	, HOUSE.N_STypes_VEs			, HOUSE.HowManySeroStatuses	) ;
		SetEqual_3D_Arrays(SumRhoKs			, AnotherStruct.SumRhoKs		, HOUSE.TotalCountries			, HOUSE.HowManyCaseCategories	, HOUSE.Num_K_Params		) ;
		SetEqual_3D_Arrays(SumRhoEffs		, AnotherStruct.SumRhoEffs		, HOUSE.HowManyCaseCategories	, HOUSE.TotalCountries			, HOUSE.HowManySeroStatuses	) ;
		//SetEqual_2D_Arrays(WaningValues		, AnotherStruct.WaningValues, 2, DATA.NoDaysOfFollowUp + 1);
		SetEqual_3D_Arrays(WaningMults		, AnotherStruct.WaningMults, HOUSE.NumAges_Waning, HOUSE.HowManySeroStatuses, DATA.NoDaysOfFollowUp + 1);
	}
	void WhichPartsUnequal(const Params_Struct &AnotherStruct, const Housekeeping_Struct &HOUSE, const DATA_struct &DATA)
	{
		bool BaseHazardsEqual		= 1;
		bool Int_BaseHazardsEqual	= 1;

		if (BaselineHazardValues != 0 && IntBaseHazLookUp != 0)
			for (int country = 0; country < HOUSE.TotalCountries; country++)
				if (DATA.NumCalendarDaysFollowUp != 0)
					for (int CalendarDay = 0; CalendarDay <= DATA.NumCalendarDaysFollowUp[country]; CalendarDay++)
					{
						if (BaselineHazardValues[country][CalendarDay] != AnotherStruct.BaselineHazardValues[country][CalendarDay]) BaseHazardsEqual		= 0;
						if (IntBaseHazLookUp	[country][CalendarDay] != AnotherStruct.IntBaseHazLookUp	[country][CalendarDay]) Int_BaseHazardsEqual	= 0;
					}

		bool ParamNosYouWillFit_Equal = true;
		if (ParamNosYouWillFit.size() != AnotherStruct.ParamNosYouWillFit.size()) ParamNosYouWillFit_Equal = false; //// check sizes equal
		else
			for (int i = 0; i < ParamNosYouWillFit.size(); i++) //// check elements equal. 
				if (ParamNosYouWillFit[i] != AnotherStruct.ParamNosYouWillFit[i])  ParamNosYouWillFit_Equal = false;


		if (	!ParamNosYouWillFit_Equal) std::cout << "ParamNosYouWillFit_Equal not equal" << endl;
		if (	!BaseHazardsEqual	) std::cout << "BaselineHazardValues not equal" << endl;
		if (	!Int_BaseHazardsEqual) std::cout << "IntBaseHazLookUp not equal" << endl;

		if (	KA_0				!= AnotherStruct.KA_0						)	std::cout << "KA_0 not equal" << endl;
		if (	KA_1				!= AnotherStruct.KA_1						)	std::cout << "KA_1 not equal" << endl;
		if (	KA_2				!= AnotherStruct.KA_2						)	std::cout << "KA_2 not equal" << endl;
		if (	KH_0				!= AnotherStruct.KH_0						)	std::cout << "KH_0 not equal" << endl;
		if (	KH_1				!= AnotherStruct.KH_1						)	std::cout << "KH_1 not equal" << endl;
		if (	KH_2				!= AnotherStruct.KH_2						)	std::cout << "KH_2 not equal" << endl;
		if (	PosEfficacy			!= AnotherStruct.PosEfficacy				)	std::cout << "PosEfficacy not equal" << endl;
		if (	NegEfficacy			!= AnotherStruct.NegEfficacy				)	std::cout << "NegEfficacy not equal" << endl;
		if (	PosWaning			!= AnotherStruct.PosWaning					)	std::cout << "PosWaning	 not equal" << endl;
		if (	NegWaning			!= AnotherStruct.NegWaning					)	std::cout << "NegWaning	 not equal" << endl;
		if (	ParamVec			!= AnotherStruct.ParamVec					)	std::cout << "ParamVec	 not equal" << endl;
		if (	NamesOfParameters	!= AnotherStruct.NamesOfParameters			)	std::cout << "NamesOfParameters not equal" << endl;

		if (	LikeSize			!= AnotherStruct.LikeSize					)	std::cout << "LikeSize not equal" << endl;
		if (	LikeFull			!= AnotherStruct.LikeFull					)	std::cout << "LikeFull not equal" << endl;
		if (	w_LikeFull			!= AnotherStruct.w_LikeFull					)	std::cout << "w_LikeFull not equal" << endl;

		if (	!ParamNumbers.IsEqualToAnother(AnotherStruct.ParamNumbers)		)	std::cout << "ParamNumbers not equal" << endl;

		if (	!IsEqual_1D_Arrays(SumRhos			, AnotherStruct.SumRhos			, HOUSE.TotalCountries								))	std::cout << "SumRhos not equal" << endl;

		if (	!IsEqual_2D_Arrays(LikeParts		, AnotherStruct.LikeParts		, HOUSE.TotalCountries	, HOUSE.LCPerC				))	std::cout << "LikeParts not equal" << endl;

		if (	!IsEqual_2D_Arrays(w_LikeParts		, AnotherStruct.w_LikeParts		, HOUSE.TotalCountries	, HOUSE.LCPerC))	std::cout << "w_LikeParts	not equal" << endl;
		if (	!IsEqual_2D_Arrays(ptr_LikeParts	, AnotherStruct.ptr_LikeParts	, HOUSE.TotalCountries	, HOUSE.LCPerC))	std::cout << "ptr_LikeParts not equal" << endl;

		if (	!IsEqual_3D_Arrays(SumRhoEffKs		, AnotherStruct.SumRhoEffKs		, HOUSE.TotalCountries			, HOUSE.HowManyCaseCategories	, HOUSE.HowManySeroStatuses)) std::cout << "SumRhoEffKs not equal" << endl;
		if (	!IsEqual_3D_Arrays(SumRhoEffs		, AnotherStruct.SumRhoEffs		, HOUSE.HowManyCaseCategories	, HOUSE.TotalCountries			, HOUSE.HowManySeroStatuses)) std::cout << "SumRhoEffs not equal" << endl;
		
		if (	!IsEqual_2D_Arrays(xKnots			, AnotherStruct.xKnots			, HOUSE.TotalCountries	, HOUSE.KnotsPerCountry		))	std::cout << "xKnots not equal" << endl;
		if (	!IsEqual_2D_Arrays(yKnots			, AnotherStruct.yKnots			, HOUSE.TotalCountries	, HOUSE.KnotsPerCountry		))	std::cout << "yKnots not equal" << endl;
		if (	!IsEqual_2D_Arrays(rhos				, AnotherStruct.rhos			, HOUSE.TotalCountries	, HOUSE.N_STypes			))	std::cout << "rhos	 not equal" << endl;
		if (	!IsEqual_3D_Arrays(WaningMults		, AnotherStruct.WaningMults		, HOUSE.NumAges_Waning	, HOUSE.HowManySeroStatuses, DATA.NoDaysOfFollowUp + 1	))	std::cout << "WaningMults	 not equal" << endl;

		if (	!IsEqual_3D_Arrays(Efficacies		, AnotherStruct.Efficacies		, HOUSE.HowManyCaseCategories	, HOUSE.N_STypes_VEs			, HOUSE.HowManySeroStatuses		))	std::cout << "Efficacies not equal" << endl;
		if (	!IsEqual_3D_Arrays(SplineCoeffs		, AnotherStruct.SplineCoeffs	, HOUSE.TotalCountries			, HOUSE.PolynomialsPerCountry	, HOUSE.MaxSplineDegree + 1		))	std::cout << "SplineCoeffs not equal" << endl;
		
		if (	!IsEqual_4D_Arrays(K_s				, AnotherStruct.K_s				, HOUSE.TotalCountries			, HOUSE.HowManyCaseCategories	, HOUSE.Num_K_Params, HOUSE.N_STypes_Ks))	std::cout << "RelRisks not equal" << endl;
		if (	!IsEqual_3D_Arrays(KplusValues		, AnotherStruct.KplusValues		, HOUSE.TotalCountries			, HOUSE.HowManyAges				, HOUSE.HowManyCaseCategories	))	std::cout << "KplusValues not equal" << endl;
		if (	!IsEqual_3D_Arrays(VacHazLikes		, AnotherStruct.VacHazLikes		, HOUSE.TotalCountries			, HOUSE.HowManySeroStatuses		, HOUSE.HowManyCaseCategories	))	std::cout << "VacHazLikes not equal" << endl;
		if (	!IsEqual_3D_Arrays(SumRhoKs			, AnotherStruct.SumRhoKs		, HOUSE.TotalCountries			, HOUSE.HowManyCaseCategories	, HOUSE.Num_K_Params			))	std::cout << "SumRhoKs	 not equal" << endl;

		if (	!IsEqual_4D_Arrays(Meta_KplusValues	, AnotherStruct.Meta_KplusValues, HOUSE.HowManyTrialArms	, HOUSE.TotalCountries, HOUSE.HowManyAges, HOUSE.HowManyCaseCategories) )	std::cout << "Meta_KplusValues	not equal" << endl;
		if (	!IsEqual_2D_Arrays(ASVE_Params		, AnotherStruct.ASVE_Params		, HOUSE.HowManySeroStatuses	, HOUSE.NumASVE_Function_Params) 										)	std::cout << "ASVE_Params		not equal" << endl;
		if (	!IsEqual_2D_Arrays(AgeEff_Mult		, AnotherStruct.AgeEff_Mult		, HOUSE.HowManySeroStatuses	, HOUSE.HowManyAges)													)	std::cout << "AgeEff_Mult		not equal" << endl;
	}
};
struct SeroPrev_Struct {

	DType ***AgeSpecificSeroPrevs	= NULL;		// indexed by i) Country; ii) index (augmented patients, non-augmented patients, or all patients) iii) age // (NumC + 3) because want CYD14, CYD15, and CYD14 and CYD15 combined. 
	DType ***AgeSpecificPopSizes	= NULL;		// indexed by i) Country; ii) index (augmented patients, non-augmented patients, or all patients) iii) age // (NumC + 3) because want CYD14, CYD15, and CYD14 and CYD15 combined. 

	DType ***Aug_Patients_MetaAgeSpecificSeroPrevs = NULL;
	DType ***All_Patients_MetaAgeSpecificSeroPrevs = NULL;

	DType **Aug_Patients_Mean_AgeSpecificSeroPrevs			= NULL;		/* indexed by: i) country; ii) age */
	DType **Aug_Patients_Median_AgeSpecificSeroPrevs		= NULL;		/* indexed by: i) country; ii) age */
	DType **Aug_Patients_LowerCI_AgeSpecificSeroPrevs		= NULL;		/* indexed by: i) country; ii) age */
	DType **Aug_Patients_UpperCI_AgeSpecificSeroPrevs		= NULL;		/* indexed by: i) country; ii) age */
															  			   								   
	DType **All_Patients_Mean_AgeSpecificSeroPrevs			= NULL;		/* indexed by: i) country; ii) age */
	DType **All_Patients_Median_AgeSpecificSeroPrevs		= NULL;		/* indexed by: i) country; ii) age */
	DType **All_Patients_LowerCI_AgeSpecificSeroPrevs		= NULL;		/* indexed by: i) country; ii) age */
	DType **All_Patients_UpperCI_AgeSpecificSeroPrevs		= NULL;		/* indexed by: i) country; ii) age */

	void AllocateMemory						(const Housekeeping_Struct &HOUSE, int NumPosteriorSamples)
	{
#ifdef PRINT_PROGRAM_PROGRESS
		std::cerr << "SeroPrev_Struct::AllocateMemory: NumPosteriorSamples " << NumPosteriorSamples;
#endif
		Allocate_3D_Array(Aug_Patients_MetaAgeSpecificSeroPrevs, NumPosteriorSamples, HOUSE.TotalCountries + 3, HOUSE.HowManyAges);
		Allocate_3D_Array(All_Patients_MetaAgeSpecificSeroPrevs, NumPosteriorSamples, HOUSE.TotalCountries + 3, HOUSE.HowManyAges);

		Allocate_3D_Array(AgeSpecificSeroPrevs, HOUSE.TotalCountries + 3, 3, HOUSE.HowManyAges);
		Allocate_3D_Array(AgeSpecificPopSizes , HOUSE.TotalCountries + 3, 3, HOUSE.HowManyAges);

		Allocate_2D_Array(Aug_Patients_Mean_AgeSpecificSeroPrevs		, HOUSE.TotalCountries + 3, HOUSE.HowManyAges);
		Allocate_2D_Array(Aug_Patients_Median_AgeSpecificSeroPrevs		, HOUSE.TotalCountries + 3, HOUSE.HowManyAges);
		Allocate_2D_Array(Aug_Patients_LowerCI_AgeSpecificSeroPrevs		, HOUSE.TotalCountries + 3, HOUSE.HowManyAges);
		Allocate_2D_Array(Aug_Patients_UpperCI_AgeSpecificSeroPrevs		, HOUSE.TotalCountries + 3, HOUSE.HowManyAges);

		Allocate_2D_Array(All_Patients_Mean_AgeSpecificSeroPrevs		, HOUSE.TotalCountries + 3, HOUSE.HowManyAges);
		Allocate_2D_Array(All_Patients_Median_AgeSpecificSeroPrevs		, HOUSE.TotalCountries + 3, HOUSE.HowManyAges);
		Allocate_2D_Array(All_Patients_LowerCI_AgeSpecificSeroPrevs		, HOUSE.TotalCountries + 3, HOUSE.HowManyAges);
		Allocate_2D_Array(All_Patients_UpperCI_AgeSpecificSeroPrevs		, HOUSE.TotalCountries + 3, HOUSE.HowManyAges);
#ifdef PRINT_PROGRAM_PROGRESS
		std::cerr << " DONE ";
#endif
	}
	void Calc_PopSizes_And_NonAugSeroPrev	(const DATA_struct &DATA, const Housekeeping_Struct &HOUSE)
	{
		int Country, AgeInYears;
		for (int patient = 0; patient < DATA.NonAugmentedIndices.size(); patient++)
		{
			AgeInYears = (int)floor(DATA.ai_s[DATA.NonAugmentedIndices[patient]]);
			AgeInYears = int(floor(AgeInYears));
		
			Country = DATA.ci_s[DATA.NonAugmentedIndices[patient]];

			if (DATA.Ii_s[DATA.NonAugmentedIndices[patient]] == 1) //// if seropos at baseline, add to AgeSpecificSeroPrevs
			{
				// add to appropriate country-seroprevalence categories. 
				/*add to country*/	AgeSpecificSeroPrevs[Country	][NonAugIndex][AgeInYears]++;		

				if (Country <= 4)	AgeSpecificSeroPrevs[NumC		][NonAugIndex][AgeInYears]++;		// CYD14, need NumC+1'th term (C++ indexing). 
				else				AgeSpecificSeroPrevs[NumC + 1	][NonAugIndex][AgeInYears]++;		// CYD15, need NumC+2'th term (C++ indexing). 

				/*add to CYD1415*/	AgeSpecificSeroPrevs[NumC + 2	][NonAugIndex][AgeInYears]++;			// CYD14 and 15 together, need NumC+3'th term (C++ indexing). 
			}
			// add to population size. 
			AgeSpecificPopSizes[Country	][NonAugIndex][AgeInYears]++;
			AgeSpecificPopSizes[NumC + 2	][NonAugIndex][AgeInYears]++;								// CYD14 and 15 together, need NumC+3'th term (C++ indexing).
		
			if (Country <= 4)	AgeSpecificPopSizes[NumC		][NonAugIndex][AgeInYears]++;			// CYD14, need NumC+1'th term (C++ indexing). 
			else				AgeSpecificPopSizes[NumC + 1	][NonAugIndex][AgeInYears]++;			// CYD15, need NumC+2'th term (C++ indexing). 
		}

		//// divide SeroPrevs by PopSize denominator if non-zero. 
		for (int country = 0; country < (NumC + 3); country++)
			for (int age = 0; age < HOUSE.HowManyAges; age++)
				if (AgeSpecificPopSizes[country][NonAugIndex][age] > 0) 
					AgeSpecificSeroPrevs[country][NonAugIndex][age] = AgeSpecificSeroPrevs[country][NonAugIndex][age] / AgeSpecificPopSizes[country][NonAugIndex][age];

		// calculate age-specific population sizes (but not seroprevalences) for augmented patients. 
		for (int patient = 0; patient < DATA.AugmentedIndices.size(); patient++)
		{
			AgeInYears = (int)floor(DATA.ai_s[DATA.AugmentedIndices[patient]]);
			AgeInYears = int	(floor(AgeInYears));

			Country		= DATA.ci_s[DATA.AugmentedIndices[patient]];

			// add to population size. 
			AgeSpecificPopSizes[Country	][AugIndex][AgeInYears]++;
			AgeSpecificPopSizes[NumC + 2	][AugIndex][AgeInYears]++;							// CYD14 and 15 together, need NumC+3'th term (NumC++ indexing).
	
			if (Country <= 4)	AgeSpecificPopSizes[NumC		][AugIndex][AgeInYears]++;		// CYD14, need NumC+1'th term (NumC++ indexing). 
			else				AgeSpecificPopSizes[NumC + 1	][AugIndex][AgeInYears]++;		// CYD15, need NumC+2'th term (NumC++ indexing). 
		}
	}

	void init							(const DATA_struct &DATA, const Housekeeping_Struct &HOUSE, int NumPosteriorSamples, bool AreWeCalculatingSeroPrevs)
	{
		if (AreWeCalculatingSeroPrevs)
		{
#ifdef PRINT_PROGRAM_PROGRESS
			std::cerr << endl << " SEROPREV init: NumPosteriorSamples " << NumPosteriorSamples << endl;
#endif
			AllocateMemory(HOUSE, NumPosteriorSamples);
			Calc_PopSizes_And_NonAugSeroPrev(DATA, HOUSE);
		}
	}
	void Calc_Aug_AgeSpecificSeroPrev	(const DATA_struct &DATA, const Housekeeping_Struct &HOUSE, int PostSampleNum)
	{
#ifdef PRINT_PROGRAM_PROGRESS
		std::cerr << "START SeroPrev, " ; 
#endif
		//// clear array from previous iterations. 
		for (int country = 0; country < (NumC + 3); country++)
		{
			for (int age = 0; age < HOUSE.HowManyAges; age++)
			{
				AgeSpecificSeroPrevs[country][AugIndex			][age] = 0;
				AgeSpecificSeroPrevs[country][AugAndNonAugIndex][age] = 0;
			}
		}

		// only recalculate age-specific seroprevalence (not population sizes) in augmented data as non-augmented prevs constant. 
		int Country, AgeInYears;

		for (int patient = 0; patient < DATA.AugmentedIndices.size(); patient++)
		{
			AgeInYears = (int)floor(DATA.ai_s[DATA.AugmentedIndices[patient]]);
			AgeInYears = int(floor(AgeInYears));

			Country = DATA.ci_s[DATA.AugmentedIndices[patient]];
			if (DATA.Ii_s[DATA.AugmentedIndices[patient]] == 1)
			{
				// add to appropriate country-seroprevalence category. 
				AgeSpecificSeroPrevs[Country	][AugIndex][AgeInYears]++;
				AgeSpecificSeroPrevs[NumC + 2	][AugIndex][AgeInYears]++;							// CYD14 and 15 together, need NumC+3'th term (C++ indexing). 
				if (Country <= 4)	AgeSpecificSeroPrevs[NumC		][AugIndex][AgeInYears]++;		// CYD14, need NumC+1'th term (C++ indexing). 
				else				AgeSpecificSeroPrevs[NumC + 1	][AugIndex][AgeInYears]++;		// CYD15, need NumC+2'th term (C++ indexing). 
			}
		}
		//// divide SeroPrevs by PopSize denominator if non-zero. 
		DType AugPrev, NonAugPrev, AugSize, NonAugSize; 
		for (int country = 0; country < (NumC + 3); country++)
			for (int age = 0; age < HOUSE.HowManyAges; age++)
				if (AgeSpecificPopSizes[country][AugIndex][age] > 0)
				{
					AugSize		= AgeSpecificPopSizes[country][AugIndex	][age];
					NonAugSize	= AgeSpecificPopSizes[country][NonAugIndex	][age];

					AugPrev		= AgeSpecificSeroPrevs	[country][AugIndex		][age] / AugSize;
					NonAugPrev	= AgeSpecificSeroPrevs	[country][NonAugIndex	][age];	//// don't need to divide by NonAugSize as done previously and has not changed. 

					AgeSpecificSeroPrevs	[country][AugIndex][age] = AugPrev; // store value. 
					//// average of augmented and non-augmented seroprevalences, weighted by population size. 
					AgeSpecificSeroPrevs	[country][AugAndNonAugIndex][age] = ((AugSize * AugPrev) + (NonAugSize * NonAugPrev) )		/	(AugSize + NonAugSize);

					Aug_Patients_MetaAgeSpecificSeroPrevs[PostSampleNum][country][age] = AgeSpecificSeroPrevs[country][AugIndex			][age]; // store in meta
					All_Patients_MetaAgeSpecificSeroPrevs[PostSampleNum][country][age] = AgeSpecificSeroPrevs[country][AugAndNonAugIndex	][age]; // store in meta
				}
#ifdef PRINT_PROGRAM_PROGRESS
		std::cerr << "END SeroPrev" << endl;
#endif
	}
	void Calc_Final_SeroPrev_Output		(const Housekeeping_Struct &HOUSE, int NumPosteriorSamples, int NumElementsOutside_CrI_Tails)
	{
#ifdef PRINT_PROGRAM_PROGRESS
		std::cerr << "Calc_Final_SeroPrev_Output: NumPosteriorSamples " << NumPosteriorSamples << " NumElementsOutside_CrI_Tails " << NumElementsOutside_CrI_Tails << endl;
#endif
		int NoElementsForMedian = NumPosteriorSamples / 2;  
		int LowerIndex = max(0, int(NumElementsOutside_CrI_Tails));
		int UpperIndex = min(int(NumPosteriorSamples - 1), int(NumPosteriorSamples - NumElementsOutside_CrI_Tails));

		if (LowerIndex >= NumPosteriorSamples) LowerIndex = NumPosteriorSamples - 1;
		if (UpperIndex >= NumPosteriorSamples) UpperIndex = NumPosteriorSamples - 1;

		if (NumPosteriorSamples > 0)
		{
			std::vector<DType> AgeSpecificSeroPrevsOverLoadsOfIterations(NumPosteriorSamples, 0);
			std::vector<DType> All_AgeSpecificSeroPrevsOverLoadsOfIterations(NumPosteriorSamples, 0);

			/// calculate mean/median and 95% credible itervals for each day post dose and country hazard group, so far as it stands.
			for (int country = 0; country < (HOUSE.TotalCountries + 3); country++)
				for (int age = 0; age < HOUSE.HowManyAges; age++)
				{
					DType MeanSeroPrev		= 0;
					DType MeanSeroPrev_All	= 0;
					for (int iter = 0; iter < NumPosteriorSamples; iter++)
					{
						AgeSpecificSeroPrevsOverLoadsOfIterations		[iter] = Aug_Patients_MetaAgeSpecificSeroPrevs[iter][country][age]; // add/amend SurvivalProbsOverLoadsOfIterations vector. 
						All_AgeSpecificSeroPrevsOverLoadsOfIterations	[iter] = All_Patients_MetaAgeSpecificSeroPrevs[iter][country][age]; // add/amend SurvivalProbsOverLoadsOfIterations vector. 
						MeanSeroPrev		+= Aug_Patients_MetaAgeSpecificSeroPrevs[iter][country][age];
						MeanSeroPrev_All	+= All_Patients_MetaAgeSpecificSeroPrevs[iter][country][age];
					}
					// sort vector 
					sort(AgeSpecificSeroPrevsOverLoadsOfIterations.begin()		, AgeSpecificSeroPrevsOverLoadsOfIterations.end());
					sort(All_AgeSpecificSeroPrevsOverLoadsOfIterations.begin()	, All_AgeSpecificSeroPrevsOverLoadsOfIterations.end());

					MeanSeroPrev		/= NumPosteriorSamples;
					MeanSeroPrev_All	/= NumPosteriorSamples;

					Aug_Patients_Mean_AgeSpecificSeroPrevs		[country][age] = MeanSeroPrev;
					Aug_Patients_Median_AgeSpecificSeroPrevs	[country][age] = AgeSpecificSeroPrevsOverLoadsOfIterations		[NoElementsForMedian];
					Aug_Patients_LowerCI_AgeSpecificSeroPrevs	[country][age] = AgeSpecificSeroPrevsOverLoadsOfIterations		[LowerIndex];
					Aug_Patients_UpperCI_AgeSpecificSeroPrevs	[country][age] = AgeSpecificSeroPrevsOverLoadsOfIterations		[UpperIndex];

					All_Patients_Mean_AgeSpecificSeroPrevs		[country][age] = MeanSeroPrev_All;
					All_Patients_Median_AgeSpecificSeroPrevs	[country][age] = All_AgeSpecificSeroPrevsOverLoadsOfIterations	[NoElementsForMedian];
					All_Patients_LowerCI_AgeSpecificSeroPrevs	[country][age] = All_AgeSpecificSeroPrevsOverLoadsOfIterations	[LowerIndex];
					All_Patients_UpperCI_AgeSpecificSeroPrevs	[country][age] = All_AgeSpecificSeroPrevsOverLoadsOfIterations	[UpperIndex];
				}
		}
	}
};
struct Chains_Struct {

	SeroPrev_Struct SEROPREV; 
	int * Final_Iis = new int[NPat](); //// used for storing actual Iis whn you need to write output. Bit hacky but keep anyway. 
	DType ** FinalLikeParts = NULL; 
	std::vector<DType> FinalParamVec; 
	DType Final_L_Full = 0; 

	DType ProportionUpdatesAccepted	= 0; //// proportion of updates  accepted
	DType AverageAcceptanceProb		= 0; //// average Acceptance probability (across all parameters)
	DType Proposed_Param, ProposedVsOldlogPosteriorRatio, AcceptOrReject_RandomnNumber, logAcceptanceProb, Temp_ThisIteration = 1, LogPrior, LogPosterior;
	std::vector<int> SimulataneousUpdateParams;
	std::vector<DType> Proposed_Param_Vector; //// Used if you need to update parameters simultaneously, but with different values (e.g. Fixed_Severe_RelRisks)
	bool LikeAndLikeAlike = 1;
	DType ** LikePartsForChecking = NULL, *** VacHazLike_CheckingArray = NULL;

	int AddtoChainEveryHowManyIterations = 100, AugmentEveryHowManyIterations = 5, NumPosteriorSamples = 0, MaxNumPosteriorSamples = 10000, NumElementsOutside_CrI_Tails = 0, HowFarDidYouGet = 0;
	int No_Iterations = NULL, BurnIn = NULL, WhenToFlipTemperatureBack = 0, thread_num = 0;
	string OutputFolder = "Output\\"; 

	// (Summary) states of chain. Even though these refer to parameter values, best to keep them with chains to avoid confusion with CurrentPARAMS and ProposedPARAMS
	std::vector<DType>	MeanPostParamVec;				//// mean posterior param vec. This is the mean of all sampled parameter values. 
	std::vector<DType>	MaxLike_ParamVec;				//// max likelihood parameters. 
	std::vector<DType>	ModalPost_ParamVec;				//// Modal Posterior parameters. 
	DType MaxLikeSoFar		= -DBL_MAX;					//// used for max likelihood. 
	DType ModalPostSoFar	= -DBL_MAX;					//// used for ModalPosterior
	int *MaxLike_Ii_s		= new int	[NPat]();		//// maximum likelihood state of augmentation chain. 
	int *ModalPost_Ii_s		= new int	[NPat]();		//// maximum likelihood state of augmentation chain. 
	DType *MeanPost_Ii_s	= new DType	[NPat]();		//// maximum likelihood state of augmentation chain. 
	DType *** rhos = NULL;		//// (for serospec efficacies & Ks) indexed by i) country,	ii) serotype, iii) iteration. proportions of each serotype in each country.
	DType alpha, *logLikeChain = NULL, **ParamChain = NULL, *AcceptArray = NULL, *AcceptArray_Aug = NULL, *log_PosteriorChain = NULL, *Total_Likelihood_Chain = NULL;
	DType *AugAccArray_SingleIter = NULL; 
	DType ** OldParamChain = NULL;
	DType *** LL_SPrev_Chains = NULL;					//// indexed by i) Posterior Sample; ii) Country - inc. (pooled) trial; iii) Serostatus (SNeg or SPos)
	DType *** LL_SPrev_Chains_ImSub = NULL;				//// indexed by i) Posterior Sample; ii) Country - inc. (pooled) trial; iii) Serostatus (SNeg or SPos)

	DType UnacceptableLogLikeDifference		= 50;
	const bool BreakifLikeDifferencesLarge	= 0;
	bool LikeDifferencesAcceptable			= 1, ShouldYouStop = 0; 
	double randomNo = 0;
	bool ParamAcceptCondition = 1;

	bool ImportOldChain_Params = 0			,		ImportOldChain_AugmentedData = 0;
	string OldChains_ParamFileName = ""		,		OldChains_AugmentedDataFileName = "";

	bool OutputParticularChainState = 0; 
	string ParticularChainState_ParamFileName, ParticularChainState_AugDataFileName;
	int IterationToOutput = NULL; 

	bool AreWeFittingParameters				= 1;
	bool AreWeAugmenting					= 1;
	bool AreWeChecking						= 0;
	bool CalculateFreshLikeLihood			= 0;
	bool CheckIndividualAugmentation		= 0;
	bool AreWeCalculatingSeroPrevs			= 1;
	bool SimulatedAnnealing					= 0;

	//// ==== //// ==== //// ==== //// ==== 
	//// ==== Survival Curve / Attack rate / Posterior hazard ratio logic Gates. 
	
	bool Calc_SCsARsHRPs_Any 				= 1;	//// do you calculate ANY Survival Curves / Attack Rates / Posterior hazard ratios, (either mean/credible intervals/maxlike/modalpost/)? If so this shoudl be true. Essentialy this will allocate memory for Curve Tails (which is wasteful if you don't need them), but will also calculate Rownames and millions of other associated quantities, so tolerate the waste.  
	bool Calc_SCsARsHRPs_MeanAndCrIs		= 1;	 
	bool Calc_SCsARsHRPs_ModalMaxLike		= 1;	//// if you want only the modal posterior (will be most useful after WBICs). 

	bool HazRatiosOnly						= false; 
	bool OutputIndividualSurvivalTables		= 0;
	
	bool Output_LL_minus_Aug = 0;
	DType * logLike_minus_Aug_Chain = NULL;


	bool CoolDuringBurnIn				= false;
	bool CoolAfterBurnIn				= false;
	bool ExponentialSchedule			= false;  /// if false then linear
	bool UseDefault_FinalTemperature	= true; /// will set FinalTemperature to 1/log(HOUSE.TotalPatients) in initialize chains functions if true. HOUSE.TotalPatients will vary if not fitting both trials together. 
	DType FinalTemperature				= 1/log(NPat);  ///// This is 1/log(31126) final value of T that multiplies ratio of posteriors/likelihoods. 

	void init(const Housekeeping_Struct &HOUSE, const DATA_struct &DATA, const Params_Struct &PARAMS)
	{
		if (!HOUSE.AreWeAugmenting	)		AreWeAugmenting = false; //// HOUSE variable is a catch all for both CHAINS and WBIC_CHAINS (or any other chains structure you may need). 
		if (AreWeChecking			)		std::cerr << "CHECKING MCMC THE SLOW WAY " << endl;
		if (AreWeChecking			)		std::cerr << "CHECKING MCMC THE SLOW WAY " << endl;
		if (AreWeChecking			)		std::cerr << "CHECKING MCMC THE SLOW WAY " << endl;
		if (AreWeChecking			)		std::cerr << "CHECKING MCMC THE SLOW WAY " << endl;
		if (AreWeChecking			)		std::cerr << "CHECKING MCMC THE SLOW WAY " << endl;
		if (CalculateFreshLikeLihood)		std::cerr << "CALCULATING MCMC THE SLOW WAY " << endl;
		if (CalculateFreshLikeLihood)		std::cerr << "CALCULATING MCMC THE SLOW WAY " << endl;
		if (CalculateFreshLikeLihood)		std::cerr << "CALCULATING MCMC THE SLOW WAY " << endl;
		if (CalculateFreshLikeLihood)		std::cerr << "CALCULATING MCMC THE SLOW WAY " << endl;
		Allocate_2D_Array(FinalLikeParts, HOUSE.TotalCountries, HOUSE.LCPerC); 

		if (Calc_SCsARsHRPs_Any == false) //// if so, then overwrite the other parts. 
		{
			Calc_SCsARsHRPs_MeanAndCrIs		= false;
			Calc_SCsARsHRPs_ModalMaxLike	= false;
		}

#ifdef PRINT_PROGRAM_PROGRESS
		std::cerr << " InitializeChains ";
#endif
		if (WhenToFlipTemperatureBack == 0) WhenToFlipTemperatureBack = No_Iterations + 1; //// i.e. if still set to default, then default is don't flip temperature back during runtime. 

		alpha = 0.05;

		NumPosteriorSamples = (No_Iterations - BurnIn) / AddtoChainEveryHowManyIterations;

		if (NumPosteriorSamples > MaxNumPosteriorSamples)
		{
			AddtoChainEveryHowManyIterations *= (NumPosteriorSamples / MaxNumPosteriorSamples); //// i.e. if you have more samples than you need, must increase no. iterations before you add to chain) add to chain less often. 

			//// redefine after resetting AddtoChainEveryHowManyIterations
			NumPosteriorSamples = (No_Iterations - BurnIn) / AddtoChainEveryHowManyIterations;
		}

		NumElementsOutside_CrI_Tails = NumPosteriorSamples * alpha * 0.5 + 1; // not quite true, but you want everything WITHIN the middle 95% CrI. If you don't add an extra element you will get the values immediately outsde the credible interval. 
		if (NumElementsOutside_CrI_Tails < 1) NumElementsOutside_CrI_Tails = 1; /// when you run proper analysis, 

		Allocate_2D_Array(ParamChain, HOUSE.No_Parameters, NumPosteriorSamples);
		logLikeChain			= new DType[NumPosteriorSamples	]();
		log_PosteriorChain		= new DType[NumPosteriorSamples	]();
		Total_Likelihood_Chain	= new DType[NumPosteriorSamples	]();
		AcceptArray				= new DType[HOUSE.No_Parameters	]();	
		AcceptArray_Aug			= new DType[HOUSE.max_threads	]();	
		AugAccArray_SingleIter	= new DType[HOUSE.max_threads	]();	
		if (HOUSE.OutputSeroPrev_LLChains)
		{
			Allocate_3D_Array(LL_SPrev_Chains		, NumPosteriorSamples, HOUSE.TotalCountries + 3, HOUSE.HowManySeroStatuses);
			Allocate_3D_Array(LL_SPrev_Chains_ImSub	, NumPosteriorSamples, HOUSE.TotalCountries + 3, HOUSE.HowManySeroStatuses);
		}

		MeanPostParamVec	 = PARAMS.ParamVec; 
		MaxLike_ParamVec	 = PARAMS.ParamVec; 
		ModalPost_ParamVec	 = PARAMS.ParamVec; 

		if (HOUSE.SeroSpecificEfficacies || HOUSE.SeroSpecific_K_values)
			if (!HOUSE.BaselinePartition)
				Allocate_3D_Array(rhos, HOUSE.TotalCountries, HOUSE.N_STypes, NumPosteriorSamples);

		if (Output_LL_minus_Aug)
			logLike_minus_Aug_Chain = new DType[NumPosteriorSamples]();

		if (UseDefault_FinalTemperature)
		{
#ifdef PRINT_PROGRAM_PROGRESS
			std::cerr << "Init Chains using default temp" << std::endl;
#endif
			FinalTemperature = (DType)1 / log(HOUSE.TotalPatients);
		}

		SEROPREV.init(DATA, HOUSE, NumPosteriorSamples, AreWeCalculatingSeroPrevs); 

#ifdef PRINT_PROGRAM_PROGRESS
		std::cerr << "NumElementsOutside_CrI_Tails " << NumElementsOutside_CrI_Tails << " NumPosteriorSamples " << NumPosteriorSamples << " AddtoChainEveryHowManyIterations " << AddtoChainEveryHowManyIterations << endl;
		std::cerr << " InitializeChains DONE" << endl;
#endif
	}
};
class SCurves {
public:

	DType **Strata_Sizes				= NULL;		//// indexed by				i)  strata; ii)  day post dose. The number of subjects in a given strata will decrease as people drop out of the study, although not many drop outs. 
	DType *** Threaded_Strata_Sizes		= NULL;		//// indexed by i) thread;	ii) strata; iii) day post dose. The number of subjects in a given strata will decrease as people drop out of the study, although not many drop outs. 

	DType *** PostSample_SurvivalTables					= NULL;		//// indexed by				i) disease (either, mild or severe);  ii) strata; iii) day post dose
	DType **** Threaded_PostSample_SurvivalTables		= NULL;		//// indexed by i) thread; ii) disease (either, mild or severe); iii) strata; iv)  day post dose
	DType **** FinalPosteriorSurvivalCurves				= NULL;		//// indexed by i) Disease Severity; ii) Statistic (i.e. Mean, LowerCrI, or UpperCrI); iii) stratum; iv) day post dose 

	DType **** MaxMin_Vals_ForTails				= NULL;		//// indexed by i) Disease Severity; ii) Index (UpperTail_Index or LowerTail_Index); iii) stratum; iv) day post dose
	int **** MaxMin_Indices_ForTails			= NULL;		//// indexed by i) Disease Severity; ii) Index (UpperTail_Index or LowerTail_Index); iii) stratum; iv) day post dose
	DType ***** CurveTails						= NULL;		//// indexed by i) Disease Severity; ii) Index (UpperTail_Index or LowerTail_Index); iii) stratum; iv) day post dose; v) index of vector //// want to bundle lower and upper tails (as well as Max/Min index etc). //// want to index these using upper or lower index, Survival Curve row (SubjectCategory), Survival Curve column (Days post dose), then tail index, i.e. the elements of the tail. 

	SCurves(int NumStrata, int NumDays, int NumDiseaseCategories, int max_threads, int NoElementsInPostDistTails)
	{
		init(NumStrata, NumDays, NumDiseaseCategories, max_threads, NoElementsInPostDistTails);
	}
	SCurves() {}; //// this default constructor does nothing. But will allow you to create array of pointers to type SCurves and allocate memory properly later. 
	void init(int NumStrata, int NumDays, int NumDiseaseCategories, int max_threads, int NoElementsInPostDistTails)
	{
		Allocate_2D_Array(Strata_Sizes													, NumStrata, NumDays);
		Allocate_3D_Array(Threaded_Strata_Sizes	, max_threads							, NumStrata, NumDays);

		//// Final survival tables (mean, lower CrI and upper CrIs)
		Allocate_4D_Array(FinalPosteriorSurvivalCurves, NumDiseaseCategories, 3	, NumStrata, NumDays);

		//// Survival table (for individual parameter sets). 
		Allocate_3D_Array(PostSample_SurvivalTables							, NumDiseaseCategories, NumStrata, NumDays);
		Allocate_4D_Array(Threaded_PostSample_SurvivalTables, max_threads	, NumDiseaseCategories, NumStrata, NumDays);

		//// for tails/95% credible intervals. Lower tails have 1 as initial value. UpperTails have 0 as initial value. Initial index where min/max value occurs is arbitrary, so setting to zero (default in allocate function) okay. 
		Allocate_4D_Array(MaxMin_Indices_ForTails	, NumDiseaseCategories, 2, NumStrata, NumDays);
		Allocate_4D_Array(MaxMin_Vals_ForTails		, NumDiseaseCategories, 2, NumStrata, NumDays);

		//// initialize the min and max values of the lower and upper tails.
		//DType InitialTailValues[2] = { 1,0 }; //// every probability will be lower than 1, and higher than 0, hence these are initial values for Lower and Upper tails. 
		//// every probability will be lower than 1, and higher than 0, hence these are initial values for Lower and Upper tails. 
		//// Update: as you'll be using the same "machinery" for Hazard ratios, these will not be bounded by 1. Set upper limit to something ludicrous like 500. Same logic should apply (probabilities less than 500 after all). Only difference may be if you run a couple of samples only, which you never do. 
		DType InitialTailValues[2] = { DBL_MAX,-DBL_MAX }; 
		for (int DiseaseSeverity = 0; DiseaseSeverity < NumDiseaseCategories; DiseaseSeverity++)
			for (int category = 0; category < NumStrata; category++)
				for (int daypostdose = 0; daypostdose < NumDays; daypostdose++)
					for (int WhichTail = 0; WhichTail < 2; WhichTail++)
						MaxMin_Vals_ForTails[DiseaseSeverity][WhichTail][category][daypostdose] = InitialTailValues[WhichTail];

		/// Initialize the tails vectors. 
		Allocate_5D_Array(CurveTails, NumDiseaseCategories, 2, NumStrata, NumDays, NoElementsInPostDistTails);
		for (int DiseaseSeverity = 0; DiseaseSeverity < NumDiseaseCategories; DiseaseSeverity++)
			for (int WhichTail = 0; WhichTail < 2; WhichTail++)
				for (int category = 0; category < NumStrata; category++)
					for (int daypostdose = 0; daypostdose < NumDays; daypostdose++)
						for (int ElementIndex = 0; ElementIndex < NoElementsInPostDistTails; ElementIndex++)
							CurveTails[DiseaseSeverity][WhichTail][category][daypostdose][ElementIndex] = MaxMin_Vals_ForTails[DiseaseSeverity][WhichTail][category][daypostdose];
	}
	bool is_initialised()
	{
		return (Strata_Sizes != NULL) & (Threaded_Strata_Sizes != NULL) & (PostSample_SurvivalTables != NULL) & (Threaded_PostSample_SurvivalTables != NULL) & (FinalPosteriorSurvivalCurves != NULL) &
			(MaxMin_Vals_ForTails != NULL) & (MaxMin_Indices_ForTails != NULL) & (CurveTails != NULL);
	}
};
struct Survival_Struct {  //// version with mild and severe added. 

	int **Strata_SCurves = NULL;		//// indexed by i) Trial Arm; ii) BaselineSerostatus. Matrix of various strata using old labelling. 

	int HowManyDiseaseSeverities	= 1; //// Either 1 if not doing mild and severe, or 3 if doing mild and severe (need either, mild or severe). Note this is diffirent to HOUSE.HowManyCaseCategories as that is a catch all term for modelling either both the active and passive phase, or modelling mild and severe disease, which in any case wouldn't contain "either". 
	int HowManyTrialPhases			= 2; //// Different than HOUSE.HowManyTrialPhases as need an option to include ActiveAndPassivePhase. 
	
	std::vector<string> DiseaseNames;
	std::vector<string> StatisticNames; /// Mean, LowerCrI, and UpperCrI
	std::vector<string> TrialPhaseNames; /// _ActivePhaseOnly, Passive_PhaseOnly, and Both (will be simply ""). 

	int PosteriorSampleCounter = 0; //// was used for debugging and to ensure that the process of thinning out posterior samples for survival curves worked. 
	int MaxNumSurvivalPosteriorSamples = 1000, NoSurvivePostSamples = 0, AddToSurvivalCurvesEveryHowManyIterations = 0, NumElementsOutside_CrI_Tails = 0;

	int HazGroupsPerCountry = 4, ExtraHazGroupsPerCountry = 5, TotalGroupsPerCountry = HazGroupsPerCountry + ExtraHazGroupsPerCountry;				// ExtraHazGroupsPerCountry = 5;  because in addition to regular hazard groups you want to combine for sets I, Ic, V, Vc and ALL.
	int TotalCountries		= NumC + 3; //// includes (pooled) trials. 
	int NoSubjectCategoriesPerAgeGroup = TotalGroupsPerCountry * TotalCountries;	// (NumC + 3) because in addition to all countries, want CYD14, CYD15, and CYD14 and CYD15 combined. 
	int NumAgeGroupCategories = 9;
	
	/*
	0		= any age group
	1,2,3	= CYD14 trial 2 to 5 years, 6 to 11 years, 12 to 14 years
	4,5		= CYD15 trial 9 to 11 years, 12 to 16 years.
	6,7		= CYD14 trial under 9 years, over 9 years
	8		= CYD14 and CYD15 trial over 9 years (no need for CYD15 under 9 years as all CYD15 participants are under 9 years, so this is just the same as CYD15 any age group).
	*/

	int NoDaysOfFollowUp, NoActiveDaysFollowUp, NoPassiveDaysFollowUp;

	int NoSubjectCategories = NumAgeGroupCategories * NoSubjectCategoriesPerAgeGroup; /// this assumes that all age groups have same number of countries - false but should just result in loads of zeros (and combined will be identical to individual trial) so though a little wasteful, should be easier to code. 
	std::vector<string> SurvivalTableRowNames, SurvivalTableColNames, AttackRateColNames, PassivePhaseSurvivalTableColNames;

	SCurves SC_WT;  //// survival curves - whole trial
	SCurves SC_PP;  //// survival curves - passive phase

	///// ///// ///// ///// /////	Hazard ratios
	/*
	Tricky calculation. You are calculating the ratio of mean hazards: i.e. the mean hazard values in the vaccine group, divided by the mean hazard values in the control group. 
	You will do this for all age groups, "countries", and serostatuses, and  serotypes, for a number of example days. 
	The idea of the calculation is as follows: 
		i)		Sum the value of the hazard on a particular day post dose (which will correspond to different calendar days). Must be threaded and then the threads must be summed. Must remember to clear threaded quantities. 
					Store results in MeanHazVals instance of SCurves class. 
		ii)		Divide the summed hazard by the size of the strata, which you are already calculating for survival curves (specifically SC_WT instance of SCurves class). 
					Store results in MeanHazVals instance of SCurves class. You don't need curve tails for this either. You only need curve tails for RATIOS. Split up existing functions accordingly. 
		iii)	Now have mean hazards in vaccine and control groups (for each age etc.). Must then calculate their ratio. Store results in HRs instance of SCurves class. 
		iv)		Now have single posterior sample for hazard ratios for all strata you care about for a selection of days post dose. Then must apply functions to add to the mean and accept or reject value in to the tails. 
				

		To use Survival curve machinery (e.g. multiple (threaded) tables by disease severity, functions that sum/reset threads etc.), you create multiple instances of SCurves for the mean hazards and the mean hazard ratios. 
		This means that you will create some things you don't need, like strata sizes (you won't add to them as calculating the exact same quantity in SC_WT instance of SCurves class), but still worth it. 
		Also you initialize the hazard RATIOS (HRs instance of SCurves class) differently: they don't need strata broken down by trial arm as they're always comparing vaccine over control. 
	*/

	SCurves * HRs_BySerotype			= NULL;		//// Hazard ratios - you don't separate these for whole trial. You will initialize these with different (but similar) values of NumStrata etc., as you don't need separate vaccine or control arms, only their ratio. 
	SCurves * MeanHazVals_BySerotype	= NULL;		//// Hazard ratios - you don't separate these for whole trial. You will use the same strata here as for SC_WT
	std::vector<string> HRsRowNames, HRsColNames; 
	std::vector<string> HRs_SeroNames; 


	int HRs_NumDaysPostDoseToCalculate = 0; //// you will add to this in HRs_DaysPostDose.push_back below. 
	std::vector<int> HRs_DaysPostDose;

	int HRs_NumHRsPerAgeGroupAndCountry = 3; // want Vac/Control for i) seronegatives; ii) seropositives; iii) any serostatus; 
	int HRs_NumAgeGroups = 8;
	int HRs_NumCountries = NumC + 3; //// Countries also refers to CYD14, CYD15, and CYD14 and CYD15 combined. 
	int HRs_NumStrata = HRs_NumAgeGroups * HRs_NumCountries * HRs_NumHRsPerAgeGroupAndCountry;
	int HRs_NumSTypes = 1; //// if not considering serotype, this will be 1. Otherwise will be HOUSE.N_STypes + 1, as you'll want hazard of combined serotypes too. Importantly, HR_serotype = 0 will refer to "combined" serotypes, not serotype 1 as you have everywhere else. Individual serotypes 1 to 4 will be labelled 1 to 4, contrary to usual Cpp indexing. 

	///// For separate Passive phase Survival curves. 
	bool CalculateSeparatePassivePhaseCurves	= 1 ;

	///// ///// ///// ///// /////	Attack Rates

	DType **** MetaAttackRates			= NULL;		//// indexed by				i) trial phase;  ii) disease severity; iii) Hazard Group/Strata  (e.g. seropositive either arm), and iv) Posterior Sample number
	DType **** Threaded_AttackRates		= NULL;		//// indexed by i) thread; ii) trial phase; iii) disease severity; iv)  Hazard Group/Strata  (e.g. seropositive either arm). NOTE: DON'T NEED posterior sample number as will reset thread specific meta attack rates with each iteration. 
	DType *** Sum_FU_Durations			= NULL;		//// For attack rate denominators: indexed by i) thread; ii) trial phase; iii) Hazard Group/Stratum. NOTE: DON'T NEED Disease Severity as follow up durations same regardless of whether considering mild, severe etc. 
	DType **** MeanModeAttackRates		= NULL;		//// indexed by				i) trial phase;  ii) disease severity; iii) Hazard Group/Strata  (e.g. seropositive either arm), and iv) Statistic (mean or modal posterior). 

	int Num_AR_SummaryStats = 3; /// 3 because you will want mean posterior, mode posterior, and max likelihood. 

	void init		(const DATA_struct &DATA, const Housekeeping_Struct &HOUSE, const Chains_Struct &CHAINS, const Chains_Struct &WBIC_CHAINS)
	{
		if (CHAINS.Calc_SCsARsHRPs_Any || WBIC_CHAINS.Calc_SCsARsHRPs_Any)
		{
#ifdef PRINT_PROGRAM_PROGRESS
			std::cerr << "InitializeSurvivalCurveDataStructure ";
#endif
			//// Initialize SCurves Strata - hard code the nasty old labelling of survival curve strata. Hard code it here and then it should work elsewhere. 
			/*		HAZARD-GROUP NUMBERS / DEFINITIONS (REALLY IRRITATING THAT THESE ARE NOT THE SAME AS OTHER STRATA NUMBERS -THE JOYS OF LEGACY CODE)
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

			Allocate_2D_Array(Strata_SCurves, HOUSE.HowManyTrialArms + 1, HOUSE.HowManySeroStatuses + 1);
			Strata_SCurves[ControlGroup		][SeroNeg			] = 1; 
			Strata_SCurves[ControlGroup		][SeroPos			] = 0;
			Strata_SCurves[ControlGroup		][EitherSeroStatus	] = 4;
			Strata_SCurves[VaccineGroup		][SeroNeg			] = 3;
			Strata_SCurves[VaccineGroup		][SeroPos			] = 2;
			Strata_SCurves[VaccineGroup		][EitherSeroStatus	] = 5;
			Strata_SCurves[EitherTrialArm	][SeroNeg			] = 7;
			Strata_SCurves[EitherTrialArm	][SeroPos			] = 6;
			Strata_SCurves[EitherTrialArm	][EitherSeroStatus	] = 8;

			///// CHOOSE Num Days of follow up
			if (HOUSE.SFU)
			{
				if (HOUSE.PASSIVE_PHASE_ONLY)
				{
					NoDaysOfFollowUp		= 336;
					NoActiveDaysFollowUp	= NoDaysOfFollowUp;
					NoPassiveDaysFollowUp	= 0;
				}
				else 	if (HOUSE.ActiveOrPassivePhase == ACTIVE_PHASE_ONLY)
				{
					NoDaysOfFollowUp		= DATA.NoDaysOfFollowUp;
					NoActiveDaysFollowUp	= NoDaysOfFollowUp;
					NoPassiveDaysFollowUp	= 0;
				}
				else	if (HOUSE.ActiveOrPassivePhase == DO_ACTIVE_AND_PASSIVE)
				{
					NoPassiveDaysFollowUp = 335;
					NoActiveDaysFollowUp = 760;

					if (!HOUSE.Include_Late_Cases)
					{
						NoDaysOfFollowUp		= 1227; 
					}
					else if (HOUSE.Include_Late_Cases && !HOUSE.FakeExtObs)
					{
						NoDaysOfFollowUp		= 1227; 
					}
					else if (HOUSE.Include_Late_Cases && HOUSE.FakeExtObs)
					{
						NoDaysOfFollowUp		= 1372;  
					}
				}
			}
			else
			{
				if (HOUSE.PASSIVE_PHASE_ONLY)
				{
					NoDaysOfFollowUp		= 336;
					NoActiveDaysFollowUp	= NoDaysOfFollowUp;
					NoPassiveDaysFollowUp	= 0;
				}
				else 	if (HOUSE.ActiveOrPassivePhase == ACTIVE_PHASE_ONLY)
				{
					NoDaysOfFollowUp		= 730;
					NoActiveDaysFollowUp	= DATA.NoActiveDaysFollowUp;
					NoPassiveDaysFollowUp	= 0;
				}
				else	if (HOUSE.ActiveOrPassivePhase == DO_ACTIVE_AND_PASSIVE)
				{
					NoActiveDaysFollowUp	= DATA.NoActiveDaysFollowUp;
					NoPassiveDaysFollowUp	= 365;

					if (!HOUSE.Include_Late_Cases)
					{
						NoDaysOfFollowUp		= DATA.NoDaysOfFollowUp; 
					}
					else if (HOUSE.Include_Late_Cases)	
					{
						NoDaysOfFollowUp		= DATA.NoDaysOfFollowUp; /// this will help with survival curve output
					}
				}
			}	

			if (HOUSE.ActiveOrPassivePhase == ACTIVE_PHASE_ONLY) CalculateSeparatePassivePhaseCurves = 0; //// guard. 

			if (HOUSE.MildAndSevere == TREATED_EQUALLY) HowManyDiseaseSeverities = 1; else HowManyDiseaseSeverities = 3; //// 1 -> {Either}, 3 -> {Either, Mild, Severe}
			if (HOUSE.MildAndSevere == TREATED_EQUALLY) DiseaseNames.push_back(""); 
			else	{									DiseaseNames.push_back("_Either");			DiseaseNames.push_back("_Mild");			DiseaseNames.push_back("_Severe");	}
			StatisticNames.push_back("Mean_");	StatisticNames.push_back("LowerCI_");		StatisticNames.push_back("UpperCI_");

			if (!HOUSE.PASSIVE_PHASE_ONLY) TrialPhaseNames.push_back("_ActivePhaseOnly"); else TrialPhaseNames.push_back(""); //// Because you relabel the passive phase as the "active" phase. Labels will look weird unless you accout for this. E.g. you would get "AttackRates_VAC_SILENT_PASSIVE_ONLY_Cs01234_PooledCountries_SFU_Dummy_ActivePhaseOnly.txt", which is confusing. 
			if (HOUSE.ActiveOrPassivePhase == DO_ACTIVE_AND_PASSIVE)
			{
				TrialPhaseNames.push_back("_PassivePhaseOnly");
				TrialPhaseNames.push_back(""); //// "" corresponds to whole trial. 
			}

			if (HOUSE.ActiveOrPassivePhase == ACTIVE_PHASE_ONLY) HowManyTrialPhases = 1; else HowManyTrialPhases = HOUSE.HowManyTrialPhases + 1;  

			NoSurvivePostSamples						= CHAINS.NumPosteriorSamples;				//// default. 
			AddToSurvivalCurvesEveryHowManyIterations	= CHAINS.AddtoChainEveryHowManyIterations;	//// default. 

			if (NoSurvivePostSamples > MaxNumSurvivalPosteriorSamples)
			{
				AddToSurvivalCurvesEveryHowManyIterations *= (CHAINS.NumPosteriorSamples / MaxNumSurvivalPosteriorSamples); //// i.e. if you have more samples than you need, must increase no. iterations before you add to chain) add to chain less often. 

				//// redefine after resetting AddtoChainEveryHowManyIterations
				NoSurvivePostSamples = (CHAINS.No_Iterations - CHAINS.BurnIn) / AddToSurvivalCurvesEveryHowManyIterations;
			}
			NumElementsOutside_CrI_Tails = (int) ((DType) NoSurvivePostSamples * CHAINS.alpha * 0.5 + 1);

			//// check
			if (AddToSurvivalCurvesEveryHowManyIterations % CHAINS.AddtoChainEveryHowManyIterations != 0) std::cerr << "Survival update iterations not multiple of Parameter update iterations" << endl << endl << endl << endl << endl << endl << endl << endl << endl << endl;
			std::cerr << "NumElementsOutside_CrI_Tails " << NumElementsOutside_CrI_Tails << " NoSurvivePostSamples " << NoSurvivePostSamples << " AddToSurvivalCurvesEveryHowManyIterations " << AddToSurvivalCurvesEveryHowManyIterations << endl;

			//// Rownames of survival tables
			for (int agegroup = 0; agegroup < NumAgeGroupCategories; agegroup++)	// rownames
				for (int country = 0; country < (HOUSE.TotalCountries + 3); country++)
					for (int hazardgroup = 0; hazardgroup < TotalGroupsPerCountry; hazardgroup++)
						SurvivalTableRowNames.push_back("A_" + std::to_string(agegroup) + " C_" + std::to_string(country) + " Hz_" + std::to_string(hazardgroup + 1));	// hazardgroup + 1 as that's how various R scripts are coded.
			// Colnames of survival tables
			for (int timepoint = 0; timepoint < NoDaysOfFollowUp + 1; timepoint++) SurvivalTableColNames.push_back("SurvivalProb_t_" + std::to_string(timepoint)); // colnames

			// Attack_Rate_Colnames
			for (int sample = 0; sample < NoSurvivePostSamples; sample++)	AttackRateColNames.push_back("Sample_" + std::to_string(sample));

			std::cerr << "Whole Trial initialization" << endl; 
			//// Initialize Whole Trial survival curve and associated quantities. 
			SC_WT.init(NoSubjectCategories, NoDaysOfFollowUp + 1, HowManyDiseaseSeverities, HOUSE.max_threads, NumElementsOutside_CrI_Tails); 

			//// //// //// //// //// //// //// //// Passive Survival Curve (same as above, but NoPassiveDays is different). 
			if (CalculateSeparatePassivePhaseCurves)
			{
				// Colnames
				for (int timepoint = 0; timepoint < NoPassiveDaysFollowUp + 1; timepoint++) 
					PassivePhaseSurvivalTableColNames.push_back("SurvivalProb_t_" + std::to_string(timepoint)); // colnames

				std::cerr << "Passive Phase initialization" << endl;
				//// initialize passive phase Survival curves and associated quantities. 
				SC_PP.init(NoSubjectCategories, NoPassiveDaysFollowUp + 1, HowManyDiseaseSeverities, HOUSE.max_threads, NumElementsOutside_CrI_Tails); 
			}

			///// ///// ///// /////				Attack Rates
			std::cerr << "Attack Rates initialization" << endl;
			Allocate_4D_Array(MetaAttackRates							, HowManyTrialPhases, HowManyDiseaseSeverities	, NoSubjectCategories, NoSurvivePostSamples);	//// HowManyTrialPhases as you want WholeTrial. 
			Allocate_4D_Array(Threaded_AttackRates	, HOUSE.max_threads	, HowManyTrialPhases, HowManyDiseaseSeverities	, NoSubjectCategories);									//// HowManyTrialPhases as you want WholeTrial. 
			Allocate_3D_Array(Sum_FU_Durations		, HOUSE.max_threads	, HowManyTrialPhases							, NoSubjectCategories); 
			Allocate_4D_Array(MeanModeAttackRates						, HowManyTrialPhases, HowManyDiseaseSeverities	, NoSubjectCategories, Num_AR_SummaryStats);	
			std::cerr << "DONE ";

			///// ///// ///// /////				Hazard Ratios initialization
			std::cerr << "Hazard Ratios initialization" << endl;
			HRs_DaysPostDose.push_back(0);		 
			HRs_DaysPostDose.push_back(365);	 
			HRs_DaysPostDose.push_back(730);	 
			HRs_DaysPostDose.push_back(760);
			HRs_DaysPostDose.push_back(910);
			HRs_DaysPostDose.push_back(1095);	 
			HRs_NumDaysPostDoseToCalculate = HRs_DaysPostDose.size(); 

			//// Rownames & colnames of HR tables
			for (int agegroup = 0; agegroup < HRs_NumAgeGroups; agegroup++)	// rownames
				for (int country = 0; country < HRs_NumCountries; country++)
					for (int BaselineSeroStatus = 0; BaselineSeroStatus < HRs_NumHRsPerAgeGroupAndCountry; BaselineSeroStatus++)
						HRsRowNames.push_back("A_" + std::to_string(agegroup) + "_C_" + std::to_string(country) + "_BS_" + std::to_string(BaselineSeroStatus));
			for (int year = 0; year < HRs_NumDaysPostDoseToCalculate; year++) 
				HRsColNames.push_back("D_" + std::to_string(HRs_DaysPostDose[year])); // colnames

			if (HOUSE.SeroSpecificEfficacies || HOUSE.SeroSpecific_K_values) HRs_NumSTypes = HOUSE.N_STypes + 1; else HRs_NumSTypes = 1; 
			
			HRs_BySerotype			= new SCurves[HRs_NumSTypes]; 
			MeanHazVals_BySerotype	= new SCurves[HRs_NumSTypes]; 
			for (int HR_serotype = 0; HR_serotype < HRs_NumSTypes; HR_serotype++) 
			{
				if (HR_serotype == 0)	HRs_SeroNames.push_back(""); else HRs_SeroNames.push_back("_sero" + std::to_string(HR_serotype));
				HRs_BySerotype			[HR_serotype].init(NoSubjectCategories, HRs_NumDaysPostDoseToCalculate, HowManyDiseaseSeverities, HOUSE.max_threads, NumElementsOutside_CrI_Tails);
				MeanHazVals_BySerotype	[HR_serotype].init(NoSubjectCategories, HRs_NumDaysPostDoseToCalculate, HowManyDiseaseSeverities, HOUSE.max_threads, NumElementsOutside_CrI_Tails);
				MeanHazVals_BySerotype	[HR_serotype].Strata_Sizes = SC_WT.Strata_Sizes; //// So that you can refer to MeanHazVals_BySerotype	[HR_serotype].Strata_Sizes but not calculate the same thing repeatedly, set the pointer. 
			}
			std::cerr << "DONE ";
		}
	}
	Survival_Struct	(const DATA_struct &DATA, const Housekeeping_Struct &HOUSE, Survival_Struct &SURVIVE, const Chains_Struct &CHAINS, const Chains_Struct &WBIC_CHAINS)
	{
		init(DATA, HOUSE, CHAINS, WBIC_CHAINS);
	}
	Survival_Struct() {}; //// this default constructor does nothing. But will allow you to create instance of structure, then initialize later. Useful e.g. when creating array of pointers to type "this structure". 
};
struct Augment_Struct {

	int NoLikeDiffsPerAugPatient = 0;
	std::vector<int> AugmentedIndices, NonAugmentedIndices;
	DType ** LikeDiffs			= NULL;		//// Differences made to individual likelihood components through augmentation. indexed by i) patient; ii) component 

	bool Aug_AllCool = 1;

	DType ** Mults				= NULL;		//// indexed by i) thread; ii) Baseline Serostatus. // augmentations basically the same: in one case you add, in the other you subtract. This gives what used to be classed as Multiplier, and -Multiplier
	DType *** iBaseHaz_Mult		= NULL;		//// indexed by i) thread; ii) PhaseSeverity; iii) Baseline Serostatus (not number of previous infections, so will include a K_seropos value fairly often). Usually this is equal to Survive Ks, but if HOUSE.ResidEffs, then you multiply the integrated baseline hazard by (1 - VE_inf), in which case it doesn't equal the "SurviveK"
	DType ** Case_Ks			= NULL;		//// indexed by i) thread; ii) Baseline Serostatus (not number of previous infections, so will include a K_seropos value fairly often). Don't need PhaseSeverity index as only one type of cases. 
	//DType **** sero_Ks			= NULL;		//// indexed by i) thread; ii) PhaseSeverity; iii) Baseline Serostatus (not number of previous infections, so will include a K_seropos value fairly often) and iv) serotype (needed only for running SS_VEs and SS_Ks together as a multiplier of the vaccine hazard (not covered with iVacHacHaz_mult as this is the "aggregate multiplier" - you need individual serotypes for quick likelihood calculation/augmentation (hence both serostatuses). ). 
	DType ** iBaseHaz			= NULL;		//// indexed by i) thread; ii) PhaseSeverity. IntBase Haz. 
	DType *** iVacHaz			= NULL;		//// IntVac Haz: indexed by i) thread; ii) PhaseSeverity; iii) Baseline Serostatus

	DType ** P_Survive					= NULL;		//// indexed by i) thread; ii) Baseline Serostatus (used for Gibbs aug). 
	DType ** HazCaseMult				= NULL;		//// indexed by i) thread; ii) Baseline Serostatus (used for Gibbs aug). 
	DType ** P_BS						= NULL;		//// Prob(baseline serostatus); indexed by i) thread; ii) Baseline Serostatus (used for Gibbs aug). And Total likelihood. 
	DType ** P_OGivenS_PS				= NULL;		//// P(0|S)P(S); indexed by i) thread; ii) Baseline Serostatus (used for Gibbs aug). And Total likelihood. 
	DType *  P_Outcome					= NULL;		//// i.e. P(0|S-)P(S-) + P(0|S+)P(S+); indexed by i) thread


	///// NOTE: You don't use iVacHaz_Mult for the Change_Aug_Data function, but you do use iVacHaz_Mult for Survival curves and for ProbSerostatusGivenOutcome. Confusing but important. Essentially iVacHaz_Mult refers to the multiplier after all serotypes taken into account (i.e. having survived all four serotypes). 
	///// Means that any effects like age or whatever have to be coded into iVacHaz_Mult (for ProbSerostatusGivenOutcome and for Survival curves) as well as coded into Change_Aug_Data function. 
	///// Bit hacky and confusing, but iVacHaz_Mult mostly refers to what it sounds like, apart from when doing SS_VEs (with or without SS_Ks), because here you divide up the likelihood into individual serotype parts. In this case iVacHaz_Mult does not multiply the integrated vaccine hazard. iVacHaz_Mult is used for Survival curves and for ProbSerostatusGivenOutcome. Confusing but important. 
	DType *** iVacHaz_Mult		= NULL;		//// indexed by i) thread; ii) PhaseSeverity; iii) Baseline Serostatus (not number of previous infections, so will include a K_seropos value fairly often). For SS_Ks and SS_VEs together, within particular hazard group, the "survivor K value" no longer same for both the integrated baseline hazard and integrate vaccine hazard. IntBase Haz has Survive K = sum(rho x K), and IntVacHaz has VAC_Survive_K = sum(rho x K x Eff). In all other cases, these terms equal (although perhaps age-specific vaccine efficacies may need attention). 

	///// For Hazard ratios
	DType **** HR_VacHaz_seroMults = NULL;		//// (needed for HazardRatios) indexed by i) thread; ii) PhaseSeverity; iii) Baseline Serostatus (not number of previous infections, so will include a K_seropos value fairly often) and iv) serotype 
	DType **** HR_BaseHaz_seroMults = NULL;		//// (needed for HazardRatios) indexed by i) thread; ii) PhaseSeverity; iii) Baseline Serostatus (not number of previous infections, so will include a K_seropos value fairly often) and iv) serotype 

	void init(Housekeeping_Struct &HOUSE, const DATA_struct &DATA)
	{
		NoLikeDiffsPerAugPatient = HOUSE.LCPerC;

		// Array of likelihood differences or each augmented data patient. Used for parallelizing MCMC. If there are X LCPerC, then there are Y=X times NumC likelihood components. Therefore there could be Y * NoAugmented entries to multi array. This is wasteful as each patient can only contribute to their country's likelihood component. Thus the array should be LCPerC + 1 rows by NoAugmented columns. plus one because you want LikeFull as well
		Allocate_2D_Array(LikeDiffs		, DATA.NoAugmented	, HOUSE.LCPerC + 1);

		Allocate_2D_Array(Mults			, HOUSE.max_threads	, HOUSE.HowManySeroStatuses);
		Allocate_2D_Array(Case_Ks		, HOUSE.max_threads	, HOUSE.HowManySeroStatuses);
		//Allocate_3D_Array(Survive_Ks	, HOUSE.max_threads	, HOUSE.HowManyCaseCategories, HOUSE.HowManySeroStatuses); //// (not number of K params).
		Allocate_3D_Array(iBaseHaz_Mult	, HOUSE.max_threads	, HOUSE.HowManyCaseCategories, HOUSE.HowManySeroStatuses); //// (not number of K params).
		Allocate_2D_Array(iBaseHaz		, HOUSE.max_threads	, HOUSE.HowManyCaseCategories); //// 
		Allocate_3D_Array(iVacHaz		, HOUSE.max_threads	, HOUSE.HowManyCaseCategories, HOUSE.HowManySeroStatuses); //// 

		Allocate_2D_Array(P_Survive					, HOUSE.max_threads	, HOUSE.HowManySeroStatuses);
		Allocate_2D_Array(HazCaseMult				, HOUSE.max_threads	, HOUSE.HowManySeroStatuses);
		Allocate_2D_Array(P_BS						, HOUSE.max_threads	, HOUSE.HowManySeroStatuses);
		Allocate_2D_Array(P_OGivenS_PS				, HOUSE.max_threads	, HOUSE.HowManySeroStatuses);
		P_Outcome = new DType[HOUSE.max_threads](); 

		Allocate_3D_Array(iVacHaz_Mult	, HOUSE.max_threads	, HOUSE.HowManyCaseCategories, HOUSE.HowManySeroStatuses); //// (not number of K params).
	
		if (HOUSE.SeroSpecificEfficacies || HOUSE.SeroSpecific_K_values) //// Previously you needed these only when both SS_VEs and SS_Ks. Now you need them for hazard ratios. 
		{
			Allocate_4D_Array(HR_VacHaz_seroMults	, HOUSE.max_threads, HOUSE.HowManyCaseCategories, HOUSE.HowManySeroStatuses, HOUSE.N_STypes);
			Allocate_4D_Array(HR_BaseHaz_seroMults	, HOUSE.max_threads, HOUSE.HowManyCaseCategories, HOUSE.HowManySeroStatuses, HOUSE.N_STypes);
		}

	}
}; 


#endif
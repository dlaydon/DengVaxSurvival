
#include "HeaderAndSwitches.h"
#include "ParamNumber.h"
#include "Splines.h"
#include "ReadAndProcessData.h"
#include "Augmentation.h"
#include "Probability.h"
#include "CalcParamsEtc.h"

int set_NCores()
{
	char * val;
	char * endchar;
	val = getenv("CCP_NUMCPUS");
	int cores = omp_get_max_threads();
	if (val != NULL) cores = strtol(val, &endchar, 10);
	omp_set_num_threads(cores);
	return cores;
}
Age_Option			Convert_AS_String			(const std::string& OptionString)
{
		 if (OptionString == "INDEPENDENT"	) return Age_Option::INDEPENDENT	;
	else if (OptionString == "HILL"			) return Age_Option::HILL			;
	else if (OptionString == "CATEGORICAL"	) return Age_Option::CATEGORICAL	;
	else if (OptionString == "SPLINE"		) return Age_Option::SPLINE			;
	else if (OptionString == "SPLINE_LINE"	) return Age_Option::SPLINE_LINE	;
	else if (OptionString == "SPLINE_STEP"	) return Age_Option::SPLINE_STEP	;
	else if (OptionString == "CUBIC"		) return Age_Option::CUBIC			;
	else std::cerr << endl << "Convert_AS_String ERROR: String not recognized" << endl; 
}
ExtImSub_Option		Convert_ExtImSub_String		(const std::string& OptionString)
{
		 if (OptionString == "IGNORED"	) return ExtImSub_Option::IGNORED		;
	else if (OptionString == "AS_DATA"	) return ExtImSub_Option::AS_DATA		;
	else if (OptionString == "AS_PROB"	) return ExtImSub_Option::AS_PROB		;
	else std::cerr << endl << "Convert_AS_String ERROR: String not recognized" << endl; 
}
EffNegWane_Option	Convert_EffNegWane_String	(const std::string& OptionString)
{
		 if (OptionString == "DEFAULT"		) return EffNegWane_Option::DEFAULT			;
	else if (OptionString == "FROM_ZERO"	) return EffNegWane_Option::FROM_ZERO		;
	else if (OptionString == "NO_WANE"		) return EffNegWane_Option::NO_WANE			;
	else std::cerr << endl << "Convert_NegEffWane_String ERROR: String not recognized" << endl;
}

template <typename TYPE> void DeAllocate_2D_Array	(TYPE **  &OBJECT, int Dim1) 
{
	for (int row = 0; row < Dim1; row++) delete[] OBJECT[row];
	delete[] OBJECT;
}
template <typename TYPE> void DeAllocate_3D_Array	(TYPE *** &OBJECT, int Dim1, int Dim2) 
{
	for (int row = 0; row < Dim1; row++) DeAllocate_2D_Array(OBJECT[row], Dim2);
	delete[] OBJECT;
}
 
void CreateAndPopulateKnots				(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	Allocate_2D_Array(PARAMS.xKnots, HOUSE.TotalCountries, HOUSE.KnotsPerCountry);
	Allocate_2D_Array(PARAMS.yKnots, HOUSE.TotalCountries, HOUSE.KnotsPerCountry);

	DType NoUnitsInYear = 1, IntervalsPerUnit = 3; 

	bool ImportCondition = HOUSE.SFU & !HOUSE.PASSIVE_PHASE_ONLY & HOUSE.ActiveOrPassivePhase == DO_ACTIVE_AND_PASSIVE;

	if (ImportCondition)
	{
		ifstream KnotInput; KnotInput.open(HOUSE.KnotsInputFilename);
		for (int country = 0; country < HOUSE.TotalCountries; country++)
			for (int knot = 0; knot < HOUSE.KnotsPerCountry; knot++)
				KnotInput >> PARAMS.yKnots[country][knot];
		KnotInput.close();
	}

	for (int country = 0; country < HOUSE.TotalCountries; country++)
		for (int knot = 0; knot < HOUSE.KnotsPerCountry; knot++)
		{
			PARAMS.xKnots[country][knot] = (NoUnitsInYear / IntervalsPerUnit) * (knot + 1);
			if (!ImportCondition)
				PARAMS.yKnots[country][knot] = 0.05;
		}
}
void CreateAndPopulateSplineCoeffs		(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	Allocate_3D_Array(PARAMS.SplineCoeffs, HOUSE.TotalCountries, HOUSE.PolynomialsPerCountry, HOUSE.MaxSplineDegree + 1);

	for (int country = 0; country < HOUSE.TotalCountries; country++)
		for (int poly = 0; poly < HOUSE.PolynomialsPerCountry; poly++)
			CoefficientsFromKnots(PARAMS.xKnots[country], PARAMS.yKnots[country], poly, PARAMS.SplineCoeffs[country][poly]);
}

void Initialize_IntVacHazardValues		(Params_Struct &PARAMS, const DATA_struct &DATA, const Housekeeping_Struct &HOUSE)
{
	/*
		tempting to make this part of Initialize Params - but don't. 
		Calc_IVH_Values_c_BS_PS function below relies on Sets (Spos/Sneg by arm and country) being properly calculated (need to know serostatus of each person)
		The initial assignment of Serostatus requires historical hazards, which are part of parameters
		However your augmented data is in the DATA structure. So DO. NOT. TOUCH
	*/
	Allocate_2D_Array(PARAMS.IVH_vals, NPat, HOUSE.HowManyCaseCategories);

	bool ConditionToSkip; 
	for (int country = 0; country < HOUSE.TotalCountries; country++)
	{
		ConditionToSkip = !(std::any_of(HOUSE.WhichCountries.begin(), HOUSE.WhichCountries.end(), [&](int j) {return j == country; }));	//// if not doing a particular country, skip that country
		if (ConditionToSkip) continue;

		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
				Calc_IVH_Values_c_BS_PS(country, BaselineSeroStatus, PhaseSeverity, PARAMS, DATA, HOUSE);
	}
}
void Initialize_IntVacHazardValues		(Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const DATA_struct &DATA, const Housekeeping_Struct &HOUSE) 
{
	Initialize_IntVacHazardValues(CurrentPARAMS , DATA, HOUSE);
	Initialize_IntVacHazardValues(ProposedPARAMS, DATA, HOUSE);
}

void Initialize_KplusValues				(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	Allocate_3D_Array(PARAMS.KplusValues, HOUSE.TotalCountries, HOUSE.HowManyAges, HOUSE.HowManyCaseCategories);

	// Allocate Meta_KplusValues. Useful as there is one scenario where the (survivor) KplusValues differ between vaccine and control arms, i.e. in K_SEROPOS variant when doing both HOUSE.SeroSpecificEfficacies && HOUSE.SeroSpecific_K_values, as here the serotype specific efficacies and relative risks need to be in the survivor K's, but obviously not for the control group.
	PARAMS.Meta_KplusValues = new DType ***[HOUSE.HowManyTrialArms]();

	PARAMS.Meta_KplusValues[No_Effs] = PARAMS.KplusValues;

	if (HOUSE.ModelVariant == K_SEROPOS && HOUSE.SeroSpecificEfficacies && HOUSE.SeroSpecific_K_values)
		Allocate_3D_Array(PARAMS.Meta_KplusValues[With_Effs], HOUSE.TotalCountries, HOUSE.HowManyAges, HOUSE.HowManyCaseCategories);
	else PARAMS.Meta_KplusValues[With_Effs] = PARAMS.KplusValues;

	if (HOUSE.ModelVariant != SIMPLE_NUMERICAL) Calc_KplusValues_All(PARAMS, HOUSE);
}
void Initialize_SumRhoKs				(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	Allocate_3D_Array(PARAMS.SumRhoKs, HOUSE.TotalCountries, HOUSE.HowManyCaseCategories, HOUSE.Num_K_Params);
	Calc_SumRhoKs_ALL(PARAMS , HOUSE);
}
void InitializeLikeArray				(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	PARAMS.LikeSize	= (HOUSE.LCPerC * HOUSE.TotalCountries);
	Allocate_2D_Array(PARAMS.LikeParts, HOUSE.TotalCountries, HOUSE.LCPerC);

	if (HOUSE.Weighting_Pass_Sev)
	{
		Allocate_2D_Array(PARAMS.w_LikeParts, HOUSE.TotalCountries, HOUSE.LCPerC);

		PARAMS.ptr_LikeParts	= PARAMS.w_LikeParts;
		PARAMS.ptr_L_Full		= &PARAMS.w_LikeFull;
	}
	else
	{
		PARAMS.ptr_LikeParts	= PARAMS.LikeParts;
		PARAMS.ptr_L_Full		= &PARAMS.LikeFull;
	}	
}
void Initialize_ParamSeroPrevs			(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, const SeroPrev_Struct &SEROPREV)
{
	///// This only does anything if (HOUSE.Empirical_SeroPrevs), i.e. if using Empirical Seroprevalences in the likelihood and augmentation, rather than exp(-ha) or (1 - exp(-ha)) as usual. 
	///// PARAMS.SeroPrevs array is used in likelihood either way, but has been allocated (and populated using historical hazards as default) in Initialize_Params function. 
	if (HOUSE.Empirical_SeroPrevs)
	{
		//// take P(SPos) from all countries and ages of SEROPREV AgeSpecificSeroPrevs. Use NonAugIndex as this refers to immuno subset, i.e. patients for whom augmentation not required. 
		for (int country = 0; country < HOUSE.TotalCountries; country++)
			for (int age = 0; age < HOUSE.HowManyAges; age++)
					 if (SEROPREV.AgeSpecificSeroPrevs[country][NonAugIndex][age] == 1) 
					PARAMS.SeroPrevs[NonLogIndex][country][SeroPos][age] = 0.99;
				else if (SEROPREV.AgeSpecificSeroPrevs[country][NonAugIndex][age] == 0)
					PARAMS.SeroPrevs[NonLogIndex][country][SeroPos][age] = 0.01;
				else
					PARAMS.SeroPrevs[NonLogIndex][country][SeroPos][age] = SEROPREV.AgeSpecificSeroPrevs[country][NonAugIndex][age];

		//// take P(SNeg) = 1 - P(SPos). 
		for (int country = 0; country < HOUSE.TotalCountries; country++)
			for (int age = 0; age < HOUSE.HowManyAges; age++)
				PARAMS.SeroPrevs[NonLogIndex][country][SeroNeg][age] = 1 - PARAMS.SeroPrevs[NonLogIndex][country][SeroPos][age];

		//// calculate log(P(SNeg)) and log(P(SPos)) for each age. 
		for (int country = 0; country < HOUSE.TotalCountries; country++)
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
				for (int age = 0; age < HOUSE.HowManyAges; age++)
					PARAMS.SeroPrevs[LogIndex][country][BaselineSeroStatus][age] = log(PARAMS.SeroPrevs[NonLogIndex][country][BaselineSeroStatus][age]);
	}
}
void Initialize_ParamSeroPrevs			(Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE, const SeroPrev_Struct &SEROPREV)
{
	Initialize_ParamSeroPrevs(CurrentPARAMS , HOUSE, SEROPREV);
	Initialize_ParamSeroPrevs(ProposedPARAMS, HOUSE, SEROPREV);
}

//// subroutines of Initialize_Params function 
void InitializeCommandLineParams		(Params_Struct &PARAMS, Params_Struct &CurrentPARAMS)
{
	PARAMS.KA_0									= CurrentPARAMS.KA_0								;
	PARAMS.KA_1									= CurrentPARAMS.KA_1								;
	PARAMS.KA_2									= CurrentPARAMS.KA_2								;
	PARAMS.KH_0									= CurrentPARAMS.KH_0								;
	PARAMS.KH_1									= CurrentPARAMS.KH_1								;
	PARAMS.KH_2									= CurrentPARAMS.KH_2								;
	PARAMS.PosEfficacy							= CurrentPARAMS.PosEfficacy							;
	PARAMS.NegEfficacy							= CurrentPARAMS.NegEfficacy							;
	PARAMS.PosWaning							= CurrentPARAMS.PosWaning							;
	PARAMS.NegWaning							= CurrentPARAMS.NegWaning							;
	PARAMS.Hosp_K_mult							= CurrentPARAMS.Hosp_K_mult							;				
	PARAMS.Initial_BS_BaseHazMult				= CurrentPARAMS.Initial_BS_BaseHazMult				; 
	PARAMS.Hosp_K_mult							= CurrentPARAMS.Hosp_K_mult							; 
	PARAMS.Initial_Halflife_Hill				= CurrentPARAMS.Initial_Halflife_Hill				; 
	PARAMS.Initial_Power_Hill					= CurrentPARAMS.Initial_Power_Hill					; 
	PARAMS.Initial_Halflife_ASVE				= CurrentPARAMS.Initial_Halflife_ASVE				; 
	PARAMS.Initial_Power_ASVE					= CurrentPARAMS.Initial_Power_ASVE					; 
	PARAMS.Initial_Prop_ASVE					= CurrentPARAMS.Initial_Prop_ASVE					; 
	PARAMS.Initial_Halflife_AS_Haz				= CurrentPARAMS.Initial_Halflife_AS_Haz				; 
	PARAMS.Initial_Power_AS_Haz					= CurrentPARAMS.Initial_Power_AS_Haz				; 
	PARAMS.Initial_Prop_AS_Haz					= CurrentPARAMS.Initial_Prop_AS_Haz					; 
	PARAMS.Initial_Halflife_AS_Waning			= CurrentPARAMS.Initial_Halflife_AS_Waning			;
	PARAMS.Initial_AS_Prime_PropIndptAge_SNeg	= CurrentPARAMS.Initial_AS_Prime_PropIndptAge_SNeg	;
	PARAMS.Initial_AS_Prime_PropIndptAge_SPos	= CurrentPARAMS.Initial_AS_Prime_PropIndptAge_SPos	;
	PARAMS.Initial_AS_Prime_SNeg_rate			= CurrentPARAMS.Initial_AS_Prime_SNeg_rate			;
	PARAMS.Initial_AS_Prime_SPos_rate			= CurrentPARAMS.Initial_AS_Prime_SPos_rate			;
}
void InitializeRelativeRisks			(int &ParamCounter, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	Allocate_2D_Array(PARAMS.Initial_Ks, HOUSE.HowManyCaseCategories, HOUSE.Num_K_Params);
	PARAMS.Initial_Ks[ActiveMild	][0] = PARAMS.KA_0;
	PARAMS.Initial_Ks[ActiveMild	][1] = PARAMS.KA_1;
	PARAMS.Initial_Ks[ActiveMild	][2] = PARAMS.KA_2;
	
	if (HOUSE.HowManyCaseCategories == 2)
	{
		PARAMS.Initial_Ks[PassiveSevere][0] = PARAMS.KH_0;
		PARAMS.Initial_Ks[PassiveSevere][1] = PARAMS.KH_1;
		PARAMS.Initial_Ks[PassiveSevere][2] = PARAMS.KH_2;
	}

	Allocate_4D_Array(PARAMS.K_s, HOUSE.TotalCountries, HOUSE.HowManyCaseCategories, HOUSE.Num_K_Params, HOUSE.N_STypes_Ks);
	string ParamNameString = "";
	std::vector<string> K_sero_strings(HOUSE.N_STypes_Ks, "");
	if (HOUSE.SeroSpecific_K_values)	for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++) K_sero_strings[serotype] = "_sero" + to_string(serotype + 1);
	else								for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++) K_sero_strings[serotype] = "";

	PARAMS.ParamNumbers.Min_K = ParamCounter;
	for (int country = 0; country < HOUSE.TotalCountries; country++)
		for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
				for (int PrevInf = 0; PrevInf < HOUSE.Num_K_Params; PrevInf++)
				{
					//// add to K's array. 
					PARAMS.K_s[country][PhaseSeverity][PrevInf][serotype] = PARAMS.Initial_Ks[PhaseSeverity][PrevInf];		//// i.e. start off all countries and serotypes with same K's

					ParamNameString = ""; 
					ParamNameString = "K" + HOUSE.PhaseSeverity_strings[PhaseSeverity] + "_" + to_string(PrevInf) + K_sero_strings[serotype];

					if (country == 0) //// when not modelling hostpitalized disease, do not add to ParamVec, NamesOfParameters for every country. 
					{
						// Add to Paramvec
						PARAMS.ParamVec.push_back(PARAMS.K_s[country][PhaseSeverity][PrevInf][serotype]);

						//// add to NamesOfParameters
						PARAMS.NamesOfParameters.push_back(ParamNameString);
						
						//// Add to Params you'll fit (if appropriate)
						if (!((PhaseSeverity == ActiveMild) && (PrevInf == 1))) //// i.e. if not baseline. don't count KA1 for all serotypes. 
							if (!HOUSE.SkipRelativeRisks)
							{
								if (HOUSE.SSKs_FitMatrix[PhaseSeverity][PrevInf] || serotype == 0) ///// i.e. Either if fitting all serotypes for this particular K_PhaseSeverity_PrevInf (works for either SSKs or !SSKs), OR. After that, do the various "skips" below. 
									if (!(HOUSE.Skip_KA_0 && PhaseSeverity == ActiveMild	&& PrevInf == 0))
									if (!(HOUSE.Skip_KA_2 && PhaseSeverity == ActiveMild	&& PrevInf == 2))
									if (!(HOUSE.Skip_KH_0 && PhaseSeverity == PassiveSevere	&& PrevInf == 0))
									if (!(HOUSE.Skip_KH_1 && PhaseSeverity == PassiveSevere	&& PrevInf == 1))
									if (!(HOUSE.Skip_KH_2 && PhaseSeverity == PassiveSevere	&& PrevInf == 2))
										if (!HOUSE.RelRisksSameBtwSerotypes || (HOUSE.Which_Sero_FitsAll == serotype)) //// either not doing Single_SNeg_Eff, or doing Single_SNeg_Eff AND the serotype you're fitting is this one from the loop.
											PARAMS.ParamNosYouWillFit.push_back(ParamCounter);
							}
						// Increment param ParamCounter
						ParamCounter++;
					}
				}
	PARAMS.ParamNumbers.Max_K = ParamCounter - 1;

	if (HOUSE.ModellingHospitalized && (HOUSE.HowManyCaseCategories == 2 || HOUSE.PASSIVE_PHASE_ONLY))
		if (HOUSE.ModelHosp_Indie_Ks)
		{
			int PassSev_PhaseSeverity = (HOUSE.PASSIVE_PHASE_ONLY) ? ActiveMild : PassiveSevere; //// again, this is because you've renamed the passive phsae as the active if fitting PASSIVE_PHASE_ONLY, but PassiveSevere still refers to 1st (i.e. 2nd in Cpp) index
			PARAMS.ParamNumbers.Min_Hosp_K = ParamCounter;
			for (int countryindex = 0; countryindex < HOUSE.CYD_15_countries.size(); countryindex++)
			{
				int Country = HOUSE.CYD_15_countries[countryindex]; 
				for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
					for (int PrevInf = 0; PrevInf < HOUSE.Num_K_Params; PrevInf++)
					{
						if (countryindex == 0) //// do not add to ParamVec, NamesOfParameters for every country.  
						{
							/////	/////	/////	Don't add to K's array, CYD-15 countries can be initialized to have the same hospital K's as CYD-14 countries. 	/////	/////	/////	/////	
							////	Add to Paramvec
							PARAMS.ParamVec.push_back(PARAMS.K_s[Country][PassSev_PhaseSeverity][PrevInf][serotype]);

							ParamNameString = "K" + HOUSE.PhaseSeverity_strings[PassSev_PhaseSeverity] + "_" + to_string(PrevInf) + K_sero_strings[serotype] + "_Hosp";

							////	Add to NamesOfParameters
							PARAMS.NamesOfParameters.push_back(ParamNameString);

							////	Add to Params you'll fit (if appropriate)
							if (!HOUSE.SkipRelativeRisks)
								if (HOUSE.SSKs_FitMatrix[PassSev_PhaseSeverity][PrevInf] || serotype == 0) ///// i.e. Either if fitting all serotypes for this particular K_PhaseSeverity_PrevInf (works for either SSKs or !SSKs), OR. After that, do the various "skips" below. 
									if (!(HOUSE.Skip_KH_0 && PrevInf == 0)) 
									if (!(HOUSE.Skip_KH_1 && PrevInf == 1))
									if (!(HOUSE.Skip_KH_2 && PrevInf == 2))
										if (!HOUSE.RelRisksSameBtwSerotypes || (HOUSE.Which_Sero_FitsAll == serotype)) //// either not doing Single_SNeg_Eff, or doing Single_SNeg_Eff AND the serotype you're fitting is this one from the loop.
											PARAMS.ParamNosYouWillFit.push_back(ParamCounter);
							////// Increment param ParamCounter
							ParamCounter++;
						}
					}
			}
			PARAMS.ParamNumbers.Max_Hosp_K = ParamCounter - 1; 
		}
		else
		{							
			Calc_Hospital_Ks(PARAMS, HOUSE); 

			////	Add to Paramvec
			PARAMS.ParamVec.push_back(PARAMS.Hosp_K_mult);

			////	Add to NamesOfParameters 
			PARAMS.NamesOfParameters.push_back("Hosp_K_mult");

			////	Add to Params you'll fit (if appropriate)
			if (!HOUSE.SkipRelativeRisks)
				if (!HOUSE.Skip_Hosp_Mult)
					PARAMS.ParamNosYouWillFit.push_back(ParamCounter);

			////// Record then Increment param ParamCounter
			PARAMS.ParamNumbers.Hosp_K_multiplier = ParamCounter++; 
		}


	////// Fixed_Severe_RelRisks
	if (HOUSE.Fixed_Severe_RelRisks && (HOUSE.HowManyCaseCategories == 2 || HOUSE.PASSIVE_PHASE_ONLY))
	{
		PARAMS.Fixed_SevereK_ratios = new DType[HOUSE.Num_K_Params](); 

		if (HOUSE.FSKs_Ratio_SetNum == 0) // default
		{
			PARAMS.Fixed_SevereK_ratios[0] = 0.25;
			PARAMS.Fixed_SevereK_ratios[1] = 1;
			PARAMS.Fixed_SevereK_ratios[2] = (DType)(PARAMS.Fixed_SevereK_ratios[0] / (DType)4);
		}
		else if (HOUSE.FSKs_Ratio_SetNum == 1)
		{
			PARAMS.Fixed_SevereK_ratios[0] = 0.125;
			PARAMS.Fixed_SevereK_ratios[1] = 1;
			PARAMS.Fixed_SevereK_ratios[2] = (DType)((DType)1 / (DType)16);
		}
		else if (HOUSE.FSKs_Ratio_SetNum == 2)
		{
			PARAMS.Fixed_SevereK_ratios[0] = 0.25;
			PARAMS.Fixed_SevereK_ratios[1] = 1;
			PARAMS.Fixed_SevereK_ratios[2] = 0.25;
		}
		else std::cerr << "Init Params Error: HOUSE.FSKs_Ratio_SetNum value not recognized" << endl;

		Update_Severe_Ks(PARAMS, HOUSE);
	}

	if (HOUSE.PS_Ks_Multiply_AM_Ks) //// i.e. if Passive Ks multiply Active Ks (KAM_i x KPS_i), then simply redefine. 
		for (int country = 0; country < HOUSE.TotalCountries; country++)
			for (int PrevInf = 0; PrevInf < HOUSE.Num_K_Params; PrevInf++)
				for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
					PARAMS.K_s[country][PassiveSevere][PrevInf][serotype] *= PARAMS.K_s[country][ActiveMild][PrevInf][serotype];

	DeAllocate_2D_Array(PARAMS.Initial_Ks, HOUSE.HowManyCaseCategories);  //// move this to after hospitalisation. 
}
void Initialize_AS_PrimingParams		(int &ParamCounter, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	if (HOUSE.ModelVariant == AS_PRIME)
	{
		Allocate_2D_Array(PARAMS.AS_Priming_Params, HOUSE.HowManySeroStatuses, HOUSE.Num_AS_Priming_ParamsPer_BS); 

		PARAMS.ParamNumbers.Min_AS_PrimingParam = ParamCounter;

		//// Populate AS_Priming_Params (and the rest). 
		PARAMS.AS_Priming_Params[SeroNeg][Prop_Index]				= PARAMS.Initial_AS_Prime_PropIndptAge_SNeg;
		PARAMS.ParamVec.push_back(PARAMS.AS_Priming_Params[SeroNeg][Prop_Index]				);
		PARAMS.NamesOfParameters.push_back("SNeg_AS_Prime_Prop");
		if (!HOUSE.Skip_AS_Prime_SNegProp)
			PARAMS.ParamNosYouWillFit.push_back(ParamCounter);
		ParamCounter++;

		PARAMS.AS_Priming_Params[SeroNeg][AS_Priming_Rate_index]	= PARAMS.Initial_AS_Prime_SNeg_rate;
		PARAMS.ParamVec.push_back(PARAMS.AS_Priming_Params[SeroNeg][AS_Priming_Rate_index]);
		PARAMS.NamesOfParameters.push_back("SNeg_AS_Prime_Dur");
		if (!HOUSE.Skip_AS_Prime_SNegRate)
			PARAMS.ParamNosYouWillFit.push_back(ParamCounter);
		ParamCounter++;

		PARAMS.AS_Priming_Params[SeroPos][Prop_Index]				= PARAMS.Initial_AS_Prime_PropIndptAge_SPos;
		PARAMS.ParamVec.push_back(PARAMS.AS_Priming_Params[SeroPos][Prop_Index]);
		PARAMS.NamesOfParameters.push_back("SPos_AS_Prime_Prop");
		if (!HOUSE.Skip_AS_Prime_SPosProp)
			PARAMS.ParamNosYouWillFit.push_back(ParamCounter);
		ParamCounter++;

		PARAMS.AS_Priming_Params[SeroPos][AS_Priming_Rate_index]	= PARAMS.Initial_AS_Prime_SPos_rate;
		PARAMS.ParamVec.push_back(PARAMS.AS_Priming_Params[SeroPos][AS_Priming_Rate_index]	);
		PARAMS.NamesOfParameters.push_back("SPos_AS_Prime_Dur");
		if (!HOUSE.Skip_AS_Prime_SPosRate)
			PARAMS.ParamNosYouWillFit.push_back(ParamCounter);
		ParamCounter++;

		PARAMS.ParamNumbers.Max_AS_PrimingParam = ParamCounter - 1; 

		//// K's and KPlus values
		Allocate_5D_Array(PARAMS.Ks_Prime				, HOUSE.TotalCountries, HOUSE.HowManyCaseCategories	, HOUSE.HowManySeroStatuses, HOUSE.N_STypes_Ks, HOUSE.HowManyAges); //// NOTE: HOUSE.HowManySeroStatuses and not  HOUSE.Num_K_Params as above. 
		Allocate_3D_Array(PARAMS.KPlusPrimeValues		, HOUSE.TotalCountries, HOUSE.HowManyAges			, HOUSE.HowManyCaseCategories);
		Allocate_3D_Array(PARAMS.SumRhoK0_SNeg_Primes	, HOUSE.TotalCountries, HOUSE.HowManyCaseCategories	, HOUSE.HowManyAges);
	}
}
void InitializeEfficacies				(int &ParamCounter, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	//// decide on skipping conditions and initializations for Efficacy parameters, if doing age specific vaccine efficacy.
	if (HOUSE.ASVE != Age_Option::INDEPENDENT & !HOUSE.SeroSpecificEfficacies)
	{
		//// In these scenarios, PosEfficacy and NegEfficacy variables (and the PARAMS.Efficacies array they go on to populate) won't be fitted, as your ASVE params will deal with efficacy. 
		//// But to save lots of if statements, easier to include them in the likelihood calculations etc. but simply fix them at 1 so they alter nothing. 
		if (HOUSE.ASVE_OnlyOneSeroStatus)
		{
				 if (HOUSE.ASVE_BS == SeroPos)	PARAMS.PosEfficacy = 1;
			else if (HOUSE.ASVE_BS == SeroNeg)	PARAMS.NegEfficacy = 1;
		}
		else
		{
			PARAMS.PosEfficacy = 1;
			PARAMS.NegEfficacy = 1;
		}
	}

	////// allocate memory for Efficacies 
	Allocate_3D_Array(PARAMS.Efficacies, HOUSE.HowManyCaseCategories, HOUSE.N_STypes_VEs, HOUSE.HowManySeroStatuses); //// Note that first dimension is NOT HOUSE.NumEffsPer_BS_And_Serotype. If not modelling separate efficacies for each PhaseSeverity, then HOUSE.NumEffsPer_BS_And_Serotype = 1, however still need to refer to whatever efficacy is in the PassiveSevere PhaseSeverity. 
	////// allocate memory for VacHazLikes (i.e. sum of IntVacHazards for all appropriate patients)
	Allocate_3D_Array(PARAMS.VacHazLikes, HOUSE.TotalCountries, HOUSE.HowManySeroStatuses, HOUSE.HowManyCaseCategories);

	//if (HOUSE.ResidEffs) //// allocate the memory regardless of whether you'll do HOUSE.ResidEffs, as then you can have (1 - VE_inf) = 1, so multipliers will be unaffected and you'll save on if statements. 
	Allocate_3D_Array(PARAMS.Inf_Effs, HOUSE.HowManyCaseCategories, HOUSE.N_STypes_VEs, HOUSE.HowManySeroStatuses); 
	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++) 
		for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
				PARAMS.Inf_Effs[PhaseSeverity][serotype][BaselineSeroStatus] = 0;

	string PhaseSeverityString_Dummy = ""; 
	PARAMS. ParamNumbers.Min_VacE = ParamCounter;
	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.NumEffsPer_BS_And_SType; PhaseSeverity++) //// Note: NOT PhaseSeverity < HOUSE.HowManyCaseCategories !!!!!!!!
		for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
		{
			//// Populate
			if (HOUSE.ASVE != Age_Option::INDEPENDENT & HOUSE.SeroSpecificEfficacies & (serotype == 0 || (HOUSE.Single_SNeg_Eff & HOUSE.Single_SPos_Eff))) 
			{
				//// need serotype 1 (0 in Cpp) to be the baseline. 
				//// If additive, need to add 0 to age efficacy profile to leave it unchanged. But if multiplicative, need to multiply age efficacy profile by 1 to leave it unchanged. 
				//// You ensure that serotype 1 (0 in Cpp) PARAMS.Efficacies are skipped later in this function (i.e. you don't add to PARAMS.ParamNosYouWillFit)
				DType Init_Identity_Value = HOUSE.SSASVE_Additive ? 0 : 1; 
				PARAMS.Efficacies[PhaseSeverity][serotype][SeroNeg] = Init_Identity_Value;
				PARAMS.Efficacies[PhaseSeverity][serotype][SeroPos] = Init_Identity_Value;
			}
			else
			{
				PARAMS.Efficacies[PhaseSeverity][serotype][SeroNeg] = PARAMS.NegEfficacy;
				PARAMS.Efficacies[PhaseSeverity][serotype][SeroPos] = PARAMS.PosEfficacy;
			}

			//// Make String
			if (HOUSE.PSVEs) PhaseSeverityString_Dummy = HOUSE.PhaseSeverity_strings[PhaseSeverity] + "_"; else PhaseSeverityString_Dummy = ""; 

			bool AgeContinueCondition = HOUSE.ASVE == Age_Option::INDEPENDENT || !HOUSE.SeroSpecificEfficacies || serotype != 0 || (HOUSE.ASVE_FitAllSero_VEs & HOUSE.Skip_All_ASVEs); //// Fit param if not considering age or SSVEs, or if seroype > 0 (because if fitting age and serotype effects need only skip first serotype). Last two conditions are for debugging. 

			//// SeroNeg 	//// do  SeroNeg first as helpful for SeroNeg to be zero when doing modulo arithmetic. 
			PARAMS.ParamVec.			push_back(PARAMS.Efficacies[PhaseSeverity][serotype][SeroNeg]);
			PARAMS.NamesOfParameters.	push_back(PhaseSeverityString_Dummy + "SNegEff_" + std::to_string(serotype + 1)); //// plus 1 so serotype names range from 1 to 4.
			if (!HOUSE.Skip_NegEfficacies)
				if (!HOUSE.Single_SNeg_Eff || (HOUSE.Which_SNeg_SeroFitsAll == serotype)) //// either not doing Single_SNeg_Eff, or doing Single_SNeg_Eff AND the serotype you're fitting is this one from the loop. 
					if (AgeContinueCondition ||  (HOUSE.ASVE_OnlyOneSeroStatus && HOUSE.ASVE_BS != SeroNeg))
						PARAMS.ParamNosYouWillFit.push_back(ParamCounter);
			ParamCounter++;

			//// SeroPos	//// do  SeroPos 2nd  as helpful for SeroPos to be 1 when doing modulo arithmetic. 
			PARAMS.ParamVec.			push_back(PARAMS.Efficacies[PhaseSeverity][serotype][SeroPos]);	//// 
			PARAMS.NamesOfParameters.	push_back(PhaseSeverityString_Dummy + "SPosEff_" + std::to_string(serotype + 1)); //// plus 1 so serotype names range from 1 to 4. 
			if (!HOUSE.Skip_PosEfficacies)
				if (!HOUSE.Single_SPos_Eff || (HOUSE.Which_SPos_SeroFitsAll == serotype))
					if (AgeContinueCondition || (HOUSE.ASVE_OnlyOneSeroStatus && HOUSE.ASVE_BS != SeroPos))
						PARAMS.ParamNosYouWillFit.push_back(ParamCounter);
			ParamCounter++;
		}
	PARAMS. ParamNumbers.Max_VacE = ParamCounter - 1; 


	////// ResidEffs. 
	if (HOUSE.ResidEffs)
	{
		PARAMS. ParamNumbers.Min_VacE_atInf = ParamCounter;
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.NumEffsPer_BS_And_SType; PhaseSeverity++) //// Note: NOT PhaseSeverity < HOUSE.HowManyCaseCategories !!!!!!!!
			for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
			{
				//// Populate
				PARAMS.Inf_Effs[PhaseSeverity][serotype][SeroNeg]	= 0; ///// i.e. start off with default model where efficacies decline to zero at infinity. 
				PARAMS.Inf_Effs[PhaseSeverity][serotype][SeroPos]	= 0; ///// i.e. start off with default model where efficacies decline to zero at infinity. 

				//// Make String
				if (HOUSE.PSVEs) PhaseSeverityString_Dummy = HOUSE.PhaseSeverity_strings[PhaseSeverity] + "_"; else PhaseSeverityString_Dummy = ""; 

				//// SeroNeg 	//// do  SeroNeg first as helpful for SeroNeg to be zero when doing modulo arithmetic. 
				PARAMS.ParamVec.			push_back(PARAMS.Inf_Effs[PhaseSeverity][serotype][SeroNeg]);
				PARAMS.NamesOfParameters.	push_back(PhaseSeverityString_Dummy + "SNegEff_atInf_" + std::to_string(serotype + 1)); //// plus 1 so serotype names range from 1 to 4.
				if (!HOUSE.Skip_NegEfficacies)
					if (!HOUSE.Single_SNeg_Eff || (HOUSE.Which_SNeg_SeroFitsAll == serotype)) //// either not doing Single_SNeg_Eff, or doing Single_SNeg_Eff AND the serotype you're fitting is this one from the loop. 
						PARAMS.ParamNosYouWillFit.push_back(ParamCounter);
				ParamCounter++;

				//// SeroPos	//// do  SeroPos 2nd  as helpful for SeroPos to be 1 when doing modulo arithmetic. 
				PARAMS.ParamVec.			push_back(PARAMS.Inf_Effs[PhaseSeverity][serotype][SeroPos]);	//// 
				PARAMS.NamesOfParameters.	push_back(PhaseSeverityString_Dummy + "SPosEff_atInf_" + std::to_string(serotype + 1)); //// plus 1 so serotype names range from 1 to 4. 
				if (!HOUSE.Skip_PosEfficacies)
					if (!HOUSE.Single_SPos_Eff || (HOUSE.Which_SPos_SeroFitsAll == serotype))
						PARAMS.ParamNosYouWillFit.push_back(ParamCounter);
				ParamCounter++; 
			}
		PARAMS. ParamNumbers.Max_VacE_atInf = ParamCounter - 1; 
	}

	if (!HOUSE.PSVEs && HOUSE.HowManyCaseCategories == 2) //// i.e. make both PhaseSeverities point to same value.
	{
		PARAMS.Efficacies	[PassiveSevere]	= PARAMS.Efficacies	[ActiveMild];
		PARAMS.Inf_Effs		[PassiveSevere]	= PARAMS.Inf_Effs	[ActiveMild];
	}

}
void InitializeWaningParamsAndValues	(int &ParamCounter, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, const DATA_struct &DATA)
{
	Allocate_2D_Array(PARAMS.WaningParams							, HOUSE.HowManySeroStatuses, HOUSE.NumWaningParamsPer_BS);
	Allocate_3D_Array(PARAMS.WaningMults , HOUSE.HowManyAges		, HOUSE.HowManySeroStatuses, DATA.N_WaningDays + 1); //// Note you don't use HOUSE.NumAges_Waning, as you will still want to refer to a waning value for a particular age, even if age irrelevant. 

	if (HOUSE.AS_Waning == Age_Option::INDEPENDENT) //// i.e. make all ages point to the same value if not doing age-specific waning.
		for (int age = 0; age < HOUSE.HowManyAges; age++)
			PARAMS.WaningMults[age] = PARAMS.WaningMults[2]; //// 2 because you will not calculate waning values for age 0, as minimum is 2. Shouldn't actually matter but keep anyway. 

	PARAMS.ParamNumbers.MinWaningParam = ParamCounter; 
	DType WaningParamValue; 
	string BS_String = "", WaningParamNameString = "";

	if (!HOUSE.HillWaning)
	{
		if (HOUSE.AS_Waning == Age_Option::SPLINE || HOUSE.AS_Waning == Age_Option::SPLINE_LINE || HOUSE.AS_Waning == Age_Option::SPLINE_STEP || HOUSE.AS_Waning == Age_Option::CUBIC)
		{
			Allocate_3D_Array(PARAMS.Age_SplineCoeffs_DurationRate	, HOUSE.HowManySeroStatuses, HOUSE.PolynomialsPerSpline_WaningDuration, HOUSE.MaxSplineDegree_WaningDuration + 1);
			Allocate_2D_Array(PARAMS.Age_DurationRate_xKnots		, HOUSE.HowManySeroStatuses, HOUSE.KnotsPerSpline_WaningDuration);
		}

		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.NumSeroStatuses_AS_Waning; BaselineSeroStatus++)
		{
			for (int WaningParamNo = 0; WaningParamNo < HOUSE.NumWaningParamsPer_BS; WaningParamNo++)
			{
				bool SkipParam = false;
				bool ParamIsADurationOrRate =	WaningParamNo == WaningRateDuration				|| 
												HOUSE.AS_Waning == Age_Option::CATEGORICAL		|| 
												HOUSE.AS_Waning == Age_Option::SPLINE			|| 
												HOUSE.AS_Waning == Age_Option::SPLINE_LINE		|| 
												HOUSE.AS_Waning == Age_Option::SPLINE_STEP		||
												HOUSE.AS_Waning == Age_Option::CUBIC			; 

				if (ParamIsADurationOrRate) //// in case of categorical durations/rates, or where we use splines, each value/knot will be a duration/rate, and not a halflife or power. 
				{
						 if (BaselineSeroStatus == SeroNeg)		WaningParamValue = PARAMS.NegWaning; 
					else if (BaselineSeroStatus == SeroPos)		WaningParamValue = PARAMS.PosWaning;
				}
				else if (WaningParamNo == Halflife_Index)		WaningParamValue = PARAMS.Initial_Halflife_AS_Waning;
				else if (WaningParamNo == Power_Index	)		WaningParamValue = PARAMS.Initial_Power_AS_Waning	;

				///// Reset yKnots (i.e. WaningParamValue) for splines to make them start at decent locations. 
				if (HOUSE.AS_Waning == Age_Option::SPLINE || HOUSE.AS_Waning == Age_Option::SPLINE_STEP || HOUSE.AS_Waning == Age_Option::SPLINE_LINE || HOUSE.AS_Waning == Age_Option::CUBIC)
				{
					//// Initialize from Default / Command line values. 
					if (BaselineSeroStatus == SeroNeg) WaningParamValue = PARAMS.NegWaning; else if (BaselineSeroStatus == SeroPos) WaningParamValue = PARAMS.PosWaning; 

					/*ORIG*/
					if (HOUSE.AS_Waning_KnotSet == 0) /// i.e. default knot set 
					{
						if (BaselineSeroStatus == SeroNeg & !HOUSE.Skip_NegWaning) //// i.e. don't overwrite WaningParamValue if not fitting parameter as then you'll be stuck with whatever you input here, rather than from command line. 
						{
								 if (WaningParamNo == 0)	WaningParamValue = 4;
							else if (WaningParamNo == 1)	WaningParamValue = 4;
							else if (WaningParamNo == 2)	WaningParamValue = 6;
							else if (WaningParamNo == 3)	WaningParamValue = 8;
						}
						else if (BaselineSeroStatus == SeroPos & !HOUSE.Skip_PosWaning)
						{
								 if (WaningParamNo == 0)	WaningParamValue = 9;
							else if (WaningParamNo == 1)	WaningParamValue = 10;
							else if (WaningParamNo == 2)	WaningParamValue = 12;
							else if (WaningParamNo == 3)	WaningParamValue = 15;
						}
					}
				}

				//// add to waning params array
				PARAMS.WaningParams[BaselineSeroStatus][WaningParamNo] = WaningParamValue;
				//// add to ParamVec
				PARAMS.ParamVec.push_back(WaningParamValue);

				//// add to ParamNames
						if (HOUSE.AS_Waning_Homogeneous		)	BS_String	= "";
				else	if (BaselineSeroStatus == SeroNeg	)	BS_String	= "SeroNeg";
				else	if (BaselineSeroStatus == SeroPos	)	BS_String	= "SeroPos";

				if (HOUSE.AS_Waning == Age_Option::HILL)
				{
							if (WaningParamNo == WaningRateDuration	)			WaningParamNameString	= "Waning"			;
					else	if (WaningParamNo == Halflife_Index		)			WaningParamNameString	= "WaningHalflife"	;
					else	if (WaningParamNo == Power_Index		)			WaningParamNameString	= "WaningPower"		;
				}
				else if (HOUSE.AS_Waning == Age_Option::CATEGORICAL	)	WaningParamNameString	= "Waning_AG_"		+ std::to_string(WaningParamNo + 1);
				else if (HOUSE.AS_Waning == Age_Option::SPLINE		)	WaningParamNameString	= "Waning_AgeKnot_" + std::to_string(WaningParamNo + 1); 
				else if (HOUSE.AS_Waning == Age_Option::SPLINE_LINE	)	WaningParamNameString	= "Waning_AgeKnot_" + std::to_string(WaningParamNo + 1); 
				else if (HOUSE.AS_Waning == Age_Option::SPLINE_STEP	)	WaningParamNameString	= "Waning_AgeKnot_" + std::to_string(WaningParamNo + 1); 
				else if (HOUSE.AS_Waning == Age_Option::CUBIC		)	WaningParamNameString	= "Waning_AgeKnot_" + std::to_string(WaningParamNo + 1); 
				else if (HOUSE.AS_Waning == Age_Option::INDEPENDENT	)	WaningParamNameString	= "Waning";
				PARAMS.NamesOfParameters.push_back(BS_String + WaningParamNameString);	
			
				//// decide on skipping conditions (update when necessary) and add to ParamNosYouWillFit if appropriate.
				if (BaselineSeroStatus == SeroNeg)
				{
						 if (ParamIsADurationOrRate 			&& (HOUSE.Skip_NegWaning	|| (HOUSE.AgeEffectsSame_Waning & WaningParamNo != 0))) SkipParam = true;	/// either you're skipping all Neg Waning, or your fixing them to be the same as each other (controlled in FindSimulataneousUpdateParams function), in which case you need only the first one. 
					else if (WaningParamNo == Halflife_Index	&& HOUSE.Skip_NegWaning_HLife	) SkipParam = true;
					else if (WaningParamNo == Power_Index		&& HOUSE.Skip_NegWaning_Power	) SkipParam = true;
				}
				else if (BaselineSeroStatus == SeroPos)
				{
						 if (ParamIsADurationOrRate 			&& (HOUSE.Skip_PosWaning	|| (HOUSE.AgeEffectsSame_Waning & WaningParamNo != 0))) SkipParam = true; /// either you're skipping all Pos Waning, or your fixing them to be the same as each other (controlled in FindSimulataneousUpdateParams function), in which case you need only the first one. 
					else if (WaningParamNo == Halflife_Index	&& HOUSE.Skip_PosWaning_HLife	) SkipParam = true;
					else if (WaningParamNo == Power_Index		&& HOUSE.Skip_PosWaning_Power	) SkipParam = true;
				}
				if (!SkipParam)	PARAMS.ParamNosYouWillFit.push_back(ParamCounter);

				//// increment param counter
				ParamCounter++;
			}

			//// Having populated yKnots (here WaningParams) Populate xKnots and calculate coefficients, if doing Age_Option::SPLINE
			if (HOUSE.AS_Waning == Age_Option::SPLINE || HOUSE.AS_Waning == Age_Option::SPLINE_STEP || HOUSE.AS_Waning == Age_Option::SPLINE_LINE || HOUSE.AS_Waning == Age_Option::CUBIC)
			{
				if (HOUSE.AS_Waning_KnotSet == 0)
				{
					PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus][0] = 2.0;
					PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus][1] = 6.0;
					PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus][2] = 12.0;
					if (HOUSE.AS_Waning != Age_Option::SPLINE_STEP) PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus][3] = 16.0;
				}
				else if (HOUSE.AS_Waning_KnotSet == 1)
				{
					PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus][0] = 2.0;
					PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus][1] = 5.9;
					PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus][2] = 6.1;
					PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus][3] = 11.9;
					PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus][4] = 12.1;
					if (HOUSE.AS_Waning != Age_Option::SPLINE_STEP) PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus][5] = 16.0;
				}
				else if (HOUSE.AS_Waning_KnotSet == 2)
				{
					PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus][0] = 2.0;
					PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus][1] = 4.8;
					PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus][2] = 7.6;
					PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus][3] = 10.4;
					PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus][4] = 13.2;
					if (HOUSE.AS_Waning != Age_Option::SPLINE_STEP) PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus][5] = 16.0;
				}
				else if (HOUSE.AS_Waning_KnotSet == 3)
				{
					PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus][0] = 2.0;
					PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus][1] = 5.5;
					PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus][2] = 9.0;
					PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus][3] = 12.5;
					if (HOUSE.AS_Waning != Age_Option::SPLINE_STEP) PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus][4] = 16.0;
				}
				else std::cerr << "InitializeWaningParamsAndValues ERROR: AS_Waning_KnotSet not recognized\nAS_Waning_KnotSet not recognized\nAS_Waning_KnotSet not recognized\nAS_Waning_KnotSet not recognized\n"; 

				//// Calculate coefficients
				if (HOUSE.AS_Waning == Age_Option::SPLINE)
					for (int poly = 0; poly < HOUSE.PolynomialsPerSpline_WaningDuration; poly++)
						CoefficientsFromKnots(PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus], PARAMS.WaningParams[BaselineSeroStatus], poly, PARAMS.Age_SplineCoeffs_DurationRate[BaselineSeroStatus][poly]);
				else if (HOUSE.AS_Waning == Age_Option::SPLINE_LINE)
					for (int poly = 0; poly < HOUSE.PolynomialsPerSpline_WaningDuration; poly++)
						CoefficientsFromKnots_Line(PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus], PARAMS.WaningParams[BaselineSeroStatus], poly, PARAMS.Age_SplineCoeffs_DurationRate[BaselineSeroStatus][poly]);
				else if (HOUSE.AS_Waning == Age_Option::CUBIC)
					for (int poly = 0; poly < HOUSE.PolynomialsPerSpline_WaningDuration; poly++)
						CoefficientsFromKnots_Cubic(PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus], PARAMS.WaningParams[BaselineSeroStatus], poly, PARAMS.Age_SplineCoeffs_DurationRate[BaselineSeroStatus][poly]);
			}
		}
	}
	else
	{
		string SeroStatusName = "";

		////// Half-lifes
		PARAMS. Hill_Halflives					= new DType[HOUSE.HowManySeroStatuses]();
		PARAMS. ParamNumbers.Min_HillHalfLife	= ParamCounter;
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
		{
			SeroStatusName = (BaselineSeroStatus == 0) ? "SNeg_" : "SPos_"; 

			PARAMS. Hill_Halflives[BaselineSeroStatus] = PARAMS.Initial_Halflife_Hill;
			
			PARAMS.ParamVec.			push_back(PARAMS.Hill_Halflives	[BaselineSeroStatus]);
			PARAMS.NamesOfParameters.push_back(SeroStatusName + "HillHalfLife");

			if (!(HOUSE.Skip_NegWaning && BaselineSeroStatus == 0)) //// i.e. include in fitted parameters if not skipping NegWaning or not correct serostatus. 
				if (!(HOUSE.Skip_PosWaning && BaselineSeroStatus == 1))
					PARAMS.ParamNosYouWillFit.push_back(ParamCounter);
			ParamCounter++;
		}
		PARAMS. ParamNumbers.Max_HillHalfLife = ParamCounter - 1;

		////// Powers	
		PARAMS.Hill_Powers					= new DType[HOUSE.HowManySeroStatuses]();
		PARAMS.ParamNumbers.Min_HillPower	= ParamCounter;
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
		{
			SeroStatusName = (BaselineSeroStatus == 0) ? "SNeg_" : "SPos_"; 

			PARAMS.Hill_Powers[BaselineSeroStatus] = PARAMS.Initial_Power_Hill;

			PARAMS.ParamVec.			push_back(PARAMS.Hill_Powers[BaselineSeroStatus]);
			PARAMS.NamesOfParameters.push_back(SeroStatusName + "HillPower");

			if (!(HOUSE.Skip_NegWaning && BaselineSeroStatus == 0)) //// i.e. include in fitted parameters if not skipping Negwaning or not correct serostatus. 
				if (!(HOUSE.Skip_PosWaning && BaselineSeroStatus == 1))
					PARAMS.ParamNosYouWillFit.push_back(ParamCounter);
			ParamCounter++;
		}
		PARAMS. ParamNumbers.Max_HillPower = ParamCounter - 1;
	}
	PARAMS.ParamNumbers.MaxWaningParam = ParamCounter - 1;

	//// if AS_Waning_Homogeneous  assign seropositive pointer of WaningParams and WaningMults to seronegative (for each age). Both serostatuses point to same values (for each age). 
	if (HOUSE.AS_Waning_Homogeneous)
	{
		PARAMS.WaningParams[SeroPos] = PARAMS.WaningParams[SeroNeg];

		for (int age = 0; age < HOUSE.NumAges_Waning; age++)
			PARAMS.WaningMults[age][SeroPos] = PARAMS.WaningMults[age][SeroNeg];
	}

	//// Calc the WANING VALUES using the parameters set up above. 
	for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.NumSeroStatuses_AS_Waning; BaselineSeroStatus++)
		Calc_WaningValues(PARAMS, HOUSE, BaselineSeroStatus, DATA.N_WaningDays);
}
void InitializeASVEParams				(int &ParamCounter, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	Allocate_2D_Array(PARAMS.AgeEff_Mult, HOUSE.HowManySeroStatuses, HOUSE.HowManyAges); 

	//// Initialize AgeEff_Mult's to 1. Need to do this firstly if HOUSE.ASVE == Age_Option::INDEPENDENT, but also if you're only considering a single serostatus etc. 
	for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
		for (int age = 0; age < HOUSE.HowManyAges; age++)
			PARAMS.AgeEff_Mult[BaselineSeroStatus][age] = 1;

	if (HOUSE.ASVE != Age_Option::INDEPENDENT)
	{
		////// Populate parameter array
		Allocate_2D_Array(PARAMS.ASVE_Params, HOUSE.HowManySeroStatuses, HOUSE.NumASVE_Function_Params); 

		if (HOUSE.ASVE == Age_Option::SPLINE || HOUSE.ASVE == Age_Option::SPLINE_LINE || HOUSE.ASVE == Age_Option::SPLINE_STEP || HOUSE.ASVE == Age_Option::CUBIC)
		{
			Allocate_3D_Array(PARAMS.Age_SplineCoeffs_Effs	, HOUSE.HowManySeroStatuses, HOUSE.PolynomialsPerSpline_EffMultiplier, HOUSE.MaxSplineDegree_EffMultiplier + 1);
			Allocate_2D_Array(PARAMS.Age_Effs_xKnots		, HOUSE.HowManySeroStatuses, HOUSE.KnotsPerSpline_EffMultiplier);
		}

		string SeroStatusName;
		PARAMS.ParamNumbers.Min_ASVE_Param = ParamCounter;

		int Start_BaselineSeroStatus, End_BaselineSeroStatus;

		//// if (HOUSE.ASVE_OnlyOneSeroStatus) need only add a set of ASVE parameters for that serostatus, and also HOUSE.AS_VE_Homogeneous has been set to false in InitializeHousekeeping function. 
		//// if (!HOUSE.ASVE_OnlyOneSeroStatus)  either  apply one hill function to both serostatuses (i.e. HOUSE.AS_VE_Homogeneous == true), or  apply one to each serostatus (i.e. HOUSE.AS_VE_Homogeneous == false). 
		//// if the former, then HOUSE.NumSeroStatuses_ASVEs = 1, so only one set of parameters added. If the latter, then NumSeroStatuses_ASVEs = HOUSE.HowManySeroStatuses.
		//// if/else statement below captures this. 
		if (HOUSE.ASVE_OnlyOneSeroStatus)	{ Start_BaselineSeroStatus = HOUSE.ASVE_BS	; End_BaselineSeroStatus = HOUSE.ASVE_BS + 1;}
		else								{ Start_BaselineSeroStatus = 0				; End_BaselineSeroStatus = HOUSE.NumSeroStatuses_ASVEs;}	

		for (int BaselineSeroStatus = Start_BaselineSeroStatus; BaselineSeroStatus < End_BaselineSeroStatus; BaselineSeroStatus++)
		{
			if (HOUSE.ASVE == Age_Option::HILL)
			{
					 if (HOUSE.AS_VE_Homogeneous)			SeroStatusName = "";
				else if (BaselineSeroStatus == SeroNeg)		SeroStatusName = "SNeg_";
				else if (BaselineSeroStatus == SeroPos)		SeroStatusName = "SPos_";

				for (int Type = 0; Type < HOUSE.NumASVE_Function_Params; Type++)
				{
					if (Type == Halflife_Index) //// doing it this way with the various if statements should render hash define definitions arbitrary. 
					{
						PARAMS.								ASVE_Params[BaselineSeroStatus][Halflife_Index]	= PARAMS.Initial_Halflife_ASVE;
						PARAMS.ParamVec.push_back(PARAMS.	ASVE_Params[BaselineSeroStatus][Halflife_Index]);
						PARAMS.NamesOfParameters.push_back(SeroStatusName + "ASVE_HalfLife");

						if (!HOUSE.Skip_All_ASVEs)
							if (!HOUSE.Skip_ASVE_Halflives)
								PARAMS.ParamNosYouWillFit.push_back(ParamCounter);
					}
					else if (Type == Power_Index)
					{
						PARAMS.								ASVE_Params[BaselineSeroStatus][Power_Index]	= PARAMS.Initial_Power_ASVE;
						PARAMS.ParamVec.push_back(PARAMS.	ASVE_Params[BaselineSeroStatus][Power_Index]);
						PARAMS.NamesOfParameters.push_back(SeroStatusName + "ASVE_Power"	);

						if (!HOUSE.Skip_All_ASVEs)
							if (!HOUSE.Skip_ASVEsPowers)
								PARAMS.ParamNosYouWillFit.push_back(ParamCounter);
					}
					else if (Type == Prop_Index)
					{
						PARAMS.								ASVE_Params[BaselineSeroStatus][Prop_Index]		= PARAMS.Initial_Prop_ASVE;
						PARAMS.ParamVec.push_back(PARAMS.	ASVE_Params[BaselineSeroStatus][Prop_Index]);
						PARAMS.NamesOfParameters.push_back(SeroStatusName + "ASVE_Prop"	);

						if (!HOUSE.Skip_All_ASVEs)
							if (!HOUSE.Skip_ASVE_Props)
								PARAMS.ParamNosYouWillFit.push_back(ParamCounter);
					}
					ParamCounter++;
				}
				
				//// now Calc_AgeVaccEffMults
				Calc_AgeVaccEffMults(PARAMS, HOUSE, BaselineSeroStatus);
			}
			else 
			{
					 if (HOUSE.AS_VE_Homogeneous)			SeroStatusName = "";
				else if (BaselineSeroStatus == SeroNeg)		SeroStatusName = "SNeg";
				else if (BaselineSeroStatus == SeroPos)		SeroStatusName = "SPos";

				for (int AgeEffParam = 0; AgeEffParam < HOUSE.NumASVE_Function_Params; AgeEffParam++)
				{
					if (BaselineSeroStatus == SeroNeg)	PARAMS.ASVE_Params[BaselineSeroStatus][AgeEffParam] = 0.27;
					if (BaselineSeroStatus == SeroPos)	
						if (HOUSE.SeroSpecificEfficacies & HOUSE.SSASVE_Additive) 
							PARAMS.ASVE_Params[BaselineSeroStatus][AgeEffParam] = 0.3;
						else
							PARAMS.ASVE_Params[BaselineSeroStatus][AgeEffParam] = 0.5;

					///// Reset yKnots (i.e. AgeEffParam) for splines to make them start at decent locations. 
					if (HOUSE.SeroSpecificEfficacies & HOUSE.ASVE_FitAllSero_VEs & HOUSE.Skip_All_ASVEs)	PARAMS.ASVE_Params[BaselineSeroStatus][AgeEffParam] = 1; 
					else if (HOUSE.ASVE == Age_Option::SPLINE || HOUSE.ASVE == Age_Option::SPLINE_STEP || HOUSE.ASVE == Age_Option::SPLINE_LINE || HOUSE.ASVE == Age_Option::CUBIC)
					{
						if (!HOUSE.ASVE_AdditionalKnots & !HOUSE.SSASVE_Additive) /// i.e. default knot set 
						{
							/*ORIG*/
							if (BaselineSeroStatus == SeroNeg)			//// i.e. default knot set 
							{
									 if (AgeEffParam == 0)	PARAMS.ASVE_Params[BaselineSeroStatus][AgeEffParam] = -0.4;
								else if (AgeEffParam == 1)	PARAMS.ASVE_Params[BaselineSeroStatus][AgeEffParam] = 0.4;
								else if (AgeEffParam == 2)	PARAMS.ASVE_Params[BaselineSeroStatus][AgeEffParam] = 0.35;
								else if (AgeEffParam == 3)	PARAMS.ASVE_Params[BaselineSeroStatus][AgeEffParam] = 0.6;
							}
							else if (BaselineSeroStatus == SeroPos)
							{
									 if (AgeEffParam == 0)	PARAMS.ASVE_Params[BaselineSeroStatus][AgeEffParam] = 0.2;
								else if (AgeEffParam == 1)	PARAMS.ASVE_Params[BaselineSeroStatus][AgeEffParam] = 0.4;
								else if (AgeEffParam == 2)	PARAMS.ASVE_Params[BaselineSeroStatus][AgeEffParam] = 0.6;
								else if (AgeEffParam == 3)	PARAMS.ASVE_Params[BaselineSeroStatus][AgeEffParam] = 0.8;
							}
						}
					}

					PARAMS.ParamVec.push_back(PARAMS.ASVE_Params[BaselineSeroStatus][AgeEffParam]);

					if (HOUSE.ASVE == Age_Option::CATEGORICAL)
						PARAMS.NamesOfParameters.push_back(SeroStatusName + "Eff_AgeGroup_" + std::to_string(AgeEffParam));
					else if (HOUSE.ASVE == Age_Option::SPLINE || HOUSE.ASVE == Age_Option::SPLINE_LINE || HOUSE.ASVE == Age_Option::SPLINE_STEP || HOUSE.ASVE == Age_Option::CUBIC)
						PARAMS.NamesOfParameters.push_back(SeroStatusName + "Eff_AgeKnot_" + std::to_string(AgeEffParam));

					int AgeEffectsSameDummyParamNo = HOUSE.ASVE == Age_Option::CATEGORICAL ? 1 : 0; 

					bool SkipParam = HOUSE.Skip_All_ASVEs || (HOUSE.Skip_SNeg_ASVEs & BaselineSeroStatus == SeroNeg) || (HOUSE.Skip_SPos_ASVEs & BaselineSeroStatus == SeroPos);

					if (AgeEffParam > 0 || HOUSE.ASVE == Age_Option::SPLINE || HOUSE.ASVE == Age_Option::SPLINE_LINE || HOUSE.ASVE == Age_Option::SPLINE_STEP || HOUSE.ASVE == Age_Option::CUBIC) 
						if (!SkipParam)
							if (!HOUSE.AgeEffectsSame_VE || AgeEffParam == AgeEffectsSameDummyParamNo) 
								PARAMS.ParamNosYouWillFit.push_back(ParamCounter);
					ParamCounter++;
				}
				
				///// Populate xKnots and manually input their equidistant values between 2 and 16. Here, ASVE_Params[BaselineSeroStatus] are the yKnots.  
				if (HOUSE.ASVE == Age_Option::SPLINE || HOUSE.ASVE == Age_Option::SPLINE_LINE || HOUSE.ASVE == Age_Option::SPLINE_STEP || HOUSE.ASVE == Age_Option::CUBIC)
				{
					if (HOUSE.ASVE_AdditionalKnots)
					{
						PARAMS.Age_Effs_xKnots[BaselineSeroStatus][0] = 2.0;
						PARAMS.Age_Effs_xKnots[BaselineSeroStatus][1] = 5.9;
						PARAMS.Age_Effs_xKnots[BaselineSeroStatus][2] = 6.1;
						PARAMS.Age_Effs_xKnots[BaselineSeroStatus][3] = 11.9;
						PARAMS.Age_Effs_xKnots[BaselineSeroStatus][4] = 12.1;
						if (HOUSE.ASVE != Age_Option::SPLINE_STEP) PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus][5] = 16.0;
					}
					else
					{
						PARAMS.Age_Effs_xKnots[BaselineSeroStatus][0] = 2.0;
						PARAMS.Age_Effs_xKnots[BaselineSeroStatus][1] = 6.0;
						PARAMS.Age_Effs_xKnots[BaselineSeroStatus][2] = 12.0;
						if (HOUSE.ASVE != Age_Option::SPLINE_STEP) PARAMS.Age_Effs_xKnots[BaselineSeroStatus][3] = 16.0;
					}

					if (HOUSE.ASVE == Age_Option::SPLINE)
						for (int poly = 0; poly < HOUSE.PolynomialsPerSpline_EffMultiplier; poly++)
							CoefficientsFromKnots(PARAMS.Age_Effs_xKnots[BaselineSeroStatus], PARAMS.ASVE_Params[BaselineSeroStatus], poly, PARAMS.Age_SplineCoeffs_Effs[BaselineSeroStatus][poly]);
					else if (HOUSE.ASVE == Age_Option::SPLINE_LINE)
						for (int poly = 0; poly < HOUSE.PolynomialsPerSpline_EffMultiplier; poly++)
							CoefficientsFromKnots_Line(PARAMS.Age_Effs_xKnots[BaselineSeroStatus], PARAMS.ASVE_Params[BaselineSeroStatus], poly, PARAMS.Age_SplineCoeffs_Effs[BaselineSeroStatus][poly]);
					else if (HOUSE.ASVE == Age_Option::CUBIC)
						for (int poly = 0; poly < HOUSE.PolynomialsPerSpline_EffMultiplier; poly++)
							CoefficientsFromKnots_Cubic(PARAMS.Age_Effs_xKnots[BaselineSeroStatus], PARAMS.ASVE_Params[BaselineSeroStatus], poly, PARAMS.Age_SplineCoeffs_Effs[BaselineSeroStatus][poly]);
				}
				
				//// once all parameters added (and coefficients calculated if doing spline) to arrays/vectors etc.  Calc_AgeVaccEffMults
				Calc_AgeVaccEffMults(PARAMS, HOUSE, BaselineSeroStatus);
			}
		}
		PARAMS.ParamNumbers.Max_ASVE_Param = ParamCounter - 1;

		//// if AS_VE_Homogeneous then assign seropositive pointer of ASVE_Params and AgeEff_Mult to seronegative. Both serostatuses point to same values. 
		if (HOUSE.AS_VE_Homogeneous)
		{
			PARAMS.ASVE_Params[SeroPos] = PARAMS.ASVE_Params[SeroNeg];
			PARAMS.AgeEff_Mult[SeroPos] = PARAMS.AgeEff_Mult[SeroNeg];
		}

		if (HOUSE.ASVE_OnlyOneSeroStatus) 
		{
			//// change AgeEff_Mult
			for (int age = 0; age < HOUSE.HowManyAges; age++)
				if (!HOUSE.AdditiveSSASVEs)
					PARAMS.AgeEff_Mult[1 - HOUSE.ASVE_BS][age] = 1;		//// note the 1 - HOUSE.ASVE_BS. Won't work if serostatus non-binary. 
				else 
					PARAMS.AgeEff_Mult[1 - HOUSE.ASVE_BS][age] = 0;		//// note the 1 - HOUSE.ASVE_BS. Won't work if serostatus non-binary. 
			//// Set proportion of VE affected by age to zero for other baselineserostatus (i.e. not HOUSE.ASVE_BS) 
			if (HOUSE.ASVE == Age_Option::HILL)
				PARAMS.ASVE_Params[1 - HOUSE.ASVE_BS][Prop_Index] = 0;	//// note the 1 - HOUSE.ASVE_BS. Won't work if serostatus non-binary. 
		}
		if (HOUSE.HowManySeroStatuses > 2)
			for (int blah = 0; blah < 80000; blah++) std::cerr << "ASVEs not coded for HOUSE.HowManySeroStatuses > 2" << endl; 
	}
}
void Initialize_AS_Haz_Params			(int &ParamCounter, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	PARAMS.AgeHaz_Mult = new DType[HOUSE.HowManyAges]();
	for (int age = 0; age < HOUSE.HowManyAges; age++)	PARAMS.AgeHaz_Mult[age] = 1;

	if (HOUSE.AS_Haz == Age_Option::SPLINE || HOUSE.AS_Haz == Age_Option::SPLINE_LINE || HOUSE.AS_Haz == Age_Option::SPLINE_STEP || HOUSE.AS_Haz == Age_Option::CUBIC)
	{
		Allocate_2D_Array(PARAMS.Age_SplineCoeffs_HazMult, HOUSE.PolynomialsPerSpline_HazMultiplier, HOUSE.MaxSplineDegree_HazMultiplier + 1);
		PARAMS.Age_HazMult_xKnots = new DType[HOUSE.KnotsPerSpline_HazMultiplier]();
	}

	if (HOUSE.AS_Haz != Age_Option::INDEPENDENT)
	{
		////// Populate parameter array

		PARAMS.ASHaz_Params = new DType[HOUSE.NumAS_Haz_Function_Params]();
		PARAMS.ParamNumbers.Min_ASHaz_Param = ParamCounter;

			 if (HOUSE.AS_Haz == Age_Option::HILL)
		{
			for (int Type = 0; Type < HOUSE.NumAS_Haz_Function_Params; Type++)
			{
				if (Type == Halflife_Index) //// doing it this way with the various if statements should render hash define definitions arbitrary. 
				{
					PARAMS.ASHaz_Params[Halflife_Index]	= PARAMS.Initial_Halflife_AS_Haz;
					PARAMS.ParamVec.push_back(PARAMS.ASHaz_Params[Halflife_Index]);
					PARAMS.NamesOfParameters.push_back("AS_Haz_Halflife");

					if (!HOUSE.Skip_All_AS_Haz)
						if (!HOUSE.Skip_AS_Haz_Halflives)
							PARAMS.ParamNosYouWillFit.push_back(ParamCounter);
				}
				else if (Type == Power_Index)
				{
					PARAMS.ASHaz_Params[Power_Index]	= PARAMS.Initial_Power_AS_Haz;
					PARAMS.ParamVec.push_back(PARAMS.ASHaz_Params[Power_Index]);
					PARAMS.NamesOfParameters.push_back("AS_Haz_Power"	);

					if (!HOUSE.Skip_All_AS_Haz)
						if (!HOUSE.Skip_AS_Haz_Powers)
							PARAMS.ParamNosYouWillFit.push_back(ParamCounter);
				}
				else if (Type == Prop_Index)
				{
					PARAMS.ASHaz_Params[Prop_Index] = PARAMS.Initial_Prop_AS_Haz;
					PARAMS.ParamVec.push_back(PARAMS.ASHaz_Params[Prop_Index]);
					PARAMS.NamesOfParameters.push_back("AS_Haz_Prop");

					if (!HOUSE.Skip_All_AS_Haz)
						if (!HOUSE.Skip_AS_Haz_Props)
							PARAMS.ParamNosYouWillFit.push_back(ParamCounter);
				}
				ParamCounter++;
			}
		}
		else if (HOUSE.AS_Haz == Age_Option::CATEGORICAL)
		{
			for (int AS_Haz_AgeMult = 0; AS_Haz_AgeMult < HOUSE.NumAS_Haz_Function_Params; AS_Haz_AgeMult++)
			{
				PARAMS.ASHaz_Params[AS_Haz_AgeMult] = 1;
				PARAMS.ParamVec.push_back(PARAMS.ASHaz_Params[AS_Haz_AgeMult]);
				PARAMS.NamesOfParameters.push_back("AS_Haz_AgeGroupMult_" + std::to_string(AS_Haz_AgeMult));

				if (AS_Haz_AgeMult > 1) //// > 1 as NumAS_Haz_Function_Params = 4 for categorical (because of Age groups 0,1,2 and 3). Age group 0 is everyone, 1 is 2-5, 2 is 6-11, 3 is 12-14. Need to set one of the age groups (not everyone) as baseline, so set Age group 1, 2-5yrs, as baseline)
					if (!HOUSE.Skip_All_AS_Haz)
						PARAMS.ParamNosYouWillFit.push_back(ParamCounter);
				ParamCounter++;
			}
		}
		else if (HOUSE.AS_Haz == Age_Option::SPLINE || HOUSE.AS_Haz == Age_Option::SPLINE_LINE || HOUSE.AS_Haz == Age_Option::SPLINE_STEP || HOUSE.AS_Haz == Age_Option::CUBIC)
		{
			for (int AS_Haz_Knot = 0; AS_Haz_Knot < HOUSE.NumAS_Haz_Function_Params; AS_Haz_Knot++)
			{
				PARAMS.ASHaz_Params[AS_Haz_Knot] = 1;
				if (AS_Haz_Knot == 1) PARAMS.ASHaz_Params[AS_Haz_Knot] = 1.549;
				if (AS_Haz_Knot == 2) PARAMS.ASHaz_Params[AS_Haz_Knot] = 1.567;
				if (AS_Haz_Knot == 3) PARAMS.ASHaz_Params[AS_Haz_Knot] = 0.8775;


				PARAMS.ParamVec.push_back(PARAMS.ASHaz_Params[AS_Haz_Knot]);
				PARAMS.NamesOfParameters.push_back("AS_Haz_Knot_" + std::to_string(AS_Haz_Knot));

				if (AS_Haz_Knot > 0) /// 
					if (!HOUSE.Skip_All_AS_Haz)
						PARAMS.ParamNosYouWillFit.push_back(ParamCounter);
				ParamCounter++;
			}

			/// Populate xKnots and manually input their equidistant values between 2 and 16. Here, ASVE_Params[BaselineSeroStatus] are the yKnots.  
			if (HOUSE.AS_Haz == Age_Option::SPLINE || HOUSE.AS_Haz == Age_Option::SPLINE_LINE || HOUSE.AS_Haz == Age_Option::SPLINE_STEP || HOUSE.AS_Haz == Age_Option::CUBIC)
			{
				if (HOUSE.AS_Haz_AdditionalKnots)
				{
					PARAMS.Age_HazMult_xKnots[0] = 2.0;
					PARAMS.Age_HazMult_xKnots[1] = 5.9;
					PARAMS.Age_HazMult_xKnots[2] = 6.1;
					PARAMS.Age_HazMult_xKnots[3] = 11.9;
					PARAMS.Age_HazMult_xKnots[4] = 12.1;
					if (HOUSE.AS_Haz != Age_Option::SPLINE_STEP) PARAMS.Age_HazMult_xKnots[5] = 16.0;
				}
				else
				{
					PARAMS.Age_HazMult_xKnots[0] = 2.0;
					PARAMS.Age_HazMult_xKnots[1] = 6.0;
					PARAMS.Age_HazMult_xKnots[2] = 12.0;
					if (HOUSE.AS_Haz != Age_Option::SPLINE_STEP) PARAMS.Age_HazMult_xKnots[3] = 16.0;
				}

				if (HOUSE.AS_Haz == Age_Option::SPLINE)
					for (int poly = 0; poly < HOUSE.PolynomialsPerSpline_EffMultiplier; poly++)
						CoefficientsFromKnots(PARAMS.Age_HazMult_xKnots, PARAMS.ASHaz_Params, poly, PARAMS.Age_SplineCoeffs_HazMult[poly]);
				else if (HOUSE.AS_Haz == Age_Option::SPLINE_LINE)
					for (int poly = 0; poly < HOUSE.PolynomialsPerSpline_EffMultiplier; poly++)
						CoefficientsFromKnots_Line(PARAMS.Age_HazMult_xKnots, PARAMS.ASHaz_Params, poly, PARAMS.Age_SplineCoeffs_HazMult[poly]);
				else if (HOUSE.AS_Haz == Age_Option::CUBIC)
					for (int poly = 0; poly < HOUSE.PolynomialsPerSpline_EffMultiplier; poly++)
						CoefficientsFromKnots_Cubic(PARAMS.Age_HazMult_xKnots, PARAMS.ASHaz_Params, poly, PARAMS.Age_SplineCoeffs_HazMult[poly]);
			}
		}
		PARAMS.ParamNumbers.Max_ASHaz_Param = ParamCounter - 1;

		//// now Calc_AgeVaccEffMults
		Calc_Age_HazMults(PARAMS, HOUSE);
	}
}
void Initialize_BS_BaseHazMults			(int &ParamCounter, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	if (HOUSE.AdjHaz)
	{
		//// allocate
		PARAMS.BS_BaseHazMults = new DType[HOUSE.HowManySeroStatuses](); 

		//// populate
		PARAMS.BS_BaseHazMults[HOUSE.Which_BS_BaseHazMult		] = PARAMS.Initial_BS_BaseHazMult;
		PARAMS.BS_BaseHazMults[1 - HOUSE.Which_BS_BaseHazMult	] = 1; //// note the 1 - HOUSE.Which_BS_BaseHazMult. Won't work if serostatus non-binary. 


		if (!HOUSE.Skip_BS_BaseHazMult)
		{
			string SeroStatusName = (HOUSE.Which_BS_BaseHazMult == SeroPos) ? "SeroPos" : "SeroNeg";
			PARAMS.ParamVec.push_back(PARAMS.BS_BaseHazMults[HOUSE.Which_BS_BaseHazMult]);
			PARAMS.NamesOfParameters.push_back(SeroStatusName + "BaseHazMult");

			PARAMS.ParamNosYouWillFit.push_back(ParamCounter);

			//// record param number. 
			PARAMS.ParamNumbers.BS_BaseHazMult = ParamCounter;

			//// increment
			ParamCounter++;
		}
		
		if (HOUSE.HowManySeroStatuses > 2)
			for (int blah = 0; blah < 80000; blah++) 
				std::cerr << "AdjHaz not coded for HOUSE.HowManySeroStatuses > 2" << endl;
	}
}
void Initialize_Rhos					(int &ParamCounter, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	////// allocate memory for rhos 
	Allocate_2D_Array(PARAMS.rhos, HOUSE.TotalCountries, HOUSE.N_STypes);

	if (HOUSE.SeroSpecificEfficacies || HOUSE.SeroSpecific_K_values) 
	{
		if (!HOUSE.BaselinePartition)
		{
			PARAMS.ParamNumbers.Min_qval = ParamCounter;
			for (int country = 0; country < HOUSE.TotalCountries; country++)
			{
				//// first do ParamVec, then do rhos. 
				for (int serotype = 0; serotype < HOUSE.N_STypes - 1; serotype++) //// note the minus 1. Proportions of 4 serotypes uniquely determined by 3 parameters. Also no need to do serotype plus 1 for names now that these don't explicitly refer to proportions. 
				{
					//// populate ParamVec and NamesOfParameters
					PARAMS.ParamVec.			push_back( (DType (1) / DType (HOUSE.N_STypes - serotype))); //// first qval will be 1/4, second = 1/3, third = 1/2. Works out that then all proportions equal. 
					PARAMS.NamesOfParameters.	push_back("qval_" + std::to_string(country) + "_" + std::to_string(serotype));

					if (!HOUSE.Skip_qvals) //// if not skipping qvalues
						if (std::any_of(HOUSE.WhichCountries.begin(), HOUSE.WhichCountries.end(), [&](int i) {return i == country; })) //// if country included in WhichCountries
							PARAMS.ParamNosYouWillFit.push_back(ParamCounter);
					ParamCounter++;
				}
				///// Calc rhos from qvalues in ParamVec
				Calc_RhosFrom_qParams(country, PARAMS, HOUSE);
			}
			PARAMS. ParamNumbers.Max_qval = ParamCounter - 1;
		}
		else
		{
			PARAMS. SumRhos = new DType[HOUSE.TotalCountries]();
			PARAMS. ParamNumbers.Min_Rho = ParamCounter;
			for (int country = 0; country < HOUSE.TotalCountries; country++)
			{
				//// first do ParamVec, then do the rhos. 
				for (int serotype = 0; serotype < HOUSE.N_STypes; serotype++) //// unlike with qval partition, here loop goes over all serotypes (just that baseline one isn't fitted). 
				{
					//// populate ParamVec and NamesOfParameters
					PARAMS.ParamVec.			push_back(DType(1)	); //// all rhos equal to one. Will also change interpretation of basline hazard (would now expect knots to be smaller). 
					PARAMS.NamesOfParameters.push_back("rho_" + std::to_string(country) + "_" + std::to_string(serotype + 1));

					PARAMS. rhos[country][serotype] = DType(1);

					if (!HOUSE.Skip_Rhos) //// if not skipping qvalues
						if (std::any_of(HOUSE.WhichCountries.begin(), HOUSE.WhichCountries.end(), [&](int i) {return i == country; })) //// if country included in WhichCountries
							if (!(serotype == HOUSE.Baseline_Serotype))	PARAMS.ParamNosYouWillFit.push_back(ParamCounter); //// keep all of them as parameters (will make various "Find blah from blah" and "IsParamBlah" functions easier, but jsut don't fit the baseline serotype proportion, e.g. p1 = 1. 
					ParamCounter++;
				}
			}
			Calc_SumRhos_ALL(PARAMS, HOUSE);
			PARAMS. ParamNumbers.Max_Rho = ParamCounter - 1;
		}
	}
	else for (int country = 0; country < HOUSE.TotalCountries; country++) PARAMS.rhos[country][0] = (DType)1; //// i.e. only considering a single serotype so equal to 1. 
}
void PopulateParamRanges				(Params_Struct &PARAMS, Housekeeping_Struct &HOUSE)
{
#ifdef PRINT_PROGRAM_PROGRESS
	std::cerr << "PopulateParamRanges ";
#endif
	
	Allocate_2D_Array(PARAMS.ParamRanges, 2, HOUSE.No_Parameters); 

	string ParamRangeFileName_Full = "ParamRanges\\" + HOUSE.ParamRangeFileName + ".txt";
	ifstream ParamRangeInput; 
	string ParamName, LowerBound_String, UpperBound_String; 
	
	if (FileExists(ParamRangeFileName_Full))
	{
		ParamRangeInput.open(ParamRangeFileName_Full);
		while (!ParamRangeInput.eof())
		{
			std::getline(ParamRangeInput, ParamName			, '\t');
			std::getline(ParamRangeInput, LowerBound_String	, '\t');
			std::getline(ParamRangeInput, UpperBound_String	, '\n');

					if (ParamName == "KAM_0" 					)	{ HOUSE.IntialParRanges.KAM_0					[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.KAM_0					[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "KAM_1" 					)	{ HOUSE.IntialParRanges.KAM_1					[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.KAM_1					[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "KAM_2" 					)	{ HOUSE.IntialParRanges.KAM_2					[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.KAM_2					[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "KPS_0" 					)	{ HOUSE.IntialParRanges.KPS_0					[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.KPS_0					[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "KPS_1" 					)	{ HOUSE.IntialParRanges.KPS_1					[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.KPS_1					[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "KPS_2" 					)	{ HOUSE.IntialParRanges.KPS_2					[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.KPS_2					[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "Hosp_K_multiplier"		)	{ HOUSE.IntialParRanges.Hosp_K_multiplier		[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.Hosp_K_multiplier		[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "SNegEff_1"				)	{ HOUSE.IntialParRanges.SNegEff_1				[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.SNegEff_1				[UpperBound] = std::stod(UpperBound_String); }
 			else	if (ParamName == "SPosEff_1"				)	{ HOUSE.IntialParRanges.SPosEff_1				[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.SPosEff_1				[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "SeroNegWaning"			)	{ HOUSE.IntialParRanges.SeroNegWaning			[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.SeroNegWaning			[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "SeroPosWaning"			)	{ HOUSE.IntialParRanges.SeroPosWaning			[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.SeroPosWaning			[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "HillHalfLife_SNeg"		)	{ HOUSE.IntialParRanges.HillHalfLife_SNeg		[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.HillHalfLife_SNeg		[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "HillHalfLife_SPos"		)	{ HOUSE.IntialParRanges.HillHalfLife_SPos		[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.HillHalfLife_SPos		[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "knots"					)	{ HOUSE.IntialParRanges.knots					[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.knots					[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "knots_logged"				)	{ HOUSE.IntialParRanges.knots_logged			[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.knots_logged			[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "Late_knots"				)	{ HOUSE.IntialParRanges.Late_knots				[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.Late_knots				[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "Late_knots_logged"		)	{ HOUSE.IntialParRanges.Late_knots_logged		[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.Late_knots_logged		[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "SeroNegWaning_Flipped"	)	{ HOUSE.IntialParRanges.SeroNegWaning_Flipped	[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.SeroNegWaning_Flipped	[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "SeroPosWaning_Flipped"	)	{ HOUSE.IntialParRanges.SeroPosWaning_Flipped	[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.SeroPosWaning_Flipped	[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "HillPower_SNeg"			)	{ HOUSE.IntialParRanges.HillPower_SNeg			[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.HillPower_SNeg			[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "HillPower_SPos"			)	{ HOUSE.IntialParRanges.HillPower_SPos			[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.HillPower_SPos			[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "HistHaz"					)	{ HOUSE.IntialParRanges.HistHaz					[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.HistHaz				[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "ASVEHalflife"				)	{ HOUSE.IntialParRanges.ASVEHalflife			[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.ASVEHalflife			[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "ASVEPower"				)	{ HOUSE.IntialParRanges.ASVEPower				[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.ASVEPower				[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "ASVEProp"					)	{ HOUSE.IntialParRanges.ASVEProp				[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.ASVEProp				[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "AS_Haz_Halflife"			)	{ HOUSE.IntialParRanges.AS_Haz_Halflife			[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.AS_Haz_Halflife		[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "AS_Haz_Power"				)	{ HOUSE.IntialParRanges.AS_Haz_Power			[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.AS_Haz_Power			[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "AS_Haz_Prop"				)	{ HOUSE.IntialParRanges.AS_Haz_Prop				[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.AS_Haz_Prop			[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "ASWaning_Halflife"		)	{ HOUSE.IntialParRanges.ASWaning_Halflife		[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.ASWaning_Halflife		[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "ASWaning_Power"			)	{ HOUSE.IntialParRanges.ASWaning_Power			[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.ASWaning_Power			[UpperBound] = std::stod(UpperBound_String); }
			else	if (ParamName == "qval"						)	{ HOUSE.IntialParRanges.qval					[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.qval					[UpperBound] = std::stod(UpperBound_String); }
  			else	if (ParamName == "rho"						)	{ HOUSE.IntialParRanges.rho						[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.rho					[UpperBound] = std::stod(UpperBound_String); }
  			else	if (ParamName == "BS_BaseHazMult"			)	{ HOUSE.IntialParRanges.BS_BaseHazMult			[LowerBound] = std::stod(LowerBound_String)		;  HOUSE.IntialParRanges.BS_BaseHazMult			[UpperBound] = std::stod(UpperBound_String); }
			else	std::cerr << "PopulateParamRanges ERROR: ParamName " << ParamName << " not recognized from " + ParamRangeFileName_Full << endl;
		}
		ParamRangeInput.close();
	}
	else std::cerr << "PopulateParamRanges ERROR: ParamRangeFileName_Full " << ParamRangeFileName_Full << " not recognized" << endl;

	for (int param = 0; param < HOUSE.No_Parameters; param++)
	{
				if (IsParamAWaningParam(param, PARAMS.ParamNumbers))
		{

			int ParamType				= Find_ParamType_From_ASWaning_Param(param, HOUSE, PARAMS.ParamNumbers); 
			int BaselineSeroStatus		= Find_SeroStatus_From_Waning_Param	(param, HOUSE, PARAMS.ParamNumbers, ParamType); 
	
			bool ParamIsADurationOrRate =	ParamType == WaningRateDuration					|| 
											HOUSE.AS_Waning == Age_Option::CATEGORICAL		|| 
											HOUSE.AS_Waning == Age_Option::SPLINE			|| 
											HOUSE.AS_Waning == Age_Option::SPLINE_LINE		|| 
											HOUSE.AS_Waning == Age_Option::SPLINE_STEP		|| 
											HOUSE.AS_Waning == Age_Option::CUBIC			; 

				 if (IsParamAHillHalflife	(param, PARAMS.ParamNumbers)	&& BaselineSeroStatus == SeroNeg)	{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.HillHalfLife_SNeg	[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.HillHalfLife_SNeg	[UpperBound];	}
			else if (IsParamAHillHalflife	(param, PARAMS.ParamNumbers)	&& BaselineSeroStatus == SeroPos)	{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.HillHalfLife_SPos	[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.HillHalfLife_SPos	[UpperBound];	}
			else if (IsParamAHillPower		(param, PARAMS.ParamNumbers)	&& BaselineSeroStatus == SeroNeg)	{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.HillPower_SNeg	[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.HillPower_SNeg	[UpperBound];	}
			else if (IsParamAHillPower		(param, PARAMS.ParamNumbers)	&& BaselineSeroStatus == SeroPos)	{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.HillPower_SPos	[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.HillPower_SPos	[UpperBound];	}
			else if (BaselineSeroStatus == SeroNeg && ParamIsADurationOrRate)									{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.SeroNegWaning		[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.SeroNegWaning		[UpperBound];	}
			else if (BaselineSeroStatus == SeroPos && ParamIsADurationOrRate)									{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.SeroPosWaning		[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.SeroPosWaning		[UpperBound];	}
			else if (ParamType == Halflife_Index)																{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.ASWaning_Halflife	[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.ASWaning_Halflife	[UpperBound]; }
			else if (ParamType == Power_Index	)																{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.ASWaning_Power	[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.ASWaning_Power	[UpperBound];	}
		}
		else	if (IsParamAKnot		(param, PARAMS.ParamNumbers))
		{	
			if (Find_Knot_FromKnotParam(param, HOUSE, PARAMS.ParamNumbers) <= 6)
				if (HOUSE.LinKnts)
				{
					PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.knots[LowerBound];
					PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.knots[UpperBound];
				}
				else
				{
					PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.knots_logged[LowerBound];
					PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.knots_logged[UpperBound];
				}
			else						//// i.e. for passive phase knots)
				if (HOUSE.LinKnts)
				{
					PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.Late_knots[LowerBound];
					PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.Late_knots[UpperBound];
				}
				else
				{
					PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.Late_knots_logged[LowerBound];
					PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.Late_knots_logged[UpperBound];
				}
		}
		else	if (IsParamAHistHaz		(param, PARAMS.ParamNumbers))											{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.HistHaz				[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.HistHaz				[UpperBound];	}
		else	if (IsParamAnEfficacy	(param, PARAMS.ParamNumbers) || IsParamAn_atInf_Efficacy(param, PARAMS.ParamNumbers))	
		{
			int BaselineSeroStatus = Find_SeroStatus_FromEffParam(param, HOUSE, PARAMS.ParamNumbers);
			if (HOUSE.ASVE != Age_Option::INDEPENDENT & HOUSE.SeroSpecificEfficacies == 1) ///// in this scenario, the parameters you've labelled as vaccine efficacies are actually multipliers or intercepts (the AgeEffMults are the real efficacies). hence why they can go higher than 1. 
			{
				if (HOUSE.SSASVE_Additive)
				{
					/*Lower Bounds*/
					// true for both seropositive and seronegative. Must allow efficacy of serotypes 2,3,4 to be LESS than serotype 1. 

					if (BaselineSeroStatus = SeroNeg & HOUSE.IntialParRanges.SNegEff_1[LowerBound] < 0)  PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.SNegEff_1[LowerBound]; 
					else PARAMS.ParamRanges[LowerBound][param] = -1;

					/*Upper Bounds*/
					PARAMS.ParamRanges[UpperBound][param] = 1;
				}
				else  // i.e. multiplicative
				{  
					if (BaselineSeroStatus = SeroNeg & HOUSE.IntialParRanges.SNegEff_1[LowerBound] < 0)  PARAMS.ParamRanges[LowerBound][param] = -5; else PARAMS.ParamRanges[LowerBound][param] = 0;
					PARAMS.ParamRanges[UpperBound][param] = 6;
				}
			}
			else
			{
					 if (BaselineSeroStatus == SeroNeg)				{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.SNegEff_1			[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.SNegEff_1			[UpperBound];	}
				else if (BaselineSeroStatus == SeroPos)				{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.SPosEff_1			[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.SPosEff_1			[UpperBound];	}
			}
		}
		else	if (IsParamA_qval		(param, PARAMS.ParamNumbers))											{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.qval					[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.qval					[UpperBound];	}
		else	if (IsParamA_rho		(param, PARAMS.ParamNumbers))											{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.rho					[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.rho					[UpperBound];	}
		else	if (IsParamARelativeRisk(param, PARAMS.ParamNumbers) || IsParamARelativeRisk_Hosp(param, PARAMS.ParamNumbers))
		{
			int PrevInf			= Find_PrevInf_From_K_Param			(param, HOUSE, PARAMS.ParamNumbers); 
			int PhaseSeverity	= Find_PhaseSeverity_From_K_Param	(param, HOUSE, PARAMS.ParamNumbers, PrevInf); 

			if (PhaseSeverity == ActiveMild)
			{
				if (PrevInf == 0) 																		{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.KAM_0				[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.KAM_0				[UpperBound];	}
				if (PrevInf == 1) 																		{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.KAM_1				[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.KAM_1				[UpperBound];	}
				if (PrevInf == 2) 																		{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.KAM_2				[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.KAM_2				[UpperBound];	}
			}
			else if (PhaseSeverity == PassiveSevere)
			{
				if (PrevInf == 0) 																		{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.KPS_0				[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.KPS_0				[UpperBound];	}
				if (PrevInf == 1) 																		{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.KPS_1				[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.KPS_1				[UpperBound];	}
				if (PrevInf == 2) 																		{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.KPS_2				[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.KPS_2				[UpperBound];	}
			}
		}
		else	if (param == PARAMS.ParamNumbers.Hosp_K_multiplier	)											{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.Hosp_K_multiplier	[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.Hosp_K_multiplier	[UpperBound];	}
		else	if (IsParamAn_ASVE_Param(param, PARAMS.ParamNumbers)) 
		{
			int Type				= Find_Type_From_ASVE_Param			(param, HOUSE, PARAMS.ParamNumbers); 
			int BaselineSeroStatus	= Find_SeroStatus_FromASVE_Param	(param, HOUSE, PARAMS.ParamNumbers, Type);

			if (HOUSE.ASVE == Age_Option::HILL)
			{
					 if (Type == Halflife_Index	)		{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.ASVEHalflife	[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.ASVEHalflife	[UpperBound];	}
				else if (Type == Power_Index	)		{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.ASVEPower		[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.ASVEPower		[UpperBound];	}
				else if (Type == Prop_Index		)		{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.ASVEProp		[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.ASVEProp		[UpperBound];	}
			}
			else if (HOUSE.ASVE != Age_Option::INDEPENDENT)
			{
				  	 if (BaselineSeroStatus == SeroNeg)	{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.SNegEff_1	[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.SNegEff_1	[UpperBound];	}
				else if (BaselineSeroStatus == SeroPos)	{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.SPosEff_1	[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.SPosEff_1	[UpperBound];	}
			}
			else std::cerr << "PopulatePARAMS.ParamRanges ERROR: ASVE param if/else wrong" << endl; 
		}
		else	if (IsParamAn_ASHaz_Param(param, PARAMS.ParamNumbers))
		{
			int Type = Find_Type_From_ASHaz_Param(param, HOUSE, PARAMS.ParamNumbers); 

			if (HOUSE.AS_Haz == Age_Option::HILL)
			{
					 if (Type == Halflife_Index	)		{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.AS_Haz_Halflife			[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.AS_Haz_Halflife			[UpperBound];	}
				else if (Type == Power_Index	)		{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.AS_Haz_Power				[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.AS_Haz_Power				[UpperBound];	}
				else if (Type == Prop_Index		)		{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.AS_Haz_Prop				[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.AS_Haz_Prop				[UpperBound];	}
			}
			else 
			{
				PARAMS.ParamRanges[LowerBound][param] = 0;
				PARAMS.ParamRanges[UpperBound][param] = 5;
			}
		}
		else	if (IsParamAn_AS_Prime_Param(param, PARAMS.ParamNumbers))
		{
			int Type				= Find_ParamType_From_AS_Prime_Param	(param, HOUSE, PARAMS.ParamNumbers);  
			int BaselineSeroStatus	= Find_SeroStatus_From_AS_Prime_Param	(param, HOUSE, PARAMS.ParamNumbers, Type);

			std::cout << "PopulatePARAMS.ParamRanges_AS_Prime: BS " << BaselineSeroStatus << " Type " << Type << endl;

				 if (Type == Prop_Index				) { PARAMS.ParamRanges[LowerBound][param] = 0	;	PARAMS.ParamRanges[UpperBound][param] = 1; }
			else if (Type == AS_Priming_Rate_index	) { PARAMS.ParamRanges[LowerBound][param] = 0.1;	PARAMS.ParamRanges[UpperBound][param] = 30; }
			else std::cout << "PopulatePARAMS.ParamRanges error: AS_Prime_Param type not recognized" << endl;

		}
		else	if (param == PARAMS.ParamNumbers.BS_BaseHazMult) 
		{
			PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.BS_BaseHazMult[LowerBound];
			PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.BS_BaseHazMult[UpperBound];
		 }
		else	std::cerr << endl << endl << "you haven't specified upper limits for every parameter. " << endl << endl;
	}

	if (HOUSE.FitWaningRate)
		for (int param = 0; param < HOUSE.No_Parameters; param++)
			if (IsParamAWaningParam(param, PARAMS.ParamNumbers))
			{
				int ParamType			= Find_ParamType_From_ASWaning_Param(param, HOUSE, PARAMS.ParamNumbers); 
				int BaselineSeroStatus	= Find_SeroStatus_From_Waning_Param	(param, HOUSE, PARAMS.ParamNumbers, ParamType); 

					 if (BaselineSeroStatus == SeroNeg && ParamType == WaningRateDuration)  {	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.SeroNegWaning_Flipped[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.SeroNegWaning_Flipped[UpperBound];	}
				else if (BaselineSeroStatus == SeroPos && ParamType == WaningRateDuration)	{	PARAMS.ParamRanges[LowerBound][param] = HOUSE.IntialParRanges.SeroPosWaning_Flipped[LowerBound];	PARAMS.ParamRanges[UpperBound][param] = HOUSE.IntialParRanges.SeroPosWaning_Flipped[UpperBound];	}
			}

#ifdef PRINT_PROGRAM_PROGRESS
	std::cerr << "DONE ";
#endif
}
void PopulateStandardDevs				(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
#ifdef PRINT_PROGRAM_PROGRESS
	std::cerr << "PopulateStandardDevs ";
#endif
	PARAMS.ProposalStandDevs = new DType[HOUSE.No_Parameters]();

	for (int param = 0; param < HOUSE.No_Parameters; param++)
	{
				if (IsParamAWaningParam(param, PARAMS.ParamNumbers))
		{
			int ParamType				= Find_ParamType_From_ASWaning_Param(param, HOUSE, PARAMS.ParamNumbers); 
			int BaselineSeroStatus		= Find_SeroStatus_From_Waning_Param	(param, HOUSE, PARAMS.ParamNumbers, ParamType); 
			bool ParamIsADurationOrRate =	ParamType == WaningRateDuration				|| 
											HOUSE.AS_Waning == Age_Option::CATEGORICAL		|| 
											HOUSE.AS_Waning == Age_Option::SPLINE			|| 
											HOUSE.AS_Waning == Age_Option::SPLINE_LINE		|| 
											HOUSE.AS_Waning == Age_Option::SPLINE_STEP		|| 
											HOUSE.AS_Waning == Age_Option::CUBIC			; 

				 if	(ParamIsADurationOrRate		)	PARAMS.ProposalStandDevs[param] = 2;
			else if (ParamType == Power_Index	) 	PARAMS.ProposalStandDevs[param] = 1;			//// power for vaccine if using Hill function 
			else if (ParamType == Halflife_Index) 	PARAMS.ProposalStandDevs[param] = 1;			//// age at which vaccine attains highest efficacy. 
		}
		else	if (IsParamAKnot		(param, PARAMS.ParamNumbers)													)	PARAMS.ProposalStandDevs[param] = 0.005	;
		else	if (IsParamAHistHaz		(param, PARAMS.ParamNumbers)													)	PARAMS.ProposalStandDevs[param] = 0.01		;
		else	if (IsParamA_qval		(param, PARAMS.ParamNumbers)													)	PARAMS.ProposalStandDevs[param] = 0.2		;
		else	if (IsParamA_rho		(param, PARAMS.ParamNumbers)													)	PARAMS.ProposalStandDevs[param] = 0.1		;
		else	if (IsParamAnEfficacy	(param, PARAMS.ParamNumbers) || IsParamAn_atInf_Efficacy(param, PARAMS.ParamNumbers	))	
		{
			if (HOUSE.ASVE != Age_Option::INDEPENDENT & HOUSE.SeroSpecificEfficacies == 1 & !HOUSE.SSASVE_Additive)	PARAMS.ProposalStandDevs[param] = 0.25; ///// in this scenario, the parameters you've labelled as vaccine efficacies are actually multipliers (the AgeEffMults are the real efficacies). 
			else																									PARAMS.ProposalStandDevs[param] = 0.1;
		}
		else	if (IsParamAHillHalflife(param, PARAMS.ParamNumbers)													)	PARAMS.ProposalStandDevs[param] = 1		;
		else	if (IsParamAHillPower	(param, PARAMS.ParamNumbers)													)	PARAMS.ProposalStandDevs[param] = 1		;
		else	if (IsParamARelativeRisk(param, PARAMS.ParamNumbers) || IsParamARelativeRisk_Hosp(param, PARAMS.ParamNumbers	))
		{
			int PrevInf		= Find_PrevInf_From_K_Param(param, HOUSE, PARAMS.ParamNumbers);;
			DType K_Step	= 0.05; 

			PARAMS.ProposalStandDevs[param] = K_Step;
		}		
		else	if (IsParamAn_ASVE_Param(param, PARAMS.ParamNumbers))
		{
			int Type = Find_Type_From_ASVE_Param(param, HOUSE, PARAMS.ParamNumbers); //// either half life or power or prop

			if (HOUSE.ASVE == Age_Option::HILL)
			{
					 if (Type == Power_Index	) 	PARAMS.ProposalStandDevs[param] = 1;			//// Power for vaccine if using Hill function 
				else if (Type == Halflife_Index	) 	PARAMS.ProposalStandDevs[param] = 1;			//// age at which vaccine attains highest efficacy. 
				else if (Type == Prop_Index		) 	PARAMS.ProposalStandDevs[param] = 0.1;			//// age at which vaccine attains highest efficacy. 
			}
			else PARAMS.ProposalStandDevs[param] = 0.1;
		}
		else	if (IsParamAn_ASHaz_Param(param, PARAMS.ParamNumbers))
		{
			int Type = Find_Type_From_ASHaz_Param(param, HOUSE, PARAMS.ParamNumbers);	//// either half life or power or prop

			if (HOUSE.AS_Haz == Age_Option::HILL)
			{
					 if (Type == Power_Index	) 	PARAMS.ProposalStandDevs[param] = 1;		//// Power for vaccine if using Hill function 
				else if (Type == Halflife_Index	) 	PARAMS.ProposalStandDevs[param] = 1;		//// age at which vaccine attains highest efficacy. 
				else if (Type == Prop_Index		) 	PARAMS.ProposalStandDevs[param] = 0.1;		//// age at which vaccine attains highest efficacy. 
			}
			else PARAMS.ProposalStandDevs[param] = 1;
		}
		else	if (param == PARAMS.ParamNumbers.Hosp_K_multiplier)												PARAMS.ProposalStandDevs[param] = 0.5;
		else	if (IsParamAn_AS_Prime_Param(param, PARAMS.ParamNumbers))
		{
			int Type				= Find_ParamType_From_AS_Prime_Param	(param, HOUSE, PARAMS.ParamNumbers); //// either half life or coeff. 
			int BaselineSeroStatus	= Find_SeroStatus_From_AS_Prime_Param	(param, HOUSE, PARAMS.ParamNumbers, Type);

				 if (Type == Prop_Index				) PARAMS.ProposalStandDevs[param] = 0.5;
			else if (Type == AS_Priming_Rate_index	) PARAMS.ProposalStandDevs[param] = 2;
			else std::cout << "PopulateStandardDevs error: AS_Prime_Param type not recognized" << endl;
		}
		else	if (param == PARAMS.ParamNumbers.BS_BaseHazMult) PARAMS.ProposalStandDevs[param] = 0.1;
		else	std::cerr << endl << endl << "PopulateStandardDevs error: Param number not recognized." << endl;
	}

	if (!HOUSE.LinKnts)
		for (int param = 0; param < HOUSE.No_Parameters; param++)
			if (IsParamAKnot(param, PARAMS.ParamNumbers)) PARAMS.ProposalStandDevs[param] = 0.5;
				
	if (HOUSE.FitWaningRate)
		for (int param = 0; param < HOUSE.No_Parameters; param++)
			if (IsParamAWaningParam(param, PARAMS.ParamNumbers))
			{
				int ParamType			= Find_ParamType_From_ASWaning_Param(param, HOUSE, PARAMS.ParamNumbers); 
				int BaselineSeroStatus	= Find_SeroStatus_From_Waning_Param	(param, HOUSE, PARAMS.ParamNumbers, ParamType); 

					 if (BaselineSeroStatus == SeroNeg && ParamType == WaningRateDuration) PARAMS.ProposalStandDevs[param] = 1;
				else if (BaselineSeroStatus == SeroPos && ParamType == WaningRateDuration) PARAMS.ProposalStandDevs[param] = 1;

			}
#ifdef PRINT_PROGRAM_PROGRESS
	std::cerr << "DONE ";
#endif
}
void Initialize_SumRhoEffNegs			(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	///// this requires rhos, Efficacies, AgeEffMults AND parameter ranges to be initialized, hence why it gets its own separate function to be called in Initialize_Params after above quantities allocated and calculated. 
	if (HOUSE.EffNegWane == EffNegWane_Option::FROM_ZERO || HOUSE.EffNegWane == EffNegWane_Option::NO_WANE)
	{
		Allocate_4D_Array(PARAMS.SumRhoEffNegs, HOUSE.HowManyCaseCategories, HOUSE.TotalCountries, HOUSE.HowManyAges, HOUSE.HowManySeroStatuses);  
		if (!HOUSE.PSVEs && HOUSE.HowManyCaseCategories == 2)  PARAMS.SumRhoEffNegs[PassiveSevere] = PARAMS.SumRhoEffNegs[ActiveMild]; //// i.e. make both PhaseSeverities point to same value.
		Calc_SumRhoEffNegs_All(PARAMS, HOUSE);
	}
}
void Initialize_Params					(Params_Struct &PARAMS, Params_Struct &CurrentPARAMS, Housekeeping_Struct &HOUSE, const DATA_struct &DATA)
{
#ifdef PRINT_PROGRAM_PROGRESS
	std::cerr << "Initialize_Params " << endl;
#endif
	InitializeLikeArray(PARAMS, HOUSE);

	////// ////// ////// ////// ////// ////// ////// ////// ////// ////// 
	////// Build ParamVecs And Names. Must add to ParamVec (the actual value), NamesOfParameters (it's name), and ParamNumbers. The latter is to try to get around hash defines and strings. 
	////// Don't want loads of string comparisons for say AmendProposedParams, but also can't use #defines due to parameter numbers varying by ModelVariant, trial phase and whether or not you're doing sero-specific efficacies. 
	////// ////// ////// ////// ////// ////// ////// ////// ////// ////// 
	int ParamCounter = 0;
	InitializeCommandLineParams(PARAMS, CurrentPARAMS); 

	////// Relative Risks
	InitializeRelativeRisks(ParamCounter, PARAMS, HOUSE); 

	////// Efficacies
	InitializeEfficacies(ParamCounter, PARAMS, HOUSE); 

	////// Waning Values
	InitializeWaningParamsAndValues(ParamCounter, PARAMS, HOUSE, DATA); 

	////// Age Efficacy multipliers
	InitializeASVEParams(ParamCounter, PARAMS, HOUSE); 

	//////// Baseline hazard age multipliers
	Initialize_AS_Haz_Params(ParamCounter, PARAMS, HOUSE); 

	//////// AS_PrimingParams
	Initialize_AS_PrimingParams(ParamCounter, PARAMS, HOUSE);

	////// Baseline hazard baseline serostatus multipliers
	Initialize_BS_BaseHazMults(ParamCounter, PARAMS, HOUSE); 

	////// rhos	
	Initialize_Rhos(ParamCounter, PARAMS, HOUSE);

	////// Knots, spline coefficients, then calculate: baseline hazard values, integrated baseline hazard values, integrated vaccine hazard values. 
	//// Allocate memory
	Allocate_2D_Array(PARAMS.IntBaseHazLookUp		, HOUSE.TotalCountries, DATA.MaxStartFollowUp_CalendarTime + DATA.NoDaysOfFollowUp + 1);
	Allocate_2D_Array(PARAMS.BaselineHazardValues	, HOUSE.TotalCountries, DATA.MaxStartFollowUp_CalendarTime + DATA.NoDaysOfFollowUp + 1);

	CreateAndPopulateKnots			(PARAMS, HOUSE);
	CreateAndPopulateSplineCoeffs	(PARAMS, HOUSE);
	for (int country = 0; country < HOUSE.TotalCountries; country++)
	{
		Calc_BaseHazValues		(country, PARAMS, DATA, HOUSE); 
		Calc_IntBaseHazValues	(country, PARAMS, DATA, HOUSE);
	}
	PARAMS. ParamNumbers.Min_Knot = ParamCounter;
	for (int country = 0; country < HOUSE.TotalCountries; country++)
		for (int knot = 0; knot < HOUSE.KnotsPerCountry; knot++) 
		{
			PARAMS.ParamVec.			push_back(PARAMS. yKnots[country][knot]); //// knots are calculated previously. 
			PARAMS.NamesOfParameters.	push_back("Knot_" + std::to_string(country) + "_" + std::to_string(knot));

			ParamCounter++;
		}
	PARAMS. ParamNumbers.Max_Knot = ParamCounter - 1;

	////// Historical Hazards	
	PARAMS. ParamNumbers.Min_HHaz = ParamCounter;
	DType *FixedHistHazrds = new DType[HOUSE.TotalCountries]();
	ifstream HistHazards; HistHazards.open(HOUSE.HH_InputFilename);
	for (int country = 0; country < HOUSE.TotalCountries; country++)
	{
		HistHazards >> FixedHistHazrds[country];	
		PARAMS.ParamVec.			push_back(FixedHistHazrds[country]); // add the historical hazards to ParamVec
		PARAMS.NamesOfParameters.	push_back("h_" + std::to_string(country));
		ParamCounter++;
	}
	delete[] FixedHistHazrds;
	HistHazards.close(); 
	PARAMS. ParamNumbers.Max_HHaz = ParamCounter - 1;


	///// Initialize ParamSeroPrevalence (Need the PARAMS.SeroPrevs array in any case, but so long as you're not using Empirical Seroprevalence in likelihood and augmentation, then PARAMS.SeroPrevs is based solely on Paramteres
	//// Parameter Seroprevalence arrays
	Allocate_4D_Array(PARAMS.SeroPrevs, 2, HOUSE.TotalCountries, HOUSE.HowManySeroStatuses, HOUSE.HowManyAges); //// 2 for log-or-not-log;

	if (!HOUSE.Empirical_SeroPrevs) //// i.e. default using historical hazards to assess seroprevalence in likelihood and augmentation. 
		for (int country = 0; country < HOUSE.TotalCountries; country++)
			Calc_ParamSeroPrevs_country(country, PARAMS, HOUSE);

	// knots
	if (!HOUSE.SkipKnots)
		for (int country_index = 0; country_index < HOUSE.WhichCountries.size(); country_index++)
			for (int knotparam = (HOUSE.WhichCountries[country_index] * HOUSE.KnotsPerCountry); knotparam < ((HOUSE.WhichCountries[country_index] + 1) * HOUSE.KnotsPerCountry); knotparam++)
				PARAMS.ParamNosYouWillFit.push_back(PARAMS.ParamNumbers.Min_Knot + knotparam);
	
	// historical hazards
	if (!HOUSE.SkipHistHazards)
		for (int country_index = 0; country_index < HOUSE.WhichCountries.size(); country_index++)
			PARAMS.ParamNosYouWillFit.push_back(PARAMS.ParamNumbers.Min_HHaz + HOUSE.WhichCountries[country_index]);

	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
	///// ///// /////  Various derivatives / look up tables associated with parameter vectors
	
	///// allocate memory for "rho efficacy products" - do this even if not considering serospecific efficacies as means fewer if statements in costly augmentation. Also useful in survival curves. 
	Allocate_3D_Array(PARAMS.SumRhoEffs, HOUSE.HowManyCaseCategories, HOUSE.TotalCountries, HOUSE.HowManySeroStatuses);

	if (!HOUSE.PSVEs && HOUSE.HowManyCaseCategories == 2) 
		PARAMS.SumRhoEffs[PassiveSevere] = PARAMS.SumRhoEffs[ActiveMild]; 
	Calc_SumRhoEffs_ALL(PARAMS, HOUSE);

	///// Kplus Values
	Initialize_KplusValues(PARAMS, HOUSE); //// Must be done after Historical hazards have been added to ParamVec. 

	//// SumRhoKs
	Allocate_3D_Array(PARAMS.SumRhoKs, HOUSE.TotalCountries, HOUSE.HowManyCaseCategories, HOUSE.Num_K_Params); //// must be after rhos, and Ks defined
	Calc_SumRhoKs_ALL(PARAMS, HOUSE);

	//// SumRhoEffKs
	if (HOUSE.SeroSpecificEfficacies && HOUSE.SeroSpecific_K_values)
		Allocate_3D_Array(PARAMS.SumRhoEffKs, HOUSE.TotalCountries, HOUSE.HowManyCaseCategories, HOUSE.HowManySeroStatuses);
	else if (HOUSE.SeroSpecific_K_values)
		PARAMS.SumRhoEffKs = PARAMS.SumRhoKs; 
	Calc_SumRhoEffKs_ALL(PARAMS, HOUSE);

	//// AS_PrimingParams
	if (HOUSE.ModelVariant == AS_PRIME)
	{
		Calc_K_Primes_ALL					(PARAMS, HOUSE);
		Calc_KplusPrimes_All			(PARAMS, HOUSE); 
		Calc_SumRhoK0_SNeg_Primes_ALL	(PARAMS, HOUSE);
	}

	HOUSE.No_Parameters = ParamCounter;

	PopulateParamRanges	(PARAMS, HOUSE); 
	PopulateStandardDevs(PARAMS, HOUSE); 
	Initialize_SumRhoEffNegs(PARAMS, HOUSE); ///// used if HOUSE.AltWaneEffNeg

	//// CalcPriors. 
	PARAMS.ParamArray_logPriors = new DType[HOUSE.No_Parameters](); //// allocate array. 
	CalculateLogPriors_all(PARAMS, HOUSE); //// calculates all log prior components and their sum. Non-knot components won't change (uniform) but values may still be important with comparing DICs between model runs (e.g. when using different param ranges / prior files). 

#ifdef PRINT_PROGRAM_PROGRESS
	std::cerr << "DONE ";
#endif
}
void Initialize_FollowUp				(DATA_struct &DATA, const Housekeeping_Struct &HOUSE)
{
	////*****////				////*****////			////*****////
	////*****////		Explanation of defnitions		////*****////
	////*****////				////*****////			////*****////

	////*****///// DATA.NoDaysOfFollowUp refers to maximum follow up duration over all patients. Will be different depending on whether modelling Active only or Active+Passive, whether using Sanofi's follow up definitions, and whether including late cases. 
	//// if ACTIVE_PHASE_ONLY, need max(Data$EndActivePhase_InDays - Data$StartFollowUpDate_InDays). If DO_ACTIVE_AND_PASSIVE, need max(Data$EndPassivePhase_InDays - Data$StartFollowUpDate_InDays) in R Script.
	//// You need waning values for days 0 to DATA.NoDaysOfFollowUp. 

	////****///// SURVIVE.NoDaysOfFollowUp need not be the same as DATA.NoDaysOfFollowUp, although it can be. This is simply the number of days post first dose for which you want to calculate survival probabilities.
	//// If calculating attack rates based on the patient specific duration of passive surveillance, can't just have 3*365 = 1095 days - need max(Data$EndPassivePhase_InDays_WithoutCases - Data$StartFollowUpDate_InDays). 
	//// Note that if including late cases then  max(Data$EndPassivePhase_InDays_WithoutCases) < max(Data$EndPassivePhase_InDays), in which case SURVIVE.NoDaysOfFollowUp != DATA.NoDaysOfFollowUp. 

	if (HOUSE.SFU)
	{
		if (HOUSE.PASSIVE_PHASE_ONLY)
		{
			DATA.	NoDaysOfFollowUp		= 337; //// Should be 335, but add in a couple more to be safe in case of memory issues. 
			DATA.	NoActiveDaysFollowUp	= DATA.		NoDaysOfFollowUp; //// under Sanofi definitions, there is no single number of active surveillance days (it's patient specific)_
		}
		else 	if (HOUSE.ActiveOrPassivePhase == ACTIVE_PHASE_ONLY)
		{
			DATA.	NoDaysOfFollowUp		= 891; //// longest active phase duration, so need waning values for all of them. 
			DATA.	NoActiveDaysFollowUp	= DATA.		NoDaysOfFollowUp; //// under Sanofi definitions, there is no single number of active surveillance days (it's patient specific)_
		}
		else	if (HOUSE.ActiveOrPassivePhase == DO_ACTIVE_AND_PASSIVE)
		{
			DATA.NoDaysOfFollowUp		= 1227; //// longest active+passive phase duration, so you need waning values for all of them. 
			DATA.NoActiveDaysFollowUp	= 891;		// using Sanofi definitions then there is no single NoActiveDaysFollowUp value - different for all patients. 
		}
		if (!HOUSE.PASSIVE_PHASE_ONLY)	DATA.N_WaningDays = DATA.NoDaysOfFollowUp;	//// In ordinary circumstances these are the same. If PASSIVE_PHASE_ONLY, you change 
		else 							DATA.N_WaningDays = 1227;					//// if (HOUSE.PASSIVE_PHASE_ONLY)	set DATA.N_WaningDays to what DATA.NoDaysOfFollowUp would have been in passive phase. 
	}
	else
	{
				if (HOUSE.PASSIVE_PHASE_ONLY)
		{
			DATA.	NoDaysOfFollowUp		= 337; //// Should be 335, but add in a couple more to be safe in case of memory issues. 
			DATA.	NoActiveDaysFollowUp	= DATA.		NoDaysOfFollowUp; //// under Sanofi definitions, there is no single number of active surveillance days (it's patient specific)_
		}
		else 	if (HOUSE.ActiveOrPassivePhase == ACTIVE_PHASE_ONLY)
		{
			DATA.	NoDaysOfFollowUp		= 730;
			DATA.	NoActiveDaysFollowUp	= DATA.NoDaysOfFollowUp;
		}
		else	if (HOUSE.ActiveOrPassivePhase == DO_ACTIVE_AND_PASSIVE)
		{
			DATA.NoDaysOfFollowUp = 1095;
			DATA.NoActiveDaysFollowUp = 730;
		}
			 if (!HOUSE.PASSIVE_PHASE_ONLY)										DATA.N_WaningDays = DATA.NoDaysOfFollowUp;	//// In ordinary circumstances these are the same. If PASSIVE_PHASE_ONLY, you change 
		else DATA.N_WaningDays = 1095;					//// if (HOUSE.PASSIVE_PHASE_ONLY)	set DATA.N_WaningDays to what DATA.NoDaysOfFollowUp would have been in passive phase. 
	}	
}
void AllocateMemoryAndPopulateSets		(DATA_struct &DATA, const Housekeeping_Struct &HOUSE)
{
#ifdef PRINT_PROGRAM_PROGRESS
	std::cerr << "AllocateMemoryAndPopulateSets ";
#endif
	int SetsPerCountry = 9; // 7 sets per country: ALL (i.e. all individuals in that country)(index 0), IVc (index 1), IcVc (index 2), IV (index 3), and IcV (index 4), I (index 5), Ic (index 6).

	for (int Country = 0; Country < HOUSE.TotalCountries; Country++)
	{
		DATA.Set_Array	[Country] = new std::vector<int>[SetsPerCountry]();
		DATA.Cases_Array[Country] = new std::vector<int>*[SetsPerCountry](); //// extra dimension in Cases_Array

		for (int set = 0; set < SetsPerCountry; set++)	DATA.Cases_Array[Country][set] = new std::vector<int>[HOUSE.HowManyCaseCategories]();
	}

	PopulateSets(DATA, HOUSE);

	//// PopulateSets function only populates sets for SeroNeg and SeroPos (not "EitherSeroStatus", as this is unaffected by augmentation). Populate these sets here. 
	int country;
	bool ConditionToSkip; 
	for (int i = 0; i < NPat; i++)
	{
		country = DATA.ci_s[i];
		ConditionToSkip = !(std::any_of(HOUSE.WhichCountries.begin(), HOUSE.WhichCountries.end(), [&](int j) {return j == country; }));	
		if (ConditionToSkip) continue;

		(DATA.Set_Array[DATA.ci_s[i]][	HOUSE.Strata[DATA.Vi_s[i]		][EitherSeroStatus]	]).push_back(i);	// Any Baseline Immunity, patient's trial arm, country c_i 
		(DATA.Set_Array[DATA.ci_s[i]][	HOUSE.Strata[EitherTrialArm		][EitherSeroStatus]	]).push_back(i);	// Any Baseline Immunity, Either Arm, country c_i 

		if (DATA.IsCase[i])
		{
			(DATA.Cases_Array[DATA.ci_s[i]][HOUSE.Strata[DATA.Vi_s[i]		][EitherSeroStatus]][DATA.Case_PhaseSeverity[i]]).push_back(i);				// Any Baseline Immunity, patient's trial arm, country c_i (cases)
			(DATA.Cases_Array[DATA.ci_s[i]][HOUSE.Strata[EitherTrialArm		][EitherSeroStatus]][DATA.Case_PhaseSeverity[i]]).push_back(i);				// Any Baseline Immunity, Either Arm, country c_i (cases)
		}
	}
#ifdef PRINT_PROGRAM_PROGRESS
	std::cerr << "DONE " << endl;
#endif
}
void DeleteData							(DATA_struct &DATA, const Housekeeping_Struct &HOUSE)
{
	delete[] DATA.ai_s;
	delete[] DATA.Ii_s;
	delete[] DATA.Vi_s;
	delete[] DATA.ci_s;
	delete[] DATA.IsCase;
	delete[] DATA.IsCase_ActiveMild;
	delete[] DATA.IsCase_PassiveSevere;
	delete[] DATA.AgeGroup1;
	delete[] DATA.AgeGroup2;
	delete[] DATA.NumCalendarDaysFollowUp;

	delete[] DATA.CaseSerotype;

	DeAllocate_2D_Array(DATA.CountryMinMaxCalendarTime	, HOUSE.TotalCountries);

	delete[] DATA.IsCase_AMandPS;
}
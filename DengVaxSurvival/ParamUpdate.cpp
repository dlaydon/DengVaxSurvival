
#include "HeaderAndSwitches.h"
#include "ParamNumber.h"
#include "Likelihood.h"
#include "Splines.h"
#include "Probability.h"
#include "CalcParamsEtc.h"
#include "ConsolePrinting.h"

/*
This script has two kinds of functions:
	
	1)	Functions that amend PARAMS structures with various precalculated quantities, by calling the relevant funtions in CalcParamsEtc.cpp (see. note at top of CalcParamsEtc.cpp). e.g. changing rho will change SumRhoEffs, and could affect stored likelihood values of vaccine hazard etc. 
	2)	Functions that amend the likelihood components AFTER parameter sturctures have been set appropriately. These are the Change_ functions below. 
*/


//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** 
//// ****			1) AmendParams and various overloads essentially call functions in CalcParamsEtc when required, e.g. if changing historical hazrd, need to change array of KPlusValues BEFORE you amend the likelihood. 
//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** 


void ChooseMaxMinPhaseSeverities(const Housekeeping_Struct &HOUSE, int PhaseSeverity)
{
	if (HOUSE.PSVEs)
	{
		HOUSE.MinMaxPhaseSeveritiesToLoopOver[0] = PhaseSeverity;
		HOUSE.MinMaxPhaseSeveritiesToLoopOver[1] = PhaseSeverity + 1;
	}
	else 
	{
		HOUSE.MinMaxPhaseSeveritiesToLoopOver[0] = 0;
		HOUSE.MinMaxPhaseSeveritiesToLoopOver[1] = HOUSE.HowManyCaseCategories;
	}
}
void PopulateAndOrSwap_Ks_Array	(std::vector<int> WhichCountries, int PhaseSeverity, int PrevInf, int serotype, DType &New_Param, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	//// Populate K_s array - do all countries first (i.e. have two "country loops") as this makes it easier to update hospital K's if not doing ModelHosp_Indie_Ks. 
	/* 
	If HOUSE.PS_Ks_Multiply_AM_Ks, then 
		either	PhaseSeverity	== PassiveSevere, in which case the new parameter KPS_i (new param) multiplies the active KAM_i. 
		OR		PhaseSeverity	== ActiveMild	, in which case order is important. Need to divide KPS_i by old KAM_i and multiply by new KAM_i (New_Param) before resetting (otherwise KPS_i will recursively include every previous KAM_i value). Must then update Ks array for ActiveMild value of KAM_i
		
	*/
	/*NEW - PASS_Ks_MULTIPLY_ACT_Ks*/
	for (int countryindex = 0; countryindex < WhichCountries.size(); countryindex++)
		if (HOUSE.PS_Ks_Multiply_AM_Ks && (PhaseSeverity == PassiveSevere))
			PARAMS.K_s[WhichCountries[countryindex]][PassiveSevere][PrevInf][serotype]		=  New_Param * PARAMS.K_s[WhichCountries[countryindex]][ActiveMild][PrevInf][serotype]; //// the new parameter KPS_i (new param) multiplies the active KAM_i. 
		else //// so either HOUSE.PS_Ks_Multiply_AM_Ks == false or PhaseSeverity == ActiveMild
		{
			if (HOUSE.PS_Ks_Multiply_AM_Ks && (PhaseSeverity == ActiveMild))
				PARAMS.K_s[WhichCountries[countryindex]][PassiveSevere][PrevInf][serotype]	*= New_Param / PARAMS.K_s[WhichCountries[countryindex]][ActiveMild][PrevInf][serotype]; //// divide KPS_i by old KAM_i and multiply by new KAM_i (New_Param)
			
			PARAMS.K_s[WhichCountries[countryindex]][PhaseSeverity][PrevInf][serotype] = New_Param; //// don't be tempted to change PhaseSeverity here to ActiveMild or PassiveSevere. If not doing HOUSE.PS_Ks_Multiply_AM_Ks needs to apply to either PhaseSeverity
		} 

}

void AmendParams				(int &param_no, DType &New_Param, Params_Struct &PARAMS, const DATA_struct &DATA, const Housekeeping_Struct &HOUSE)
{
#ifdef PRINT_WITHIN_CHANGE_FUNCTIONS
#ifdef PRINT_WITHIN_CHANGE_KNOT_FUNCTIONS
	std::cout << "AmendParams p" << param_no << " " << PARAMS.NamesOfParameters[param_no] << " value " << New_Param << endl;
#else 
	if (!IsParamAKnot(param_no, PARAMS.ParamNumbers)) std::cout << "AmendParams p" << param_no << " " << PARAMS.NamesOfParameters[param_no] << " value " << New_Param;
#endif 
#endif
	PARAMS.ParamVec[param_no] = New_Param; //// All parameters require this change.

	//// **** PRIORS
	AmendLogPrior(param_no, PARAMS, HOUSE);  //// Changes LogPriorFull and PARAMS.ParamArray_logPriors[param_no].

	//// **** Parameter values and associated quantities (e.g. relative risks / precalculated integrated hazards etc.)
	if (IsParamARelativeRisk(param_no, PARAMS.ParamNumbers) || IsParamARelativeRisk_Hosp(param_no, PARAMS.ParamNumbers))
	{
		int PrevInf			= Find_PrevInf_From_K_Param			(param_no, HOUSE, PARAMS.ParamNumbers); 
		int PhaseSeverity	= Find_PhaseSeverity_From_K_Param	(param_no, HOUSE, PARAMS.ParamNumbers, PrevInf); 
		int serotype		= Find_Serotype_From_K_Param		(param_no, HOUSE, PARAMS.ParamNumbers, PrevInf, PhaseSeverity);;

		std::vector<int> WhichCountries;
		if (IsParamARelativeRisk(param_no, PARAMS.ParamNumbers))
			if (HOUSE.ModellingHospitalized && PhaseSeverity == PassiveSevere && HOUSE.ModelHosp_Indie_Ks)	WhichCountries = HOUSE.CYD_14_countries;	//// only do CYD-14 countries. 
			else																							WhichCountries = HOUSE.WhichCountries;		//// otherwise do all countries
		else if (IsParamARelativeRisk_Hosp(param_no, PARAMS.ParamNumbers))									WhichCountries = HOUSE.CYD_15_countries;

		if (HOUSE.SSKs_FitMatrix[PhaseSeverity][PrevInf]) //// i.e. if fitting all serotypes for this K_PhaseSeverity_PrevInf, then do the below one serotype (this one) at a time. Automatically true for !SSKs too as default values of matrix all true (and changed only if SSKs). 
			PopulateAndOrSwap_Ks_Array(WhichCountries, PhaseSeverity, PrevInf, serotype, New_Param, PARAMS, HOUSE);
		else //// so SSKs == true & HOUSE.SSKs_FitMatrix[PhaseSeverity][PrevInf] == false. i.e. this K_PhaseSeverity_PrevInf is not serotype specific. Hence you will amend all serotype Ks for this PhaseSeverity and PrevInf at the same time. 
			for (int serotypeDummy = 0; serotypeDummy < HOUSE.N_STypes_Ks; serotypeDummy++)
				PopulateAndOrSwap_Ks_Array(WhichCountries, PhaseSeverity, PrevInf, serotypeDummy, New_Param, PARAMS, HOUSE);

		int Min_PhaseSeverity = PhaseSeverity, Max_PhaseSeverity = PhaseSeverity + 1; //// i.e. default is that you do one phase severity at a time in loop below. 
		if (HOUSE.PS_Ks_Multiply_AM_Ks && PhaseSeverity == ActiveMild) { Min_PhaseSeverity = ActiveMild; Max_PhaseSeverity = HOUSE.HowManyCaseCategories; } //// i.e. if changing Passive only this doesn't affect active, so default definition okay. 

		for (int PhaseSeverity = Min_PhaseSeverity; PhaseSeverity < Max_PhaseSeverity; PhaseSeverity++)
		{
			//// then update the hospital K's (as Hosp_K_CYD15 = regular hosp K times Hosp K multiplier. 
			if (HOUSE.ModellingHospitalized && PhaseSeverity == PassiveSevere && !HOUSE.ModelHosp_Indie_Ks) Calc_Hospital_Ks(PARAMS, HOUSE);

			if (HOUSE.ModelVariant == AS_PRIME && PrevInf > 0)
			{
				//// Ks_Prime
				for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++) //// as PrevInf > 0, looks like you'll only need SeroPos, but K1 and K2 affect K0prime,- and K1prime,+
					for (int countryindex = 0; countryindex < WhichCountries.size(); countryindex++)
						Calc_K_Primes_BS_PhaseSeverity_country_serotype(PARAMS, HOUSE, BaselineSeroStatus, PhaseSeverity, WhichCountries[countryindex], serotype);

				//// KPlusPrimeValues and SumRhoK0_SNeg_Primes: Both consider SS_Ks, therefore has rhos. Therefore need to do all countries for this K, even if Ks same between countries. Note  don't store arrays of KplusValues or KplusPrimes for single/case serotypes, only "survivors". 
				for (int countryindex = 0; countryindex < WhichCountries.size(); countryindex++) //// have to do all countries here as
				{
					Calc_KplusPrimes_country_PhaseSeverity			(PARAMS, HOUSE, WhichCountries[countryindex], PhaseSeverity);
					Calc_SumRhoK0_SNeg_Primes_country_PhaseSeverity(PARAMS, HOUSE, WhichCountries[countryindex], PhaseSeverity); 
				}
			}

			for (int countryindex = 0; countryindex < WhichCountries.size(); countryindex++)
			{
				//// recalculate SumRhoKProds for that PrevInf and PhaseSeverity for all countries in WhichCountries. 
				Calc_SumRhoK_country_PhaseSeverity_PrevInf(PARAMS, HOUSE, WhichCountries[countryindex], PhaseSeverity, PrevInf);
			
				if (HOUSE.AdditiveSSASVEs && HOUSE.SeroSpecific_K_values &&												//// In this scenario, sero Ks combine with sero VEs and age VEs, and are therefore included in the #pragma loop for l_Int_Vac_Haz. Hence change in K1 or K2 requires recalculation. 
					(PrevInf > 0 /*|| HOUSE.ModelVariant == SIMPLE_NUMERICAL || HOUSE.ModelVariant == K_SEROPOS*/))		//// If SIMPLE_NUMERICAL or K_SEROPOS, K0 affects vaccine group. For VAC_SILENT it doesn't. AS_PRIME dealt with below. 
				{
					//// recalculate Kplus value
					Calc_KplusValues_country_PhaseSeverity	(PARAMS, HOUSE, WhichCountries[countryindex], PhaseSeverity); //// this also covers K_SEROPOS SS_Ks and SS_VEs together
					
					CalcAndStoreSumIntVacHazOverPatients	(PrevInf - 1, WhichCountries[countryindex], DATA, PARAMS, HOUSE, PhaseSeverity); //// Note PrevInf - 1. Because K2 => SeroPositives (1); K1 => Seronegatives (0). 
				}
				else if (!(HOUSE.ModelVariant == SIMPLE_NUMERICAL))
				{
					//// recalculate Kplus value
					Calc_KplusValues_country_PhaseSeverity(PARAMS, HOUSE, WhichCountries[countryindex], PhaseSeverity); //// this also covers K_SEROPOS SS_Ks and SS_VEs together

					//// and therefore VacHazLikes if necessary
					if (HOUSE.ModelVariant == K_SEROPOS || HOUSE.ModelVariant == AS_PRIME)
						CalcAndStoreSumIntVacHazOverPatients(SeroPos, WhichCountries[countryindex], DATA, PARAMS, HOUSE, PhaseSeverity);
					if (HOUSE.ModelVariant == AS_PRIME)
						CalcAndStoreSumIntVacHazOverPatients(SeroNeg, WhichCountries[countryindex], DATA, PARAMS, HOUSE, PhaseSeverity);
				}
			}
		}
	}
	else	if (param_no == PARAMS.ParamNumbers.Hosp_K_multiplier) //// will only be true if (HOUSE.ModellingHospitalized !HOUSE.ModelHosp_Indie_Ks) 
	{
		if (HOUSE.ModelVariant == AS_PRIME) std::cerr << "AS_PRIME not coded for non-indie hosp Ks" << endl << "AS_PRIME not coded for non-indie hosp Ks" << endl;

		PARAMS.Hosp_K_mult = New_Param; 
		Calc_Hospital_Ks(PARAMS, HOUSE);

		std::vector<int> WhichCountries = HOUSE.CYD_15_countries; 
		int PhaseSeverity = PassiveSevere; 
		for (int countryindex = 0; countryindex < WhichCountries.size(); countryindex++)
			for (int PrevInf = 0; PrevInf < HOUSE.Num_K_Params; PrevInf++)
			{
				//// recalculate SumRhoKProds for that PrevInf and PhaseSeverity for all countries in WhichCountries. 
				Calc_SumRhoK_country_PhaseSeverity_PrevInf(PARAMS, HOUSE, WhichCountries[countryindex], PhaseSeverity, PrevInf);
			
				//// recalculate Kplus value
				if (!(HOUSE.ModelVariant == SIMPLE_NUMERICAL))	
					if (PrevInf > 0)
					{
						Calc_KplusValues_country_PhaseSeverity(PARAMS, HOUSE, WhichCountries[countryindex], PhaseSeverity);  //// this also covers K_SEROPOS SS_Ks and SS_VEs together

						///// and therefore VacHazLikes if necessary
						if (HOUSE.ModelVariant == K_SEROPOS && HOUSE.SeroSpecificEfficacies)
							CalcAndStoreSumIntVacHazOverPatients(SeroPos, WhichCountries[countryindex], DATA, PARAMS, HOUSE, PhaseSeverity);
					}
			}
	}
	else	if (IsParamAWaningParam	(param_no, PARAMS.ParamNumbers))
	{
		int ParamType			= Find_ParamType_From_ASWaning_Param(param_no, HOUSE, PARAMS.ParamNumbers); 
		int BaselineSeroStatus	= Find_SeroStatus_From_Waning_Param	(param_no, HOUSE, PARAMS.ParamNumbers, ParamType); 

		if (HOUSE.AS_Waning_Homogeneous && BaselineSeroStatus != 0)
			std::cerr << "AmendPARAMS error: HOUSE.AS_Waning_Homogeneous && BaselineSeroStatus != 0" << endl;
		
		if (HOUSE.HillWaning)
		{
				 if (IsParamAHillHalflife	(param_no, PARAMS.ParamNumbers)) PARAMS.Hill_Halflives	[BaselineSeroStatus] = New_Param;
			else if (IsParamAHillPower		(param_no, PARAMS.ParamNumbers)) PARAMS.Hill_Powers	[BaselineSeroStatus] = New_Param;
		}
		else PARAMS.WaningParams[BaselineSeroStatus][ParamType] = New_Param;

		//// recalculate coefficients if spline
		if (HOUSE.AS_Waning == Age_Option::SPLINE)
			for (int poly = 0; poly < HOUSE.PolynomialsPerSpline_WaningDuration; poly++)
				CoefficientsFromKnots(PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus], PARAMS.WaningParams[BaselineSeroStatus], poly, PARAMS.Age_SplineCoeffs_DurationRate[BaselineSeroStatus][poly]);
		else if (HOUSE.AS_Waning == Age_Option::SPLINE_LINE)  ///// no need for anything for SPLINE_STEP as spline calculated purely from knots and no coefficients needed. 
			for (int poly = 0; poly < HOUSE.PolynomialsPerSpline_WaningDuration; poly++)
				CoefficientsFromKnots_Line(PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus], PARAMS.WaningParams[BaselineSeroStatus], poly, PARAMS.Age_SplineCoeffs_DurationRate[BaselineSeroStatus][poly]);
		else if (HOUSE.AS_Waning == Age_Option::CUBIC)  ///// no need for anything for SPLINE_STEP as spline calculated purely from knots and no coefficients needed. 
			for (int poly = 0; poly < HOUSE.PolynomialsPerSpline_WaningDuration; poly++)
				CoefficientsFromKnots_Cubic(PARAMS.Age_DurationRate_xKnots[BaselineSeroStatus], PARAMS.WaningParams[BaselineSeroStatus], poly, PARAMS.Age_SplineCoeffs_DurationRate[BaselineSeroStatus][poly]);

		Calc_WaningValues(PARAMS, HOUSE, BaselineSeroStatus, DATA.N_WaningDays);

		//// if (HOUSE.AS_Waning_Homogeneous), Pointers of WaningParams and WaningMults should be the same for SeroNeg and SeroPos. 
		//// so even thought preceding code will only change derived quantities for 0 (== SeroNeg), in effect change happening for both serostatuses. 
		//// Howevere, you would still need to change the vaccine hazards for both serostatuses. 

		//// Amend individual vaccine hazards (for each baseline serostatus and trial phase
		for (int countryindex = 0; countryindex < HOUSE.NoCountriesToFit; countryindex++)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
				if (HOUSE.AS_Waning_Homogeneous)
					for (int BaseSero = 0; BaseSero < HOUSE.HowManySeroStatuses; BaseSero++)
						Calc_IVH_Values_c_BS_PS(HOUSE.WhichCountries[countryindex], BaseSero, PhaseSeverity, PARAMS, DATA, HOUSE);
				else 
					Calc_IVH_Values_c_BS_PS(HOUSE.WhichCountries[countryindex], BaselineSeroStatus, PhaseSeverity, PARAMS, DATA, HOUSE);

		//// recalculate summary integrated vaccine hazards. 
		for (int countryindex = 0; countryindex < HOUSE.NoCountriesToFit; countryindex++)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++) //// this takes account of hospital severe as well if you're doing that. 
				if (HOUSE.AS_Waning_Homogeneous)
					for (int BaseSero = 0; BaseSero < HOUSE.HowManySeroStatuses; BaseSero++)
						CalcAndStoreSumIntVacHazOverPatients(BaseSero, HOUSE.WhichCountries[countryindex], DATA, PARAMS, HOUSE, PhaseSeverity);
				else
					CalcAndStoreSumIntVacHazOverPatients(BaselineSeroStatus, HOUSE.WhichCountries[countryindex], DATA, PARAMS, HOUSE, PhaseSeverity);
	}
	else	if (IsParamAn_ASVE_Param(param_no, PARAMS.ParamNumbers))
	{
		int Type				= Find_Type_From_ASVE_Param			(param_no, HOUSE, PARAMS.ParamNumbers); //// either half life or power or prop. (Or if !HOUSE.Hill_ASVE, then an age group multiplier)
		int BaselineSeroStatus	= Find_SeroStatus_FromASVE_Param	(param_no, HOUSE, PARAMS.ParamNumbers, Type);

		if (HOUSE.ASVE_OnlyOneSeroStatus) BaselineSeroStatus = HOUSE.ASVE_BS; 

		if (HOUSE.AS_VE_Homogeneous && BaselineSeroStatus != 0)
			std::cerr << "AmendPARAMS error: HOUSE.AS_VE_Homogeneous && BaselineSeroStatus != 0" << endl; 

		//// add to ASVE_Params
		PARAMS.ASVE_Params[BaselineSeroStatus][Type] = New_Param;

		//// recalculate coefficients if spline
		if (HOUSE.ASVE == Age_Option::SPLINE)
			for (int poly = 0; poly < HOUSE.PolynomialsPerSpline_EffMultiplier; poly++)
				CoefficientsFromKnots(PARAMS.Age_Effs_xKnots[BaselineSeroStatus], PARAMS.ASVE_Params[BaselineSeroStatus], poly, PARAMS.Age_SplineCoeffs_Effs[BaselineSeroStatus][poly]);
		else if (HOUSE.ASVE == Age_Option::SPLINE_LINE)
			for (int poly = 0; poly < HOUSE.PolynomialsPerSpline_EffMultiplier; poly++)
				CoefficientsFromKnots_Line(PARAMS.Age_Effs_xKnots[BaselineSeroStatus], PARAMS.ASVE_Params[BaselineSeroStatus], poly, PARAMS.Age_SplineCoeffs_Effs[BaselineSeroStatus][poly]);
		else if (HOUSE.ASVE == Age_Option::CUBIC)
			for (int poly = 0; poly < HOUSE.PolynomialsPerSpline_EffMultiplier; poly++)
				CoefficientsFromKnots_Cubic(PARAMS.Age_Effs_xKnots[BaselineSeroStatus], PARAMS.ASVE_Params[BaselineSeroStatus], poly, PARAMS.Age_SplineCoeffs_Effs[BaselineSeroStatus][poly]);

		//// Change AgeVaccEffMults
		Calc_AgeVaccEffMults(PARAMS, HOUSE, BaselineSeroStatus);

		//// SumRhoEffNegs - must be done AFTER AgeVaccEffMults
		if (HOUSE.EffNegWane == EffNegWane_Option::FROM_ZERO || HOUSE.EffNegWane == EffNegWane_Option::NO_WANE)
			Calc_SumRhoEffNegs_BaselineSeroStatus(BaselineSeroStatus, PARAMS, HOUSE); 

		//// Amend VacHazLikes
		for (int countryindex = 0; countryindex < HOUSE.NoCountriesToFit; countryindex++)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++) //// this takes account of hospital severe as well if you're doing that. 
				if (HOUSE.AS_VE_Homogeneous)
					for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
						CalcAndStoreSumIntVacHazOverPatients(BaselineSeroStatus, HOUSE.WhichCountries[countryindex], DATA, PARAMS, HOUSE, PhaseSeverity);
				else
					CalcAndStoreSumIntVacHazOverPatients(BaselineSeroStatus, HOUSE.WhichCountries[countryindex], DATA, PARAMS, HOUSE, PhaseSeverity);
	}
	else	if (IsParamAn_ASHaz_Param(param_no, PARAMS.ParamNumbers))
	{
		int Type				= Find_Type_From_ASHaz_Param			(param_no, HOUSE, PARAMS.ParamNumbers); //// either half life or power or prop. 

		//// add to ASHaz_Params
		PARAMS.ASHaz_Params[Type] = New_Param; 

		//// recalculate coefficients if spline
		if (HOUSE.AS_Haz == Age_Option::SPLINE)
			for (int poly = 0; poly < HOUSE.PolynomialsPerSpline_HazMultiplier; poly++)
				CoefficientsFromKnots(PARAMS.Age_HazMult_xKnots, PARAMS.ASHaz_Params, poly, PARAMS.Age_SplineCoeffs_HazMult[poly]);
		else if (HOUSE.AS_Haz == Age_Option::SPLINE_LINE)
			for (int poly = 0; poly < HOUSE.PolynomialsPerSpline_HazMultiplier; poly++)
				CoefficientsFromKnots_Line(PARAMS.Age_HazMult_xKnots, PARAMS.ASHaz_Params, poly, PARAMS.Age_SplineCoeffs_HazMult[poly]);
		else if (HOUSE.AS_Haz == Age_Option::CUBIC)
			for (int poly = 0; poly < HOUSE.PolynomialsPerSpline_HazMultiplier; poly++)
				CoefficientsFromKnots_Cubic(PARAMS.Age_HazMult_xKnots, PARAMS.ASHaz_Params, poly, PARAMS.Age_SplineCoeffs_HazMult[poly]);

		//// Change Age_HazMults
		Calc_Age_HazMults(PARAMS, HOUSE);

		//// Amend VacHazLikes (for both serostatuses)
		for (int countryindex = 0; countryindex < HOUSE.NoCountriesToFit; countryindex++)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++) //// this takes account of hospital severe as well if you're doing that. 
				for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
					CalcAndStoreSumIntVacHazOverPatients(BaselineSeroStatus, HOUSE.WhichCountries[countryindex], DATA, PARAMS, HOUSE, PhaseSeverity);
	}
	else	if (IsParamAKnot		(param_no, PARAMS.ParamNumbers)) 
	{
		int knot	= Find_Knot_FromKnotParam	(param_no, HOUSE, PARAMS.ParamNumbers);
		int country = Find_Country_FromKnotParam(param_no, HOUSE, PARAMS.ParamNumbers, knot	);

		PARAMS.yKnots[country][knot] = New_Param;

		//// Amend SplineCoeffs
		for (int poly = 0; poly < HOUSE.PolynomialsPerCountry; poly++)
			CoefficientsFromKnots(	PARAMS.xKnots[country], PARAMS.yKnots[country], poly, PARAMS.SplineCoeffs[country][poly]);

		//// Change Baseline Hazard Values. 
		Calc_BaseHazValues(country, PARAMS, DATA, HOUSE); 

		//// Change Integrated Hazard Values
		Calc_IntBaseHazValues(country, PARAMS, DATA, HOUSE); 

		//// Amend individual vaccine hazards (for each baseline serostatus and trial phase
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
				Calc_IVH_Values_c_BS_PS(country, BaselineSeroStatus, PhaseSeverity, PARAMS, DATA, HOUSE);

		//// Amend summary level vaccine hazards (for each baseline serostatus and trial phase
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
				CalcAndStoreSumIntVacHazOverPatients(BaselineSeroStatus, country, DATA, PARAMS, HOUSE, PhaseSeverity);
	}
	else	if (IsParamAHistHaz		(param_no, PARAMS.ParamNumbers)) 
	{
		int country = Find_Country_FromHistHazParam(param_no, PARAMS.ParamNumbers);

		if (HOUSE.ModelVariant != SIMPLE_NUMERICAL)
		{
			Calc_KplusValues_country(PARAMS, HOUSE, country);
			if (HOUSE.ModelVariant == AS_PRIME)
			{
				for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
					for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
						Calc_K_Primes_BS_PhaseSeverity_country_serotype(PARAMS, HOUSE, SeroPos, PhaseSeverity, country, serotype); 

				Calc_KplusPrimes_country(PARAMS, HOUSE, country);
			}
		}

		if (HOUSE.ModelVariant == K_SEROPOS || HOUSE.ModelVariant == AS_PRIME)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
				CalcAndStoreSumIntVacHazOverPatients(SeroPos, country, DATA, PARAMS, HOUSE, PhaseSeverity);

		if (!HOUSE.Empirical_SeroPrevs)
			Calc_ParamSeroPrevs_country(country, PARAMS, HOUSE);
	}
	else	if (IsParamA_qval		(param_no, PARAMS.ParamNumbers))
	{
		int qval	= Find_qval_From_qParam		(param_no, HOUSE, PARAMS.ParamNumbers);
		int country = Find_Country_From_qParam(param_no, HOUSE, PARAMS.ParamNumbers, qval);

		Calc_RhosFrom_qParams	(country, PARAMS, HOUSE);															   
		Calc_SumRhoEffs_country(PARAMS, HOUSE, country);

		//// SumRhoEffNegs - must be done AFTER Rhos
		if (HOUSE.EffNegWane == EffNegWane_Option::FROM_ZERO || HOUSE.EffNegWane == EffNegWane_Option::NO_WANE)
			Calc_SumRhoEffNegs_country(country, PARAMS, HOUSE);

		//// recalculate SumRhoKProds for that PrevInf and PhaseSeverity for all countries. 
		if (HOUSE.SeroSpecific_K_values)	//// recalculate Kplus value
		{
			Calc_SumRhoKs_country		(PARAMS, HOUSE, country);
			if (HOUSE.ModelVariant != SIMPLE_NUMERICAL)
				Calc_KplusValues_country	(PARAMS, HOUSE, country); //// this also covers K_SEROPOS SS_Ks and SS_VEs together

			if (HOUSE.ModelVariant == K_SEROPOS)
				for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
					CalcAndStoreSumIntVacHazOverPatients(SeroPos, country, DATA, PARAMS, HOUSE, PhaseSeverity);  //// this also covers K_SEROPOS SS_Ks and SS_VEs together
		}
		if (HOUSE.SeroSpecificEfficacies && HOUSE.ASVE != Age_Option::INDEPENDENT && HOUSE.SSASVE_Additive)
			for (int BS = 0; BS < HOUSE.HowManySeroStatuses; BS++)
				for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
					CalcAndStoreSumIntVacHazOverPatients(BS, country, DATA, PARAMS, HOUSE, PhaseSeverity);
	}
	else	if (IsParamA_rho		(param_no, PARAMS.ParamNumbers))
	{
		int serotype	= Find_rho_From_rhoParam	(param_no, PARAMS.ParamNumbers.Min_Rho, HOUSE.N_STypes);
		int country		= Find_Country_From_rhoParam(param_no, PARAMS.ParamNumbers.Min_Rho, HOUSE.N_STypes, serotype);

		PARAMS.rhos[country][serotype] = New_Param;

		Calc_RhosFrom_qParams(country, PARAMS, HOUSE);
		Calc_SumRhoEffs_country(PARAMS, HOUSE, country);	//// used even when not doing sero-specific efficacies (in augmentation and survival). 
																
		//// SumRhoEffNegs - must be done AFTER Rhos
		if (HOUSE.EffNegWane == EffNegWane_Option::FROM_ZERO || HOUSE.EffNegWane == EffNegWane_Option::NO_WANE)
			Calc_SumRhoEffNegs_country(country, PARAMS, HOUSE);

		//// recalculate SumRhoKProds for that PrevInf and PhaseSeverity for all countries. 
		if (HOUSE.SeroSpecific_K_values)	//// recalculate Kplus value
		{
			Calc_SumRhoKs_country(PARAMS, HOUSE, country);
			if (HOUSE.ModelVariant != SIMPLE_NUMERICAL)
				Calc_KplusValues_country(PARAMS, HOUSE, country); //// this also covers K_SEROPOS SS_Ks and SS_VEs together

			if (HOUSE.ModelVariant == K_SEROPOS)
				for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
					CalcAndStoreSumIntVacHazOverPatients(SeroPos, country, DATA, PARAMS, HOUSE, PhaseSeverity); //// this also covers K_SEROPOS SS_Ks and SS_VEs together
		}
		if (HOUSE.SeroSpecificEfficacies && HOUSE.ASVE != Age_Option::INDEPENDENT && HOUSE.SSASVE_Additive)
			for (int BS = 0; BS < HOUSE.HowManySeroStatuses; BS++)
				for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
					CalcAndStoreSumIntVacHazOverPatients(BS, country, DATA, PARAMS, HOUSE, PhaseSeverity);
	}
	else	if (IsParamAnEfficacy	(param_no, PARAMS.ParamNumbers) || IsParamAn_atInf_Efficacy(param_no, PARAMS.ParamNumbers))
	{
		int BaselineSeroStatus	= Find_SeroStatus_FromEffParam		(param_no, HOUSE, PARAMS.ParamNumbers);
		int SeroType			= Find_Serotype_FromEffParam		(param_no, HOUSE, PARAMS.ParamNumbers, BaselineSeroStatus);
		int PhaseSeverity		= Find_PhaseSeverity_From_EffParam	(param_no, HOUSE, PARAMS.ParamNumbers, BaselineSeroStatus, SeroType);

		if (!HOUSE.ResidEffs) //// PARAMS.Efficacies is what it describes (i.e. default)
			PARAMS.Efficacies[PhaseSeverity][SeroType][BaselineSeroStatus] = New_Param;
		else  //// PARAMS.Efficacies actually refers to (VE_initial - VE_atInf), so you need to be more careful. 
		{
			if (IsParamAnEfficacy(param_no, PARAMS.ParamNumbers)) //// proposed parameter is VE_initial. 
				PARAMS.Efficacies[PhaseSeverity][SeroType][BaselineSeroStatus] = New_Param - PARAMS.Inf_Effs[PhaseSeverity][SeroType][BaselineSeroStatus];
			else if (IsParamAn_atInf_Efficacy(param_no, PARAMS.ParamNumbers)) //// Init_Effs not stored, so need to be a little more clever. Your proposed parameter is VE_atInf 
			{
				//// order is important here. 
				PARAMS.Efficacies	[PhaseSeverity][SeroType][BaselineSeroStatus] += PARAMS.Inf_Effs[PhaseSeverity][SeroType][BaselineSeroStatus]; //// i.e. Efficacy now equal to VE_initial. 
				PARAMS.Inf_Effs		[PhaseSeverity][SeroType][BaselineSeroStatus] = New_Param;
				PARAMS.Efficacies	[PhaseSeverity][SeroType][BaselineSeroStatus] -= PARAMS.Inf_Effs[PhaseSeverity][SeroType][BaselineSeroStatus]; //// Efficacyrefers to (VE_initial - VE_atInf) again. 
			}
		}
		Calc_SumRhoEffs_BaselineSeroStatus	(PARAMS, HOUSE, BaselineSeroStatus);
		ChooseMaxMinPhaseSeverities			(HOUSE, PhaseSeverity); //// changes HOUSE.MinMaxPhaseSeveritiesToLoopOver (potentially)

		//// SumRhoEffNegs - must be done AFTER Efficacies
		if (HOUSE.EffNegWane == EffNegWane_Option::FROM_ZERO || HOUSE.EffNegWane == EffNegWane_Option::NO_WANE)
			for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++)
				Calc_SumRhoEffNegs_c_BS_PS(HOUSE.WhichCountries[countryindex], PhaseSeverity, BaselineSeroStatus, PARAMS, HOUSE);

		bool Scenario_1 = HOUSE.SeroSpecificEfficacies && (HOUSE.ModelVariant == K_SEROPOS) && HOUSE.SeroSpecific_K_values && (BaselineSeroStatus == SeroPos);
		bool Scenario_2 = HOUSE.AdditiveSSASVEs; //// true for either SSKs or not. 

		if (Scenario_1 || Scenario_2) //// in all other scenarios, efficacy doesn't affect VacHazLikes
		{
			if (Scenario_1) Calc_KplusValues_All(PARAMS, HOUSE);
			for (int PhaseSeverity = HOUSE.MinMaxPhaseSeveritiesToLoopOver[0]; PhaseSeverity < HOUSE.MinMaxPhaseSeveritiesToLoopOver[1]; PhaseSeverity++)
				for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++)
					CalcAndStoreSumIntVacHazOverPatients(BaselineSeroStatus, HOUSE.WhichCountries[countryindex], DATA, PARAMS, HOUSE, PhaseSeverity);
		}
	}
	else	if (IsParamAn_AS_Prime_Param(param_no, PARAMS.ParamNumbers))
	{
		int Type				= Find_ParamType_From_AS_Prime_Param	(param_no, HOUSE, PARAMS.ParamNumbers); //// either half life or coeff. 
		int BaselineSeroStatus	= Find_SeroStatus_From_AS_Prime_Param	(param_no, HOUSE, PARAMS.ParamNumbers, Type);

		PARAMS.AS_Priming_Params[BaselineSeroStatus][Type] = New_Param; 

		//// Ks_Prime
		Calc_K_Primes_BS(PARAMS, HOUSE, BaselineSeroStatus);

		//// KPlusPrimeValues and SumRhoK0_SNeg_Primes: Both consider SS_Ks, therefore has rhos. 
		if (BaselineSeroStatus == SeroNeg)
			Calc_SumRhoK0_SNeg_Primes_ALL(PARAMS, HOUSE);
		else if (BaselineSeroStatus == SeroPos)
			Calc_KplusPrimes_All(PARAMS, HOUSE);

		//// VacHazLikes, all countries and phase severities for that serostatus. 
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++)
				CalcAndStoreSumIntVacHazOverPatients(BaselineSeroStatus, HOUSE.WhichCountries[countryindex], DATA, PARAMS, HOUSE, PhaseSeverity);
	}	
	else	if (param_no == PARAMS.ParamNumbers.BS_BaseHazMult)
	{
		if (HOUSE.AdjHaz && !HOUSE.Skip_BS_BaseHazMult)
		{
			//// populate
			PARAMS.BS_BaseHazMults[HOUSE.Which_BS_BaseHazMult] = New_Param;

			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
				for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++)
					CalcAndStoreSumIntVacHazOverPatients(HOUSE.Which_BS_BaseHazMult, HOUSE.WhichCountries[countryindex], DATA, PARAMS, HOUSE, PhaseSeverity); // only recalculate for the appropriate serostatus.
		}
	}
	else	std::cerr << "AmendPARAMS error: param_no not recognized" << endl;

#ifdef PRINT_WITHIN_CHANGE_FUNCTIONS
#ifdef PRINT_WITHIN_CHANGE_KNOT_FUNCTIONS
	std::cout << "AmendParams p" << param_no << " " << PARAMS.NamesOfParameters[param_no] << " value " << New_Param << endl;
#else 
	if (!IsParamAKnot(param_no, PARAMS.ParamNumbers)) std::cout << " DONE" << endl;
#endif 
#endif
}
void AmendParams				(std::vector<DType> ParamVec, Params_Struct &PARAMS, const DATA_struct &DATA, const Housekeeping_Struct &HOUSE)  //// does all parameters (used for mean and modal posterior say, or input from previous chain. )
{
	//// overload of single param AmendParams. Assumes that i) ParamVec.size() == HOUSE.No_Parameters. 

	if (ParamVec.size() != HOUSE.No_Parameters) 
		for (int error = 0; error < 20; error++) 
			std::cerr << "AmendParams overload 1 error: ParamVec.size() = " << ParamVec.size() << " HOUSE.No_Parameters = " << HOUSE.No_Parameters << endl; 

	//// apply AmendParams repeatedly. 
	for (int param_no = 0; param_no < HOUSE.No_Parameters; param_no++)
		AmendParams(param_no, ParamVec[param_no], PARAMS, DATA, HOUSE); 
}
void AmendParams				(std::vector<DType> ParamVec, std::vector<int> PVecIndices, Params_Struct &PARAMS, const DATA_struct &DATA, const Housekeeping_Struct &HOUSE) //// does a selction of parameters (used for simulataneous updates)
{
	//// overload of single param AmendParams. Assumes that i) ParamVec.size() == PVecIndices.size(). 
	if (ParamVec.size() != HOUSE.No_Parameters) 
		for (int error = 0; error < 20; error++) 
			std::cerr << "AmendParams overload 2 error: ParamVec.size() = " << ParamVec.size() << " HOUSE.No_Parameters = " << HOUSE.No_Parameters << endl;
	//// apply AmendParams repeatedly. 
	for (int param_no_index = 0; param_no_index < PVecIndices.size(); param_no_index++)
		AmendParams(PVecIndices[param_no_index], ParamVec[param_no_index], PARAMS, DATA, HOUSE);
}
void UpdateOrReset_ParamsEtc	(int &param_no, const DATA_struct &DATA, Params_Struct &ParamsToUpdate, const Params_Struct &CorrectPARAMS, const Housekeeping_Struct &HOUSE) ///  //// because things are now a function of PARAMS as a whole, you need to reset more things than when objects like lookup tables were floating around unstructured. 
{
#ifdef PRINT_WITHIN_CHANGE_FUNCTIONS
#ifdef PRINT_WITHIN_CHANGE_KNOT_FUNCTIONS
	std::cout << "UpdateOrReset_ParamsEtc p" << param_no << " " << PARAMS.NamesOfParameters[param_no] << " value " << New_Param << endl;
#else 
	if (!IsParamAKnot(param_no, ParamsToUpdate.ParamNumbers)) std::cout << "UpdateOrReset_ParamsEtc p" << param_no << " " << ParamsToUpdate.NamesOfParameters[param_no] << endl;
#endif 
#endif

#ifdef COPY_ENTIRE_PARAM_STRUCTS
	ParamsToUpdate.SetEqualToAnother(CorrectPARAMS, HOUSE, DATA);
#else 
	ParamsToUpdate.ParamVec[param_no] = CorrectPARAMS.ParamVec[param_no]; //// All parameters require this change.

	//// Copy Likelihood array. 
	for (int country = 0; country < HOUSE.TotalCountries; country++)
		for (int component = 0; component < HOUSE.LCPerC; component++)
			ParamsToUpdate.LikeParts[country][component] = CorrectPARAMS.LikeParts[country][component];
	
	//// Copy Full Likelihood (and posterior). And Priors. 
	ParamsToUpdate.LikeFull							= CorrectPARAMS.LikeFull;
	ParamsToUpdate.LogPosterior						= CorrectPARAMS.LogPosterior;
	ParamsToUpdate.LogPriorFull						= CorrectPARAMS.LogPriorFull;
	ParamsToUpdate.Total_logLike					= CorrectPARAMS.Total_logLike;
	ParamsToUpdate.ParamArray_logPriors[param_no]	= CorrectPARAMS.ParamArray_logPriors[param_no];

	//// **** //// **** //// **** //// **** //// **** //// **** //// **** 
	//// **** Parameter values and associated quantities (e.g. relative risks / precalculated integrated hazards etc.)

			if (IsParamARelativeRisk(param_no, CorrectPARAMS.ParamNumbers) || IsParamARelativeRisk_Hosp(param_no, CorrectPARAMS.ParamNumbers))
	{
		int PrevInf			= Find_PrevInf_From_K_Param			(param_no, HOUSE, CorrectPARAMS.ParamNumbers);
		int PhaseSeverity	= Find_PhaseSeverity_From_K_Param	(param_no, HOUSE, CorrectPARAMS.ParamNumbers, PrevInf);
		int serotype		= Find_Serotype_From_K_Param		(param_no, HOUSE, CorrectPARAMS.ParamNumbers, PrevInf, PhaseSeverity);;
		std::vector<int> WhichCountries; 

		if (IsParamARelativeRisk(param_no, CorrectPARAMS.ParamNumbers))
			if (HOUSE.ModellingHospitalized && PhaseSeverity == PassiveSevere && HOUSE.ModelHosp_Indie_Ks)	WhichCountries = HOUSE.CYD_14_countries;	//// only do CYD-14 countries. 
			else																							WhichCountries = HOUSE.WhichCountries;		//// otherwise do all countries
		else if (IsParamARelativeRisk_Hosp(param_no, CorrectPARAMS.ParamNumbers))											WhichCountries = HOUSE.CYD_15_countries;

		int Min_PhaseSeverity = PhaseSeverity, Max_PhaseSeverity = PhaseSeverity + 1; //// i.e. default is that you do one phase severity at a time in loop below. 
		if (HOUSE.PS_Ks_Multiply_AM_Ks && PhaseSeverity == ActiveMild) { Min_PhaseSeverity = ActiveMild; Max_PhaseSeverity = HOUSE.HowManyCaseCategories; } //// i.e. if changing Passive only this doesn't affect active, so default definition okay. 

		for (int PhaseSeverity = Min_PhaseSeverity; PhaseSeverity < Max_PhaseSeverity; PhaseSeverity++)
		{
			//// K_s array. 
			for (int countryindex = 0; countryindex < WhichCountries.size(); countryindex++)
				if (HOUSE.SSKs_FitMatrix[PhaseSeverity][PrevInf]) //// i.e. if fitting all serotypes for this K_PhaseSeverity_PrevInf, then do the below one serotype (this one) at a time. Automatically true for !SSKs too as default values of matrix all true (and changed only if SSKs). 
					ParamsToUpdate.K_s[WhichCountries[countryindex]][PhaseSeverity][PrevInf][serotype] = CorrectPARAMS.K_s[WhichCountries[countryindex]][PhaseSeverity][PrevInf][serotype];
				else //// so SSKs == true & HOUSE.SSKs_FitMatrix[PhaseSeverity][PrevInf] == false. i.e. this K_PhaseSeverity_PrevInf is not serotype specific. Hence you will amend all serotype Ks for this PhaseSeverity and PrevInf at the same time. 
					for (int serotypeDummy = 0; serotypeDummy < HOUSE.N_STypes_Ks; serotypeDummy++)
						ParamsToUpdate.K_s[WhichCountries[countryindex]][PhaseSeverity][PrevInf][serotypeDummy] = CorrectPARAMS.K_s[WhichCountries[countryindex]][PhaseSeverity][PrevInf][serotypeDummy];

			//// then update the hospital K's (as Hosp_K_CYD15 = regular hosp K times Hosp K multiplier. 
			if (HOUSE.ModellingHospitalized && PhaseSeverity == PassiveSevere && !HOUSE.ModelHosp_Indie_Ks) Calc_Hospital_Ks(ParamsToUpdate, HOUSE);

			if (HOUSE.ModelVariant == AS_PRIME && PrevInf > 0)
			{
				//// Ks_Prime
				for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++) //// as PrevInf > 0, looks like you'll only need SeroPos, but K1 and K2 affect K0prime,- and K1prime,+
					for (int Age = 0; Age < HOUSE.HowManyAges; Age++)
						//ParamsToUpdate.Ks_Prime[WhichCountries[0]][PhaseSeverity][BaselineSeroStatus][serotype][Age] = CorrectPARAMS.Ks_Prime[WhichCountries[0]][PhaseSeverity][BaselineSeroStatus][serotype][Age];
						for (int countryindex = 0; countryindex < WhichCountries.size(); countryindex++)
							ParamsToUpdate.Ks_Prime[WhichCountries[countryindex]][PhaseSeverity][BaselineSeroStatus][serotype][Age] = CorrectPARAMS.Ks_Prime[WhichCountries[countryindex]][PhaseSeverity][BaselineSeroStatus][serotype][Age];

				//// KPlusPrimeValues and SumRhoK0_SNeg_Primes: Both consider SS_Ks, therefore has rhos. Therefore need to do all countries for this K, even if Ks same between countries. Note you don't store arrays of KplusValues or KplusPrimes for single/case serotypes, only "survivors". 
				for (int countryindex = 0; countryindex < WhichCountries.size(); countryindex++) //// have to do all countries here as
					for (int Age = 0; Age < HOUSE.HowManyAges; Age++)
					{
						ParamsToUpdate.KPlusPrimeValues		[WhichCountries[countryindex]][Age][PhaseSeverity]	= CorrectPARAMS.KPlusPrimeValues	[WhichCountries[countryindex]][Age][PhaseSeverity];
						ParamsToUpdate.SumRhoK0_SNeg_Primes	[WhichCountries[countryindex]][PhaseSeverity][Age]	= CorrectPARAMS.SumRhoK0_SNeg_Primes[WhichCountries[countryindex]][PhaseSeverity][Age];
					}
			}

			for (int countryindex = 0; countryindex < WhichCountries.size(); countryindex++)
			{
				//// update Kplus values 
				if (!(HOUSE.ModelVariant == SIMPLE_NUMERICAL))	
					if (PrevInf > 0)
					{
						for (int age = 0; age < HOUSE.HowManyAges; age++)
							ParamsToUpdate.KplusValues[WhichCountries[countryindex]][age][PhaseSeverity] = CorrectPARAMS.KplusValues[WhichCountries[countryindex]][age][PhaseSeverity];

						///// and therefore VacHazLikes if necessary
						if (HOUSE.ModelVariant == K_SEROPOS || HOUSE.ModelVariant == AS_PRIME)
							ParamsToUpdate.VacHazLikes[WhichCountries[countryindex]][SeroPos][PhaseSeverity] = CorrectPARAMS.VacHazLikes[WhichCountries[countryindex]][SeroPos][PhaseSeverity];
						if (HOUSE.ModelVariant == AS_PRIME)
							ParamsToUpdate.VacHazLikes[WhichCountries[countryindex]][SeroNeg][PhaseSeverity] = CorrectPARAMS.VacHazLikes[WhichCountries[countryindex]][SeroNeg][PhaseSeverity];

						if ((HOUSE.ModelVariant == K_SEROPOS) && HOUSE.SeroSpecific_K_values && HOUSE.SeroSpecificEfficacies) //// Do the same for KplusValues with Effs. Don't need last condition in if statement as PrevInf > 0; 
							for (int age = 0; age < HOUSE.HowManyAges; age++)
								ParamsToUpdate.Meta_KplusValues[With_Effs][WhichCountries[countryindex]][age][PhaseSeverity] = CorrectPARAMS.Meta_KplusValues[With_Effs][WhichCountries[countryindex]][age][PhaseSeverity];
					}
				//// update SumRhoKs for all countries
				ParamsToUpdate.SumRhoKs[WhichCountries[countryindex]][PhaseSeverity][PrevInf] = CorrectPARAMS.SumRhoKs[WhichCountries[countryindex]][PhaseSeverity][PrevInf];
			}
		}
	}
	else	if (param_no == CorrectPARAMS.ParamNumbers.Hosp_K_multiplier) //// will only be true if (HOUSE.ModellingHospitalized !HOUSE.ModelHosp_Indie_Ks) 
	{
		ParamsToUpdate.Hosp_K_mult = CorrectPARAMS.Hosp_K_mult;
	
		Calc_Hospital_Ks(ParamsToUpdate, HOUSE);

		std::vector<int> WhichCountries = HOUSE.CYD_15_countries; 
		int PhaseSeverity = PassiveSevere; 
		for (int countryindex = 0; countryindex < WhichCountries.size(); countryindex++)
			for (int PrevInf = 0; PrevInf < HOUSE.Num_K_Params; PrevInf++)
			{
				//// update SumRhoKs for all countries
				ParamsToUpdate.SumRhoKs[WhichCountries[countryindex]][PhaseSeverity][PrevInf] = CorrectPARAMS.SumRhoKs[WhichCountries[countryindex]][PhaseSeverity][PrevInf];

				//// recalculate Kplus value
				if (!(HOUSE.ModelVariant == SIMPLE_NUMERICAL))	
					if (PrevInf > 0)
					{
						for (int age = 0; age < HOUSE.HowManyAges; age++)
							ParamsToUpdate.KplusValues[WhichCountries[countryindex]][age][PhaseSeverity] = CorrectPARAMS.KplusValues[WhichCountries[countryindex]][age][PhaseSeverity];

						///// and therefore VacHazLikes if necessary
						if ((HOUSE.ModelVariant == K_SEROPOS) && HOUSE.SeroSpecificEfficacies)
							ParamsToUpdate.VacHazLikes[WhichCountries[countryindex]][SeroPos][PhaseSeverity] = CorrectPARAMS.VacHazLikes[WhichCountries[countryindex]][SeroPos][PhaseSeverity];

						if ((HOUSE.ModelVariant == K_SEROPOS) && HOUSE.SeroSpecific_K_values && HOUSE.SeroSpecificEfficacies /*&& BaselineSeroStatus == SeroPos*/) //// Do the same for KplusValues with Effs. Don't need last condition in if statement as PrevInf > 0; 
							for (int age = 0; age < HOUSE.HowManyAges; age++)
								ParamsToUpdate.Meta_KplusValues[With_Effs][WhichCountries[countryindex]][age][PhaseSeverity] = CorrectPARAMS.Meta_KplusValues[With_Effs][WhichCountries[countryindex]][age][PhaseSeverity];
					}
			}
	}
	else	if (IsParamAWaningParam	(param_no, CorrectPARAMS.ParamNumbers)	)
	{
		int ParamType			= Find_ParamType_From_ASWaning_Param(param_no, HOUSE, CorrectPARAMS.ParamNumbers); 
		int BaselineSeroStatus	= Find_SeroStatus_From_Waning_Param	(param_no, HOUSE, CorrectPARAMS.ParamNumbers, ParamType); 

		if (HOUSE.AS_Waning_Homogeneous && BaselineSeroStatus != 0)
			std::cerr << "UpdateOrReset_ParamsEtc error: HOUSE.AS_Waning_Homogeneous && BaselineSeroStatus != 0" << endl;

		if (HOUSE.HillWaning)
		{
				 if (IsParamAHillHalflife	(param_no, CorrectPARAMS.ParamNumbers)) 	ParamsToUpdate.Hill_Halflives	[BaselineSeroStatus] = CorrectPARAMS.Hill_Halflives	[BaselineSeroStatus];
			else if (IsParamAHillPower		(param_no, CorrectPARAMS.ParamNumbers)) 	ParamsToUpdate.Hill_Powers		[BaselineSeroStatus] = CorrectPARAMS.Hill_Powers	[BaselineSeroStatus];
		}
		else	ParamsToUpdate.WaningParams[BaselineSeroStatus][ParamType] = CorrectPARAMS.WaningParams[BaselineSeroStatus][ParamType];

		//// recalculate coefficients (xKnots don't change, and ASVE_Params are the yKnots)

		if (HOUSE.AS_Waning == Age_Option::SPLINE || HOUSE.AS_Waning == Age_Option::SPLINE_LINE || HOUSE.AS_Waning == Age_Option::CUBIC) ///// no need for anything for SPLINE_STEP as spline calculated purely from knots and no cofficents needed. 
			for (int poly = 0; poly < HOUSE.PolynomialsPerSpline_WaningDuration; poly++)
				for (int coeff = 0; coeff < HOUSE.MaxSplineDegree_WaningDuration; coeff++)
					ParamsToUpdate.Age_SplineCoeffs_DurationRate[BaselineSeroStatus][poly][coeff] = CorrectPARAMS.Age_SplineCoeffs_DurationRate[BaselineSeroStatus][poly][coeff];

		//// reset waning rates
		for (int AgeInYears = 0; AgeInYears < HOUSE.NumAges_Waning; AgeInYears++)
			for (int daypostdose = 0; daypostdose < DATA.N_WaningDays; daypostdose++)
				ParamsToUpdate.WaningMults[AgeInYears][BaselineSeroStatus][daypostdose] = CorrectPARAMS.WaningMults[AgeInYears][BaselineSeroStatus][daypostdose];

		//// Amend individual vaccine hazards (for each baseline serostatus and trial phase) and reset VacHazLikes
		int Stratum_Index = NULL, Max_BaseSero = NULL;
		if (HOUSE.AS_Waning_Homogeneous) Max_BaseSero = HOUSE.HowManySeroStatuses; else Max_BaseSero = BaselineSeroStatus + 1; //// i.e. if AS_Waning_Homogeneous do all serostatuses, otherwise do only one. 
		for (int BaseSero = BaselineSeroStatus; BaseSero < Max_BaseSero; BaseSero++)
		{
			Stratum_Index = HOUSE.Strata[VaccineGroup][BaseSero];
			for (int countryindex = 0; countryindex < HOUSE.NoCountriesToFit; countryindex++)
				for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
					for (int i = 0; i < DATA.Set_Array[HOUSE.WhichCountries[countryindex]][Stratum_Index].size(); i++)
					{
						int person_i = DATA.Set_Array[HOUSE.WhichCountries[countryindex]][Stratum_Index][i];
						ParamsToUpdate.IVH_vals[person_i][PhaseSeverity] = CorrectPARAMS.IVH_vals[person_i][PhaseSeverity]; 
					}

			//// reset VacHazLikes
			for (int country = 0; country < HOUSE.TotalCountries; country++)
				for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
					ParamsToUpdate.VacHazLikes[country][BaseSero][PhaseSeverity] = CorrectPARAMS.VacHazLikes[country][BaseSero][PhaseSeverity];
		}
	}
	else	if (IsParamAn_ASVE_Param(param_no, CorrectPARAMS.ParamNumbers))
	{
		int Type				= Find_Type_From_ASVE_Param			(param_no, HOUSE, CorrectPARAMS.ParamNumbers); //// either half life or power or prop. 
		int BaselineSeroStatus	= Find_SeroStatus_FromASVE_Param	(param_no, HOUSE, CorrectPARAMS.ParamNumbers, Type);
		/*std::cout << CorrectPARAMS.NamesOfParameters[param_no] << " BaselineSeroStatus " << BaselineSeroStatus << " Type " << Type  << endl;*/

		if (HOUSE.ASVE_OnlyOneSeroStatus) BaselineSeroStatus = HOUSE.ASVE_BS;

		if (HOUSE.AS_VE_Homogeneous && BaselineSeroStatus != 0)
			std::cerr << "UpdateOrReset_ParamsEtc error: HOUSE.AS_VE_Homogeneous && BaselineSeroStatus != 0" << endl; 

		//// add to ASVE_Params
		ParamsToUpdate.ASVE_Params[BaselineSeroStatus][Type] = CorrectPARAMS.ASVE_Params[BaselineSeroStatus][Type];

		//// Change AgeVaccEffMults
		for (int age = 0; age < HOUSE.HowManyAges; age++)
			ParamsToUpdate.AgeEff_Mult[BaselineSeroStatus][age] = CorrectPARAMS.AgeEff_Mult[BaselineSeroStatus][age];

		//// recalculate coefficients (xKnots don't change, and ASVE_Params are the yKnots)
		if (HOUSE.ASVE == Age_Option::SPLINE || HOUSE.ASVE == Age_Option::SPLINE_LINE || HOUSE.ASVE == Age_Option::CUBIC)  ///// no need for anything for SPLINE_STEP OR CATEGORICAL as spline calculated purely from knots and no cofficents needed. 
			for (int poly = 0; poly < HOUSE.PolynomialsPerSpline_EffMultiplier; poly++)
				for (int coeff = 0; coeff < HOUSE.MaxSplineDegree_EffMultiplier; coeff++)
					ParamsToUpdate.Age_SplineCoeffs_Effs[BaselineSeroStatus][poly][coeff] = CorrectPARAMS.Age_SplineCoeffs_Effs[BaselineSeroStatus][poly][coeff];

		//// SumRhoEffNegs
		if (HOUSE.EffNegWane == EffNegWane_Option::FROM_ZERO || HOUSE.EffNegWane == EffNegWane_Option::NO_WANE)
			for (int country = 0; country < HOUSE.TotalCountries; country++)
				for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++) //// this takes account of hospital severe as well if you're doing that. 
					for (int age = 0; age < HOUSE.HowManyAges; age++)
						ParamsToUpdate.SumRhoEffNegs[PhaseSeverity][country][age][BaselineSeroStatus] = CorrectPARAMS.SumRhoEffNegs[PhaseSeverity][country][age][BaselineSeroStatus];

		//// must reset VacHazLikes
		for (int country = 0; country < HOUSE.TotalCountries; country++)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++) //// this takes account of hospital severe as well if you're doing that. 
				if (HOUSE.AS_VE_Homogeneous)
					for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
						ParamsToUpdate.VacHazLikes[country][BaselineSeroStatus][PhaseSeverity] = CorrectPARAMS.VacHazLikes[country][BaselineSeroStatus][PhaseSeverity];
				else
					ParamsToUpdate.VacHazLikes[country][BaselineSeroStatus][PhaseSeverity] = CorrectPARAMS.VacHazLikes[country][BaselineSeroStatus][PhaseSeverity];
	}
	else	if (IsParamAn_ASHaz_Param(param_no, CorrectPARAMS.ParamNumbers))
	{
		int Type				= Find_Type_From_ASHaz_Param			(param_no, HOUSE, CorrectPARAMS.ParamNumbers); //// either half life or power or prop. 

		//// add to ASHaz
		ParamsToUpdate.ASHaz_Params[Type] = CorrectPARAMS.ASHaz_Params[Type];

		//// Change AgeHaz_Mult
		for (int age = 0; age < HOUSE.HowManyAges; age++)
			ParamsToUpdate.AgeHaz_Mult[age] = CorrectPARAMS.AgeHaz_Mult[age];

		//// recalculate coefficients (xKnots don't change, and AS_Haz are the yKnots)
		if (HOUSE.AS_Haz == Age_Option::SPLINE || HOUSE.AS_Haz == Age_Option::SPLINE_LINE || HOUSE.AS_Haz == Age_Option::CUBIC)
			for (int poly = 0; poly < HOUSE.PolynomialsPerSpline_HazMultiplier; poly++)
				for (int coeff = 0; coeff < HOUSE.MaxSplineDegree_HazMultiplier; coeff++)
					ParamsToUpdate.Age_SplineCoeffs_HazMult[poly][coeff] = CorrectPARAMS.Age_SplineCoeffs_HazMult[poly][coeff];

		//// must reset VacHazLikes
		for (int country = 0; country < HOUSE.TotalCountries; country++)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++) //// this takes account of hospital severe as well if you're doing that. 
				for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
					ParamsToUpdate.VacHazLikes[country][BaselineSeroStatus][PhaseSeverity] = CorrectPARAMS.VacHazLikes[country][BaselineSeroStatus][PhaseSeverity];
	}
	else	if (IsParamAKnot		(param_no, CorrectPARAMS.ParamNumbers)) //// for knots, must change yKnots as well as param vec. 
	{
		int knot	= Find_Knot_FromKnotParam	(param_no, HOUSE, CorrectPARAMS.ParamNumbers);
		int country = Find_Country_FromKnotParam(param_no, HOUSE, CorrectPARAMS.ParamNumbers, knot);

		//// reset knots
		ParamsToUpdate.yKnots[country][knot] = CorrectPARAMS.yKnots[country][knot];
		
		//// reset coeffs
		for (int poly = 0; poly < HOUSE.PolynomialsPerCountry; poly++)
			for (int coeff = 0; coeff <= HOUSE.MaxSplineDegree; coeff++)	ParamsToUpdate.SplineCoeffs[country][poly][coeff] = CorrectPARAMS.SplineCoeffs[country][poly][coeff];

		//// reset IVH_vals. affects both trial phases and serostatuses, but only this country. 
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
		{
			int Stratum_Index = HOUSE.Strata[VaccineGroup][BaselineSeroStatus];
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
				for (int i = 0; i < DATA.Set_Array[country][Stratum_Index].size(); i++)
				{
					int person_i = DATA.Set_Array[country][Stratum_Index][i];
					ParamsToUpdate.IVH_vals[person_i][PhaseSeverity] = CorrectPARAMS.IVH_vals[person_i][PhaseSeverity];
				}
		}

		//// reset VacHazLikes. affects both trial phases and serostatuses, but only this country. 
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
				ParamsToUpdate.VacHazLikes[country][BaselineSeroStatus][PhaseSeverity] = CorrectPARAMS.VacHazLikes[country][BaselineSeroStatus][PhaseSeverity];

		//// reset (integrated) hazards
		if (!(HOUSE.ModelVariant == SIMPLE_ANALYTICAL))	
		{
			for (int CalendarDay = 0; CalendarDay <= DATA.NumCalendarDaysFollowUp[country]; CalendarDay++)
			{
				ParamsToUpdate.BaselineHazardValues	[country][CalendarDay] = CorrectPARAMS.BaselineHazardValues[country][CalendarDay];
				ParamsToUpdate.IntBaseHazLookUp		[country][CalendarDay] = CorrectPARAMS.IntBaseHazLookUp	[country][CalendarDay];
			}
		}
		else 
		{
			ParamsToUpdate.yKnots[country][knot] = CorrectPARAMS.yKnots[country][knot];
			for (int poly = 0; poly < HOUSE.PolynomialsPerCountry; poly++)
			{
				ParamsToUpdate.AreasUnderFullPolynomials[country][poly] = CorrectPARAMS.AreasUnderFullPolynomials[country][poly];
				for (int coeff = 0; coeff <= HOUSE.MaxSplineDegree; coeff++)		ParamsToUpdate.SplineCoeffs[country][poly][coeff] = CorrectPARAMS.SplineCoeffs[country][poly][coeff];

				(*ParamsToUpdate.RootsArray[country][poly]) = (*CorrectPARAMS.RootsArray[country][poly]);
			}
		}
	}
	else	if (IsParamAHistHaz		(param_no, CorrectPARAMS.ParamNumbers)) //// NumC because of hist hazards
	{
		int country = Find_Country_FromHistHazParam(param_no, CorrectPARAMS.ParamNumbers);

		/// reset Active and Hospital KplusValues for all ages in this country. 
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int age = 0; age < HOUSE.HowManyAges; age++)
				ParamsToUpdate.KplusValues[country][age][PhaseSeverity] = CorrectPARAMS.KplusValues[country][age][PhaseSeverity];

		if (HOUSE.ModelVariant == K_SEROPOS && HOUSE.SeroSpecific_K_values && HOUSE.SeroSpecificEfficacies) //// Do the same for KplusValues with Effs. 
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
				for (int age = 0; age < HOUSE.HowManyAges; age++)
				ParamsToUpdate.Meta_KplusValues[With_Effs][country][age][PhaseSeverity] = CorrectPARAMS.Meta_KplusValues[With_Effs][country][age][PhaseSeverity];

		if (HOUSE.ModelVariant == AS_PRIME)
		{
			//// Ks_Prime
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
				for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
					for (int age = 0; age < HOUSE.HowManyAges; age++)
						ParamsToUpdate.Ks_Prime[country][PhaseSeverity][SeroPos][serotype][age] = CorrectPARAMS.Ks_Prime[country][PhaseSeverity][SeroPos][serotype][age];

			//// KPlusPrimeValues
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
				for (int age = 0; age < HOUSE.HowManyAges; age++)
					ParamsToUpdate.KPlusPrimeValues[country][age][PhaseSeverity] = CorrectPARAMS.KPlusPrimeValues[country][age][PhaseSeverity];
		}
		if (HOUSE.ModelVariant == K_SEROPOS || HOUSE.ModelVariant == AS_PRIME)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++) //// must reset VacHazLikes too. Affects both trial phases, but only seropositives in this country. 
				ParamsToUpdate.VacHazLikes[country][SeroPos][PhaseSeverity] = CorrectPARAMS.VacHazLikes[country][SeroPos][PhaseSeverity];

		if (!HOUSE.Empirical_SeroPrevs)
			for (int log_or_non_log_index = 0; log_or_non_log_index < 2; log_or_non_log_index++)
				for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
					for (int age = 0; age < HOUSE.HowManyAges; age++)
						ParamsToUpdate.SeroPrevs[log_or_non_log_index][country][BaselineSeroStatus][age] = CorrectPARAMS.SeroPrevs[log_or_non_log_index][country][BaselineSeroStatus][age]; 
	}
	else	if (IsParamA_qval		(param_no, CorrectPARAMS.ParamNumbers))
	{
		int qval	= Find_qval_From_qParam		(param_no, HOUSE, CorrectPARAMS.ParamNumbers);
		int country = Find_Country_From_qParam	(param_no, HOUSE, CorrectPARAMS.ParamNumbers, qval);

		//// not necessarilly all rhos/proportions affected by change in a given qvalue, but being conservative won't hurt. 
		for (int serotype = 0; serotype < HOUSE.N_STypes; serotype++)	
			ParamsToUpdate.rhos[country][serotype] = CorrectPARAMS.rhos[country][serotype];

		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
				ParamsToUpdate.SumRhoEffs[PhaseSeverity][country][BaselineSeroStatus] = CorrectPARAMS.SumRhoEffs[PhaseSeverity][country][BaselineSeroStatus];

		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int PrevInf = 0; PrevInf < HOUSE.Num_K_Params; PrevInf++)
				ParamsToUpdate.SumRhoKs[country][PhaseSeverity][PrevInf] = CorrectPARAMS.SumRhoKs[country][PhaseSeverity][PrevInf]; //// SumRhoKs affected by 

		if (HOUSE.EffNegWane == EffNegWane_Option::FROM_ZERO || HOUSE.EffNegWane == EffNegWane_Option::NO_WANE)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
				for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
					for (int country = 0; country < HOUSE.TotalCountries; country++)
						for (int age = 0; age < HOUSE.HowManyAges; age++)
							ParamsToUpdate.SumRhoEffNegs[PhaseSeverity][country][age][BaselineSeroStatus] = CorrectPARAMS.SumRhoEffNegs[PhaseSeverity][country][age][BaselineSeroStatus];

		if (HOUSE.SeroSpecific_K_values)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
				for (int country = 0; country < HOUSE.TotalCountries; country++)
				{
					for (int age = 0; age < HOUSE.HowManyAges; age++)
						ParamsToUpdate.KplusValues[country][age][PhaseSeverity] = CorrectPARAMS.KplusValues[country][age][PhaseSeverity];

					if (HOUSE.ModelVariant == K_SEROPOS) //// qval affects rhos affects KplusValues affects VacHazLikes if (HOUSE.ModelVariant == K_SEROPOS)
					{
						if (HOUSE.SeroSpecificEfficacies) //// Do the same for KplusValues with Effs. Don't need SS_Ks or K_SEROPOS conditions in if statement as HOUSE.SeroSpecific_K_values && HOUSE.ModelVariant == K_SEROPOS
							for (int age = 0; age < HOUSE.HowManyAges; age++)
								ParamsToUpdate.Meta_KplusValues[With_Effs][country][age][PhaseSeverity] = CorrectPARAMS.Meta_KplusValues[With_Effs][country][age][PhaseSeverity];

						for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
							ParamsToUpdate.VacHazLikes[country][BaselineSeroStatus][PhaseSeverity] = CorrectPARAMS.VacHazLikes[country][BaselineSeroStatus][PhaseSeverity];
					}

					if (HOUSE.SeroSpecificEfficacies && HOUSE.ASVE != Age_Option::INDEPENDENT && HOUSE.SSASVE_Additive)
						for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
							ParamsToUpdate.VacHazLikes[country][BaselineSeroStatus][PhaseSeverity] = CorrectPARAMS.VacHazLikes[country][BaselineSeroStatus][PhaseSeverity];
				}
	}
	else	if (IsParamA_rho		(param_no, CorrectPARAMS.ParamNumbers))
	{
		int serotype	= Find_rho_From_rhoParam	(param_no, CorrectPARAMS.ParamNumbers.Min_Rho, HOUSE.N_STypes);
		int country		= Find_Country_From_rhoParam(param_no, CorrectPARAMS.ParamNumbers.Min_Rho, HOUSE.N_STypes, serotype);

		ParamsToUpdate.rhos		[country][serotype]	= CorrectPARAMS.rhos	[country][serotype];
		ParamsToUpdate.SumRhos	[country]			= CorrectPARAMS.SumRhos	[country];

		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
				ParamsToUpdate.SumRhoEffs[PhaseSeverity][country][BaselineSeroStatus] = CorrectPARAMS.SumRhoEffs[PhaseSeverity][country][BaselineSeroStatus];

		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int PrevInf = 0; PrevInf < HOUSE.Num_K_Params; PrevInf++)
				ParamsToUpdate.SumRhoKs[country][PhaseSeverity][PrevInf] = CorrectPARAMS.SumRhoKs[country][PhaseSeverity][PrevInf]; 

		if (HOUSE.EffNegWane == EffNegWane_Option::FROM_ZERO || HOUSE.EffNegWane == EffNegWane_Option::NO_WANE)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
				for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
					for (int country = 0; country < HOUSE.TotalCountries; country++)
						for (int age = 0; age < HOUSE.HowManyAges; age++)
							ParamsToUpdate.SumRhoEffNegs[PhaseSeverity][country][age][BaselineSeroStatus] = CorrectPARAMS.SumRhoEffNegs[PhaseSeverity][country][age][BaselineSeroStatus];

		if (HOUSE.SeroSpecific_K_values)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
				for (int country = 0; country < HOUSE.TotalCountries; country++)
				{
					for (int age = 0; age < HOUSE.HowManyAges; age++)
						ParamsToUpdate.KplusValues[country][age][PhaseSeverity] = CorrectPARAMS.KplusValues[country][age][PhaseSeverity];

					if (HOUSE.ModelVariant == K_SEROPOS) //// qval affects rhos affects KplusValues affects VacHazLikes if (HOUSE.ModelVariant == K_SEROPOS)
					{
						if (HOUSE.SeroSpecificEfficacies) //// Do the same for KplusValues with Effs. Don't need SS_Ks or K_SEROPOS conditions in if statement as HOUSE.SeroSpecific_K_values && HOUSE.ModelVariant == K_SEROPOS
							for (int age = 0; age < HOUSE.HowManyAges; age++)
								ParamsToUpdate.Meta_KplusValues[With_Effs][country][age][PhaseSeverity] = CorrectPARAMS.Meta_KplusValues[With_Effs][country][age][PhaseSeverity];

						for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
							ParamsToUpdate.VacHazLikes[country][BaselineSeroStatus][PhaseSeverity] = CorrectPARAMS.VacHazLikes[country][BaselineSeroStatus][PhaseSeverity];
					}

					if (HOUSE.SeroSpecificEfficacies && HOUSE.ASVE != Age_Option::INDEPENDENT && HOUSE.SSASVE_Additive)
						for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
							ParamsToUpdate.VacHazLikes[country][BaselineSeroStatus][PhaseSeverity] = CorrectPARAMS.VacHazLikes[country][BaselineSeroStatus][PhaseSeverity];
				}
	}
	else	if (IsParamAnEfficacy	(param_no, CorrectPARAMS.ParamNumbers) || IsParamAn_atInf_Efficacy(param_no, CorrectPARAMS.ParamNumbers))
	{
		int BaselineSeroStatus	= Find_SeroStatus_FromEffParam		(param_no, HOUSE, CorrectPARAMS.ParamNumbers);
		int SeroType			= Find_Serotype_FromEffParam		(param_no, HOUSE, CorrectPARAMS.ParamNumbers, BaselineSeroStatus);
		int PhaseSeverity		= Find_PhaseSeverity_From_EffParam	(param_no, HOUSE, CorrectPARAMS.ParamNumbers, BaselineSeroStatus, SeroType);

		ParamsToUpdate.Efficacies[PhaseSeverity][SeroType][BaselineSeroStatus] = CorrectPARAMS.Efficacies[PhaseSeverity][SeroType][BaselineSeroStatus];

		if (HOUSE.ResidEffs)
			ParamsToUpdate.Inf_Effs[PhaseSeverity][SeroType][BaselineSeroStatus] = CorrectPARAMS.Inf_Effs[PhaseSeverity][SeroType][BaselineSeroStatus];

		for (int country = 0; country < HOUSE.TotalCountries; country++)
			ParamsToUpdate.SumRhoEffs[PhaseSeverity][country][BaselineSeroStatus] = CorrectPARAMS.SumRhoEffs[PhaseSeverity][country][BaselineSeroStatus];

		if (HOUSE.EffNegWane == EffNegWane_Option::FROM_ZERO || HOUSE.EffNegWane == EffNegWane_Option::NO_WANE)
			for (int country = 0; country < HOUSE.TotalCountries; country++)
				for (int age = 0; age < HOUSE.HowManyAges; age++)
					ParamsToUpdate.SumRhoEffNegs[PhaseSeverity][country][age][BaselineSeroStatus] = CorrectPARAMS.SumRhoEffNegs[PhaseSeverity][country][age][BaselineSeroStatus];

		ChooseMaxMinPhaseSeverities(HOUSE, PhaseSeverity);  //// changes HOUSE.MinMaxPhaseSeveritiesToLoopOver (potentially)

		bool Scenario_1 = HOUSE.SeroSpecificEfficacies && (HOUSE.ModelVariant == K_SEROPOS)		&& HOUSE.SeroSpecific_K_values && (BaselineSeroStatus == SeroPos);
		bool Scenario_2 = HOUSE.SeroSpecificEfficacies && HOUSE.ASVE != Age_Option::INDEPENDENT && HOUSE.SSASVE_Additive; 

		if (Scenario_1 || Scenario_2) //// in all other scenarios, efficacy doesn't affect VacHazLikes
			for (int PhaseSeverity = HOUSE.MinMaxPhaseSeveritiesToLoopOver[0]; PhaseSeverity < HOUSE.MinMaxPhaseSeveritiesToLoopOver[1]; PhaseSeverity++)
				for (int country = 0; country < HOUSE.TotalCountries; country++)
				{
					ParamsToUpdate.VacHazLikes[country][BaselineSeroStatus][PhaseSeverity] = CorrectPARAMS.VacHazLikes[country][BaselineSeroStatus][PhaseSeverity];

					if (Scenario_1)
						for (int age = 0; age < HOUSE.HowManyAges; age++)
							ParamsToUpdate.Meta_KplusValues[With_Effs][country][age][PhaseSeverity] = CorrectPARAMS.Meta_KplusValues[With_Effs][country][age][PhaseSeverity];
				}
	}
	else	if (IsParamAn_AS_Prime_Param(param_no, CorrectPARAMS.ParamNumbers))
	{
		int Type				= Find_ParamType_From_AS_Prime_Param	(param_no, HOUSE, CorrectPARAMS.ParamNumbers); //// either half life or coeff. 
		int BaselineSeroStatus	= Find_SeroStatus_From_AS_Prime_Param	(param_no, HOUSE, CorrectPARAMS.ParamNumbers, Type);

		ParamsToUpdate.AS_Priming_Params[BaselineSeroStatus][Type] = CorrectPARAMS.AS_Priming_Params[BaselineSeroStatus][Type];

		//// Ks_Prime
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++) //// as PrevInf > 0, looks like you'll only need SeroPos, but K1 and K2 affect K0prime,- and K1prime,+
				for (int Age = 0; Age < HOUSE.HowManyAges; Age++)
					for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++)
						ParamsToUpdate.Ks_Prime[HOUSE.WhichCountries[countryindex]][PhaseSeverity][BaselineSeroStatus][serotype][Age] = CorrectPARAMS.Ks_Prime[HOUSE.WhichCountries[countryindex]][PhaseSeverity][BaselineSeroStatus][serotype][Age];

		//// KPlusPrimeValues and SumRhoK0_SNeg_Primes: Both consider SS_Ks, therefore has rhos. Therefore need to do all countries for this K, even if Ks same between countries. Note you don't store arrays of KplusValues or KplusPrimes for single/case serotypes, only "survivors". 
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++) 
				for (int Age = 0; Age < HOUSE.HowManyAges; Age++)
				{
						 if (BaselineSeroStatus == SeroNeg)	ParamsToUpdate.SumRhoK0_SNeg_Primes	[HOUSE.WhichCountries[countryindex]][PhaseSeverity][Age] = CorrectPARAMS.SumRhoK0_SNeg_Primes[HOUSE.WhichCountries[countryindex]][PhaseSeverity][Age];
					else if (BaselineSeroStatus == SeroPos)	ParamsToUpdate.KPlusPrimeValues		[HOUSE.WhichCountries[countryindex]][Age][PhaseSeverity] = CorrectPARAMS.KPlusPrimeValues	[HOUSE.WhichCountries[countryindex]][Age][PhaseSeverity];
				}
		//// VacHazLikes, all countries and phase severities for that serostatus. 
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++)
				ParamsToUpdate.VacHazLikes[HOUSE.WhichCountries[countryindex]][BaselineSeroStatus][PhaseSeverity] = CorrectPARAMS.VacHazLikes[HOUSE.WhichCountries[countryindex]][BaselineSeroStatus][PhaseSeverity];
	}
	else	if (param_no == CorrectPARAMS.ParamNumbers.BS_BaseHazMult)
	{
		if (HOUSE.AdjHaz && !HOUSE.Skip_BS_BaseHazMult)
		{
			//// populate
			ParamsToUpdate.BS_BaseHazMults[HOUSE.Which_BS_BaseHazMult] = CorrectPARAMS.BS_BaseHazMults[HOUSE.Which_BS_BaseHazMult];

			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
				for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++)
					ParamsToUpdate.VacHazLikes[HOUSE.WhichCountries[countryindex]][HOUSE.Which_BS_BaseHazMult][PhaseSeverity] = 
					CorrectPARAMS. VacHazLikes[HOUSE.WhichCountries[countryindex]][HOUSE.Which_BS_BaseHazMult][PhaseSeverity]; // only necessary for the appropriate serostatus.
		}
	}

	
	else	std::cerr << "UpdateOrReset_ParamsEtc error: param_no not recognized" << endl;

#endif
}


//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** 
//// ****			2) Likelihood amendments
//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** 

void Change_K			(int PhaseSeverity, int PrevInf, int serotype, const DATA_struct &DATA, const Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE, const std::vector<int> &WhichCountries)
{
#ifdef PRINT_WITHIN_CHANGE_FUNCTIONS
	std::cout << "Change_K serotype " << serotype << " PhaseSeverity " << PhaseSeverity << " PrevInf " << PrevInf << " Value " << ProposedPARAMS.K_s[WhichCountries[0]][PhaseSeverity][PrevInf][serotype] << endl;
#endif
	//// want to have essentially the same code in the following four scenarios: i) neither SS_VEs nor SS_Ks (i.e. default); ii) SS_VEs only; iii) SS_Ks only; and iv) both SS_Ks & SS_VEs
	//// Want to avoid allocating memory, i.e. without allocating a container WhichVacHazSerotypesToUpdate 
	//// i)		neither SS_VEs nor SS_Ks:	serotype function argument always zero, HOUSE.N_STypes_VEs = 1;												WhichVacHazSerotypesToUpdate = 0; 
	//// ii)	SS_VEs only:				serotype function argument always zero, HOUSE.N_STypes_VEs = 4;	Must update IntVacHaz for all serotypes.	WhichVacHazSerotypesToUpdate = 0,1,2,3; 
	//// iii)	SS_Ks only:					serotype function argument varies, HOUSE.N_STypes_VEs = 1; 													WhichVacHazSerotypesToUpdate = 0; 
	//// iv)	SS_Ks & SS_VEs:				serotype function argument varies, HOUSE.N_STypes_VEs = 4, but only update one IntVacHaz at a time.;		WhichVacHazSerotypesToUpdate = serotype (i.e. serotype in function argument); 

	//// As per the above, definitions below state that if SSKs and SSVEs, want to consider one serotype at a time. Otherwise do all of them. 
	//// Note thought that you may redefine this if PrevInf > 0. 
	int Min_Vac_serotype = (HOUSE.SeroSpecificEfficacies && HOUSE.SeroSpecific_K_values) ? serotype			: 0; //// If both SS_Ks & SS_VEs, serotype function argument should refer to both efficacies and relative risks.  If doing SS_Ks and not SS_VEs, only one Eff value per baseline serostatus	, which affects L_Indices etc., hence why you need to be careful here. 
	int Max_Vac_Serotype = (HOUSE.SeroSpecificEfficacies && HOUSE.SeroSpecific_K_values) ? (serotype + 1)	: HOUSE.N_STypes_VEs;


	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
	///// Change (log) K cases (This loop true for all serotypes, PhaseSeverity's, PrevInfs and ModelVariants). 
	int Country;
	for (int countryindex = 0; countryindex < WhichCountries.size(); countryindex++)
	{
		Country = WhichCountries[countryindex];

		if (	HOUSE.ModelVariant == VAC_SILENT															|| 
				((HOUSE.ModelVariant == K_SEROPOS || HOUSE.ModelVariant == AS_PRIME)	&& ( PrevInf == 0))	||
				((HOUSE.ModelVariant == SIMPLE_NUMERICAL)								&& ( PrevInf < 2))	)
		ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.Ks[serotype][PhaseSeverity][PrevInf]]	=
		CurrentPARAMS. LikeParts[Country][HOUSE.L_Indices.Ks[serotype][PhaseSeverity][PrevInf]]	* log(ProposedPARAMS.K_s[Country][PhaseSeverity][PrevInf][serotype]) / log(CurrentPARAMS.K_s[Country][PhaseSeverity][PrevInf][serotype]);

		if (PrevInf > 0)
		{
			if (HOUSE.ModelVariant == VAC_SILENT || HOUSE.ModelVariant == K_SEROPOS || HOUSE.ModelVariant == AS_PRIME)
				ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.Kplus[serotype][PhaseSeverity]] = l_Kplus_case(DATA, ProposedPARAMS, HOUSE, Country, PhaseSeverity, serotype);

			if (HOUSE.ModelVariant == AS_PRIME)
			{
				ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.KPrimes	[serotype][PhaseSeverity]] = l_K0SNegPrime_case(DATA, ProposedPARAMS, HOUSE, Country, PhaseSeverity, serotype);
				ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.KplusPrime[serotype][PhaseSeverity]] = l_KPlusPrime_case(DATA, ProposedPARAMS, HOUSE, Country, PhaseSeverity, serotype);
			}
		}
	}

	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
	///// Change integrated hazards (baseline and vaccine). 

	//// K_0 is hazard multiplier for control group only for VAC_SILENT, otherwise for both trial arms. Remember, loops will be < MaxTrialArm_IBH, so for VAC_SILENT and AS_PRIME, will be < VaccineGroup, i.e. only the control group. 
	//// Confusingly, this is different to HOUSE.NumTrialArms_IBH, because that refers to number of trial arms where an integrated baseline hazard calculation is required
	//// MaxTrialArm_IBH as defined here refers to number of trial arms affected by a change in K, AND THAT REQUIRE THE SAME KIND OF CHANGE (either full recalculation or simple multiplier change). 
	//// DO NOT change this definition. 
	int MaxTrialArm_IBH = (HOUSE.ModelVariant == VAC_SILENT || HOUSE.ModelVariant == AS_PRIME) ? VaccineGroup : HOUSE.NumTrialArms_IBH;

	if (PrevInf == 0)
	{
		int BaselineSeroStatus = SeroNeg; 

		for (int countryindex = 0; countryindex < WhichCountries.size(); countryindex++)
		{
			Country = WhichCountries[countryindex]; 

			///// Change IntBase Haz
			for (int TrialArm = 0; TrialArm < MaxTrialArm_IBH; TrialArm++)
				ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][SeroNeg][TrialArm]]	=
				CurrentPARAMS. LikeParts[Country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][SeroNeg][TrialArm]]	* ProposedPARAMS.SumRhoKs[Country][PhaseSeverity][PrevInf]		/ CurrentPARAMS.SumRhoKs[Country][PhaseSeverity][PrevInf];
		
			///// Change IntVac Haz (not for VAC_SILENT or AS_PRIME). 
			if ((HOUSE.ModelVariant == SIMPLE_NUMERICAL) || (HOUSE.ModelVariant == K_SEROPOS))
				for (int Vac_serotype = Min_Vac_serotype; Vac_serotype < Max_Vac_Serotype; Vac_serotype++) //// 
					if (HOUSE.SeroSpecificEfficacies && HOUSE.SeroSpecific_K_values)	//// must use serotype specific K 
						ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][Vac_serotype]] =
						CurrentPARAMS. LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][Vac_serotype]] * ProposedPARAMS.K_s[Country][PhaseSeverity][PrevInf][Vac_serotype] / CurrentPARAMS.K_s[Country][PhaseSeverity][PrevInf][Vac_serotype];
					else
						ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][Vac_serotype]] =
						CurrentPARAMS. LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][Vac_serotype]] * ProposedPARAMS.SumRhoKs[Country][PhaseSeverity][PrevInf]		/ CurrentPARAMS.SumRhoKs[Country][PhaseSeverity][PrevInf];
		}
	}
	else
	{
		//// in control arm: K1 maps to seropositive
		//// in vaccine arm: relative risk of K1 (baseline if ActiveMild) in vaccine arm will be seropositive if SIMPLE_NUMERICAL or K_SEROPOS, otherwise seronegative. 
		//// In vaccine arm, for AS_PRIME, both K1 and K2 affect both Seropositive and seronegative, but for control arm only seropositive. 
		int BaselineSeroStatus_VaccineArm = (HOUSE.ModelVariant == VAC_SILENT && PrevInf == 1) ? SeroNeg : SeroPos;

		//// REDEFINE Min_Vac_serotype and Max_Vac_Serotype if necessary. Essentially, if you don't divide up likecomponents by serotype, then a change in one seroK affects a change in the summed vaccine hazard. Not true if you do divide up. 
		if ((HOUSE.ModelVariant == K_SEROPOS || HOUSE.AdditiveSSASVEs) && 
			(HOUSE.SeroSpecificEfficacies && HOUSE.SeroSpecific_K_values))	
		{	
			Min_Vac_serotype = 0; 
			Max_Vac_Serotype = HOUSE.N_STypes_VEs;
		}  ///// ///// For K_SEROPOS seropositive SSKs and SSVEs, because there are no common factors to take out, then you don't have 4 separate likelihood components for the integrated vaccine hazards. Hence a change in the K for any one serotype affects the amalgamated IntVacHaz component. So do all to be safe. 

		for (int countryindex = 0; countryindex < WhichCountries.size(); countryindex++)
		{
			Country = WhichCountries[countryindex];

			//// IntBaseHaz
			for (int TrialArm = 0; TrialArm < MaxTrialArm_IBH; TrialArm++)
				if ((HOUSE.ModelVariant == VAC_SILENT) || (HOUSE.ModelVariant == K_SEROPOS || HOUSE.ModelVariant == AS_PRIME))			/// l_Int_Base_Haz must be fully recalculated
					ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][SeroPos][TrialArm]]	= l_Int_Base_Haz(Country, DATA, ProposedPARAMS, HOUSE, PhaseSeverity, SeroPos, TrialArm);
				else																					/// l_Int_Base_Haz can be quickly recalculated
					ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][SeroPos][TrialArm]]	=
					CurrentPARAMS. LikeParts[Country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][SeroPos][TrialArm]]	* ProposedPARAMS.SumRhoKs[Country][PhaseSeverity][PrevInf] / CurrentPARAMS.SumRhoKs[Country][PhaseSeverity][PrevInf]; //// can recalculate l_Int_Base_Haz_Neg_Vac quickly

			//// IntVacHaz
			for (int Vac_serotype = Min_Vac_serotype; Vac_serotype < Max_Vac_Serotype; Vac_serotype++)
			{
				if (HOUSE.ModelVariant == K_SEROPOS || (HOUSE.AdditiveSSASVEs && HOUSE.SeroSpecific_K_values))										/// l_Int_Vac_Haz must be recalculated. VacHazLikes done in AmendProposedParams function (because of K+ multipliers for Seropositives). This function changes just the efficacies (and rhos)
						ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus_VaccineArm][Vac_serotype]] = l_Int_Vac_Haz_sero_i(ProposedPARAMS, HOUSE, Country, Vac_serotype, BaselineSeroStatus_VaccineArm, PhaseSeverity);
				else if ((HOUSE.ModelVariant == VAC_SILENT) || (HOUSE.ModelVariant == SIMPLE_NUMERICAL))		/// l_Int_Vac_Haz can be quickly recalculated
					if (HOUSE.SeroSpecificEfficacies && HOUSE.SeroSpecific_K_values)	//// must use serotype specific K 
						ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus_VaccineArm][Vac_serotype]]	=
						CurrentPARAMS. LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus_VaccineArm][Vac_serotype]]	* ProposedPARAMS.K_s[Country][PhaseSeverity][PrevInf][Vac_serotype] / CurrentPARAMS.K_s[Country][PhaseSeverity][PrevInf][Vac_serotype];
					else																//// all IntVacHaz components have same "survive" K, whether that is with rho's or not. 
						ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus_VaccineArm][Vac_serotype]]	=
						CurrentPARAMS. LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus_VaccineArm][Vac_serotype]]	* ProposedPARAMS.SumRhoKs[Country][PhaseSeverity][PrevInf]			/ CurrentPARAMS.SumRhoKs[Country][PhaseSeverity][PrevInf];

				else if (HOUSE.ModelVariant == AS_PRIME)																				/// l_Int_Vac_Haz must be recalculated. VacHazLikes done in AmendProposedParams function (because of K+prime multipliers for Seropositives, and K0Prime multipliers for seronegatives). This function changes just the efficacies (and rhos)
					for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
						ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][Vac_serotype]] = l_Int_Vac_Haz_sero_i(ProposedPARAMS, HOUSE, Country, Vac_serotype, BaselineSeroStatus, PhaseSeverity);
			}

			//// IntBaseHaz for vaccine arm
			if (HOUSE.ModelVariant == AS_PRIME)
				for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
					ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus][VaccineGroup]]			= l_Int_Base_Haz(Country, DATA, ProposedPARAMS, HOUSE, PhaseSeverity, BaselineSeroStatus			, VaccineGroup);
			else if (HOUSE.ModelVariant == K_SEROPOS || (HOUSE.AdditiveSSASVEs && HOUSE.SeroSpecific_K_values) || ((HOUSE.EffNegWane == EffNegWane_Option::NO_WANE && HOUSE.SeroSpecific_K_values)))
				ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus_VaccineArm][VaccineGroup]]	= l_Int_Base_Haz(Country, DATA, ProposedPARAMS, HOUSE, PhaseSeverity, BaselineSeroStatus_VaccineArm	, VaccineGroup); 
			else if (HOUSE.ModelVariant == VAC_SILENT)
				ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus_VaccineArm][VaccineGroup]]	=
				CurrentPARAMS. LikeParts[Country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus_VaccineArm][VaccineGroup]]	* ProposedPARAMS.SumRhoKs[Country][PhaseSeverity][PrevInf] / CurrentPARAMS.SumRhoKs[Country][PhaseSeverity][PrevInf]; //// can recalculate l_Int_Base_Haz_Neg_Vac quickly
		}
	}

#ifdef PRINT_WITHIN_CHANGE_FUNCTIONS
	std::cout << "Finished Change_K " <<  endl;
#endif

}
void Change_K			(int PhaseSeverity, int PrevInf, int serotype, const DATA_struct &DATA, const Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE)
{
	Change_K(PhaseSeverity, PrevInf, serotype, DATA, CurrentPARAMS, ProposedPARAMS, HOUSE, HOUSE.WhichCountries);
}
void Change_Efficacy	(const DATA_struct &DATA, const Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE, char BaselineSeroStatus, int &serotype, int PhaseSeverity)
{
#ifdef PRINT_WITHIN_CHANGE_FUNCTIONS
	std::cout << "Change_Efficacy PhaseSeverity " << PhaseSeverity << " BaselineSeroStatus " << int(BaselineSeroStatus) << " serotype " << serotype << endl;
#endif

	int Country; 

	// choose likelihood vector indices. 
	DType Current_Efficacy, Proposed_Efficacy;

	bool Scenario_1 = HOUSE.SeroSpecificEfficacies && (HOUSE.ModelVariant == K_SEROPOS)		&& HOUSE.SeroSpecific_K_values && (BaselineSeroStatus == SeroPos);
	bool Scenario_2 = HOUSE.SeroSpecificEfficacies && (HOUSE.ASVE != Age_Option::INDEPENDENT) && HOUSE.SSASVE_Additive; 

	for (int countryindex = 0; countryindex < HOUSE.NoCountriesToFit; countryindex++)
    {
		Country = HOUSE.WhichCountries[countryindex];

		for (int PhaseSeverity = HOUSE.MinMaxPhaseSeveritiesToLoopOver[0]; PhaseSeverity < HOUSE.MinMaxPhaseSeveritiesToLoopOver[1]; PhaseSeverity++) 
		{
			Current_Efficacy	= CurrentPARAMS. Efficacies[PhaseSeverity][serotype][BaselineSeroStatus];
			Proposed_Efficacy	= ProposedPARAMS.Efficacies[PhaseSeverity][serotype][BaselineSeroStatus];

			if (HOUSE.EffNegWane == EffNegWane_Option::FROM_ZERO) 
			{
				Current_Efficacy	= abs(Current_Efficacy	); 
				Proposed_Efficacy	= abs(Proposed_Efficacy	); 
			}

			if (Current_Efficacy == 0) ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] = l_Int_Vac_Haz_sero_i_long_way(BaselineSeroStatus, Country, DATA, ProposedPARAMS, HOUSE, PhaseSeverity, serotype);
			else
			{
#ifdef VACHAZLIKE_SLOW_WAY
				ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] = l_Int_Vac_Haz_sero_i_long_way(BaselineSeroStatus, Country, DATA, ProposedPARAMS, HOUSE, PhaseSeverity, serotype);
#else			
				if (Scenario_1 || Scenario_2)  ///// For K_SEROPOS seropositive SSKs and SSVEs, because there are no common factors to take out, then you don't have 4 separate likelihood components for the integrated vaccine hazards. Hence a change in the K for any one serotype affects the amalgamated IntVacHaz component. So do all to be safe. 
					for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
						ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] = l_Int_Vac_Haz_sero_i(ProposedPARAMS, HOUSE, Country, serotype, BaselineSeroStatus, PhaseSeverity);
				else
					ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] =
					CurrentPARAMS. LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] * Proposed_Efficacy / Current_Efficacy;
#endif		
			}

			ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.WaningEffs[serotype][BaselineSeroStatus][PhaseSeverity]] = l_WaningEfficacy(BaselineSeroStatus, Country, DATA, ProposedPARAMS, HOUSE, serotype, PhaseSeverity); /// not require in TWO_CASE_CATEGORIES if statement as can only be a case in at most one trial phase. 
			
			if (HOUSE.ResidEffs) //// i.e. recalculate log(1 - VE_inf) part of likelihood. Update: not needed as covered in l_WaningEfficacy above
			{
				DType Current_InfEff	= CurrentPARAMS. Inf_Effs[PhaseSeverity][serotype][BaselineSeroStatus]; 
				DType Proposed_InfEff	= ProposedPARAMS.Inf_Effs[PhaseSeverity][serotype][BaselineSeroStatus]; 

				ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus][VaccineGroup]] =
				CurrentPARAMS. LikeParts[Country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus][VaccineGroup]] * (1 - Proposed_InfEff) / (1 - Current_InfEff);
			}

			if (HOUSE.EffNegWane == EffNegWane_Option::FROM_ZERO || HOUSE.EffNegWane == EffNegWane_Option::NO_WANE)
				if (CanEfficacyBeNegative(BaselineSeroStatus, HOUSE)) 
					for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++)
						for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
							ProposedPARAMS.LikeParts[HOUSE.WhichCountries[countryindex]][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus][VaccineGroup]] = l_Int_Base_Haz(HOUSE.WhichCountries[countryindex], DATA, ProposedPARAMS, HOUSE, PhaseSeverity, BaselineSeroStatus, VaccineGroup);
		}
    }
}
void Change_Waning		(const DATA_struct &DATA, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE, char BaselineSeroStatus)
{
#ifdef PRINT_WITHIN_CHANGE_FUNCTIONS
	std::cout << "Change_Waning BaselineSeroStatus " << int(BaselineSeroStatus) << endl;
#endif
	int Country; 
	for (int countryindex = 0; countryindex < HOUSE.NoCountriesToFit; countryindex++)
	{
		Country = HOUSE.WhichCountries[countryindex];

		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++) 
			for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
			{
#ifdef VACHAZLIKE_SLOW_WAY
				ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]]	= l_Int_Vac_Haz_sero_i_long_way(BaselineSeroStatus, Country, DATA, ProposedPARAMS, HOUSE, PhaseSeverity, serotype);
#else
				ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]]	= l_Int_Vac_Haz_sero_i(ProposedPARAMS, HOUSE, Country, serotype, BaselineSeroStatus, PhaseSeverity);
#endif
				ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.WaningEffs[serotype][BaselineSeroStatus][PhaseSeverity]]	= l_WaningEfficacy(BaselineSeroStatus, Country, DATA, ProposedPARAMS, HOUSE, serotype, PhaseSeverity);
			}
	}
}
void Change_ASVE_Param	(const DATA_struct &DATA, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE, char BaselineSeroStatus)
{
	//// This function is essentially a wrapper for change vaccine hazard likelihood components. 
	//// However if you're doing HOUSE.AltWaneEffNeg, then you also need to change int base hazard for vaccine arm (due to waning function now being given by -exp(-t/tau) + 1 if "efficacy" negative. 
	Change_Waning(DATA, ProposedPARAMS, HOUSE, BaselineSeroStatus); 

	if (HOUSE.EffNegWane == EffNegWane_Option::FROM_ZERO || HOUSE.EffNegWane == EffNegWane_Option::NO_WANE)
		if (CanEfficacyBeNegative(BaselineSeroStatus, HOUSE)) 
			for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++)
				for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++) 
					ProposedPARAMS.LikeParts[HOUSE.WhichCountries[countryindex]][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus][VaccineGroup]] = l_Int_Base_Haz(HOUSE.WhichCountries[countryindex], DATA, ProposedPARAMS, HOUSE, PhaseSeverity, BaselineSeroStatus, VaccineGroup);
}

void Change_AS_Haz		(const DATA_struct &DATA, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE)
{
#ifdef PRINT_WITHIN_CHANGE_FUNCTIONS
	std::cout << "Change_AS_Haz " << endl;
#endif
	int Country;
	for (int countryindex = 0; countryindex < HOUSE.NoCountriesToFit; countryindex++)
	{
		Country = HOUSE.WhichCountries[countryindex];

		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++) 
		{
			//// AS_Haz_Cases
			ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.AgeHazMult[PhaseSeverity]] = l_AS_HazMult(Country, DATA, ProposedPARAMS, HOUSE, PhaseSeverity);

			//// IntBaseHaz
			for (int TrialArm = 0; TrialArm < HOUSE.HowManyTrialArms; TrialArm++)
				for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
					ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus][TrialArm]] = l_Int_Base_Haz(Country, DATA, ProposedPARAMS, HOUSE, PhaseSeverity, BaselineSeroStatus, TrialArm);

			//// IntVacHaz
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
				for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
#ifdef VACHAZLIKE_SLOW_WAY
					ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]]	= l_Int_Vac_Haz_sero_i_long_way(BaselineSeroStatus, Country, DATA, ProposedPARAMS, HOUSE, PhaseSeverity, serotype);
#else
					ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]]	= l_Int_Vac_Haz_sero_i(ProposedPARAMS, HOUSE, Country, serotype, BaselineSeroStatus, PhaseSeverity);
#endif
		}
	}
}
void Change_Knot		(int &country, int WhichKnot, int param_no, const DATA_struct &DATA, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE)
{
#ifdef PRINT_WITHIN_CHANGE_KNOT_FUNCTIONS
	std::cout << "Change_Knot country " << country << " knot " << WhichKnot << " param_no " << param_no << endl;
#endif
	if (!(ProposedPARAMS.yKnots[country][WhichKnot] == ProposedPARAMS.ParamVec[param_no])) std::cerr << "ChangeKnot error: ProposedPARAMS.yKnots[country][WhichKnot] NOT EQUAL TO ProposedPARAMS.ParamVec[param_no]" << endl;
	if (ProposedPARAMS.yKnots[country][WhichKnot] < 0) std::cerr << "ChangeKnot error: NewYKnotValue < 0" << endl;
	
	//// Change case Baseline Hazard
	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
		ProposedPARAMS.LikeParts[country][HOUSE.L_Indices.BaseHaz[PhaseSeverity]] = l_BaseHaz		(country, DATA, ProposedPARAMS, HOUSE, PhaseSeverity);								///// cases				haz 1,...4. 

	//// Change Integrated Baseline Hazard values for each PhaseSeverity, BaselineSerostatus and TrialArm. 
	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			for (int TrialArm = 0; TrialArm < HOUSE.NumTrialArms_IBH; TrialArm++)
				ProposedPARAMS.LikeParts[country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus][TrialArm]] = l_Int_Base_Haz(country, DATA, ProposedPARAMS, HOUSE, PhaseSeverity, BaselineSeroStatus, TrialArm);
	
	//// Change Integrated Vaccine Hazard values for each PhaseSeverity and BaselineSerostatus 
	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++) //// this takes account of hospital severe as well if you're doing that. Keep loop separate from below - you want serotype to be last in nested loop. For loop below don't want to needlessly repeat calculations. 
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
#ifdef VACHAZLIKE_SLOW_WAY
				ProposedPARAMS.LikeParts[country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] = l_Int_Vac_Haz_sero_i_long_way(Set_Array, BaselineSeroStatus, country, DATA, ProposedPARAMS, HOUSE, PhaseSeverity, serotype);
#else
				ProposedPARAMS.LikeParts[country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] = l_Int_Vac_Haz_sero_i(ProposedPARAMS, HOUSE, country, serotype, BaselineSeroStatus, PhaseSeverity);	///// survivors
#endif
}
void Change_HistHazard	(int &country, const DATA_struct &DATA, const Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE)  
{
#ifdef PRINT_WITHIN_CHANGE_HH_FUNCTIONS
	std::cout << "Change_HistHazard country " << country << endl;
#endif
	DType Old_hci = CurrentPARAMS. ParamVec[CurrentPARAMS. ParamNumbers.Min_HHaz + country];
	DType New_hci = ProposedPARAMS.ParamVec[ProposedPARAMS.ParamNumbers.Min_HHaz + country];

	if (!HOUSE.Empirical_SeroPrevs)
	{
		ProposedPARAMS.LikeParts[country][HOUSE.L_Indices.BSs[SeroNeg]] = New_hci * CurrentPARAMS.LikeParts[country][HOUSE.L_Indices.BSs[SeroNeg]] / Old_hci;
		ProposedPARAMS.LikeParts[country][HOUSE.L_Indices.BSs[SeroPos]] = l_Aug(country, SeroPos, DATA, ProposedPARAMS, HOUSE);		// note you are not looping through countries, so do not need to specify country part of LikelihoodVec index with House.WhichCountries[country] as you do for global parameters. 
	}

	////// then K plus values. 
	if (HOUSE.ModelVariant == VAC_SILENT || HOUSE.ModelVariant == K_SEROPOS || HOUSE.ModelVariant == AS_PRIME)
		for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
				ProposedPARAMS.LikeParts[country][HOUSE.L_Indices.Kplus[serotype][PhaseSeverity]] = l_Kplus_case(DATA, ProposedPARAMS, HOUSE, country, PhaseSeverity, serotype);

	//// HH's affect K+ (and K+prime) terms. 
	//// For VAC_SILENT, HHs affects SeroPos Control Group. 
	//// For K_SEROPOS, HH affects SeroPos either arm 
	//// For AS_PRIME, HH affects SeroPos either arm (K+ Control Group, K+prime Vaccine Group)
	int MaxTrialArm_IBH = (HOUSE.ModelVariant == VAC_SILENT) ? VaccineGroup : HOUSE.NumTrialArms_IBH; //// i.e. for VAC_SILENT HH doesn't affect Seropositive vaccinees, only seropositive controls. Think also true of AS_PRIME but need to check. 
	//// Change Integrated Baseline Hazard values for each PhaseSeverity, BaselineSerostatus and TrialArm. 
	if (HOUSE.ModelVariant == VAC_SILENT || HOUSE.ModelVariant == K_SEROPOS || HOUSE.ModelVariant == AS_PRIME)
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int TrialArm = 0; TrialArm < MaxTrialArm_IBH; TrialArm++)
				ProposedPARAMS.LikeParts[country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][SeroPos][TrialArm]] = l_Int_Base_Haz(country, DATA, ProposedPARAMS, HOUSE, PhaseSeverity, SeroPos, TrialArm);

	if (HOUSE.ModelVariant == K_SEROPOS || HOUSE.ModelVariant == AS_PRIME) 	//// Change Vaccine hazards if necessary.  
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
				ProposedPARAMS.LikeParts[country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][SeroPos][serotype]] = l_Int_Vac_Haz_sero_i(ProposedPARAMS, HOUSE, country, serotype, SeroPos, PhaseSeverity);

	if (HOUSE.ModelVariant == AS_PRIME)
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
				ProposedPARAMS.LikeParts[country][HOUSE.L_Indices.KplusPrime[serotype][PhaseSeverity]] = l_KPlusPrime_case(DATA, ProposedPARAMS, HOUSE, country, PhaseSeverity, serotype);

}
void Change_Rho			(int &country, int &serotype, const DATA_struct &DATA, const Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE)
{
#ifdef PRINT_WITHIN_CHANGE_FUNCTIONS
	std::cerr << "Change_Rho country " << country << " serotype " << serotype << endl; 
#endif

	//// For cases - need to change log(rho_c_d)
	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
		if (!HOUSE.BaselinePartition)
				ProposedPARAMS.LikeParts[country][HOUSE.L_Indices.rhos[serotype][PhaseSeverity]] =
				CurrentPARAMS. LikeParts[country][HOUSE.L_Indices.rhos[serotype][PhaseSeverity]] * log(ProposedPARAMS.rhos[country][serotype]) / log(CurrentPARAMS.rhos[country][serotype]);
		else
		{
			if (CurrentPARAMS.rhos[country][serotype] == 1)		ProposedPARAMS.LikeParts[country][HOUSE.L_Indices.rhos[serotype][PhaseSeverity]] = l_rho_case(DATA, ProposedPARAMS, HOUSE, country, serotype, PhaseSeverity);
			else
			{
				ProposedPARAMS.LikeParts[country][HOUSE.L_Indices.rhos[serotype][PhaseSeverity]] =
				CurrentPARAMS. LikeParts[country][HOUSE.L_Indices.rhos[serotype][PhaseSeverity]] * log(ProposedPARAMS.rhos[country][serotype]) / log(CurrentPARAMS.rhos[country][serotype]);
			} 	 
		}

	//// change IntVacHaz (and possibly IntBaseHaz) quantities. 
	if (HOUSE.SeroSpecificEfficacies) 
	{
		//// change IntVacHaz quantities. 
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			{
				bool Scenario_1 = /*HOUSE.SeroSpecificEfficacies && - already counted */ ((HOUSE.ModelVariant == K_SEROPOS) || (HOUSE.ModelVariant == AS_PRIME)) && HOUSE.SeroSpecific_K_values && (BaselineSeroStatus == SeroPos);
				bool Scenario_2 = HOUSE.AdditiveSSASVEs; //// true for SSKs or not. 

				//// IntVacHaz
#ifdef VACHAZLIKE_SLOW_WAY
				ProposedPARAMS.LikeParts[CLkInd + HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] = l_Int_Vac_Haz_sero_i_long_way(Set_Array, BaselineSeroStatus, country, DATA, ProposedPARAMS, HOUSE, PhaseSeverity, serotype);
#else			
				if (Scenario_1 || Scenario_2)
					for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
						ProposedPARAMS.LikeParts[country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] = l_Int_Vac_Haz_sero_i(ProposedPARAMS, HOUSE, country, serotype, BaselineSeroStatus, PhaseSeverity);
				else
						ProposedPARAMS.LikeParts[country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] =
						CurrentPARAMS. LikeParts[country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] * ProposedPARAMS.rhos[country][serotype] / CurrentPARAMS.rhos[country][serotype];
#endif	
				if (HOUSE.EffNegWane == EffNegWane_Option::FROM_ZERO || HOUSE.EffNegWane == EffNegWane_Option::NO_WANE)
					if (CanEfficacyBeNegative(BaselineSeroStatus, HOUSE)) 
						ProposedPARAMS.LikeParts[country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus][VaccineGroup]] = l_Int_Base_Haz(country, DATA, ProposedPARAMS, HOUSE, PhaseSeverity, BaselineSeroStatus, VaccineGroup);
			}

		if (HOUSE.BaselinePartition) 
			if (!HOUSE.SeroSpecific_K_values)
			{
				DType SumRho_Ratio = ProposedPARAMS.SumRhos[country] / CurrentPARAMS.SumRhos[country];
				for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
					for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
						for (int TrialArm = 0; TrialArm < HOUSE.NumTrialArms_IBH; TrialArm++)
							ProposedPARAMS.LikeParts[country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus][TrialArm]] = 
							CurrentPARAMS. LikeParts[country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus][TrialArm]] * SumRho_Ratio;
			}
	}

	if (HOUSE.SeroSpecific_K_values)
	{
		///// Change Integrated baseline hazards. 
		///// first change SumRhoK's and KplusValues in ProposedPARAMS (done in Amend function). Then change common factors when you can and recaculate when necessary (depending on ModelVariant). 
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
				for (int TrialArm = 0; TrialArm < HOUSE.NumTrialArms_IBH; TrialArm++)
				{
					if (	(HOUSE.ModelVariant == VAC_SILENT	&& BaselineSeroStatus == SeroPos && TrialArm == ControlGroup)		|| /// For SS_Ks, change in rho changes KPlusValues, so have to recalculate l_Int_Base_Haz fully if doing Seropositive Controls (not for SIMPLE_NUMERICAL though). 
							(HOUSE.ModelVariant == K_SEROPOS	&& BaselineSeroStatus == SeroPos)									|| 
							(HOUSE.ModelVariant == AS_PRIME		&& (BaselineSeroStatus == SeroPos || TrialArm == VaccineGroup))		||
							((HOUSE.EffNegWane == EffNegWane_Option::FROM_ZERO || HOUSE.EffNegWane == EffNegWane_Option::NO_WANE) && TrialArm == VaccineGroup)
						) // For AS_PRIME, need to recalculate IntBaseHaz fully for each serostatus in VaccineGroup, but only for Seropositive controls. 
						ProposedPARAMS.LikeParts[country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus][TrialArm]] = l_Int_Base_Haz(country, DATA, ProposedPARAMS, HOUSE, PhaseSeverity, BaselineSeroStatus, TrialArm);
					else
					{
						if (HOUSE.ModelVariant == VAC_SILENT || HOUSE.ModelVariant == AS_PRIME)
							ProposedPARAMS.LikeParts[country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus][TrialArm]] = 
							CurrentPARAMS. LikeParts[country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus][TrialArm]] * ProposedPARAMS.SumRhoKs[country][PhaseSeverity][BaselineSeroStatus + TrialArm] / CurrentPARAMS.SumRhoKs[country][PhaseSeverity][BaselineSeroStatus + TrialArm]; //// BaselineSeroStatus + TrialArm as Vaccine is Silent Infection. 
						else if (HOUSE.ModelVariant == K_SEROPOS)
							ProposedPARAMS.LikeParts[country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus][TrialArm]] = 
							CurrentPARAMS. LikeParts[country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus][TrialArm]] * ProposedPARAMS.SumRhoKs[country][PhaseSeverity][0] / CurrentPARAMS.SumRhoKs[country][PhaseSeverity][0]; //// BaselineSeroStatus + TrialArm as Vaccine is Silent Infection. 
						else std::cerr << "Change_Rho ERROR: ModelVariant not recognized" << endl; 
					}
				}

		///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
		///// Change Integrated Vaccine hazards. 
		int K_mult_Index; 

		if (!HOUSE.SeroSpecificEfficacies) //// your SeroSpecificEfficacies change above
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			{
				K_mult_Index = BaselineSeroStatus; 
				if (HOUSE.ModelVariant == VAC_SILENT) K_mult_Index++;  //// this pretty much defines the VAC_SILENT model variant. 

				for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
					for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
					{
						if (HOUSE.ModelVariant == AS_PRIME || (HOUSE.ModelVariant == K_SEROPOS && BaselineSeroStatus == SeroPos)) //// must recalculate IntVac_Haz parts becuase of K_seropos multiplier. 
							ProposedPARAMS.LikeParts[country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] = l_Int_Vac_Haz_sero_i(ProposedPARAMS, HOUSE, country, serotype, BaselineSeroStatus, PhaseSeverity);
						else
							ProposedPARAMS.LikeParts[country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] =
							CurrentPARAMS. LikeParts[country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] * ProposedPARAMS.SumRhoKs[country][PhaseSeverity][K_mult_Index] / CurrentPARAMS.SumRhoKs[country][PhaseSeverity][K_mult_Index];
					}
			}
	}
}
void Change_qVal		(int &qval, int &country, const DATA_struct &DATA, const Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE)
{
#ifdef PRINT_WITHIN_CHANGE_FUNCTIONS
	std::cerr << "Change_qVal country " << country << " qval " << qval << endl;
#endif
	/*
		change in qc0 affects		ALL serotypes		in country c. (i.e. 0,1,2,3 in cpp indexing)
		change in qc1 affects		serotypes 2,3,4		in country c. (i.e.   1,2,3 in cpp indexing) 
		change in qc2 affects		serotypes 3,4		in country c. (i.e.     2,3 in cpp indexing)
	*/
	for (int serotype = qval; serotype < HOUSE.N_STypes; serotype++)	
		Change_Rho(country, serotype, DATA, CurrentPARAMS, ProposedPARAMS, HOUSE);
}
void Change_AS_Prime	(int BaselineSeroStatus, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE, const DATA_struct &DATA)
{
	if (HOUSE.ModelVariant != AS_PRIME) std::cerr << "Change_AS_Prime ERROR: HOUSE.ModelVariant != AS_PRIME" << endl << "Change_AS_Prime ERROR: HOUSE.ModelVariant != AS_PRIME" << endl << "Change_AS_Prime ERROR: HOUSE.ModelVariant != AS_PRIME" << endl;
#ifdef PRINT_WITHIN_CHANGE_FUNCTIONS
	std::cout << "Change_AS_Prime BaselineSeroStatus " << int(BaselineSeroStatus) << endl;
#endif
	int Country; 
	for (int countryindex = 0; countryindex < HOUSE.NoCountriesToFit; countryindex++)
	{
		Country = HOUSE.WhichCountries[countryindex];

		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
		{
			//// KPrimes Cases
			if (BaselineSeroStatus == SeroNeg)
				for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
					ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.KPrimes[serotype][PhaseSeverity]] = l_K0SNegPrime_case(DATA, ProposedPARAMS, HOUSE, Country, PhaseSeverity, serotype);
			//// KplusPrime Cases
			else if (BaselineSeroStatus == SeroPos)
				for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
					ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.KplusPrime[serotype][PhaseSeverity]] = l_KPlusPrime_case(DATA, ProposedPARAMS, HOUSE, Country, PhaseSeverity, serotype);

			//// IntBaseHaz
			ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus][VaccineGroup]] = l_Int_Base_Haz(Country, DATA, ProposedPARAMS, HOUSE, PhaseSeverity, BaselineSeroStatus, VaccineGroup);

			//// IntVacHaz
			for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
#ifdef VACHAZLIKE_SLOW_WAY
				ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] = l_Int_Vac_Haz_sero_i_long_way(BaselineSeroStatus, Country, DATA, ProposedPARAMS, HOUSE, PhaseSeverity, serotype);
#else
				ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] = l_Int_Vac_Haz_sero_i(ProposedPARAMS, HOUSE, Country, serotype, BaselineSeroStatus, PhaseSeverity);
#endif
		}
	}
}
void Change_BS_BaseHazMult(Params_Struct& ProposedPARAMS, const Housekeeping_Struct& HOUSE, const DATA_struct& DATA)
{
#ifdef PRINT_WITHIN_CHANGE_FUNCTIONS
	std::cout << "Change_BS_BaseHazMult: Prop value = " << ProposedPARAMS.BS_BaseHazMults[HOUSE.Which_BS_BaseHazMult] << endl;
#endif
	int Country = 0; 
	int BS		= HOUSE.Which_BS_BaseHazMult; 
	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
		for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++)
		{
			Country = HOUSE.WhichCountries[countryindex]; 
			ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.l_BS_BaseHaz_Mults[PhaseSeverity][BS]] = l_BS_BaseHaz_Mult(Country, DATA, ProposedPARAMS, HOUSE, PhaseSeverity, BS);	///// cases	

			/////	Integrated Baseline Hazards
			for (int TrialArm = 0; TrialArm < HOUSE.NumTrialArms_IBH; TrialArm++)
				ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BS][TrialArm]] = l_Int_Base_Haz(Country, DATA, ProposedPARAMS, HOUSE, PhaseSeverity, BS, TrialArm);		///// survivors	

			/////	Integrated Vaccine Hazards
			for (int serotype = 0; serotype < HOUSE.N_STypes; serotype++)
				ProposedPARAMS.LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BS][serotype]] = l_Int_Vac_Haz_sero_i(ProposedPARAMS, HOUSE, Country, serotype, BS, PhaseSeverity);	///// survivors
		}
}

void Change_Param		(int param_no, DType &New_Param, const DATA_struct &DATA, const Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE)
{
	if (IsParamARelativeRisk(param_no, CurrentPARAMS.ParamNumbers) || IsParamARelativeRisk_Hosp(param_no, CurrentPARAMS.ParamNumbers))
	{
		int	PrevInf			= Find_PrevInf_From_K_Param			(param_no, HOUSE, CurrentPARAMS.ParamNumbers);
		int PhaseSeverity	= Find_PhaseSeverity_From_K_Param	(param_no, HOUSE, CurrentPARAMS.ParamNumbers, PrevInf);
		int serotype		= Find_Serotype_From_K_Param		(param_no, HOUSE, CurrentPARAMS.ParamNumbers, PrevInf, PhaseSeverity);

		std::vector<int> WhichCountries; 
		if (IsParamARelativeRisk(param_no, CurrentPARAMS.ParamNumbers))
			if (HOUSE.ModellingHospitalized && PhaseSeverity == PassiveSevere && HOUSE.ModelHosp_Indie_Ks)	WhichCountries = HOUSE.CYD_14_countries;	//// only do CYD-14 countries. 
			else																							WhichCountries = HOUSE.WhichCountries;		//// otherwise do all countries
		else if (IsParamARelativeRisk_Hosp(param_no, CurrentPARAMS.ParamNumbers))							WhichCountries = HOUSE.CYD_15_countries;	//// only do CYD-15 countries. 

		int Min_PhaseSeverity = PhaseSeverity, Max_PhaseSeverity = PhaseSeverity + 1; //// i.e. default is that you do one phase severity at a time in loop below. 
		if (HOUSE.PS_Ks_Multiply_AM_Ks && PhaseSeverity == ActiveMild) { Min_PhaseSeverity = ActiveMild; Max_PhaseSeverity = HOUSE.HowManyCaseCategories; } //// i.e. if changing Passive only this doesn't affect active, so default definition okay. 

		/*Order of if statement here really important. Don't be tempted to change. e.g. if SS_KAM_2 == false but SS_KPS_2 == true*/
		if (HOUSE.SSKs_FitMatrix[PhaseSeverity][PrevInf])	//// i.e. if fitting all serotypes for this K_PhaseSeverity_PrevInf, then do the below one serotype (this one) at a time. Automatically true for !SSKs too as default values of matrix all true (and changed only if SSKs). 
			for (int PhaseSeverityDummy = Min_PhaseSeverity; PhaseSeverityDummy < Max_PhaseSeverity; PhaseSeverityDummy++)
					Change_K(PhaseSeverityDummy, PrevInf, serotype		, DATA, CurrentPARAMS, ProposedPARAMS, HOUSE, WhichCountries);
		else												//// so SSKs == true & HOUSE.SSKs_FitMatrix[PhaseSeverity][PrevInf] == false. i.e. this K_PhaseSeverity_PrevInf is not serotype specific. Hence you will amend all serotype Ks for this PhaseSeverity and PrevInf at the same time. 
			for (int PhaseSeverityDummy = Min_PhaseSeverity; PhaseSeverityDummy < Max_PhaseSeverity; PhaseSeverityDummy++)
				for (int serotypeDummy = 0; serotypeDummy < HOUSE.N_STypes_Ks; serotypeDummy++)
					Change_K(PhaseSeverityDummy, PrevInf, serotypeDummy	, DATA, CurrentPARAMS, ProposedPARAMS, HOUSE, WhichCountries);
	}
	else	if (param_no == CurrentPARAMS.ParamNumbers.Hosp_K_multiplier) //// will only be true if (HOUSE.ModellingHospitalized !HOUSE.ModelHosp_Indie_Ks) 
	{
#ifdef PRINT_WITHIN_CHANGE_FUNCTIONS
		std::cout << "Change_Hosp_K_multiplier " << endl;
#endif
		for (int PrevInf = 0; PrevInf < HOUSE.Num_K_Params; PrevInf++)
			for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
				Change_K(PassiveSevere, PrevInf, serotype, DATA, CurrentPARAMS, ProposedPARAMS, HOUSE, HOUSE.CYD_15_countries); //// now do CYD-15 countries. 
	}
	else	if (IsParamAKnot		(param_no, CurrentPARAMS.ParamNumbers)	) 	//// if changing a knot
	{
		int knot	= Find_Knot_FromKnotParam	(param_no, HOUSE, CurrentPARAMS.ParamNumbers);
		int country = Find_Country_FromKnotParam(param_no, HOUSE, CurrentPARAMS.ParamNumbers, knot	);

		Change_Knot(country, knot, param_no, DATA, ProposedPARAMS, HOUSE);
	}
	else	if (IsParamAHistHaz		(param_no, CurrentPARAMS.ParamNumbers)	)
	{
		int country = Find_Country_FromHistHazParam(param_no, CurrentPARAMS.ParamNumbers);

		Change_HistHazard(country, DATA, CurrentPARAMS, ProposedPARAMS, HOUSE);
	}
	else	if (IsParamA_qval		(param_no, CurrentPARAMS.ParamNumbers)	)
	{
		int qval	= Find_qval_From_qParam		(param_no, HOUSE, CurrentPARAMS.ParamNumbers);
		int country = Find_Country_From_qParam	(param_no, HOUSE, CurrentPARAMS.ParamNumbers, qval);

		Change_qVal(qval, country, DATA, CurrentPARAMS, ProposedPARAMS, HOUSE);
	}
	else	if (IsParamA_rho		(param_no, CurrentPARAMS.ParamNumbers)	)
	{
		int serotype = Find_rho_From_rhoParam		(param_no, CurrentPARAMS.ParamNumbers.Min_Rho, HOUSE.N_STypes);
		int country	 = Find_Country_From_rhoParam	(param_no, CurrentPARAMS.ParamNumbers.Min_Rho, HOUSE.N_STypes, serotype);

		Change_Rho(country, serotype, DATA, CurrentPARAMS, ProposedPARAMS, HOUSE);
	}
	else	if (IsParamAnEfficacy	(param_no, CurrentPARAMS.ParamNumbers) || IsParamAn_atInf_Efficacy(param_no, CurrentPARAMS.ParamNumbers))
	{
		int BaselineSeroStatus	= Find_SeroStatus_FromEffParam		(param_no, HOUSE, CurrentPARAMS.ParamNumbers);
		int serotype			= Find_Serotype_FromEffParam		(param_no, HOUSE, CurrentPARAMS.ParamNumbers, BaselineSeroStatus);
		int PhaseSeverity		= Find_PhaseSeverity_From_EffParam	(param_no, HOUSE, CurrentPARAMS.ParamNumbers, BaselineSeroStatus, serotype);

		Change_Efficacy(DATA, CurrentPARAMS, ProposedPARAMS, HOUSE, BaselineSeroStatus, serotype, PhaseSeverity);
	}
	else	if (IsParamAWaningParam	(param_no, CurrentPARAMS.ParamNumbers))
	{
		if (HOUSE.AS_Waning_Homogeneous) //// change vaccine hazard and waning 
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
				Change_Waning(DATA, ProposedPARAMS, HOUSE, BaselineSeroStatus);
		else
		{
			int ParamType			= Find_ParamType_From_ASWaning_Param(param_no, HOUSE, CurrentPARAMS.ParamNumbers); 
			int BaselineSeroStatus	= Find_SeroStatus_From_Waning_Param	(param_no, HOUSE, CurrentPARAMS.ParamNumbers, ParamType); 

			Change_Waning(DATA, ProposedPARAMS, HOUSE, BaselineSeroStatus);
		}
	}
	else	if (IsParamAn_ASVE_Param(param_no, CurrentPARAMS.ParamNumbers))
	{
		int Type				= Find_Type_From_ASVE_Param		(param_no, HOUSE, CurrentPARAMS.ParamNumbers); //// either half life or coeff. 
		int BaselineSeroStatus	= Find_SeroStatus_FromASVE_Param(param_no, HOUSE, CurrentPARAMS.ParamNumbers, Type);

		if (HOUSE.ASVE_OnlyOneSeroStatus) BaselineSeroStatus = HOUSE.ASVE_BS;
		
		if (HOUSE.AS_VE_Homogeneous) //// change vaccine hazard and waning 
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
				Change_ASVE_Param(DATA, ProposedPARAMS, HOUSE, BaselineSeroStatus);
		else 
			Change_ASVE_Param(DATA, ProposedPARAMS, HOUSE, BaselineSeroStatus);
	}
	else	if (IsParamAn_ASHaz_Param(param_no, CurrentPARAMS.ParamNumbers))	Change_AS_Haz(DATA, ProposedPARAMS, HOUSE);
	else	if (IsParamAn_AS_Prime_Param(param_no, CurrentPARAMS.ParamNumbers))
	{
		int Type				= Find_ParamType_From_AS_Prime_Param	(param_no, HOUSE, CurrentPARAMS.ParamNumbers); //// either half life or coeff. 
		int BaselineSeroStatus	= Find_SeroStatus_From_AS_Prime_Param	(param_no, HOUSE, CurrentPARAMS.ParamNumbers, Type);

		Change_AS_Prime(BaselineSeroStatus, ProposedPARAMS, HOUSE, DATA);
	}
	else	if (param_no == CurrentPARAMS.ParamNumbers.BS_BaseHazMult)
	{
		Change_BS_BaseHazMult(ProposedPARAMS, HOUSE, DATA); 
	}
	else	std::cerr << "ChangeParam error: param_no not recognized " << endl;
}




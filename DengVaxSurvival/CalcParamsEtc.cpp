#include "Splines.h"
#include "Probability.h"
#include "CalcParamsEtc.h"

/*
	This script contains functions that calculate associated quantities of parameters. For example, "Efficacies" will affect SumRhoEffs etc; knots will affect (integrated) baseline hazard values. 
	In addition, these functions store various quantities that can be precalculated to speed computation. e.g. KPlusValues is an array calculated once for each unique value then used in Likelihood calculations as necessary. 
*/



bool CanEfficacyBeNegative				(int BaselineSeroStatus, const Housekeeping_Struct &HOUSE)
{
	 	 if (BaselineSeroStatus == SeroNeg && HOUSE.IntialParRanges.SNegEff_1[LowerBound] >= 0) return false;
	else if (BaselineSeroStatus == SeroPos && HOUSE.IntialParRanges.SPosEff_1[LowerBound] >= 0)	return false;
	else																						return true	; 
}
bool IsEfficacyNegative					(int BaselineSeroStatus, int PhaseSeverity, int AgeInYears, int serotype, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE) //// really asking IS the initial value of "fixed-time immunity" negative, whether there are age, phase/severity, or serotype, serostatus effects. Not asking whether it can in principle be negative. That is covered in function "CanEfficacyBeNegative" above. 
{
	DType EfficacyDummy = NULL; ///// need to know if "Efficacy" is greater or less than zero, which depends on a number of factors. Are there Age Effects? Are there Serotype effects as well? If so do they combine multiplicatively or additively?
	bool AdditiveScenario = HOUSE.SeroSpecificEfficacies && HOUSE.ASVE != Age_Option::INDEPENDENT && HOUSE.SSASVE_Additive;
	//// determine value of "efficacy". 
	if (AdditiveScenario)	EfficacyDummy = PARAMS.AgeEff_Mult[BaselineSeroStatus][AgeInYears] + PARAMS.Efficacies[PhaseSeverity][serotype][BaselineSeroStatus];
	else					EfficacyDummy = PARAMS.AgeEff_Mult[BaselineSeroStatus][AgeInYears] * PARAMS.Efficacies[PhaseSeverity][serotype][BaselineSeroStatus];

	if (EfficacyDummy >= 0)		return false; 
	else						return true; 
}
DType Calc_AggTransImmunity				(int BaselineSeroStatus, int serotype,	int AgeInYears, int PhaseSeverity,									const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	//// Guard against calling this function for serotype > 0 when not doing SeroSpecificEfficacies. Could happen with HOUSE.SeroSpecific_K_values. 
	if (!HOUSE.SeroSpecificEfficacies) serotype = 0;

	//// This provides the aggregate efficacy (for a single serotype, i.e. not when summing over serotypes e.g. for IntBaseHaz) WITHOUT asking for the absolute value.
	//// The overload below (and one that is used most often), has the absolute value bit essentially tacked onto this. You need the latter when doing integrated hazards / augmentation etc.. 
	//// The reason you separate it out is that in likelihood function l_WaningEfficacy (and associated bit in Augmentation), need to know whether "Efficacy" was negative as it affects whether you add 1 to the waning multiplier 
	//// (remember when HOUSE.AltWaneEffNeg your Waning multiplier is exp(-t/tau) if VE>=0 but is exp(-t/tau) - 1 if VE<0. 
	DType AggregateEff		= 0; 
	DType BS_SerotypePart	= PARAMS.Efficacies[PhaseSeverity][serotype][BaselineSeroStatus];
	DType AgePart			= PARAMS.AgeEff_Mult[BaselineSeroStatus][AgeInYears];

	if (HOUSE.AdditiveSSASVEs)	AggregateEff = BS_SerotypePart + AgePart;
	else						AggregateEff = BS_SerotypePart * AgePart;

	return(AggregateEff);
}
DType Calc_AggTransImmunity				(int BaselineSeroStatus, int serotype,	int AgeInYears, int PhaseSeverity, EffNegWane_Option ENW_Option,	const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	//// See note in overload of this function why you separate out calculating the aggregate transient immunity and whether or not to take the absolute value. This function calls the former and performs the latter. 
	DType AggregateEff = Calc_AggTransImmunity(BaselineSeroStatus, serotype, AgeInYears, PhaseSeverity, PARAMS, HOUSE);
		 if (ENW_Option == EffNegWane_Option::FROM_ZERO	)						AggregateEff = abs(AggregateEff);
	else if (ENW_Option == EffNegWane_Option::NO_WANE && AggregateEff < 0)		AggregateEff = 0;					
	return(AggregateEff);
}
DType Calc_SumRhoAggregateTransImmunity	(int BaselineSeroStatus, int country,	int AgeInYears, int PhaseSeverity, EffNegWane_Option ENW_Option,	const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	//// Only ever used if HOUSE.AdditiveSSASVEs and applied to INTEGRATED vaccine hazard. 
	DType AggregateSummedEff = 0;
	for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
		AggregateSummedEff += PARAMS.rhos[country][serotype] * Calc_AggTransImmunity(BaselineSeroStatus, serotype, AgeInYears, PhaseSeverity, ENW_Option, PARAMS, HOUSE);

	return AggregateSummedEff; 
}
DType Calc_SumRho_Ks_AggTransImmunity	(int BaselineSeroStatus, int country,	int AgeInYears, int PhaseSeverity, EffNegWane_Option ENW_Option,	const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	//// Only ever used if HOUSE.AdditiveSSASVEs AND HOUSE.Sero and applied to INTEGRATED vaccine hazard. 
	DType Aggregate_Sum_Rho_Age_Sero_Eff_K = 0;
	for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
		Aggregate_Sum_Rho_Age_Sero_Eff_K +=	PARAMS.rhos[country][serotype]																	* //// rho 
								Calc_AggTransImmunity(BaselineSeroStatus, serotype, AgeInYears, PhaseSeverity, ENW_Option, PARAMS, HOUSE)	* //// TI_bd(age): i.e. transient immunity as an (additive or multiplicative) function of age and serotype. 
								Choose_K(VaccineGroup, BaselineSeroStatus, country, AgeInYears, PhaseSeverity, serotype, PARAMS, HOUSE)		; //// serotype specific K. 

	return Aggregate_Sum_Rho_Age_Sero_Eff_K;
}
DType DetermineAltWaneEffBaseHazMult	(int BaselineSeroStatus, int PhaseSeverity, int AgeInYears, int serotype, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	//// **** //// **** //// **** USED for HOUSE.AltWaneEffNeg
	//// Determines whether to multiply basline hazard by (1 - VE). 
	DType EfficacyDummy = 0, OneMinusEffDummy = 0; ///// need to know if "Efficacy" is greater or less than zero, which depends on a number of factors. Are there Age Effects? Are there Serotype effects as well? If so do they combine multiplicatively or additively?
	EfficacyDummy		= Calc_AggTransImmunity(BaselineSeroStatus, serotype, AgeInYears, PhaseSeverity, PARAMS, HOUSE);
	if (EfficacyDummy >= 0)		OneMinusEffDummy = 1;
	else						OneMinusEffDummy = (1 - EfficacyDummy);
	return OneMinusEffDummy;
}
void Calc_SumRhoEffNegs_c_BS_PS			(int country, int BaselineSeroStatus, int PhaseSeverity, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	///// NOTE: This depends on AgeEffs, Efficacies and rhos etc. being properly calculated already. 

	if ((HOUSE.EffNegWane == EffNegWane_Option::FROM_ZERO || HOUSE.EffNegWane == EffNegWane_Option::NO_WANE))
	{
		//// If efficacy can in principle be negative be negative, need to run summation of rho_cd^*
		if (CanEfficacyBeNegative(BaselineSeroStatus, HOUSE) & !HOUSE.BaselinePartition) //// if BaselinePartition then sum(rho_cd) != 1. 
		{
			DType RhoCDStar_Dummy = 0;
			DType EfficacyDummy = NULL; ///// need to know if "Efficacy" is greater or less than zero, which depends on a number of factors. Are there Age Effects? Are there Serotype effects as well? If so do they combine multiplicatively or additively?

			for (int AgeInYears = 0; AgeInYears < HOUSE.HowManyAges; AgeInYears++)
			{
				PARAMS.SumRhoEffNegs[PhaseSeverity][country][AgeInYears][BaselineSeroStatus] = 0; //// reset to zero
				for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
				{
					RhoCDStar_Dummy = PARAMS.rhos[country][serotype] * DetermineAltWaneEffBaseHazMult(BaselineSeroStatus, PhaseSeverity, AgeInYears, serotype, PARAMS, HOUSE); //// will be either rho * 1	OR	 rho * (1 - VE_bd(alpha))
					PARAMS.SumRhoEffNegs[PhaseSeverity][country][AgeInYears][BaselineSeroStatus] += RhoCDStar_Dummy;
				}
			}
		}
		else //// If efficacy cannot be negative, this should always be sum(rho_cd) = 1. 
			for (int AgeInYears = 0; AgeInYears < HOUSE.HowManyAges; AgeInYears++)
				PARAMS.SumRhoEffNegs[PhaseSeverity][country][AgeInYears][BaselineSeroStatus] = 1;
	}
}

DType Choose_K					(int TrialArm, int BaselineSeroStatus, int country, int Age, int PhaseSeverity, int serotype, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	//// Guard against calling this function for serotype > 0 when not doing SeroSpecific_K_values. Could happen with HOUSE.SeroSpecificEfficacies. 
	if (!HOUSE.SeroSpecific_K_values) serotype = 0; 

	DType K		= 0;
	if (HOUSE.ModelVariant == VAC_SILENT)
		if (TrialArm == ControlGroup)
			if (BaselineSeroStatus == SeroNeg)			K = PARAMS.K_s[country][PhaseSeverity][0][serotype];
			else 										K = K_seropos(PARAMS.ParamVec[PARAMS.ParamNumbers.Min_HHaz + country], Age, PARAMS.K_s[country][PhaseSeverity][1][serotype], PARAMS.K_s[country][PhaseSeverity][2][serotype]);
		else											K = PARAMS.K_s[country][PhaseSeverity][BaselineSeroStatus + 1][serotype];
	else if (HOUSE.ModelVariant == K_SEROPOS)
			if (BaselineSeroStatus == SeroNeg)			K = PARAMS.K_s[country][PhaseSeverity][0][serotype];
			else 										K = K_seropos(PARAMS.ParamVec[PARAMS.ParamNumbers.Min_HHaz + country], Age, PARAMS.K_s[country][PhaseSeverity][1][serotype], PARAMS.K_s[country][PhaseSeverity][2][serotype]);
	else if (HOUSE.ModelVariant == AS_PRIME)
		if (TrialArm == ControlGroup)
		{
			if (BaselineSeroStatus == SeroNeg)			K = PARAMS.K_s[country][PhaseSeverity][0][serotype];
			else 										K = K_seropos(PARAMS.ParamVec[PARAMS.ParamNumbers.Min_HHaz + country], Age, PARAMS.K_s[country][PhaseSeverity][1][serotype], PARAMS.K_s[country][PhaseSeverity][2][serotype]);
		}
		else /// i.e. vaccine group
		{
			if (BaselineSeroStatus == SeroNeg)			K = PARAMS.Ks_Prime[country][PhaseSeverity][SeroNeg][serotype][Age];
			else 										K = K_seropos(PARAMS.ParamVec[PARAMS.ParamNumbers.Min_HHaz + country], Age, PARAMS.K_s[country][PhaseSeverity][1][serotype], PARAMS.K_s[country][PhaseSeverity][2][serotype]);
		}
	else if (HOUSE.ModelVariant == SIMPLE_NUMERICAL)	K = PARAMS.K_s[country][PhaseSeverity][BaselineSeroStatus][serotype]; //// for SIMPLE_NUMERICAL K multipliers the same between trial arms. And there is no K2 so can just use BaselineSeroStatus as the same as PrevInf. 

	return K; 
}
DType Choose_Sum_Rho_K			(int TrialArm, int BaselineSeroStatus, int country, int Age, int PhaseSeverity, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	DType SumRhoK = 0;
	for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++) 
		SumRhoK += PARAMS.rhos[country][serotype] * Choose_K(TrialArm, BaselineSeroStatus, country, Age, PhaseSeverity, serotype, PARAMS, HOUSE); 
	return SumRhoK; 
}
DType Choose_Survive_K			(int TrialArm, int BaselineSeroStatus, int country, int Age, int PhaseSeverity, const Params_Struct &PARAMS, char ModelVariant)
{
	DType SurviveK = 0;
	if (ModelVariant == VAC_SILENT)
		if (TrialArm == ControlGroup)
			if (BaselineSeroStatus == SeroNeg)		SurviveK = PARAMS.SumRhoKs[country][PhaseSeverity][0];
			else 									SurviveK = PARAMS.Meta_KplusValues[No_Effs][country][Age][PhaseSeverity];
		else										SurviveK = PARAMS.SumRhoKs[country][PhaseSeverity][BaselineSeroStatus + 1];
	else if (ModelVariant == K_SEROPOS)	
			if (BaselineSeroStatus == SeroNeg)		SurviveK = PARAMS.SumRhoKs[country][PhaseSeverity][0];			//// for K_SEROPOS K multipliers the same between trial arms. 
			else 									SurviveK = PARAMS.KplusValues[country][Age][PhaseSeverity];		//// for K_SEROPOS K multipliers the same between trial arms. 
	else if (ModelVariant == AS_PRIME)
		if (TrialArm == ControlGroup)
		{
			if (BaselineSeroStatus == SeroNeg)		SurviveK = PARAMS.SumRhoKs[country][PhaseSeverity][0];
			else 									SurviveK = PARAMS.Meta_KplusValues[No_Effs][country][Age][PhaseSeverity];
		}
		else /// i.e. vaccine group
		{
			if (BaselineSeroStatus == SeroNeg)		SurviveK  = PARAMS.SumRhoK0_SNeg_Primes[country][PhaseSeverity][Age];
			else 									SurviveK  = PARAMS.KPlusPrimeValues[country][Age][PhaseSeverity];
		}
	else if (ModelVariant == SIMPLE_NUMERICAL)		SurviveK = PARAMS.SumRhoKs[country][PhaseSeverity][BaselineSeroStatus]; //// for SIMPLE_NUMERICAL K multipliers the same between trial arms. And there is no K2 so you can just use BaselineSeroStatus as the same as PrevInf. 


	return SurviveK;
}
DType Choose_BaseHazMult		(DType K_singlesero,	int TrialArm,	int BaselineSeroStatus, int country, int Age, int PhaseSeverity, int serotype,	const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	//// When using common factors (usually Ks), can call this function within likelihood #pragma loops with K_singlesero = 1 so that usual machinery will work. 
	DType BaseHazMult_SingleSerotype = K_singlesero;

	if (HOUSE.AdjHaz) 
		BaseHazMult_SingleSerotype *= PARAMS.BS_BaseHazMults[BaselineSeroStatus];
	if (HOUSE.AS_Haz != Age_Option::INDEPENDENT)
		BaseHazMult_SingleSerotype *= PARAMS.AgeHaz_Mult[Age];
	if (HOUSE.SeroSpecificEfficacies || HOUSE.SeroSpecific_K_values) 
		BaseHazMult_SingleSerotype *= PARAMS.rhos[country][serotype];
	if ((HOUSE.EffNegWane == EffNegWane_Option::FROM_ZERO || HOUSE.EffNegWane == EffNegWane_Option::NO_WANE) && TrialArm == VaccineGroup)
		BaseHazMult_SingleSerotype *= DetermineAltWaneEffBaseHazMult(BaselineSeroStatus, PhaseSeverity, Age, serotype, PARAMS, HOUSE); //// will be either rho * 1	OR	 rho * (1 - VE_bd(alpha))

	return BaseHazMult_SingleSerotype; 
}
DType Choose_BaseHazMult		(						int TrialArm,	int BaselineSeroStatus, int country, int Age, int PhaseSeverity, int serotype,	const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	///// Overload of the above, which assumes K_singlesero already calculated
	DType K_singlesero					= Choose_K			(				TrialArm, BaselineSeroStatus, country, Age, PhaseSeverity, serotype, PARAMS, HOUSE); 
	DType BaseHazMult_SingleSerotype	= Choose_BaseHazMult(K_singlesero,	TrialArm, BaselineSeroStatus, country, Age, PhaseSeverity, serotype, PARAMS, HOUSE); 
	return BaseHazMult_SingleSerotype;
}
DType Choose_SumBaseHazMult		(int TrialArm,							int BaselineSeroStatus, int country, int Age, int PhaseSeverity,				const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	DType SumRhoBaseHazMult = 0;
	for (int serotype = 0; serotype < HOUSE.N_STypes; serotype++)
		SumRhoBaseHazMult += Choose_BaseHazMult(TrialArm, BaselineSeroStatus, country, Age, PhaseSeverity, serotype, PARAMS, HOUSE); //// includes rhos
	return SumRhoBaseHazMult;
}

DType Choose_VacHazMult			(DType K_singlesero,					int BaselineSeroStatus, int country, int Age, int PhaseSeverity, int serotype,	const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	//// When using common factors (usually Ks), can call this function within likelihood #pragma loops with K_singlesero = 1 so that usual machinery will work. 
	DType VacHazMult_SingleSerotype = K_singlesero;
	VacHazMult_SingleSerotype		*= Calc_AggTransImmunity(BaselineSeroStatus, serotype, Age, PhaseSeverity, HOUSE.EffNegWane, PARAMS, HOUSE); //// here you use ENW_Option, (will either not amend aggregate transient immunity (DEFAULT), will return abs value (FROM_ZERO), or will return zero (NO_WANE). In the latter case there is no vaccine hazard as there is no waning). 

	if (HOUSE.AdjHaz)
		VacHazMult_SingleSerotype	*= PARAMS.BS_BaseHazMults[BaselineSeroStatus];
	if (HOUSE.AS_Haz != Age_Option::INDEPENDENT)
		VacHazMult_SingleSerotype	*= PARAMS.AgeHaz_Mult[Age];
	if (HOUSE.SeroSpecificEfficacies || HOUSE.SeroSpecific_K_values)
		VacHazMult_SingleSerotype	*= PARAMS.rhos[country][serotype];

	return VacHazMult_SingleSerotype;
}
DType Choose_SumVacHazMult		(DType K_singlesero,					int BaselineSeroStatus, int country, int Age, int PhaseSeverity,				const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	//// For the sum of integrated vaccine hazards, can often use common factors that are not specific to serotype, age (i.e. a person) or both. 
	//// For example, if SSKs but not SSVEs, there are no age effects in the Ks (if the model is VAC_SILENT), and hence this can be taken out of the l_Int_Vac_Haz #pragma loop. 
	//// Or if both SSVEs and SSKs, but no additive age effects, then can also take out common factors of Ks and VEs from #pragma loop. 
	//// In this case however (and others like it), apply the common factors to multipe serotype-specific likelihood components, which various change and amend functions rely upon. 
	//// Hence this solution. By allowing this function to input it's own value of K_singlesero = 1, can still use common factors and not ruin the rest of time-saving code. 
	//// Essentially, unless doing AdditiveSSASVEs and SSKs, for the LIKELIHOOD (and only the likelihood), call this function with K_singlesero = 1, and apply various K terms elsewhere. 

	DType SumVacHazMult = 0;
	for (int serotype = 0; serotype < HOUSE.N_STypes; serotype++)
		SumVacHazMult += Choose_VacHazMult(K_singlesero,	BaselineSeroStatus, country, Age, PhaseSeverity, serotype, PARAMS, HOUSE); //// includes rhos so don't include here. 
	return SumVacHazMult;
}
DType Choose_VacHazMult			(										int BaselineSeroStatus, int country, int Age, int PhaseSeverity, int serotype,	const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	///// Overload of the above, which assumes K_singlesero already calculated (and thus might save recalculation). 
	DType K_singlesero					= Choose_K			(VaccineGroup, BaselineSeroStatus, country, Age, PhaseSeverity, serotype, PARAMS, HOUSE); //// obviously multiplier of vaccine hazard only needed for VaccineGroup. 
	DType VacHazMult_SingleSerotype		= Choose_VacHazMult	(K_singlesero, BaselineSeroStatus, country, Age, PhaseSeverity, serotype, PARAMS, HOUSE);
	
	return VacHazMult_SingleSerotype;
}
DType Choose_SumVacHazMult		(										int BaselineSeroStatus, int country, int Age, int PhaseSeverity,				const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	DType SumVacHazMult = 0;
	for (int serotype = 0; serotype < HOUSE.N_STypes; serotype++)
		SumVacHazMult += Choose_VacHazMult(					BaselineSeroStatus, country, Age, PhaseSeverity, serotype, PARAMS, HOUSE); //// includes rhos so don't include here. 
	return SumVacHazMult;
}




void Calc_ParamSeroPrevs_country					(int &country, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	//// calculate P(SNeg) and P(SPos) for each age. 
	for (int age = 0; age < HOUSE.HowManyAges; age++)
	{
		PARAMS.SeroPrevs[NonLogIndex][country][SeroNeg][age] = exp(-PARAMS.ParamVec[PARAMS.ParamNumbers.Min_HHaz + country] * (DType)age);			//// P(SNeg)
		PARAMS.SeroPrevs[NonLogIndex][country][SeroPos][age] = 1 - PARAMS.SeroPrevs[NonLogIndex][country][SeroNeg][age];							//// P(SPos) = 1 - P(SNeg)
	}
	//// calculate log(P(SNeg)) and log(P(SPos)) for each age. 
	for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
		for (int age = 0; age < HOUSE.HowManyAges; age++)
			PARAMS.SeroPrevs[LogIndex][country][BaselineSeroStatus][age] = log(PARAMS.SeroPrevs[NonLogIndex][country][BaselineSeroStatus][age]);
}
void Calc_KplusValues_country_PhaseSeverity			(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int country, int PhaseSeverity)
{
	DType AgeDummy, hci = PARAMS.ParamVec[PARAMS.ParamNumbers.Min_HHaz + country];
	DType K1, K2;

	DType AmountToAdd = 0;
	if (HOUSE.ModelVariant != SIMPLE_NUMERICAL)
		for (int age = 0; age < HOUSE.HowManyAges; age++)
		{
			AgeDummy = age;

			//// reset value to zero, then sum over serotypes. 
			PARAMS.KplusValues[country][age][PhaseSeverity] = 0; 

			if (HOUSE.ModelVariant == K_SEROPOS && HOUSE.SeroSpecific_K_values && HOUSE.SeroSpecificEfficacies)
				PARAMS.Meta_KplusValues[With_Effs][country][age][PhaseSeverity] = 0;

			for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
			{
				K1 = PARAMS.K_s[country][PhaseSeverity][1][serotype]; 
				K2 = PARAMS.K_s[country][PhaseSeverity][2][serotype];

				AmountToAdd = K_seropos(hci, AgeDummy, K1, K2);

				if (HOUSE.SeroSpecific_K_values)
					AmountToAdd *= PARAMS.rhos[country][serotype];

				PARAMS.KplusValues[country][age][PhaseSeverity] += AmountToAdd;

				if (HOUSE.ModelVariant == K_SEROPOS && HOUSE.SeroSpecific_K_values && HOUSE.SeroSpecificEfficacies)
					PARAMS.Meta_KplusValues[With_Effs][country][age][PhaseSeverity] += AmountToAdd * PARAMS.Efficacies[PhaseSeverity][serotype][SeroPos];
			}
		}
}
void Calc_KplusPrimes_country_PhaseSeverity			(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int country, int PhaseSeverity)
{
	if (HOUSE.ModelVariant == AS_PRIME)
	{
		DType AgeDummy, K1, K2, AmountToAdd = 0, hci = PARAMS.ParamVec[PARAMS.ParamNumbers.Min_HHaz + country];
		for (int age = 0; age < HOUSE.HowManyAges; age++)
		{
			AgeDummy = age;

			//// reset value to zero, then sum over serotypes. 
			PARAMS.KPlusPrimeValues[country][age][PhaseSeverity] = 0;

			for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
			{
				K1 = PARAMS.Ks_Prime[country][PhaseSeverity][1][serotype][age];
				K2 = PARAMS.K_s		[country][PhaseSeverity][2][serotype];

				AmountToAdd = K_seropos(hci, AgeDummy, K1, K2);

				if (HOUSE.SeroSpecific_K_values)
					AmountToAdd *= PARAMS.rhos[country][serotype];

				PARAMS.KPlusPrimeValues[country][age][PhaseSeverity] += AmountToAdd;
			}
		}
	}
}
void Calc_K_Primes_BS_PhaseSeverity_country_serotype(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int BaselineSeroStatus, int PhaseSeverity, int country, int serotype)
{
	if (BaselineSeroStatus == SeroNeg)
		for (int age = 0; age < HOUSE.HowManyAges; age++)
			PARAMS.Ks_Prime[country][PhaseSeverity][SeroNeg][serotype][age] =
				K0_Prime(age,	PARAMS.K_s[country][PhaseSeverity][1][serotype]				, 
								PARAMS.K_s[country][PhaseSeverity][2][serotype]				, 
								PARAMS.AS_Priming_Params[SeroNeg][Prop_Index]				,
								PARAMS.AS_Priming_Params[SeroNeg][AS_Priming_Rate_index]	);
	else if (BaselineSeroStatus == SeroPos)
		for (int age = 0; age < HOUSE.HowManyAges; age++)
			PARAMS.Ks_Prime[country][PhaseSeverity][SeroPos][serotype][age] =
				K1_Prime(age,	PARAMS.K_s[country][PhaseSeverity][1][serotype]				, 
								PARAMS.K_s[country][PhaseSeverity][2][serotype]				, 
								PARAMS.AS_Priming_Params[SeroPos][Prop_Index]				,
								PARAMS.AS_Priming_Params[SeroPos][AS_Priming_Rate_index]	);
}
void Calc_SumRhos_country							(int &country, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	PARAMS.SumRhos[country] = 0;
	for (int serotype = 0; serotype < HOUSE.N_STypes; serotype++) PARAMS.SumRhos[country] += PARAMS.rhos[country][serotype];
}
void Calc_SumRhoK0_SNeg_Primes_country_PhaseSeverity(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int country, int PhaseSeverity)
{
	DType AmountToAdd; 

	for (int Age = 0; Age < HOUSE.HowManyAges; Age++)
	{
		//// reset to zero. Then add afresh. 
		AmountToAdd = 0;

		PARAMS.SumRhoK0_SNeg_Primes[country][PhaseSeverity][Age] = 0;
		for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
		{
			AmountToAdd = PARAMS.Ks_Prime[country][PhaseSeverity][SeroNeg][serotype][Age];

			if (HOUSE.SeroSpecific_K_values) AmountToAdd *= PARAMS.rhos[country][serotype];

			PARAMS.SumRhoK0_SNeg_Primes[country][PhaseSeverity][Age] += AmountToAdd;
		}
	}
}
void Calc_RhosFrom_qParams							(int &country, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	DType qc0 = PARAMS.ParamVec[PARAMS.ParamNumbers.Min_qval + (HOUSE.N_STypes - 1) * country		]; 
	DType qc1 = PARAMS.ParamVec[PARAMS.ParamNumbers.Min_qval + (HOUSE.N_STypes - 1) * country + 1	];
	DType qc2 = PARAMS.ParamVec[PARAMS.ParamNumbers.Min_qval + (HOUSE.N_STypes - 1) * country + 2	];

	PARAMS.rhos[country][Serotype_1] =								qc0		;	//// change in qc0 affects all serotypes	in country c. 
	PARAMS.rhos[country][Serotype_2] =					qc1		* (1 - qc0)	;	//// change in qc1 affects serotypes 2,3,4	in country c. 
	PARAMS.rhos[country][Serotype_3] =	qc2			* (1 - qc1)	* (1 - qc0)	;	//// change in qc2 affects serotypes 3,4	in country c. 
	PARAMS.rhos[country][Serotype_4] = (1 - qc2)	* (1 - qc1)	* (1 - qc0)	;	//// change in qc2 affects serotypes 3,4	in country c.  

	////// check. 
	//DType SumRhos = 0; 
	//for (int serotype = 0; serotype < HOUSE.N_STypes; serotype++) SumRhos += PARAMS.rhos[country][serotype]; 
	//if (!(SumRhos == 1)) std::cerr << "CalcRhosFrom_qParams error: sum != 1: diff = " << 1 - SumRhos << endl; 
}
void Calc_SumRhoEffs_c_BS_PS						(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int country, int BaselineSeroStatus, int PhaseSeverity) //// must call this (or related) when changing any efficacy (inc. when not serospecific) or rho parameter
{
	DType AmountToAdd = 0; 
	PARAMS.SumRhoEffs[PhaseSeverity][country][BaselineSeroStatus] = 0; //// first reset. 
	for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
	{
		AmountToAdd = PARAMS.Efficacies[PhaseSeverity][serotype][BaselineSeroStatus]; 
		if (HOUSE.SeroSpecificEfficacies) AmountToAdd *= PARAMS.rhos[country][serotype]; 

		PARAMS.SumRhoEffs[PhaseSeverity][country][BaselineSeroStatus] += AmountToAdd;
	}
}
void Calc_SumRhoK_country_PhaseSeverity_PrevInf		(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int &country, int PhaseSeverity, int PrevInf)
{
	//// reset to zero. Then add afresh. 
	PARAMS.SumRhoKs[country][PhaseSeverity][PrevInf] = 0;
	DType AmountToAdd = 0;

	for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
	{
		AmountToAdd = PARAMS.K_s[country][PhaseSeverity][PrevInf][serotype]; 

		if (HOUSE.SeroSpecific_K_values) AmountToAdd *= PARAMS.rhos[country][serotype];
		
		PARAMS.SumRhoKs[country][PhaseSeverity][PrevInf] += AmountToAdd;
	}
}
void Calc_SumRhoEffKs_c_BS_PS						(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int &country, int PhaseSeverity, int BaselineSeroStatus)
{
	if (HOUSE.SeroSpecificEfficacies && HOUSE.SeroSpecific_K_values)
	{
		//// reset to zero. Then add afresh. 
		PARAMS.SumRhoEffKs[country][PhaseSeverity][BaselineSeroStatus] = 0;

		int PrevInf = BaselineSeroStatus;
		if (HOUSE.ModelVariant == VAC_SILENT) PrevInf++;

		for (int serotype = 0; serotype < HOUSE.N_STypes; serotype++)
			PARAMS.SumRhoEffKs[country][PhaseSeverity][BaselineSeroStatus] += 
			PARAMS.rhos[country][serotype] * PARAMS.Efficacies[PhaseSeverity][serotype][BaselineSeroStatus] * PARAMS.K_s[country][PhaseSeverity][PrevInf][serotype];
	}
}
void Calc_BaseHazValues								(int &country, Params_Struct &PARAMS, const DATA_struct &DATA, const Housekeeping_Struct &HOUSE)
{
	//// Change Baseline Hazard Values. 
#pragma omp parallel for schedule(static,1) 
	// Can parallelize over  BaselineHazard array (each entry independent of the other), but you are incrementing Cumulative and then storing it before all threads finish. 
	for (int thread_no = 0; thread_no < HOUSE.max_threads; thread_no++)
		for (int CalendarDay = thread_no; CalendarDay <= DATA.NumCalendarDaysFollowUp[country]; CalendarDay += HOUSE.max_threads) //// Need only amend integrated hazard values that are after this value so this is wasteful but can get tricky so err on side of caution. 
			PARAMS.BaselineHazardValues[country][CalendarDay] = BaselineHazard((CalendarDay * HOUSE.TimeInterval), country, PARAMS, HOUSE);
}
void Calc_IntBaseHazValues							(int &country, Params_Struct &PARAMS, const DATA_struct &DATA, const Housekeeping_Struct &HOUSE)
{
	// Can parallelize over the BaselineHazard array (each entry independent of the other), but you are incrementing Cumulative and then storing it before all threads finish. 
	DType IntBaseHazCumulative = 0;
	for (int CalendarDay = 0; CalendarDay <= DATA.NumCalendarDaysFollowUp[country]; CalendarDay++)
	{
		// add to cumulative / approximate integrated hazard. 
		PARAMS.IntBaseHazLookUp[country][CalendarDay]	= IntBaseHazCumulative;
		IntBaseHazCumulative							+= PARAMS.BaselineHazardValues[country][CalendarDay] * HOUSE.TimeInterval;	// tempting to multiply by TimeInterval once at end as common factor, but as TimeInterval interval required for each day don't do this.
	}
}
void Calc_IVH_Values_c_BS_PS						(int country, int BaselineSeroStatus, int PhaseSeverity, Params_Struct &PARAMS, const DATA_struct &DATA, const Housekeeping_Struct &HOUSE)  //// Tempting to put this in Initialize_Params, but put it in main as otherwise you interupt following flow: Initialize params, use historical hazards to populate strata sets, calculate IntVacHaz array. . 
{
	//// Note this is different from CalcIntBaseHazardValues function. That function calculates a look up table of (integrated) baseline Hazard values, from which it is easy to calculate a particular patient's integrated baseline hazard
	//// It is not so easy to do this for the vaccine hazard, and so this function calculates a patient's vaccine hazards (i.e. int_0^t(lambda(s) * Waning(s)))
	int Stratum_index = HOUSE.Strata[VaccineGroup][BaselineSeroStatus]; 

#pragma omp parallel for schedule(static,1) 
	for (int thread_no = 0; thread_no < HOUSE.max_threads; thread_no++)
	{
		DType AmountToAdd;
		int person_i, FirstDayVacHazard;
		for (int i = thread_no; i < DATA.Set_Array[country][Stratum_index].size(); i += HOUSE.max_threads)
		{
			person_i = DATA.Set_Array[country][Stratum_index][i];
			
			//// First day usually start of this trial phase (or start of trial if HOUSE.MildAndSevere == TREATED_SEPARATELY)
			//// However, if all doses are required for this serostatus before vaccine has any effect, only integrate the vaccine hazard from the third dose onwards (baseline hazard in the vaccine group, calculated in l_Int_Base_Haz function, remains unchanged)
			//// Only thing left is to be careful, as you would still integrate from the start of the passive phase, but when HOUSE.MildAndSevere == TREATED_SEPARATELY, you always integrate from start of follow up, hence FirstDayVacHazard would be independent of PhaseSeverity. 
			FirstDayVacHazard = DATA.FollowUp[START][PhaseSeverity][person_i];

			if (HOUSE.AllDosesRequired_AG_BS[DATA.AgeGroup1[person_i]][BaselineSeroStatus])
				if (HOUSE.MildAndSevere == TREATED_SEPARATELY || PhaseSeverity == ActiveMild)
					FirstDayVacHazard = DATA.ThirdDose[person_i];

			PARAMS.IVH_vals[person_i][PhaseSeverity] = 
				IntVacHaz(	FirstDayVacHazard									,
							DATA.FollowUp[END	][PhaseSeverity	]	[person_i]	,
							DATA.FollowUp[START	][ActiveMild	]	[person_i]	,	// Note that Dose_1_Day in IntVacHaz has to be start of active phase, NEVER DATA.Start_PassiveSevere
							DATA.SecondDose							[person_i]	,
							DATA.ThirdDose							[person_i]	,
							DATA.ai_s								[person_i]	,
							BaselineSeroStatus, country, PARAMS, HOUSE)			;	
		}
	}
}


///// The functions below are all wrappers for those above
void Calc_KplusValues_PhaseSeverity			(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int PhaseSeverity)
{
	for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++)
		Calc_KplusValues_country_PhaseSeverity(PARAMS, HOUSE, HOUSE.WhichCountries[countryindex], PhaseSeverity);
}
void Calc_KplusValues_country				(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int &country)
{
	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)		
		Calc_KplusValues_country_PhaseSeverity(PARAMS, HOUSE, country, PhaseSeverity);
}
void Calc_KplusValues_All					(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++)
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)	
			Calc_KplusValues_country_PhaseSeverity(PARAMS, HOUSE, HOUSE.WhichCountries[countryindex], PhaseSeverity);
}
void Calc_KplusPrimes_PhaseSeverity			(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int PhaseSeverity)
{
	for (int country = 0; country < HOUSE.TotalCountries; country++)								
		Calc_KplusPrimes_country_PhaseSeverity(PARAMS, HOUSE, country, PhaseSeverity);
}
void Calc_KplusPrimes_country				(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int country)
{
	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)		
		Calc_KplusPrimes_country_PhaseSeverity(PARAMS, HOUSE, country, PhaseSeverity);
}
void Calc_KplusPrimes_All					(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++)
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)	
			Calc_KplusPrimes_country_PhaseSeverity(PARAMS, HOUSE, HOUSE.WhichCountries[countryindex], PhaseSeverity);
}
void Calc_K_Primes_ALL						(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
				for (int country = 0; country < HOUSE.TotalCountries; country++) //// loop over all countries, not WhichCountries. Pointers set so that pointer to country 1,2,3,4 same as to country 0 etc.
					Calc_K_Primes_BS_PhaseSeverity_country_serotype(PARAMS, HOUSE, BaselineSeroStatus, PhaseSeverity, country, serotype);
}
void Calc_K_Primes_BS						(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int BaselineSeroStatus)
{
	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
		for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++)
			for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
				Calc_K_Primes_BS_PhaseSeverity_country_serotype(PARAMS, HOUSE, BaselineSeroStatus, PhaseSeverity, HOUSE.WhichCountries[countryindex], serotype);
}
void Calc_K_Primes_serotype					(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int serotype)
{
	for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++)
				for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
					Calc_K_Primes_BS_PhaseSeverity_country_serotype(PARAMS, HOUSE, BaselineSeroStatus, PhaseSeverity, HOUSE.WhichCountries[countryindex], serotype);
}
void Calc_K_Primes_PhaseSeverity_serotype	(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int PhaseSeverity, int serotype)
{
	for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
		for (int country = 0; country < HOUSE.TotalCountries; country++)
			for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
				Calc_K_Primes_BS_PhaseSeverity_country_serotype(PARAMS, HOUSE, BaselineSeroStatus, PhaseSeverity, country, serotype);
}
void Calc_SumRhoK0_SNeg_Primes_ALL			(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++)
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			Calc_SumRhoK0_SNeg_Primes_country_PhaseSeverity(PARAMS, HOUSE, HOUSE.WhichCountries[countryindex], PhaseSeverity);
}
void Calc_SumRhos_ALL						(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	for (int country = 0; country < HOUSE.TotalCountries; country++)		Calc_SumRhos_country(country, PARAMS, HOUSE);
}
void Calc_SumRhoEffs_country							(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int &country)
{
	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)	//// Note: NOT  <HOUSE.NumEffsPer_BS_And_Serotype. Still need to refer to whatever efficacy is in the PassiveSevere PhaseSeverity. Accounted for using pointers already
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			Calc_SumRhoEffs_c_BS_PS(PARAMS, HOUSE, country, BaselineSeroStatus, PhaseSeverity);
}
void Calc_SumRhoEffs_BaselineSeroStatus					(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int BaselineSeroStatus) //// must call this (or related) when changing any efficacy (inc. when not serospecific) or rho parameter
{
	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++) 	//// Note: NOT  <HOUSE.NumEffsPer_BS_And_Serotype. Still need to refer to whatever efficacy is in the PassiveSevere PhaseSeverity. Accounted for using pointers already
		for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++)
			Calc_SumRhoEffs_c_BS_PS(PARAMS, HOUSE, HOUSE.WhichCountries[countryindex], BaselineSeroStatus, PhaseSeverity);
}
void Calc_SumRhoEffs_PhaseSeverity						(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int PhaseSeverity) //// must call this (or related) when changing any efficacy (inc. when not serospecific) or rho parameter
{
	for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++)
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			Calc_SumRhoEffs_c_BS_PS(PARAMS, HOUSE, HOUSE.WhichCountries[countryindex], BaselineSeroStatus, PhaseSeverity);
}
void Calc_SumRhoEffs_country_BaselineSeroStatus			(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int country, int BaselineSeroStatus)
{
	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)	//// Note: NOT  <HOUSE.NumEffsPer_BS_And_Serotype. Still need to refer to whatever efficacy is in the PassiveSevere PhaseSeverity. Accounted for using pointers already
		Calc_SumRhoEffs_c_BS_PS(PARAMS, HOUSE, country, BaselineSeroStatus, PhaseSeverity);
}
void Calc_SumRhoEffs_country_PhaseSeverity				(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int country, int PhaseSeverity)
{
	for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
		Calc_SumRhoEffs_c_BS_PS(PARAMS, HOUSE, country, BaselineSeroStatus, PhaseSeverity);
}
void Calc_SumRhoEffs_BaselineSeroStatus_PhaseSeverity	(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int BaselineSeroStatus, int PhaseSeverity)
{
	for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++)
		Calc_SumRhoEffs_c_BS_PS(PARAMS, HOUSE, HOUSE.WhichCountries[countryindex], BaselineSeroStatus, PhaseSeverity);
}
void Calc_SumRhoEffs_ALL								(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++)
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++) 	//// Note: NOT  <HOUSE.NumEffsPer_BS_And_Serotype. Still need to refer to whatever efficacy is in the PassiveSevere PhaseSeverity. Accounted for using pointers already
				Calc_SumRhoEffs_c_BS_PS(PARAMS, HOUSE, HOUSE.WhichCountries[countryindex], BaselineSeroStatus, PhaseSeverity);
}
void Calc_SumRhoEffs_ALL								(Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE)
{
	Calc_SumRhoEffs_ALL	(CurrentPARAMS	, HOUSE);		/// must be done before the threading. 
	Calc_SumRhoEffs_ALL	(ProposedPARAMS , HOUSE);		/// must be done before the threading. 
}
void Calc_SumRhoKs_country_PhaseSeverity		(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int &country, int PhaseSeverity)
{
	for (int PrevInf = 0; PrevInf < HOUSE.Num_K_Params; PrevInf++) 
		Calc_SumRhoK_country_PhaseSeverity_PrevInf(PARAMS, HOUSE, country, PhaseSeverity, PrevInf);
}
void Calc_SumRhoKs_PhaseSeverity_PrevInf		(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int PhaseSeverity, int PrevInf)
{
	for (int country = 0; country < HOUSE.TotalCountries; country++)
		Calc_SumRhoK_country_PhaseSeverity_PrevInf(PARAMS, HOUSE, country, PhaseSeverity, PrevInf);
}
void Calc_SumRhoKs_country_PrevInf				(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int &country, int PrevInf)
{
	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
		Calc_SumRhoK_country_PhaseSeverity_PrevInf(PARAMS, HOUSE, country, PhaseSeverity, PrevInf);
}
void Calc_SumRhoKs_country						(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int &country)
{
	for (int PrevInf = 0; PrevInf < HOUSE.Num_K_Params; PrevInf++)
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			Calc_SumRhoK_country_PhaseSeverity_PrevInf(PARAMS, HOUSE, country, PhaseSeverity, PrevInf);
}
void Calc_SumRhoKs_PhaseSeverity				(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int PhaseSeverity)
{
	for (int country = 0; country < HOUSE.TotalCountries; country++)
		for (int PrevInf = 0; PrevInf < HOUSE.Num_K_Params; PrevInf++)
			Calc_SumRhoK_country_PhaseSeverity_PrevInf(PARAMS, HOUSE, country, PhaseSeverity, PrevInf);
}
void Calc_SumRhoKs_PrevInf						(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int PrevInf)
{
	for (int country = 0; country < HOUSE.TotalCountries; country++)
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			Calc_SumRhoK_country_PhaseSeverity_PrevInf(PARAMS, HOUSE, country, PhaseSeverity, PrevInf);
}
void Calc_SumRhoKs_ALL							(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	for (int country = 0; country < HOUSE.TotalCountries; country++)
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int PrevInf = 0; PrevInf < HOUSE.Num_K_Params; PrevInf++)
				Calc_SumRhoK_country_PhaseSeverity_PrevInf(PARAMS, HOUSE, country, PhaseSeverity, PrevInf);
}
void Calc_SumRhoKs_ALL							(Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE)
{
	Calc_SumRhoKs_ALL(CurrentPARAMS , HOUSE);		/// must be done before the threading. 
	Calc_SumRhoKs_ALL(ProposedPARAMS, HOUSE);		/// must be done before the threading. 
}
void Calc_SumRhoEffKs_country_PhaseSeverity				(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int &country, int PhaseSeverity)
{
	for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
		Calc_SumRhoEffKs_c_BS_PS(PARAMS, HOUSE, country, PhaseSeverity, BaselineSeroStatus);
}
void Calc_SumRhoEffKs_PhaseSeverity_BaselineSeroStatus	(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int PhaseSeverity, int BaselineSeroStatus)
{
	for (int country = 0; country < HOUSE.TotalCountries; country++)
		Calc_SumRhoEffKs_c_BS_PS(PARAMS, HOUSE, country, PhaseSeverity, BaselineSeroStatus);
}
void Calc_SumRhoEffKs_country_BaselineSeroStatus		(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int country, int BaselineSeroStatus)
{
	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
		Calc_SumRhoEffKs_c_BS_PS(PARAMS, HOUSE, country, PhaseSeverity, BaselineSeroStatus);
}
void Calc_SumRhoEffKs_country							(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int &country)
{
	for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			Calc_SumRhoEffKs_c_BS_PS(PARAMS, HOUSE, country, PhaseSeverity, BaselineSeroStatus);
}
void Calc_SumRhoEffKs_PhaseSeverity						(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int PhaseSeverity)
{
	for (int country = 0; country < HOUSE.TotalCountries; country++)
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			Calc_SumRhoEffKs_c_BS_PS(PARAMS, HOUSE, country, PhaseSeverity, BaselineSeroStatus);
}
void Calc_SumRhoEffKs_BaselineSeroStatus				(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int BaselineSeroStatus)
{
	for (int country = 0; country < HOUSE.TotalCountries; country++)
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			Calc_SumRhoEffKs_c_BS_PS(PARAMS, HOUSE, country, PhaseSeverity, BaselineSeroStatus);
}
void Calc_SumRhoEffKs_ALL								(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	for (int country = 0; country < HOUSE.TotalCountries; country++)
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
				Calc_SumRhoEffKs_c_BS_PS(PARAMS, HOUSE, country, PhaseSeverity, BaselineSeroStatus);
}
void Calc_SumRhoEffKs_ALL								(Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE)
{
	Calc_SumRhoEffKs_ALL(CurrentPARAMS , HOUSE);		/// must be done before the threading. 
	Calc_SumRhoEffKs_ALL(ProposedPARAMS, HOUSE);		/// must be done before the threading. 
}
void Calc_Hospital_Ks_country		(int country, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
		for (int PrevInf = 0; PrevInf < HOUSE.Num_K_Params; PrevInf++)
			PARAMS.K_s[country][PassiveSevere][PrevInf][serotype] = PARAMS.K_s[0][PassiveSevere][PrevInf][serotype] * PARAMS.Hosp_K_mult; //// i.e. multiply regular K by Hosp_K_mult. Choice of country 0 is arbitrary - just needs to be not in CYD-15 countries (or whichever countries you're modelling different hospital K's for). 
}
void Calc_Hospital_Ks				(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	//// See note in Update_Severe_Ks function as to how this function differs from that. 
	for (int countryindex = 0; countryindex < HOUSE.CYD_15_countries.size(); countryindex++) //// do only CYD15 countries. 
		for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
			for (int PrevInf = 0; PrevInf < HOUSE.Num_K_Params; PrevInf++)
				PARAMS.K_s[HOUSE.CYD_15_countries[countryindex]][PassiveSevere][PrevInf][serotype] = PARAMS.K_s[0][PassiveSevere][PrevInf][serotype] * PARAMS.Hosp_K_mult; //// i.e. multiply regular K by Hosp_K_mult. Choice of country 0 is arbitrary - just needs to be not in CYD-15 countries (or whichever countries modelling different hospital K's for). 
}
void Update_Severe_Ks					(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	//// used for HOUSE.Fixed_Severe_RelRisks, not the same as Calc_Hospital_Ks function when HOUSE.ModelHosp_Indie_Ks == false (default is true). This latter function uses a fitted multiplier of hospital Ks, rather then completely independent PassiveSevere K parameters for CYD15 trial. 
	if (HOUSE.Fixed_Severe_RelRisks)
		for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++) ////  do all countries being fitted.
			for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
				for (int PrevInf = 0; PrevInf < HOUSE.Num_K_Params; PrevInf++)
					PARAMS.K_s[HOUSE.WhichCountries[countryindex]][PassiveSevere][PrevInf][serotype] =
					PARAMS.K_s[HOUSE.WhichCountries[countryindex]][PassiveSevere][1		 ][serotype] * PARAMS.Fixed_SevereK_ratios[PrevInf]; //// i.e. multiply KH_0 and KH_2 (and also KH_1) by their fixed ratios (the ratio for KH1 set to 1).  
}
void Calc_SumRhoEffNegs_country				(int country, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			Calc_SumRhoEffNegs_c_BS_PS(country, BaselineSeroStatus, PhaseSeverity, PARAMS, HOUSE);
}
void Calc_SumRhoEffNegs_PhaseSeverity		(int PhaseSeverity, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
		for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++)
			Calc_SumRhoEffNegs_c_BS_PS(HOUSE.WhichCountries[countryindex], BaselineSeroStatus, PhaseSeverity, PARAMS, HOUSE);
}
void Calc_SumRhoEffNegs_BaselineSeroStatus	(int BaselineSeroStatus, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
		for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++)
			Calc_SumRhoEffNegs_c_BS_PS(HOUSE.WhichCountries[countryindex], BaselineSeroStatus, PhaseSeverity, PARAMS, HOUSE);
}
void Calc_SumRhoEffNegs_All					(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	for (int countryindex = 0; countryindex < HOUSE.WhichCountries.size(); countryindex++)
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
				Calc_SumRhoEffNegs_c_BS_PS(HOUSE.WhichCountries[countryindex], BaselineSeroStatus, PhaseSeverity, PARAMS, HOUSE);
}
///// END wrapper functions

DType ExponentialWaning					(DType WaningValue, DType daypostdose, const Housekeeping_Struct &HOUSE)
{
	return(exp(-(HOUSE.TimeInterval * daypostdose) / WaningValue));
}
void Calc_ExpWaningMults					(DType *WaningArray, const Housekeeping_Struct &HOUSE, const int &NoDaysOfFollowUp, DType WaningValue)
{
#pragma omp parallel for schedule(static,1) 
		for (int thread_no = 0; thread_no < HOUSE.max_threads; thread_no++)
			for (int daypostdose = thread_no; daypostdose < NoDaysOfFollowUp; daypostdose += HOUSE.max_threads)
				WaningArray[daypostdose] = exp(-(HOUSE.TimeInterval * (DType)daypostdose) / WaningValue);
}
void Calc_WaningValues_BS_Age		(int BaselineSeroStatus, int AgeInYears, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, const int &NoDaysOfFollowUp)
{
	//// This function calcultes the Waning multipliers over applied to vaccine efficacy (age-specific or not) over time. 
	//// Function only called if assuming exponential waning. 
	//// Function calculates the waning value, which is always a duration (even if we fit the rate), and then used this in the exponential function #pragma loop. 
	//// Function calculates duration based on value of Age_Option. In short, WaningValue is either simply a parameter, or the duration is calculated by other parameters. 

	DType WaningValue; 
	if (HOUSE.AS_Waning == Age_Option::INDEPENDENT)	
		WaningValue = PARAMS.WaningParams[BaselineSeroStatus][WaningRateDuration];																																		//// tau_b
	else if (HOUSE.AS_Waning_OnlyOneSeroStatus && HOUSE.AS_Waning_OneSeroBS != BaselineSeroStatus) //// true regardless of how age specific waning is modelled
		WaningValue = PARAMS.WaningParams[BaselineSeroStatus][WaningRateDuration]; 
	else if (HOUSE.AS_Waning == Age_Option::HILL) //// i.e. if using a hill function to model the age differences in waning rate/duration
		WaningValue = PARAMS.WaningParams[BaselineSeroStatus][WaningRateDuration] * Hill(AgeInYears, PARAMS.WaningParams[BaselineSeroStatus][Power_Index], PARAMS.WaningParams[BaselineSeroStatus][Halflife_Index]);	//// tau_b * Hill(age, Power, Halflife)
	else if (HOUSE.AS_Waning == Age_Option::CATEGORICAL)
	{
			 if (AgeInYears <= 5)	WaningValue = PARAMS.WaningParams[BaselineSeroStatus][0];
		else if (AgeInYears <= 11)	WaningValue = PARAMS.WaningParams[BaselineSeroStatus][1];
		else						WaningValue = PARAMS.WaningParams[BaselineSeroStatus][2];
	}
	else if (HOUSE.AS_Waning == Age_Option::SPLINE || HOUSE.AS_Waning == Age_Option::SPLINE_LINE || HOUSE.AS_Waning == Age_Option::SPLINE_STEP || HOUSE.AS_Waning == Age_Option::CUBIC)
		WaningValue = WaningDurationSpline((DType)AgeInYears, BaselineSeroStatus, PARAMS, HOUSE);
	else	std::cerr << "CalcWaningValues_BS_Age ERROR: if/else tree wrong" << endl;

#pragma omp parallel for schedule(static,1) 
	for (int thread_no = 0; thread_no < HOUSE.max_threads; thread_no++)
		for (int daypostdose = thread_no; daypostdose < NoDaysOfFollowUp; daypostdose += HOUSE.max_threads)
			PARAMS.WaningMults[AgeInYears][BaselineSeroStatus][daypostdose] = exp(-(HOUSE.TimeInterval * (DType)daypostdose) / WaningValue);
}
void Calc_WaningValues				(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int BaselineSeroStatus, const int &NoDaysOfFollowUp)
{
	// ** // Wrapper function that calculates waning values for a particular serostatus, and all ages if doing HOUSE.AS_Waning. 
	if (!HOUSE.HillWaning)
		for (int AgeInYears = 0; AgeInYears < HOUSE.NumAges_Waning; AgeInYears++)
			Calc_WaningValues_BS_Age(BaselineSeroStatus, AgeInYears, PARAMS, HOUSE, NoDaysOfFollowUp);
	else
#pragma omp parallel for schedule(static,1) 
		for (int thread_no = 0; thread_no < HOUSE.max_threads; thread_no++)
			for (int daypostdose = thread_no; daypostdose < NoDaysOfFollowUp; daypostdose += HOUSE.max_threads)
				PARAMS.WaningMults[0][BaselineSeroStatus][daypostdose] = 1 - Hill(HOUSE.TimeInterval * DType(daypostdose), PARAMS.Hill_Powers[BaselineSeroStatus], PARAMS.Hill_Halflives[BaselineSeroStatus]);
}
void Calc_AgeVaccEffMults			(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int BaselineSeroStatus)
{
	if (HOUSE.ASVE == Age_Option::HILL)
	{
		for (int age = 0; age < HOUSE.HowManyAges; age++)
			PARAMS.AgeEff_Mult[BaselineSeroStatus][age] =
			(1 - PARAMS.ASVE_Params[BaselineSeroStatus][Prop_Index]) + PARAMS.ASVE_Params[BaselineSeroStatus][Prop_Index] * Hill(DType(age), PARAMS.ASVE_Params[BaselineSeroStatus][Power_Index], PARAMS.ASVE_Params[BaselineSeroStatus][Halflife_Index]);
	}
	else	if (HOUSE.ASVE == Age_Option::CATEGORICAL) //// will work only for CYD-14 only. 
	{
		for (int age = 0; age < HOUSE.HowManyAges; age++)
		{
				 if (age <= 5)	PARAMS.AgeEff_Mult[BaselineSeroStatus][age] = PARAMS.ASVE_Params[BaselineSeroStatus][1];
			else if (age <= 11) PARAMS.AgeEff_Mult[BaselineSeroStatus][age] = PARAMS.ASVE_Params[BaselineSeroStatus][2];
			else				PARAMS.AgeEff_Mult[BaselineSeroStatus][age] = PARAMS.ASVE_Params[BaselineSeroStatus][3];
		}
	}
	else	//// do this for all splines
		for (int age = 0; age < HOUSE.HowManyAges; age++)
			PARAMS.AgeEff_Mult[BaselineSeroStatus][age] = EfficacyMultiplierSpline(DType(age), BaselineSeroStatus, PARAMS, HOUSE);
}
void Calc_Age_HazMults				(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	if (HOUSE.AS_Haz == Age_Option::HILL)
		for (int age = 0; age < HOUSE.HowManyAges; age++)
			PARAMS.AgeHaz_Mult[age] =
			(1 - PARAMS.ASHaz_Params[Prop_Index]) + PARAMS.ASHaz_Params[Prop_Index] * Hill(DType(age), PARAMS.ASHaz_Params[Power_Index], PARAMS.ASHaz_Params[Halflife_Index]);
	else if (HOUSE.AS_Haz == Age_Option::CATEGORICAL)
		for (int age = 0; age < HOUSE.HowManyAges; age++)
		{
				 if (age <= 5)	PARAMS.AgeHaz_Mult[age] = PARAMS.ASHaz_Params[1];  
			else if (age <= 11) PARAMS.AgeHaz_Mult[age] = PARAMS.ASHaz_Params[2];
			else				PARAMS.AgeHaz_Mult[age] = PARAMS.ASHaz_Params[3];
		}
	else	//// do this for all splines
		for (int age = 0; age < HOUSE.HowManyAges; age++)
			PARAMS.AgeHaz_Mult[age] = AS_HazMultSpline(age, PARAMS, HOUSE); 
}

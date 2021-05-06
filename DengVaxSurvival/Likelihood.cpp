
#include "HeaderAndSwitches.h"
#include "Probability.h"
#include "Splines.h"
#include "ConsolePrinting.h"
#include "Augmentation.h"
#include "ParamUpdate.h"
#include "CalcParamsEtc.h"


//// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// 
//// CASE LIKELIHOOD FUNCTIONS - i.e. parts of instantaneous hazard

DType l_K_case			(const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, const int &country, char PhaseSeverity, int PrevInf, int &serotype)
{
#ifdef PRINT_WITHIN_LIKE_FUNCTIONS
	std::cout << "l_K_case: c " << country << " K_" << (int)PhaseSeverity << "_" << PrevInf << " serotype" << (int)serotype << endl;
#endif

	if (HOUSE.ModelVariant == AS_PRIME && PrevInf != 0)	
		std::cerr << "l_K_case error: HOUSE.ModelVariant == AS_PRIME && PrevInf != 0" << endl; 

	if ((PhaseSeverity == ActiveMild) && (PrevInf == 1))  //// baseline 
	{
		return 0;
	}
	else
	{
		int Set_Index = NULL; 
			 if (PrevInf == 0)	Set_Index = (HOUSE.ModelVariant == VAC_SILENT || HOUSE.ModelVariant == AS_PRIME) ? HOUSE.Strata[ControlGroup][SeroNeg] : HOUSE.Strata[EitherTrialArm][SeroNeg];
		else if (PrevInf == 1)	Set_Index = (HOUSE.ModelVariant == VAC_SILENT) ? HOUSE.Strata[VaccineGroup][SeroNeg] : HOUSE.Strata[EitherTrialArm][SeroPos];	//// should not call for K_SEROPOS
		else if (PrevInf == 2)	Set_Index = HOUSE.Strata[VaccineGroup][SeroPos];																				//// only have K2 for VAC_SILENT. 

		int NoCasesInHazGroup = 0;		//// will change between augmentation steps. 
#pragma omp parallel for schedule(static,1) reduction(+:NoCasesInHazGroup)
		for (int thread_no = 0; thread_no < HOUSE.max_threads; thread_no++)
			for (int i = thread_no; i < DATA.Cases_Array[country][Set_Index][PhaseSeverity].size(); i += HOUSE.max_threads)
				if ((!HOUSE.SeroSpecific_K_values) || (DATA.CaseSerotype[DATA.Cases_Array[country][Set_Index][PhaseSeverity][i]] == serotype)) //// i.e. if not doing SeroSpecific_K_values, only need to consider if patient was a case using if statment from line above. If doing SeroSpecific_K_values, serotype also needs to match. Tempting to simply record all non sero-specific cases as serotype 1 (i.e. 0 in Cpp), but this will ruin everything if you do SS_VEs but not SS_Ks, so don't be tempted. 
					NoCasesInHazGroup++; 

		DType l_K_case = NoCasesInHazGroup * log(PARAMS.K_s[country][PhaseSeverity][PrevInf][serotype]);

		return l_K_case;
	}
}
DType l_K0SNegPrime_case(const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, const int &country, char PhaseSeverity, int &serotype)
{
#ifdef PRINT_WITHIN_LIKE_FUNCTIONS
	std::cout << "l_K0SNegPrime_case: c " << country << " PS" << (int)PhaseSeverity << " serotype" << (int)serotype << endl;
#endif

	if (HOUSE.ModelVariant != AS_PRIME) std::cerr << "l_KPlusPrime_case error: HOUSE.ModelVariant != AS_PRIME" << endl; 
	DType l_10 = 0;
	int Set_Index = HOUSE.Strata[VaccineGroup][SeroNeg];
	
#pragma omp parallel for schedule(static,1) reduction(+:l_10)
	for (int thread_no = 0; thread_no < HOUSE.max_threads; thread_no++)
		for (int i = thread_no; i < DATA.Cases_Array[country][Set_Index][PhaseSeverity].size(); i += HOUSE.max_threads) 
			if ((!HOUSE.SeroSpecific_K_values) || (DATA.CaseSerotype[DATA.Cases_Array[country][Set_Index][PhaseSeverity][i]] == serotype))
				l_10 += log(PARAMS.Ks_Prime[country][PhaseSeverity][SeroNeg][serotype][DATA.ai_s[DATA.Cases_Array[country][Set_Index][PhaseSeverity][i]]]); 	
	return l_10;
}
DType l_Kplus_case		(const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, const int &country, char PhaseSeverity, int &serotype) 
{
#ifdef PRINT_WITHIN_LIKE_FUNCTIONS
	std::cout << "l_Kplus_case: c " << country << " PS" << (int)PhaseSeverity << " serotype" << (int)serotype << endl;
#endif
	if (HOUSE.ModelVariant == SIMPLE_NUMERICAL) std::cerr << "l_Kplus_case error: HOUSE.ModelVariant == SIMPLE_NUMERICAL" << endl;

	DType l_10 = 0;
	int Set_Index = (HOUSE.ModelVariant == VAC_SILENT || HOUSE.ModelVariant == AS_PRIME) ? HOUSE.Strata[ControlGroup][SeroPos] : HOUSE.Strata[EitherTrialArm][SeroPos]; //// In VAC_SILENT & AS_PRIME, the multiplier of hazard for seropositives is not the same for vaccine and control arms. In all other variants hazard multiplier the same across arms. 

#pragma omp parallel for schedule(static,1) reduction(+:l_10)
	for (int thread_no = 0; thread_no < HOUSE.max_threads; thread_no++)
		for (int i = thread_no; i < DATA.Cases_Array[country][Set_Index][PhaseSeverity].size(); i += HOUSE.max_threads) 
			if ((!HOUSE.SeroSpecific_K_values) || (DATA.CaseSerotype[DATA.Cases_Array[country][Set_Index][PhaseSeverity][i]] == serotype))
				l_10 += log(	K_seropos(		PARAMS.ParamVec[PARAMS.ParamNumbers.Min_HHaz + country]				,	//// hci
												DATA.ai_s[DATA.Cases_Array[country][Set_Index][PhaseSeverity][i]]	,	//// age
												PARAMS.K_s[country][PhaseSeverity][1][serotype]						,	//// K1
												PARAMS.K_s[country][PhaseSeverity][2][serotype]						)	//// K2
					);
	return l_10;
}
DType l_KPlusPrime_case	(const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, const int &country, char PhaseSeverity, int &serotype)  
{
#ifdef PRINT_WITHIN_LIKE_FUNCTIONS
	std::cout << "l_KPlusPrime_case: c " << country << " PS" << (int)PhaseSeverity << " serotype" << (int)serotype << endl;
#endif

	//// Only difference between this and l_Kplus_case is that you use Ks_Prime rather than K_s for K1
	if (HOUSE.ModelVariant != AS_PRIME) std::cerr << "l_KPlusPrime_case error: HOUSE.ModelVariant != AS_PRIME" << endl; //// take this out when all working nicely. 

	DType l_10 = 0;
	int Set_Index = HOUSE.Strata[VaccineGroup][SeroPos];

#pragma omp parallel for schedule(static,1) reduction(+:l_10)
	for (int thread_no = 0; thread_no < HOUSE.max_threads; thread_no++)
		for (int i = thread_no; i < DATA.Cases_Array[country][Set_Index][PhaseSeverity].size(); i += HOUSE.max_threads)
			if ((!HOUSE.SeroSpecific_K_values) || (DATA.CaseSerotype[DATA.Cases_Array[country][Set_Index][PhaseSeverity][i]] == serotype))
				l_10 += log(	K_seropos(		PARAMS.ParamVec[PARAMS.ParamNumbers.Min_HHaz + country]	,	//// hci
												DATA.ai_s[DATA.Cases_Array[country][Set_Index][PhaseSeverity][i]],	//// age
												PARAMS.Ks_Prime[country][PhaseSeverity][SeroPos][serotype][DATA.ai_s[DATA.Cases_Array[country][Set_Index][PhaseSeverity][i]]],	//// K1, +, prime
												PARAMS.K_s[country][PhaseSeverity][2][serotype]			)	//// K2, +, prime = K2
					);
	return l_10;
}
DType l_BaseHaz			(const int &country, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, char PhaseSeverity)
{
#ifdef PRINT_WITHIN_LIKE_FUNCTIONS
	std::cout << "l_BaseHaz: c " << country << " PS" << (int)PhaseSeverity << endl;
#endif
	DType l_1 = 0;
	int Set_Index = HOUSE.Strata[EitherTrialArm][EitherSeroStatus];

#pragma omp parallel for schedule(static,1) reduction(+:l_1)
	for (int thread_no = 0; thread_no < HOUSE.max_threads; thread_no++)
		for (int i = thread_no; i < DATA.Cases_Array[country][Set_Index][PhaseSeverity].size(); i += HOUSE.max_threads)
			l_1 += log(PARAMS.BaselineHazardValues[country][DATA.FollowUp[END][PassiveSevere][DATA.Cases_Array[country][Set_Index][PhaseSeverity][i]]]);
	return l_1;
}
DType WaningEfficacy	(char BaselineSeroStatus, int &serotype, char PhaseSeverity, int AgeInYears, DType WaningMult, EffNegWane_Option ENW_Option, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, bool NaturalLog)
{
	DType WaningEffBaseHazMult	= 1; //// (natural log of) multiplier to (not-integrated) baseline hazard, which accounts for efficacy (inc BaselineSeroStatus, serotype and age) and waning. 

	/*	NOTE: You've taken out InfEffs stuff from here as you never used it.	*/
	/*
		if TI_bd(age) >= 0			OR		if (HOUSE.EffNegWane == EffNegWane_Option::DEFAULT)

			want (1 -  (TI_bd(age) x Waning(t_final))		
	
		if TI_bd(age) < 0

			if (HOUSE.EffNegWane == EffNegWane_Option::FROM_ZERO)	, want (1 -  ( abs(TI_bd(age)) x (Waning(t) - 1)	) 
			if (HOUSE.EffNegWane == EffNegWane_Option::NO_WANE	)	, want (1 -  ( TI_bd(age) )

		If calculating for Log likelihoods, will take natural log of above quantities. If calculating Total_Likelihoods (or gibbs augmentation), will not need natural logs of above quantities. 
	*/

	//// get transient immunity for this baseline serostatus (b), dengue serotype (d), and age in years. 
	DType TI_bd_age					= Calc_AggTransImmunity(BaselineSeroStatus, serotype, AgeInYears, PhaseSeverity, PARAMS, HOUSE); //// Haven't used overload with ENW_Option as need to know what to add to waning mult (0 if AggEff >=0 or !HOUSE.AltWaneEffNeg, 1 otherwise). 

		 if (TI_bd_age >= 0 || HOUSE.EffNegWane == EffNegWane_Option::DEFAULT)		WaningEffBaseHazMult = 	1 - (	TI_bd_age	*	WaningMult		)	; 
	else if (HOUSE.EffNegWane == EffNegWane_Option::FROM_ZERO)						WaningEffBaseHazMult = 	1 - (abs(TI_bd_age) * (WaningMult - 1)	)	;
	else if (HOUSE.EffNegWane == EffNegWane_Option::NO_WANE)						WaningEffBaseHazMult = 	1 -				TI_bd_age					;

	if (NaturalLog) WaningEffBaseHazMult = log(WaningEffBaseHazMult); 

	return WaningEffBaseHazMult;
}
DType l_WaningEfficacy	(char BaselineSeroStatus, const int &country, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int &serotype, char PhaseSeverity)
{
#ifdef PRINT_WITHIN_LIKE_FUNCTIONS
	std::cout << "l_WaningEfficacy: c " << country << " BS" << (int)BaselineSeroStatus << " PS" << (int)PhaseSeverity << " serotype" << (int)serotype << endl;
#endif
	DType l_c7or9 = 0;
	int Set_Index = HOUSE.Strata[VaccineGroup][BaselineSeroStatus];

#pragma omp parallel for schedule(static,1) reduction(+:l_c7or9)
	for (int thread_no = 0; thread_no < HOUSE.max_threads; thread_no++)
		for (int i = thread_no; i < DATA.Cases_Array[country][Set_Index][PhaseSeverity].size(); i += HOUSE.max_threads)
			if ((!HOUSE.SeroSpecificEfficacies) || (DATA.CaseSerotype[DATA.Cases_Array[country][Set_Index][PhaseSeverity][i]] == serotype))
			{
				int person_i			= DATA.Cases_Array[country][Set_Index][PhaseSeverity][i];											
				int AgeInYears_person_i = DATA.ai_s[person_i];																				
				DType WaningMult		= PARAMS.WaningMults[AgeInYears_person_i][BaselineSeroStatus][DATA.TimePost_Final_Dose[person_i]];	

				l_c7or9 += WaningEfficacy(BaselineSeroStatus, serotype, PhaseSeverity, AgeInYears_person_i, 
					PARAMS.WaningMults[AgeInYears_person_i][BaselineSeroStatus][DATA.TimePost_Final_Dose[person_i]], HOUSE.EffNegWane, PARAMS, HOUSE, /*NaturalLog =*/ true);
			}
	return l_c7or9;
}
DType l_rho_case		(const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, const int &country, int serotype, char PhaseSeverity)
{
#ifdef PRINT_WITHIN_LIKE_FUNCTIONS
	std::cout << "l_WaningEfficacy: c " << country << " PS" << (int)PhaseSeverity << " serotype" << (int)serotype << endl;
#endif
	int Set_Index = HOUSE.Strata[EitherTrialArm][EitherSeroStatus];

	int NoCasesInHazGroup = 0; //// will change between augmentation steps. 
#pragma omp parallel for schedule(static,1) reduction(+:NoCasesInHazGroup)
	for (int thread_no = 0; thread_no < HOUSE.max_threads; thread_no++)
		for (int i = thread_no; i < DATA.Cases_Array[country][Set_Index][PhaseSeverity].size(); i += HOUSE.max_threads)
			if (DATA.CaseSerotype[DATA.Cases_Array[country][Set_Index][PhaseSeverity][i]] == serotype) NoCasesInHazGroup++;

	DType l_log_rho_cd	= NoCasesInHazGroup * log(PARAMS.rhos[country][serotype]);

	return l_log_rho_cd;
}
DType l_AS_HazMult		(const int &country, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, char PhaseSeverity)
{
#ifdef PRINT_WITHIN_LIKE_FUNCTIONS
	std::cout << "l_AS_HazMult: c " << country <<  " PS" << (int)PhaseSeverity <<  endl;
#endif
	DType l_AS_HazMult = 0;
	int Set_Index = HOUSE.Strata[EitherTrialArm][EitherSeroStatus];

#pragma omp parallel for schedule(static,1) reduction(+:l_AS_HazMult)
	for (int thread_no = 0; thread_no < HOUSE.max_threads; thread_no++)
		for (int i = thread_no; i < DATA.Cases_Array[country][Set_Index][PhaseSeverity].size(); i += HOUSE.max_threads)
			l_AS_HazMult += log(PARAMS.AgeHaz_Mult[DATA.ai_s[DATA.Cases_Array[country][Set_Index][PhaseSeverity][i]]]);
	return l_AS_HazMult;
}
DType l_BS_BaseHaz_Mult	(const int &country, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, char PhaseSeverity, char BaselineSeroStatus)
{
#ifdef PRINT_WITHIN_LIKE_FUNCTIONS
	std::cout << "l_WaningEfficacy: c " << country << " BS" << (int)BaselineSeroStatus << " PS" << (int)PhaseSeverity <<  endl;
#endif
	DType l_BS_BaseHaz_Mult = 0;

	if (HOUSE.Which_BS_BaseHazMult == BaselineSeroStatus)
	{
		int Set_Index = HOUSE.Strata[EitherTrialArm][BaselineSeroStatus];

		int NoCasesInHazGroup = 0; //// will change between augmentation steps. 
#pragma omp parallel for schedule(static,1) reduction(+:NoCasesInHazGroup)
		for (int thread_no = 0; thread_no < HOUSE.max_threads; thread_no++)
			for (int i = thread_no; i < DATA.Cases_Array[country][Set_Index][PhaseSeverity].size(); i += HOUSE.max_threads) NoCasesInHazGroup++; 

		l_BS_BaseHaz_Mult = NoCasesInHazGroup * log(PARAMS.BS_BaseHazMults[BaselineSeroStatus]); 

	} /// else do nothing (i.e. leave l_BS_BaseHaz_Mult equal to 0;)
	return l_BS_BaseHaz_Mult;
}


//// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// 
//// AUGMENTATION LIKELIHOOD FUNCTION - used for both cases and non-cases.

DType l_Aug				(const int &country, int BaselineSeroStatus, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, bool ImSubOnly)
{
#ifdef PRINT_WITHIN_LIKE_FUNCTIONS
	std::cout << "l_Aug: c " << country << " BS " << BaselineSeroStatus << endl;
#endif
	DType Sum_LL_BS = 0, hci = PARAMS.ParamVec[PARAMS.ParamNumbers.Min_HHaz + country];;
	int Set_Index = HOUSE.Strata[EitherTrialArm][BaselineSeroStatus];
#pragma omp parallel for schedule(static,1) reduction(+:Sum_LL_BS)
	for (int thread_no = 0; thread_no < HOUSE.max_threads; thread_no++)
		for (int i = thread_no; i < DATA.Set_Array[country][Set_Index].size(); i += HOUSE.max_threads)
			if (!ImSubOnly | DATA.InImmSub[DATA.Set_Array[country][Set_Index][i]])
				Sum_LL_BS += PARAMS.SeroPrevs[LogIndex][country][BaselineSeroStatus][DATA.ai_s[DATA.Set_Array[country][Set_Index][i]]];
	return Sum_LL_BS;
}
DType l_Aug				(const int &country, int BaselineSeroStatus, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	//// overload of l_Aug function above used to apply default: ImSubOnly == false. This default will be called most often
	DType l_BS_Default = l_Aug(country, BaselineSeroStatus, DATA, PARAMS, HOUSE, false);
	return l_BS_Default; 
}

//// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// 
//// SURVIVOR LIKELIHOOD FUNCTIONS - i.e. integrated hazards needed for cases and non-cases. 

DType l_Int_Base_Haz			(const int &country, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, char PhaseSeverity, char BaselineSeroStatus, char TrialArm) 
{
#ifdef PRINT_WITHIN_LIKE_FUNCTIONS
	std::cout << "l_Int_Base_Haz: Arm " << (int)TrialArm << " Kmult" << Kmult << " BS" << int(BaselineSeroStatus) << " c" << country << " PS" << (int)PhaseSeverity << endl;
#endif
	DType LikeComponent = 0; 
	int Set_Index = HOUSE.Strata[TrialArm][BaselineSeroStatus];

#pragma omp parallel for schedule(static,1) reduction(+:LikeComponent)
	for (int thread_no = 0; thread_no < HOUSE.max_threads; thread_no++) 
	{
		DType AmountToSubtract;
		int person_i; 
		for (int i = thread_no; i < DATA.Set_Array[country][Set_Index].size(); i += HOUSE.max_threads)
		{	
			person_i			= DATA.Set_Array[country][Set_Index][i]; 
			AmountToSubtract	= IntBaseHaz(DATA.FollowUp[START][PhaseSeverity][person_i], DATA.FollowUp[END][PhaseSeverity][person_i], country, PARAMS);
			AmountToSubtract	*= Choose_SumBaseHazMult(TrialArm, BaselineSeroStatus, country, DATA.ai_s[person_i], PhaseSeverity, PARAMS, HOUSE);

			LikeComponent		-= AmountToSubtract;
		}
	}
	return LikeComponent;
}
DType l_Int_Vac_Haz				(char BaselineSeroStatus, const int &country, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, char PhaseSeverity)	 
{
	DType LikeComponent = 0;
	int Set_Index = HOUSE.Strata[VaccineGroup][BaselineSeroStatus];

	//// Choose Kmult: i.e. will equal 1 if cannot take out K as common factor (e.g. if K_Seropos variant or if SSKs and SSVEs together). 
	//// The below applies to: i) no serotype effects; ii) SS_Ks only; iii) SS_VEs only; iv) SSASVEs default i.e. multiplicative. 

	DType Kmult = 1; /// Used for K_SEROPOS, where K+ values required in #pragma loop; and for SS_Ks and SS_VEs together, where each serotype efficacy has a corresponding serotype RR, i.e. not a common factor K. Note that if not SSKs, SumRhoKs are simply Ks, hence rhos not taken into account. That comes in l_Int_Vac_Haz_sero_i. 
	if (BaselineSeroStatus == SeroPos)
	{
			 if (HOUSE.ModelVariant == VAC_SILENT		)										Kmult = PARAMS.SumRhoKs[country][PhaseSeverity][2];
		else if (HOUSE.ModelVariant == SIMPLE_NUMERICAL	)										Kmult = PARAMS.SumRhoKs[country][PhaseSeverity][1]; //// no need for K_SEROPOS as Kmult already equals 1 and you have KplusValues in loop
		// else don't do anything: K_SEROPOS and AS_PRIME don't have common factor
	}
	else	if (BaselineSeroStatus == SeroNeg) 
	{ 
			 if (HOUSE.ModelVariant == VAC_SILENT)												Kmult = PARAMS.SumRhoKs[country][PhaseSeverity][1];
		else if ((HOUSE.ModelVariant == K_SEROPOS) || (HOUSE.ModelVariant == SIMPLE_NUMERICAL)) Kmult = PARAMS.SumRhoKs[country][PhaseSeverity][0];
		// else don't do anything: AS_PRIME doesn't have common factor
	}
	if (HOUSE.SeroSpecificEfficacies && HOUSE.SeroSpecific_K_values)							Kmult = 1;		/// again not true, but likelihood is the same up to a constant scaling factor for the various serotypes. each serotype efficacy has a corresponding serotype RR, i.e. not a common factor K. 

#ifdef PRINT_WITHIN_LIKE_FUNCTIONS
	std::cout << "l_Int_Vac_Haz: K " << Kmult << " SeroStatus " << int(BaselineSeroStatus) << " country " << country << endl;
#endif

#pragma omp parallel for schedule(static,1) reduction(+:LikeComponent)	
	for (int thread_no = 0; thread_no < HOUSE.max_threads; thread_no++)
	{
		DType AmountToAdd		= NULL;
		DType AgeEffMultiplier	= NULL; //// unless HOUSE.ASVEs, this will equal one. 
		int person_i;
		for (int i = thread_no; i < DATA.Set_Array[country][Set_Index].size(); i += HOUSE.max_threads)
		{
			person_i		= DATA.Set_Array[country][Set_Index][i];
			AmountToAdd		= PARAMS.IVH_vals[person_i][PhaseSeverity];

			if (HOUSE.AdditiveSSASVEs)
				if (HOUSE.SeroSpecific_K_values || (HOUSE.ModelVariant == K_SEROPOS && BaselineSeroStatus == SeroPos)) //// See note in overload of Choose_SumVacHazMult. Essentially if both SSVEs SSKs and Additive age effects  need entire Vaccine hazard multiplier, summed over all serotypes, in this #prgama loop, as no common factors can be taken out....
					AmountToAdd *= Choose_SumVacHazMult(			BaselineSeroStatus, country, DATA.ai_s[person_i], PhaseSeverity, PARAMS, HOUSE);
				else							///// .... however, if NOT doing SSKs, then  can take out a common factor at the end (which also helps with the rest of the code, inc. Change_ functions). K=1 here so that multiplying through by common multiplier doesn't "double count". 			
					AmountToAdd *= Choose_SumVacHazMult((DType)1,	BaselineSeroStatus, country, DATA.ai_s[person_i], PhaseSeverity, PARAMS, HOUSE);					//// these two should be equivalent. 
			else
			{
				if (HOUSE.EffNegWane == EffNegWane_Option::FROM_ZERO)
					AgeEffMultiplier = abs(PARAMS.AgeEff_Mult[BaselineSeroStatus][DATA.ai_s[person_i]]);	//// When SSASVEs multiplicative, you take "efficacies" out of patient sum - so you can't use them here. Hence CalcSumRhoAggregateTransImmunity function above isn't good. Also see note whereby you use absolute value if HOUSE.AltWaneEffNeg, but not otherwise (previous and default). 
				else
					AgeEffMultiplier = PARAMS.AgeEff_Mult[BaselineSeroStatus][DATA.ai_s[person_i]]; //// i.e. either MultSSASVES and Regular waning or nonSSVEs or non-ASVEs (in which case PARAMS.AgeEff_Mult[BaselineSeroStatus][DATA.ai_s[person_i]] = 1). 
				
				AmountToAdd *= AgeEffMultiplier; 
				AmountToAdd *= PARAMS.AgeHaz_Mult[DATA.ai_s[person_i]];	//// unless HOUSE.AS_Hazards, this will equal one. 
				
				if (HOUSE.ModelVariant == K_SEROPOS && BaselineSeroStatus == SeroPos)
					AmountToAdd *= PARAMS.Meta_KplusValues[With_Effs][country][DATA.ai_s[person_i]][PhaseSeverity]; /// this cannot be taken out as a common factor for all patients. 

				if (HOUSE.ModelVariant == AS_PRIME)
					if (BaselineSeroStatus == SeroPos)
						AmountToAdd *= PARAMS.KPlusPrimeValues[country][DATA.ai_s[person_i]][PhaseSeverity];
					else if (BaselineSeroStatus == SeroNeg)
						AmountToAdd *= PARAMS.SumRhoK0_SNeg_Primes[country][PhaseSeverity][DATA.ai_s[person_i]];
			}
			LikeComponent += AmountToAdd;
		}
	}
	// Multiply by common factors where you can. 
	LikeComponent *= Kmult;

	if (HOUSE.AdjHaz & !HOUSE.AdditiveSSASVEs) //// as if HOUSE.AdditiveSSASVEs then AdjHaz already accounted for. 
		LikeComponent *= PARAMS.BS_BaseHazMults[BaselineSeroStatus];

	return LikeComponent;
}
void CalcAndStoreSumIntVacHazOverPatients	(char BaselineSeroStatus, const int &country, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, char PhaseSeverity, DType *** &VacHazLikes)
{
	///// if doing SSVEs etc. the different likelihood components all use sum_i(\int(IntVacHaz(person_i)). If common factors can be taken out (i.e. do not depend on patient or i, e.g. Ks (sometimes)), then they are. Otherwise, they're left in. 
	VacHazLikes[country][BaselineSeroStatus][PhaseSeverity] = l_Int_Vac_Haz(BaselineSeroStatus, country, DATA, PARAMS, HOUSE, PhaseSeverity);
}
void CalcAndStoreSumIntVacHazOverPatients	(char BaselineSeroStatus, const int &country, const DATA_struct &DATA, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, char PhaseSeverity)
{
	CalcAndStoreSumIntVacHazOverPatients(BaselineSeroStatus, country, DATA, PARAMS, HOUSE, PhaseSeverity, PARAMS.VacHazLikes);
}
DType l_Int_Vac_Haz_sero_i			(const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, const int &country, int serotype, char BaselineSeroStatus, char PhaseSeverity, DType *** &VacHazLike_CheckingArray)
{
	///// For K_SEROPOS seropositive SSKs and SSVEs, because there are no common factors to take out, don't have 4 separate likelihood components for the integrated vaccine hazards. 
	///// However  still have HOUSE.NSTypes_VEs = 4
	//// Solution: for serotype == HOUSE.KS_SSVEsKs_BaselineSerotype, have this function return VacHazLike_CheckingArray[country][BaselineSeroStatus][PhaseSeverity] (having already accounted for efficacy or K, and hence not needing to here). 
	//// For all other serotypes, return 0.
	//// Therefore not required to amend Change_K and Change_Efficacy code - as they will apply K*/K multipliers to 0, which still works. 

	DType sumIntVacHaz_wAge_Efficacy_KsEtc = 0;

	if (HOUSE.ModelVariant == K_SEROPOS && HOUSE.SeroSpecificEfficacies && HOUSE.SeroSpecific_K_values && BaselineSeroStatus == SeroPos) //// this is the only scenario in which efficacy included in l_Int_Vac_Haz (either in #pragma loop or rest of function). 
		if (serotype == HOUSE.KS_SSVEsKs_BaselineSerotype)		sumIntVacHaz_wAge_Efficacy_KsEtc = VacHazLike_CheckingArray[country][BaselineSeroStatus][PhaseSeverity]; ///// Ks , rhos and efficacies already considered. 
		else													sumIntVacHaz_wAge_Efficacy_KsEtc = 0;
	else if (HOUSE.AdditiveSSASVEs)
	{
		if (serotype == HOUSE.SSASVEs_BaselineSerotype)			sumIntVacHaz_wAge_Efficacy_KsEtc = VacHazLike_CheckingArray[country][BaselineSeroStatus][PhaseSeverity]; ///// Ks , rhos and efficacies already considered (true for either SSKs or !SSKs). 
		else													sumIntVacHaz_wAge_Efficacy_KsEtc = 0;
	}
	else
	{
		sumIntVacHaz_wAge_Efficacy_KsEtc = VacHazLike_CheckingArray[country][BaselineSeroStatus][PhaseSeverity]; //// this is your starting point. Will contain Ks if either HOUSE.SS_Ks or HOUSE.SS_VEs, but not both. 
		
		///// Apart from above scenarios, you always multiply by efficacy
		DType EffDummy = PARAMS.Efficacies[PhaseSeverity][serotype][BaselineSeroStatus];
		if (HOUSE.EffNegWane == EffNegWane_Option::FROM_ZERO) EffDummy = abs(EffDummy);

		sumIntVacHaz_wAge_Efficacy_KsEtc *= EffDummy;

		//// If just SS_Ks, don't need to multiply by K as this is covered in l_Int_Vac_Haz (either in #pragma loop for K_SEROPOS variant, or as Kmult in rest of function for VAC_SILENT and SIMPLE_NUMERICAL variants)
		if (HOUSE.SeroSpecificEfficacies) sumIntVacHaz_wAge_Efficacy_KsEtc *= PARAMS.rhos[country][serotype]; //// SS_Ks also have rhos, but they are covered in Kmult of l_Int_Vac_Haz function. 

		if (HOUSE.SeroSpecificEfficacies && HOUSE.SeroSpecific_K_values) //// must multiply by sero-specific K value...
			if (!(HOUSE.ModelVariant == K_SEROPOS && BaselineSeroStatus == SeroPos)) //// ... but not for seropositive K_SEROPOS as is alread included in the l_Int_Vac_Haz #pragma loop (i.e. is patient specific). 
			{
				int PrevInf = (HOUSE.ModelVariant == VAC_SILENT) ? BaselineSeroStatus + 1 : BaselineSeroStatus; //// For VAC_SILENT need K1 for SeroNeg (i.e. for BaselineSeroStatus = 0) and K2 for SeroPos (i.e. BaselineSeroStatus = 1). For SIMPLE_NUMERICAL, need simply either K1 or K2. For K_SEROPOS, doesn't apply as K's in earlier #pragma loop. 
				sumIntVacHaz_wAge_Efficacy_KsEtc *= PARAMS.K_s[country][PhaseSeverity][PrevInf][serotype]; // each serotype efficacy has a corresponding serotype RR, i.e. not a common factor K.
			}
	}

	return sumIntVacHaz_wAge_Efficacy_KsEtc;  
}
DType l_Int_Vac_Haz_sero_i			(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, const int &country, int serotype, char BaselineSeroStatus, char PhaseSeverity)
{
	DType Blah = l_Int_Vac_Haz_sero_i(PARAMS, HOUSE, country, serotype, BaselineSeroStatus, PhaseSeverity, PARAMS.VacHazLikes);
	return Blah;
}
DType l_Int_Vac_Haz_sero_i_long_way	(char BaselineSeroStatus, const int &country, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, char PhaseSeverity, int serotype)	  
{
	DType Blah;

	if (HOUSE.ModelVariant == K_SEROPOS && HOUSE.SeroSpecificEfficacies && HOUSE.SeroSpecific_K_values && BaselineSeroStatus == SeroPos)
		if (serotype == HOUSE.KS_SSVEsKs_BaselineSerotype)		 Blah = l_Int_Vac_Haz(BaselineSeroStatus, country, DATA, PARAMS, HOUSE, PhaseSeverity); else Blah = 0;
	else
	{
		Blah = l_Int_Vac_Haz(BaselineSeroStatus, country, DATA, PARAMS, HOUSE, PhaseSeverity) * PARAMS.Efficacies[PhaseSeverity][serotype][BaselineSeroStatus]; //// not appropriate for K_SEROPOS when doing SSKs and SSVEs, otherwise fine. 
		if (HOUSE.SeroSpecificEfficacies) Blah *= PARAMS.rhos[country][serotype];

		if (HOUSE.SeroSpecificEfficacies && HOUSE.SeroSpecific_K_values) //// must multiply by sero-specific K value. 
																		 //if (HOUSE.ModelVariant != K_SEROPOS) //// but not for this model variant as is alread included in the l_Int_Vac_Haz #pragma loop (i.e. is patient specific). 
		{
			int PrevInf = (HOUSE.ModelVariant == VAC_SILENT) ? BaselineSeroStatus + 1 : BaselineSeroStatus; //// For VAC_SILENT need K1 for SeroNeg (i.e. for BaselineSeroStatus = 0) and K2 for SeroPos (i.e. BaselineSeroStatus = 1). For SIMPLE_NUMERICAL, need simply either K1 or K2. For K_SEROPOS, doesn't apply as K's in earlier #pragma loop. 
			Blah *= PARAMS.K_s[country][PhaseSeverity][PrevInf][serotype]; // each serotype efficacy has a corresponding serotype RR, i.e. not a common factor K.
		}
	}
	return Blah;
}

//// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// 
//// FULL LIKELIHOOD FUNCTIONS  

DType l_full				(const DATA_struct &DATA, DType ** LikeArray, const Housekeeping_Struct &HOUSE)
{
	DType l_full = 0;
	for (int country = 0; country < HOUSE.TotalCountries; country++)
		for (int component = 0; component < HOUSE.LCPerC; component++)
			l_full += LikeArray[country][component];
	l_full += log(HOUSE.TimeInterval) * DATA.NoInfected;
	return l_full;
}
DType l_full				(const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	return l_full(DATA, PARAMS.LikeParts, HOUSE);	///// this calculates full Likelihood from unweighted likelihood array. 
}
DType l_full_Minus_Aug		(const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	DType l_full_minus_aug = 0;
	for (int country = 0; country < HOUSE.TotalCountries; country++)
		for (int component = 0; component < HOUSE.LCPerC; component++)
		{
			if (component == HOUSE.L_Indices.BSs[SeroNeg]) continue;
			if (component == HOUSE.L_Indices.BSs[SeroPos]) continue;

			l_full_minus_aug += PARAMS.LikeParts[country][component];
		}
	l_full_minus_aug += log(HOUSE.TimeInterval) * DATA.NoInfected;
	return l_full_minus_aug;
}

void LikelihoodFromScratch	(const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, DType **LikeParts, DType *** &VacHazLike_CheckingArray)
{
	int Country; 
	for (int countryindex = 0; countryindex < HOUSE.NoCountriesToFit; countryindex++)
    {
		Country		= HOUSE.WhichCountries[countryindex];

		//// Prob serostatus components
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			LikeParts[Country][HOUSE.L_Indices.BSs[BaselineSeroStatus]] = l_Aug(Country, BaselineSeroStatus, DATA, PARAMS, HOUSE);	// sum(log(exp(-ha)))  OR sum(log(1-exp(-ha)))

		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
		{
			LikeParts[Country][HOUSE.L_Indices.BaseHaz[PhaseSeverity]] = l_BaseHaz(Country, DATA, PARAMS, HOUSE, PhaseSeverity);	// Baseline hazard for cases	any arm, any serostatus. 

			if (HOUSE.AdjHaz)
				for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
					LikeParts[Country][HOUSE.L_Indices.l_BS_BaseHaz_Mults[PhaseSeverity][BaselineSeroStatus]] = l_BS_BaseHaz_Mult(Country, DATA, PARAMS, HOUSE, PhaseSeverity, BaselineSeroStatus);	///// cases	

			/////	Integrated Baseline Hazards
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
				for (int TrialArm = 0; TrialArm < HOUSE.NumTrialArms_IBH; TrialArm++)
					LikeParts[Country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus][TrialArm]] = l_Int_Base_Haz(Country, DATA, PARAMS, HOUSE, PhaseSeverity, BaselineSeroStatus, TrialArm);		///// survivors	

			//// K's
			for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
				for (int PrevInf = 0; PrevInf < HOUSE.Num_K_Likes; PrevInf++)
					LikeParts[Country][HOUSE.L_Indices.Ks[serotype][PhaseSeverity][PrevInf]] = l_K_case(DATA, PARAMS, HOUSE, Country, PhaseSeverity, PrevInf, serotype);

			//////  K plus values. 
			if (HOUSE.ModelVariant == VAC_SILENT || HOUSE.ModelVariant == K_SEROPOS || HOUSE.ModelVariant == AS_PRIME)
				for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
					LikeParts[Country][HOUSE.L_Indices.Kplus[serotype][PhaseSeverity]] = l_Kplus_case(DATA, PARAMS, HOUSE, Country, PhaseSeverity, serotype);

			if (HOUSE.ModelVariant == AS_PRIME)
				for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
				{
					LikeParts[Country][HOUSE.L_Indices.KPrimes		[serotype][PhaseSeverity]] = l_K0SNegPrime_case	(DATA, PARAMS, HOUSE, Country, PhaseSeverity, serotype);
					LikeParts[Country][HOUSE.L_Indices.KplusPrime	[serotype][PhaseSeverity]] = l_KPlusPrime_case	(DATA, PARAMS, HOUSE, Country, PhaseSeverity, serotype);
				}

			/////	Integrated Vaccine Hazards and Waning Efficacies
	
			for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
				for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
					LikeParts[Country][HOUSE.L_Indices.WaningEffs[serotype][BaselineSeroStatus][PhaseSeverity]] = l_WaningEfficacy(BaselineSeroStatus, Country, DATA, PARAMS, HOUSE, serotype, PhaseSeverity);
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			{
#ifdef VACHAZLIKE_SLOW_WAY
				for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
					LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] = l_Int_Vac_Haz_sero_i_long_way(Set_Array, BaselineSeroStatus, Country, DATA, PARAMS, HOUSE, PhaseSeverity, serotype);
#else
				CalcAndStoreSumIntVacHazOverPatients(BaselineSeroStatus, Country, DATA, PARAMS, HOUSE, PhaseSeverity, VacHazLike_CheckingArray);
				for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
					LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] = l_Int_Vac_Haz_sero_i(PARAMS, HOUSE, Country, serotype, BaselineSeroStatus, PhaseSeverity, VacHazLike_CheckingArray);	///// survivors
#endif
			}
			///// rho cases. 
			if ((HOUSE.SeroSpecificEfficacies) || (HOUSE.SeroSpecific_K_values))
				for (int serotype = 0; serotype < HOUSE.N_STypes; serotype++)
					LikeParts[Country][HOUSE.L_Indices.rhos[serotype][PhaseSeverity]] = l_rho_case(DATA, PARAMS, HOUSE, Country, serotype, PhaseSeverity);

			///// AS_Hazards cases. 
			if (HOUSE.AS_Haz != Age_Option::INDEPENDENT)
				LikeParts[Country][HOUSE.L_Indices.AgeHazMult[PhaseSeverity]] = l_AS_HazMult(Country, DATA, PARAMS, HOUSE, PhaseSeverity);
		}
    }
}
void LikelihoodFromScratch	(const DATA_struct &DATA, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, DType **LikeParts)
{
	LikelihoodFromScratch(DATA, PARAMS, HOUSE, LikeParts, PARAMS.VacHazLikes);
}

DType LL_Total	(const DATA_struct &DATA, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, Augment_Struct &AUG) 
{
	/*
		distinct from l_full function, which calculates log likelihood by summing log like components.
		This calculates "Total likelihood", where non-immune subset patients are not assumend to have a particular serostatus - i.e. we don't use augmented data! For non-immunongenicity subset patients, we use P(O).
		Format of function:

		two #pragma loops, 
			first calculates Sum_i^N(	log(P(O))) = Product_i^N(	log(	P(O|S-)P(S-) + P(O|S+)P(S+)	)	) for augmented patients
			second calculates Sum_i^N(	log(P(O|S)P(S))) for immunogenicity subset (non augmented patients); 
	*/

	DType Total_logLike = 0;

	DType * Total_logLikeArray = new DType[HOUSE.max_threads](); /// for parrelization. 

	/*
		1) Augmented data patients (Non-ImmunogenicitySubset)
	*/
#pragma omp parallel for schedule(static,1) 
	for (int thread_no = 0; thread_no < HOUSE.max_threads; thread_no++)
	{
		DType DummyValue = 0; //// this is a value that does nothing - Your ProbSerostatusGivenOutcome function returns a DType, but also calculates quantities in AUG structure (e.g. AUG.P_OGivenS_PS below)
		for (int AugPatient = thread_no; AugPatient < DATA.NoAugmented; AugPatient += HOUSE.max_threads)
		{
			///// Choose/Propose new value for augmented data patient (choose in case of Gibbs sampling, propose in case of MCMC). 
			Adjust_Aug_thread(thread_no, DATA.AugmentedIndices[AugPatient], DATA, PARAMS, HOUSE, AUG);
			DummyValue						= ProbSerostatusGivenOutcome(DATA.AugmentedIndices[AugPatient], DATA, PARAMS, HOUSE, AUG, thread_no);
			Total_logLikeArray[thread_no]	+= log(		AUG.P_OGivenS_PS[thread_no][		DATA.Ii_s[DATA.AugmentedIndices[AugPatient]]	]			);
		}
	}

	/*
		2) ImmunogenicitySubset (Non-Augmented data patients)
	*/
#pragma omp parallel for schedule(static,1) 
	for (int thread_no = 0; thread_no < HOUSE.max_threads; thread_no++)
	{
		DType DummyValue = 0; //// this is a value that does nothing - Your ProbSerostatusGivenOutcome function returns a DType, but also calculates quantities in AUG structure (e.g. AUG.P_Outcome below)
		for (int NonAugPatient = thread_no; NonAugPatient < DATA.NonAugmentedIndices.size(); NonAugPatient += HOUSE.max_threads)
		{
			///// Choose/Propose new value for augmented data patient (choose in case of Gibbs sampling, propose in case of MCMC). 
			Adjust_Aug_thread(thread_no, DATA.NonAugmentedIndices[NonAugPatient], DATA, PARAMS, HOUSE, AUG);

			DummyValue						= ProbSerostatusGivenOutcome(DATA.NonAugmentedIndices[NonAugPatient], DATA, PARAMS, HOUSE, AUG, thread_no);
			Total_logLikeArray[thread_no]	+= log(AUG.P_Outcome[thread_no]);
		}
	}

	//// multiply all threads together
	for (int thread_no = 0; thread_no < HOUSE.max_threads; thread_no++)
		Total_logLike += Total_logLikeArray[thread_no];

	delete Total_logLikeArray;

	return Total_logLike;
}

bool CompareLikelihoods		(string ParamsOrAugment, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, DType ** LikePartsForChecking, DType *** &VacHazLike_CheckingArray)
{
#ifdef PRINT_WITHIN_CHANGE_FUNCTIONS
	std::cout << "CompareLikelihoods " << ParamsOrAugment << endl;
#endif

	bool LikeAndLikeAlike = 1;

	//// clear to zero
	for (int country = 0; country < HOUSE.TotalCountries; country++)
		for (int component = 0; component < HOUSE.LCPerC; component++)
			LikePartsForChecking[country][component] = 0;

	LikelihoodFromScratch(DATA, PARAMS, HOUSE, LikePartsForChecking, VacHazLike_CheckingArray);

#ifdef PRINT_WITHIN_CHANGE_FUNCTIONS
	std::cout << "CompareLikelihoods finished  LikelihoodFromScratch"  << endl;
#endif
	DType AbsDiff_threshold = 0.0001;

	for (int country = 0; country < HOUSE.TotalCountries; country++)
		for (int component = 0; component < HOUSE.LCPerC; component++)
		{
			if (abs(PARAMS.LikeParts[country][component] - LikePartsForChecking[country][component]) > AbsDiff_threshold)  // i.e. must differ by more than floating point error
			{
				LikeAndLikeAlike = 0;

				std::cerr << endl << ParamsOrAugment << " check discrepancy: comp = " << component << " country " << country << endl;
				std::cerr << "1st is: " << PARAMS.LikeParts[country][component] << ", 2nd is: " << LikePartsForChecking[country][component] << endl;
			}
		}

	string ComponentString = "";
	DType AbsDiff = 0;
	if (LikeAndLikeAlike == 0)
		for (int country = 0; country < HOUSE.TotalCountries; country++)
			for (int component = 0; component < HOUSE.LCPerC; component++)
			{
				ComponentString = WhichLikeComp(component, HOUSE);

				AbsDiff = abs(PARAMS.LikeParts[country][component] - LikePartsForChecking[country][component]);

				if (AbsDiff > AbsDiff_threshold)
				{
					std::cerr << " country " << country << ", lc" << component << " " << ComponentString <<
						" Fast " << PARAMS.LikeParts[country][component] << " slow " << LikePartsForChecking[country][component] <<
						" %Diff " << AbsDiff / LikePartsForChecking[country][component] << endl;
			}
		}

	///// calculate and output full likelihood. 
	if (LikeAndLikeAlike == 0)
	{
		DType L_Full_quick = 0, L_Full_Check = 0;

		std::cout << "PARAMS.LikeSize " << PARAMS.LikeSize << endl;

		for (int country = 0; country < HOUSE.TotalCountries; country++)
			for (int component = 0; component < HOUSE.LCPerC; component++)
				L_Full_quick += PARAMS.LikeParts[country][component];
		
		for (int country = 0; country < HOUSE.TotalCountries; country++)
			for (int component = 0; component < HOUSE.LCPerC; component++)
				L_Full_Check += LikePartsForChecking[country][component];

		L_Full_quick += log(HOUSE.TimeInterval) * DATA.NoInfected;
		L_Full_Check += log(HOUSE.TimeInterval) * DATA.NoInfected;

		std::cout << "L_Full_quick " << L_Full_quick << endl;
		std::cout << "L_Full_Check " << L_Full_Check << endl;
	}

#ifdef PRINT_WITHIN_CHANGE_FUNCTIONS
	std::cout << "END CompareLikelihoods " << endl;
#endif

	return LikeAndLikeAlike;
}

void CalcWeightedLike		(const Housekeeping_Struct &HOUSE, DType **w_LikeParts, DType **LikeParts, int LikeSize)
{
	if (HOUSE.HowManyCaseCategories == 2)
	{
		///// Format of this is to copy LikeParts to w_LikeParts, then multiply each PassiveSevere component of w_LikeParts by HOUSE.PS_Weight. 
		for (int country = 0; country < HOUSE.TotalCountries; country++)
			for (int component = 0; component < HOUSE.LCPerC; component++)
				w_LikeParts[country][component] = LikeParts[country][component]; 

		int Country, PhaseSeverity = PassiveSevere;
		for (int countryindex = 0; countryindex < HOUSE.NoCountriesToFit; countryindex++) 
		{
			Country			= HOUSE.WhichCountries[countryindex];

			w_LikeParts[Country][HOUSE.L_Indices.BaseHaz[PhaseSeverity]] = HOUSE.PS_Weight * 
			w_LikeParts[Country][HOUSE.L_Indices.BaseHaz[PhaseSeverity]];	

			/////	Integrated Baseline Hazards
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
				for (int TrialArm = 0; TrialArm < HOUSE.NumTrialArms_IBH; TrialArm++)
					w_LikeParts[Country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus][TrialArm]] = HOUSE.PS_Weight *
					w_LikeParts[Country][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus][TrialArm]];		///// survivors	

			//// K's (int base haz already covered)
			for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
				for (int PrevInf = 0; PrevInf < HOUSE.Num_K_Likes; PrevInf++)
					w_LikeParts[Country][HOUSE.L_Indices.Ks[serotype][PhaseSeverity][PrevInf]] = HOUSE.PS_Weight *
					w_LikeParts[Country][HOUSE.L_Indices.Ks[serotype][PhaseSeverity][PrevInf]];

			////// then K plus values. 
			if ((HOUSE.ModelVariant == VAC_SILENT) || (HOUSE.ModelVariant == K_SEROPOS))
				for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
					w_LikeParts[Country][HOUSE.L_Indices.Kplus[serotype][PhaseSeverity]] = HOUSE.PS_Weight *
					w_LikeParts[Country][HOUSE.L_Indices.Kplus[serotype][PhaseSeverity]];

			/////	Integrated Vaccine Hazards and Waning Efficacies
			for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
				for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
					w_LikeParts[Country][HOUSE.L_Indices.WaningEffs[serotype][BaselineSeroStatus][PhaseSeverity]] = HOUSE.PS_Weight *
					w_LikeParts[Country][HOUSE.L_Indices.WaningEffs[serotype][BaselineSeroStatus][PhaseSeverity]];

			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
				for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
					w_LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] = HOUSE.PS_Weight *
					w_LikeParts[Country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]];	///// survivors

			///// rho cases. 
			if ((HOUSE.SeroSpecificEfficacies) || (HOUSE.SeroSpecific_K_values))
				for (int serotype = 0; serotype < HOUSE.N_STypes; serotype++)
					w_LikeParts[Country][HOUSE.L_Indices.rhos[serotype][PhaseSeverity]] = HOUSE.PS_Weight *
					w_LikeParts[Country][HOUSE.L_Indices.rhos[serotype][PhaseSeverity]];
		}
	}
}
void CalcWeightedLike		(const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	CalcWeightedLike(HOUSE, PARAMS.w_LikeParts, PARAMS.LikeParts, PARAMS.LikeSize);
}
void CalcWeightedLikes		(const DATA_struct &DATA, Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE)
{
	CalcWeightedLike(CurrentPARAMS , HOUSE);
	CalcWeightedLike(ProposedPARAMS, HOUSE);

	CurrentPARAMS. w_LikeFull	= l_full(DATA, CurrentPARAMS. w_LikeParts, HOUSE);
	ProposedPARAMS.w_LikeFull	= l_full(DATA, ProposedPARAMS.w_LikeParts, HOUSE);
}

//// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// 
//// CHECKING FUNCTIONS

bool Likelihood_Okay		(Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS)
{
	bool AllGood = 1; 

	if ((ProposedPARAMS.LikeFull > 0) || (CurrentPARAMS.LikeFull > 0) || isnan(CurrentPARAMS.LikeFull)) 
	{
		std::cerr << "Breaking MCMC - either NaNs or LL > 0" << endl;
		std::cerr << "ProposedPARAMS.LikeFull " << ProposedPARAMS.LikeFull << endl;
		std::cerr << "CurrentPARAMS.LikeFull " << CurrentPARAMS.LikeFull << endl;
		std::cerr << endl;
		AllGood = 0; 
	}
	return AllGood; 
}
bool VacHazLikes_Okay		(const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	bool All_Equal = 1;
#ifdef TEST_VACCINE_LIKELIHOODS
	DType AbsDiff					= 0.0001;
	DType VacHazLikeLongWay			= 0;
	DType VacHazLikeLongWay_sero_i	= 0;
	DType LLVector_entry = 0; 
	DType VacHazLike_w_efficacies	= 0; 

	///// check for all countries, serostatuses, and Phases/Severities.
	for (int country = 0; country < HOUSE.TotalCountries; country++)
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			{
				if (All_Equal == 0) break; 

				VacHazLikeLongWay = 0;
				for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
				{
					VacHazLikeLongWay_sero_i	= 0; ///// reset
					VacHazLikeLongWay_sero_i	= l_Int_Vac_Haz_sero_i_long_way(BaselineSeroStatus, country, DATA, PARAMS, HOUSE, PhaseSeverity, serotype); ///  store

					////// check equal to what is in like vector
					LLVector_entry = PARAMS.LikeParts[country][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]];

					if (abs(VacHazLikeLongWay_sero_i - LLVector_entry) > AbsDiff)
					{
						std::cerr << "VacHazLikeLongWay_sero_i NOT EQUAL TO PARAMS.LikeParts entry " << endl;
						std::cerr << "VacHazLikeLongWay_sero_i " << VacHazLikeLongWay_sero_i << ", PARAMS.LikeParts entry " << LLVector_entry << endl;
						std::cerr << " country " << country << ", BaselineSeroStatus " << BaselineSeroStatus << ", PhaseSeverity " << PhaseSeverity << ", serotype " << serotype << endl;
						All_Equal = 0;
					}

					/// increment. 
					VacHazLikeLongWay			+= VacHazLikeLongWay_sero_i; 
				}
				if (All_Equal == 0) break;

				if (HOUSE.Single_SNeg_Eff && HOUSE.Single_SPos_Eff)
				{
					VacHazLike_w_efficacies = PARAMS.VacHazLikes[country][BaselineSeroStatus][PhaseSeverity] * PARAMS.Efficacies[PhaseSeverity][2][BaselineSeroStatus]; 

					///// sums will only be equal if efficacies equal between serotype
					if (abs(VacHazLikeLongWay - VacHazLike_w_efficacies) > AbsDiff)
					{
						std::cerr << "VacHazLikeLongWay NOT EQUAL TO PARAMS.VacHazLikes " << endl;
						std::cerr << "VacHazLikeLongWay " << VacHazLikeLongWay << ", CurrentWay " << PARAMS.VacHazLikes[country][BaselineSeroStatus][PhaseSeverity] << endl;
						std::cerr << " country " << country << ", BaselineSeroStatus " << BaselineSeroStatus << ", PhaseSeverity " << PhaseSeverity << endl;

						All_Equal = 0;
						break;
					}
				}
			}
	if (!All_Equal) std::cout << "VacHazLikes check failed: " << endl;

#endif
	return All_Equal; 
}
bool Rhos_Okay				(const Params_Struct &CurrentPARAMS, const Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE)
{
	bool AllCool = 1;
#ifdef Rhos_Okay
	if (HOUSE.SeroSpecificEfficacies || HOUSE.SeroSpecific_K_values)
		if (!HOUSE.BaselinePartition) //// should only sum to one if this condition satisfied i.e. not baseline partition
		{
			DType AbsDiff = 0.0000001, One_dummy = 1; //// latter as you need a DType not an int. 
			for (int country = 0; country < HOUSE.TotalCountries; country++)
			{
				DType SumRho_ThisCountry_Current = 0;
				DType SumRho_ThisCountry_Proposed = 0;
				for (int serotype = 0; serotype < HOUSE.N_STypes; serotype++)
				{
					SumRho_ThisCountry_Current  += CurrentPARAMS. rhos[country][serotype];
					SumRho_ThisCountry_Proposed += ProposedPARAMS.rhos[country][serotype];
				}
				if ((abs(SumRho_ThisCountry_Current - One_dummy) > AbsDiff) || (abs(SumRho_ThisCountry_Proposed - One_dummy) > AbsDiff))
				{
					if ((abs(SumRho_ThisCountry_Current - One_dummy) > AbsDiff)) std::cerr << "SumRho_Current_country  " << country << " NOT EQUAL TO 1" << endl;
					if ((abs(SumRho_ThisCountry_Proposed - One_dummy) > AbsDiff)) std::cerr << "SumRho_Proposed_country " << country << " NOT EQUAL TO 1" << endl;

					AllCool = 0;
					break;
				}
			}
		}
		else
		{
			for (int country = 0; country < HOUSE.TotalCountries; country++)
			{
				for (int serotype = 0; serotype < HOUSE.N_STypes; serotype++)
					if (CurrentPARAMS.rhos[country][serotype] < 0 || ProposedPARAMS.rhos[country][serotype] < 0)
					{
						if (CurrentPARAMS. rhos[country][serotype] < 0) std::cout << "CurrentPARAMS.rhos["  << country << " ][" << serotype << "] " << CurrentPARAMS.rhos[country][serotype] << endl;
						if (ProposedPARAMS.rhos[country][serotype] < 0) std::cout << "ProposedPARAMS.rhos[" << country << " ][" << serotype << "] " << ProposedPARAMS.rhos[country][serotype] << endl;
							
						AllCool = 0; 
						break; 
					}
				
				if (CurrentPARAMS.SumRhos[country] < 1 || ProposedPARAMS.SumRhos[country] < 1)
				{
					std::cout << "CurrentPARAMS. SumRhos[" << country << "] " << CurrentPARAMS.SumRhos[country] << endl;
					std::cout << "ProposedPARAMS.SumRhos[" << country << "] " << ProposedPARAMS.SumRhos[country] << endl;
					AllCool = 0;
					break;
				}
			}
		}
	
	if (!AllCool)	std::cout << "RhosCool check failed:" << endl;
#endif
	return AllCool; 
}
bool SumRhoKs_Okay			(const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, string TestingProposedOrCurrent)
{
	bool AllCool = 1;
#ifdef CHECK_SERO_Ks

	DType AbsDiff = 0.0001;
	DType SumRhoKDummy = 0;
	DType  KtoCompare		= 0;
	string KtoCompare_NAME	= "";

	for (int country = 0; country < HOUSE.TotalCountries; country++)
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int PrevInf = 0; PrevInf < HOUSE.Num_K_Params; PrevInf++)
				{
					SumRhoKDummy = 0;
					
					if (HOUSE.SeroSpecific_K_values)	for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++) SumRhoKDummy += PARAMS.rhos[country][serotype] * PARAMS.K_s[country][PhaseSeverity][PrevInf][serotype];
					else								for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++) SumRhoKDummy += PARAMS.K_s[country][PhaseSeverity][PrevInf][serotype];

					if (abs(PARAMS.SumRhoKs[country][PhaseSeverity][PrevInf] - SumRhoKDummy) > AbsDiff)
					{
						std::cout << TestingProposedOrCurrent << " SumRhoKs NOT RIGHT" << endl;
						std::cerr << " country " << country << ", PhaseSeverity " << PhaseSeverity << ", PrevInf " << PrevInf << 
							", SumRhoK " << PARAMS.SumRhoKs[country][PhaseSeverity][PrevInf] << " SumRhoKDummy " << SumRhoKDummy << endl;

						AllCool = 0;
						break;
					}
				}
#endif
	return AllCool;
}
bool SumRhoKs_Okay			(const Params_Struct &CurrentPARAMS, const Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE)
{
	return (SumRhoKs_Okay(CurrentPARAMS , HOUSE, "Current ") & SumRhoKs_Okay(ProposedPARAMS, HOUSE, "Proposed")); 
}
bool Everything_Okay		(int param_no, const DATA_struct &DATA, Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE, string ParamOrAugment)
{
	bool All_Good = 1;

	if (!Rhos_Okay(CurrentPARAMS, ProposedPARAMS, HOUSE))				All_Good = 0;
	if (!VacHazLikes_Okay(DATA, CurrentPARAMS, HOUSE))					All_Good = 0;
	if (!SumRhoKs_Okay(CurrentPARAMS, ProposedPARAMS, HOUSE))			All_Good = 0;

#ifdef CHECK_STRUCTURES_EQUAL
	if (!CurrentPARAMS.IsEqualToAnother(ProposedPARAMS, HOUSE, DATA))
	{
		std::cout << "ParamSTRUCTs NOT EQUAL! before update of " << ParamOrAugment << endl;
		CurrentPARAMS.WhichPartsUnequal(ProposedPARAMS, HOUSE, DATA);
		ProposedPARAMS.SetEqualToAnother(CurrentPARAMS, HOUSE, DATA);
		All_Good = 0;
	}
#endif
	return All_Good;
}

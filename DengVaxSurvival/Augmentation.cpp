
#include "HeaderAndSwitches.h"
#include "Probability.h"
#include "Splines.h"
#include "CalcParamsEtc.h"
#include "Likelihood.h"
#include "ReadAndProcessData.h"


void PopulateSets			(DATA_struct &DATA, const Housekeeping_Struct &HOUSE) 
{
	bool ConditionToSkip; 
	for (int i = 0; i < NPat; i++)
	{
		ConditionToSkip = !(std::any_of(HOUSE.WhichCountries.begin(), HOUSE.WhichCountries.end(), [&](int j) {return j == DATA.ci_s[i]; }));//// if not doing a particular country, skip that country

		if (ConditionToSkip) continue;

		//// Add to strata sets (their trial arm and serostatus, and "either arm" and serostatus) Don't need to populate (or clear) EitherSeroStatus (as will be unchanged by augmentation). 
		(DATA.Set_Array[DATA.ci_s[i]][	HOUSE.Strata[DATA.Vi_s[i]		][DATA.Ii_s[i]]		]).push_back(i);
		(DATA.Set_Array[DATA.ci_s[i]][	HOUSE.Strata[EitherTrialArm		][DATA.Ii_s[i]]		]).push_back(i);

		//// Do the same for Cases_Array
		if (DATA.IsCase[i])
		{
			(DATA.Cases_Array[DATA.ci_s[i]][	HOUSE.Strata[DATA.Vi_s[i]		][DATA.Ii_s[i]]	]	[DATA.Case_PhaseSeverity[i]]).push_back(i);
			(DATA.Cases_Array[DATA.ci_s[i]][	HOUSE.Strata[EitherTrialArm		][DATA.Ii_s[i]]	]	[DATA.Case_PhaseSeverity[i]]).push_back(i);
		}
	}
}
void ClearSets				(DATA_struct &DATA, const Housekeeping_Struct &HOUSE)
{
	bool ConditionToSkip;			
	for (int Country = 0; Country < HOUSE.TotalCountries; Country++)
	{
		ConditionToSkip = !(std::any_of(HOUSE.WhichCountries.begin(), HOUSE.WhichCountries.end(), [&](int i){return i == Country; }));//// if not doing a particular country, skip that country
		if (ConditionToSkip) continue;

		for (int TrialArm = 0; TrialArm <= HOUSE.HowManyTrialArms; TrialArm++)	//// note the <= as want "either arm". 
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)	
			{
				(DATA.Set_Array[Country][	HOUSE.Strata[TrialArm][BaselineSeroStatus]	]).clear();
				for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
					(DATA.Cases_Array[Country][	HOUSE.Strata[TrialArm][BaselineSeroStatus]	][PhaseSeverity]).clear();	
			}
	}
}
void UpdateOrResetAugData	(int * IisToChange, int * SourceIis, DATA_struct &DATA, const Housekeeping_Struct &HOUSE)
{
	SetEqual_1D_Arrays(IisToChange, SourceIis, NPat); 
	//// These two commands only relevant if IisToChange = DATA.Iis, but include here anyway. 
	ClearSets	(DATA, HOUSE); 
	PopulateSets(DATA, HOUSE);
}
void InitAugDataValues		(DATA_struct &DATA, const Params_Struct &CurrentPARAMS, const Housekeeping_Struct &HOUSE)
{
	DType ProbNaive, hci, RandomNumber = NULL;

	string Old_AugState_FileName = "Output\\" + HOUSE.OldChainFileNamePrefix_Aug + HOUSE.OutputStringForOldChainInput_Aug + ".txt"; 
	if (HOUSE.StartFromPreviousChain_Aug && FileExists(Old_AugState_FileName))
	{
		std::cerr << "Initializing Augmented data using " << HOUSE.OldChainFileNamePrefix_Aug << " from previous chain: " << Old_AugState_FileName << endl;
		std::ifstream Ii_sFile_input;
		Ii_sFile_input.open(Old_AugState_FileName);
		for (int patient = 0; patient < NPat; patient++)	Ii_sFile_input >> DATA.Ii_s[patient];
		Ii_sFile_input.close();
	}
	else
	{
#ifdef PRINT_PROGRAM_PROGRESS
		std::cerr << "Initializing Augmented data using random numbers and historical hazards" << endl;
#endif
		for (int patient = 0; patient < DATA.NoAugmented; patient++)
		{
			hci = CurrentPARAMS.ParamVec[CurrentPARAMS.ParamNumbers.Min_HHaz + DATA.ci_s[DATA.AugmentedIndices[patient]]];
			ProbNaive = PNaiveBaseline(hci, DATA.ai_s[DATA.AugmentedIndices[patient]]);
			RandomNumber = ranf_mt(0);
			DATA.Ii_s[DATA.AugmentedIndices[patient]] = (RandomNumber < ProbNaive) ? 0 : 1;
		}
	}
}
void Choose_iHazMults		(int th, int person_i, int BaselineSeroStatus,	const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, Augment_Struct &AUG)
{
	///// Format of this function is, for various comibinations of model variant, serostatus, and trial arm etc.: 
	//// i)  Choose multiplier of baseline hazard for each PhaseSeverity
	//// ii) Choose multiplier of vaccine hazard, i.e. int_0^t(\lambda(t) & W(t)), for each PhaseSeverity. 
	
	int country = DATA.ci_s[person_i];

	///// iBaseHaz_Mult. 
	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
		AUG.iBaseHaz_Mult[th][PhaseSeverity][BaselineSeroStatus] = Choose_SumBaseHazMult(DATA.Vi_s[person_i], BaselineSeroStatus, country, DATA.ai_s[person_i], PhaseSeverity, PARAMS, HOUSE);

	///// iVacHaz_Mult. 
	if (DATA.Vi_s[person_i] == VaccineGroup)
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			AUG.iVacHaz_Mult[th][PhaseSeverity][BaselineSeroStatus] = Choose_SumVacHazMult(BaselineSeroStatus, country, DATA.ai_s[person_i], PhaseSeverity, PARAMS, HOUSE);
}
void Choose_iHazMults		(int th, int person_i,							const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, Augment_Struct &AUG)
{
	//// Overload that chooses (and stores in relevant thread) multipliers of integrated baseline and vaccine hazards for that patient for both serostatus. Need both BSs when augmenting or doing non-aug likelihood. Need only one (current value of DATA.Iis) fo survival curves. 
	for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
		Choose_iHazMults(th, person_i, BaselineSeroStatus, DATA, PARAMS, HOUSE, AUG);
}

void Choose_Case_Ks					(int th, int person_i, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, Augment_Struct &AUG)
{
	int Case_K_serotype		= (HOUSE.SeroSpecific_K_values)		? DATA.CaseSerotype[person_i] : 0; //// will be set to MDValue (= 9999) if not set properly. If doing	SS_VEs		and not		SS_Ks	, then still only one K value	per prior exposure		, which affects L_Indices etc., hence why you need to be careful here. 

	if (DATA.IsCase[person_i])
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			AUG.Case_Ks[th][BaselineSeroStatus] =
				Choose_K(DATA.Vi_s[person_i], BaselineSeroStatus, DATA.ci_s[person_i], DATA.ai_s[person_i], DATA.Case_PhaseSeverity[person_i], Case_K_serotype, PARAMS, HOUSE);
}
void Choose_HazMults_SingleSero		(int th, int person_i, int BaselineSeroStatus,	const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, Augment_Struct &AUG)
{
	//// Called when doing Hazard ratios - only when either HOUSE.SeroSpecific_K_values OR HOUSE.SeroSpecificEfficacies
	int country		= DATA.ci_s[person_i];
	int Age			= DATA.ai_s[person_i];

	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
		for (int serotype = 0; serotype < HOUSE.N_STypes; serotype++)
		{
			//// baseline hazard multipliers
			AUG.HR_BaseHaz_seroMults[th][PhaseSeverity][BaselineSeroStatus][serotype] =
				Choose_BaseHazMult(DATA.Vi_s[person_i], BaselineSeroStatus, country, Age, PhaseSeverity, serotype, PARAMS, HOUSE);

			//// vaccine hazard multipliers
			if (DATA.Vi_s[person_i] == VaccineGroup)	
					AUG.HR_VacHaz_seroMults[th][PhaseSeverity][BaselineSeroStatus][serotype] = Choose_VacHazMult(BaselineSeroStatus, country, Age, PhaseSeverity, serotype, PARAMS, HOUSE);
			else	AUG.HR_VacHaz_seroMults[th][PhaseSeverity][BaselineSeroStatus][serotype] = 0;
		}
}
void Choose_HazMults_SingleSero		(int th, int person_i,							const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, Augment_Struct &AUG)
{
	//// Overload that chooses (and stores in relevant thread) multipliers of baseline and vaccine hazards for that patient for both serostatus.
	for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
		Choose_HazMults_SingleSero(th, person_i, BaselineSeroStatus, DATA, PARAMS, HOUSE, AUG);
}
void Adjust_Aug_thread				(int th, int person_i, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, Augment_Struct &AUG)
{
	///// Function changes various quantities in Augment_Struct &AUG that you use in both Change_Aug_Data function (for likelihood) and in ProbSerostatusGivenOutcome (which calculates probabilities for Gibb's augmentation). 

	int country				= DATA.ci_s[person_i];
	int Case_K_serotype		= (HOUSE.SeroSpecific_K_values)		? DATA.CaseSerotype[person_i] : 0; //// will be set to MDValue (= 9999) if not set properly. If doing	SS_VEs		and not		SS_Ks	, then still only one K value	per prior exposure		, which affects L_Indices etc., hence why you need to be careful here. 
	int Case_Eff_serotype	= (HOUSE.SeroSpecificEfficacies)	? DATA.CaseSerotype[person_i] : 0; //// will be set to MDValue (= 9999) if not set properly. If doing	SS_Ks		and not		SS_VEs	, then still only one Eff value per baseline serostatus	, which affects L_Indices etc., hence why you need to be careful here. 

	DType hci = PARAMS.ParamVec[PARAMS.ParamNumbers.Min_HHaz + country];

	///// ///// ///// ///// ///// ///// ///// ///// 
	//// Augmentation multipliers
	
	// if Seronegative (and augmenting to Seropositive),	SUBTRACT	all Seronegative terms and		ADD				all the Seropositive, 
	// if Seropositive (and augmenting to Seropositive),	ADD			all Seronegative terms and		SUBTRACT		all the Seropositive 

	if (DATA.Ii_s[person_i] == 0)
	{
		AUG.Mults[th][SeroPos] = 1;
		AUG.Mults[th][SeroNeg] = -1;
	}
	else
	{
		AUG.Mults[th][SeroPos] = -1;
		AUG.Mults[th][SeroNeg] =  1;
	}

	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
	///// Calculate integrated baseline hazard and integrated vaccine hazard for this patient (for both phases/severities if necessary)
	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
	
	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
	{
		//// IntBaseHaz
		AUG.iBaseHaz[th][PhaseSeverity] = IntBaseHaz(DATA.FollowUp[START][PhaseSeverity][person_i], DATA.FollowUp[END][PhaseSeverity][person_i], country, PARAMS);

		//// if baseline partition (i.e. rhos not proportions, but relative multipliers of baseline hazard). Not really doing this anymore. 
		if (HOUSE.SeroSpecificEfficacies || HOUSE.SeroSpecific_K_values)
			if (HOUSE.BaselinePartition) AUG.iBaseHaz[th][PhaseSeverity] *= PARAMS.SumRhos[country];

		//// IntVacHaz
		if (DATA.Vi_s[person_i])
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			{
				if (BaselineSeroStatus == DATA.Ii_s[person_i]) AUG.iVacHaz[th][PhaseSeverity][BaselineSeroStatus] = PARAMS.IVH_vals[person_i][PhaseSeverity]; //// i.e. just use the pre-calculated one for the right BaselineSeroStatus
				else
				{
					int FirstDayVacHazard = DATA.FollowUp[START][PhaseSeverity][person_i];

					if (HOUSE.AllDosesRequired_AG_BS[DATA.AgeGroup1[person_i]][BaselineSeroStatus])
						if (HOUSE.MildAndSevere == TREATED_SEPARATELY || PhaseSeverity == ActiveMild)
							FirstDayVacHazard = DATA.ThirdDose[person_i];

					AUG.iVacHaz[th][PhaseSeverity][BaselineSeroStatus] =	IntVacHaz(	FirstDayVacHazard,
																						DATA.FollowUp	[END	][PhaseSeverity	]	[person_i], 
																						DATA.FollowUp	[START	][ActiveMild	]	[person_i], //// i.e. dose 1 date
																						DATA.SecondDose								[person_i], 
																						DATA.ThirdDose								[person_i], 
																						DATA.ai_s									[person_i], 
																						BaselineSeroStatus, country, PARAMS, HOUSE			);
				}

			}
		else for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)	AUG.iVacHaz[th][PhaseSeverity][BaselineSeroStatus] = 0;
	}


	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
	///// Choose K values			- Augmentations moves subject from one hazard group to another (but arm stays the same). Hence number of K values per trial phase/disease type limited to two. Will obviously be different for SeroSpecific_K_Values. 
	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 

	Choose_iHazMults(th, person_i, DATA, PARAMS, HOUSE, AUG); 

	///// ///// ///// ///// ///// ///// ///// ///// 
	///// Case K's and BaselineHazard...

	Choose_Case_Ks(th, person_i, DATA, PARAMS, HOUSE, AUG); 

}
void Change_Aug_Data_patient		(int person_i, int pi_augindex, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, Augment_Struct &AUG, int th)
{
	//// th is thread_number; 

	//////////////// Augmented Data
	// In Parallel Version, instead of returning the New_LikelihoodComponentVec, return DIFFERENCES made to each component of the vector and the full likelihood (former needed to carry on the MCMC for subsequent parameters in subsequent iterations, latter needed for MCMC accept/reject step. 
	// This function doesn't write new augmented values, just returns the component differences. Sum/implement the differences later

	int country				= DATA.ci_s[person_i];
	int AgeInYears_person_i = DATA.ai_s[person_i];

	int Case_K_serotype		= (HOUSE.SeroSpecific_K_values)		? DATA.CaseSerotype[person_i] : 0; //// will be set to MDValue (= 9999) if not set properly. If doing	SS_VEs		and not		SS_Ks	, then still only one K value	per prior exposure		, which affects L_Indices etc., hence why you need to be careful here. 
	int Case_Eff_serotype	= (HOUSE.SeroSpecificEfficacies)	? DATA.CaseSerotype[person_i] : 0; //// will be set to MDValue (= 9999) if not set properly. If doing	SS_Ks		and not		SS_VEs	, then still only one Eff value per baseline serostatus	, which affects L_Indices etc., hence why you need to be careful here. 

	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
	///// Change Prob(SeroStatus) components - same for all model variants/phases/seroefficacies etc. 
	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 

	AUG.LikeDiffs[pi_augindex][HOUSE.L_Indices.BSs[SeroNeg]] = AUG.Mults[th][SeroNeg] * PARAMS.SeroPrevs[LogIndex][country][SeroNeg][DATA.ai_s[person_i]];	//// Multiplier * subtract																										
	AUG.LikeDiffs[pi_augindex][HOUSE.L_Indices.BSs[SeroPos]] = AUG.Mults[th][SeroPos] * PARAMS.SeroPrevs[LogIndex][country][SeroPos][DATA.ai_s[person_i]];	//// Multiplier * add																											

	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
	///// change waning efficacy					- same for all model variants 
	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 

	if (DATA.IsCase[person_i] && DATA.Vi_s[person_i])
	{
		DType AggregateEff_singleStype			= NULL; 
		DType ZeroOrOneToSubtractFromWaningMult = NULL; 
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			AUG.LikeDiffs[pi_augindex][HOUSE.L_Indices.WaningEffs[Case_Eff_serotype][BaselineSeroStatus][DATA.Case_PhaseSeverity[person_i]]] = AUG.Mults[th][BaselineSeroStatus] * 
				WaningEfficacy(BaselineSeroStatus, Case_Eff_serotype, DATA.Case_PhaseSeverity[person_i], AgeInYears_person_i,
				PARAMS.WaningMults[AgeInYears_person_i][BaselineSeroStatus][DATA.TimePost_Final_Dose[person_i]], HOUSE.EffNegWane, PARAMS, HOUSE, /*NaturalLog =*/ true);
	}

	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
	//// change integrated baseline hazards
	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
	
	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			AUG.LikeDiffs[pi_augindex][HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus][DATA.Vi_s[person_i]]] = 
				-AUG.Mults[th][BaselineSeroStatus] * AUG.iBaseHaz_Mult[th][PhaseSeverity][BaselineSeroStatus] * AUG.iBaseHaz[th][PhaseSeverity];		// Multiplier * adding		a minus hazard  	

	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
	//// change integrated vaccine hazards - same for all model variants - only difference is if sero efficacies, dealt with in values WaningEffInd_SPos and SPosVE_Cases etc. 
	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
	
	if (DATA.Vi_s[person_i])
		for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
				for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
					if (HOUSE.AdditiveSSASVEs)
					{
						if (serotype == HOUSE.SSASVEs_BaselineSerotype)
							AUG.LikeDiffs[pi_augindex][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] =
								AUG.Mults[th][BaselineSeroStatus] * //// Multiplier * add		
								Choose_SumVacHazMult(BaselineSeroStatus, country, DATA.ai_s[person_i], PhaseSeverity, PARAMS, HOUSE) *
								AUG.iVacHaz[th]		[PhaseSeverity][BaselineSeroStatus] ; 
						else AUG.LikeDiffs[pi_augindex][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] = 0;
					}
					else
					{
						//// i.e. (add Multiplier) * AgeHaz_Mult * IntVachaz
						AUG.LikeDiffs[pi_augindex][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]]	 = 
							AUG.Mults[th][BaselineSeroStatus] * 
							PARAMS.AgeHaz_Mult	[DATA.ai_s[person_i]] *
							AUG.iVacHaz[th]		[PhaseSeverity]				[BaselineSeroStatus];											//// Multiplier * add		

						//// Multiply iVacHaz by "Efficacy", for which you (may) need both Age and serotype/serotstatus components. If doing HOUSE.AltWaneEffNeg  need the absolute value of "efficacy" (see notes). 
						DType EffSeroMultDummy	= HOUSE.EffNegWane == EffNegWane_Option::FROM_ZERO ? abs(PARAMS.Efficacies[PhaseSeverity][serotype][BaselineSeroStatus]	) : PARAMS.Efficacies[PhaseSeverity][serotype][BaselineSeroStatus];
						DType EffAgeMultDummy	= HOUSE.EffNegWane == EffNegWane_Option::FROM_ZERO ? abs(PARAMS.AgeEff_Mult[BaselineSeroStatus][DATA.ai_s[person_i]]	) : PARAMS.AgeEff_Mult[BaselineSeroStatus][DATA.ai_s[person_i]];
						//// now multiply...
						AUG.LikeDiffs[pi_augindex][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] *= EffSeroMultDummy	;
						AUG.LikeDiffs[pi_augindex][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] *= EffAgeMultDummy	;

						//// Multiplier by correct K value (basically always the Survive K but not if doing both SS_K and SS_VEs)
						if (HOUSE.SeroSpecific_K_values & !HOUSE.SeroSpecificEfficacies)
							AUG.LikeDiffs[pi_augindex][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] *= Choose_Sum_Rho_K(VaccineGroup, BaselineSeroStatus, country, DATA.ai_s[person_i], PhaseSeverity, PARAMS, HOUSE);
						else
							AUG.LikeDiffs[pi_augindex][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] *= Choose_K(VaccineGroup, BaselineSeroStatus, country, DATA.ai_s[person_i], PhaseSeverity, serotype, PARAMS, HOUSE);

						//// multiply by rho if necessary. 
						if (HOUSE.SeroSpecificEfficacies) 
							AUG.LikeDiffs[pi_augindex][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] *= PARAMS.rhos[country][serotype]; 

						//// multiply by BS_BaseHazMult if necessary. 
						if (HOUSE.AdjHaz) //// included in BaseHazMult but in AUG.Survive_Ks, so must multiply here if doing. 
							AUG.LikeDiffs[pi_augindex][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype]] *= PARAMS.BS_BaseHazMults[BaselineSeroStatus];
					}
	
	if (HOUSE.ModelVariant == K_SEROPOS && HOUSE.SeroSpecific_K_values && HOUSE.SeroSpecificEfficacies) //// need to re-do seropositive vaccinees. Seronegatives are fine. 
		if (DATA.Vi_s[person_i])
			for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
				for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
					if (serotype == HOUSE.KS_SSVEsKs_BaselineSerotype)
						AUG.LikeDiffs[pi_augindex][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][SeroPos][serotype]] =
							AUG.Mults[th][SeroPos] * //// Multiplier * add		
							AUG.iVacHaz[th][PhaseSeverity][SeroPos] *
							PARAMS.Meta_KplusValues[With_Effs][country][DATA.ai_s[person_i]][PhaseSeverity]; //// includes efficacies and rhos. 
					else AUG.LikeDiffs[pi_augindex][HOUSE.L_Indices.IntVacHaz[PhaseSeverity][SeroPos][serotype]] = 0;

	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
	//// K multipliers
	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
	
	int K_0_ind = NULL, K_1_ind = NULL, K_2_ind = NULL, K_plus_ind = NULL, SNeg_K_Prime_ind = NULL, SPos_Kplus_Prime = NULL; //// Either phase. Will only be a case in single PhaseSeverity

	if (DATA.IsCase[person_i]) //// need this condition otherwise DATA.CaseSerotype[person_i] will give memory access violation. 
	{
		K_0_ind = HOUSE.L_Indices.Ks[Case_K_serotype][DATA.Case_PhaseSeverity[person_i]][0];
		if (HOUSE.ModelVariant == VAC_SILENT || HOUSE.ModelVariant == SIMPLE_NUMERICAL)								K_1_ind		= HOUSE.L_Indices.Ks	[Case_K_serotype][DATA.Case_PhaseSeverity[person_i]][1];
		if (HOUSE.ModelVariant == VAC_SILENT)																		K_2_ind		= HOUSE.L_Indices.Ks	[Case_K_serotype][DATA.Case_PhaseSeverity[person_i]][2];
		if (HOUSE.ModelVariant == VAC_SILENT || HOUSE.ModelVariant == K_SEROPOS || HOUSE.ModelVariant == AS_PRIME)	K_plus_ind	= HOUSE.L_Indices.Kplus	[Case_K_serotype][DATA.Case_PhaseSeverity[person_i]];
		if (HOUSE.ModelVariant == AS_PRIME)
		{
			SNeg_K_Prime_ind = HOUSE.L_Indices.KPrimes		[Case_K_serotype][DATA.Case_PhaseSeverity[person_i]];
			SPos_Kplus_Prime = HOUSE.L_Indices.KplusPrime	[Case_K_serotype][DATA.Case_PhaseSeverity[person_i]];
		}
	}

	if (	((HOUSE.ModelVariant == VAC_SILENT || HOUSE.ModelVariant == AS_PRIME) && (DATA.Vi_s[person_i] == 0))	||
			 (HOUSE.ModelVariant == K_SEROPOS)																		||
			 (HOUSE.ModelVariant == SIMPLE_NUMERICAL)																)
		if (DATA.IsCase[person_i])
		{
			AUG.LikeDiffs[pi_augindex][K_0_ind	] = AUG.Mults[th][SeroNeg] * log(AUG.Case_Ks[th][SeroNeg]);	//// Multiplier * subtract
			if (HOUSE.ModelVariant == VAC_SILENT || HOUSE.ModelVariant == K_SEROPOS || HOUSE.ModelVariant == AS_PRIME)
				AUG.LikeDiffs[pi_augindex][K_plus_ind	] = AUG.Mults[th][SeroPos] * log(AUG.Case_Ks[th][SeroPos]); //// Multiplier * add
			if (HOUSE.ModelVariant == SIMPLE_NUMERICAL)										
				AUG.LikeDiffs[pi_augindex][K_1_ind	] = AUG.Mults[th][SeroPos] * log(AUG.Case_Ks[th][SeroPos]); //// Multiplier * add
		}
	
	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
	//// Extra changes needed for VAC_SILENT AS_PRIME vaccine group
	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
	
	if (HOUSE.ModelVariant == VAC_SILENT)
	{
		if (DATA.Vi_s[person_i])
			//// change K multipliers
			if (DATA.IsCase[person_i])
			{
				AUG.LikeDiffs[pi_augindex][K_2_ind] = AUG.Mults[th][SeroPos] * log(AUG.Case_Ks[th][SeroPos]);
				AUG.LikeDiffs[pi_augindex][K_1_ind] = AUG.Mults[th][SeroNeg] * log(AUG.Case_Ks[th][SeroNeg]);
			}
	} 
	else if (HOUSE.ModelVariant == AS_PRIME)
	{
		if (DATA.Vi_s[person_i])
			//// change Kprime/KplusPrime 
			if (DATA.IsCase[person_i])
			{
				AUG.LikeDiffs[pi_augindex][SNeg_K_Prime_ind] = AUG.Mults[th][SeroNeg] * log(AUG.Case_Ks[th][SeroNeg]);
				AUG.LikeDiffs[pi_augindex][SPos_Kplus_Prime] = AUG.Mults[th][SeroPos] * log(AUG.Case_Ks[th][SeroPos]);
			}
	}

	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
	//// AdjHaz
	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 

	if (HOUSE.AdjHaz)
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			if (DATA.IsCase_AMandPS[PhaseSeverity][person_i])
				for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
					AUG.LikeDiffs[pi_augindex][HOUSE.L_Indices.l_BS_BaseHaz_Mults[PhaseSeverity][BaselineSeroStatus]] = AUG.Mults[th][BaselineSeroStatus] * log(PARAMS.BS_BaseHazMults[BaselineSeroStatus]); 


	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
	//// Calculate full likelihood. 
	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 

	DType L_full = PARAMS.LikeFull; /// careful: DO NOT change CurrentlogL_Full itself as function passes by reference and calling this function in parallel. 
	for (int component = 0; component < AUG.NoLikeDiffsPerAugPatient; component++)	L_full += AUG.LikeDiffs[pi_augindex][component];
	AUG.LikeDiffs[pi_augindex][AUG.NoLikeDiffsPerAugPatient] = L_full; /// note the lack of minus 1 at the end. L_Full is not a like diff. 
}
DType ProbSerostatusGivenOutcome	(int person_i, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, Augment_Struct &AUG, int th)
{
	///// As well as returning a probability, this function changes various quantities in Augment_Struct &AUG
	DType P_Spos_Given_Outcome = NULL; 
	
	///// Must compute: i) P(S+|C+), or ii) P(S+|C-). 
	///// C refers to whether subject a case (+) or not (-). S refers to whether you are seropositive (+), or seronegative (-) at baseline. P(S+) = 1 - P(S-) and P(C+) = 1 - P(C-). 

	///// Format of gibbs aug main code is. Compute probability of seropositivity given outcome (C+ or C-). If AugRandomNumber > P(S+|Outcome) then assign augmented data Ii = 1, otherwise Ii = 0 
	///// If the assigned value is DIFFERENT from existing value of augmented data point, compute the changes this will make to the likelihood. Do all of the above in parallel, recording the changes made to the likelihood for ALL patients where a different aug value is assigned (accepted), then as before compute the full likelihood as the sum of the differences made. 

	int country				= DATA.ci_s[person_i]; 
	int AgeInYears_person_i = DATA.ai_s[person_i];
	int Case_K_serotype		= (HOUSE.SeroSpecific_K_values)		? DATA.CaseSerotype[person_i] : 0; //// will be set to MDValue (= 9999) if not set properly. If doing	SS_VEs		and not		SS_Ks	, then still only one K value	per prior exposure		, which affects L_Indices etc., hence need to be careful. 
	int Case_Eff_serotype	= (HOUSE.SeroSpecificEfficacies)	? DATA.CaseSerotype[person_i] : 0; //// will be set to MDValue (= 9999) if not set properly. If doing	SS_Ks		and not		SS_VEs	, then still only one Eff value per baseline serostatus	, which affects L_Indices etc., hence need to be careful. 

	//// Compute P(S+) and P(S-). 
	if (HOUSE.ExtImSub == ExtImSub_Option::AS_PROB && DATA.Imputed_ProbSPos[person_i] != MDvalue)
	{
		AUG.P_BS[th][SeroPos] = DATA.Imputed_ProbSPos[person_i];
		AUG.P_BS[th][SeroNeg] = 1 - AUG.P_BS[th][SeroPos];
	}
	else //// i.e. if HOUSE.ExtImSub == ExtImSub_Option::IGNORED or HOUSE.ExtImSub == ExtImSub_Option::AS_DATA (in the latter case not augmenting "extended Immune subset" patients). 
	{
		AUG.P_BS[th][SeroPos] = PARAMS.SeroPrevs[NonLogIndex][country][SeroPos][DATA.ai_s[person_i]];
		AUG.P_BS[th][SeroNeg] = PARAMS.SeroPrevs[NonLogIndex][country][SeroNeg][DATA.ai_s[person_i]];
	}

	//// BaselineHazard
	DType BaseHazValue	= (DATA.IsCase[person_i]) ? PARAMS.BaselineHazardValues[DATA.ci_s[person_i]][DATA.FollowUp[END][PassiveSevere][person_i]] : -10000;				//// // //// In data for all active phase cases, End_PassiveSevere = Start_PassiveSevere = End_ActiveMild. Set to something impossible to cause a crash if incorrectly set or not set at all. 

	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
	///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// 
	///// Now compute: i) P(SEROSTATUS|OUTCOME)

	///// Reset P_Survive to zero and HazCaseMult to 1. 
	for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
	{
		AUG.P_Survive	[th][BaselineSeroStatus] = 0;
		AUG.HazCaseMult	[th][BaselineSeroStatus] = 1;
	}

	if (DATA.Vi_s[person_i] == 0) /// subject in control group
	{
		//// format is: i) add to int hazards; ii) exponentiate; iii) multiply by instantaneous hazard if a case. 
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
		{
			//// add integrated hazards from all trial phases / disease severities. (actually subtract since you need exponential of negative integrated hazard). 
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++) 
				AUG.P_Survive[th][BaselineSeroStatus] -= AUG.iBaseHaz[th][PhaseSeverity] * AUG.iBaseHaz_Mult[th][PhaseSeverity][BaselineSeroStatus];
			
			//// exponentiate to get survival prob. 
			AUG.P_Survive[th][BaselineSeroStatus] = exp(AUG.P_Survive[th][BaselineSeroStatus]); 

			///// if a case then need to multiply survival probability by instantaneous hazard to get density, then by time interval to get probability. 
			if (DATA.IsCase[person_i])
			{
				DType rho_c_d = 1; //// default is to set rho_c_d equal to 1, unless doing serospecific efficacies.  
				if (HOUSE.SeroSpecificEfficacies || HOUSE.SeroSpecific_K_values)	rho_c_d = PARAMS.rhos[country][DATA.CaseSerotype[person_i]];

				AUG.HazCaseMult[th][BaselineSeroStatus] = rho_c_d * HOUSE.TimeInterval * BaseHazValue * AUG.Case_Ks[th][BaselineSeroStatus];  
				if (HOUSE.AdjHaz) AUG.HazCaseMult[th][BaselineSeroStatus] *= PARAMS.BS_BaseHazMults[BaselineSeroStatus];
			}
		}
	}
	else /// subject in vaccine group
	{
		//// format is: i) add to int hazards; ii) exponentiate; iii) multiply by instantaneous hazard if a case. 
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
		{
			//// add integrated hazards from all trial phases / disease severities. (actually subtract since you need exponential of negative integrated hazard). 
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++) //// add integrated hazards from all trial phases / disease severities. (actually subtract coz you need exponential of negative integrated hazard). 
				AUG.P_Survive[th][BaselineSeroStatus] -= (  (AUG.iBaseHaz_Mult	[th][PhaseSeverity][BaselineSeroStatus] * AUG.iBaseHaz[th][PhaseSeverity]) -						////			(IntBaseHaz		* Multiplier (K and may or may not include rhos		) )  minus... 
															(AUG.iVacHaz_Mult	[th][PhaseSeverity][BaselineSeroStatus] * AUG.iVacHaz [th][PhaseSeverity][BaselineSeroStatus]));	//// ...minus	(Int_Vac_Haz	* Multiplier (K and Effs may or may not include rhos) ). No need for age multiplier here (if relevant), as would be included in AUG.iVacHaz
	
			//// exponentiate to get survival prob. 
			AUG.P_Survive[th][BaselineSeroStatus] = exp(AUG.P_Survive[th][BaselineSeroStatus]);

			///// if a case then need to multiply survival probability by instantaneous hazard to get density, then by time interval to get probability. 
			if (DATA.IsCase[person_i])
			{
				DType rho_c_d = (HOUSE.SeroSpecificEfficacies || HOUSE.SeroSpecific_K_values) ? PARAMS.rhos[country][DATA.CaseSerotype[person_i]] : 1; //// default is to set rho_c_d equal to 1, unless doing serospecific efficacies or K's.  
				
				AUG.HazCaseMult[th][BaselineSeroStatus] = rho_c_d * HOUSE.TimeInterval * BaseHazValue * AUG.Case_Ks[th][BaselineSeroStatus] * PARAMS.AgeHaz_Mult[DATA.ai_s[person_i]];
				if (HOUSE.AdjHaz) AUG.HazCaseMult[th][BaselineSeroStatus] *= PARAMS.BS_BaseHazMults[BaselineSeroStatus];
				AUG.HazCaseMult[th][BaselineSeroStatus] *= WaningEfficacy(BaselineSeroStatus, Case_Eff_serotype, DATA.Case_PhaseSeverity[person_i], AgeInYears_person_i,
					PARAMS.WaningMults[AgeInYears_person_i][BaselineSeroStatus][DATA.TimePost_Final_Dose[person_i]], HOUSE.EffNegWane, PARAMS, HOUSE, /*NaturalLog =*/ false);
			}
		}
	}
	//// Calculate P(O|S) x P(S)
	for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
		AUG.P_OGivenS_PS[th][BaselineSeroStatus] = AUG.P_Survive[th][BaselineSeroStatus] * AUG.HazCaseMult[th][BaselineSeroStatus] * AUG.P_BS[th][BaselineSeroStatus];

	AUG.P_Outcome[th]		= AUG.P_OGivenS_PS[th][SeroNeg] + AUG.P_OGivenS_PS[th][SeroPos];
	P_Spos_Given_Outcome	= AUG.P_OGivenS_PS[th][SeroPos] / AUG.P_Outcome[th];	//// i.e.   P(Outcome|S+)P(S+) / P(Outcome)		{ = P(Outcome|S+)P(S+) / P(Outcome|S+)P(S+) + P(Outcome|S-)P(S-)	} 

	/*if (isnan(P_Spos_Given_Outcome) || (P_Spos_Given_Outcome < 0) || (P_Spos_Given_Outcome > 1))
	{
		std::fflush(stderr);

		std::cerr << endl << endl << "P_Spos_Given_Outcome WRONG P th" << th << " " << P_Spos_Given_Outcome << " Numerator " << AUG.P_OGivenS_PS[th][SeroPos] << " Denominator " << AUG.P_Outcome[th] << endl;
		std::cerr << "person_i th" << th << " " << person_i << " serostatus " << DATA.Ii_s[person_i] << " country " << country << " IsCase_AM " << DATA.IsCase_ActiveMild[person_i] << " IsCase_PS " << DATA.IsCase_PassiveSevere[person_i] << " TrialARM " << DATA.Vi_s[person_i] << endl;

		std::cerr << "th" << th << endl; 
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
		{
			std::cerr << "P_Survive_Haz BS" << BaselineSeroStatus << " " << AUG.P_Survive[th][BaselineSeroStatus] << endl;
			std::cerr << "HazCaseMult BS" << BaselineSeroStatus << " " << AUG.HazCaseMult[th][BaselineSeroStatus] << endl;
		}
		std::cerr << "Age Haz Mult " << PARAMS.AgeHaz_Mult[DATA.ai_s[person_i]] << endl;

		for (int countryDummy = 0; countryDummy < 10; countryDummy++)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
				for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
					std::cerr << "PARAMS.SumRhoEffKs country " << countryDummy << " PS " << PhaseSeverity << " BS " << BaselineSeroStatus << " " << PARAMS.SumRhoEffKs[countryDummy][PhaseSeverity][BaselineSeroStatus] << endl;

		for (int countryDummy = 0; countryDummy < 10; countryDummy ++)
			for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
				std::cerr << "PARAMS.rhos countryDummy " << countryDummy << " serotype" << serotype << "  " << PARAMS.rhos[countryDummy][serotype] << endl;

		std::cerr << "Ks" << endl;
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int PrevInf = 0; PrevInf < 3; PrevInf++)
				for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
					std::cerr << "K Phase" << PhaseSeverity << " PrevInf " << PrevInf << " serotype " << serotype << " " << PARAMS.K_s[country][PhaseSeverity][PrevInf][serotype] << endl; 

		std::cerr << "Efficacies (multipliers)" << endl;
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
				for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
					std::cerr << "PARAMS.Effs PS" << PhaseSeverity << " BS " << BaselineSeroStatus << " serotype" << serotype << " " << PARAMS.Efficacies[PhaseSeverity][serotype][BaselineSeroStatus] << endl;

		std::cerr << "Efficacies Ages" << endl;
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			std::cerr << "PARAMS.AgeEff_Mult BS " << BaselineSeroStatus << " age " << DATA.ai_s[person_i] << " " << PARAMS.AgeEff_Mult[BaselineSeroStatus][DATA.ai_s[person_i]] << endl;

		std::cerr << "iBaseHaz" << endl;
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			std::cerr << "AUG.iBaseHaz PS " << PhaseSeverity << " " << AUG.iBaseHaz[th][PhaseSeverity] << endl;

		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			{
				std::cerr << "BaselineSeroStatus " << BaselineSeroStatus << " PhaseSeverity " << PhaseSeverity << endl;
				std::cerr << "iVacHaz " << AUG.iVacHaz[th][PhaseSeverity][BaselineSeroStatus] << endl;
				std::cerr << "iVacHaz_Mult " << AUG.iVacHaz_Mult[th][PhaseSeverity][BaselineSeroStatus] << endl;
				std::cerr << "iBaseHaz_Mult " << AUG.iBaseHaz_Mult[th][PhaseSeverity][BaselineSeroStatus] << endl;
			}
		std::cerr << "P_Survive_Haz_SPos th" << th << " " << AUG.P_Survive[th][SeroPos]	<< " P_Survive_Haz_SNeg "			<< AUG.P_Survive[th][SeroNeg]	<< endl;
		std::cerr << "AUG.iBaseHaz[th][ActiveMild] th" << th << " " << AUG.iBaseHaz[th][ActiveMild] << " AUG.iBaseHaz[th][PassiveSevere] " << AUG.iBaseHaz[th][PassiveSevere] << endl;
		std::cerr << "AUG.iVacHaz[th][ActiveMild][SeroPos] th" << th << " " << AUG.iVacHaz[th][ActiveMild][SeroPos] << " AUG.iVacHaz[th][ActiveMild][SeroNeg]  " << AUG.iVacHaz[th][ActiveMild][SeroNeg] << endl;
		std::cerr << "AUG.iVacHaz[th][PassiveSevere][SeroPos] th" << th << " " << AUG.iVacHaz[th][PassiveSevere][SeroPos] << " AUG.iVacHaz[th][PassiveSevere][SeroNeg]  " << AUG.iVacHaz[th][PassiveSevere][SeroNeg] << endl;
		std::cerr << "AUG.iVacHaz_Mult[th][ActiveMild][SeroPos] th" << th << " " << AUG.iVacHaz_Mult[th][ActiveMild][SeroPos] << " AUG.iVacHaz_Mult[th][ActiveMild][SeroNeg]  " << AUG.iVacHaz_Mult[th][ActiveMild][SeroNeg] << endl;
		std::cerr << "AUG.iVacHaz_Mult[th][PassiveSevere][SeroPos] th" << th << " " << AUG.iVacHaz_Mult[th][PassiveSevere][SeroPos] << " AUG.iVacHaz_Mult[th][PassiveSevere][SeroNeg]  " << AUG.iVacHaz_Mult[th][PassiveSevere][SeroNeg] << endl;
		std::cerr << "AUG.iBaseHaz_Mult[th][ActiveMild][SeroPos] th" << th << " " << AUG.iBaseHaz_Mult[th][ActiveMild][SeroPos] << " AUG.iBaseHaz_Mult[th][ActiveMild][SeroNeg]  " << AUG.iBaseHaz_Mult[th][ActiveMild][SeroNeg] << endl;
		std::cerr << "AUG.iBaseHaz_Mult[th][PassiveSevere][SeroPos] th" << th << " " << AUG.iBaseHaz_Mult[th][PassiveSevere][SeroPos] << " AUG.iBaseHaz_Mult[th][PassiveSevere][SeroNeg]  " << AUG.iBaseHaz_Mult[th][PassiveSevere][SeroNeg] << endl;
		std::cerr << "AUG.HazCaseMult[th][SeroPos] th" << th << " " << AUG.HazCaseMult[th][SeroPos] << " AUG.HazCaseMult[th][SeroNeg] " << AUG.HazCaseMult[th][SeroNeg] << endl;
		std::cerr << "PARAMS.Eff_SPos th" << th << " " << PARAMS.Efficacies[DATA.Case_PhaseSeverity[person_i]][Case_Eff_serotype][SeroPos] << " PARAMS.Eff_SNeg " << PARAMS.Efficacies[DATA.Case_PhaseSeverity[person_i]][Case_Eff_serotype][SeroNeg] << endl;
		std::cerr << "PARAMS.AgeEff_Mult_SPos th" << th << " " << PARAMS.AgeEff_Mult[SeroPos][DATA.ai_s[person_i]] << " PARAMS.AgeEff_Mult_SNeg " << PARAMS.AgeEff_Mult[SeroNeg][DATA.ai_s[person_i]] << endl;
		std::cerr << "P_Spos " << AUG.P_BS[th][SeroPos] << " P_Sneg " << AUG.P_BS[th][SeroNeg] << endl;
		std::cerr << "age " << DATA.ai_s[person_i] << endl;
		std::cerr << "BaseHazValue " << BaseHazValue << " Waning_SPos " << PARAMS.WaningMults[DATA.ai_s[person_i]][SeroPos][DATA.TimePost_Final_Dose[person_i]] << " Waning_SNeg " << PARAMS.WaningMults[DATA.ai_s[person_i]][SeroNeg][DATA.TimePost_Final_Dose[person_i]] << endl;
		std::cerr << "K_SPos_case " << AUG.Case_Ks[th][SeroPos] << " K_SNeg_case " << AUG.Case_Ks[th][SeroNeg] << endl;
	
		if (HOUSE.SeroSpecificEfficacies || HOUSE.SeroSpecific_K_values) if (DATA.IsCase[person_i])
		{
			std::cerr << " CaseSerotype "	<< DATA.CaseSerotype[person_i] << endl;
			std::cerr << " rho "			<< PARAMS.rhos[country][DATA.CaseSerotype[person_i]] << endl;
		}
		std::fflush(stderr);
	}*/

	return (P_Spos_Given_Outcome); 
}
void AugmentData					(Augment_Struct &AUG, DATA_struct &DATA, Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE, Chains_Struct &CHAINS)
{
	Calc_SumRhoEffs_ALL		(CurrentPARAMS	, ProposedPARAMS, HOUSE);		/// must be done before the threading. 
	Calc_SumRhoKs_ALL		(CurrentPARAMS	, ProposedPARAMS, HOUSE);		/// must be done before the threading. 
	Calc_SumRhoEffKs_ALL	(CurrentPARAMS	, ProposedPARAMS, HOUSE);		/// must be done before the threading. 

	//// reset Aug Acceptance array for this iteration
	for (int thread_no = 0; thread_no < HOUSE.max_threads; thread_no++)		CHAINS.AugAccArray_SingleIter[thread_no] = 0; 
	//// loop over augmented data.
#pragma omp parallel for schedule(static,1) 
	for (int thread_no = 0; thread_no < HOUSE.max_threads; thread_no++)
	{
		int ProposedAugmentedValue;
		DType AcceptOrReject_RandomnNumber_Aug, logAcceptanceProb_Aug, ProposedVsOldlogPosteriorRatio_Aug, r_Aug = 0, AugRandomNumber = 0, ProbSposGivenOutcome = NULL, GibbsAugRandomNumber = NULL /*, Proposed_PriorComponent_Aug, Proposed_logPrior_Aug, Proposed_logPosterior_Aug*/;

		for (int patient = thread_no; patient < DATA.NoAugmented; patient += HOUSE.max_threads)
		{
			if (AUG.Aug_AllCool == 0) break;

			Adjust_Aug_thread(thread_no, DATA.AugmentedIndices[patient], DATA, CurrentPARAMS, HOUSE, AUG);

			///// Choose/Propose new value for augmented data patient (choose in case of Gibbs sampling, propose in case of MCMC). 
			if (HOUSE.Aug_MH_or_Gibbs == GIBBS_AUG)
			{
				ProbSposGivenOutcome	= ProbSerostatusGivenOutcome(DATA.AugmentedIndices[patient], DATA, CurrentPARAMS, HOUSE, AUG, thread_no);
				GibbsAugRandomNumber	= RandomNumber_Mt(thread_no); 
				ProposedAugmentedValue	= (GibbsAugRandomNumber < ProbSposGivenOutcome) ? 1 : 0;
			}
			else if (HOUSE.Aug_MH_or_Gibbs == MH_AUG)
			{
				AugRandomNumber			= RandomNumber_Mt(thread_no);
				ProposedAugmentedValue	= (AugRandomNumber > 0.5) ? 1 : 0;
			}
			else std::cerr << "HOUSE.AugMCMC_or_Gibbs not recognized" << endl;

			if (AUG.Aug_AllCool == 0) patient = DATA.NoAugmented + 1;

			if ((ProposedAugmentedValue != DATA.Ii_s[DATA.AugmentedIndices[patient]]))  // i.e. if not proposing the same augmented value as is currently stored. 
			{
				Change_Aug_Data_patient(DATA.AugmentedIndices[patient], patient, DATA, CurrentPARAMS, HOUSE, AUG, thread_no); ///// for gibbs sampling, do not change value of Ii until after this is called, as multiplier of differences (i.e. add or subtract from various components) depends on original value. Must change, without accept reject step, afterwards. 

				if (HOUSE.Aug_MH_or_Gibbs == MH_AUG)
				{
					////// Calculate posterior prob of proposed parameter. 
					ProposedVsOldlogPosteriorRatio_Aug	= AUG.LikeDiffs[patient][AUG.NoLikeDiffsPerAugPatient] - CurrentPARAMS.LikeFull;
					logAcceptanceProb_Aug				= ProposedVsOldlogPosteriorRatio_Aug / CHAINS.Temp_ThisIteration;
					r_Aug								= RandomNumber_Mt(thread_no);
					AcceptOrReject_RandomnNumber_Aug = log(r_Aug);

					if (AcceptOrReject_RandomnNumber_Aug < logAcceptanceProb_Aug) // i.e. if you accept...
					{
						(CHAINS.AcceptArray_Aug[thread_no])++;
						(CHAINS.AugAccArray_SingleIter[thread_no])++; 
						DATA.Ii_s[DATA.AugmentedIndices[patient]] = ProposedAugmentedValue;

						//// Reset the integrated vaccine hazard for that patient (for both Current & Proposed Params)
						for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
						{
							CurrentPARAMS. IVH_vals[DATA.AugmentedIndices[patient]][PhaseSeverity] = AUG.iVacHaz[thread_no][PhaseSeverity][DATA.Ii_s[DATA.AugmentedIndices[patient]]];
							ProposedPARAMS.IVH_vals[DATA.AugmentedIndices[patient]][PhaseSeverity] = AUG.iVacHaz[thread_no][PhaseSeverity][DATA.Ii_s[DATA.AugmentedIndices[patient]]];
						}
					}
					else for (int diff = 0; diff < AUG.NoLikeDiffsPerAugPatient + 1; diff++) AUG.LikeDiffs[patient][diff] = 0;  //// need to reset differences if you reject
				}
				else if (HOUSE.Aug_MH_or_Gibbs == GIBBS_AUG)
				{
					DATA.Ii_s[DATA.AugmentedIndices[patient]] = ProposedAugmentedValue;	
					(CHAINS.AcceptArray_Aug[thread_no])++;								
					(CHAINS.AugAccArray_SingleIter[thread_no])++;

					//// Reset the integrated vaccine hazard for that patient (for both Current & Proposed Params)
					for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
					{
						CurrentPARAMS. IVH_vals[DATA.AugmentedIndices[patient]][PhaseSeverity] = AUG.iVacHaz[thread_no][PhaseSeverity][DATA.Ii_s[DATA.AugmentedIndices[patient]]];
						ProposedPARAMS.IVH_vals[DATA.AugmentedIndices[patient]][PhaseSeverity] = AUG.iVacHaz[thread_no][PhaseSeverity][DATA.Ii_s[DATA.AugmentedIndices[patient]]];
					}
				}
				else std::cerr << "HOUSE.AugMCMC_or_Gibbs not recognized" << endl;
			}
		}
	}

	for (int patient = 0; patient < DATA.NoAugmented; patient++)
		for (int component = 0; component < HOUSE.LCPerC; component++)
		{
			CurrentPARAMS.LikeParts[DATA.ci_s[DATA.AugmentedIndices[patient]]][component] += AUG.LikeDiffs[patient][component];
			AUG.LikeDiffs[patient][component] = 0; // reset difference to zero
		} 

	CurrentPARAMS.LikeFull	= l_full(DATA, CurrentPARAMS, HOUSE);
	for (int country = 0; country < HOUSE.TotalCountries; country++)
		for (int component = 0; component < HOUSE.LCPerC; component++)
			ProposedPARAMS.LikeParts[country][component] = CurrentPARAMS.LikeParts[country][component]; //// reset proposed likelihood too for future iterations
	ProposedPARAMS.LikeFull = CurrentPARAMS.LikeFull;

	if (HOUSE.Weighting_Pass_Sev)		CalcWeightedLikes(DATA, CurrentPARAMS, ProposedPARAMS, HOUSE);

	//// Must clear and recalculate sets if augmenting.  
	ClearSets	(DATA, HOUSE);		
	PopulateSets(DATA, HOUSE);		
	for (int countryindex = 0; countryindex < HOUSE.NoCountriesToFit; countryindex++)
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			{ 
				CalcAndStoreSumIntVacHazOverPatients(BaselineSeroStatus, HOUSE.WhichCountries[countryindex], DATA, CurrentPARAMS, HOUSE, PhaseSeverity);																		//// recalculate for CurrentPARAMS. 
				ProposedPARAMS.VacHazLikes[HOUSE.WhichCountries[countryindex]][BaselineSeroStatus][PhaseSeverity] = 
				CurrentPARAMS. VacHazLikes[HOUSE.WhichCountries[countryindex]][BaselineSeroStatus][PhaseSeverity];	//// copy to ProposedPARAMS
			}
	Calc_SumRhoEffKs_ALL(CurrentPARAMS , HOUSE);		
	Calc_SumRhoEffKs_ALL(ProposedPARAMS, HOUSE);		
	if (CHAINS.AreWeChecking)		// Check likelihood is the same whether you calculate efficiently or from scratch 
		CHAINS.LikeAndLikeAlike = CompareLikelihoods("Augmentation", DATA, CurrentPARAMS, HOUSE, CHAINS.LikePartsForChecking, CHAINS.VacHazLike_CheckingArray);
}
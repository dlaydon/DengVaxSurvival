
#include "StructureDefs.h"
#include "ParamNumber.h"
#include "ParamUpdate.h"
#include "Likelihood.h"
#include "Probability.h"
#include "ConsolePrinting.h"
#include "CalcParamsEtc.h"
#include "Augmentation.h"
#include "SurvivalCurves.h"

std::vector<int>	FindSimulataneousUpdateParams	(int &param_no, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE) 
{
	//// function returns vector of size 1 with entry param_no if not doing simultaneous update, i.e. vast majority of the time. 
	
	std::vector<int> SimulataneousUpdateParams; 

	if (IsParamAnEfficacy(param_no, PARAMS.ParamNumbers))
	{
		int BaselineSeroStatus = Find_SeroStatus_FromEffParam(param_no, HOUSE, PARAMS.ParamNumbers);

		if ((HOUSE.Single_SNeg_Eff && BaselineSeroStatus == SeroNeg) || ((HOUSE.Single_SPos_Eff && BaselineSeroStatus == SeroPos)))
		{
			for (int VE_param = PARAMS.ParamNumbers.Min_VacE; VE_param <= PARAMS.ParamNumbers.Max_VacE; VE_param++) 
				if (VE_param % HOUSE.HowManySeroStatuses == BaselineSeroStatus) SimulataneousUpdateParams.push_back(VE_param);
		}
		else if (HOUSE.SingleEff)
		{
			for (int VE_param = PARAMS.ParamNumbers.Min_VacE; VE_param <= PARAMS.ParamNumbers.Max_VacE; VE_param++) 
				SimulataneousUpdateParams.push_back(VE_param);
		}	
		else SimulataneousUpdateParams.push_back(param_no);
	}
	else if (IsParamARelativeRisk(param_no, PARAMS.ParamNumbers) || IsParamARelativeRisk_Hosp(param_no, PARAMS.ParamNumbers))
	{
		int PrevInf_of_ParamNo			= Find_PrevInf_From_K_Param			(param_no, HOUSE, PARAMS.ParamNumbers);
		int PhaseSeverity_of_ParamNo	= Find_PhaseSeverity_From_K_Param	(param_no, HOUSE, PARAMS.ParamNumbers, PrevInf_of_ParamNo);

		if (HOUSE.RelRisksSameBtwSerotypes)
		{
			int PrevInf, PhaseSeverity;
			///// cycle through all K parameters, add ones corresponding to PrevInf to SimulataneousUpdateParams
			for (int K_Param = PARAMS.ParamNumbers.Min_K; K_Param <= PARAMS.ParamNumbers.Max_K; K_Param++)
			{
				PrevInf			= Find_PrevInf_From_K_Param			(K_Param, HOUSE, PARAMS.ParamNumbers);
				PhaseSeverity	= Find_PhaseSeverity_From_K_Param	(K_Param, HOUSE, PARAMS.ParamNumbers, PrevInf);

				if (PrevInf == PrevInf_of_ParamNo && PhaseSeverity == PhaseSeverity_of_ParamNo) SimulataneousUpdateParams.push_back(K_Param);
			}
		}
		else if (HOUSE.Fixed_Severe_RelRisks)
		{
			if (PrevInf_of_ParamNo == 1 && PhaseSeverity_of_ParamNo == PassiveSevere) //// if param is KPS_1 of any serotype.
			{
				SimulataneousUpdateParams.push_back(param_no - 1);  //// i.e. KPS_0 of any serotype
				SimulataneousUpdateParams.push_back(param_no);		//// i.e. KPS_1 of any serotype
				SimulataneousUpdateParams.push_back(param_no + 1);	//// i.e. KPS_2 of any serotype
			}
			else SimulataneousUpdateParams.push_back(param_no);
		}
		else SimulataneousUpdateParams.push_back(param_no);
	}
	else if (IsParamAWaningParam(param_no, PARAMS.ParamNumbers))
	{
		if (HOUSE.AgeEffectsSame_Waning)
		{
			int BaselineSeroStatus_OrigParam = Find_SeroStatus_From_Waning_Param(param_no, HOUSE, PARAMS.ParamNumbers);
			for (int WaningParam_param = PARAMS.ParamNumbers.MinWaningParam; WaningParam_param <= PARAMS.ParamNumbers.MaxWaningParam; WaningParam_param++) 
				if (Find_SeroStatus_From_Waning_Param(WaningParam_param, HOUSE, PARAMS.ParamNumbers) == BaselineSeroStatus_OrigParam) SimulataneousUpdateParams.push_back(WaningParam_param);
		}	
		else SimulataneousUpdateParams.push_back(param_no);
	}
	else if (IsParamAn_ASVE_Param(param_no, PARAMS.ParamNumbers))
	{
		if (HOUSE.AgeEffectsSame_VE)
		{
			int BaselineSeroStatus_OrigParam = Find_SeroStatus_FromASVE_Param(param_no, HOUSE, PARAMS.ParamNumbers);
			for (int ASVE_param = PARAMS.ParamNumbers.Min_ASVE_Param; ASVE_param <= PARAMS.ParamNumbers.Max_ASVE_Param; ASVE_param++) //// note <= equality sign. 
				if (Find_SeroStatus_FromASVE_Param(ASVE_param, HOUSE, PARAMS.ParamNumbers) == BaselineSeroStatus_OrigParam) SimulataneousUpdateParams.push_back(ASVE_param);
		}	
		else SimulataneousUpdateParams.push_back(param_no);
	}
	else		SimulataneousUpdateParams.push_back(param_no);

	return SimulataneousUpdateParams; 
}
std::vector<DType>	Populate_Proposed_Param_Vec		(std::vector<int> &SimulataneousUpdateParams, int &param_no, DType Proposed_Param, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	std::vector<DType> Proposed_Param_Vector(SimulataneousUpdateParams.size(), Proposed_Param); 
	
	//// Amend Proposed_Param_Vector as necessary. 
	if (HOUSE.Fixed_Severe_RelRisks)
		if (IsParamARelativeRisk(param_no, PARAMS.ParamNumbers) || IsParamARelativeRisk_Hosp(param_no, PARAMS.ParamNumbers))
		{
			int PrevInf_of_ParamNo			= Find_PrevInf_From_K_Param			(param_no, HOUSE, PARAMS.ParamNumbers);
			int PhaseSeverity_of_ParamNo	= Find_PhaseSeverity_From_K_Param	(param_no, HOUSE, PARAMS.ParamNumbers, PrevInf_of_ParamNo);

			if (PrevInf_of_ParamNo == 1 && PhaseSeverity_of_ParamNo == PassiveSevere) //// if param is KPS_1 of any serotype.
				for (int sim_param = 0; sim_param < SimulataneousUpdateParams.size(); sim_param++)
					Proposed_Param_Vector[sim_param] *= PARAMS.Fixed_SevereK_ratios[sim_param];
		}
	return Proposed_Param_Vector;
}
std::vector<DType>	Calc_RowMeans					(DType ** ParamChain, int NumRows, int NumCols) 
{
	std::vector<DType> MeanPostParamVec(NumRows, (DType)0); // initialize to zero

	for (int row = 0; row < NumRows; row++)
	{
		for (int col = 0; col < NumCols; col++)
			MeanPostParamVec[row] += ParamChain[row][col]; 
		MeanPostParamVec[row] /= NumCols;
	}
	return MeanPostParamVec; 
}
DType				ChooseTemp						(int &iteration, int &BurnIn, int &No_Iterations, bool CoolDuringBurnIn, bool CoolAfterBurnIn, DType FinalTemperature, /*DType ExpRate, DType Slope, */bool ExponentialTrue_LinearFalse)
{
	/*
	TempAfterBurnIn = set to 1 if you want to run a regular chain after (hopefully) climbing to region of high probability mass. Set to w_n = 1 / log (NPat) (approx 1/10) for WBIC. Set to something smaller for colder chain (e.g. 1/100)
	
	*/

	DType Temp = 1; 
	if (iteration < BurnIn)
	{
		if (CoolDuringBurnIn)
		{
			if (ExponentialTrue_LinearFalse)
			{
				DType ExpRate	= log(FinalTemperature) / BurnIn; // leave out minus signs here and in Temp calculation as they cancel
				Temp			= exp(ExpRate * iteration);
			}
			else
			{
				DType Slope		= (FinalTemperature - 1) / BurnIn; //// i.e. delta(y) / delta (x)
				Temp			= (Slope * iteration) + 1;
			}
		}
	}
	else
	{
		if (CoolAfterBurnIn)
		{
			if (ExponentialTrue_LinearFalse)
			{
				DType ExpRate	= log(FinalTemperature) /  (No_Iterations - BurnIn); // leave out minus signs here and in Temp calculation as they cancel
				Temp			= exp(ExpRate * (iteration - BurnIn));
			}
			else
			{
				DType Slope		= (FinalTemperature - 1) / (No_Iterations - BurnIn); //// i.e. delta(y) / delta (x)
				Temp			= (Slope * (iteration - BurnIn)) + 1; //// assuming y-intercept always 1. 
			}
		}
		else Temp = FinalTemperature;
	}
	return Temp; 
}
DType				ProposeNewParam					(int &param_no, int &thread_num, const Params_Struct &CurrentPARAMS, const Housekeeping_Struct &HOUSE)
{
#ifdef PRINT_WITHIN_CHANGE_FUNCTIONS
#ifdef PRINT_WITHIN_CHANGE_KNOT_FUNCTIONS
	std::cout << " ProposeNewParam p" << param_no << " " << CurrentPARAMS.NamesOfParameters[param_no] << endl;
#else 
	if (!IsParamAKnot(param_no, CurrentPARAMS.ParamNumbers)) std::cout << " ProposeNewParam p" << param_no << " " << CurrentPARAMS.NamesOfParameters[param_no] << endl;
#endif 
#endif

	DType Proposed_Param = CurrentPARAMS.ParamVec[param_no];	//// set Proposed equal to current

	if (!HOUSE.LinKnts) if (IsParamAKnot(param_no, CurrentPARAMS.ParamNumbers)) Proposed_Param = log10(Proposed_Param); // take log of knots if desired. 
	if (HOUSE.FitWaningRate && IsParamAWaningParam(param_no, CurrentPARAMS.ParamNumbers))
		if (Find_ParamType_From_ASWaning_Param(param_no, HOUSE, CurrentPARAMS.ParamNumbers) == WaningRateDuration) 
			Proposed_Param = 1 / Proposed_Param;

#ifdef USE_RANDOM_NUMBERS
	Proposed_Param = Proposed_Param + (CurrentPARAMS.ProposalStandDevs[param_no] * snorm_mt(thread_num)); //// choose new param based on normal propoosal and standard devs. note for knots, proposals still normal, just on log scale. Check your notes on this. Essentially you've got a metric transform (of 1/x if natural log, 1/(ln(10) * x) if base 10) on your uniform prior and normal proposal. These values cancel hence your acceptance probabilities are calculated in the same way as below. Note however than you are proposing values on a log normal distribution, which would normally require a correction factor. The fact that you don't have one is due to there being a prior which cancels out said correction factor. 
#else
	Proposed_Param = Proposed_Param + (CurrentPARAMS.ProposalStandDevs[param_no] * 1); //// choose new param based on normal propoosal and standard devs. note for knots, proposals still normal, just on log scale. Check your notes on this. Essentially you've got a metric transform (of 1/x if natural log, 1/(ln(10) * x) if base 10) on your uniform prior and normal proposal. These values cancel hence your acceptance probabilities are calculated in the same way as below. Note however than you are proposing values on a log normal distribution, which would normally require a correction factor. The fact that you don't have one is due to there being a prior which cancels out said correction factor. 
#endif
	//// "wrapping method"  
	while ((Proposed_Param < CurrentPARAMS.ParamRanges[LowerBound][param_no]) || (Proposed_Param > CurrentPARAMS.ParamRanges[UpperBound][param_no])) // i.e. if proposed parameter is either lower than lower bound or higher than higher bound, choose a new parameter. 
	{
		if ((Proposed_Param < CurrentPARAMS.ParamRanges[LowerBound][param_no])) Proposed_Param = CurrentPARAMS.ParamRanges[UpperBound][param_no] - (CurrentPARAMS.ParamRanges[LowerBound][param_no] - Proposed_Param);  // 
		if ((Proposed_Param > CurrentPARAMS.ParamRanges[UpperBound][param_no])) Proposed_Param = CurrentPARAMS.ParamRanges[LowerBound][param_no] + (Proposed_Param - CurrentPARAMS.ParamRanges[UpperBound][param_no]);
	}
	if (!HOUSE.LinKnts) if (IsParamAKnot(param_no, CurrentPARAMS.ParamNumbers)) Proposed_Param = pow(10, Proposed_Param); // if fitting log of knots, if param is a knot, raise 10 to power Proposed_Param. 
	if (HOUSE.FitWaningRate && IsParamAWaningParam(param_no, CurrentPARAMS.ParamNumbers))
		if (Find_ParamType_From_ASWaning_Param(param_no, HOUSE, CurrentPARAMS.ParamNumbers) == WaningRateDuration)
			Proposed_Param = 1 / Proposed_Param;

	return Proposed_Param; 
}

bool MCMC_step_param		(int param_no, const DATA_struct &DATA, Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE, Chains_Struct &CHAINS)
{
	bool All_Good = true; 

	/*Checks*/	if (!Everything_Okay(param_no, DATA, CurrentPARAMS, ProposedPARAMS, HOUSE, CurrentPARAMS.NamesOfParameters[param_no])) All_Good = false;
	CHAINS.Proposed_Param				= ProposeNewParam(param_no, CHAINS.thread_num, CurrentPARAMS, HOUSE);
	CHAINS.SimulataneousUpdateParams	= FindSimulataneousUpdateParams(param_no, ProposedPARAMS, HOUSE); //// assumes that parameters updated simultaneously will have same value. 
	CHAINS.Proposed_Param_Vector		= Populate_Proposed_Param_Vec(CHAINS.SimulataneousUpdateParams, param_no, CHAINS.Proposed_Param, ProposedPARAMS, HOUSE);

	for (int sim_update_param = 0; sim_update_param < CHAINS.SimulataneousUpdateParams.size(); sim_update_param++)
	{
		AmendParams		(CHAINS.SimulataneousUpdateParams[sim_update_param], CHAINS.Proposed_Param_Vector[sim_update_param], ProposedPARAMS, DATA, HOUSE);

		if (HOUSE.ASVE != Age_Option::INDEPENDENT & HOUSE.SeroSpecificEfficacies) //// Ensure efficacy and multiplier combination does not result in >100% efficacy for all ages, serostatuses and serotypes. 
			if (		IsParamAnEfficacy	(CHAINS.SimulataneousUpdateParams[sim_update_param], ProposedPARAMS.ParamNumbers)	| 
						IsParamAn_ASVE_Param(param_no, ProposedPARAMS.ParamNumbers)												|
						IsParamA_rho		(param_no, ProposedPARAMS.ParamNumbers)												|
						IsParamA_qval		(param_no, ProposedPARAMS.ParamNumbers)												)
				for (int BS = 0; BS < HOUSE.HowManySeroStatuses; BS++)
					for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
						for (int age = 0; age < HOUSE.HowManyAges; age++)
						{
							DType UpperBound_TI = (BS == SeroNeg) ? HOUSE.IntialParRanges.SNegEff_1[UpperBound] : HOUSE.IntialParRanges.SPosEff_1[UpperBound];
							DType LowerBound_TI = (BS == SeroNeg) ? HOUSE.IntialParRanges.SNegEff_1[LowerBound] : HOUSE.IntialParRanges.SPosEff_1[LowerBound];
							
							// Aggregate Fixed-time immunity
							DType Agg_TI = HOUSE.AdditiveSSASVEs ? 
											ProposedPARAMS.Efficacies[ActiveMild][serotype][BS] + ProposedPARAMS.AgeEff_Mult[BS][age]	: 
											ProposedPARAMS.Efficacies[ActiveMild][serotype][BS] * ProposedPARAMS.AgeEff_Mult[BS][age]	;
							
							if (Agg_TI > UpperBound_TI || Agg_TI < LowerBound_TI)
							{
								CHAINS.logAcceptanceProb = log(0);
								goto AcceptOrReject_LABEL;
							}
						}
		if (CHAINS.CalculateFreshLikeLihood)	LikelihoodFromScratch	(DATA, ProposedPARAMS, HOUSE, ProposedPARAMS.LikeParts, ProposedPARAMS.VacHazLikes);
		else									Change_Param			(CHAINS.SimulataneousUpdateParams[sim_update_param], CHAINS.Proposed_Param, DATA, CurrentPARAMS, ProposedPARAMS, HOUSE);
	}
	ProposedPARAMS.LikeFull			= l_full(DATA, ProposedPARAMS, HOUSE);
	ProposedPARAMS.LogPosterior		= ProposedPARAMS.LikeFull + ProposedPARAMS.LogPriorFull; 
	if (HOUSE.Weighting_Pass_Sev)	CalcWeightedLikes(DATA, CurrentPARAMS, ProposedPARAMS, HOUSE);
	if (!Likelihood_Okay(CurrentPARAMS, ProposedPARAMS)) All_Good = false;
	if (All_Good == false) PrintFittedParamNumbers(HOUSE, CurrentPARAMS);
	if (CHAINS.AreWeChecking)	// Checks that likelihood is the same whether calculated efficiently or from scratch each time. 
	{
		CHAINS.LikeAndLikeAlike = CompareLikelihoods("Parameter", DATA, ProposedPARAMS, HOUSE, CHAINS.LikePartsForChecking, CHAINS.VacHazLike_CheckingArray);
		if (CHAINS.LikeAndLikeAlike == 0)	All_Good = false; 
	}

	// Calculate posterior prob of proposal 
	CHAINS.ProposedVsOldlogPosteriorRatio	= *ProposedPARAMS.ptr_L_Full - *CurrentPARAMS.ptr_L_Full; 
	CHAINS.logAcceptanceProb				= CHAINS.ProposedVsOldlogPosteriorRatio / CHAINS.Temp_ThisIteration;
	if (CHAINS.logAcceptanceProb > 0) CHAINS.logAcceptanceProb = (DType)0; 

AcceptOrReject_LABEL:
	CHAINS.randomNo						= RandomNumber_Mt(CHAINS.thread_num);
	CHAINS.AcceptOrReject_RandomnNumber	= log(CHAINS.randomNo);
	CHAINS.ParamAcceptCondition			= CHAINS.AcceptOrReject_RandomnNumber < CHAINS.logAcceptanceProb;

	MaybePrint(param_no, CHAINS.Proposed_Param, CurrentPARAMS, ProposedPARAMS, HOUSE, CHAINS.ParamAcceptCondition); 

	if (CHAINS.ParamAcceptCondition)  
	{
		++(CHAINS.AcceptArray[param_no]);
		for (int sim_update_param = 0; sim_update_param < CHAINS.SimulataneousUpdateParams.size(); sim_update_param++)
			UpdateOrReset_ParamsEtc(CHAINS.SimulataneousUpdateParams[sim_update_param], DATA, CurrentPARAMS, ProposedPARAMS, HOUSE);

		for (int sim_update_param = 0; sim_update_param < CHAINS.SimulataneousUpdateParams.size(); sim_update_param++)
		{
			CHAINS.ProportionUpdatesAccepted++;
			CHAINS.AverageAcceptanceProb += exp(CHAINS.logAcceptanceProb); 
		}
	}
	else for (int sim_update_param = 0; sim_update_param < CHAINS.SimulataneousUpdateParams.size(); sim_update_param++)
			UpdateOrReset_ParamsEtc(CHAINS.SimulataneousUpdateParams[sim_update_param], DATA, ProposedPARAMS, CurrentPARAMS, HOUSE);	/*if you reject (i.e. same function with arguments reversed)*/

	return All_Good; 
}
void AddTo_LL_SPrev_Chains	(int PostSampleNum, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, Chains_Struct &CHAINS, DType *** Which_LL_SPrevChain, bool TotalPopulation_True_ImSubOnly_False)
{
#ifdef PRINT_PROGRAM_PROGRESS
	std::cerr << "AddTo_LL_SPrev_Chains " << endl;
#endif
	// Which_LL_SPrevChain is either CHAINS.LL_SPrev_Chains or CHAINS.LL_SPrev_Chains_ImSub
	// Similar in structure to l_full_Minus_Aug function in Likelihood.cpp. But this outputs likelihood components associated with serostatus (either sum(log(exp(-ha))) or sum(log(1-exp(-ha))) )
	DType LL_Value = 0; 
	for (int country = 0; country < HOUSE.TotalCountries; country++)
	{
		std::vector<int> CountriesToAddTo; 
		CountriesToAddTo.push_back(country);													// country 
		if (country < 5) CountriesToAddTo.push_back(10); else CountriesToAddTo.push_back(11);	// trial
		CountriesToAddTo.push_back(12);															// combined CYD114 and CYD15 pooled trials

		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
		{
			 //calculate or retreive LL component
			if (TotalPopulation_True_ImSubOnly_False)	LL_Value = PARAMS.LikeParts[country][HOUSE.L_Indices.BSs[BaselineSeroStatus]];
			else										LL_Value = l_Aug(country, BaselineSeroStatus, DATA, PARAMS, HOUSE, true);  // i.e. ImSubOnly = true; 
			// add LL component to appropriate "countries" (for that Posterior Sample). 
			for (int countryindex = 0; countryindex < CountriesToAddTo.size(); countryindex++)
				Which_LL_SPrevChain[PostSampleNum][CountriesToAddTo[countryindex]][BaselineSeroStatus] += LL_Value;
		}
	}
}
void AddToChains			(int iter, Chains_Struct &CHAINS, const DATA_struct &DATA, Params_Struct &CurrentPARAMS, const Housekeeping_Struct &HOUSE, Survival_Struct &SURVIVE, Survival_Struct &SURVIVE_ImSub, Augment_Struct &AUG)
{
	//// Parameter chains and Seroprevalence. Wrapper function that adds to chains and calculates SeroPrevs
	if (((iter + 1) % CHAINS.AddtoChainEveryHowManyIterations == 0) && (iter >= CHAINS.BurnIn))
	{
#ifdef PRINT_PROGRAM_PROGRESS
		std::cerr << "AddTo_Param_SeroPrev_SCurve_Chains " << endl;
#endif
		int PostSampleNum = (iter - CHAINS.BurnIn) / CHAINS.AddtoChainEveryHowManyIterations;

		CurrentPARAMS.Total_logLike = LL_Total(DATA, CurrentPARAMS, HOUSE, AUG);

		if (PostSampleNum < CHAINS.NumPosteriorSamples) // to ensure no memory access violations
		{
			if (CHAINS.AreWeCalculatingSeroPrevs)
				CHAINS.SEROPREV.Calc_Aug_AgeSpecificSeroPrev(DATA, HOUSE, PostSampleNum); // SeroPrevs

			// log likelihood chain(s). 
			CHAINS.logLikeChain			[PostSampleNum]		= CurrentPARAMS.LikeFull;		// log Likelihood 
			CHAINS.log_PosteriorChain	[PostSampleNum]		= CurrentPARAMS.LogPosterior;	// log Posterior 
			CHAINS.Total_Likelihood_Chain[PostSampleNum]	= CurrentPARAMS.Total_logLike;	// Total Likelihood that doesn't assume serostatus of non-immune subset. 

			if (HOUSE.OutputSeroPrev_LLChains)
			{
				AddTo_LL_SPrev_Chains(PostSampleNum, DATA, CurrentPARAMS, HOUSE, CHAINS, CHAINS.LL_SPrev_Chains			, true);
				AddTo_LL_SPrev_Chains(PostSampleNum, DATA, CurrentPARAMS, HOUSE, CHAINS, CHAINS.LL_SPrev_Chains_ImSub	, false);
			}

			//// Parameter chain:
			for (int param = 0; param < HOUSE.No_Parameters; param++)
				CHAINS.ParamChain[param][PostSampleNum] = CurrentPARAMS.ParamVec[param];

			/// Add to rho chains. 
			if (HOUSE.SeroSpecificEfficacies || HOUSE.SeroSpecific_K_values)
				if (!HOUSE.BaselinePartition) //// because rhos already outputted as parametes in Baseline partition 
					if (CHAINS.NumPosteriorSamples > 0) 
						for (int country = 0; country < HOUSE.TotalCountries; country++)
							for (int serotype = 0; serotype < HOUSE.N_STypes; serotype++)
								CHAINS.rhos[country][serotype][PostSampleNum] = CurrentPARAMS.rhos[country][serotype];

			//// log likelihood (minus augmentation parts) chain if doing. 
			if (CHAINS.Output_LL_minus_Aug)		CHAINS.logLike_minus_Aug_Chain[PostSampleNum] = l_full_Minus_Aug(DATA, CurrentPARAMS, HOUSE);

			//// add to augmented data mean posterior. 
			for (int patient = 0; patient < NPat; patient++) CHAINS.MeanPost_Ii_s[patient] += DATA.Ii_s[patient];

			//// update current max likelihood. 
			if (CurrentPARAMS.LikeFull > CHAINS.MaxLikeSoFar)
			{
				//// reset MaxLikeSoFar
				CHAINS.MaxLikeSoFar				= CurrentPARAMS.LikeFull;
				//// reset Reset MaxLike ParamVec
				CHAINS.MaxLike_ParamVec			= CurrentPARAMS.ParamVec;
				//// reset MaxLike Augmented data
				for (int subject = 0; subject < NPat; subject++)	CHAINS.MaxLike_Ii_s[subject] = DATA.Ii_s[subject]; 
			}

			//// update current Modal posterior.
			if (CurrentPARAMS.LogPosterior > CHAINS.ModalPostSoFar)
			{
				//// reset ModalPostSoFar
				CHAINS.ModalPostSoFar			= CurrentPARAMS.LogPosterior;
				//// reset Reset ModalPost ParamVec
				CHAINS.ModalPost_ParamVec		= CurrentPARAMS.ParamVec;
				//// reset ModalPost Augmented data
				for (int subject = 0; subject < NPat; subject++)	CHAINS.ModalPost_Ii_s[subject] = DATA.Ii_s[subject];
			}
		}
	}

	/////////////// /////////////// /////////////// /////////////// /////////////// /////////////// /////////////// /////////////// /////////////// /////////////// /////////////// /////////////// /////////////// /////////////// /////////////// /////////////// 
	/////////////// SURVIVAL Curves / Attack rates/ Hazard Ratios
	if (CHAINS.Calc_SCsARsHRPs_MeanAndCrIs)
		if (((iter + 1) % SURVIVE.AddToSurvivalCurvesEveryHowManyIterations == 0) && ((iter + 1) > CHAINS.BurnIn))
		{
			int SurvivePostSampleNum = (iter - CHAINS.BurnIn) / SURVIVE.AddToSurvivalCurvesEveryHowManyIterations;

			if (SurvivePostSampleNum < SURVIVE.NoSurvivePostSamples) //// to guard against memory access violations
			{
				Calc_SumRhoEffs_ALL(CurrentPARAMS, HOUSE); /// must be done before threading. 

				GenerateSurvivalCurves(DATA, CurrentPARAMS, HOUSE, SURVIVE		, iter, CHAINS, AUG, true	);
				GenerateSurvivalCurves(DATA, CurrentPARAMS, HOUSE, SURVIVE_ImSub, iter, CHAINS, AUG, false	);
			}
		}
}
void MCMC_All				(int FirstIteration, Chains_Struct &CHAINS, DATA_struct &DATA, Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, Augment_Struct &AUG, const Housekeeping_Struct &HOUSE, Survival_Struct &SURVIVE, Survival_Struct &SURVIVE_ImSub)
{
	int param_no = 0; 

	//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** 
	//// **** //// **** //// **** ////			PARAMETERS
	for (int iter = FirstIteration; iter < CHAINS.No_Iterations; iter++) // For regular chains, FirstIteration = 1 not 0 as initial guesses calculated already. For WBIC chains, set to zero. 
	{
		if (CHAINS.ShouldYouStop == 1) break;

		Print_MCMC_Progress(iter, CHAINS, CurrentPARAMS); 
		CHAINS.HowFarDidYouGet = iter; std::fflush(stderr);
	
		if (CHAINS.SimulatedAnnealing)	CHAINS.Temp_ThisIteration = ChooseTemp(iter, CHAINS.BurnIn, CHAINS.No_Iterations, CHAINS.CoolDuringBurnIn, CHAINS.CoolAfterBurnIn, CHAINS.FinalTemperature, CHAINS.ExponentialSchedule); 
		else							CHAINS.Temp_ThisIteration = 1; 

		if (CHAINS.AreWeFittingParameters)
		{
			// loop over parameters
			for (int param_no_dummy = 0; param_no_dummy < CurrentPARAMS.ParamNosYouWillFit.size(); param_no_dummy++)
			{
				param_no = CurrentPARAMS.ParamNosYouWillFit[param_no_dummy];
				bool Continue_MCMC = MCMC_step_param(param_no, DATA, CurrentPARAMS, ProposedPARAMS, HOUSE, CHAINS); 
				if (!Continue_MCMC)
				{
					CHAINS.HowFarDidYouGet = iter; 
					std::cerr << "breaking at p" << param_no << " " << CurrentPARAMS.NamesOfParameters[param_no] << " Cur = " << CurrentPARAMS.ParamVec[param_no] << " Prop = " << ProposedPARAMS.ParamVec[param_no] << endl;
					goto EndMCMC_LABEL;
				}
			}
			if (CHAINS.LikeAndLikeAlike == 0) break;
		}

		// at this point, structures should be equal. 
		if (!Everything_Okay(param_no, DATA, CurrentPARAMS, ProposedPARAMS, HOUSE, "augmentation")) goto EndMCMC_LABEL;

		//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** 
		//// **** //// **** //// **** ////			AUGMENTATION
		if (CHAINS.AreWeAugmenting) if ((iter+1) % CHAINS.AugmentEveryHowManyIterations == 0)
		{
			AugmentData(AUG, DATA, CurrentPARAMS, ProposedPARAMS, HOUSE, CHAINS); 
			if (CHAINS.LikeAndLikeAlike == 0) { CHAINS.HowFarDidYouGet = iter;	break; }
		}	

		if (AUG.Aug_AllCool == 0) goto EndMCMC_LABEL;

		//// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// 
		//// Add to chains / SURVIVAL Curves / Attack rates/ Hazard Ratios
		AddToChains(iter, CHAINS, DATA, CurrentPARAMS, HOUSE, SURVIVE, SURVIVE_ImSub, AUG);

		if (!Likelihood_Okay(CurrentPARAMS, ProposedPARAMS)) goto EndMCMC_LABEL;

		if (CHAINS.OutputParticularChainState)
			if (iter == CHAINS.WhenToFlipTemperatureBack)
			{
				// Params
				std::ofstream ParamStateFile;
				ParamStateFile.open(CHAINS.ParticularChainState_ParamFileName);
				for (int param_number = 0; param_number < HOUSE.No_Parameters; param_number++)	ParamStateFile << CurrentPARAMS.ParamVec[param_number] << endl;
				ParamStateFile.close();

				// Augmented Data
				std::ofstream Ii_sFile;
				Ii_sFile.open(CHAINS.ParticularChainState_AugDataFileName);
				for (int row = 0; row < NPat; row++)	Ii_sFile << DATA.Ii_s[row] << endl;
				Ii_sFile.close();
			}
	}	// END MCMC Loop

EndMCMC_LABEL:

	// Store current state of augmented data (DATA.Iis), parameters and likelihood.
	UpdateOrResetAugData(CHAINS.Final_Iis, DATA.Ii_s, DATA, HOUSE);
	CHAINS.FinalParamVec	= CurrentPARAMS.ParamVec; 
	CHAINS.Final_L_Full		= CurrentPARAMS.LikeFull; 
	SetEqual_2D_Arrays(CHAINS.FinalLikeParts, CurrentPARAMS.LikeParts, HOUSE.TotalCountries, HOUSE.LCPerC); 

	std::cerr << endl << "END MCMC" << endl;
#ifdef PRINT_PROGRAM_PROGRESS
	std::cerr << "SURVIVE.PosteriorSampleCounter " << SURVIVE.PosteriorSampleCounter << std::endl;
	std::cerr << "SURVIVE.NoSurvivePostSamples " << SURVIVE.NoSurvivePostSamples << std::endl;
	std::cerr << "SURVIVE.AddToSurvivalCurvesEveryHowManyIterations " << SURVIVE.AddToSurvivalCurvesEveryHowManyIterations << std::endl;
#endif	
	std::cerr << endl << "Number of iterations completed: " << CHAINS.HowFarDidYouGet + 1 << " out of " << CHAINS.No_Iterations << endl; std::fflush(stderr);
	std::cerr << endl << "CurrentPARAMS.LikeFull = " << CurrentPARAMS.LikeFull << endl << "CurrentPARAMS.Total_logLike = " << CurrentPARAMS.Total_logLike << endl << endl;
}


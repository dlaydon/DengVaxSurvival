
#include "HeaderAndSwitches.h"
#include "ParamNumber.h"

DType RandomNumber_Mt			(int thread_no)
{
	//// wrapper of ranf_mt function - used to clean ifdefs etc out of main code 
	DType RandomNo = 0;
#if defined(USE_RANDOM_NUMBERS)
	RandomNo = ranf_mt(thread_no); /// this will change depending on thread_no (obviously), and thus answer will change with max_threads.
	while (RandomNo == 0) RandomNo = ranf_mt(thread_no);
#else
	RandomNo = NON_RANDOM_NUBMER_USED_AS_DUMMY;
#endif
	return RandomNo;
}
DType PNaiveBaseline			(const DType &HistHazard, const DType &Age)
{
	return exp(-HistHazard * Age);
}
DType PImmuneBaseline			(const DType &HistHazard, const DType &Age)
{
	return (1 - PNaiveBaseline(HistHazard, Age));
}
DType logPNaiveBaseline			(const DType &HistHazard, const DType &Age)		
{
	return (-HistHazard * Age);
}
DType logPImmuneBaseline		(const DType &HistHazard, const DType &Age)		
{
	return log(1 - PNaiveBaseline(HistHazard, Age));
}
DType K_seropos					(const DType &HistHazard, const DType &Age, const DType &K_1, const DType &K_2) 
{
	DType P0 = exp(-HistHazard * Age); /// naive
	DType P_Survive_single_serotype = exp(-HistHazard * Age / 4); 
	DType P1 = 4 * (1 - P_Survive_single_serotype) * (P_Survive_single_serotype*P_Survive_single_serotype*P_Survive_single_serotype); /// infected exactly once. 
	DType Quantity = P1 / (1 - P0); //// Proportion people infected exactly once out of seropositives. 
	DType K_plus =	(		Quantity  * K_1) +
					((1 -	Quantity) * K_2) ;
	return K_plus;
}
DType Prob_SNeg_PrimeOne		(int Age, DType Prop_Indpt_Age, DType SNegRate) //// For HOUSE.ModelVariant == AS_PRIME. Check functional form, but should get smaller as age increases
{
	return Prop_Indpt_Age + (1 - Prop_Indpt_Age) * exp(-(DType)Age * SNegRate);
}
DType Prob_SPos_NoPriming		(int Age, DType Prop_Indpt_Age, DType SPosRate) //// For HOUSE.ModelVariant == AS_PRIME
{
	return (1 - Prop_Indpt_Age) * exp(-(DType)Age * SPosRate);
}
DType K0_Prime					(int Age, DType K1, DType K2, DType Prop_Indpt_Age, DType SNegRate)
{
	DType Prob_Move_One_age = Prob_SNeg_PrimeOne(Age, Prop_Indpt_Age, SNegRate);

	return Prob_Move_One_age * K1 + (1 - Prob_Move_One_age) * K2;
}
DType K1_Prime					(int Age, DType K1, DType K2, DType Prop_Indpt_Age, DType SPosRate)
{
	DType Prob_SPos_NoPrime = Prob_SPos_NoPriming(Age, Prop_Indpt_Age, SPosRate);

	return Prob_SPos_NoPrime * K1 + (1 - Prob_SPos_NoPrime) * K2;
}
DType UniformPrior				(DType LowerLim, DType UpperLim)
{
	return 1 / (UpperLim - LowerLim);
}
DType Log_UniformPrior			(DType LowerLim, DType UpperLim)
{
	//// not the same as prior used for logged knots. This is logged because you consider the log likelihhood, log posterior, and log prior. 
	return log(UniformPrior(LowerLim, UpperLim)); 
}
DType UniformPrior_Log10		(DType x, DType LowerLim_log, DType UpperLim_log)
{
	///// NOTE: this is the prior probability of the parameter ON THE NON-LOG scale. The log10 is here as a conversion factor. 
	return 1 / (x * log(10) * (UpperLim_log - LowerLim_log)); 
}
DType Log_UniformPrior_Log10	(DType x, DType LowerLim, DType UpperLim)
{
	//// Double logged because, by fitting log10 of knots, you have a (1/x*ln(10)) prior. You then need the (natural) log of this prior, as calculating log likelihhood, log posterior etc. 
	return log(UniformPrior_Log10(x, LowerLim, UpperLim));
}

void CalculateLogPriors_single	(int param_no, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	if (PARAMS.ParamRanges[UpperBound][param_no] - PARAMS.ParamRanges[LowerBound][param_no] == 0) PARAMS.ParamArray_logPriors[param_no] = 1; 
	else if (IsParamAKnot(param_no, PARAMS.ParamNumbers))
		PARAMS.ParamArray_logPriors[param_no] = Log_UniformPrior_Log10	(PARAMS.ParamVec[param_no], PARAMS.ParamRanges[LowerBound][param_no], PARAMS.ParamRanges[UpperBound][param_no]);
	else	
		PARAMS.ParamArray_logPriors[param_no] = Log_UniformPrior		(							PARAMS.ParamRanges[LowerBound][param_no], PARAMS.ParamRanges[UpperBound][param_no]);
}
void CalculateLogPriors_all		(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	//// reset log prior. 
	PARAMS.LogPriorFull = 0; 
	for (int param_no = 0; param_no < HOUSE.No_Parameters; param_no++)
	{
		//// change component of PARAMS.ParamArray_logPriors
		CalculateLogPriors_single(param_no, PARAMS, HOUSE);

		//// sum log prior (multiply non-log prior). 
		PARAMS.LogPriorFull += PARAMS.ParamArray_logPriors[param_no]; 
	}
}
void CalculateTotalLogPrior		(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	DType TotalLogPrior = 0;  //// i.e. Total (non-logged) prior  = 1
	for (int param_no = 0; param_no < HOUSE.No_Parameters; param_no++) TotalLogPrior += PARAMS.ParamArray_logPriors[param_no]; 
	PARAMS.LogPriorFull = TotalLogPrior; 
}
void AmendLogPrior				(int param_no, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	//// Assumes that PARAMS.ParamVec[param_no] has been changed previously. 
	if (IsParamAKnot(param_no, PARAMS.ParamNumbers)) ///// currently Knots are the only parameters with non-uniform priors. Hence you only need to amend the prior vector and (log) product for knots. 
	{
		DType Old_LogParameterPrior = PARAMS.ParamArray_logPriors[param_no];
		//// change component of PARAMS.ParamArray_logPriors
		CalculateLogPriors_single(param_no, PARAMS, HOUSE); //// changes ProposedPARAMS.ParamArray_logPriors[param_no];
		//// change Full prior
		PARAMS.LogPriorFull = PARAMS.LogPriorFull - Old_LogParameterPrior + PARAMS.ParamArray_logPriors[param_no]; //// i.e. divide by old, multiply by new on non-logged scale. 
	}
}


DType DIC						(int NumPosteriorSamples, DType LogLike_Value, DType * Chain)
{
	DType MeanLogLike		= 0;
	/// calculate mean
	for (int sample = 0; sample < NumPosteriorSamples; sample++) MeanLogLike += Chain[sample];
	MeanLogLike /= NumPosteriorSamples;
	/// calculate variance
	DType Var_Chain = 0;
	for (int sample = 0; sample < NumPosteriorSamples; sample++) Var_Chain += (Chain[sample] - MeanLogLike) * (Chain[sample] - MeanLogLike); /// i.e. squared, but faster to write it this way. 
	Var_Chain /= NumPosteriorSamples;

	//// calc DIC = LogLike(Data | MeanOrModeParameters) - 
	return -2 * (LogLike_Value - 2 * Var_Chain);
}

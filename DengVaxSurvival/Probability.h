#ifndef PROBABILTIY_FUNCTIONS_HEADER_INCLUDED
#define PROBABILTIY_FUNCTIONS_HEADER_INCLUDED

#include "StructureDefs.h"

DType RandomNumber_Mt			(int thread_no);
DType logPNaiveBaseline			(const DType &HistHazard, const DType &Age);
DType logPImmuneBaseline		(const DType &HistHazard, const DType &Age);
DType PNaiveBaseline			(const DType &HistHazard, const DType &Age);
DType PImmuneBaseline			(const DType &HistHazard, const DType &Age);
DType K_seropos					(const DType &HistHazard, const DType &Age, const DType &K_1, const DType &K_2);
DType Prob_SNeg_PrimeOne		(int Age, DType Prop_Indpt_Age, DType SNegRate);						//// For HOUSE.ModelVariant == AS_PRIME. 
DType Prob_SPos_NoPriming		(int Age, DType Prop_Indpt_Age, DType SPosRate);						//// For HOUSE.ModelVariant == AS_PRIME
DType K0_Prime					(int Age, DType K1, DType K2, DType Prop_Indpt_Age, DType SNegRate);
DType K1_Prime					(int Age, DType K1, DType K2, DType Prop_Indpt_Age, DType SPosRate);
DType UniformPrior				(DType LowerLim, DType UpperLim);
DType UniformPrior_Log10		(DType x, DType LowerLim_log, DType UpperLim_log);
DType Log_UniformPrior			(DType LowerLim, DType UpperLim); 
DType Log_UniformPrior_Log10	(DType x, DType LowerLim, DType UpperLim); 
void CalculateLogPriors_single	(int param_no, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE); 
void CalculateLogPriors_all		(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE); 
void CalculateTotalLogPrior		(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE); 
void AmendLogPrior				(int param_no, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE); 
DType DIC						(int NumPosteriorSamples, DType LogLike_Value, DType * Chain);

#endif

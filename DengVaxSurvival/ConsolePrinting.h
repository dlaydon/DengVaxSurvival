#pragma once

#include "HeaderAndSwitches.h"

string WhichLikeComp(int ModuloComponent, const Housekeeping_Struct &HOUSE);

void PrintKsArray(int country, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
void printVector		(std::vector<int> VEC, string VecName);
void printVector		(std::vector<DType> VEC, string VecName);
void printVector		(std::vector<string> VEC, string VecName);
void print_LikeArray	(DType ** &LikeParts, const Housekeeping_Struct &HOUSE, int country);

void PrintNaNComponents			(const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
void PrintUnequalComponents		(DType ** &CurrentLikeParts, DType ** &ProposedLikeParts, const Housekeeping_Struct &HOUSE);
void Print_MCMC_Progress		(int iter, Chains_Struct &CHAINS, Params_Struct &CurrentPARAMS); 
void MaybePrint					(int param_no, DType Proposed_Param, Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE, bool ParamAcceptCondition);
void PrintVacHazLikeDiscrep		(int param_no, Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, DType *** VacHazLike_CheckingArray);
void TestAndPrintParamNumberFucntions	(Housekeeping_Struct &HOUSE, Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS);
void PrintFittedParamNumbers			(const Housekeeping_Struct &HOUSE, Params_Struct &CurrentPARAMS);
void Print_Like_Indices					(Housekeeping_Struct &HOUSE, Params_Struct &CurrentPARAMS);
void Check_Age_ParamsEtc(Age_Option AS_VE_Wane_Haz, const Housekeeping_Struct &HOUSE, DType * ParamContainer, DType * MultContainer, DType ** CoeffContainer, int NumFunctionParams, int PolysPerSpline, int MaxSplineDegree, Params_Struct &PARAMS);
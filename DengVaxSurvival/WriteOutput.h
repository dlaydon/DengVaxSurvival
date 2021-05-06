#pragma once

#include "StructureDefs.h"


void WriteOutput					(string FILENAME, DType **Object, int NRows, int NCols); 
void WriteParameterChainOutput		(string FILENAME, std::vector<string> NamesOfParameters, const Housekeeping_Struct &HOUSE, Chains_Struct &CHAINS);
void WriteAcceptanceArray			(string FILENAME, DType * AcceptanceArray, DType * AcceptanceArray_Aug, int max_threads, int NumPosteriorSamples, int No_Parameters, int NoAugmented, int No_Iterations, int AugmentEveryHowManyIterations, std::vector<string> NamesOfParameters); 
void WriteSeroPrevOutput			(string FILENAME, DType **AugOnlySeroPrevs, DType **AllSeroPrevs, DType ***AgeSpecificSeroPrevs, int HowManyAges, int NoCountries);
void WriteSingleChain				(string FILENAME, DType *Chain, int NumPosteriorSamples, string Colname);
void WriteAllOutput					(const DATA_struct &DATA, Survival_Struct &SURVIVE, Survival_Struct &SURVIVE_ImSub, const Params_Struct &CurrentPARAMS, const Housekeeping_Struct &HOUSE, Chains_Struct &CHAINS);
void WriteSurvivalTable				(string FILENAME, DType **Object, int NRows, int NCols, std::vector<string> RowNames, std::vector<string> ColNames);
void WriteAugDataOutput				(int * SerostatusArray, string FileName);
void WriteAugDataOutput				(DType * SerostatusArray, string FileName);

void WriteEverythingForParticularChain(DATA_struct &DATA, Survival_Struct &SURVIVE, Survival_Struct &SURVIVE_ImSub,
	Params_Struct &CurrentPARAMS, const Housekeeping_Struct &HOUSE, Chains_Struct &CHAINS, Augment_Struct &AUG); 

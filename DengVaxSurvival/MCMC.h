#pragma once

#include "StructureDefs.h"

std::vector<int>	FindSimulataneousUpdateParams	(int &param_no, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
std::vector<DType>	Populate_Proposed_Param_Vec		(std::vector<int> &SimulataneousUpdateParams, int &param_no, DType Proposed_Param, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
std::vector<DType>	Calc_RowMeans					(DType ** ParamChain, int NumRows, int NumCols); 

DType ProposeNewParam	(int &param_no, int &thread_num, const Params_Struct &CurrentPARAMS, const Housekeeping_Struct &HOUSE);
bool MCMC_step_param	(int param_no, const DATA_struct &DATA, Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE, Chains_Struct &CHAINS);
void AddToChains		(int iter, Chains_Struct &CHAINS, const DATA_struct &DATA, Params_Struct &CurrentPARAMS, const Housekeeping_Struct &HOUSE, Survival_Struct &SURVIVE, Survival_Struct &SURVIVE_ImSub, Augment_Struct &AUG);
void MCMC_All			(int FirstIteration, Chains_Struct &CHAINS, DATA_struct &DATA, Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, Augment_Struct &AUG, const Housekeeping_Struct &HOUSE, Survival_Struct &SURVIVE, Survival_Struct &SURVIVE_ImSub);


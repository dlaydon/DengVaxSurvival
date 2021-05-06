#ifndef PARAM_UPDATE_FUNCTIONS_INCLUDED
#define PARAM_UPDATE_FUNCTIONS_INCLUDED

#include "StructureDefs.h"

/*

This script (and associated header) has two kinds of functions:
	
	1)	Functions that amend your PARAMS structures with various precalculated quantities by calling the relevant funtions in CalcParamsEtc.cpp (see. note at top of CalcParamsEtc.cpp). e.g. changing rho will change SumRhoEffs, and could affect stored likelihood values of vaccine hazard etc. 
	2)	Functions that amend the likelihood components AFTER parameter sturctures have been set appropriately. These are the Change_ functions below. 

*/


//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** 
//// ****			1) AmendParams and various overloads essentially call functions in CalcParamsEtc when required, e.g. if changing historical hazrd, need to change array of KPlusValues BEFORE you amend the likelihood. 
//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** 

void AmendParams				(int &param_no, DType &New_Param, Params_Struct &PARAMS, const DATA_struct &DATA, const Housekeeping_Struct &HOUSE);
void AmendParams				(std::vector<DType> ParamVec, std::vector<int> PVecIndices, Params_Struct &PARAMS, const DATA_struct &DATA, const Housekeeping_Struct &HOUSE);
void AmendParams				(std::vector<DType> ParamVec, Params_Struct &PARAMS, const DATA_struct &DATA, const Housekeeping_Struct &HOUSE);
void UpdateOrReset_ParamsEtc	(int &param_no, const DATA_struct &DATA, Params_Struct &ParamsToUpdate, const Params_Struct &CorrectPARAMS, const Housekeeping_Struct &HOUSE); ///  //// because things are now a function of PARAMS as a whole, you need to reset more things than when objects like lookup tables were floating around unstructured. 

//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** 
//// ****			2) Likelihood amendments
//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** 

void Change_K			(int PhaseSeverity, int PrevInf, int serotype, const DATA_struct &DATA, const Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE, const std::vector<int> &WhichCountries);
void Change_K			(int PhaseSeverity, int PrevInf, int serotype, const DATA_struct &DATA, const Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE);
void Change_Efficacy	(const DATA_struct &DATA, const Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE, char BaselineSeroStatus, int &serotype, int PhaseSeverity);
void Change_Waning		(const DATA_struct &DATA, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE, char BaselineSeroStatus);
void Change_Knot		(int &country, int WhichKnot, int param_no, const DATA_struct &DATA, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE);
void Change_HistHazard	(int &country, const DATA_struct &DATA, const Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE);
void Change_Rho			(int &country, int &serotype, const DATA_struct &DATA, const Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE);
void Change_qVal		(int &qval, int &country, const DATA_struct &DATA, const Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE);
void Change_Param		(int param_no, DType &New_Param, const DATA_struct &DATA, const Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE);

#endif
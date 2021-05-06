
#ifndef PARAM_NUMBER_FUNCTIONS_HEADER_INCLUDED
#define PARAM_NUMBER_FUNCTIONS_HEADER_INCLUDED

#include "StructureDefs.h"

bool IsParamAKnot				(const int &param_no, const ParamNumbers_Struct &ParamNumbers);
bool IsParamAHistHaz			(const int &param_no, const ParamNumbers_Struct &ParamNumbers);
bool IsParamA_qval				(const int &param_no, const ParamNumbers_Struct &ParamNumbers);
bool IsParamA_rho				(const int &param_no, const ParamNumbers_Struct &ParamNumbers);
bool IsParamAnEfficacy			(const int &param_no, const ParamNumbers_Struct &ParamNumbers);
bool IsParamAn_atInf_Efficacy	(const int &param_no, const ParamNumbers_Struct &ParamNumbers);
bool IsParamARelativeRisk		(const int &param_no, const ParamNumbers_Struct &ParamNumbers);
bool IsParamARelativeRisk_Hosp	(const int &param_no, const ParamNumbers_Struct &ParamNumbers);
bool IsParamAHillHalflife		(const int &param_no, const ParamNumbers_Struct &ParamNumbers); 
bool IsParamAHillPower			(const int &param_no, const ParamNumbers_Struct &ParamNumbers);
bool IsParamAWaningParam		(const int &param_no, const ParamNumbers_Struct &ParamNumbers); 
bool IsParamAn_ASVE_Param		(const int &param_no, const ParamNumbers_Struct &ParamNumbers);
bool IsParamAn_ASHaz_Param		(const int &param_no, const ParamNumbers_Struct &ParamNumbers);
bool IsParamAn_AS_Prime_Param	(const int &param_no, const ParamNumbers_Struct &ParamNumbers);


int Find_Knot_FromKnotParam					(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers);
int Find_Country_FromKnotParam				(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers, const int &knot);
int Find_Country_FromHistHazParam			(const int &param_no, const ParamNumbers_Struct &ParamNumbers);
int Find_qval_From_qParam					(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers);
int Find_Country_From_qParam				(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers, const int &qval);
int Find_rho_From_rhoParam					(const int &param_no, const int &Min_rho_no, const int &No_STypes);
int Find_Country_From_rhoParam				(const int &param_no, const int &Min_rho_no, const int &No_STypes, const int &rho);
int Find_SeroStatus_FromEffParam			(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers);
int Find_Serotype_FromEffParam				(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers, const int &SeroStatus);
int Find_PhaseSeverity_From_EffParam		(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers, const int &SeroStatus, const int &serotype); 
int Find_PrevInf_From_K_Param				(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers);
int Find_PhaseSeverity_From_K_Param			(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers, const int &PrevInf);
int Find_Serotype_From_K_Param				(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers, const int &PrevInf, const int &PhaseSeverity); 
int Find_SeroStatus_From_HalfLife			(const int &param_no, const int &Min_HillHalfLife);
int Find_SeroStatus_From_HillCoeff			(const int &param_no, const int &Min_HillPower);
int Find_SeroStatus_FromASVE_Param			(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers, const int &type);
int Find_SeroStatus_FromASVE_Param			(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers);
int Find_Type_From_ASVE_Param				(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers);
int Find_Type_From_ASHaz_Param				(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers);
int Find_ParamType_From_ASWaning_Param		(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers); 
int Find_SeroStatus_From_Waning_Param		(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers, const int &ParamType); 
int Find_SeroStatus_From_Waning_Param		(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers);
int Find_ParamType_From_AS_Prime_Param		(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers);
int Find_SeroStatus_From_AS_Prime_Param		(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers, const int &ParamType);

#endif
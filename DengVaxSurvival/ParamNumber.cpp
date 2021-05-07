
#include "HeaderAndSwitches.h"

///// These are all the Same function but helpful to have them named. So too with Find functions below. 

bool IsParamAKnot				(const int &param_no, const ParamNumbers_Struct &ParamNumbers)
{
	return (param_no >= ParamNumbers.Min_Knot) && (param_no <= ParamNumbers.Max_Knot); 
}
bool IsParamAHistHaz			(const int &param_no, const ParamNumbers_Struct &ParamNumbers)
{
	return (param_no >= ParamNumbers.Min_HHaz) && (param_no <= ParamNumbers.Max_HHaz);
}
bool IsParamA_qval				(const int &param_no, const ParamNumbers_Struct &ParamNumbers)
{
	return (param_no >= ParamNumbers.Min_qval) && (param_no <= ParamNumbers.Max_qval);
}
bool IsParamA_rho				(const int &param_no, const ParamNumbers_Struct &ParamNumbers)
{
	return (param_no >= ParamNumbers.Min_Rho) && (param_no <= ParamNumbers.Max_Rho);
}
bool IsParamAnEfficacy			(const int &param_no, const ParamNumbers_Struct &ParamNumbers)
{
	return (param_no >= ParamNumbers.Min_VacE) && (param_no <= ParamNumbers.Max_VacE);
}
bool IsParamAn_atInf_Efficacy	(const int &param_no, const ParamNumbers_Struct &ParamNumbers)
{
	return (param_no >= ParamNumbers.Min_VacE_atInf) && (param_no <= ParamNumbers.Max_VacE_atInf);
}
bool IsParamARelativeRisk		(const int &param_no, const ParamNumbers_Struct &ParamNumbers)
{
	return (param_no >= ParamNumbers.Min_K) && (param_no <= ParamNumbers.Max_K);
}
bool IsParamARelativeRisk_Hosp	(const int &param_no, const ParamNumbers_Struct &ParamNumbers)
{
	return (param_no >= ParamNumbers.Min_Hosp_K) && (param_no <= ParamNumbers.Max_Hosp_K);
}
bool IsParamAHillHalflife		(const int &param_no, const ParamNumbers_Struct &ParamNumbers)
{
	return (param_no >= ParamNumbers.Min_HillHalfLife) && (param_no <= ParamNumbers.Max_HillHalfLife);
}
bool IsParamAHillPower			(const int &param_no, const ParamNumbers_Struct &ParamNumbers)
{
	return (param_no >= ParamNumbers.Min_HillPower) && (param_no <= ParamNumbers.Max_HillPower);
}
bool IsParamAWaningParam		(const int &param_no, const ParamNumbers_Struct &ParamNumbers)
{
	return( (param_no >= ParamNumbers.MinWaningParam && param_no <= ParamNumbers.MaxWaningParam)	|| 
			IsParamAHillHalflife(param_no, ParamNumbers) || IsParamAHillPower(param_no, ParamNumbers));
}
bool IsParamAn_ASVE_Param		(const int &param_no, const ParamNumbers_Struct &ParamNumbers)
{
	return (param_no >= ParamNumbers.Min_ASVE_Param) && (param_no <= ParamNumbers.Max_ASVE_Param);
}
bool IsParamAn_ASHaz_Param		(const int &param_no, const ParamNumbers_Struct &ParamNumbers)
{
	return (param_no >= ParamNumbers.Min_ASHaz_Param) && (param_no <= ParamNumbers.Max_ASHaz_Param);
}
bool IsParamAn_AS_Prime_Param	(const int &param_no, const ParamNumbers_Struct &ParamNumbers)
{
	return (param_no >= ParamNumbers.Min_AS_PrimingParam) && (param_no <= ParamNumbers.Max_AS_PrimingParam);
}

int Find_Final_Index_Array			(const int &param_no, const int &Offset, const int &FinalIndexDim)
{
	return ((param_no				- Offset) % FinalIndexDim);
}
int Find_Penultimate_ArrayIndex		(const int &param_no, const int &Offset, const int &FinalIndexDim, const int &Final_Index) //// don't delete this yet as Find_Country_From_rhoParam can't take the overloaded version yet. 
{
	return ((param_no - Final_Index - Offset) / FinalIndexDim);
}
int Find_Penultimate_ArrayIndex		(const int &param_no, const int &Offset, const int &FinalIndexDim, const int &Final_Index, const int &PenultimateIndex_Dim)
{
	return(((param_no - Offset - Final_Index) / FinalIndexDim) % PenultimateIndex_Dim); 
}
int Find_3rdLast_ArrayIndex			(const int &param_no, const int &Offset, const int &FinalIndexDim, const int &Final_Index, const int &PenultimateIndexDim, const int &Penultimate_Index, const int& ThirdLastIndexDim)
{
	return (		(((param_no - Offset - Final_Index) / FinalIndexDim) - Penultimate_Index) / PenultimateIndexDim				) % ThirdLastIndexDim;
}

int Find_Knot_FromKnotParam					(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers)
{
	return (Find_Final_Index_Array(param_no, ParamNumbers.Min_Knot, HOUSE.KnotsPerCountry));
}
int Find_Country_FromKnotParam				(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers, const int &knot)
{
	return (Find_Penultimate_ArrayIndex(param_no, ParamNumbers.Min_Knot, HOUSE.KnotsPerCountry, knot, HOUSE.TotalCountries));
}
int Find_Country_FromHistHazParam			(const int &param_no, const ParamNumbers_Struct &ParamNumbers)
{
	return (param_no - ParamNumbers.Min_HHaz);
}
int Find_qval_From_qParam					(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers)
{
	return (Find_Final_Index_Array(param_no, ParamNumbers.Min_qval, (HOUSE.N_STypes - 1)));
}
int Find_Country_From_qParam				(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers, const int &qval)
{
	return (Find_Penultimate_ArrayIndex(param_no, ParamNumbers.Min_qval, (HOUSE.N_STypes - 1), qval, HOUSE.TotalCountries));
}
int Find_rho_From_rhoParam					(const int &param_no, const int &Min_rho_no	, const int &No_STypes)
{
	return (Find_Final_Index_Array(param_no, Min_rho_no, No_STypes));
}
int Find_Country_From_rhoParam				(const int &param_no, const int &Min_rho_no	, const int &No_STypes			, const int &rho)
{
	return (Find_Penultimate_ArrayIndex(param_no, Min_rho_no, No_STypes, rho));
}
int Find_SeroStatus_FromEffParam			(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers)
{
		 if (IsParamAnEfficacy			(param_no, ParamNumbers))	return (Find_Final_Index_Array(param_no, ParamNumbers.Min_VacE		, HOUSE.HowManySeroStatuses));
	else if (IsParamAn_atInf_Efficacy	(param_no, ParamNumbers))	return (Find_Final_Index_Array(param_no, ParamNumbers.Min_VacE_atInf, HOUSE.HowManySeroStatuses));
	else return 0; // Should never happen.
}
int Find_Serotype_FromEffParam				(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers, const int &SeroStatus)
{
		 if (IsParamAnEfficacy			(param_no, ParamNumbers))	return (Find_Penultimate_ArrayIndex(param_no, ParamNumbers.Min_VacE			, HOUSE.HowManySeroStatuses, SeroStatus, HOUSE.N_STypes_VEs));
	else if (IsParamAn_atInf_Efficacy	(param_no, ParamNumbers))	return (Find_Penultimate_ArrayIndex(param_no, ParamNumbers.Min_VacE_atInf	, HOUSE.HowManySeroStatuses, SeroStatus, HOUSE.N_STypes_VEs));
	else return 0; // Should never happen.
}
int Find_PhaseSeverity_From_EffParam		(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers, const int &SeroStatus, const int &serotype)
{
		 if (IsParamAnEfficacy			(param_no, ParamNumbers))	return (Find_3rdLast_ArrayIndex(param_no, ParamNumbers.Min_VacE			, HOUSE.HowManySeroStatuses, SeroStatus, HOUSE.N_STypes_VEs, serotype, HOUSE.HowManyCaseCategories));
	else if (IsParamAn_atInf_Efficacy	(param_no, ParamNumbers))	return (Find_3rdLast_ArrayIndex(param_no, ParamNumbers.Min_VacE_atInf	, HOUSE.HowManySeroStatuses, SeroStatus, HOUSE.N_STypes_VEs, serotype, HOUSE.HowManyCaseCategories));
	else return 0; // Should never happen.
}
int Find_PrevInf_From_K_Param				(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers)
{
		 if (IsParamARelativeRisk		(param_no, ParamNumbers))	return (Find_Final_Index_Array(param_no, ParamNumbers.Min_K		, HOUSE.Num_K_Params));
	else if (IsParamARelativeRisk_Hosp	(param_no, ParamNumbers))	return (Find_Final_Index_Array(param_no, ParamNumbers.Min_Hosp_K, HOUSE.Num_K_Params));
	else return 0; // Should never happen.
}
int Find_PhaseSeverity_From_K_Param			(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers, const int &PrevInf)
{
	  	 if (IsParamARelativeRisk		(param_no, ParamNumbers))								return (Find_Penultimate_ArrayIndex(param_no, ParamNumbers.Min_K, HOUSE.Num_K_Params, PrevInf, HOUSE.HowManyCaseCategories)); 
	else if (IsParamARelativeRisk_Hosp	(param_no, ParamNumbers) && !HOUSE.PASSIVE_PHASE_ONLY)	return (	PassiveSevere	);
	else if (IsParamARelativeRisk_Hosp	(param_no, ParamNumbers) &&  HOUSE.PASSIVE_PHASE_ONLY)	return (	ActiveMild		); /// because relabelled passive phase to be "active phase" to make code run more easily.
	else return 0; // Should never happen.
}
int Find_Serotype_From_K_Param				(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers, const int &PrevInf, const int &PhaseSeverity)
{
		 if (IsParamARelativeRisk		(param_no, ParamNumbers))	return (	Find_3rdLast_ArrayIndex(param_no, ParamNumbers.Min_K, HOUSE.Num_K_Params, PrevInf, HOUSE.HowManyCaseCategories, PhaseSeverity, HOUSE.N_STypes_Ks)	);
	else if (IsParamARelativeRisk_Hosp	(param_no, ParamNumbers))	return (	  (param_no - ParamNumbers.Min_Hosp_K	- PrevInf) / HOUSE.Num_K_Params);
	else return 0; // Should never happen.
}
int Find_Type_From_ASVE_Param				(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers)
{
	return ((param_no - ParamNumbers.Min_ASVE_Param) % HOUSE.NumASVE_Function_Params);
}
int Find_SeroStatus_FromASVE_Param			(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers, const int &type)
{
	return ((param_no - type - ParamNumbers.Min_ASVE_Param) / HOUSE.NumASVE_Function_Params);
}
int Find_SeroStatus_FromASVE_Param			(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers)
{
	int type				= Find_Type_From_ASVE_Param		(param_no, HOUSE, ParamNumbers); 
	int BaselineSeroStatus	= Find_SeroStatus_FromASVE_Param(param_no, HOUSE, ParamNumbers, type);
	return BaselineSeroStatus;
}
int Find_Type_From_ASHaz_Param				(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers)
{
	return ((param_no - ParamNumbers.Min_ASHaz_Param) % HOUSE.NumAS_Haz_Function_Params);
}
int Find_SeroStatus_From_HalfLife			(const int &param_no, const int &Min_HillHalfLife)
{
	return (param_no - Min_HillHalfLife);
}
int Find_SeroStatus_From_HillCoeff			(const int &param_no, const int &Min_HillPower)
{
	return (param_no - Min_HillPower);
}

int Find_ParamType_From_ASWaning_Param		(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers)
{
		 if (HOUSE.HillWaning	)	return 0;
	else if (HOUSE.AS_Waning == Age_Option::INDEPENDENT	)	return 0;
	else return ((param_no - ParamNumbers.MinWaningParam) % HOUSE.NumWaningParamsPer_BS);
}
int Find_SeroStatus_From_Waning_Param		(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers, const int &ParamType)
{
	int BaselineSeroStatus = -8000; //// chosen to cause crash if not set properly. 

	if (IsParamAWaningParam(param_no, ParamNumbers))
	{
		if (IsParamAHillHalflife	(param_no, ParamNumbers))	BaselineSeroStatus = Find_SeroStatus_From_HalfLife	(param_no, ParamNumbers.Min_HillHalfLife);
		else if (IsParamAHillPower		(param_no, ParamNumbers))	BaselineSeroStatus = Find_SeroStatus_From_HillCoeff	(param_no, ParamNumbers.Min_HillPower	);
		else														BaselineSeroStatus = (param_no - ParamType - ParamNumbers.MinWaningParam) / HOUSE.NumWaningParamsPer_BS;
	}
	return BaselineSeroStatus;
}
int Find_SeroStatus_From_Waning_Param		(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers)
{
	int ParamType			= Find_ParamType_From_ASWaning_Param(param_no, HOUSE, ParamNumbers); 
	int BaselineSeroStatus	= Find_SeroStatus_From_Waning_Param	(param_no, HOUSE, ParamNumbers, ParamType); 
	return BaselineSeroStatus;
}
int Find_ParamType_From_AS_Prime_Param		(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers)
{
	return ((param_no - ParamNumbers.Min_AS_PrimingParam) % HOUSE.Num_AS_Priming_ParamsPer_BS);
}
int Find_SeroStatus_From_AS_Prime_Param		(const int &param_no, const Housekeeping_Struct &HOUSE, const ParamNumbers_Struct &ParamNumbers, const int &ParamType)
{
	return (param_no - ParamType - ParamNumbers.Min_AS_PrimingParam) / HOUSE.Num_AS_Priming_ParamsPer_BS;;
}
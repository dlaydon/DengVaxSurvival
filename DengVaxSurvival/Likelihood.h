
#ifndef LIKELIHOOD_FUNCTIONS_HEADER_INCLUDED
#define LIKELIHOOD_FUNCTIONS_HEADER_INCLUDED

#include "StructureDefs.h"

DType l_K_case					(const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, const int &country, char PhaseSeverity, int PrevInf, int &serotype);
DType l_Kplus_case				(const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, const int &country, char PhaseSeverity, int &serotype);
DType l_K0SNegPrime_case		(const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, const int &country, char PhaseSeverity, int &serotype);
DType l_KPlusPrime_case			(const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, const int &country, char PhaseSeverity, int &serotype);

DType l_BaseHaz						(const int &country, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, char PhaseSeverity);
DType l_Aug							(const int &country, int BaselineSeroStatus, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
DType l_Aug							(const int &country, int BaselineSeroStatus, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, bool ImSubOnly);

DType l_Int_Base_Haz				(const int &country, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, char PhaseSeverity, char BaselineSeroStatus, char TrialArm);
DType l_Int_Vac_Haz					(char BaselineSeroStatus, const int &country, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, char PhaseSeverity);
DType l_Int_Vac_Haz_sero_i			(const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, const int &country, int serotype, char BaselineSeroStatus, char PhaseSeverity, DType *** &VacHazLikes);
DType l_Int_Vac_Haz_sero_i			(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, const int &country, int serotype, char BaselineSeroStatus, char PhaseSeverity);
DType l_Int_Vac_Haz_sero_i_long_way	(char BaselineSeroStatus, const int &country, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, char PhaseSeverity, int serotype);
void CalcAndStoreSumIntVacHazOverPatients					(char BaselineSeroStatus, const int &country, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, char PhaseSeverity, DType *** &VacHazLikes);
void CalcAndStoreSumIntVacHazOverPatients					(char BaselineSeroStatus, const int &country, const DATA_struct &DATA, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, char PhaseSeverity);

DType WaningEfficacy				(char BaselineSeroStatus, int &serotype, char PhaseSeverity, int AgeInYears, DType WaningMult, EffNegWane_Option ENW_Option, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, bool NaturalLog);
DType l_WaningEfficacy				(char BaselineSeroStatus, const int &country, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int &serotype, char PhaseSeverity);
DType l_rho_case					(const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, const int &country, int serotype, char PhaseSeverity);
DType l_AS_HazMult					(const int &country, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, char PhaseSeverity);
DType l_BS_BaseHaz_Mult(const int& country, const DATA_struct& DATA, const Params_Struct& PARAMS, const Housekeeping_Struct& HOUSE, char PhaseSeverity, char BaselineSeroStatus);
DType l_full						(const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
DType l_full						(const DATA_struct &DATA, DType **LikeArray, const Housekeeping_Struct &HOUSE);
DType l_full_Minus_Aug				(const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);

void LikelihoodFromScratch			(const DATA_struct &DATA, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, DType **LikeParts);  //// in this version, provide the address of an array to populate. While PARAMS has one such array, if you are checking then you will need another. 
void LikelihoodFromScratch			(const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, DType **LikeParts, DType *** &VacHazLike_CheckingArray);

DType LL_Total						(const DATA_struct &DATA, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, Augment_Struct &AUG);

void CalcWeightedLikes			(const DATA_struct &DATA, Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE);
bool Everything_Okay			(int param_no, const DATA_struct &DATA, Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE, string ParamOrAugment);


bool CompareLikelihoods			(string ParamsOrAugment, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, DType **LikePartsForChecking, DType *** &VacHazLike_CheckingArray);
bool Likelihood_Okay			(Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS);



#endif


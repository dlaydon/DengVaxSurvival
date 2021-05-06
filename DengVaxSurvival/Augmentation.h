#ifndef AUGMENTATION_FUNCTIONS_INCLUDED
#define AUGMENTATION_FUNCTIONS_INCLUDED

#include "StructureDefs.h"

void PopulateSets			(DATA_struct &DATA, const Housekeeping_Struct &HOUSE);
void ClearSets				(DATA_struct &DATA, const Housekeeping_Struct &HOUSE);
void UpdateOrResetAugData	(int * IisToChange, int * SourceIis, DATA_struct &DATA, const Housekeeping_Struct &HOUSE);
void InitAugDataValues		(DATA_struct &DATA, const Params_Struct &CurrentPARAMS, const Housekeeping_Struct &HOUSE);

void Choose_HazMults_SingleSero		(int th, int person_i, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, Augment_Struct &AUG);
void Choose_iHazMults				(int th, int person_i, int BaselineSeroStatus,	const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, Augment_Struct &AUG);
void Choose_iHazMults				(int th, int person_i,							const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, Augment_Struct &AUG);
void Adjust_Aug_thread				(int th, int person_i, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, Augment_Struct &AUG);
DType ProbSerostatusGivenOutcome	(int person_i, const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, Augment_Struct &AUG, int th);
void AugmentData					(Augment_Struct &AUG, DATA_struct &DATA, Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE, Chains_Struct &CHAINS); 

#endif


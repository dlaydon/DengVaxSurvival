#ifndef INITIALIZE_FUNCTIONS_INCLUDED
#define INITIALIZE_FUNCTIONS_INCLUDED

#include "StructureDefs.h"

int					set_NCores				();
Age_Option			Convert_AS_String		(const std::string& OptionString);
ExtImSub_Option		Convert_ExtImSub_String	(const std::string& OptionString);
EffNegWane_Option	Convert_EffNegWane_String(const std::string& OptionString);

template <typename TYPE> void DeAllocate_2D_Array	(TYPE **   &OBJECT, int Dim1);
template <typename TYPE> void DeAllocate_3D_Array	(TYPE ***  &OBJECT, int Dim1, int Dim2); 

void Initialize_Params				(Params_Struct &PARAMS, Params_Struct &CurrentPARAMS, Housekeeping_Struct &HOUSE, const DATA_struct &DATA);
void Initialize_IntVacHazardValues	(Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const DATA_struct &DATA, const Housekeeping_Struct &HOUSE);
void Initialize_ParamSeroPrevs		(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, const SeroPrev_Struct &SEROPREV);
void Initialize_ParamSeroPrevs		(Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE, const SeroPrev_Struct &SEROPREV);
void Initialize_FollowUp			(DATA_struct &DATA, const Housekeeping_Struct &HOUSE);
void AllocateMemoryAndPopulateSets	(DATA_struct &DATA, const Housekeeping_Struct &HOUSE); 

void DeleteData(DATA_struct &DATA, const Housekeeping_Struct &HOUSE);

#endif
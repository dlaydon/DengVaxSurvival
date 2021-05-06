#ifndef READ_AND_PROCESS_DATA_HEADER_INCLUDED
#define READ_AND_PROCESS_DATA_HEADER_INCLUDED

#include "StructureDefs.h"

void ReadInParams						(Housekeeping_Struct &HOUSE, Chains_Struct &CHAINS, Chains_Struct &WBIC_CHAINS, Params_Struct &CurrentPARAMS, string pParamFileName);
void InitializeDATApointers				(DATA_struct &DATA, const Housekeeping_Struct &HOUSE); 
void ReadInData							(DATA_struct &DATA, Housekeeping_Struct &HOUSE);
void ProcessData						(DATA_struct &DATA, const Housekeeping_Struct &HOUSE);
bool CheckCaseDefinitions				(const DATA_struct &DATA, const Housekeeping_Struct &HOUSE);
bool FileExists							(string FileName); 
std::vector<DType> ReadOldParamChain	(string ParamChainFileName, int NoParams, bool PrintEverything);
void ProcessOldParamChain				(DATA_struct &DATA, const Housekeeping_Struct &HOUSE, Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS);
void ImportOldParamChains				(const Housekeeping_Struct &HOUSE, Chains_Struct &CHAINS, Chains_Struct &WBIC_CHAINS); 

#endif


#include "HeaderAndSwitches.h"
#include "Initialize.h"
#include "ParamUpdate.h"
#include "Augmentation.h"
#include "ConsolePrinting.h"



bool FileExists				(string FileName)
{
	std::ifstream infile(FileName);
	bool FileOkay = infile.good();
	infile.close(); 
	return FileOkay;
}
void ReadInParams			(Housekeeping_Struct &HOUSE, Chains_Struct &CHAINS, Chains_Struct &WBIC_CHAINS, Params_Struct &CurrentPARAMS, string pParamFileName)
{
#ifdef USE_COMMAND_LINE

	ifstream ParamsEtc; 
	std::cerr << "Reading Param file: " << pParamFileName ; 
	
	string param_name, param_value_string, ModVariantStringDummy, SingleOrMultiDoseStringDummy, ActiveOrHospitalStringDummy, MildAndSevereStringDummy, FittingHazardsStringDummy, AugMH_or_GibbsDummy;

	if (FileExists(pParamFileName))
	{
		ParamsEtc.open(pParamFileName);
		while (!ParamsEtc.eof())
		{
			//// Format is get the parameter name, then move to next tab, then get the parameter value, then move to next line.
			std::getline(ParamsEtc, param_name			, '\t');
			std::getline(ParamsEtc, param_value_string	, '\n');

			//// these contribute to the output string. 
					if (param_name == "ModelVariant"					)	ModVariantStringDummy			= param_value_string;				//// string
			else	if (param_name == "SingleOrMultiDose"				)	SingleOrMultiDoseStringDummy	= param_value_string;				//// string
			else	if (param_name == "ActiveOrPassivePhase"			)	ActiveOrHospitalStringDummy		= param_value_string;				//// string
			else	if (param_name == "MildAndSevere"					)	MildAndSevereStringDummy		= param_value_string;				//// string
			else	if (param_name == "AugMH_or_Gibbs"					)	AugMH_or_GibbsDummy				= param_value_string;				//// string
			else	if (param_name == "LinKnts"							)	HOUSE.LinKnts					= std::stoi(param_value_string);	//// string
			else	if (param_name == "EffNegWane"						)	HOUSE.EffNegWane				= Convert_EffNegWane_String(param_value_string);	//// string
			
			else	if (param_name == "FitWaningRate"					)	HOUSE.FitWaningRate				= std::stoi(param_value_string);
			else	if (param_name == "HillWaning"						)	HOUSE.HillWaning				= std::stoi(param_value_string);	//// bool
			else	if (param_name == "ResidEffs"						)	HOUSE.ResidEffs					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "SSASVE_Additive"					)	HOUSE.SSASVE_Additive			= std::stoi(param_value_string);	//// bool
			
			else	if (param_name == "SS_VEs"							)	HOUSE.SeroSpecificEfficacies	= std::stoi(param_value_string);	//// bool
			else	if (param_name == "SS_Ks"							)	HOUSE.SeroSpecific_K_values		= std::stoi(param_value_string);	//// bool

			else	if (param_name == "SS_KAM_0"						)	HOUSE.SS_KAM_0					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "SS_KAM_2"						)	HOUSE.SS_KAM_2					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "SS_KPS_0"						)	HOUSE.SS_KPS_0					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "SS_KPS_1"						)	HOUSE.SS_KPS_1					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "SS_KPS_2"						)	HOUSE.SS_KPS_2					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "PS_Ks_Multiply_AM_Ks"			)	HOUSE.PS_Ks_Multiply_AM_Ks		= std::stoi(param_value_string);	//// bool
			else	if (param_name == "BaselinePartition"				)	HOUSE.BaselinePartition			= std::stoi(param_value_string);	//// bool
			else	if (param_name == "PASSIVE_PHASE_ONLY"				)	HOUSE.PASSIVE_PHASE_ONLY		= std::stoi(param_value_string);	//// bool
			
			else	if (param_name == "SingleEff"						)	HOUSE.SingleEff					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Single_SNeg_Eff"					)	HOUSE.Single_SNeg_Eff			= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Single_SPos_Eff"					)	HOUSE.Single_SPos_Eff			= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Which_SNeg_SeroFitsAll"			)	HOUSE.Which_SNeg_SeroFitsAll	= std::stoi(param_value_string);	//// int
			else	if (param_name == "Which_SPos_SeroFitsAll"			)	HOUSE.Which_SPos_SeroFitsAll	= std::stoi(param_value_string);	//// int
			else	if (param_name == "RelRisksSameBtwSerotypes"		)	HOUSE.RelRisksSameBtwSerotypes	= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Which_Sero_FitsAll"				)	HOUSE.Which_Sero_FitsAll		= std::stoi(param_value_string);	//// int
			
			else	if (param_name == "LTFU_SurvivalCurves"				)	HOUSE.LTFU_SurvivalCurves				= std::stoi(param_value_string);	//// bool
			else	if (param_name == "OutputStringExtra"				)	HOUSE.OutputStringExtra					= param_value_string;				//// string
			else	if (param_name == "OutputStringExtra_PrevChain"		)	HOUSE.OutputStringExtra_PrevChain		= param_value_string;				//// string
			else	if (param_name == "OutputStringForOldChainInput_Aug")	HOUSE.OutputStringForOldChainInput_Aug	= param_value_string;				//// string
			else	if (param_name == "OldChainFileNamePrefix_Aug"		)	HOUSE.OldChainFileNamePrefix_Aug		= param_value_string;				//// string

			else	if (param_name == "StartFromPreviousChain"			)	HOUSE.StartFromPreviousChain			= std::stoi(param_value_string);	//// bool
			else	if (param_name == "StartFromPreviousChain_Aug"		)	HOUSE.StartFromPreviousChain_Aug		= std::stoi(param_value_string);	//// bool
			else	if (param_name == "AddOutputStringExtraToOldChains"	)	HOUSE.AddOutputStringExtraToOldChains	= std::stoi(param_value_string);	//// bool
			
			else	if (param_name == "ParamRangeFileName"				)	HOUSE.ParamRangeFileName				= param_value_string;				//// string
			else	if (param_name == "WBIC_RunOnly"					)	HOUSE.WBIC_RunOnly						= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Do_WBIC_Runs"					)	HOUSE.Do_WBIC_Runs						= std::stoi(param_value_string);	//// bool
			
			
			else	if (param_name == "ModellingHospitalized"			)	HOUSE.ModellingHospitalized					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "ModelHosp_Indie_Ks"				)	HOUSE.ModelHosp_Indie_Ks					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Fixed_Severe_RelRisks"			)	HOUSE.Fixed_Severe_RelRisks					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "FSKs_Ratio_SetNum"				)	HOUSE.FSKs_Ratio_SetNum						= std::stoi(param_value_string);	//// int
			
			else	if (param_name == "SerotypeString"					)	HOUSE.SerotypeString						= param_value_string;				//// string
			else	if (param_name == "SFU"								)	HOUSE.SFU									= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Include_Late_Cases"				)	HOUSE.Include_Late_Cases					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "FakeExtObs"						)	HOUSE.FakeExtObs							= std::stoi(param_value_string);	//// bool
			else	if (param_name == "PooledCountries"					)	HOUSE.PooledCountries						= std::stoi(param_value_string);	//// bool
			else	if (param_name == "PooledTrials"					)	HOUSE.PooledTrials							= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Empirical_SeroPrevs"				)	HOUSE.Empirical_SeroPrevs					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Weighting_Pass_Sev"				)	HOUSE.Weighting_Pass_Sev					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "PS_Weight"						)	HOUSE.PS_Weight								= std::stod(param_value_string);	//// double 
			else	if (param_name == "ASVE"							)	HOUSE.ASVE									= Convert_AS_String(param_value_string);	//// Age_Option from string
			else	if (param_name == "ASVE_AdditionalKnots"			)	HOUSE.ASVE_AdditionalKnots					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "AS_Haz"							)	HOUSE.AS_Haz								= Convert_AS_String(param_value_string);	//// Age_Option from string
			else	if (param_name == "AS_Haz_AdditionalKnots"			)	HOUSE.AS_Haz_AdditionalKnots				= std::stoi(param_value_string);	//// bool
			else	if (param_name == "AS_Waning"						)	HOUSE.AS_Waning								= Convert_AS_String(param_value_string);	//// Age_Option from string
			else	if (param_name == "AS_Waning_OnlyOneSeroStatus"		)	HOUSE.AS_Waning_OnlyOneSeroStatus			= std::stoi(param_value_string);	//// bool
			else	if (param_name == "AS_Waning_OneSeroBS"				)	HOUSE.AS_Waning_OneSeroBS					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "AS_Waning_KnotSet"				)	HOUSE.AS_Waning_KnotSet						= std::stoi(param_value_string);	//// int
			else	if (param_name == "AS_Waning_Homogeneous"			)	HOUSE.AS_Waning_Homogeneous					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "AgeEffectsSame_Waning"			)	HOUSE.AgeEffectsSame_Waning					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "AgeEffectsSame_VE"				)	HOUSE.AgeEffectsSame_VE						= std::stoi(param_value_string);	//// bool
			
			else	if (param_name == "PSVEs"							)	HOUSE.PSVEs									= std::stoi(param_value_string);	//// bool
			else	if (param_name == "AS_VE_Homogeneous"				)	HOUSE.AS_VE_Homogeneous						= std::stoi(param_value_string);	//// bool
			else	if (param_name == "ASVE_OnlyOneSeroStatus"			)	HOUSE.ASVE_OnlyOneSeroStatus				= std::stoi(param_value_string);	//// bool
			else	if (param_name == "ASVE_BS"							)	HOUSE.ASVE_BS								= std::stoi(param_value_string);	//// int 
			else	if (param_name == "AdjHaz"							)	HOUSE.AdjHaz								= std::stoi(param_value_string);	//// bool 
			else	if (param_name == "Which_BS_BaseHazMult"			)	HOUSE.Which_BS_BaseHazMult					= std::stoi(param_value_string);	//// int 
			
			else	if (param_name == "AllDosesRequired_SNeg"			)	HOUSE.AllDosesRequired_SNeg					= std::stoi(param_value_string);	//// bool 
			else	if (param_name == "AllDosesRequired_SPos"			)	HOUSE.AllDosesRequired_SPos					= std::stoi(param_value_string);	//// bool 
			else	if (param_name == "AllDosesRequired_SNeg_AgeGroup1"	)	HOUSE.AllDosesRequired_SNeg_AgeGroup1		= std::stoi(param_value_string);	//// bool 
			else	if (param_name == "AllDosesRequired_SNeg_AgeGroup2"	)	HOUSE.AllDosesRequired_SNeg_AgeGroup2		= std::stoi(param_value_string);	//// bool 
			else	if (param_name == "AllDosesRequired_SNeg_AgeGroup3"	)	HOUSE.AllDosesRequired_SNeg_AgeGroup3		= std::stoi(param_value_string);	//// bool 
			else	if (param_name == "AllDosesRequired_SNeg_AgeGroup4"	)	HOUSE.AllDosesRequired_SNeg_AgeGroup4		= std::stoi(param_value_string);	//// bool 
			else	if (param_name == "AllDosesRequired_SNeg_AgeGroup5"	)	HOUSE.AllDosesRequired_SNeg_AgeGroup5		= std::stoi(param_value_string);	//// bool 
			else	if (param_name == "AllDosesRequired_SPos_AgeGroup1"	)	HOUSE.AllDosesRequired_SPos_AgeGroup1		= std::stoi(param_value_string);	//// bool 
			else	if (param_name == "AllDosesRequired_SPos_AgeGroup2"	)	HOUSE.AllDosesRequired_SPos_AgeGroup2		= std::stoi(param_value_string);	//// bool 
			else	if (param_name == "AllDosesRequired_SPos_AgeGroup3"	)	HOUSE.AllDosesRequired_SPos_AgeGroup3		= std::stoi(param_value_string);	//// bool 
			else	if (param_name == "AllDosesRequired_SPos_AgeGroup4"	)	HOUSE.AllDosesRequired_SPos_AgeGroup4		= std::stoi(param_value_string);	//// bool 
			else	if (param_name == "AllDosesRequired_SPos_AgeGroup5"	)	HOUSE.AllDosesRequired_SPos_AgeGroup5		= std::stoi(param_value_string);	//// bool 
			else	if (param_name == "ExtImSub"						)	HOUSE.ExtImSub								= Convert_ExtImSub_String(param_value_string);	//// bool 
			else	if (param_name == "AreWeAugmenting"					)	HOUSE.AreWeAugmenting						= std::stoi(param_value_string);	//// bool

			else	if (param_name == "BurnIn"							)	CHAINS.BurnIn									= std::stoi(param_value_string);	//// int
			else	if (param_name == "No_Iterations"					)	CHAINS.No_Iterations							= std::stoi(param_value_string);	//// int
			else	if (param_name == "WBIC_BurnIn"						)	WBIC_CHAINS.BurnIn								= std::stoi(param_value_string);	//// int
			else	if (param_name == "WBIC_No_Iterations"				)	WBIC_CHAINS.No_Iterations						= std::stoi(param_value_string);	//// int
			else	if (param_name == "AugmentEveryHowManyIterations"	)	CHAINS.AugmentEveryHowManyIterations			= std::stoi(param_value_string);	//// int
			else	if (param_name == "AddtoChainEveryHowManyIterations")	CHAINS.AddtoChainEveryHowManyIterations			= std::stoi(param_value_string);	//// int
			else	if (param_name == "WhenToFlipTemperatureBack"		)	CHAINS.WhenToFlipTemperatureBack				= std::stoi(param_value_string);	//// int
			else	if (param_name == "SimulatedAnnealing"				)	CHAINS.SimulatedAnnealing						= std::stoi(param_value_string);	//// bool
			else	if (param_name == "CoolDuringBurnIn"				)	CHAINS.CoolDuringBurnIn							= std::stoi(param_value_string);	//// bool
			else	if (param_name == "FinalTemperature"				)	CHAINS.FinalTemperature							= std::stod(param_value_string);	//// double
			else	if (param_name == "WBIC_FinalTemperature"			)	WBIC_CHAINS.FinalTemperature					= std::stod(param_value_string);	//// double
			else	if (param_name == "WBIC_UseDefault_FinalTemperature")	WBIC_CHAINS.UseDefault_FinalTemperature			= std::stod(param_value_string);	//// double
			
			else	if (param_name == "Run_MCMC"							)	HOUSE.Run_MCMC								= std::stoi(param_value_string);	//// bool
			else	if (param_name == "AreWeFittingParameters"				)	CHAINS.AreWeFittingParameters				= std::stoi(param_value_string);	//// bool
			else	if (param_name == "AreWeChecking"						)	CHAINS.AreWeChecking						= std::stoi(param_value_string);	//// bool
			else	if (param_name == "CalculateFreshLikeLihood"			)	CHAINS.CalculateFreshLikeLihood				= std::stoi(param_value_string);	//// bool
			else	if (param_name == "CheckIndividualAugmentation"			)	CHAINS.CheckIndividualAugmentation			= std::stoi(param_value_string);	//// bool
			else	if (param_name == "AreWeCalculatingSeroPrevs"			)	CHAINS.AreWeCalculatingSeroPrevs			= std::stoi(param_value_string);	//// bool
			else	if (param_name == "HazRatiosOnly"						)	CHAINS.HazRatiosOnly						= std::stoi(param_value_string);	//// bool
			else	if (param_name == "OutputIndividualSurvivalTables"		)	CHAINS.OutputIndividualSurvivalTables		= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Calc_SCsARsHRPs_Any"					)	CHAINS.Calc_SCsARsHRPs_Any					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Calc_SCsARsHRPs_MeanAndCrIs"			)	CHAINS.Calc_SCsARsHRPs_MeanAndCrIs			= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Calc_SCsARsHRPs_ModalMaxLike"		)	CHAINS.Calc_SCsARsHRPs_ModalMaxLike			= std::stoi(param_value_string);	//// bool
			else	if (param_name == "WBIC_Calc_SCsARsHRPs_Any"			)	WBIC_CHAINS.Calc_SCsARsHRPs_Any				= std::stoi(param_value_string);	//// bool
			else	if (param_name == "WBIC_Calc_SCsARsHRPs_MeanAndCrIs"	)	WBIC_CHAINS.Calc_SCsARsHRPs_MeanAndCrIs		= std::stoi(param_value_string);	//// bool
			else	if (param_name == "WBIC_Calc_SCsARsHRPs_ModalMaxLike"	)	WBIC_CHAINS.Calc_SCsARsHRPs_ModalMaxLike	= std::stoi(param_value_string);	//// bool
			
			else	if (param_name == "seed1"							)	HOUSE.seed1							= std::stoi(param_value_string);	//// int
			else	if (param_name == "seed2"							)	HOUSE.seed2							= std::stoi(param_value_string);	//// int
			else	if (param_name == "DataFilename"					)	HOUSE.DataFilename					= param_value_string;				//// string
			else	if (param_name == "UseSyntheticData"				)	HOUSE.UseSyntheticData				= std::stoi(param_value_string);	//// bool
			else	if (param_name == "RandomCases"						)	HOUSE.RandomCases					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "RandCase_Seed"					)	HOUSE.RandCase_Seed					= std::stoi(param_value_string);	//// int
			
			else	if (param_name == "Output_LL_minus_Aug"				)	CHAINS.Output_LL_minus_Aug				= std::stoi(param_value_string);	//// bool
			else	if (param_name == "AugmentEveryHowManyIterations"	)	CHAINS.AugmentEveryHowManyIterations	= std::stoi(param_value_string);	//// int
			else	if (param_name == "AddtoChainEveryHowManyIterations")	CHAINS.AddtoChainEveryHowManyIterations	= std::stoi(param_value_string);	//// int
			
			//// ///// **** **** //// ///// **** **** //// ///// **** **** //// ///// **** **** //// ///// **** **** //// ///// **** **** //// /////  **** **** //// ///// 
			//// //// (initial) model Parameter input
			
			else	if (param_name == "KA_0"								)	CurrentPARAMS.KA_0									= std::stod(param_value_string);	//// double
			else	if (param_name == "KA_1"								)	CurrentPARAMS.KA_1									= std::stod(param_value_string);	//// double
			else	if (param_name == "KA_2"								)	CurrentPARAMS.KA_2									= std::stod(param_value_string);	//// double
			else	if (param_name == "KH_0"								)	CurrentPARAMS.KH_0									= std::stod(param_value_string);	//// double
			else	if (param_name == "KH_1"								)	CurrentPARAMS.KH_1									= std::stod(param_value_string);	//// double
			else	if (param_name == "KH_2"								)	CurrentPARAMS.KH_2									= std::stod(param_value_string);	//// double
			else	if (param_name == "PosEfficacy"							)	CurrentPARAMS.PosEfficacy							= std::stod(param_value_string);	//// double
			else	if (param_name == "NegEfficacy"							)	CurrentPARAMS.NegEfficacy							= std::stod(param_value_string);	//// double
			else	if (param_name == "PosWaning"							)	CurrentPARAMS.PosWaning								= std::stod(param_value_string);	//// double
			else	if (param_name == "NegWaning"							)	CurrentPARAMS.NegWaning								= std::stod(param_value_string);	//// double
			else	if (param_name == "Initial_BS_BaseHazMult"				)	CurrentPARAMS.Initial_BS_BaseHazMult				= std::stod(param_value_string);	//// double
			else	if (param_name == "Hosp_K_mult"							)	CurrentPARAMS.Hosp_K_mult							= std::stod(param_value_string);	//// double
			else	if (param_name == "Initial_Halflife_Hill"				)	CurrentPARAMS.Initial_Halflife_Hill					= std::stod(param_value_string);	//// double
			else	if (param_name == "Initial_Power_Hill"					)	CurrentPARAMS.Initial_Power_Hill					= std::stod(param_value_string);	//// double
			else	if (param_name == "Initial_Halflife_ASVE"				)	CurrentPARAMS.Initial_Halflife_ASVE					= std::stod(param_value_string);	//// double
			else	if (param_name == "Initial_Power_ASVE"					)	CurrentPARAMS.Initial_Power_ASVE					= std::stod(param_value_string);	//// double
			else	if (param_name == "Initial_Prop_ASVE"					)	CurrentPARAMS.Initial_Prop_ASVE						= std::stod(param_value_string);	//// double
			else	if (param_name == "Initial_Halflife_AS_Haz"				)	CurrentPARAMS.Initial_Halflife_AS_Haz				= std::stod(param_value_string);	//// double
			else	if (param_name == "Initial_Power_AS_Haz"				)	CurrentPARAMS.Initial_Power_AS_Haz					= std::stod(param_value_string);	//// double
			else	if (param_name == "Initial_Prop_AS_Haz"					)	CurrentPARAMS.Initial_Prop_AS_Haz					= std::stod(param_value_string);	//// double
			else	if (param_name == "Initial_Halflife_AS_Waning"			)	CurrentPARAMS.Initial_Halflife_AS_Waning			= std::stod(param_value_string);	//// double
			else	if (param_name == "Initial_Power_AS_Waning"				)	CurrentPARAMS.Initial_Power_AS_Waning				= std::stod(param_value_string);	//// double
			else	if (param_name == "Initial_AS_Prime_PropIndptAge_SNeg"	)	CurrentPARAMS.Initial_AS_Prime_PropIndptAge_SNeg	= std::stod(param_value_string);	//// double
			else	if (param_name == "Initial_AS_Prime_PropIndptAge_SPos"	)	CurrentPARAMS.Initial_AS_Prime_PropIndptAge_SPos	= std::stod(param_value_string);	//// double
			else	if (param_name == "Initial_AS_Prime_SNeg_rate"			)	CurrentPARAMS.Initial_AS_Prime_SNeg_rate			= std::stod(param_value_string);	//// double
			else	if (param_name == "Initial_AS_Prime_SPos_rate"			)	CurrentPARAMS.Initial_AS_Prime_SPos_rate			= std::stod(param_value_string);	//// double

			//// ///// **** **** //// ///// **** **** //// ///// **** **** //// ///// **** **** //// ///// **** **** //// ///// **** **** //// /////  **** **** //// ///// 
			//// //// Fit or skip various parameters
			
			else	if (param_name == "FitAllCountries"					)	HOUSE.FitAllCountries					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Fit_c0"							)	HOUSE.Fit_c0							= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Fit_c1"							)	HOUSE.Fit_c1							= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Fit_c2"							)	HOUSE.Fit_c2							= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Fit_c3"							)	HOUSE.Fit_c3							= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Fit_c4"							)	HOUSE.Fit_c4							= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Fit_c5"							)	HOUSE.Fit_c5							= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Fit_c6"							)	HOUSE.Fit_c6							= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Fit_c7"							)	HOUSE.Fit_c7							= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Fit_c8"							)	HOUSE.Fit_c8							= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Fit_c9"							)	HOUSE.Fit_c9							= std::stoi(param_value_string);	//// bool
			else	if (param_name == "FitAllGlobalParams"				)	HOUSE.FitAllGlobalParams				= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_KA_0"						)	HOUSE.Skip_KA_0							= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_NegWaning"					)	HOUSE.Skip_NegWaning					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_NegWaning_HLife"			)	HOUSE.Skip_NegWaning_HLife				= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_NegWaning_Power"			)	HOUSE.Skip_NegWaning_Power				= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_PosWaning"					)	HOUSE.Skip_PosWaning					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_PosWaning_HLife"			)	HOUSE.Skip_PosWaning_HLife				= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_PosWaning_Power"			)	HOUSE.Skip_PosWaning_Power				= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_KA_2"						)	HOUSE.Skip_KA_2							= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_KH_0"						)	HOUSE.Skip_KH_0							= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_KH_1"						)	HOUSE.Skip_KH_1							= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_KH_2"						)	HOUSE.Skip_KH_2							= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_PosEfficacies"				)	HOUSE.Skip_PosEfficacies				= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_NegEfficacies"				)	HOUSE.Skip_NegEfficacies				= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_qvals"						)	HOUSE.Skip_qvals						= std::stoi(param_value_string);	//// bool
			else	if (param_name == "SkipRelativeRisks"				)	HOUSE.SkipRelativeRisks					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "SkipKnots"						)	HOUSE.SkipKnots							= std::stoi(param_value_string);	//// bool
			else	if (param_name == "SkipHistHazards"					)	HOUSE.SkipHistHazards					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_Hosp_Mult"					)	HOUSE.Skip_Hosp_Mult					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_Rhos"						)	HOUSE.Skip_Rhos							= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_All_ASVEs"					)	HOUSE.Skip_All_ASVEs					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_SNeg_ASVEs"					)	HOUSE.Skip_SNeg_ASVEs					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_SPos_ASVEs"					)	HOUSE.Skip_SPos_ASVEs					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_ASVEsPowers"				)	HOUSE.Skip_ASVEsPowers					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_ASVE_Halflives"				)	HOUSE.Skip_ASVE_Halflives				= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_ASVE_Props"					)	HOUSE.Skip_ASVE_Props					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_All_AS_Haz"					)	HOUSE.Skip_All_AS_Haz					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_AS_Haz_Powers"				)	HOUSE.Skip_AS_Haz_Powers				= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_AS_Haz_Halflives"			)	HOUSE.Skip_AS_Haz_Halflives				= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_AS_Haz_Props"				)	HOUSE.Skip_AS_Haz_Props					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_BS_BaseHazMult"				)	HOUSE.Skip_BS_BaseHazMult				= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_AS_Prime_All"				)	HOUSE.Skip_AS_Prime_All					= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_AS_Prime_SNegRate"			)	HOUSE.Skip_AS_Prime_SNegRate			= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_AS_Prime_SPosRate"			)	HOUSE.Skip_AS_Prime_SPosRate			= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_AS_Prime_SNegProp"			)	HOUSE.Skip_AS_Prime_SNegProp			= std::stoi(param_value_string);	//// bool
			else	if (param_name == "Skip_AS_Prime_SPosProp"			)	HOUSE.Skip_AS_Prime_SPosProp			= std::stoi(param_value_string);	//// bool
			else	std::cerr << " ReadInParams ERROR: param_name " << param_name << " not recognized. param_value_string = " << param_value_string << endl;
		}
		ParamsEtc.close();
	}
		 if (ModVariantStringDummy == "SIMPLE_ANALYTICAL"			) HOUSE.ModelVariant			= SIMPLE_ANALYTICAL			;
	else if (ModVariantStringDummy == "SIMPLE_NUMERICAL"			) HOUSE.ModelVariant			= SIMPLE_NUMERICAL			;
	else if (ModVariantStringDummy == "K_SEROPOS"					) HOUSE.ModelVariant			= K_SEROPOS					;
	else if (ModVariantStringDummy == "DROP_K"						) HOUSE.ModelVariant			= DROP_K					;
	else if (ModVariantStringDummy == "VAC_SILENT"					) HOUSE.ModelVariant			= VAC_SILENT				;
	else if (ModVariantStringDummy == "AS_PRIME"					) HOUSE.ModelVariant			= AS_PRIME					;
	else std::cerr << "ModelVariant variable incorrectly specified: using default VAC_SILENT" << endl;

		 if (SingleOrMultiDoseStringDummy == "SINGLE_DOSE"			) HOUSE.SingleOrMultiDose		= SINGLE_DOSE				;
	else if (SingleOrMultiDoseStringDummy == "MULTI_DOSE"			) HOUSE.SingleOrMultiDose		= MULTI_DOSE				;
	else std::cerr << "SingleOrMultiDose variable incorrectly specified: using default MULTI_DOSE" << endl; 

		 if (ActiveOrHospitalStringDummy == "ACTIVE_PHASE_ONLY"		) HOUSE.ActiveOrPassivePhase	= ACTIVE_PHASE_ONLY			;
	else if (ActiveOrHospitalStringDummy == "DO_ACTIVE_AND_PASSIVE"	) HOUSE.ActiveOrPassivePhase	= DO_ACTIVE_AND_PASSIVE		;
	else std::cerr << "ActiveOrPassivePhase variable incorrectly specified: using default ACTIVE_PHASE_ONLY" << endl; 

		 if (MildAndSevereStringDummy == "TREATED_EQUALLY"			) HOUSE.MildAndSevere			= TREATED_EQUALLY			;
	else if (MildAndSevereStringDummy == "TREATED_SEPARATELY"		) HOUSE.MildAndSevere			= TREATED_SEPARATELY		;
	else std::cerr << "MildAndSevere variable incorrectly specified: using default TREATED_EQUALLY" << endl; 

		 if (AugMH_or_GibbsDummy == "GIBBS_AUG"						) HOUSE.Aug_MH_or_Gibbs			= GIBBS_AUG					;
	else if (AugMH_or_GibbsDummy == "MH_AUG"						) HOUSE.Aug_MH_or_Gibbs			= MH_AUG					;
	else std::cerr << "AugMH_or_GibbsDummy variable incorrectly specified: using default GIBBS_AUG" << endl; 

	std::cerr << " DONE " << endl; std::fflush(stderr); 

#endif
}
void InitializeDATApointers	(DATA_struct &DATA, const Housekeeping_Struct &HOUSE)
{
	DATA.FollowUp							= new int  **[2]();	//// start and finish. 
	DATA.FollowUp[START]					= new int  * [2]();	
	DATA.FollowUp[END	]					= new int  * [2]();
	DATA.IsCase_AMandPS						= new bool * [2]();

	DATA.FollowUp[START]	[ActiveMild	]	= new int	[NPat]();
	DATA.FollowUp[END	]	[ActiveMild	]	= new int	[NPat]();

	DATA.IsCase_AMandPS		[ActiveMild	]		= DATA.IsCase_ActiveMild		;

	DATA.FollowUp[START]	[PassiveSevere	]		= new int	[NPat]();
	DATA.FollowUp[END	]	[PassiveSevere	]		= new int	[NPat]();

	DATA.IsCase_AMandPS		[PassiveSevere	]		= DATA.IsCase_PassiveSevere		;
}
void ReadInData				(DATA_struct &DATA, Housekeeping_Struct &HOUSE)
{
	ifstream inFile;
	inFile.open(HOUSE.DataFilename);

	/// Need minimum and maximum calendar times period for each country to prevent  calculating and storing (integrated) hazard values and waning values for calendar times before first dose of first patient in that country and after latest end of follow up in that country. 
	for (int country = 0; country < NumC; country++)	
	{	
		DATA.CountryMinMaxCalendarTime[country]					= new DType[2]();		
		DATA.CountryMinMaxCalendarTime[country][Min_TimeIndex]	= 60000; 		// index 0 = minimum. Pick arbitrarily large number definitely greater than minimum calendar time. 
		DATA.CountryMinMaxCalendarTime[country][Max_TimeIndex]	= 0;			// index 1 = maximum. 
	}

	bool *ptr_IsCase = NULL; //// if only doing active phase, want cases in active. Otherwise want cases in active and passive. All other parameters the same. 
			if (HOUSE.PASSIVE_PHASE_ONLY)			ptr_IsCase		= DATA.IsCase_PassiveSevere	;
	else	if (HOUSE.HowManyCaseCategories == 1)	ptr_IsCase		= DATA.IsCase_ActiveMild	;
	else	if (HOUSE.HowManyCaseCategories == 2)	ptr_IsCase		= DATA.IsCase				;
	
	if (HOUSE.PASSIVE_PHASE_ONLY & HOUSE.MildAndSevere == TREATED_SEPARATELY)
		for (int line = 0; line < 100; line++) std::cerr << "PASSIVE_PHASE_ONLY & MILDSEVERE - NOT CODED UP!!!" << endl;

	int *ptr_Iis = (HOUSE.ExtImSub == ExtImSub_Option::AS_DATA) ? DATA.Imputed_Ii_s : DATA.Ii_s;
	int *ptr_EndFollowUp = NULL;
			if (HOUSE.PASSIVE_PHASE_ONLY)								ptr_EndFollowUp = DATA.FollowUp[END][PassiveSevere	];
	else	if (HOUSE.ActiveOrPassivePhase == ACTIVE_PHASE_ONLY		)	ptr_EndFollowUp = DATA.FollowUp[END][ActiveMild		];	
	else	if (HOUSE.ActiveOrPassivePhase == DO_ACTIVE_AND_PASSIVE	)	ptr_EndFollowUp = DATA.FollowUp[END][PassiveSevere	];	
	int *ptr_StartFollowUp = (HOUSE.PASSIVE_PHASE_ONLY) ? DATA.FollowUp[START][PassiveSevere] : DATA.FollowUp[START][ActiveMild];

	DATA.CaseSerotype = new int[NPat](); //// even when not considering serotype fill this with zeros
	if ((HOUSE.SeroSpecificEfficacies) || (HOUSE.SeroSpecific_K_values)) DATA.NoDataVariables++;		

	std::vector<string> DataVariableNames(DATA.NoDataVariables);
	for (int variable = 0; variable < DATA.NoDataVariables; variable++) inFile >> DataVariableNames[variable]; /*printVector(DataVariableNames, "DataVariableNames"); */

	for (int a = 0; a < NPat; a++)
	{
		//// add to Patient index number. 
		DATA.PatientIndexNumbers.push_back(a); 

		inFile >>	DATA.FollowUp[START][ActiveMild][a] >> DATA.SecondDose[a] >> DATA.ThirdDose[a] >> DATA.FollowUp[END][ActiveMild][a] >> 
					DATA.FollowUp[START][PassiveSevere][a] >> DATA.FollowUp[END][PassiveSevere][a] >> DATA.TimePost_1st_Dose[a] >> DATA.TimePost_Final_Dose[a] >>
					DATA.Ii_s[a] >> DATA.Vi_s[a] >> DATA.ai_s[a] >> DATA.AgeGroup1[a] >> DATA.AgeGroup2[a] >> DATA.ci_s[a] >> 
			DATA.IsCase[a] >> DATA.IsCase_ActiveMild[a] >> DATA.IsCase_PassiveSevere[a] >> DATA.Imputed_Ii_s[a] >> DATA.Imputed_ProbSPos[a];
		
		inFile >> DATA.EndActiveSurveillance	[a];
		inFile >> DATA.EndPassiveSurveillance	[a];
		inFile >> DATA.LTFU						[a];
		if ((HOUSE.SeroSpecificEfficacies) || (HOUSE.SeroSpecific_K_values)) inFile >> DATA.CaseSerotype[a];

		///// If aggregating countries, change data here. Must be done before CountryMinMaxCalendarTime and lambda stuff below 
		if (HOUSE.PooledCountries)
			if (DATA.ci_s[a] <= 4)		DATA.ci_s[a] = HOUSE.Combined_CYD_14_country;
			else						DATA.ci_s[a] = HOUSE.Combined_CYD_15_country;
		else if (HOUSE.PooledTrials)	DATA.ci_s[a] = HOUSE.Combined_CYD_14_15_country;

		if (std::any_of(HOUSE.WhichCountries.begin(), HOUSE.WhichCountries.end(), [&](int i){return i == DATA.ci_s[a]; }))
		{
			if (ptr_Iis[a] == MDvalue)	{	DATA.AugmentedIndices.push_back(a)		; DATA.NoAugmented++;	}	// store indices of		augmented data patients. And add to the number of augmented patients. 
			else 						{	DATA.NonAugmentedIndices.push_back(a)	; DATA.InImmSub[a] = 1;	}	// store indices of non-augmented data patients. And record that this patient is in the ImmuneSubset. 
			if (ptr_IsCase[a] == 1) DATA.NoInfected++;
			HOUSE.TotalPatients++; //// add to number of people in this model run. 
		}

		if (HOUSE.PASSIVE_PHASE_ONLY) //// when setting the minimum calendar time that a patient entered passive phase, need to be careful not to accidentally include active phase cases or LTFUs, as their "Start (and end) Passive Phase" has been set to the last day of their active phase, to make sure any integrated hazards are zero. But here would erroneously result in early passivee phase. 
		{
			if (ptr_StartFollowUp[a]	< DATA.CountryMinMaxCalendarTime[DATA.ci_s[a]][Min_TimeIndex] & !DATA.IsCase_ActiveMild[a] & !DATA.LTFU[a] ) DATA.CountryMinMaxCalendarTime[DATA.ci_s[a]][Min_TimeIndex] = ptr_StartFollowUp[a];
		}
		else
		{
			if (ptr_StartFollowUp[a]	< DATA.CountryMinMaxCalendarTime[DATA.ci_s[a]][Min_TimeIndex]) DATA.CountryMinMaxCalendarTime[DATA.ci_s[a]][Min_TimeIndex] = ptr_StartFollowUp[a];
		}
		if (ptr_EndFollowUp		[a]	> DATA.CountryMinMaxCalendarTime[DATA.ci_s[a]][Max_TimeIndex]) DATA.CountryMinMaxCalendarTime[DATA.ci_s[a]][Max_TimeIndex] = ptr_EndFollowUp	[a];

		if (ptr_StartFollowUp	[a] > DATA.MaxStartFollowUp_CalendarTime) DATA.MaxStartFollowUp_CalendarTime	= ptr_StartFollowUp[a];
		if (ptr_EndFollowUp		[a]	> DATA.MaxEndFollowUp_CalendarTime	) DATA.MaxEndFollowUp_CalendarTime		= ptr_EndFollowUp[a];
	}
	std::cerr << "DATA.NoInfected " << DATA.NoInfected << " HOUSE.TotalPatients " << HOUSE.TotalPatients << endl;


	if (HOUSE.LTFU_SurvivalCurves && HOUSE.Use_WithoutCases_FU_Defns)
	{
		DATA.EndActiveSurveillance	= DATA.FollowUp[END][ActiveMild		];
		DATA.EndPassiveSurveillance = DATA.FollowUp[END][PassiveSevere	];


		std::cerr << " Use_WithoutCases_FU_Defns! " << endl;
		std::cerr << " Use_WithoutCases_FU_Defns! " << endl;
		std::cerr << " Use_WithoutCases_FU_Defns! " << endl;
		std::cerr << " Use_WithoutCases_FU_Defns! " << endl;
		std::cerr << " Use_WithoutCases_FU_Defns! " << endl;
		std::cerr << " Use_WithoutCases_FU_Defns! " << endl;
	}

	std::cerr << " MaxEndFollowUp_CalendarTime " << DATA.MaxEndFollowUp_CalendarTime << endl;
	inFile.close();
}
void ProcessData			(DATA_struct &DATA, const Housekeeping_Struct &HOUSE)
{
	if ((HOUSE.SingleOrMultiDose == SINGLE_DOSE) || (HOUSE.ModelVariant == SIMPLE_ANALYTICAL)) for (int a = 0; a < NPat; a++) DATA.TimePost_Final_Dose[a] = DATA.TimePost_1st_Dose[a];  //// Easier to simply redefine TimePostDose data than to have loads of if statements. You do this because otherwise Integrated Vaccine hazards wouldn't agree with case probability densities in SIMPLE_ANALYTICAL

	//// redefine follow up dates for PASSIVE_PHASE_ONLY
	if (HOUSE.PASSIVE_PHASE_ONLY)
	{
		//// Need to find the earliest calendar date that somebody started their passive phase (don't include LTFU patients or active phase cases)
		int EarliestPassivePhaseStart_CalendarDay = 90000; //// i.e. something huge
		for (int a = 0; a < NPat; a++)
			if (!DATA.IsCase_ActiveMild[a] && !DATA.LTFU[a])
				if (DATA.FollowUp[START][PassiveSevere][a] < EarliestPassivePhaseStart_CalendarDay)
					EarliestPassivePhaseStart_CalendarDay = DATA.FollowUp[START][PassiveSevere][a]; 

		//// redefine passive phase calendar times and all derived quantities (effectively resetting "Day_0" that you use in R script)
		for (int a = 0; a < NPat; a++)	
		{
			DATA.FollowUp[START	][PassiveSevere][a] -= EarliestPassivePhaseStart_CalendarDay;
			DATA.FollowUp[END	][PassiveSevere][a] -= EarliestPassivePhaseStart_CalendarDay;
			DATA.EndPassiveSurveillance			[a]	-= EarliestPassivePhaseStart_CalendarDay;
			DATA.SecondDose						[a] -= EarliestPassivePhaseStart_CalendarDay; //// need to reset second and third dose date so that the time differences between values of t and dose date are the same. 
			DATA.ThirdDose						[a] -= EarliestPassivePhaseStart_CalendarDay;
		}

		for (int country = 0; country < HOUSE.TotalCountries; country++)
		{
			DATA.CountryMinMaxCalendarTime[country][Min_TimeIndex] -= EarliestPassivePhaseStart_CalendarDay; 
			DATA.CountryMinMaxCalendarTime[country][Max_TimeIndex] -= EarliestPassivePhaseStart_CalendarDay;
		}

		DATA.MaxStartFollowUp_CalendarTime	-= EarliestPassivePhaseStart_CalendarDay;
		DATA.MaxEndFollowUp_CalendarTime	-= EarliestPassivePhaseStart_CalendarDay;

		for (int a = 0; a < NPat; a++)	
		{
			if (!DATA.IsCase_ActiveMild[a] && !DATA.LTFU[a]) //// For all patients who where in passive phase, redefine "active" phase to equal passive. 
			{
				DATA.FollowUp[START	][ActiveMild][a]	= DATA.FollowUp[START	][PassiveSevere][a];
				DATA.FollowUp[END	][ActiveMild][a]	= DATA.FollowUp[END		][PassiveSevere][a];
				DATA.EndActiveSurveillance[a]			= DATA.EndPassiveSurveillance[a];
			}
			else //// for all active phase cases and LTFUs, redefine active phase to start and end at zero
			{
				DATA.FollowUp[START	][ActiveMild][a]	= 0;
				DATA.FollowUp[END	][ActiveMild][a]	= 0;
				DATA.EndActiveSurveillance[a]			= 0;
			}
		}

		///// Now redefine / relabel.  cases.
		for (int a = 0; a < NPat; a++)
		{
			DATA.IsCase[a]					= DATA.IsCase_PassiveSevere[a]; 
			DATA.IsCase_ActiveMild[a]		= DATA.IsCase_PassiveSevere[a];
			DATA.IsCase_PassiveSevere[a]	= 0; 
		}
	}
	else if (HOUSE.HowManyCaseCategories == 1)	 //// i.e. if HOUSE.ActiveOrPassivePhase == ACTIVE_PHASE_ONLY	AND		HOUSE.MildAndSevere == TREATED_EQUALLY. 
		for (int a = 0; a < NPat; a++)
		{
			DATA.IsCase					[a] = DATA.IsCase_ActiveMild[a];
			DATA.IsCase_PassiveSevere	[a] = 0; 
		}

	for (int patient = 0; patient < NPat; patient++)
	{
			 if (DATA.IsCase_AMandPS[ActiveMild		][patient]) DATA.Case_PhaseSeverity[patient] = ActiveMild	; 
		else if (DATA.IsCase_AMandPS[PassiveSevere	][patient]) DATA.Case_PhaseSeverity[patient] = PassiveSevere;
		else													DATA.Case_PhaseSeverity[patient] = MDvalue		; //// intended to cause a crash if you use DATA.Case_PhaseSeverity for non-cases. 
	}

	if (HOUSE.MildAndSevere == TREATED_SEPARATELY)
	{
		for (int patient = 0; patient < NPat; patient++) //// here you consider the entire trial, but need to evaluate two integrated hazards. Easiest to simply change the data. 
			DATA.FollowUp[START][PassiveSevere][patient] = DATA.FollowUp[START][ActiveMild][patient];
		
		if (HOUSE.ActiveOrPassivePhase == ACTIVE_PHASE_ONLY)
			for (int patient = 0; patient < NPat; patient++) //// here you consider ACTIVE_PHASE_ONLY, but need still to evaluate two integrated hazards. Easiest to simply change the data. 
				DATA.FollowUp[END][PassiveSevere][patient] = DATA.FollowUp[END][ActiveMild][patient];
	}
	 
	///// **************** ///// If doing random serotypes for cases (sanity check), randomly assign serotype. Must be done after IsCase adjustment above. 
	srand(HOUSE.RandCase_Seed);
	if (HOUSE.SeroSpecificEfficacies || HOUSE.SeroSpecific_K_values)
		if (HOUSE.UseSyntheticData && HOUSE.RandomCases)
			for (int patient = 0; patient < NPat; patient++)
				if (DATA.IsCase[patient] == 1) DATA.CaseSerotype[patient] = rand() % HOUSE.N_STypes;

	for (int country = 0; country < NumC; country++)
	{
		DATA.CountryMinMaxCalendarTime[country][Min_TimeIndex] = 0;
		DATA.CountryMinMaxCalendarTime[country][Max_TimeIndex] = DATA.MaxEndFollowUp_CalendarTime;
	}
	for (int country = 0; country < NumC; country++) 
	{ 
		DATA.CountryMinMaxCalendarTime[country][Min_TimeIndex] /= 365; 
		DATA.CountryMinMaxCalendarTime[country][Max_TimeIndex] /= 365; 
	}

	for (int country = 0; country < NumC; country++) DATA.NumCalendarDaysFollowUp[country] = (int)round((HOUSE.TimeIntervalReciprocal) * (DATA.CountryMinMaxCalendarTime[country][Max_TimeIndex] - DATA.CountryMinMaxCalendarTime[country][Min_TimeIndex]));

	//// Initialize durations (only need to do this once, rather than every survival curve). 
	Allocate_2D_Array(DATA.FU_Duration_days	, HOUSE.HowManyTrialPhases + 1, NPat);
	Allocate_2D_Array(DATA.FU_Duration_years, HOUSE.HowManyTrialPhases + 1, NPat);

	for (int patient = 0; patient < NPat; patient++)
	{
		DATA.FU_Duration_days	[ActivePhase	][patient]	= DATA.EndActiveSurveillance			[patient] - DATA.FollowUp[START][ActiveMild][patient];
		DATA.FU_Duration_years	[ActivePhase	][patient]	= DATA.FU_Duration_days[ActivePhase	]	[patient] * HOUSE.TimeInterval; 

		if (HOUSE.ActiveOrPassivePhase == DO_ACTIVE_AND_PASSIVE)
		{
			DATA.FU_Duration_days	[PassivePhase	][patient] = DATA.EndPassiveSurveillance		[patient] - DATA.EndActiveSurveillance[patient];
			DATA.FU_Duration_years	[PassivePhase	][patient] = DATA.FU_Duration_days[PassivePhase][patient] * HOUSE.TimeInterval;
			DATA.FU_Duration_days	[WholeTrial		][patient] = DATA.EndPassiveSurveillance		[patient] - DATA.FollowUp[START][ActiveMild][patient];
			DATA.FU_Duration_years	[WholeTrial		][patient] = DATA.FU_Duration_days[WholeTrial]	[patient] * HOUSE.TimeInterval;
		}
	}
#ifdef PRINT_PROGRAM_PROGRESS
	std::cout << " ProcessData DONE" << endl;
#endif

	if (HOUSE.ExtImSub == ExtImSub_Option::AS_DATA) DATA.Ii_s = DATA.Imputed_Ii_s; //// i.e. just change the pointer
}
bool CheckCaseDefinitions	(const DATA_struct &DATA, const Housekeeping_Struct &HOUSE)
{
	bool All_Good = 1;

	if (HOUSE.FitAllCountries == 1)
	{
		int TotCases_def1 = 0;	for (int patient = 0; patient < NPat; patient++) TotCases_def1 += DATA.IsCase[patient];
		int TotCases_def2 = 0;	for (int patient = 0; patient < NPat; patient++) TotCases_def2 += DATA.IsCase_ActiveMild[patient] + DATA.IsCase_PassiveSevere[patient];
		int TotCases_def3 = 0;	for (int patient = 0; patient < NPat; patient++) TotCases_def3 += DATA.IsCase_AMandPS[ActiveMild][patient] + DATA.IsCase_AMandPS[PassiveSevere][patient];

		bool CaseDefs_Okay_1 = 1;
		bool CaseDefs_Okay_2 = 1;
		bool CaseDefs_Okay_3 = 1;

		if (TotCases_def1 != DATA.NoInfected) { std::cerr << "TotCases_def1 != DATA.NoInfected" << endl; CaseDefs_Okay_1 = 0; }
		if (TotCases_def2 != DATA.NoInfected) { std::cerr << "TotCases_def2 != DATA.NoInfected" << endl; CaseDefs_Okay_2 = 0; }
		if (TotCases_def3 != DATA.NoInfected) { std::cerr << "TotCases_def3 != DATA.NoInfected" << endl; CaseDefs_Okay_3 = 0; }

		bool CaseDefs_Consistent_1 = 1;
		bool CaseDefs_Consistent_2 = 1;
		bool CaseDefs_Consistent_3 = 1;

		for (int patient = 0; patient < NPat; patient++) //// this check same as below just using different pointers. Should be identical. 
			if (DATA.IsCase_ActiveMild[patient] || DATA.IsCase_PassiveSevere[patient]) //// i.e. patient should be a case
				if (DATA.IsCase[patient] == 0) /* i.e. if patient not recorded as a case*/ { std::cerr << "Cases definitions inconsistent: patient " << patient << endl; CaseDefs_Consistent_1 = 0; }
		for (int patient = 0; patient < NPat; patient++) //// this check same as above just using different pointers. Should be identical. 
			if (DATA.IsCase_AMandPS[ActiveMild][patient] || DATA.IsCase_AMandPS[PassiveSevere][patient]) //// i.e. patient should be a case
				if (DATA.IsCase[patient] == 0) /* i.e. if patient not recorded as a case*/ { std::cerr << "Cases definitions inconsistent: patient " << patient << endl; CaseDefs_Consistent_2 = 0; }
		if (HOUSE.SeroSpecificEfficacies || HOUSE.SeroSpecific_K_values)
			for (int patient = 0; patient < NPat; patient++)
				if (DATA.IsCase_AMandPS[ActiveMild][patient] || DATA.IsCase_AMandPS[PassiveSevere][patient])
					if (DATA.CaseSerotype[patient] == MDvalue) { std::cerr << "Case Serotype missing for case: patient " << patient << endl; CaseDefs_Consistent_3 = 0; }

		if (!CaseDefs_Okay_1 || !CaseDefs_Okay_1 || !CaseDefs_Okay_1 || !CaseDefs_Consistent_1 || !CaseDefs_Consistent_2 || !CaseDefs_Consistent_3) All_Good = 0;
	}
	return All_Good;
}

std::vector<DType>			ReadOldParamChain(string ParamChainFileName, int NoParams, bool PrintEverything, bool JustFinalChain, DType ** ChainContainerToAddTo)
{
	//// structure of this function is to read entire old parameter chain, and spit out just the final posterior sample as "new" ParamVec. 
	std::ifstream OldParamChain;
	OldParamChain.open(ParamChainFileName);

	if (PrintEverything) std::cout << "ParamChainFileName " << ParamChainFileName << endl;

	string LL_Name_string = "";
	std::vector<string> ParamNames				(NoParams);
	std::vector<string> ParamVec_stringVersion	(NoParams);

	if (PrintEverything) std::cerr << endl << "Testing ReadOldParams" << endl;

	//////// NAMES 
	//// Get LL "name" first
	std::getline(OldParamChain, LL_Name_string, '\t');
	//// Get other param names first
	for (int param_no = 0; param_no < NoParams; param_no++)
		if (param_no == (NoParams - 1)) std::getline(OldParamChain, ParamNames[param_no], '\n'); else std::getline(OldParamChain, ParamNames[param_no], '\t');
	if (PrintEverything) for (int param_no = 0; param_no < NoParams; param_no++) std::cout << " param_name " << ParamNames[param_no] << endl;

	if (PrintEverything) std::cout << "getting values " << endl;

	string Old_LL_Value_stringVersion;
	int postsample = 0;
	while (!OldParamChain.eof())
	{
		if (PrintEverything) std::cerr << " postsample " << postsample << endl;
		//// Get LL first
		std::getline(OldParamChain, Old_LL_Value_stringVersion, '\t');
		if (PrintEverything) std::cout << "Old_LL_Value_stringVersion " << Old_LL_Value_stringVersion << endl;
		for (int param_no = 0; param_no < NoParams; param_no++)
		{
			if (param_no == (NoParams - 1)) std::getline(OldParamChain, ParamVec_stringVersion[param_no], '\n'); else std::getline(OldParamChain, ParamVec_stringVersion[param_no], '\t');
			if (!JustFinalChain)
				ChainContainerToAddTo[param_no][postsample] = std::stod(ParamVec_stringVersion[param_no]);
		}
		postsample++;
		if (PrintEverything) if (OldParamChain.eof()) std::cout << endl << "End of ParamChain" << endl;
	}
	OldParamChain.close();

	std::vector<DType>	ParamVec(NoParams);
	for (int param_no = 0; param_no < NoParams; param_no++)
		ParamVec[param_no] = std::stod(ParamVec_stringVersion[param_no]);
	if (PrintEverything)	std::cout << endl << "ReadOldParamChain FINISHED" << endl;

	return ParamVec;
}
std::vector<DType>			ReadOldParamChain(string ParamChainFileName, int NoParams, bool PrintEverything)
{
	DType ** DummyPointerThatDoesNothingAsYoureApplyingDefaults = NULL; 
	std::vector<DType>	ParamVec = ReadOldParamChain(ParamChainFileName, NoParams, PrintEverything, /*JustFinalChain =*/ true, DummyPointerThatDoesNothingAsYoureApplyingDefaults);
	return ParamVec; 
}

void ProcessOldParamChain	(DATA_struct &DATA, const Housekeeping_Struct &HOUSE, Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS)
{
	if (HOUSE.StartFromPreviousChain)
	{
		string OldParamChainFileName		= "Output\\ParameterChainOutput"		+ HOUSE.OutputStringForOldChainInput + ".txt";
		string WBIC_OldParamChainFileName	= "Output\\WBIC_ParameterChainOutput"	+ HOUSE.OutputStringForOldChainInput + ".txt";

		if (FileExists(WBIC_OldParamChainFileName) & HOUSE.StartFrom_WBIC_Chain) OldParamChainFileName = WBIC_OldParamChainFileName;
		if (FileExists(OldParamChainFileName))
		{
			std::cerr << "Amending Params with " << OldParamChainFileName << endl;
			std::vector<DType> VectorToSeedParamsWith = ReadOldParamChain(OldParamChainFileName, HOUSE.No_Parameters, /*PrintEverything = */false);
			if (VectorToSeedParamsWith.size() == 0) std::cerr << "ProcessOldParamChain error: VectorToSeedParamsWith empty" << std::endl;
			AmendParams(VectorToSeedParamsWith, CurrentPARAMS , DATA, HOUSE);
			AmendParams(VectorToSeedParamsWith, ProposedPARAMS, DATA, HOUSE);
			if (!FileExists("Output\\FinalState_AugData" + HOUSE.OutputStringForOldChainInput + ".txt")) //// InitializeAugmentedData already takes account of old chain. So this loop only relevant if FinalState_AugData doesn't exist.  
			{
				////// Your ParamVec etc. will should now be correct. However, your Augmented data, your Set_Array etc. are based on old Historical hazard parameters. Therefore need to do again, having reset your sets. Not true if final state of the augmented data already imported from an old chain. 
				InitAugDataValues	(DATA, CurrentPARAMS, HOUSE);
				ClearSets			(DATA, HOUSE);		//// Clear Set_Array in Data
				PopulateSets		(DATA, HOUSE);		//// Populate Set_Array in Data
				AmendParams			(VectorToSeedParamsWith, CurrentPARAMS , DATA, HOUSE);
				AmendParams			(VectorToSeedParamsWith, ProposedPARAMS, DATA, HOUSE);
			}
			else std::cerr << "Output\\FinalState_AugData" + HOUSE.OutputStringForOldChainInput + ".txt" << " does not exist" << endl;
		}
		else std::cerr << " OldParamChainFileName: " << OldParamChainFileName << " does not exist" << endl; 
	}
}

void ImportOldParamChains(const Housekeeping_Struct &HOUSE, Chains_Struct &CHAINS, Chains_Struct &WBIC_CHAINS)
{
	std::cerr << "ImportOldParamChains ";
	string OldParamChainFileName		= "Output\\ParameterChainOutput"		+ HOUSE.OutputStringForOldChainInput + ".txt";
	string WBIC_OldParamChainFileName	= "Output\\WBIC_ParameterChainOutput"	+ HOUSE.OutputStringForOldChainInput + ".txt";

	///// differs from the above in that you keep full chains. 
	if (FileExists(OldParamChainFileName		))	
		std::vector<DType> DummyOldParamsVec		= ReadOldParamChain(OldParamChainFileName		, HOUSE.No_Parameters, /*PrintEverything*/ false, /*JustFinalChain*/ false, CHAINS.		ParamChain);
	if (FileExists(WBIC_OldParamChainFileName	))	
		std::vector<DType> DummyOldWBICParamsVec	= ReadOldParamChain(WBIC_OldParamChainFileName	, HOUSE.No_Parameters, /*PrintEverything*/ false, /*JustFinalChain*/ false, WBIC_CHAINS.ParamChain);

	std::cerr << "DONE " << endl;
}


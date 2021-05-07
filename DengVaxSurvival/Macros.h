#ifndef MACROS_HEADER
#define MACROS_HEADER

using namespace std;

//#define USE_CLUSTER			//////////////// REMEMBER TO SET COMPILER TO INTEL
//#define USE_COMMAND_LINE


/// FLAGS AND INDICES.
#define MDvalue		9999			// Number ascribed to missing data. 
#define NPat		31126			// Number of study participants. 
#define NumC		10				// Number of countries. 

#define SIMPLE_ANALYTICAL				0	// numbers arbitrary - but must be different to each other and less than or equal to 255 as ModelVariant is a char. 
#define SIMPLE_NUMERICAL				1
#define K_SEROPOS						2
#define DROP_K							3	
#define VAC_SILENT						4
#define AS_PRIME						5

#define MULTI_DOSE						0
#define SINGLE_DOSE						1

#define TREATED_EQUALLY					0
#define TREATED_SEPARATELY				1

#define ACTIVE_PHASE_ONLY				0	//// these are distinct from Active and Hospital hash-defined below. These refer to whether or not you'll even consider the hospital phase. Active and Hospital (in lower case) refer to switches in your survivor and case functions, when you are including the hospital phase (if you include hospital you still include active)
#define DO_ACTIVE_AND_PASSIVE			1	//// This mostly refers to duration of trial, but also refers to whether or not you calculate certain hazards. 

#define START							0
#define END								1

#define ActiveMild						0	//// Numbers not aribtrary. You will always do ActiveMild but not always PassiveSevere - so when looping ActiveMild has to come first.
#define PassiveSevere					1	

#define ActivePhase						0	//// Could use ActiveMild and PassiveSevere, but for survival curve best not to conflate trial phase with disease severity.  
#define PassivePhase					1
#define WholeTrial						2	//// For follow up durations used in survival curves

#define EITHER							0	//// For survival curve calculation
#define MILD_NON_HOSP					1	//// For survival curve calculation
#define SEVERE_HOSP						2	//// For survival curve calculation

#define MEAN_POST						0
#define LOWER_CrI_POST					1
#define UPPER_CrI_POST					2
#define MODAL_POST						1		//// attack rate posterior summaries will have Mean posterior (i.e. ARs calculated using param vector which is the means of all parameters). Tempting to use this MODE_POST macro in the final S Curves, but that won't be as easy to implement as simply outputting the mean as an individual post sample of the scurves. 
#define MAX_LIKE						2		//// attack rate posterior summaries will have Mean posterior (i.e. ARs calculated using param vector which is the means of all parameters). Tempting to use this MODE_POST macro in the final S Curves, but that won't be as easy to implement as simply outputting the mean as an individual post sample of the scurves. 

#define NonAugIndex						0		//// for seroprev calculations. 
#define AugIndex						1
#define AugAndNonAugIndex				2

#define NonLogIndex						0 	//// for seroprev parameter arrays. Will usually be governed by historical hazards, otherwise by empirical seroprevalence. 
#define LogIndex						1

#define SeroNeg							0	/// NUMBERS Not arbitrary, must correspond to baseline serostatus in Data.Ii_s. 
#define SeroPos							1	/// NUMBERS Not arbitrary, must correspond to baseline serostatus in Data.Ii_s. 
#define EitherSeroStatus				2	/// NUMBERS Not arbitrary, must correspond to baseline serostatus in Data.Ii_s. 

#define ControlGroup					0	/// NUMBERS Not arbitrary, must correspond to Trial Arm in Data.Vi_s. 
#define VaccineGroup					1	/// NUMBERS Not arbitrary, must correspond to Trial Arm in Data.Vi_s. 
#define EitherTrialArm					2	

#define No_Effs							0	//// used for MetaKplusvalues. See explanation in Initialize_KplusValues function, or in Params_Struct definition. Tempting to just use ControlGroup and VaccineGroup macros above but then would need the ControlGroup macro for patients in the vaccine group so best to reinvent the wheel. 
#define With_Effs						1

#define WaningRateDuration				0	//// Important that this one remains zero - as if not fitting age specific waning, then only one waning param per serostatus, hence must refer to 0th parameter in arrays. 
#define Halflife_Index					1	//// Indexes are not arbitrary - you add to ParamVec in a specific order so don't be tempted to change them, or at least not without great care. 
#define Power_Index						2
#define Prop_Index						0	//// Proportion of Vaccine Efficacy that is dependent upon age. Also used for HOUSE.ModelVariant == AS_PRIME so keep as zero. 
#define AS_Priming_Rate_index			1	//// Used for HOUSE.ModelVariant == AS_PRIME. 

#define LowerTail_Index					0	//// for survival curve tails. Not arbitrary. 
#define UpperTail_Index					1

#define LowerBound						0	///// for ParamRanges
#define UpperBound						1

#define Min_TimeIndex					0	///// for CountryMinMaxCalendarTime (in DATA structure). 
#define Max_TimeIndex					1

#define GIBBS_AUG						0	//// gibbs sampler
#define MH_AUG							1	//// metropolis-hastings sampler

#define Serotype_1						0	//// useful for indexing
#define Serotype_2						1
#define Serotype_3						2
#define Serotype_4						3

///// Survival Curve Strata numbers. Chosen to match legacy code.
#define SCurve_Seropositive_Control		0
#define SCurve_Seronegative_Control 	1
#define SCurve_Seropositive_Vaccine 	2
#define SCurve_Seronegative_Vaccine 	3
#define SCurve_Either_Control 			4
#define SCurve_Either_Vaccine 			5
#define SCurve_Seropositive_Either 		6
#define SCurve_Seronegative_Either 		7
#define SCurve_Either_Either 			8



/// VARIOUS CHECKS AND CONSOLE PRINTING FLAGS

//#define CHECK_STRUCTURES_EQUAL
#ifdef USE_CLUSTER				//// to avoid horrible situation where you accidentally leave off random numbers while running on cluster.
#define USE_RANDOM_NUMBERS
#else
#define USE_RANDOM_NUMBERS		//// can be helpful to do without random numbers when debugging. Obviously nothing will work without random numbers when running properly 
#ifndef USE_RANDOM_NUMBERS
#define NON_RANDOM_NUBMER_USED_AS_DUMMY		0.5 //// setting equal to 1 => parameters always accepted; setting equal to 0 (or 0.00001 to avoid log(0)) parameters always rejected
#endif
#endif
//#define VACHAZLIKE_SLOW_WAY
//#define COPY_ENTIRE_PARAM_STRUCTS
//#define PRINT_WITHIN_LIKE_FUNCTIONS
//#define PRINT_WITHIN_CHANGE_FUNCTIONS
#ifdef PRINT_WITHIN_CHANGE_FUNCTIONS
//#define PRINT_WITHIN_CHANGE_KNOT_FUNCTIONS
//#define PRINT_WITHIN_CHANGE_HH_FUNCTIONS
#endif

//#define PRINT_PROGRAM_PROGRESS
#define MaxParamNoToPrint		100
//#define OUTPUT_L_INDICES
//#define OUTPUT_FITTED_PARAM_NOS
//#define PRINT_LIKE_AT_EVERY_ITERATION
//#define PRINT_LOADS
#define ITERATION_TO_PRINT_AFTER	0

//#define PRINT_PARAM_CONDITION true //// i.e. print all parameters
//#define PRINT_PARAM_CONDITION CurrentPARAMS.NamesOfParameters[param_no] == "KAM_0" 
//#define PRINT_PARAM_CONDITION CurrentPARAMS.NamesOfParameters[param_no] == "SNeg_AS_Prime_Dur" 
//#define PRINT_PARAM_CONDITION CurrentPARAMS.NamesOfParameters[param_no] == "KAM_2_sero1" 
//#define PRINT_PARAM_CONDITION CurrentPARAMS.NamesOfParameters[param_no] == "KAM_2_sero1" || CurrentPARAMS.NamesOfParameters[param_no] == "KAM_2_sero2" || CurrentPARAMS.NamesOfParameters[param_no] == "KAM_2_sero3" || CurrentPARAMS.NamesOfParameters[param_no] == "KAM_2_sero4"
//#define PRINT_PARAM_CONDITION CurrentPARAMS.NamesOfParameters[param_no] == "AM_SNegEff_1" || CurrentPARAMS.NamesOfParameters[param_no] == "PS_SNegEff_1" || CurrentPARAMS.NamesOfParameters[param_no] == "AM_SPosEff_1" || CurrentPARAMS.NamesOfParameters[param_no] == "PS_SPosEff_1"
//#define PRINT_PARAM_CONDITION CurrentPARAMS.NamesOfParameters[param_no] == "PS_SNegEff_1" || CurrentPARAMS.NamesOfParameters[param_no] == "PS_SPosEff_1"
//#define PRINT_PARAM_CONDITION CurrentPARAMS.NamesOfParameters[param_no] == "SNegEff_1" || CurrentPARAMS.NamesOfParameters[param_no] == "SNegEff_2" || CurrentPARAMS.NamesOfParameters[param_no] == "SNegEff_3" || CurrentPARAMS.NamesOfParameters[param_no] == "SNegEff_4" || CurrentPARAMS.NamesOfParameters[param_no] == "SNegEff_AgeKnot_0" || CurrentPARAMS.NamesOfParameters[param_no] == "SNegEff_AgeKnot_1" || CurrentPARAMS.NamesOfParameters[param_no] == "SNegEff_AgeKnot_2" || CurrentPARAMS.NamesOfParameters[param_no] == "SNegEff_AgeKnot_3" 
//#define PRINT_PARAM_CONDITION_2 CurrentPARAMS.NamesOfParameters[sim_update_param] == "SNegEff_1" || CurrentPARAMS.NamesOfParameters[sim_update_param] == "SNegEff_2" || CurrentPARAMS.NamesOfParameters[sim_update_param] == "SNegEff_3" || CurrentPARAMS.NamesOfParameters[sim_update_param] == "SNegEff_4"
//#define PRINT_PARAM_CONDITION CurrentPARAMS.NamesOfParameters[param_no] == "SNegEff_AgeKnot_0" || CurrentPARAMS.NamesOfParameters[param_no] == "SNegEff_AgeKnot_1" || CurrentPARAMS.NamesOfParameters[param_no] == "SNegEff_AgeKnot_2" || CurrentPARAMS.NamesOfParameters[param_no] == "SNegEff_AgeKnot_3" 
#define PRINT_PARAM_CONDITION CurrentPARAMS.NamesOfParameters[param_no] == "SNegEff_AgeGroup_1" || CurrentPARAMS.NamesOfParameters[param_no] == "SNegEff_AgeGroup_2" || CurrentPARAMS.NamesOfParameters[param_no] == "SNegEff_AgeGroup_3"
#define PRINT_PARAM_CONDITION_2 CurrentPARAMS.NamesOfParameters[param_no] == "SNegEff_AgeKnot_0" || CurrentPARAMS.NamesOfParameters[param_no] == "SNegEff_AgeKnot_1" || CurrentPARAMS.NamesOfParameters[param_no] == "SNegEff_AgeKnot_2" || CurrentPARAMS.NamesOfParameters[param_no] == "SNegEff_AgeKnot_3" 
//#define PRINT_PARAM_CONDITION CurrentPARAMS.NamesOfParameters[param_no] == "KAM_2_sero1"
//#define PRINT_PARAM_CONDITION CurrentPARAMS.NamesOfParameters[param_no] == "Knot_0_8"
//#define PRINT_PARAM_CONDITION param_no <= MaxParamNoToPrint


//#define TESTING_Is_and_Find_Param_FUNCTIONS 
//#define TEST_VACCINE_LIKELIHOODS //// will only work if HOUSE.Single_SPos_Eff etc. is set properly. 


//// Loop macros - would save many lines of code but Intellisense hates it so leave. 
#define BaselineSeroStatus_Loop		for (int BaselineSeroStatus = 0	; BaselineSeroStatus	< HOUSE.HowManySeroStatuses		; BaselineSeroStatus++	)		//// Cycle through each serostatus
#define PhaseSeverity_Loop 			for	(int PhaseSeverity = 0		; PhaseSeverity			< HOUSE.HowManyCaseCategories	; PhaseSeverity++		)		//// Cycle through trial phases OR disease severities
#define Age_Loop					for (int age = 0				; age					< HOUSE.HowManyAges				; age++					)		//// Cycle through ages
#define All_Countries_Loop			for (int country = 0			; country				< HOUSE.TotalCountries			; country++				)		//// Cycle through all countries (whether you're fitting all countries or not)
#define Some_Countries_Loop			for (int countryindex = 0		; countryindex			< HOUSE.NoCountriesToFit		; countryindex++		)		//// Cycle through all countries you are fitting
#define Serotype_Loop_KsOrEffs		for (int serotype = 0			; serotype				< HOUSE.N_STypes				; serotype++			)		//// Cycle through serotypes (Serotype_Loop_KsOrEffs
#define Serotype_Loop_Ks			for (int serotype = 0			; serotype				< HOUSE.N_STypes_Ks				; serotype++			)		//// Cycle through serotypes (Ks)
#define Serotype_Loop_VEs			for (int serotype = 0			; serotype				< HOUSE.N_STypes_VEs			; serotype++			)		//// Cycle through serotypes (Effs)
#define PrevInf_Loop				for (int PrevInf = 0			; PrevInf				< HOUSE.Num_K_Params			; PrevInf++				)		//// Cycle through previous infections (K0, K1, K2)
#define Knot_Loop					for (int knot = 0				; knot					< HOUSE.KnotsPerCountry			; knot++				)		//// Cycle through knots


#endif
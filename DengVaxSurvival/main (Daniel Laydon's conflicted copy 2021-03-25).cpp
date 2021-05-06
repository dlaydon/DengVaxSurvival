
#include "HeaderAndSwitches.h"
#include "Initialize.h"
#include "Likelihood.h"
#include "ReadAndProcessData.h"
#include "ParamNumber.h"
#include "ParamUpdate.h"
#include "Augmentation.h"
#include "ConsolePrinting.h"
#include "WriteOutput.h"
#include "MCMC.h"


int main(int argc, char **argv)
{
	// initialize Data, Parameter, Housekeeping, Survival, SeroPrev, and Chains data structures.
	DATA_struct				DATA							;
	Housekeeping_Struct		HOUSE							;	
	Params_Struct			CurrentPARAMS, ProposedPARAMS	;
	Survival_Struct			SURVIVE, SURVIVE_ImSub			;
	Chains_Struct			CHAINS, WBIC_CHAINS				;
	Augment_Struct			AUG								;

#ifdef USE_CLUSTER
	HOUSE.max_threads = set_NCores(); // Open MP. 
#else
	HOUSE.max_threads = 6; 
	omp_set_num_threads(HOUSE.max_threads);
#endif
#ifdef USE_COMMAND_LINE
	string pParamFileName = argv[1]; 
	ReadInParams(HOUSE, CHAINS, WBIC_CHAINS, CurrentPARAMS, pParamFileName);
#else
	CHAINS.No_Iterations	= 20000; 
	CHAINS.BurnIn			=  0;

	HOUSE.ModelVariant	= VAC_SILENT		;
	HOUSE.Do_WBIC_Runs	= false; 
	HOUSE.AdjHaz		= true; 

	HOUSE.DataFilename = "Data\\SimData.txt";
	//HOUSE.FitAllCountries = false; 
	//HOUSE.Fit_c0 = false;
	//////HOUSE.Fit_c1 = false;
	//////HOUSE.Fit_c2 = false;
	//////HOUSE.Fit_c3 = false;
	//////HOUSE.Fit_c4 = false;
	//////HOUSE.Fit_c5 = false;
	//////HOUSE.Fit_c6 = false;
	//HOUSE.Fit_c7 = false;
	//////HOUSE.Fit_c8 = false;
	//////HOUSE.Fit_c9 = false; 

	HOUSE.SSASVE_Additive			= true; 
	HOUSE.SeroSpecificEfficacies	= 1;
	HOUSE.ASVE						= Age_Option::CATEGORICAL;
	HOUSE.AS_Haz					= Age_Option::CATEGORICAL;
	HOUSE.ParamRangeFileName		= "prs1_2";
	HOUSE.EffNegWane				= EffNegWane_Option::NO_WANE;
	HOUSE.PS_Ks_Multiply_AM_Ks		= true; 

	CHAINS.Calc_SCsARsHRPs_Any				= 0;
	CHAINS.AddtoChainEveryHowManyIterations = 1;
#endif
	// Set seeds. 
#if defined(USE_RANDOM_NUMBERS)
	initSeeds(HOUSE.seed1, HOUSE.seed2, HOUSE.max_threads);
#endif
	if (HOUSE.Aug_MH_or_Gibbs == MH_AUG && HOUSE.Weighting_Pass_Sev) 
		for (int i = 0; i < 500; i++) std::cerr << " Metropolis-Hastings Augmentation not coded for Weighting_Pass_Sev" << endl;  

	HOUSE.init();
	Initialize_FollowUp		(DATA, HOUSE);
	InitializeDATApointers	(DATA, HOUSE); 

	//// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// 
	// 1. Import and process data
	//// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// 

	ReadInData	(DATA, HOUSE);
	ProcessData	(DATA, HOUSE);
#ifndef USE_CLUSTER
	if (HOUSE.MakeEveryoneAControl	) for (int pat = 0; pat < NPat; pat++) DATA.Vi_s[pat] = ControlGroup;
	if (HOUSE.MakeEveryoneAVaccinee	) for (int pat = 0; pat < NPat; pat++) DATA.Vi_s[pat] = VaccineGroup;
#endif

	////// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// 
	//// 2. Define Parameters (initial values)
	////// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// 

	Initialize_Params(CurrentPARAMS , CurrentPARAMS, HOUSE, DATA);	 
	Initialize_Params(ProposedPARAMS, CurrentPARAMS, HOUSE, DATA);	 
	// Initial Augmentated Data 
	InitAugDataValues				(DATA, CurrentPARAMS, HOUSE);		// Initializes Augmented Data based on Historical hazards or a previous chain state. 
	AllocateMemoryAndPopulateSets	(DATA, HOUSE);						// Define Set_Array in Data, based on DATA.Ii values. 
	Initialize_IntVacHazardValues	(CurrentPARAMS, ProposedPARAMS, DATA, HOUSE);
	ProcessOldParamChain			(DATA, HOUSE, CurrentPARAMS, ProposedPARAMS);
#ifdef USE_CLUSTER
	std::cerr << "HOUSE.LCPerC "			<< HOUSE.LCPerC			<< endl;
	std::cerr << "Fitting params: "; for (int param_no = 0; param_no < CurrentPARAMS.ParamNosYouWillFit.size(); param_no++) std::cerr << CurrentPARAMS.NamesOfParameters[CurrentPARAMS.ParamNosYouWillFit[param_no]] << ", "; std::cerr << endl;
#endif

	////// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// 
	//// 4. Initialize Structures
	////// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// 

	CHAINS.init			(HOUSE, DATA, CurrentPARAMS); 
	WBIC_CHAINS.init	(HOUSE, DATA, CurrentPARAMS); WBIC_CHAINS.OutputFolder = "Output\\WBIC_";
	Initialize_ParamSeroPrevs(CurrentPARAMS, ProposedPARAMS, HOUSE, CHAINS.SEROPREV);
	SURVIVE.		init		(DATA, HOUSE, CHAINS, WBIC_CHAINS);
	SURVIVE_ImSub.	init		(DATA, HOUSE, CHAINS, WBIC_CHAINS);
	AUG.			init		(HOUSE, DATA); 
	if (!CurrentPARAMS.IsEqualToAnother(ProposedPARAMS, HOUSE, DATA)) 		CurrentPARAMS.WhichPartsUnequal(ProposedPARAMS, HOUSE, DATA);

	LikelihoodFromScratch(DATA, CurrentPARAMS , HOUSE, CurrentPARAMS. LikeParts);
	LikelihoodFromScratch(DATA, ProposedPARAMS, HOUSE, ProposedPARAMS.LikeParts);  

	CurrentPARAMS. LikeFull = l_full(DATA, CurrentPARAMS, HOUSE);	
	ProposedPARAMS.LikeFull = l_full(DATA, ProposedPARAMS, HOUSE);
	CurrentPARAMS. LogPosterior = CurrentPARAMS. LikeFull + CurrentPARAMS. LogPriorFull; 
	ProposedPARAMS.LogPosterior = ProposedPARAMS.LikeFull + ProposedPARAMS.LogPriorFull; 
	std::cerr << "CurrentPARAMS.LikeFull "	<< CurrentPARAMS.LikeFull << " ProposedPARAMS.LikeFull  " << ProposedPARAMS.LikeFull << std::endl << endl;

	//for (int country = 0; country < HOUSE.WhichCountries.size(); country++)
	//{
	//	std::cerr << "\n"; 
	//	print_LikeArray(CurrentPARAMS.LikeParts, HOUSE, HOUSE.WhichCountries[country]);
	//}

	if (HOUSE.Weighting_Pass_Sev)
	{
		CalcWeightedLikes(DATA, CurrentPARAMS, ProposedPARAMS, HOUSE);

		CurrentPARAMS. w_LikeFull	= l_full(DATA, CurrentPARAMS. w_LikeParts, HOUSE);
		ProposedPARAMS.w_LikeFull	= l_full(DATA, ProposedPARAMS.w_LikeParts, HOUSE);
		std::cerr << "CurrentPARAMS.w_LikeFull " << CurrentPARAMS.w_LikeFull << " ProposedPARAMS.w_LikeFull  " << ProposedPARAMS.w_LikeFull << std::endl << endl;
	}

	if (CHAINS.AreWeChecking || CHAINS.CheckIndividualAugmentation)
	{
		Allocate_2D_Array(CHAINS.LikePartsForChecking, HOUSE.TotalCountries, HOUSE.LCPerC); 
		Allocate_3D_Array(CHAINS.VacHazLike_CheckingArray, HOUSE.TotalCountries, HOUSE.HowManySeroStatuses, HOUSE.HowManyCaseCategories);

		LikelihoodFromScratch(DATA, CurrentPARAMS, HOUSE, CHAINS.LikePartsForChecking, CHAINS.VacHazLike_CheckingArray);
	}
	AddToChains(0, CHAINS, DATA, CurrentPARAMS, HOUSE, SURVIVE, SURVIVE_ImSub, AUG);

	//// Checks (whether these functions do anything are controlled by #defines in Macros.h file). 
	Print_Like_Indices				(HOUSE, CurrentPARAMS);
	PrintFittedParamNumbers			(HOUSE, CurrentPARAMS);
	TestAndPrintParamNumberFucntions(HOUSE, CurrentPARAMS, ProposedPARAMS); 

	//////// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// 
	////// 5. MCMC
	//////// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// //// 

	if (HOUSE.Run_MCMC)
	{
		std::cerr << "Start MCMC " << endl;

		//////// //// //// //// //// //// //// //////// //// //// //// //// //// //// 
		//////// First chain

		if (!HOUSE.WBIC_RunOnly) 
			MCMC_All(1, CHAINS, DATA, CurrentPARAMS, ProposedPARAMS, AUG, HOUSE, SURVIVE, SURVIVE_ImSub); 
		if (!HOUSE.WBIC_RunOnly) //// Don't write anything from regular chains if you're only running WBICs. 
			WriteEverythingForParticularChain(DATA, SURVIVE, SURVIVE_ImSub, CurrentPARAMS, HOUSE, CHAINS, AUG);

		//////// //// //// //// //// //// //// //////// //// //// //// //// //// //// 
		//////// Second (WBIC) chain

		if (HOUSE.Do_WBIC_Runs)
		{
			//// reset values in chains - don't output seroprevs or survival curves / attack rates. 
			if (WBIC_CHAINS.Calc_SCsARsHRPs_Any == true) //// i.e. only do this if default unchanged. 
			{
				WBIC_CHAINS.Calc_SCsARsHRPs_MeanAndCrIs = false;
				WBIC_CHAINS.Calc_SCsARsHRPs_ModalMaxLike = true;
			}
			WBIC_CHAINS.AreWeFittingParameters	= true;
			WBIC_CHAINS.SimulatedAnnealing		= true; //// having restarted chains - don't vary temperature (keep it fixed at 1 / log(NPat))
			WBIC_CHAINS.CoolDuringBurnIn		= true;
			WBIC_CHAINS.CoolAfterBurnIn			= false;

			std::cerr << "Start MCMC: WBIC run " << endl;	fflush(stderr);

			MCMC_All(0, WBIC_CHAINS, DATA, CurrentPARAMS, ProposedPARAMS, AUG, HOUSE, SURVIVE, SURVIVE_ImSub);
		}
		if (HOUSE.Do_WBIC_Runs)
			WriteEverythingForParticularChain(DATA, SURVIVE, SURVIVE_ImSub, CurrentPARAMS, HOUSE, WBIC_CHAINS	, AUG);
	}
	else
	{
		ImportOldParamChains(HOUSE, CHAINS, WBIC_CHAINS);
		std::vector<DType>	DummyParamVec(HOUSE.No_Parameters);
		//// compute Total Log Likelihood for param chain values. 
		std::cerr << "compute Total Log Likelihood for old param chain values " << endl;

		for (int PostSample = 0; PostSample < CHAINS.NumPosteriorSamples; PostSample++)
		{
			std::cerr << "i" << PostSample << " ";
			if (PostSample % 100 == 0) std::cerr << endl;
			//// extract param vec for this posterior Sample. 
			for (int param_no = 0; param_no < HOUSE.No_Parameters; param_no++)
				CurrentPARAMS.ParamVec	[param_no] = CHAINS.ParamChain[param_no][PostSample];

			AmendParams(CurrentPARAMS.ParamVec, CurrentPARAMS, DATA, HOUSE); 
			CurrentPARAMS.Total_logLike = LL_Total(DATA, CurrentPARAMS, HOUSE, AUG); 
			CHAINS.Total_Likelihood_Chain[PostSample] = CurrentPARAMS.Total_logLike; 
		}
		WriteSingleChain(CHAINS.OutputFolder + "TotalLogLikeChain" + HOUSE.OutputString + ".txt", CHAINS.Total_Likelihood_Chain, CHAINS.NumPosteriorSamples, "TotalLogLike");

		//// compute Total Log Likelihood for WBIC param chain values. 
		std::cerr << "compute Total Log Likelihood for old WBIC param chain values " << endl;
		for (int PostSample = 0; PostSample < WBIC_CHAINS.NumPosteriorSamples; PostSample++)
		{
			std::cerr << "i" << PostSample << " ";
			if (PostSample % 100 == 0) std::cerr << endl;
			//// extract param vec for this posterior Sample. 
			for (int param_no = 0; param_no < HOUSE.No_Parameters; param_no++)
				CurrentPARAMS.ParamVec	[param_no] = WBIC_CHAINS.ParamChain[param_no][PostSample];

			AmendParams(CurrentPARAMS.ParamVec, CurrentPARAMS, DATA, HOUSE); 
			CurrentPARAMS.Total_logLike = LL_Total(DATA, CurrentPARAMS, HOUSE, AUG); 
			WBIC_CHAINS.Total_Likelihood_Chain[PostSample] = CurrentPARAMS.Total_logLike;
		}
		WriteSingleChain(WBIC_CHAINS.OutputFolder + "TotalLogLikeChain" + HOUSE.OutputString + ".txt", WBIC_CHAINS.Total_Likelihood_Chain, WBIC_CHAINS.NumPosteriorSamples, "TotalLogLike");
	}

	//// calculate average acceptance probabilities for regular and WBIC chains
	CHAINS.		ProportionUpdatesAccepted	/= (CHAINS.		HowFarDidYouGet * CurrentPARAMS.ParamNosYouWillFit.size());
	WBIC_CHAINS.ProportionUpdatesAccepted	/= (WBIC_CHAINS.HowFarDidYouGet * CurrentPARAMS.ParamNosYouWillFit.size());
	CHAINS.		AverageAcceptanceProb		/= (CHAINS.		HowFarDidYouGet * CurrentPARAMS.ParamNosYouWillFit.size());
	WBIC_CHAINS.AverageAcceptanceProb		/= (WBIC_CHAINS.HowFarDidYouGet * CurrentPARAMS.ParamNosYouWillFit.size());

	std::cerr << " Regular Chains ProportionUpdatesAccepted " << CHAINS.ProportionUpdatesAccepted << std::endl;
	std::cerr << " WBIC Chains ProportionUpdatesAccepted " << WBIC_CHAINS.ProportionUpdatesAccepted << std::endl;

	std::cerr << " Regular Chains AverageAcceptanceProb " << CHAINS.AverageAcceptanceProb << std::endl;
	std::cerr << " WBIC Chains AverageAcceptanceProb " << WBIC_CHAINS.AverageAcceptanceProb << std::endl;

	DeleteData(DATA, HOUSE); 

	std::cerr << "main finished: return 0 next line" << endl; std::fflush(stderr);
	return 0;
}

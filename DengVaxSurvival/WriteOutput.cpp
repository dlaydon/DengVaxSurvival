
#include "HeaderAndSwitches.h"
#include "Likelihood.h"
#include "SurvivalCurves.h"
#include "ParamUpdate.h"
#include "Augmentation.h"
#include "Probability.h"
#include "MCMC.h"


void WriteOutput						(string FILENAME, DType **Object, int NRows, int NCols)
{
	ofstream File;
	File.open(FILENAME);
	for (int row = 0; row < NRows; row++)
	{
		for (int col = 0; col < NCols - 1; col++) File << Object[row][col] << "\t";
		File << Object[row][NCols - 1];
		File << endl;
	}
	File.close();
}
void WriteParameterChainOutput			(string FILENAME, std::vector<string> NamesOfParameters, const Housekeeping_Struct &HOUSE, Chains_Struct &CHAINS)
{
	ofstream outputfile;
	outputfile.open(FILENAME);
	
	if ((HOUSE.SeroSpecificEfficacies || HOUSE.SeroSpecific_K_values) && (!HOUSE.BaselinePartition))
	{
		outputfile << "logLikelihood" << "\t"; 	for (int param_no = 0; param_no < HOUSE.No_Parameters; param_no++)	outputfile << NamesOfParameters[param_no] << "\t";	
		for (int country = 0; country < HOUSE.TotalCountries; country++)
			for (int serotype = 0; serotype < HOUSE.N_STypes; serotype++)
			{
				if ((country == (HOUSE.TotalCountries - 1)) && (serotype == (HOUSE.N_STypes - 1)))	outputfile << "rho_" + std::to_string(country) + "_" + std::to_string(serotype + 1);			
				else 																				outputfile << "rho_" + std::to_string(country) + "_" + std::to_string(serotype + 1) << "\t";
			}
		outputfile << endl;
	}
	else
	{
		outputfile << "logLikelihood" << "\t"; 	for (int param_no = 0; param_no < HOUSE.No_Parameters - 1; param_no++)	outputfile << NamesOfParameters[param_no] << "\t";	outputfile << NamesOfParameters[HOUSE.No_Parameters - 1];
		outputfile << endl;
	}

	for (int index = 0; index < CHAINS.NumPosteriorSamples; index++)
	{
		if ((HOUSE.SeroSpecificEfficacies || HOUSE.SeroSpecific_K_values) && (!HOUSE.BaselinePartition))
		{
			outputfile << CHAINS.logLikeChain[index] << "\t"; 	for (int param_no = 0; param_no < HOUSE.No_Parameters; param_no++)	outputfile << CHAINS.ParamChain[param_no][index] << "\t";
			for (int country = 0; country < HOUSE.TotalCountries; country++)
				for (int serotype = 0; serotype < HOUSE.N_STypes; serotype++)
				{
					if	((country == (HOUSE.TotalCountries - 1)) && (serotype == (HOUSE.N_STypes - 1))) outputfile << CHAINS.rhos[country][serotype][index];			
					else 																				outputfile << CHAINS.rhos[country][serotype][index] << "\t";
				}
		}
		else
		{
			outputfile << CHAINS.logLikeChain[index] << "\t"; 	for (int param_no = 0; param_no < HOUSE.No_Parameters - 1; param_no++)	outputfile << CHAINS.ParamChain[param_no][index] << "\t";
			outputfile << CHAINS.ParamChain[HOUSE.No_Parameters - 1][index];									
		}
		outputfile << endl;
	}
	outputfile.close();
}
void WriteSingleChain					(string FILENAME, DType *Chain, int NumPosteriorSamples, string Colname)
{
	ofstream outputfile;
	outputfile.open(FILENAME);
	outputfile << Colname << std::endl;
	for (int index = 0; index < NumPosteriorSamples; index++) 
		outputfile << Chain[index] << std::endl;
	outputfile.close();
}
void Write_LL_minus_Aug_Chain			(string FILENAME, Chains_Struct &CHAINS)
{
	ofstream outputfile;
	outputfile.open(FILENAME);
	outputfile << "logLikelihood_minus_aug" << endl;
	for (int index = 0; index < CHAINS.NumPosteriorSamples; index++)	outputfile << CHAINS.logLike_minus_Aug_Chain[index] << endl;
	outputfile.close();
}
void Write_SPrev_LL_Chain(string FILENAME, DType *** Which_LL_SPrevChain, const Chains_Struct &CHAINS, const Housekeeping_Struct &HOUSE)
{
	//// Create Colnames. 
	std::vector<string> Colnames;
	string BS_String = ""; 
	for (int country = 0; country < HOUSE.TotalCountries + 3; country++)
		for (int BS = 0; BS < HOUSE.HowManySeroStatuses; BS++)
		{
			if (BS == SeroNeg) BS_String = "SeroNeg"; else if (BS == SeroPos) BS_String = "SeroPos"; else std::cerr << "Write_SPrev_LL_Chain error: BS not recognized " << endl; 
			Colnames.push_back("c" + std::to_string(country) + "_" + BS_String);
		}
	ofstream outputfile;
	outputfile.open(FILENAME);
	//// Write Colnames; 
	for (int col = 0; col < Colnames.size() - 1; col++)	
		outputfile << Colnames[col] << "\t"; 	
	outputfile << Colnames[Colnames.size() - 1] << std::endl;

	//// Write values
	for (int index = 0; index < CHAINS.NumPosteriorSamples; index++)		
		for (int country = 0; country < HOUSE.TotalCountries + 3; country++)
			for (int BS = 0; BS < HOUSE.HowManySeroStatuses; BS++)
			{
				outputfile << Which_LL_SPrevChain[index][country][BS];
				if (country == HOUSE.TotalCountries + 2 & BS == HOUSE.HowManySeroStatuses - 1) outputfile << "\n";  else outputfile << "\t";
			}
	outputfile.close();
}
void WriteAcceptanceArray				(string FILENAME, DType * AcceptanceArray, DType * AcceptanceArray_Aug, int max_threads, int NumPosteriorSamples, int No_Parameters, int NoAugmented, int No_Iterations, int AugmentEveryHowManyIterations, std::vector<string> NamesOfParameters)
{
	DType Total_Aug_Acceptances = 0;
	for (int thread_no = 0; thread_no < max_threads; thread_no++) Total_Aug_Acceptances += AcceptanceArray_Aug[thread_no];
	DType Prop_Aug_Acceptances = Total_Aug_Acceptances / (NumPosteriorSamples * NoAugmented / AugmentEveryHowManyIterations);
	
	ofstream outputfile_AcceptanceArray;
	outputfile_AcceptanceArray.open(FILENAME);
	for (int param_no = 0; param_no < No_Parameters; param_no++)	outputfile_AcceptanceArray << NamesOfParameters[param_no]					<< "\t";	outputfile_AcceptanceArray << "AugmentedData";
	outputfile_AcceptanceArray << endl;
	for (int param_no = 0; param_no < No_Parameters; param_no++)	outputfile_AcceptanceArray << (AcceptanceArray[param_no] / No_Iterations)	<< "\t"	;	outputfile_AcceptanceArray << Prop_Aug_Acceptances;
	outputfile_AcceptanceArray << endl;
	outputfile_AcceptanceArray.close();
}

void WriteSurvivalTable					(string FILENAME, DType **Object, int NRows, int NCols, std::vector<string> RowNames, std::vector<string> ColNames)
{
	ofstream File;
	File.open(FILENAME);

	//// colnames
	File << "Age_Country_Arm_Hazard_Group" << "\t";	
	for (int timepoint = 0; timepoint < NCols - 1; timepoint++) File << ColNames[timepoint] << "\t";
	File << ColNames[NCols - 1] << endl;

	// values
	for (int row = 0; row < NRows; row++)
	{
		// rownames
		File << RowNames[row] << "\t";
		//values 
		for (int timepoint = 0; timepoint < NCols - 1; timepoint++) File << Object[row][timepoint] << "\t";
		File << Object[row][NCols - 1] << endl;
	}
	File.close();
}
void WriteSurivalOutput					(string OutputFolder, Survival_Struct &SURVIVE, Chains_Struct &CHAINS, const Housekeeping_Struct &HOUSE, string ImSubOrNotString)
{
	//// function is a wrapper of WriteSurvivalTable function. Writes lots of Survival tables, and attack rates, and hazard ratios. 
	if (CHAINS.Calc_SCsARsHRPs_MeanAndCrIs)
	{
		/////// Calculate and Output Mean/Median and Credible Intervals for survival curves (and hazard ratios). 
		CalculateSurvivalCurveOutput(SURVIVE);

		//// **** //// **** //// ****  write Survival tables.
		//// Whole trial survival tables
		if (!CHAINS.HazRatiosOnly)
			for (int DiseaseSeverity = 0; DiseaseSeverity < SURVIVE.HowManyDiseaseSeverities; DiseaseSeverity++)
				for (int Statistic = 0; Statistic < 3; Statistic++)
					WriteSurvivalTable(OutputFolder + SURVIVE.StatisticNames[Statistic] + "SurvivalTable" + ImSubOrNotString + HOUSE.OutputString + SURVIVE.DiseaseNames[DiseaseSeverity] + ".txt",
						SURVIVE.SC_WT.FinalPosteriorSurvivalCurves[DiseaseSeverity][Statistic],
						SURVIVE.NoSubjectCategories, SURVIVE.NoDaysOfFollowUp + 1, SURVIVE.SurvivalTableRowNames, SURVIVE.SurvivalTableColNames);

		//// passive phase survival tables. 
		if (SURVIVE.CalculateSeparatePassivePhaseCurves & !CHAINS.HazRatiosOnly)
			for (int DiseaseSeverity = 0; DiseaseSeverity < SURVIVE.HowManyDiseaseSeverities; DiseaseSeverity++)
				for (int Statistic = 0; Statistic < 3; Statistic++)
					WriteSurvivalTable(OutputFolder + SURVIVE.StatisticNames[Statistic] + "PassiveSurvivalTable" + ImSubOrNotString + HOUSE.OutputString + SURVIVE.DiseaseNames[DiseaseSeverity] + ".txt",
						SURVIVE.SC_PP.FinalPosteriorSurvivalCurves[DiseaseSeverity][Statistic],
						SURVIVE.NoSubjectCategories, SURVIVE.NoPassiveDaysFollowUp + 1, SURVIVE.SurvivalTableRowNames, SURVIVE.PassivePhaseSurvivalTableColNames);

		//// attack rates (individual posterior samples)
		if (SURVIVE.NoSurvivePostSamples > 0 & !CHAINS.HazRatiosOnly)
			for (int TrialPhase = 0; TrialPhase < SURVIVE.HowManyTrialPhases; TrialPhase++)
				for (int DiseaseSeverity = 0; DiseaseSeverity < SURVIVE.HowManyDiseaseSeverities; DiseaseSeverity++)
					WriteSurvivalTable(OutputFolder + "AttackRates" + ImSubOrNotString + HOUSE.OutputString + SURVIVE.TrialPhaseNames[TrialPhase] + SURVIVE.DiseaseNames[DiseaseSeverity] + ".txt",	//// filename
										SURVIVE.MetaAttackRates[TrialPhase][DiseaseSeverity],																						//// point to correct table. 
										SURVIVE.NoSubjectCategories, SURVIVE.NoSurvivePostSamples, SURVIVE.SurvivalTableRowNames, SURVIVE.AttackRateColNames);

		//// hazard ratios. 
		if (SURVIVE.NoSurvivePostSamples > 0)
			for (int DiseaseSeverity = 0; DiseaseSeverity < SURVIVE.HowManyDiseaseSeverities; DiseaseSeverity++)
				for (int Statistic = 0; Statistic < 3; Statistic++)
					for (int HR_serotype = 0; HR_serotype < SURVIVE.HRs_NumSTypes; HR_serotype++)
						WriteSurvivalTable(OutputFolder + SURVIVE.StatisticNames[Statistic] + "HazRatios" + SURVIVE.HRs_SeroNames[HR_serotype] + ImSubOrNotString + HOUSE.OutputString + SURVIVE.DiseaseNames[DiseaseSeverity] + ".txt",
							SURVIVE.HRs_BySerotype[HR_serotype].FinalPosteriorSurvivalCurves[DiseaseSeverity][Statistic],
							SURVIVE.HRs_NumStrata, SURVIVE.HRs_NumDaysPostDoseToCalculate, SURVIVE.HRsRowNames, SURVIVE.HRsColNames);
	}
}
void WriteSeroPrevOutput				(string FILENAME, DType **AugOnlySeroPrevs, DType **AllSeroPrevs, DType ***AgeSpecificSeroPrevs, int HowManyAges, int NoCountries) 
{
	ofstream File;
	File.open(FILENAME);
	// colnames
	File << "\t";
	for (int age = 0; age < HowManyAges; age++)	File << "age " << std::to_string(age) << "\t";
	File << endl;
	// values
	for (int country = 0; country < (NoCountries + 3); country++)
	{
		File << "Country_" + std::to_string(country) + "_NonAugmented" << "\t"; // rowname
		for (int age = 0; age < HowManyAges; age++)	File << AgeSpecificSeroPrevs[country][NonAugIndex][age] << "\t";// values 
		File << endl;
		File << "Country_" + std::to_string(country) + "_Augmented" << "\t"; // rowname
		for (int age = 0; age < HowManyAges; age++)	File << AugOnlySeroPrevs[country][age] << "\t";// values 
		File << endl;
		File << "Country_" + std::to_string(country) + "_All" << "\t"; // rowname
		for (int age = 0; age < HowManyAges; age++)	File << AllSeroPrevs[country][age] << "\t";// values 
		File << endl;
	}
	File.close();
}
void WriteAllOutput						(const DATA_struct &DATA, Survival_Struct &SURVIVE, Survival_Struct &SURVIVE_ImSub, 
	const Params_Struct &CurrentPARAMS, const Housekeeping_Struct &HOUSE, Chains_Struct &CHAINS) /// wrapper functions of the above. Just to make main shorter. 
{
	if (CHAINS.NumPosteriorSamples > 0)
	{
		string OutputFolder = CHAINS.OutputFolder; // "Output\\" + WBIC_Or_Not_String; //// 

		WriteParameterChainOutput	(OutputFolder + "ParameterChainOutput"		+ HOUSE.OutputString + ".txt", CurrentPARAMS.NamesOfParameters, HOUSE, CHAINS);
		WriteSingleChain			(OutputFolder + "LogPosteriorChain"			+ HOUSE.OutputString + ".txt", CHAINS.log_PosteriorChain	, CHAINS.NumPosteriorSamples, "logPosterior");
		WriteSingleChain			(OutputFolder + "TotalLogLikeChain"			+ HOUSE.OutputString + ".txt", CHAINS.Total_Likelihood_Chain, CHAINS.NumPosteriorSamples, "TotalLogLike");
		WriteAcceptanceArray		(OutputFolder + "AcceptanceArray"			+ HOUSE.OutputString + ".txt", CHAINS.AcceptArray	, CHAINS.AcceptArray_Aug, HOUSE.max_threads, CHAINS.NumPosteriorSamples, HOUSE.No_Parameters, DATA.NoAugmented, CHAINS.No_Iterations, CHAINS.AugmentEveryHowManyIterations, CurrentPARAMS.NamesOfParameters);
		
		if (HOUSE.OutputSeroPrev_LLChains)
		{
			Write_SPrev_LL_Chain(CHAINS.OutputFolder + "SPrev_LL_Chain"			+ HOUSE.OutputString + ".txt", CHAINS.LL_SPrev_Chains		, CHAINS, HOUSE);
			Write_SPrev_LL_Chain(CHAINS.OutputFolder + "SPrev_LL_Chain_ImSub"	+ HOUSE.OutputString + ".txt", CHAINS.LL_SPrev_Chains_ImSub	, CHAINS, HOUSE);
		}
		if (CHAINS.Output_LL_minus_Aug)		Write_LL_minus_Aug_Chain(OutputFolder + "LL_No_Aug" + HOUSE.OutputString + ".txt", CHAINS);

		if (CHAINS.AreWeCalculatingSeroPrevs)
		{
			/////// Calculate and Output Mean/Median and Credible Intervals for age-specific sero-prevalences at baseline. 
			CHAINS.SEROPREV.Calc_Final_SeroPrev_Output(HOUSE, CHAINS.NumPosteriorSamples, CHAINS.NumElementsOutside_CrI_Tails);

			WriteSeroPrevOutput(OutputFolder + "Mean_SeroPrevOutput"	+ HOUSE.OutputString + ".txt", CHAINS.SEROPREV.Aug_Patients_Mean_AgeSpecificSeroPrevs		, CHAINS.SEROPREV.All_Patients_Mean_AgeSpecificSeroPrevs	, CHAINS.SEROPREV.AgeSpecificSeroPrevs, HOUSE.HowManyAges, HOUSE.TotalCountries);
			WriteSeroPrevOutput(OutputFolder + "Median_SeroPrevOutput"	+ HOUSE.OutputString + ".txt", CHAINS.SEROPREV.Aug_Patients_Median_AgeSpecificSeroPrevs		, CHAINS.SEROPREV.All_Patients_Median_AgeSpecificSeroPrevs	, CHAINS.SEROPREV.AgeSpecificSeroPrevs, HOUSE.HowManyAges, HOUSE.TotalCountries);
			WriteSeroPrevOutput(OutputFolder + "LowerCI_SeroPrevOutput" + HOUSE.OutputString + ".txt", CHAINS.SEROPREV.Aug_Patients_LowerCI_AgeSpecificSeroPrevs	, CHAINS.SEROPREV.All_Patients_LowerCI_AgeSpecificSeroPrevs	, CHAINS.SEROPREV.AgeSpecificSeroPrevs, HOUSE.HowManyAges, HOUSE.TotalCountries);
			WriteSeroPrevOutput(OutputFolder + "UpperCI_SeroPrevOutput" + HOUSE.OutputString + ".txt", CHAINS.SEROPREV.Aug_Patients_UpperCI_AgeSpecificSeroPrevs	, CHAINS.SEROPREV.All_Patients_UpperCI_AgeSpecificSeroPrevs	, CHAINS.SEROPREV.AgeSpecificSeroPrevs, HOUSE.HowManyAges, HOUSE.TotalCountries);
		}

		if (CHAINS.Calc_SCsARsHRPs_MeanAndCrIs)
		{
			//// write Survival tables.
			WriteSurivalOutput(OutputFolder, SURVIVE		, CHAINS, HOUSE, ""			);  //// entire trial
			WriteSurivalOutput(OutputFolder, SURVIVE_ImSub	, CHAINS, HOUSE, "_ImSub"	);  //// Immune subset
		}
	}
	std::cerr << "WriteAllOutput function finished: " << endl; std::fflush(stderr);
}

void WriteAugDataOutput					(int * SerostatusArray, string FileName)
{
	std::ofstream Ii_sFile_output;
	Ii_sFile_output.open(FileName);
	for (int row = 0; row < NPat; row++)	Ii_sFile_output << SerostatusArray[row] << endl;
	Ii_sFile_output.close();
}
void WriteAugDataOutput					(DType * SerostatusArray, string FileName)
{
	std::ofstream Ii_sFile_output;
	Ii_sFile_output.open(FileName);
	for (int row = 0; row < NPat; row++)	Ii_sFile_output << SerostatusArray[row] << endl;
	Ii_sFile_output.close();
}
void Write_DIC_Output					(DATA_struct &DATA, Params_Struct &CurrentPARAMS, const Housekeeping_Struct &HOUSE, Chains_Struct &CHAINS)
{
	//// DIC given by either: 
	/// i)  DIC = logLike (MeanPosterior) - 2 Var (LogLikeChain   ); OR....
	/// ii) DIC = logLike (MeanPosterior) - 2 Var (PosteriorChain );  //// depending on whether you consider knots or logged knots as fitted parameters. 

	DType LogLikeMeanPosterior; 
	//// amend parameters with MeanPostParamVec
	////// Amend Augmented data with MaxLike_Ii_s (CHAINS.MeanPost_Ii_s are probabilities so  won't work). 
	UpdateOrResetAugData(DATA.Ii_s, CHAINS.MaxLike_Ii_s, DATA, HOUSE);
	//// amend parameters with MeanPostParamVec
	AmendParams(CHAINS.MeanPostParamVec, CurrentPARAMS, DATA, HOUSE);
	LikelihoodFromScratch(DATA, CurrentPARAMS, HOUSE, CurrentPARAMS.LikeParts); 
	CurrentPARAMS.LikeFull = l_full(DATA, CurrentPARAMS.LikeParts, HOUSE); 
	LogLikeMeanPosterior = CurrentPARAMS.LikeFull; 

	DType Var_LogPostChain_v1 = 0;
	DType Var_LogPostChain_v2 = 0;
	DType MeanLogLike_v1 = 0; 
	/// calculate means
	for (int sample = 0; sample < CHAINS.NumPosteriorSamples; sample++) MeanLogLike_v1 += CHAINS.logLikeChain[sample]; 
	MeanLogLike_v1 /= CHAINS.NumPosteriorSamples; 
	DType MeanLogLike_v2 = 0;
	for (int sample = 0; sample < CHAINS.NumPosteriorSamples; sample++) MeanLogLike_v2 += CHAINS.log_PosteriorChain[sample];
	MeanLogLike_v2 /= CHAINS.NumPosteriorSamples;
	/// calculate variances
	for (int sample = 0; sample < CHAINS.NumPosteriorSamples; sample++) Var_LogPostChain_v1 += (CHAINS.logLikeChain			[sample] - MeanLogLike_v1) * (CHAINS.logLikeChain		[sample] - MeanLogLike_v1); /// i.e. squared
	for (int sample = 0; sample < CHAINS.NumPosteriorSamples; sample++) Var_LogPostChain_v2 += (CHAINS.log_PosteriorChain	[sample] - MeanLogLike_v2) * (CHAINS.log_PosteriorChain	[sample] - MeanLogLike_v2); /// i.e. squared
	Var_LogPostChain_v1 /= CHAINS.NumPosteriorSamples; 
	Var_LogPostChain_v2	/= CHAINS.NumPosteriorSamples; 
	DType DIC_defn1 = DIC(CHAINS.NumPosteriorSamples, LogLikeMeanPosterior	, CHAINS.logLikeChain);
	DType DIC_defn2 = DIC(CHAINS.NumPosteriorSamples, LogLikeMeanPosterior	, CHAINS.log_PosteriorChain);
	DType DIC_defn3 = DIC(CHAINS.NumPosteriorSamples, CHAINS.MaxLikeSoFar	, CHAINS.logLikeChain);
	DType DIC_defn4 = DIC(CHAINS.NumPosteriorSamples, CHAINS.MaxLikeSoFar	, CHAINS.log_PosteriorChain);
	DType DIC_defn5 = DIC(CHAINS.NumPosteriorSamples, CHAINS.ModalPostSoFar , CHAINS.logLikeChain);
	DType DIC_defn6 = DIC(CHAINS.NumPosteriorSamples, CHAINS.ModalPostSoFar	, CHAINS.log_PosteriorChain);

	std::cerr << "DICs: defn 1 " << DIC_defn1 << "  defn 2 " << DIC_defn2 << std::endl;
	std::ofstream DIC_output;
	DIC_output.open(CHAINS.OutputFolder + "DIC" + HOUSE.OutputString + ".txt");
	DIC_output << "Quantity"							<< "\t" << "Value"					<< std::endl; //// colnames
	DIC_output << "DIC_defn1_LLChain"					<< "\t" << DIC_defn1				<< std::endl;
	DIC_output << "DIC_defn2_LogPostChain"				<< "\t" << DIC_defn2				<< std::endl;
	DIC_output << "DIC_defn3_MaxLike_LLChain"			<< "\t" << DIC_defn3				<< std::endl;
	DIC_output << "DIC_defn4_MaxLike_LogPostChain"		<< "\t" << DIC_defn4				<< std::endl;
	DIC_output << "DIC_defn5_ModalPost_LLChain"			<< "\t" << DIC_defn5				<< std::endl;
	DIC_output << "DIC_defn6_ModalPost_LogPostChain"	<< "\t" << DIC_defn6				<< std::endl;
	DIC_output << "MeanLogLike_v1"						<< "\t" << MeanLogLike_v1			<< std::endl;
	DIC_output << "MeanLogLike_v2"						<< "\t" << MeanLogLike_v2			<< std::endl;
	DIC_output << "Var_LogPostChain_v1"					<< "\t" << Var_LogPostChain_v1		<< std::endl;
	DIC_output << "Var_LogPostChain_v2"					<< "\t" << Var_LogPostChain_v2		<< std::endl;
	DIC_output << "LogLikeMeanPosterior"				<< "\t" << LogLikeMeanPosterior		<< std::endl;
	DIC_output << "MaxLikeValue"						<< "\t" << CHAINS.MaxLikeSoFar		<< std::endl;
	DIC_output << "ModalPostValue"						<< "\t" << CHAINS.ModalPostSoFar	<< std::endl;
	DIC_output.close();
}

void WriteEverythingForParticularChain (DATA_struct &DATA, Survival_Struct &SURVIVE, Survival_Struct &SURVIVE_ImSub,
	Params_Struct &CurrentPARAMS, const Housekeeping_Struct &HOUSE, Chains_Struct &CHAINS, Augment_Struct &AUG)
{
	/*
		Format of this function is (obviously) to output. However also to calculate various summary quantities - DIC, Modal posterior survival curves and attack rates etc. 
		For the latter, need to alter Parameters and augmented Data, i.e. CurrentPARAMS, and DATA.Ii_s (including strata "Sets"). 
		Must alter augmented data first, as some quantities in CurrentPARAMS depend on assignment of patients' serostatus (e.g. vaccine hazards, which would affect likelihood for DIC).
		Must record state of the Param and augmented data chains BEFORE they're amended. This is done at end of iterations loop in MCMC_All function in ParamUpdate.cpp
	*/
	//// final state of Augmented Data (must be done before you change DATA.Ii_s)
	WriteAugDataOutput(CHAINS.Final_Iis, CHAINS.OutputFolder + "FinalState_AugData" + HOUSE.OutputString + ".txt"); //// CHAINS.Final Iis calculated at end of iterations loop in MCMC_All function in ParamUpdate.cpp

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// // //			OUTPUT Parameters likelihood chains AcceptanceArray and AcceptanceArray_aug 
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	WriteAllOutput(DATA, SURVIVE, SURVIVE_ImSub, CurrentPARAMS, HOUSE, CHAINS);

	//// Get Mean Posterior vector
	CHAINS.MeanPostParamVec = Calc_RowMeans(CHAINS.ParamChain, HOUSE.No_Parameters, CHAINS.NumPosteriorSamples); 


	///// **** ///// ///// **** ///// ///// **** ///// ///// **** ///// ///// **** ///// ///// **** ///// 
	///// **** ///// ///// **** /////		OUTPUT AUGMENTED DATA (and related quantities).

	//// Output AugmentedData (and process MeanPost_Ii_s). 
	//// Divide augmented data mean posterior, then write. 
	for (int patient = 0; patient < NPat; patient++) CHAINS.MeanPost_Ii_s[patient] /= CHAINS.NumPosteriorSamples;
	WriteAugDataOutput(CHAINS.MeanPost_Ii_s, CHAINS.OutputFolder + "MeanPost_AugData" + HOUSE.OutputString + ".txt");
	//// max like / modal post state of Augmented Data
	WriteAugDataOutput(CHAINS.MaxLike_Ii_s, CHAINS.OutputFolder + "ModalPost_AugData" + HOUSE.OutputString + ".txt");

	///// **** ///// ///// **** ///// ///// **** ///// ///// **** ///// ///// **** ///// ///// **** ///// 
	///// **** ///// ///// **** /////		CALCULATE AND OUTPUT DIC

	Write_DIC_Output(DATA, CurrentPARAMS, HOUSE, CHAINS); 


	//// Output Max Like and modal parameters (guard against their being empty)
	if (CHAINS.MaxLike_ParamVec.size() > 0 & CHAINS.ModalPost_ParamVec.size() > 0)
	{
		std::ofstream MaxLike_ModalPost_ParamChain;
		MaxLike_ModalPost_ParamChain.open(CHAINS.OutputFolder + "ParamsMaxLikeModalPost" + HOUSE.OutputString + ".txt");
		/// Colnames
		MaxLike_ModalPost_ParamChain << "ParamName " << "\t" << "MaxLike " << "\t" << "ModalPost " << "\t" << "MeanPost" << "\t" << std::endl;
		for (int param_no_dummy = 0; param_no_dummy < HOUSE.No_Parameters; param_no_dummy++)
			MaxLike_ModalPost_ParamChain << CurrentPARAMS.NamesOfParameters[param_no_dummy] << "\t" << CHAINS.MaxLike_ParamVec[param_no_dummy] << "\t" << CHAINS.ModalPost_ParamVec[param_no_dummy] << "\t" << CHAINS.MeanPostParamVec[param_no_dummy] << std::endl;
		MaxLike_ModalPost_ParamChain.close();
	}

	///// **** ///// ///// **** ///// ///// **** ///// ///// **** ///// ///// **** ///// ///// **** ///// 
	///// **** ///// ///// **** /////		OUTPUT MEAN and MODAL POSTERIOR SURVIVAL CURVES AND ATTACK RATES (must be done after other output written. Much of SURVIVE structures are cumulative - applying GenerateSurvivalCurves function, even just to output a single sample, adds to cumulative means etc.)

	if (CHAINS.Calc_SCsARsHRPs_ModalMaxLike & !CHAINS.HazRatiosOnly)
	{
		//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** 
		//// **** modal posteroir survival curves and attack rates. (this is the modal posterior only if you consider non-logged knots as the parameters)
		//// In the three calls to GenerateSurvivalCurves, you output Survival Curve tables for Modal_SurvivalTable and MaxLike_SurvivalTable only
		//// However in each call you add to Attack rates, which you then output at the end. 
		//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** 

		////// Amend Augmented data with ModalPost_Ii_s
		UpdateOrResetAugData(DATA.Ii_s, CHAINS.ModalPost_Ii_s, DATA, HOUSE);
		//// amend parameters with ModalPost_ParamVec
		AmendParams(CHAINS.ModalPost_ParamVec, CurrentPARAMS, DATA, HOUSE);
		//// generate modal posterior curves and attack rates and output s curves
		GenerateSurvivalCurves(DATA, CurrentPARAMS, HOUSE, SURVIVE, MODAL_POST	, CHAINS, AUG, 
			/*AllPatients = */true, /*OutputIndividualSurvivalTables = */ true, /*OutputIndividual_Passive_SurvivalTables = */ false,
			/*DType **** AttackRateContainer = */ SURVIVE.MeanModeAttackRates, /*ParamSet_SCurve_RootFilename = */"Modal_SurvivalTable"	, /*ParamSet_PassiveSCurve_RootFilename = */ "Modal_PassiveSurvivalTable"		, 
			/*Use_Default_AR_Post_Sample_No = */ false, /*CleanFirstDay_IndividualCurves =*/ true);

		//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** 
		//// **** MaxLike survival curves and attack rates.  (this is the modal posterior if you consider logged knots as the parameters)
		//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** 

		////// Amend Augmented data with MaxLike_Ii_s
		UpdateOrResetAugData(DATA.Ii_s, CHAINS.MaxLike_Ii_s, DATA, HOUSE);
		//// amend parameters with ModalPost_ParamVec
		AmendParams(CHAINS.MaxLike_ParamVec, CurrentPARAMS, DATA, HOUSE);
		//// generate MaxLike curves and attack rates and output s curves
		GenerateSurvivalCurves(DATA, CurrentPARAMS, HOUSE, SURVIVE, MAX_LIKE	, CHAINS, AUG, 
			/*AllPatients = */ true, /*OutputIndividualSurvivalTables = */ true	, /*OutputIndividual_Passive_SurvivalTables = */ false, 
			/*DType **** AttackRateContainer = */ SURVIVE.MeanModeAttackRates, /*ParamSet_SCurve_RootFilename = */ "MaxLike_SurvivalTable"	, /*ParamSet_PassiveSCurve_RootFilename = */ "MaxLike_PassiveSurvivalTable"	,
			/*Use_Default_AR_Post_Sample_No = */ false, /*CleanFirstDay_IndividualCurves =*/ true);
	
		////// Amend Augmented data with MaxLike_Ii_s (duplicate this line of code so below will still work even if you change stuff around). 
		UpdateOrResetAugData(DATA.Ii_s, CHAINS.MaxLike_Ii_s, DATA, HOUSE);
		//// amend parameters with MeanPostParamVec (with MaxLike_Ii_s / CHAINS.MeanPost_Ii_s are probabilities so they won't work). 
		AmendParams(CHAINS.MeanPostParamVec, CurrentPARAMS, DATA, HOUSE);
		//// generate curves, calculate attack rates. Don't output survival curves. Will still add to MeanModeAttackRates[][][][MEAN_POST] though which will be outputted using calls to  WriteSurvivalTable below. 
		GenerateSurvivalCurves(DATA, CurrentPARAMS, HOUSE, SURVIVE, MEAN_POST	, CHAINS, AUG, 
			/*AllPatients = */ true, /*OutputIndividualSurvivalTables = */ false, /* OutputIndividual_Passive_SurvivalTables = */ false,
			/*DType **** AttackRateContainer = */ SURVIVE.MeanModeAttackRates, /*ParamSet_SCurve_RootFilename = */"MeanPost_SurvivalTable ", /*ParamSet_PassiveSCurve_RootFilename = */"MeanPost_PassiveSurvivalTable", 
			/*Use_Default_AR_Post_Sample_No = */ false, /*Use_Default_AR_Post_Sample_No = */ true); 

		//// Write attack rates (summary statistics - mean and modal posterior estimates). Have to do this independently of WriteAllOutput function above as the above calls to GenerateSurvivalCurves function (to calculate mean and modal posterior S Curves and attack rates) also changes cumulative mean posterior curves.
		std::vector<string> DummyColNames(SURVIVE.Num_AR_SummaryStats, "");
		DummyColNames[MEAN_POST] = "Mean_Post"; DummyColNames[MODAL_POST] = "Modal_Post"; DummyColNames[MAX_LIKE] = "Max_Like"; 
		for (int TrialPhase = 0; TrialPhase < SURVIVE.HowManyTrialPhases; TrialPhase++)
			for (int DiseaseSeverity = 0; DiseaseSeverity < SURVIVE.HowManyDiseaseSeverities; DiseaseSeverity++)
				WriteSurvivalTable(CHAINS.OutputFolder + "AttackRates" + "_MeanMode" + HOUSE.OutputString + SURVIVE.TrialPhaseNames[TrialPhase] + SURVIVE.DiseaseNames[DiseaseSeverity] + ".txt",	//// filename
									SURVIVE.MeanModeAttackRates[TrialPhase][DiseaseSeverity],																					//// point to correct table (note MeanModeAttackRates rather than MetaAttackRates). 
									SURVIVE.NoSubjectCategories, SURVIVE.Num_AR_SummaryStats, SURVIVE.SurvivalTableRowNames, DummyColNames); 
	}

	//// Change everything back - necessary for restarting chains from where they were.
	UpdateOrResetAugData(DATA.Ii_s, CHAINS.Final_Iis, DATA, HOUSE);
	CurrentPARAMS.ParamVec = CHAINS.FinalParamVec;
	AmendParams(CurrentPARAMS.ParamVec, CurrentPARAMS, DATA, HOUSE);
	CurrentPARAMS.LikeFull = CHAINS.Final_L_Full;
	SetEqual_2D_Arrays(CurrentPARAMS.LikeParts, CHAINS.FinalLikeParts, HOUSE.TotalCountries, HOUSE.LCPerC);
	//// Reset Likelihood
	LikelihoodFromScratch(DATA, CurrentPARAMS, HOUSE, CurrentPARAMS.LikeParts);
}


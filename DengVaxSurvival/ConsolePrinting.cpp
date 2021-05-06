
#include "HeaderAndSwitches.h"
#include "ParamNumber.h"
#include "Splines.h"


string WhichLikeComp		(int ModuloComponent, const Housekeeping_Struct &HOUSE)
{
	string ComponentString = ""; 

			if (ModuloComponent >= HOUSE.LCPerC	) std::cerr << "ModuloComponent >= HOUSE.LCPerC: remmember to use Modulo Components"	<< endl; 
	else	if (ModuloComponent < 0				) std::cerr << "ModuloComponent < 0"				<< endl;
	else
	{
				if (ModuloComponent == HOUSE.L_Indices.BSs[SeroNeg])	ComponentString = "Aug_SNeg";
		else	if (ModuloComponent == HOUSE.L_Indices.BSs[SeroPos])	ComponentString = "Aug_SPos";
		else
		{
			//// rhos
			if (HOUSE.L_Indices.rhos != NULL)
				for (int serotype = 0; serotype < HOUSE.N_STypes; serotype++)
					for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
						if (ModuloComponent == HOUSE.L_Indices.rhos[serotype][PhaseSeverity])
							ComponentString = "rho_sero_" + to_string(serotype) + "_Phase_" + to_string(PhaseSeverity);

			//// Base_Haz
			if (HOUSE.L_Indices.BaseHaz != NULL)
				for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
					if (ModuloComponent == HOUSE.L_Indices.BaseHaz[PhaseSeverity])
						ComponentString = "BaseHaz_Phase" + to_string(PhaseSeverity);

			//// Age-specific hazard case multiplier
			if (HOUSE.L_Indices.AgeHazMult != NULL)
				for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
					if (ModuloComponent == HOUSE.L_Indices.AgeHazMult[PhaseSeverity])
						ComponentString = "AgeHazMult_Phase" + to_string(PhaseSeverity);

			//// Waningeffs
			if (HOUSE.L_Indices.WaningEffs != NULL)
				for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
					for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
						for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
							if (ModuloComponent == HOUSE.L_Indices.WaningEffs[serotype][BaselineSeroStatus][PhaseSeverity])
								ComponentString = "WanEff_sero_" + to_string(serotype) + "_SStatus_" + to_string(BaselineSeroStatus) + "_Phase_" + to_string(PhaseSeverity);

			//// IntBaseHaz
			if (HOUSE.L_Indices.IntBaseHaz != NULL)
				for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
					for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
						for (int TrialArm = 0; TrialArm < HOUSE.NumTrialArms_IBH; TrialArm++)
							if (ModuloComponent == HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus][TrialArm])
							ComponentString = "IntBaseHaz_Phase_" + to_string(PhaseSeverity) + "_SStatus_" + to_string(BaselineSeroStatus) + "_Arm_" + to_string(TrialArm);

			//// IntVacHaz
			if (HOUSE.L_Indices.IntVacHaz != NULL)
				for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
					for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
						for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
							if (ModuloComponent == HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype])
								ComponentString = "IntVacHaz_Phase_" + to_string(PhaseSeverity) + "_SStatus_" + to_string(BaselineSeroStatus) + "_sero_" + to_string(serotype);

			//// SeroSpec_Ks
			if (HOUSE.L_Indices.Ks != NULL)
				for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
					for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
						for (int K_like = 0; K_like < HOUSE.Num_K_Likes; K_like++)
							if (ModuloComponent == HOUSE.L_Indices.Ks[serotype][PhaseSeverity][K_like])
								ComponentString = "Ks_sero_" + to_string(serotype) + "_PS_" + to_string(PhaseSeverity) + "_Klike_" + to_string(K_like);;	

			//// Kplus
			if (HOUSE.L_Indices.Kplus != NULL)
				for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
					for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
						if (ModuloComponent == HOUSE.L_Indices.Kplus[serotype][PhaseSeverity])
							ComponentString = "Kplus_sero_" + to_string(serotype) + "_Phase_" + to_string(PhaseSeverity);

			//// l_BS_BaseHaz_Mults
			if (HOUSE.L_Indices.l_BS_BaseHaz_Mults != NULL)
				for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
					for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
						if (ModuloComponent == HOUSE.L_Indices.l_BS_BaseHaz_Mults[PhaseSeverity][BaselineSeroStatus])
							ComponentString = "l_BS_BaseHaz_Mult_" + to_string(PhaseSeverity) + "_SStatus_" + to_string(BaselineSeroStatus);

			//// KPrime_likes
			if (HOUSE.L_Indices.KPrimes != NULL)
				for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
					for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
						if (ModuloComponent == HOUSE.L_Indices.KPrimes[serotype][PhaseSeverity])
							ComponentString = "KPrimes_PS" + to_string(PhaseSeverity) +  "_sero_" + to_string(serotype);

			//// K_plus_Prime_likes
			if (HOUSE.L_Indices.KplusPrime != NULL)
				for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
					for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
						if (ModuloComponent == HOUSE.L_Indices.KplusPrime[serotype][PhaseSeverity])
							ComponentString = "KplusPrime_PS" + to_string(PhaseSeverity) + "_sero_" + to_string(serotype);
		}
	}
	if (ComponentString == "") std::cout << "WhichLikeComp error: ComponentString not determined" << endl; 

	return ComponentString; 
}

void PrintKsArray(int country, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
	{
		for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
		{
			for (int PrevInf = 0; PrevInf < HOUSE.Num_K_Params; PrevInf++)
				std::cerr << /*"C" << country <<*/ "K" + HOUSE.PhaseSeverity_strings[PhaseSeverity] + "_" + to_string(PrevInf) + "_s" + to_string(serotype) << ": " << PARAMS.K_s[country][PhaseSeverity][PrevInf][serotype] << ", ";
			std::cerr << endl;
		}
		std::cerr << endl;
	}
}
void printVector			(std::vector<int> VEC,						string VecName)
{
	std::cerr << VecName << " contains: ";	for (std::vector<int>::iterator it = VEC.begin(); it != VEC.end(); ++it) std::cerr << ' ' << *it; std::cerr << '\n';
}
void printVector			(std::vector<DType> VEC,					string VecName)
{
	std::cerr << VecName << " contains: "; for (std::vector<DType>::iterator it = VEC.begin(); it != VEC.end(); ++it) std::cerr << ' ' << *it; std::cerr << '\n';
}
void printVector			(std::vector<string>	VEC,				string VecName)
{
	std::cerr << VecName << " contains: "; for (std::vector<string>::iterator it = VEC.begin(); it != VEC.end(); ++it) std::cerr << ' ' << *it; std::cerr << '\n';
}
void print_LikeArray		(DType ** &LikeParts, const Housekeeping_Struct &HOUSE, int country)
{
	string ComponentString = "";
	for (int component = 0; component < HOUSE.LCPerC; component++)
	{
		ComponentString = WhichLikeComp(component, HOUSE);
		std::cout << "c" << country << " comp " << component << " " << ComponentString << " " << LikeParts[country][component] << endl;
	}
}
void PrintNaNComponents		(const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE)
{
	string ComponentString = ""; 
	for (int country = 0; country < HOUSE.TotalCountries; country++)
		for (int component = 0; component < HOUSE.LCPerC; component++)
		{
			ComponentString = WhichLikeComp(component, HOUSE);

			if (isnan(PARAMS.LikeParts[country][component]))		std::cerr << "lc" << component << " " << ComponentString << endl;
		}
}
void PrintUnequalComponents	(DType ** &CurrentLikeParts, DType ** &ProposedLikeParts, const Housekeeping_Struct &HOUSE)
{
	string ComponentString = "";
	DType AbsDiff_threshold = 0.0001;

	for (int country = 0; country < HOUSE.TotalCountries; country++)
		for (int component = 0; component < HOUSE.LCPerC; component++)
			if (	abs(CurrentLikeParts[country][component] - ProposedLikeParts[country][component]) > AbsDiff_threshold	)
			{
				ComponentString = WhichLikeComp(component, HOUSE);
				std::cout << country << "_" << component << " " << ComponentString << 
					" 1st " << CurrentLikeParts [country][component] << 
					" 2nd " << ProposedLikeParts[country][component] << endl;
			}
}
void Print_MCMC_Progress	(int iter, Chains_Struct &CHAINS, Params_Struct &CurrentPARAMS)
{
#ifdef USE_CLUSTER
	if ((iter + 1) % CHAINS.AddtoChainEveryHowManyIterations == 0)	std::cerr << " i" << iter << ", " << "l=" << CurrentPARAMS.LikeFull << endl;
#else
	//std::cerr << "i" << iter << " ";
	std::cerr << " i" << iter << ", " << "l=" << CurrentPARAMS.LikeFull << std::endl;
	//std::cerr << std::endl << " i" << iter << ", " << "l=" << CurrentPARAMS.LikeFull << " Temp " << CHAINS.Temp_ThisIteration;
	//std::cerr << " i" << iter << " BurnIn " << CHAINS.BurnIn << " Temp " << CHAINS.Temp_ThisIteration << /*" Slope " << Slope << " Intercept " << Intercept <<*/ endl;
#ifdef PRINT_LIKE_AT_EVERY_ITERATION
	std::cerr << "l=" << CurrentPARAMS.LikeFull << "   ";
#endif
#endif
}
void MaybePrint				(int param_no, DType Proposed_Param, Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE, bool ParamAcceptCondition)
{
#ifdef PRINT_LOADS
#ifndef USE_CLUSTER
	if (/*iter >= ITERATION_TO_PRINT_AFTER & */PRINT_PARAM_CONDITION)
	{
		std::cerr << /*"p" << param_no << " " << */CurrentPARAMS.NamesOfParameters[param_no] << endl;
		std::cerr << "C " << CurrentPARAMS.ParamVec[param_no] << " P " << Proposed_Param;
		std::cerr << " C_LL " << CurrentPARAMS.LikeFull << " P_LL " << ProposedPARAMS.LikeFull;

		DType AcceptanceProb = exp(ProposedPARAMS.LikeFull - CurrentPARAMS.LikeFull);
		std::cerr << " AccProb " << AcceptanceProb << endl;

		PrintNaNComponents(ProposedPARAMS, HOUSE);
		if (ParamAcceptCondition)  std::cerr << " ACCEPT" << endl; else std::cerr << " REJECT" << endl;
	}
#endif 
#endif
}
void PrintVacHazLikeDiscrep	(int param_no, Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, DType *** VacHazLike_CheckingArray)
{
	std::cerr	<< "p" << param_no << " " 
				<< CurrentPARAMS.NamesOfParameters	[param_no] << " current "
				<< CurrentPARAMS.ParamVec			[param_no] << " prop "		<< ProposedPARAMS.ParamVec[param_no] << endl;
	int country_dummy = 0, PhaseDummy = 0, SerostatusDummy = 1;
	std::cout <<	"e.g. P.VacHazLikes " << ProposedPARAMS.VacHazLikes[country_dummy][SerostatusDummy][PhaseDummy] << " " <<
					"e.g. C.VacHazLikes " << CurrentPARAMS. VacHazLikes[country_dummy][SerostatusDummy][PhaseDummy] << 
					" Checking " << VacHazLike_CheckingArray[country_dummy][SerostatusDummy][PhaseDummy]
					<< " %Diff " << 
					( ProposedPARAMS.VacHazLikes[country_dummy][SerostatusDummy][PhaseDummy] - VacHazLike_CheckingArray[country_dummy][SerostatusDummy][PhaseDummy]) 
					/ ProposedPARAMS.VacHazLikes[country_dummy][SerostatusDummy][PhaseDummy] << endl; 
}
void TestAndPrintParamNumberFucntions	(Housekeeping_Struct &HOUSE, Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS)
{
#ifdef TESTING_Is_and_Find_Param_FUNCTIONS

	std::cerr << endl << "Testing param_no_functions" << endl << endl;
	bool CurrPC_knot, PropPC_knot, C_HH, P_HH, C_q, P_q, C_Eff, P_Eff, C_RR, P_RR, C_rho, P_rho;
	for (int param_no = 0; param_no < HOUSE.No_Parameters; param_no++)
	{
		CurrPC_knot = IsParamAKnot(param_no, CurrentPARAMS.ParamNumbers);
		PropPC_knot = IsParamAKnot(param_no, ProposedPARAMS.ParamNumbers);
		C_HH = IsParamAHistHaz(param_no, CurrentPARAMS.ParamNumbers);
		P_HH = IsParamAHistHaz(param_no, ProposedPARAMS.ParamNumbers);
		C_q = IsParamA_qval(param_no, CurrentPARAMS.ParamNumbers);
		P_q = IsParamA_qval(param_no, ProposedPARAMS.ParamNumbers);
		C_rho = IsParamA_rho(param_no, CurrentPARAMS.ParamNumbers);
		P_rho = IsParamA_rho(param_no, ProposedPARAMS.ParamNumbers);
		C_Eff = IsParamAnEfficacy(param_no, CurrentPARAMS.ParamNumbers);
		P_Eff = IsParamAnEfficacy(param_no, ProposedPARAMS.ParamNumbers);
		C_RR = IsParamARelativeRisk(param_no, CurrentPARAMS.ParamNumbers);
		P_RR = IsParamARelativeRisk(param_no, ProposedPARAMS.ParamNumbers);

		std::cerr << "p" << param_no << " " << CurrentPARAMS.NamesOfParameters[param_no];

		if (CurrPC_knot)
		{
			if (CurrPC_knot != PropPC_knot) cerr << "KNOT ERROR";
			cerr << " ParamIsAKnot: country ";
			int knot = Find_Knot_FromKnotParam(param_no, CurrentPARAMS.ParamNumbers.Min_Knot, HOUSE.KnotsPerCountry);
			int country = Find_Country_FromKnotParam(param_no, CurrentPARAMS.ParamNumbers.Min_Knot, HOUSE.KnotsPerCountry, knot);
			cerr << country << " knot " << knot << endl;
		}
		if (C_HH)
		{
			if (C_HH != P_HH) cerr << "HH ERROR";
			cerr << " ParamIsAnHH: country ";
			int country = Find_Country_FromHistHazParam(param_no, CurrentPARAMS.ParamNumbers.Min_HHaz);
			cerr << country << endl;
		}
		if (C_q)
		{
			if (C_q != P_q) cerr << "q ERROR";
			cerr << " ParamIsAqval: country ";
			int qval = Find_qval_From_qParam(param_no, CurrentPARAMS.ParamNumbers.Min_qval, HOUSE.N_STypes);
			int country = Find_Country_From_qParam(param_no, CurrentPARAMS.ParamNumbers.Min_qval, HOUSE.N_STypes, qval);
			cerr << country << " qval " << qval << endl;
		}
		if (C_Eff)
		{
			if (C_Eff != P_Eff) cerr << "_Eff ERROR";
			cerr << " ParamIsAnEff: BaselineSeroStatus ";
			int BaselineSeroStatus = Find_SeroStatus_FromEffParam(param_no, CurrentPARAMS.ParamNumbers.Min_VacE, HOUSE.HowManySeroStatuses);
			int serotype = Find_Serotype_FromEffParam(param_no, CurrentPARAMS.ParamNumbers.Min_VacE, HOUSE.HowManySeroStatuses, BaselineSeroStatus);
			cerr << BaselineSeroStatus << " serotype " << serotype << endl;
		}
		if (C_RR)
		{
			if (C_RR != P_RR) cerr << "_RR ERROR";
			cerr << " ParamIsAnRR: PhaseSeverity ";
			int PrevInf = Find_PrevInf_From_K_Param(param_no, CurrentPARAMS.ParamNumbers.Min_K, HOUSE.Num_K_Params);
			int PhaseSeverity = Find_PhaseSeverity_From_K_Param(param_no, CurrentPARAMS.ParamNumbers.Min_K, HOUSE.Num_K_Params, HOUSE.HowManyCaseCategories, PrevInf);
			int serotype = Find_Serotype_From_K_Param(param_no, CurrentPARAMS.ParamNumbers.Min_K, HOUSE.Num_K_Params, HOUSE.HowManyCaseCategories, PrevInf, PhaseSeverity);
			cerr << PhaseSeverity << " PrevInf " << PrevInf << " serotype " << serotype << endl;
		}
		if (C_rho)
		{
			if (C_rho != P_rho) cerr << "_rho ERROR";
			cerr << " ParamIsARho: country ";
			int rho = Find_rho_From_rhoParam(param_no, CurrentPARAMS.ParamNumbers.Min_Rho, HOUSE.N_STypes);
			int country = Find_Country_From_rhoParam(param_no, CurrentPARAMS.ParamNumbers.Min_Rho, HOUSE.N_STypes, rho);
			cerr << country << " rho " << rho << endl;
		}
	}
#endif
}
void PrintFittedParamNumbers			(const Housekeeping_Struct &HOUSE, Params_Struct &CurrentPARAMS)
{
#ifdef OUTPUT_FITTED_PARAM_NOS
	std::cerr << endl << "ParamNosYouWillFit Size " << CurrentPARAMS.ParamNosYouWillFit.size() << endl;
	for (int param_no = 0; param_no < CurrentPARAMS.ParamNosYouWillFit.size(); param_no++)
		if (!IsParamAKnot(CurrentPARAMS.ParamNosYouWillFit[param_no], CurrentPARAMS.ParamNumbers))
			if (!IsParamAHistHaz(CurrentPARAMS.ParamNosYouWillFit[param_no], CurrentPARAMS.ParamNumbers))
				std::cerr << "p" << CurrentPARAMS.ParamNosYouWillFit[param_no] << " " << CurrentPARAMS.NamesOfParameters[CurrentPARAMS.ParamNosYouWillFit[param_no]] << " " << CurrentPARAMS.ParamVec[CurrentPARAMS.ParamNosYouWillFit[param_no]]
				<< " range " << CurrentPARAMS.ParamRanges[LowerBound][CurrentPARAMS.ParamNosYouWillFit[param_no]] << " to " << CurrentPARAMS.ParamRanges[UpperBound][CurrentPARAMS.ParamNosYouWillFit[param_no]] << endl;
#endif
}
void Print_Like_Indices					(Housekeeping_Struct &HOUSE, Params_Struct &CurrentPARAMS)
{
#ifdef OUTPUT_L_INDICES
	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
		std::cerr << "BaseHaz_Phase_" << PhaseSeverity << "  " << HOUSE.L_Indices.BaseHaz[PhaseSeverity] << endl;

	std::cerr << "Aug_SPos				" << HOUSE.L_Indices.Aug_SPos << endl;
	std::cerr << "Aug_SNeg				" << HOUSE.L_Indices.Aug_SNeg << endl;

	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			for (int TrialArm = 0; TrialArm < HOUSE.NumTrialArms_IBH; TrialArm++)
				std::cerr << "IntBaseHaz_Phase_" + to_string(PhaseSeverity) + "_SStatus_" + to_string(BaselineSeroStatus) + "_TrialArm_" + to_string(TrialArm) << "   " << HOUSE.L_Indices.IntBaseHaz[PhaseSeverity][BaselineSeroStatus][TrialArm] << endl; ;

	for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int PrevInf = 0; PrevInf < HOUSE.Num_K_Likes; PrevInf++)
				std::cerr << "SeroSpec_Ks_sero_" << serotype << "_Phase_" << PhaseSeverity << "_PrevInf_" << PrevInf << "     " << HOUSE.L_Indices.Ks[serotype][PhaseSeverity][PrevInf] << endl;;

	if (!(HOUSE.ModelVariant == SIMPLE_NUMERICAL))
		for (int serotype = 0; serotype < HOUSE.N_STypes_Ks; serotype++)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
				std::cerr << "Kplus_sero_" << serotype << "_Phase_" << PhaseSeverity << "     " << HOUSE.L_Indices.Kplus[serotype][PhaseSeverity] << endl;;

	for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
		for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
			for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
				std::cerr << "IntVacHaz_Phase_" << PhaseSeverity << "_SStatus_" << BaselineSeroStatus << "_sero_" << serotype << "     " << HOUSE.L_Indices.IntVacHaz[PhaseSeverity][BaselineSeroStatus][serotype] << endl;

	for (int serotype = 0; serotype < HOUSE.N_STypes_VEs; serotype++)
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
				std::cerr << "WaningEffs_sero_" << serotype << "_SStatus_" << BaselineSeroStatus << "_Phase_" << PhaseSeverity << "   " << HOUSE.L_Indices.WaningEffs[serotype][BaselineSeroStatus][PhaseSeverity] << endl;

	if (HOUSE.SeroSpecificEfficacies || HOUSE.SeroSpecific_K_values)
		for (int serotype = 0; serotype < HOUSE.N_STypes; serotype++)
			for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
				std::cerr << "rho_sero_" << serotype + "_PhaseSeverity_" + to_string(PhaseSeverity) << "     " << HOUSE.L_Indices.rhos[serotype][PhaseSeverity] << endl;

	if (HOUSE.AdjHaz)
		for (int PhaseSeverity = 0; PhaseSeverity < HOUSE.HowManyCaseCategories; PhaseSeverity++)
			for (int BaselineSeroStatus = 0; BaselineSeroStatus < HOUSE.HowManySeroStatuses; BaselineSeroStatus++)
				std::cerr << "l_BS_BaseHaz_Mults_Phase_" << PhaseSeverity << "_SStatus_" << BaselineSeroStatus << "   " << HOUSE.L_Indices.l_BS_BaseHaz_Mults[PhaseSeverity][BaselineSeroStatus] << endl;

	std::cerr << endl << "LCPerC " << HOUSE.LCPerC << " LikeSize " << CurrentPARAMS.LikeSize << endl;
#endif
}
void Check_Age_ParamsEtc				(Age_Option AS_VE_Wane_Haz, const Housekeeping_Struct &HOUSE, DType * ParamContainer, DType * MultContainer, DType ** CoeffContainer, int NumFunctionParams, int PolysPerSpline, int MaxSplineDegree, Params_Struct &PARAMS)
{
	if (AS_VE_Wane_Haz != Age_Option::INDEPENDENT)
	{
		//// print Params
		for (int AParam = 0; AParam < NumFunctionParams; AParam++)
			std::cout << "Age_Param_" << AParam << " " << ParamContainer[AParam] << endl;
		//// print Mults (note: for waning, multipliers have a time part as well - you don't store the mean durations / rates in a container, you calculate afresh then reset and store the waning multipliers by age and daypostdose. Hence this printing will be garbage)
		for (int age = 0; age < HOUSE.HowManyAges; age++)
			std::cout << "Age_Mult_age_" << age << " " << MultContainer[age] << endl;
		//// print Spline coefficients and a couple of spline / cubic values
		if (AS_VE_Wane_Haz != Age_Option::CATEGORICAL)
			for (int poly = 0; poly < PolysPerSpline; poly++)
				for (int coeff = 0; coeff <= MaxSplineDegree; coeff++)
					//// coefficients
					std::cout << "Age_poly_" << poly << "_Coeff_" << coeff << " " << CoeffContainer[poly][coeff] << endl;
		//// print couple of spline / cubic values (bit dodgy this - too lazy to add in proper function arguments, so just using values for VE, as they're usually same as Waning, AS_Haz etc). Won't be true if doing baseline hazards
		if (AS_VE_Wane_Haz == Age_Option::SPLINE)
			for (int xval = 4; xval < 13; xval++)
				std::cout << " y " << Spline(xval, PARAMS.Age_Effs_xKnots[0], ParamContainer, CoeffContainer, HOUSE.KnotsPerSpline_EffMultiplier, PolysPerSpline);
		else if (AS_VE_Wane_Haz == Age_Option::CUBIC)
			for (int xval = 0; xval < 13; xval++)
				std::cout << " y " << Cubic(xval, PARAMS.Age_Effs_xKnots[0], ParamContainer, CoeffContainer, HOUSE.KnotsPerSpline_EffMultiplier, PolysPerSpline);

		std::cout << endl;
	}
}
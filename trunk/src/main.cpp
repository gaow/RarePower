// =====================================================================================
//
//       Filename:  main.cpp
//
//    Description:  simulations, sampling, and power calculations for gene based association studies.
//
//        Version:  1.0
//        Created:  01/04/2011 10:52:07 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Gao Wang (gw), wangow@gmail.com
//                  Baylor College of Medicine, Texas, USA
//        License:  GNU General Public License <http://www.gnu.org/licenses/>
//                  Copyright (c) 2011, Gao Wang
//
// =====================================================================================

#include <cstdlib>
#include <fstream>
#include <cmath>
#include <numeric>

#include "gw_utilities.h"
#include "person.h"
#include "assocsims.h"
#include "Argument_helper.h"
#include "main.h"
using namespace std;
using namespace gpow;

int main(int argc, const char * argv[])
{

	string program_name = string(EXENAME);
	bool verbose = false, quiet = false;

	//////
	// Parameters
	//////

	unsigned seed = 10086;
	string gFile("None");
	bool shouldUseGenPool = false;
	double boundary = 0.01;
	double neutral_cutoff = 0.0;

	vector<double> propFunctionalRv(2);
	propFunctionalRv[0] = 1.0;
	//!- proportion of effective deleterious variant (vs. non-causal)
	propFunctionalRv[1] = 1.0;
	//!- proportion of effective protective variant (vs. non-causal)
	string strmoi = "A";
	char moi = strmoi[0];
	char tmoi = moi;

	string simulationTask = "None";
	//!- options: 1.dichot-odds, 3.dichot-par, 4.qt, 5.dichot-qt, 2.pop-odds, 6.mendelian

	vector<double> oddsRatios(5);
	oddsRatios[0] = 0.0;
	oddsRatios[1] = 2.0;
	//!- odds ratio for deleterious variants
	oddsRatios[2] = 0.0;
	oddsRatios[3] = 1.0;
	//!- odds ratio for protective variants
	oddsRatios[4] = 1.0;
	//!- odds ratio for common variants
	double baselinef = 0.01;
	//!- base-line penetrance ~ prevalence

	vector<double> pars(2);
	pars[0] = 0.05;
	pars[1] = 0.0;
	bool isParVariable = false;
	//!- population attributable risk model

	vector<double> qtcoefs(3);
	//!- Locus specific effects w.r.t. standard deviation
	qtcoefs[0] = 0.0;
	qtcoefs[1] = 0.2;
	//!- effect per rare variant
	qtcoefs[2] = 0.0;
	//!- effect per common variant
	vector<double> qtcuts(2);
	qtcuts[0] = 0.1;
	qtcuts[1] = 0.9;
	bool shouldMarkBin = false;

	//!- For Mendelian simulations
	double percentageCausal = 0.1;
	bool isMendelAlleleFixed = false;
	double propHeterCases = 0;

	vector<double> propMissingData(4);
	propMissingData[0] = 0.0;
	propMissingData[1] = 0.0;
	propMissingData[2] = 0.0;
	propMissingData[3] = 0.0;
	double missingLowMaf = 0.0;
	bool shouldMarkMissing = false;

	unsigned nPopulation = 2000;
	unsigned nCases = 500;
	unsigned nCtrls = 500;
	unsigned nUnphenotyped = 0;
	vector<double> pedInfos(6, 0.0);

	bool isPedWritten = false;

	/*** Write simulated data ***/
	bool isSynoKept = false;
	bool isCvTrimmed = false;
	string projectName = "proj" + n2s(time(NULL));

	/*** settings for analysis ***/
	//!-Trim data parameters
	string test = "CMC-one-midP";
	double mafLower = 0.0;
	double mafUpper = 0.01;

	double alpha = 0.05;
	unsigned nPermutations = 2000;
	unsigned adaptive = 5000;
	unsigned nReplicates = 1000;


	//////
	// argsparser
	//////

	bool noinfo = false;
	for (int i = 0; i < argc; ++i) {
		if (strcmp(argv[i], "--help") != 0 && strcmp(argv[i], "-h") == 0) {
			noinfo = true;
			break;
		}
	}

	dsr::Argument_helper ah;
	// the only required argument
	ah.new_string("task", args_dsc("task", noinfo), simulationTask);
	ah.new_string("gdata", args_dsc("gdata", noinfo), gFile);
	ah.new_optional_string("pname", args_dsc("pname", noinfo), projectName);
	// named arguments
	ah.new_named_double('f', "define_rare", "<frequency>", args_dsc("f", noinfo), boundary);
	ah.new_named_double('q', "prop_func_deleterious", "<fraction>", args_dsc("q", noinfo), propFunctionalRv[0]);
	ah.new_named_double('p', "prop_func_protective", "<fraction>", args_dsc("p", noinfo), propFunctionalRv[1]);
	ah.new_named_string('g', "mode_of_inheritance", "<moi>", args_dsc("g", noinfo), strmoi);
	ah.new_named_double('A', "OR_deleterious_min", "<effect_size>", args_dsc("A", noinfo), oddsRatios[0]);
	ah.new_named_double('B', "OR_deleterious_max", "<effect_size>", args_dsc("B", noinfo), oddsRatios[1]);
	ah.new_named_double('C', "OR_protective_min", "<effect_size>", args_dsc("C", noinfo), oddsRatios[2]);
	ah.new_named_double('D', "OR_protective_max", "<effect_size>", args_dsc("D", noinfo), oddsRatios[3]);
	ah.new_named_double('E', "OR_common", "<effect_size>", args_dsc("E", noinfo), oddsRatios[4]);
	ah.new_named_double('F', "prevalence", "<fraction>", args_dsc("F", noinfo), baselinef);
	ah.new_named_double('G', "PAR_deleterious", "<fraction>", args_dsc("G", noinfo), pars[0]);
	ah.new_named_double('H', "PAR_protective", "<fraction>", args_dsc("H", noinfo), pars[1]);
	ah.new_flag('I', "PAR-variable", args_dsc("I", noinfo), isParVariable);
	ah.new_named_double('J', "QT_effect_min", "<multiplier>", args_dsc("J", noinfo), qtcoefs[0]);
	ah.new_named_double('K', "QT_effect_max", "<multiplier>", args_dsc("K", noinfo), qtcoefs[1]);
	ah.new_named_double('L', "QT_effect_common", "<multiplier>", args_dsc("L", noinfo), qtcoefs[2]);
	ah.new_named_double('M', "QT_lower_percentile", "<fraction>", args_dsc("M", noinfo), qtcuts[0]);
	ah.new_named_double('N', "QT_upper_percentile", "<fraction>", args_dsc("N", noinfo), qtcuts[1]);
	ah.new_flag('O', "QT-binary", args_dsc("O", noinfo), shouldMarkBin);
	ah.new_named_double('P', "Mendelian_causal", "<fraction>", args_dsc("P", noinfo), percentageCausal);
	ah.new_named_double('Q', "Mendelian_heterogeneity", "<fraction>", args_dsc("Q", noinfo), propHeterCases);
	ah.new_flag('R', "fixed_Mendelian_variant", args_dsc("R", noinfo), isMendelAlleleFixed);
	ah.new_named_unsigned_int('W', "num_all_samples", "<#samples>", args_dsc("W", noinfo), nPopulation);
	ah.new_named_unsigned_int('X', "num_cases", "<#cases>", args_dsc("X", noinfo), nCases);
	ah.new_named_unsigned_int('Y', "num_ctrls", "<#ctrls>", args_dsc("Y", noinfo), nCtrls);
	ah.new_named_unsigned_int('Z', "num_cohort_ctrls", "<#cohort_ctrls>", args_dsc("Z", noinfo), nUnphenotyped);
	ah.new_flag('U', "use_haplotype_pool", args_dsc("U", noinfo), shouldUseGenPool);
	ah.new_named_double('a', "prop_missing_deleterious", "<fraction>", args_dsc("a", noinfo), propMissingData[0]);
	ah.new_named_double('b', "prop_missing_protective", "<fraction>", args_dsc("b", noinfo), propMissingData[1]);
	ah.new_named_double('c', "prop_missing_non_causal", "<fraction>", args_dsc("c", noinfo), propMissingData[2]);
	ah.new_named_double('d', "prop_missing_synonymous", "<fraction>", args_dsc("d", noinfo), propMissingData[3]);
	ah.new_named_double('e', "missing_low_maf", "<frequency>", args_dsc("e", noinfo), missingLowMaf);
	ah.new_flag('k', "recode_missing", args_dsc("k", noinfo), shouldMarkMissing);
	ah.new_flag('i', "keep_synonymous", args_dsc("i", noinfo), isSynoKept);
	ah.new_flag('j', "remove_common_loci", args_dsc("j", noinfo), isCvTrimmed);
	ah.new_named_double('l', "maf_lower", "<frequency>", args_dsc("l", noinfo), mafLower);
	ah.new_named_double('m', "maf_upper", "<frequency>", args_dsc("m", noinfo), mafUpper);
	ah.new_named_double('n', "define_neutral", "<annotation_cutoff>", args_dsc("n", noinfo), neutral_cutoff);
	ah.new_named_string('t', "test", "<association_test>", args_dsc("t", noinfo), test);
	ah.new_named_double('s', "significance", "<alpha_level>", args_dsc("s", noinfo), alpha);
	ah.new_named_unsigned_int('r', "replicates", "<#replicates>", args_dsc("r", noinfo), nReplicates);
	ah.new_named_unsigned_int('u', "permutations", "<#permutations>", args_dsc("u", noinfo), nPermutations);
	ah.new_named_unsigned_int('y', "rng_seed", "<long_integer>", args_dsc("y", noinfo), seed);
	ah.new_flag('z', "simulation_only", args_dsc("z", noinfo), isPedWritten);
	ah.new_flag('v', "maximal_output", args_dsc("v", noinfo), verbose);
	ah.new_flag('x', "minimal_output", args_dsc("x", noinfo), quiet);

	// program information
	ah.set_name(program_name.c_str());
	ah.set_description(banner.c_str());
	ah.set_version((atof(VERSION) > 0.0) ? VERSION : SVN_REV);
	ah.set_author("Gao Wang <wangow@gmail.com>");
	ah.set_build_date(COMPILE_DATE);

	ah.process(argc, argv);

	if (quiet) verbose = false;
	if (verbose) ah.write_usage(std::clog, 1);

	//////
	// Check options and generate command in effect.
	//////
	strmoi = pystring::upper(strmoi);
	string cmdcurrent = check_options(program_name, projectName, gFile, boundary, neutral_cutoff,
		propFunctionalRv,  strmoi,  simulationTask, oddsRatios,  baselinef, pars,
		isParVariable, qtcoefs,  qtcuts,  shouldMarkBin,  percentageCausal,
		isMendelAlleleFixed, propMissingData,  missingLowMaf,  shouldMarkMissing,  nCases,
		nCtrls,  nPopulation,  nUnphenotyped,  propHeterCases,  isSynoKept,
		isCvTrimmed,  isPedWritten,  test,  mafLower, mafUpper,  alpha,  nPermutations,
		nReplicates,  verbose, quiet, seed, shouldUseGenPool);
	// adjust moi
	moi = strmoi[0];
	tmoi = (strmoi.size() == 2) ? strmoi[1] : moi;
	if (!quiet) {
		std::clog << "INFO: Genetic model/data: " << gFile << std::endl;
		if (!isPedWritten) {
			std::clog << "INFO: Association method: " << test << std::endl;
		}
		std::clog << "INFO: Running command $" << pystring::replace(cmdcurrent, "PLACEHOLDER", "\"" + test + "\"") << std::endl;
	}


	//////
	// Power calculations
	//////

	RNG rng;
	gsl_rng * gslr;
	if (seed == 0) {
		gslr = rng.get();
	} else{
		gslr = rng.get(seed);
	}

	double yl = gsl_cdf_ugaussian_Pinv(qtcuts[0]);
	double yh = gsl_cdf_ugaussian_Pinv(qtcuts[1]);

	vector2F mafDat;
	vector2F annDat;
	vector2UI posDat;
	scan_vector2F(gFile + ".ann", annDat);
	scan_vector2F(gFile + ".maf", mafDat);
	scan_vector2UI(gFile + ".pos", posDat);

	string pvalueFileName = projectName + ".pvalues";
	ofstream outPvalue;
	if (!isPedWritten && verbose) {
		outPvalue.open(pvalueFileName.c_str(), ios::app);
		if (!is_file_empty(pvalueFileName)) {
			std::cerr << "WARNING: project summary files for [ " << projectName << " ] already exist; information from this simulation will be appended to the end of existing files!" << std::endl;
		}
		outPvalue << "$ " << pystring::replace(cmdcurrent, "PLACEHOLDER", "\"" + test + "\"") << std::endl;
	}

	unsigned iReplicate = 0;
	vector<string> tests;
	pystring::split(test, tests);
	vectorUI pcounts(tests.size(), 0);


	while (iReplicate != nReplicates) {

		/*** Choose MAF data ***/
		unsigned dataIdx = gsl_rng_uniform_int(gslr, mafDat.size());
		vector<double> & mafs = mafDat[dataIdx];
		vector<double> & fnctAnnotations = annDat[dataIdx];
		vector<unsigned> & positions = posDat[dataIdx];
		vector2UI poolDat;
		if (shouldUseGenPool) {
			scan_vector2UI(gFile + "_hap/" + gFile + ".hap" + n2s(dataIdx + 1), poolDat);
			assert(poolDat.front().size() == mafs.size());
		} else {
			poolDat.resize(0);
		}

		vector2F genoFreqs(mafs.size());

		for (unsigned i = 0; i != mafs.size(); ++i) {
			genoFreqs[i].push_back( (1 - mafs[i]) * (1 - mafs[i]) );
			genoFreqs[i].push_back(mafs[i] * mafs[i]);
		}

		/*** Generate population data ***/
		gwSimulator * simulator = new gwSimulator(boundary, neutral_cutoff, pedInfos, mafs, genoFreqs, fnctAnnotations, positions);

		if (simulationTask == "6")
			simulator->createGenotypeMendelianTraitsAssociations(percentageCausal, (!isMendelAlleleFixed),
				moi, nCases, propHeterCases, nCtrls, gslr);  // "mendelian"

		else if (simulationTask == "7")
			simulator->createGenotypeComplexTraitsAssociations(propFunctionalRv, baselinef, oddsRatios,
				moi, nCases, gslr, poolDat);  // "sibpairs"

		else if (simulationTask == "1")
			simulator->createGenotypeComplexTraitsAssociations(propFunctionalRv, baselinef, oddsRatios,
				moi, nCases, nCtrls, nUnphenotyped, gslr, verbose, projectName, poolDat);  // "dichot-odds"

		else if (simulationTask == "2")
			simulator->createGenotypeComplexTraitsAssociations(propFunctionalRv, baselinef, oddsRatios,
				moi, nPopulation, gslr, verbose, projectName, poolDat);  // "pop-odds"

		else if (simulationTask == "3")
			simulator->createGenotypeComplexTraitsAssociations(propFunctionalRv, pars, (!isParVariable), moi,
				nCases, nCtrls, nUnphenotyped, gslr);  // "dichot-par"

		else if (simulationTask == "8")
			simulator->createGenotypeComplexTraitsAssociations(propFunctionalRv, baselinef, pars, (!isParVariable),
				moi, nCases, nCtrls, nUnphenotyped, gslr, verbose, projectName);  // "dichot-par-odds"

		else if (simulationTask == "4")
			simulator->createGenotypeComplexTraitsAssociations(propFunctionalRv, qtcoefs, moi, nPopulation, gslr, poolDat);  // "qt"

		else if (simulationTask == "5")
			simulator->createGenotypeComplexTraitsAssociations(propFunctionalRv, qtcoefs, moi, qtcuts, shouldMarkBin,
				nPopulation, nCases, nCtrls, nUnphenotyped, gslr, poolDat);  // "dichot-qt"
		else {
			exit(-1);
		}


		/*** Exclude some variants (a mimic genotyping procedure) ***/
		simulator->mimicGenotyping(propMissingData, missingLowMaf, shouldMarkMissing, gslr);
		simulator->createPedfileMatrix((!isSynoKept), isCvTrimmed, isPedWritten, projectName, simulationTask);

		/*** Create data matrix, write to a ped file if necessary ***/
		if (isPedWritten) {
			break;
		}

		//    if (verbose) {
		//      simulator->calcVariantsPars(test, mafLower, mafUpper);
		//    }

		vector2F genotypes = simulator->getGenotypes();
		vectorF phenotypes = simulator->getPhenotypes();
		vectorF observedMafs = simulator->getObservedMafs();

		delete simulator;

		/***Do p-value calculations***/
		unsigned pt = 0;
		while (true) {
			if (pt == tests.size()) {
				if (pt == 0) {
					exit(-1);
				}else {
					break;
				}
			}

			// set data
			gwAssocdata assocdat;
			assocdat.setXdat(genotypes);
			assocdat.setYdat(phenotypes);
			assocdat.setMafs(observedMafs, mafLower, mafUpper);
			if (simulationTask == "7") {
				// mask out loci where case/ctrl are concordent
				assocdat.markwildSibpairloci();
				// set permutator for sibpairs permutation
				assocdat.setPermutator(1);
			}
			assocdat.trimXdat();
			// set test
			gwBaseTest * atest = testFactory(tests[pt]);
			if (!atest) {
				std::clog << "WARNING: Invalid test method [ " << tests[pt] << " ]" << endl;
				tests.erase(tests.begin() + pt);
				pcounts.erase(pcounts.begin() + pt);
				delete atest;
				continue;
			}
			atest->setA(alpha);
			atest->setMOI(tmoi);
			atest->setAlternative((tests[pt].find("one") < tests[pt].size()) ? 1 : 2);
			atest->setPermutations(nPermutations, adaptive);
			if (tests[pt].find("ExtremeQT") < tests[pt].size()) {
				atest->setVar("yh", yh);
				atest->setVar("yl", yl);
			}
			if (tests[pt].find("midP") < tests[pt].size()) {
				atest->useMidP();
			}
			atest->applyMOI(assocdat);
			double pvalue = atest->apply(assocdat);
			delete atest;

			if (pvalue <= alpha) ++pcounts[pt];
			if (!isPedWritten && verbose) {
				outPvalue << tests[pt] << "\t" << pvalue << "\n";
				if (iReplicate + 1 == nReplicates) {
					outPvalue << std::endl;
					outPvalue.close();
				}
			}
			++pt;
		}

		++iReplicate;
		if (!quiet) {
			progress_bar(iReplicate, nReplicates, "Replicate");
		}
	}

	if (isPedWritten) {
		if (!quiet) std::clog << "INFO: Simulated data written to files " << projectName + ".geno " << projectName + ".phen " << projectName + ".map " << std::endl;
		std::ofstream lout((projectName + ".log").c_str(), ios::app);
		lout << "cmd = " << cmdcurrent << std::endl;
		lout.close();
	} else {
		if (!quiet) std::clog << std::endl << "INFO: output format header [ METHOD|POWER|STANDARD.ERROR|PROJECT_ID|COMMAND ]" << std::endl;
		for (unsigned i = 0; i != tests.size(); ++i) {
			double power = 1.0 * pcounts[i] / (1.0 * nReplicates);
			double pse = sqrt(power * (1.0 - power) / (1.0 * nReplicates));
			gw_round(power, 1e-4);
			gw_round(pse, 1e-6);
			std::cout << tests[i] << "|" << power << "|" << pse << "|" << projectName << "|" << pystring::replace(cmdcurrent, "PLACEHOLDER", tests[i]) << std::endl;
		}
	}
	return 0;
}


/////////////////
/////////////////


std::string check_options(std::string prog_name, std::string & projectName, std::string & gFile, double & boundary, double & neutral_cutoff,
                          vectorF & propFunctionalRv, std::string & strmoi, std::string & simulationTask, vectorF & oddsRatios, double & baselinef, vectorF & pars,
                          bool & isParVariable, vectorF & qtcoefs, vectorF & qtcuts, bool & shouldMarkBin, double & percentageCausal,
                          bool & isMendelAlleleFixed, vectorF & propMissingData, double & missingLowMaf, bool & shouldMarkMissing, unsigned & nCases,
                          unsigned & nCtrls, unsigned & nPopulation, unsigned & nUnphenotyped, double & propHeterCases, bool & isSynoKept,
                          bool & isCvTrimmed, bool & isPedWritten, std::string & test, double & mafLower,
                          double & mafUpper, double & alpha, unsigned & nPermutations, unsigned & nReplicates, bool & verbose, bool & quiet, unsigned & seed,
                          bool & shouldUseGenPool)
{
	std::string cmd = prog_name;

	vector<string> tests;
	pystring::split(test, tests);
	// check if tests are compatible under the same simulation setting
	bool isQT = false, isCC = false, need_perm = false;
	for (unsigned i = 0; i != tests.size(); ++i) {
		if (tests[i].find("QT") < tests[i].size()) {
			isQT = true;
		}else {
			isCC = true;
		}
		if (tests[i] != "CMC" && tests[i] != "WSS" && tests[i] != "RVE"
		    && tests[i] != "CMC-one" && tests[i] != "CMC-one-midP" && tests[i] != "WSS-one" 
            && tests[i] != "RVE-one" && tests[i] != "RVE-one-midP"
		    && tests[i] != "MZQT" && tests[i] != "MZQT-one") {
			need_perm = true;
		}
	}

	if (isQT && isCC) {
		std::cerr << "ERROR: Methods [ " << test << " ] not compatible under the same simulation setting" << std::endl;
		std::cerr << "ERROR: Quit now" << std::endl;
		exit(-1);
	}

	//////
	// check options
	//////
	if (strmoi.size() > 2) {
		strmoi = pystring::slice(strmoi, 0, 2);
		std::cerr << "WARNING: Invalid -g/--mode_of_inheritance input. Set it to " << strmoi << std::endl;
	}

	if (!(simulationTask == "1" || simulationTask == "2" || simulationTask == "3" || simulationTask == "4" || simulationTask == "5" || simulationTask == "6" || simulationTask == "7" || simulationTask == "8")) {
		std::cerr << "ERROR: Invalid simulation task input [ " << simulationTask << " ] . Choose task = 1~7 (see --help)" << std::endl;
		std::cerr << "ERROR: Quit now" << std::endl;
		exit(-1);
	}
	if (simulationTask == "6") { // Mendelian simulations
		if (isCvTrimmed == false || mafUpper < 1.0) {
			std::clog << "WARNING: options \"-m 1 -j\"(or \"--maf_upper 1 --remove_common_loci\") should be used for analysis of Mendelian diseases associations" << std::endl;
			std::clog << "WARNING: setting parameters \"-m 1 -j\"" << std::endl;
			isCvTrimmed = true;
			mafUpper = 1.0;
		}
		if (strmoi[0] == 'A') {
			strmoi[0] = 'D';
			std::clog << "WARNING: setting default mode of inheritance to \"dominant\" ('D')." << std::endl;
		}
		if (strmoi[0] == 'C' && isMendelAlleleFixed) {
			std::clog << "WARNING: using Compound Recessive model without allelic heterogeneity. This is essentially the Recessive model without allelic heterogeneity" << std::endl;
		}

	}
	if (simulationTask == "5") {
		if (nPopulation > 0 && (nCases > 0 || nCtrls > 0)) {
			std::clog << "WARNING: sampling extreme QT values from finite cohort size = " << nPopulation << ". Setting samples having values > " << qtcuts[1] * 100 << " percentile as cases, < " << qtcuts[0] * 100 << " percentile as ctrls" << std::endl;
			std::clog << "WARNING: to sample from population for fixed #case/ctrls, please set \"-W 0 (or --num_all_samples 0)\" and specify number of cases/ctrls" << std::endl;
			nCases = 0;
			nCtrls = 0;
		}
		if (!isPedWritten && test.find("QT") > test.size() && shouldMarkBin != true) {
			std::clog << "WARNING: recoding extreme QT values to binary variables in order to be analyzed by method [ " << test << " ]" << std::endl;
			shouldMarkBin = true;
		}
	}
	if (simulationTask == "4") {
		for (size_t i = 0; i < tests.size(); ++i) {
			if (tests[i].find("QT") > tests[i].size()) {
				std::clog << "WARNING: test method [ " << test << " ] is not valid for quantitative traits analysis. Setting test to [ MZQT ]" << std::endl;
				tests[i] = "MZQT";
				tests.resize(1);
				break;
			}
		}
	}
	if (simulationTask == "1" || simulationTask == "2" || simulationTask == "7") {
		if (!((oddsRatios[0] >= 1.0 || oddsRatios[0] == 0.0)
		      && oddsRatios[0] < oddsRatios[1] && oddsRatios[1] >= 1.0 && (oddsRatios[2] < 1.0 && oddsRatios[2] >= 0.0)
		      && oddsRatios[2] < oddsRatios[3] && oddsRatios[3] <= 1.0 && oddsRatios[4] > 0.0)) {
			std::cerr << "ERROR: Invalid odds ratio input" << std::endl;
			std::cerr << "ERROR: Quit now" << std::endl;
			exit(-1);
		}
	}
	if (simulationTask == "3") {
		if (!(pars[0] >= 0.0 && pars[1] >= 0.0 && pars[0] < 1.0 && pars[1] < 1.0)) {
			std::cerr << "ERROR: Invalid population attributable risks input" << std::endl;
			std::cerr << "ERROR: Quit now" << std::endl;
			exit(-1);
		}
	}
	if (simulationTask == "7") {
		if (nCases != nCtrls) {
			std::clog << "WARNING: setting number of unaffeceted sibling to be the same as affected, ie, " << nCases << " case/ctrls." << std::endl;
			nCtrls = nCases;
		}
	}

	if (isPedWritten && !is_file_empty(projectName + ".ped")) {
		std::cerr << "WARNING: data files [ " << projectName << ".ped/pos/log ] already exist; data from this simulation will be appended to the end of existing files!" << std::endl;
	}

	if (shouldUseGenPool && (simulationTask == "3" || simulationTask == "6")) {
		std::clog << "WARNING: Cannot use haplotype pool for given task code [ " << simulationTask << " ]. Disabling the use of external pool." << std::endl;
		shouldUseGenPool = false;
	}

	// reformat command in case the input has extra unwanted spaces
	test = pystring::join(" ", tests);

	//////
	// write command
	//////

	// common commands, in lower case
	cmd += " " + simulationTask + " " + gFile + " " + projectName
	       + " -f " + n2s(boundary) + " -n " + n2s(neutral_cutoff)
	       + " -g " + strmoi;

	if (simulationTask == "1" || simulationTask == "7") {
		cmd += " -q " + n2s(propFunctionalRv[0]) + " -p " + n2s(propFunctionalRv[1]);
		cmd += " -B " + n2s(oddsRatios[1]) + " -D " + n2s(oddsRatios[3]);
		if (!fEqual(oddsRatios[0], 0.0)) cmd += " -A " + n2s(oddsRatios[0]);
		if (!fEqual(oddsRatios[2], 0.0)) cmd += " -C " + n2s(oddsRatios[2]);
		cmd += " -F " + n2s(baselinef);
		cmd += " -X " + n2s(nCases) + " -Y " + n2s(nCtrls);
		if (simulationTask == "1") cmd += " -Z " + n2s(nUnphenotyped);
		if (shouldUseGenPool) cmd += " -U";
	}
	if (simulationTask == "2") {
		cmd += " -q " + n2s(propFunctionalRv[0]) + " -p " + n2s(propFunctionalRv[1]);
		cmd += " -B " + n2s(oddsRatios[1]) + " -D " + n2s(oddsRatios[3]);
		if (!fEqual(oddsRatios[0], 0.0)) cmd += " -A " + n2s(oddsRatios[0]);
		if (!fEqual(oddsRatios[2], 0.0)) cmd += " -C " + n2s(oddsRatios[2]);
		cmd += " -F " + n2s(baselinef);
		cmd += " -W " + n2s(nPopulation);
		if (shouldUseGenPool) cmd += " -U";
	}
	if (simulationTask == "3") {
		cmd += " -q " + n2s(propFunctionalRv[0]) + " -p " + n2s(propFunctionalRv[1]);
		cmd += " -G " + n2s(pars[0]) + " -H " + n2s(pars[1]);
		if (isParVariable) {
			cmd += " -I";
		}
		cmd += " -X " + n2s(nCases) + " -Y " + n2s(nCtrls) + " -Z " + n2s(nUnphenotyped);
	}
	if (simulationTask == "8") {
		cmd += " -q " + n2s(propFunctionalRv[0]) + " -p " + n2s(propFunctionalRv[1]);
		cmd += " -G " + n2s(pars[0]) + " -H " + n2s(pars[1]);
		if (isParVariable) {
			cmd += " -I";
		}
		cmd += " -F " + n2s(baselinef);
		cmd += " -X " + n2s(nCases) + " -Y " + n2s(nCtrls) + " -Z " + n2s(nUnphenotyped);
	}
	if (simulationTask == "4") {
		cmd += " -q " + n2s(propFunctionalRv[0]) + " -p " + n2s(propFunctionalRv[1]);
		if (!fEqual(qtcoefs[0], 0.0)) cmd += " -J " + n2s(qtcoefs[0]);
		cmd += " -K " + n2s(qtcoefs[1]) + " -L " + n2s(qtcoefs[2]);
		cmd += " -W " + n2s(nPopulation);
		if (shouldUseGenPool) cmd += " -U";
	}
	if (simulationTask == "5") {
		cmd += " -q " + n2s(propFunctionalRv[0]) + " -p " + n2s(propFunctionalRv[1]);
		if (!fEqual(qtcoefs[0], 0.0)) cmd += " -J " + n2s(qtcoefs[0]);
		cmd += " -K " + n2s(qtcoefs[1]) + " -L " + n2s(qtcoefs[2]);
		cmd += " -M " + n2s(qtcuts[0]) + " -N " + n2s(qtcuts[1]);
		if (shouldMarkBin) cmd += " -O";
		if (nPopulation > 0) cmd += " -W " + n2s(nPopulation);
		else cmd += " -X " + n2s(nCases) + " -Y " + n2s(nCtrls);
		if (shouldUseGenPool) cmd += " -U";
	}
	if (simulationTask == "6") {
		cmd += " -P " + n2s(percentageCausal) + " -Q " + n2s(propHeterCases);
		if (isMendelAlleleFixed) cmd += " -R";
		cmd += " -X " + n2s(nCases) + " -Y " + n2s(nCtrls);
	}
	if (!fEqual(propMissingData[0], 0.0)) cmd += " -a " + n2s(propMissingData[0]);
	if (!fEqual(propMissingData[1], 0.0)) cmd += " -b " + n2s(propMissingData[1]);
	if (!fEqual(propMissingData[2], 0.0)) cmd += " -c " + n2s(propMissingData[2]);
	if (!fEqual(propMissingData[3], 0.0)) cmd += " -d " + n2s(propMissingData[3]);
	if (!fEqual(missingLowMaf, 0.0)) cmd += " -e " + n2s(missingLowMaf);
	if (shouldMarkMissing) cmd += " -k";
	if (isSynoKept) cmd += " -i";
	if (isCvTrimmed) cmd += " -j";
	if (isPedWritten) {
		cmd += " -z";
		return cmd;
	}
	cmd += " -l " + n2s(mafLower) + " -m " + n2s(mafUpper);
	cmd += " -t " + n2s("PLACEHOLDER");
	cmd += " -s " + n2s(alpha);

	if (need_perm) {
		cmd += " -u " + n2s(nPermutations);
	}

	cmd += " -r " + n2s(nReplicates);

	if (seed != 10086) cmd += " -y " + n2s(seed);

	if (verbose && !quiet) cmd += " -v";
	if (quiet) cmd += " -x";

	return cmd;
}


const char * args_dsc(const std::string name, bool empty)
{

	std::string res;
	std::string formattedtestlist = "CMC, CMC-one, CMC-one-midP, RVE, RVE-one, RVE-one-midP, CMCST, CMCST-one \n\t| WSS, WSS-one, WSSPM, WSSPM-one, MZ, MZ-one, KBAC, KBAC-one, KBACST, KBACST-one, VT, VT-one, VTfisher, VTfisher-one \n\t| aSum, RBT, RBT-one, calpha, calpha-one, RareCover, RareCover-one, WF, WF-one, WF-midP, WF-one-midP, SKAT\n\t| CMCQT, CMCQT-one, MZQT, MZQT-one, MZQTPM, MZQTPM-one, ExtremeQT, ExtremeQT-one";

	if (empty) {
		res = "";
	}else if (name == "task") {
		res = "\n\tAnalysis task. \n\t| Type values \"1~7\"\n\t| 1) Case-ctrl samples given odds ratio and prevalence \n\t| 2) Population samples given odds ratio and prevalence \n\t| 3) Case-ctrl samples given population attributable risk \n\t| 4) Quantitative traits samples \n\t| 5) Extreme quantitative traits samples \n\t| 6) Mendelian traits samples \n\t| 7) Affected/unaffected sib-pairs \n\t";
	}else if (name == "gdata") {
		res = "\n\tGenetic data files for the simulation to be based on. \n\t| STRING\n\t| Proper gdata.maf, gdata.ann and gdata.pos files need to be provided to the program (gdata.hap file will be needed if --use_haplotype_pool is envoked). \n\t";
	}else if (name == "pname") {
		res = "\n\tProject name. \n\t| STRING \n\t| set output file names. \n\t| with \"-v\" option, the program will generate summary files of power calculation \n\t| with \"-z\" option, the program will generate simulated data only\n\t";
	}else if (name == "f") {
		res = "\n\tDefinition of rare variants. \n\t| FLOAT\n\t| variants having MAF <= frequency will be considered as \"rare\" variants\n\t";
	}else if (name == "q") {
		res = "\n\tProportion of FUNCTIONAL deleterious variants. \n\t| FLOAT\n\t| 0.0 <= fraction <= 1.0\n\t| (1-proportion)x100\% is the proportion of non-causal deleterious variants (noise) \n\t| This is NOT the proportion of deleterious variants, which should have been defined in <gdata.ann> file\n\t";
	}else if (name == "p") {
		res = "\n\tProportion of FUNCTIONAL protective variants. \n\t| FLOAT\n\t| 0.0 <= fraction <= 1.0\n\t| (1-proportion)x100\% is the proportion of non-causal protective variants (noise) \n\t| This is NOT the proportion of protective variants, which should have been defined in <gdata.ann> file\n\t";
	}else if (name == "g") {
		res = "\n\tMode of inheritance under which the phenotype data is simulated/analyzed.\n\t| STRING (size = 2)\n\t| \"A\" (additive), \"D\" (dominant), \"R\" (recessive), \"M\" (multiplicative), \n\t| \"C\" (compound dominant for non-mendelian traits, or compound recessive for mendelian traits)\n\t| Input should be either one character like \"X\" (for both simulation and analysis) or two characters like \"XY\" (X for simulation, Y for analysis)\n\t";
	}else if (name == "A") {
		res = "\n\t[task=1,2] Minimum odds ratio for deleterious variants. \n\t| FLOAT\n\t| 0.0 for fixed effect size model; >=1.0 for variable effect sizes model\n\t";
	}else if (name == "B") {
		res = "\n\t[task=1,2] Maximum odds ratio for deleterious variants. \n\t| FLOAT\n\t| >=1.0 AND > $OR_deleterious_min\n\t| will be the odds ratio for deleterious variants in fixed effect size model, or the maximum odds ratio in variable effect sizes model\n\t";
	}else if (name == "C") {
		res = "\n\t[task=1,2] Minimum odds ratio for protective variants. \n\t| FLOAT\n\t| 0.0 for fixed effect size model; <1.0 for variable effect sizes model\n\t";
	}else if (name == "D") {
		res = "\n\t[task=1,2] Maximum odds ratio for protective variants. \n\t| FLOAT\n\t| <=1.0 AND > $OR_protective_min\n\t| will be the odds ratio for protective variants in fixed effect size model, or the maximum odds ratio in variable effect sizes model\n\t";
	}else if (name == "E") {
		res = "\n\t[task=1,2] Odds ratio for common variants. \n\t| FLOAT\n\t| =1.0 for neutral, >1.0 for deleterious, <1.0 for protective \n\t| this is the odds ratio for all variants having MAF > $define_rare, i.e., the \"common\" variants \n\t| assuming all common variants have fixed effect size\n\t";
	}else if (name == "F") {
		res = "\n\t[task=1,2] Disease prevalence. \n\t| FLOAT \n\t| will be used as baseline penetrance of a gene (baseline penetrance ~= disease prevalence)\n\t";
	}else if (name == "G") {
		res = "\n\t[task=3] Total population attributable risk for deleterious variants. \n\t|FLOAT \n\t| 0.0 <= fraction <= 1.0\n\t";
	}else if (name == "H") {
		res = "\n\t[task=3] Total population attributable risk for protective variants. \n\t|FLOAT \n\t| 0.0 <= fraction <= 1.0\n\t";
	}else if (name == "I") {
		res = "\n\t[task=3] Locus specific population attributable risk is inversely proportional to its MAF (rather than uniformly distributed).\n\t| BOOLEAN \n\t";
	}else if (name == "J") {
		res = "\n\t[task=4,5] Minimum mean value shift per variant. \n\t| FLOAT\n\t| 0.0 for fixed effect size model; >0.0 for variable effect sizes model\n\t";
	}else if (name == "K") {
		res = "\n\t[task=4,5] Maximum mean value shift per variant. \n\t| FLOAT\n\t| >=0.0 AND > $QT_effect_min \n\t| will be locus effect to quantitative trait for fixed effect size model, or maximum effect for variable effect sizes model\n\t| the mean of the quantitative trait will be shifted by (multiplier*sigma) where sigma is standard deviation of the quantitative trait\n\t| will automatically shift trait values up for deleterious variants, down for protective variants\n\t";
	}else if (name == "L") {
		res = "\n\t[task=4,5] Mean value shift for common variants. \n\t| FLOAT\n\t| >=0.0 \n\t| this is the effect for all variants having MAF > $define_rare, i.e., the \"common\" variants \n\t| assuming all common variants have fixed effect size\n\t";
	}else if (name == "M") {
		res = "\n\t[task=5] Lower percentile cutoff for quantitative traits in extreme QT sampling. \n\t| FLOAT\n\t| 0.0 <= fraction <= $QT_upper_percentile\n\t";
	}else if (name == "N") {
		res = "\n\t[task=5] Upper percentile cutoff for quantitative traits in extreme QT sampling. \n\t| FLOAT\n\t| $QT_lower_percentile <= fraction <= 1.0\n\t";
	}else if (name == "O") {
		res = "\n\t[task=5] Re-code extreme quantitative traits using binary codings.\n\t| BOOLEAN\n\t| if envoked, will convert extreme quantitative traits into binary traits\n\t";
	}else if (name == "P") {
		res = "\n\t[task=6] Percentage of rare variants being causal in Mendelian traits simulation. \n\t| FLOAT\n\t| 0.0 <= fraction <= 1.0\n\t| any of the variants having the top (fraction)x100\% smallest MAF will be contributable to a Mendelian trait \n\t| this assumes the disease is allelic heterogeneous \n\t";
	}else if (name == "Q") {
		res = "\n\t[task=6] Porportion of cases that do not carry the disease allele at the gene region under study. \n\t| FLOAT\n\t| 0.0 <= fraction <= 1.0\n\t";
	}else if (name == "R") {
		res = "\n\t[task=6] No allelic heterogeneity for Mendelian traits. \n\t| BOOLEAN\n\t| if envoked, will fix the causal variants of the Mendelian trait to the one that has the (fraction)x100\%-th smallest MAF \n\t| rather than using all variants having the top (fraction)x100\% smallest MAF\n\t";
	}else if (name == "W") {
		res = "\n\t[task=2,4,5] Number of total samples. \n\t| INT (>=0)\n\t| for extreme QT simulations (task=5) $num_all_samples>0 will have extreme QT samples collected from this finite cohort.\n\t";
	}else if (name == "X") {
		res = "\n\t[task=1,3,5,6] Number of cases. \n\t| INT (>0)\n\t| for extreme QT simulations (task=5) it will be #samples having high QT values from the population when $num_all_samples is set to 0\n\t";
	}else if (name == "Y") {
		res = "\n\t[task=1,3,5,6] Number of ctrls. \n\t| INT (>0)\n\t| for extreme QT simulations (task=5) it will be #samples having low QT values from the population when $num_all_samples is set to 0\n\t";
	}else if (name == "Z") {
		res = "\n\t[task=1,3,6] Number of unphenotyped cohort ctrls. \n\t| INT (>0)\n\t";
	}else if (name == "U") {
		res = "\n\t[task=1,2,4,5] Randomly sample haplotypes from haplotype pool file $gdata.hap, rather than generating haplotypes on the fly. \n\t| BOOLEAN \n\t";
	}else if (name == "a") {
		res = "\n\tProportion of deleterious variants missing data. \n\t| FLOAT \n\t| missing genotypes will be coded as wildtype genotype by default\n\t";
	}else if (name == "b") {
		res = "\n\tProportion of protective variants missing data. \n\t| FLOAT \n\t| missing genotypes will be coded as wildtype genotype by default\n\t";
	}else if (name == "c") {
		res = "\n\tProportion of non-causal variants missing data. \n\t| FLOAT \n\t| missing genotypes will be coded as wildtype genotype by default\n\t";
	}else if (name == "d") {
		res = "\n\tProportion of synonymous variants missing data. \n\t| FLOAT \n\t| missing genotypes will be coded as wildtype genotype by default\n\t";
	}else if (name == "e") {
		res = "\n\tVariants having MAF < $missing_low_maf will be marked as missing data. \n\t| FLOAT \n\t| note that $missing_low_maf is compared against the haplotype pool, not the sample\n\t";
	}else if (name == "k") {
		res = "\n\tRe-code missing data. \n\t| BOOLEAN\n\t| if envoked, will re-code missing data from wildtype genotype to \"-9\", indicating missingness\n\t";
	}else if (name == "i") {
		res = "\n\tKeep synonymous variants from analysis. \n\t| BOOLEAN\n\t";
	}else if (name == "j") {
		res = "\n\tRemove common variant sites from analysis. \n\t| BOOLEAN\n\t| the \"common\" loci refers to variants in the haplotype pool having MAF > $define_rare \n\t";
	}else if (name == "l") {
		res = "\n\tLower bound of observed sample minor allele frequency. \n\t| FLOAT\n\t| loci having observed MAF < $maf_lower will not be analyzed \n\t";
	}else if (name == "m") {
		res = "\n\tUpper bound of observed sample minor allele frequency. \n\t| FLOAT\n\t| loci having observed MAF > $maf_upper will not be analyzed \n\t";
	}else if (name == "n") {
		res = "\n\tAnnotation value cut-off that defines a variant to be \"neutral\" in evolution (either synonymous or non-coding). \n\t| FLOAT\n\t| loci having annotation value (in gdata.ann) between (-$annotation_cutoff, +$annotation_cutoff) will be regarded \"neutral\" and will not contribute to phenotype \n\t";
	}else if (name == "t") {
		res = "\n\tAssociation test method. \n\t| STRING \n\t| " + formattedtestlist + "\n\t";
	}else if (name == "s") {
		res = "\n\tSignificance level at which power will be evaluated. \n\t| FLOAT\n\t";
	}else if (name == "r") {
		res = "\n\tNumber of replicates for power evaluation. \n\t| INT (>0)\n\t";
	}else if (name == "u") {
		res = "\n\tNumber of permutations, only applicable to permutation based methods. \n\t| INT (>0)\n\t";
	}else if (name == "y") {
		res = "\n\tSeed for random number generator. \n\t| INT (>=0) \n\t| =0 is to use a random seed (seed = system time + process ID)\n\t";
	}else if (name == "z") {
		res = "\n\tWrite out simulated genotype-phenotype data to file $pname, rather than calculating power. \n\t| BOOLEAN \n\t| Genotype file - one snv per row, first column: snv id, subsequent columns: genotype values (0/1/2/NA) of each sample. \n\t| Phenotype file - one subject per row, first column: subject id, second column: quantitative/binary phenotypes. \n\t| Mapping file - each row maps a snv to a gene, first column: gene id, second column: snv id, third column: external MAF, fourth column: observed MAF.\n\t";
	}else if (name == "v") {
		res = "\n\t Maximal screen output information as well as file output for intermediate statistic such as p-values, etc.\n\t";
	}else if (name == "x") {
		res = "\n\t Only output result to screen. No file output.\n\t";
	}else {
		std::cerr << "Error input in args_dsc" << std::endl;
		exit(1);
	}

	char * cstr = new char [res.size() + 1];
	strcpy(cstr, res.c_str());
	return cstr;
}



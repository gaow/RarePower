// =====================================================================================
//
//       Filename:  main.h
//
//    Description:  simulations, sampling, and power calculations for gene based association studies.
//
//        Version:  1.0
//        Created:  01/04/2011 10:52:55 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Gao Wang (gw), wangow@gmail.com
//                  Baylor College of Medicine, Texas, USA
//        License:  GNU General Public License <http://www.gnu.org/licenses/>
//                  Copyright (c) 2011, Gao Wang
//
// =====================================================================================

#include "assoctests.h"

std::string check_options(std::string prog_name, std::string & projectName, std::string & gFile, double & boundary, double & neutral_cutoff,
	vectorF & propFunctionalRv, char & moi, std::string & simulationTask, vectorF & oddsRatios, double & baselinef, vectorF & pars,
	bool & isParConst, vectorF & qtcoefs, vectorF & qtcuts, bool & shouldMarkBin, double & percentageCausal,
	bool & isAllelicHeterogeneous, vectorF & propMissingData, double & missingLowMaf, bool & shouldMarkMissing, unsigned & nCases,
	unsigned & nCtrls, unsigned & nPopulation, unsigned & nUnphenotyped, double & propHeterCases, bool & isSynoTrimmed,
	bool & isCvTrimmed, bool & isPedWritten, std::string & test, double & mafLower,
	double & mafUpper, double & alpha, unsigned & nPermutations, unsigned & nReplicates, bool & verbose, bool & quiet, unsigned & seed,
	bool & shouldUseGenPool);

gwBaseTest * testFactory(std::string const & classname)
{
	if (classname == "CMC" || classname == "CMC-one") {
		return new CmcfisherP();
	} else if (classname == "WSS" || classname == "WSS-one") {
		return new WssRankPA();
	} else if (classname == "RVE" || classname == "RVE-one") {
		return new RvefisherP();
	} else if (classname == "CMCST" || classname == "CMCST-one") {
		return new CmcstP();
	} else if (classname == "MZ" || classname == "MZ-one") {
		return new AnrvstP();
	} else if (classname == "CMCPM") {
		return new CmcchiP();
	} else if (classname == "WSSPM" || classname == "WSSPM-one") {
		return new WssRankP();
	} else if (classname == "KBAC-one" || classname == "KBAC") {
		return new KbacP();
	} else if (classname == "KBACST-one" || classname == "KBACST") {
		return new KbacstP();
	} else if (classname == "VT" || classname == "VT-one") {
		return new VtP();
	} else if (classname == "VTfisher" || classname == "VTfisher-one") {
		return new VtFisherP();
	} else if (classname == "aSum") {
		return new AsumP();
	} else if (classname == "CMCQT" || classname == "CMCQT-one") {
		return new CmcqtP();
	} else if (classname == "MZQT" || classname == "MZQT-one") {
		return new AnrvqtP();
	} else if (classname == "MZQTPM" || classname == "MZQTPM-one") {
		return new AnrvqtPermP();
	} else if (classname == "ExtremeQT" || classname == "ExtremeQT-one") {
		return new AnrvqtCondP();
	} else if (classname == "RBT" || classname == "RBT-one") {
		return new TestRareP();
	} else if (classname == "calpha" || classname == "calpha-one") {
		return new CalphaP();
	} else if (classname == "RareCover" || classname == "RareCover-one") {
		return new RareCoverP();
	} else if (classname.find("WF") < classname.size()) {
		return new WsFisherP();
	} else if (classname == "SKAT") {
		return new SkatP();
	} else {
		return NULL;
	}
}



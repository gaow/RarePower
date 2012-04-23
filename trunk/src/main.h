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

gpow::gwBaseTest * testFactory(std::string const & classname)
{
	if (classname == "CMC" || classname == "CMC-one") {
		return new gpow::CmcfisherP();
	} else if (classname == "WSS" || classname == "WSS-one") {
		return new gpow::WssRankPA();
	} else if (classname == "RVE" || classname == "RVE-one") {
		return new gpow::RvefisherP();
	} else if (classname == "CMCST" || classname == "CMCST-one") {
		return new gpow::CmcstP();
	} else if (classname == "MZ" || classname == "MZ-one") {
		return new gpow::AnrvstP();
	} else if (classname == "CMCPM") {
		return new gpow::CmcchiP();
	} else if (classname == "WSSPM" || classname == "WSSPM-one") {
		return new gpow::WssRankP();
	} else if (classname == "KBAC-one" || classname == "KBAC") {
		return new gpow::KbacP();
	} else if (classname == "KBACST-one" || classname == "KBACST") {
		return new gpow::KbacstP();
	} else if (classname == "VT" || classname == "VT-one") {
		return new gpow::VtP();
	} else if (classname == "VTfisher" || classname == "VTfisher-one") {
		return new gpow::VtFisherP();
	} else if (classname == "aSum") {
		return new gpow::AsumP();
	} else if (classname == "CMCQT" || classname == "CMCQT-one") {
		return new gpow::CmcqtP();
	} else if (classname == "MZQT" || classname == "MZQT-one") {
		return new gpow::AnrvqtP();
	} else if (classname == "MZQTPM" || classname == "MZQTPM-one") {
		return new gpow::AnrvqtPermP();
	} else if (classname == "ExtremeQT" || classname == "ExtremeQT-one") {
		return new gpow::AnrvqtCondP();
	} else if (classname == "RBT" || classname == "RBT-one") {
		return new gpow::TestRareP();
	} else if (classname == "calpha" || classname == "calpha-one") {
		return new gpow::CalphaP();
	} else if (classname == "RareCover" || classname == "RareCover-one") {
		return new gpow::RareCoverP();
	} else if (classname.find("WF") < classname.size()) {
		return new gpow::WsFisherP();
	} else if (classname == "SKAT") {
		return new gpow::SkatP();
	} else {
		return NULL;
	}
}

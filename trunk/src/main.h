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


//! \fn progress_bar()
#include <sstream>
template <typename T>
std::string n2s(T Number)
{
	std::stringstream ss;

	ss << Number;
	return ss.str();
}


void progress_bar(unsigned int x, unsigned int N);

std::string check_options(std::string prog_name, std::string & projectName, std::string & gFile, double & boundary, double & neutral_cutoff,
	vectorF & propFunctionalRv, char & moi, std::string & simulationTask, vectorF & oddsRatios, double & baselinef, vectorF & pars,
	bool & isParConst, vectorF & qtcoefs, vectorF & qtcuts, bool & shouldMarkBin, double & percentageCausal,
	bool & isAllelicHeterogeneous, vectorF & propMissingData, double & missingLowMaf, bool & shouldMarkMissing, unsigned & nCases,
	unsigned & nCtrls, unsigned & nPopulation, unsigned & nUnphenotyped, double & propHeterCases, bool & isSynoTrimmed,
	bool & isCvTrimmed, bool & isPedWritten, std::string & test, double & mafLower,
	double & mafUpper, double & alpha, unsigned & nPermutations, unsigned & nReplicates, bool & verbose, bool & quiet, unsigned & seed,
	bool & shouldUseGenPool);


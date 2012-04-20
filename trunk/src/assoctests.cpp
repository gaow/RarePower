// =====================================================================================
//
//       Filename:  assoctests.cpp
//
//    Description:  Statistical tests for association mapping
//
//        Version:  1.0
//        Created:  01/04/2011 10:48:03 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Gao Wang (gw), wangow@gmail.com
//                  Baylor College of Medicine, Texas, USA
//        License:  GNU General Public License <http://www.gnu.org/licenses/>
//                  Copyright (c) 2011, Gao Wang
//
// =====================================================================================


#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <list>
#include <numeric>
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_vector_double.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_errno.h"
#include "gw_maths.h"
#include "gw_utilities.h"
#include "assoctests.h"
#define SKAT_GSL

namespace {
const double AFFECTED = 2.0, UNAFFECTED = 1.0, HOMO_ALLELE = 2.0,
             MINOR_ALLELE = 1.0, MAJOR_ALLELE = 0.0, MISSING_ALLELE = -9.0;
}

bool fsame(double a, double b)
{
	return fabs(a - b) < 1.0e-6;
}


gwAssociations::gwAssociations(const vectorF & observedMafs, const vectorF & ydat,
	const vector2F & xdat, double mafLower, double mafUpper, double alpha)
{
	assert(xdat.size() != 0 && ydat.size() == xdat.size()
		&& observedMafs.size() == xdat[0].size() && mafLower >= 0.0
		&& mafUpper <= 1.0 && mafUpper > mafLower && alpha > 0.0 && alpha < 1.0);

	__xdat = xdat;
	__ydat = ydat;
	__observedMafs = observedMafs;
	__mafLower = mafLower;
	__mafUpper = mafUpper;
	__alpha = alpha;

	__isDebug = false;

	if (__isDebug) {
		std::cout << __observedMafs << std::endl;
		std::cout << __ydat << std::endl;
	}
}


gwAssociations::~gwAssociations()
{
}


void gwAssociations::setVerbose(int debuglevel)
{
	if (debuglevel != 0)
		__isDebug = true;
	else __isDebug = false;
	return;
}


void gwAssociations::debug(int showWhat) const
{
	m_printPunches(60);
	std::cout.precision(9);

	switch (showWhat) {
	case 1:
		std::cout << "__observedMafs:\n" << __observedMafs << "\n" << std::endl;
		break;
	case 2:
		std::cout << "__ydat:\n" << __ydat << "\n" << std::endl;
		break;
	case 3:
	{
		std::cout << "__xdat[0]:\n" << __xdat[0] << "\n" << std::endl;
		std::cout << "__xdat[__xdat.size()]:\n" << __xdat[__xdat.size()] << "\n" << std::endl;
	}
	break;
	default:
	{
		std::cout << "__observedMafs:\n" << __observedMafs << "\n" << std::endl;
		std::cout << "__ydat:\n" << __ydat << "\n" << std::endl;
		std::cout << "__xdat[0]:\n" << __xdat[0] << "\n" << std::endl;
		std::cout << "__xdat[__xdat.size() - 1]:\n" << __xdat[__xdat.size() - 1] << "\n" << std::endl;
	}
	break;
	}

	m_printPunches(60);
	return;
}


//!\brief Collapsing variants simple LR score test statistic

double gwAssociations::calcCmcstP(unsigned sided, unsigned nPermutations, unsigned adaptive)
{
	/*! Li&Leal 2008 collapsing. Statistical test is the score test for simple logistic regression model.
	 */

	//!- trim data by mafs upper-lower bounds
	m_trimXdat();

	//!- CMC scoring
	vectorF regressors = m_indicateRegionalVariants();

	if (__isDebug)
		std::cout << regressors << std::endl;

	//!- logistic regression score statistic + permutation
	double xbar = gw_mean(regressors);
	unsigned nCases = 0;
	for (unsigned i = 0; i != __ydat.size(); ++i)
		if (__ydat[i] == AFFECTED)
			++nCases;

	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= nPermutations) {
		double statistic = m_calcSimpleLogitRegScore(regressors, __ydat, xbar, nCases, sided);
		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (adaptive != 0)
				pvalue = m_checkAdaptivePvalue(permcount1, permcount2, iPermutation, adaptive, sided);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		random_shuffle(__ydat.begin(), __ydat.end());
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		if (sided == 1 || sided == 0)
			pvalue = (1.0 + permcount1) / (1.0 + nPermutations);
		else {
			double permcount = gw_dmin(permcount1, permcount2);
			pvalue = (2.0 * permcount + 2.0) / (1.0 + nPermutations);
		}
	}
	return pvalue;
}


//!\brief Counts of variants simple LR score test statistic

double gwAssociations::calcAnrvstP(unsigned sided, unsigned nPermutations, unsigned adaptive)
{
	/*! Andrew P. Morris 2009 collapsing method. Statistical test is the score test for simple logistic regression model.
	 */

	//!- trim data by mafs upper-lower bounds
	m_trimXdat();

	//!- ANRV scoring
	vectorF regressors = m_countRegionalVariants();
	if (__isDebug)
		std::cout << regressors << std::endl;
	//!- logistic regression score statistic + permutation
	double xbar = gw_mean(regressors);
	unsigned nCases = 0;
	for (unsigned i = 0; i != __ydat.size(); ++i)
		if (__ydat[i] == AFFECTED)
			++nCases;


	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= nPermutations) {
		double statistic = m_calcSimpleLogitRegScore(regressors, __ydat, xbar, nCases, sided);
		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (adaptive != 0)
				pvalue = m_checkAdaptivePvalue(permcount1, permcount2, iPermutation, adaptive, sided);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		random_shuffle(__ydat.begin(), __ydat.end());
		++iPermutation;
	}

	if (__isDebug)
		std::cout << permcount1 << " " << permcount2 << std::endl;

	if (pvalue <= 1.0) ;
	else {
		if (sided == 1 || sided == 0)
			pvalue = (1.0 + permcount1) / (1.0 + nPermutations);
		else {
			double permcount = gw_dmin(permcount1, permcount2);
			pvalue = (2.0 * permcount + 2.0) / (1.0 + nPermutations);
		}
	}
	return pvalue;
}


double gwAssociations::calcCmcchiP(unsigned sided, unsigned nPermutations, unsigned adaptive)
{
	/*!* By Li&Leal 2008 collapsing <br>
	 * * The original CMC test combine the collapsed rare variants and uncollapsed common variants in a Hotelling's T test. <br>
	 * * Here assume common variants are not analyzed -- the CMC test is now the "Collapsing test" <br>
	 * * Common variants not properly handled here but may be dealt with elsewhere (i.e. regression framework, etc) <br> <br>
	 * Implementation:
	 */

	//!- trim data by mafs upper-lower bounds
	m_trimXdat();

	//!- CMC scoring
	vectorF regressors = m_indicateRegionalVariants();
	if (__isDebug)
		std::cout << regressors << std::endl;

	//!- 2 by 2 Chisq test

	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= nPermutations) {
		double statistic = m_calc2X2Chisq(regressors, __ydat);
		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (adaptive != 0)
				pvalue = m_checkAdaptivePvalue(permcount1, permcount2, iPermutation, adaptive, sided);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		random_shuffle(__ydat.begin(), __ydat.end());
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		if (sided == 1 || sided == 0)
			pvalue = (1.0 + permcount1) / (1.0 + nPermutations);
		else {
			double permcount = gw_dmin(permcount1, permcount2);
			pvalue = (2.0 * permcount + 2.0) / (1.0 + nPermutations);
		}
	}
	return pvalue;
}


double gwAssociations::calcCmcfisherP(unsigned sided)
{
	/*!* Similar with the chi-squared test based CMC but uses the Fisher's exact test, no permutation <br> <br>
	 * Implementation:
	 */

	//!- trim data by mafs upper-lower bounds
	m_trimXdat();

	//!- CMC scoring
	vectorF regressors = m_indicateRegionalVariants();
	if (__isDebug)
		std::cout << regressors << std::endl;

	//!- 2 by 2 Fisher's exact test, one or two sided
	double pvalue = m_calc2X2Fisher(regressors, __ydat, sided);

	return pvalue;
}


double gwAssociations::calcRvefisherP(unsigned sided)
{
	/*!* Carrier frequencies are compared for variants that are only present in cases to those that are only observed controls <br> <br>
	 * Implementation:
	 */

	//!- trim data by mafs upper-lower bounds
	m_trimXdat();

	//!- RVE scoring
	vectorF regressors = m_indicateRegionalUniqueVariants();
	if (__isDebug)
		std::cout << regressors << std::endl;

	//!- 2 by 2 Fisher's exact test, one or two sided
	double pvalue = m_calc2X2Fisher(regressors, __ydat, sided);
	return pvalue;
}


double gwAssociations::calcWssRankP(const char moi, unsigned nCtrls, unsigned sided, unsigned nPermutations, unsigned adaptive)
{
	/*! * Implement Browning 2009 Wss paper <br><br>
	 * Implementation:
	 */

	//!- trim data by mafs upper-lower bounds
	m_trimXdat();

	unsigned nCases = __ydat.size() - nCtrls;


	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= nPermutations) {

		//!- Calc MAF in controls
		vectorF weights(0);
		for (unsigned j = 0; j != __xdat[0].size(); ++j) {

			unsigned nVariants = 0;
			//! - Count number of mutations in controls
			for (unsigned i = 0; i != __xdat.size(); ++i) {
				// Control only
				if (__ydat[i] == UNAFFECTED) {
					if (__xdat[i][j] != MISSING_ALLELE)
						nVariants += (unsigned)__xdat[i][j];
					else {
						std::cerr << "Input data problem in gwAssociations::calcWssRankP(). Now Quit." << std::endl;
						exit(-1);
					}
				}else
					continue;
			}
			//! - Compute the "q" for the locus
			weights.push_back((nVariants + 1.0) / (2.0 * nCtrls + 2.0));
		}

		//!- Calc the Weights for each locus
		for (unsigned j = 0; j != __xdat[0].size(); ++j)
			weights[j] = sqrt(2.0 * __xdat.size() * weights[j] * (1.0 - weights[j]));

		vectorF scores(0);

		for (unsigned i = 0; i != __xdat.size(); ++i) {
			double score = 0.0;
			//! - Define genetic score of the individual

			for (unsigned j = 0; j != __xdat[i].size(); ++j) {
				// scan all loci of an individual
				double tmp = 0.0;
				// Dominant model
				if (moi == 'D') {
					if (__xdat[i][j] == HOMO_ALLELE)
						tmp = 1.0;
				}
				// recessive model
				else if (moi == 'R') {
					if (__xdat[i][j] != MAJOR_ALLELE)
						tmp = 1.0;
				}
				// additive and multiplicative models
				else ;

				//! - Parameterize mode of inheritance I_ij = {0, 1, 2} in Browning paper
				tmp = __xdat[i][j] - tmp;
				score += tmp / weights[j];
			}

			//!- Calculate and store the genetic score
			scores.push_back(score);
		}

		if (__isDebug)
			std::cout << scores << std::endl;

		//! - Compute rank test statistic, via <b> Mann_Whitneyu() </b>
		double caseScores[nCases], ctrlScores[nCtrls];
		int tmpa = 0, tmpu = 0;
		for (unsigned i = 0; i != __ydat.size(); ++i) {
			if (__ydat[i] == AFFECTED) {
				caseScores[tmpa] = scores[i];
				++tmpa;
			}else {
				ctrlScores[tmpu] = scores[i];
				++tmpu;
			}
		}

		double statistic = Mann_Whitneyu(caseScores, nCases, ctrlScores, nCtrls);

		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (adaptive != 0)
				pvalue = m_checkAdaptivePvalue(permcount1, permcount2, iPermutation, adaptive, sided);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		random_shuffle(__ydat.begin(), __ydat.end());
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		if (sided == 1 || sided == 0)
			pvalue = (1.0 + permcount1) / (1.0 + nPermutations);
		else {
			double permcount = gw_dmin(permcount1, permcount2);
			pvalue = (2.0 * permcount + 2.0) / (1.0 + nPermutations);
		}
	}
	return pvalue;
}


double gwAssociations::calcWssRankP(const char moi, unsigned nCtrls, unsigned sided)
{
	/*! * Implement Browning 2009 Wss paper normal approximation <br><br>
	 * Implementation:
	 */

	//!- trim data by mafs upper-lower bounds
	m_trimXdat();

	int nCases = __ydat.size() - nCtrls;

	unsigned iPermutation = 0;
	double observedStatistic = 0.0;
	vectorF statistics(0);

	while (iPermutation <= 1000) {

		//!- Calc MAF in controls
		vectorF weights(0);
		for (unsigned j = 0; j != __xdat[0].size(); ++j) {

			int nVariants = 0;
			//! - Count number of mutations in controls
			for (unsigned i = 0; i != __xdat.size(); ++i) {
				// Control only; consider additive codings only
				if (__ydat[i] == UNAFFECTED) {
					if (__xdat[i][j] != MISSING_ALLELE)
						nVariants += (int)__xdat[i][j];
					else {
						std::cerr << "Input data problem in gwAssociations::calcWssRankP. Now Quit." << std::endl;
						exit(-1);
					}
				}else
					continue;
			}
			//! - Compute the "q" for the locus
			weights.push_back((nVariants + 1.0) / (2.0 * nCtrls + 2.0));
		}

		//!- Calc the Weights for each locus
		for (unsigned j = 0; j != __xdat[0].size(); ++j)
			weights[j] = sqrt(2.0 * __xdat.size() * weights[j] * (1.0 - weights[j]));

		vectorF scores(0);

		for (unsigned i = 0; i != __xdat.size(); ++i) {
			double score = 0.0;
			//! - Define genetic score of the individual

			for (unsigned j = 0; j != __xdat[i].size(); ++j) {
				// scan all loci of an individual
				double tmp = 0.0;
				// Dominant model
				if (moi == 'D') {
					if (__xdat[i][j] == HOMO_ALLELE)
						tmp = 1.0;
				}
				// recessive model
				else if (moi == 'R') {
					if (__xdat[i][j] != MAJOR_ALLELE)
						tmp = 1.0;
				}
				// additive and multiplicative models
				else ;

				//! - Parameterize mode of inheritance I_ij = {0, 1, 2} in Browning paper
				tmp = __xdat[i][j] - tmp;
				score += tmp / weights[j];
			}

			//!- Calculate and store the genetic score

			scores.push_back(score);
		}

		if (__isDebug)
			std::cout << scores << std::endl;

		//! - Compute rank test statistic, via <b> Mann_Whitneyu() </b>
		double caseScores[nCases], ctrlScores[nCtrls];
		int tmpa = 0, tmpu = 0;
		for (unsigned i = 0; i != __ydat.size(); ++i) {
			if (__ydat[i] == AFFECTED) {
				caseScores[tmpa] = scores[i];
				++tmpa;
			}else {
				ctrlScores[tmpu] = scores[i];
				++tmpu;
			}
		}

		double statistic = Mann_Whitneyu(caseScores, nCases, ctrlScores, nCtrls);

		if (iPermutation == 0)
			observedStatistic = statistic;
		else
			statistics.push_back(statistic);

		//!- Permutation
		random_shuffle(__ydat.begin(), __ydat.end());
		++iPermutation;
	}

	double mean = gw_mean(statistics);
	double sd = gw_sd(statistics);
	if (sd == 0.0)
		sd = 1.0e-6;

	double statisticstd = (observedStatistic - mean) / sd;
	double pvalue = 9.0;
	if (sided == 1) {
		pvalue = gsl_cdf_ugaussian_Q(statisticstd);
	}else if (sided == 0) {
		statisticstd = statisticstd * statisticstd;
		pvalue = gsl_cdf_chisq_Q(statisticstd, 1.0);
	}else {
		std::cerr << "ERROR: \"sided\" should be either 0 (default two-sided test) or 1 (one-sided test)" << std::endl;
		exit(-1);
	}
	return pvalue;
}


double gwAssociations::calcKbacP(unsigned sided, unsigned nPermutations, unsigned adaptive)
{
	/*! * the KBAC Statistic: sum of genotype pattern frequencies differences, weighted by hypergeometric kernel. <br>
	 * * It is a permutation based one/two-sided test.  <br>
	 * * See <em> Liu DJ 2010 PLoS Genet. </em> <br><br>
	 * * Implementation:
	 */
	if (sided == 0 || sided > 2)
		sided = 1;

	//!- trim data by mafs upper-lower bounds
	m_trimXdat();


	unsigned sampleSize = __ydat.size();
	// sample size
	unsigned regionLen = __xdat[0].size();
	// candidate region length
	unsigned nCases = 0;
	// case size

	for (unsigned i = 0; i != __ydat.size(); ++i) {
		if (__ydat[i] == AFFECTED)
			++nCases;
	}
	unsigned nCtrls = sampleSize - nCases;

	vectorF genotypeId(sampleSize);
	//!-Compute unique genotype patterns (string) as ID scores (double)

	for (unsigned i = 0; i != sampleSize; ++i) {

		double vntIdL = 0.0;
		double vntIdR = 0.0;
		const double ixiix = pow(9.0, 10.0);
		unsigned lastCnt = 0;
		unsigned tmpCnt = 0;

		for (unsigned j = 0; j != regionLen; ++j) {

			if (__xdat[i][j] != MISSING_ALLELE && __xdat[i][j] != MAJOR_ALLELE)
				vntIdR += pow(3.0, 1.0 * (j - lastCnt)) * __xdat[i][j];
			else
				continue;
			if (vntIdR >= ixiix) {
				vntIdL = vntIdL + 1.0;
				vntIdR = vntIdR - ixiix;
				lastCnt = lastCnt + tmpCnt + 1;
				tmpCnt = 0;
				continue;
			}else {
				++tmpCnt;
				continue;
			}
		}

		genotypeId[i] = vntIdL + vntIdR * 1e-10;
		// one-to-one "ID number" for a genotype pattern
	}

	// unique genotype patterns
	// use "list" and some STL algorithms

	std::list<double> uniqueId(genotypeId.begin(), genotypeId.end());
	uniqueId.remove(MAJOR_ALLELE);

	if (uniqueId.size() == 0)
		return 1.0;

	uniqueId.sort();
	uniqueId.unique();
	// remove wildtype and get unique genotype patterns
	//!- Unique genotype patterns that occur in the sample


	vectorF uniquePattern(uniqueId.size());
	copy(uniqueId.begin(), uniqueId.end(), uniquePattern.begin());
	uniqueId.clear();

	//  std::cout << uniquePattern << std::endl;

	// count number of sample individuals for each genotype pattern
	unsigned uniquePatternCounts[uniquePattern.size()];
	for (unsigned u = 0; u != uniquePattern.size(); ++u)
		uniquePatternCounts[u] = 0;

	for (unsigned i = 0; i != sampleSize; ++i) {
		// for each sample, identify/count its genotype pattern

		for (unsigned u = 0; u != uniquePattern.size(); ++u) {

			if (genotypeId[i] == uniquePattern[u]) {
				// genotype pattern identified
				++uniquePatternCounts[u];
				// count this genotype pattern
				break;
			}else ;
			// genotype pattern not found -- move on to next pattern
		}
	}


	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;
	while (iPermutation <= nPermutations) {

		// the KBAC statistic. Will be of length 1 or 2
		vectorF kbacStatistics(0);
		// two models
		for (unsigned s = 0; s != sided; ++s) {

			//!- count number of sample cases (for the 1st model, or ctrls for the 2nd model) for each genotype pattern
			unsigned uniquePatternCountsSub[uniquePattern.size()];
			for (unsigned u = 0; u != uniquePattern.size(); ++u)
				uniquePatternCountsSub[u] = 0;
			// genotype pattern counts in cases (for the 1st model, or ctrls for the 2nd model)

			for (unsigned i = 0; i != sampleSize; ++i) {
				if (__ydat[i] == (AFFECTED - 1.0 * s) ) {
					// for each "case (for the 1st model, or ctrls for 2nd model)", identify/count its genotype pattern
					for (unsigned u = 0; u != uniquePattern.size(); ++u) {
						if (genotypeId[i] == uniquePattern[u]) {
							// genotype pattern identified in cases (for the 1st model, or ctrls for 2nd model)
							++uniquePatternCountsSub[u];
							// count this genotype pattern
							break;
						}else ;
						// genotype pattern not found -- move on to next pattern
					}
				}else ;
			}

			//!- KBAC weights
			double uniquePatternWeights[uniquePattern.size()];
			// genotype pattern weights, the hypergeometric distribution cmf
			for (unsigned u = 0; u != uniquePattern.size(); ++u)
				uniquePatternWeights[u] = 0.0;

			for (unsigned u = 0; u != uniquePattern.size(); ++u) {
				if (s == 0)
					uniquePatternWeights[u] = gsl_cdf_hypergeometric_P(uniquePatternCountsSub[u], uniquePatternCounts[u], sampleSize - uniquePatternCounts[u], nCases);
				//uniquePatternWeights[u] = gw_hypergeometric_cmf(uniquePatternCountsSub[u], uniquePatternCounts[u], sampleSize - uniquePatternCounts[u], nCases);
				else
					uniquePatternWeights[u] = gsl_cdf_hypergeometric_P(uniquePatternCountsSub[u], uniquePatternCounts[u], sampleSize - uniquePatternCounts[u], nCtrls);
				//uniquePatternWeights[u] = gw_hypergeometric_cmf(uniquePatternCountsSub[u], uniquePatternCounts[u], sampleSize - uniquePatternCounts[u], nCtrls);
			}

			//!- KBAC statistic: sum of genotype pattern frequencies differences in cases vs. controls, weighted by the hypergeometric distribution kernel
			double kbac = 0.0;
			for (unsigned u = 0; u != uniquePattern.size(); ++u) {
				if (s == 0)
					kbac = kbac + ( (1.0 * uniquePatternCountsSub[u]) / (1.0 * nCases) - (1.0 * (uniquePatternCounts[u] - uniquePatternCountsSub[u])) / (1.0 * nCtrls) ) * uniquePatternWeights[u];
				else
					kbac = kbac + ( (1.0 * uniquePatternCountsSub[u]) / (1.0 * nCtrls) - (1.0 * (uniquePatternCounts[u] - uniquePatternCountsSub[u])) / (1.0 * nCases) ) * uniquePatternWeights[u];
			}

			if (__isDebug)
				std::cout << kbac << std::endl;

			//FIXME
			//gw_round(kbac, 1E-3);
			kbacStatistics.push_back(kbac);
		}

		double statistic = 0.0;
		//!- one model statistic
		if (kbacStatistics.size() == 1) {
			statistic = kbacStatistics[0];
		}
		//!- two model statistic
		else if (kbacStatistics.size() == 2) {
			statistic = fmax(kbacStatistics[0], kbacStatistics[1]);
		}else {
			std::cerr << "ERROR: KBAC statistic error" << std::endl;
			exit(-1);
		}

		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (adaptive != 0)
				pvalue = m_checkAdaptivePvalue(permcount1, permcount2, iPermutation, adaptive, 0);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		random_shuffle(__ydat.begin(), __ydat.end());
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		pvalue = (1.0 + permcount1) / (1.0 + nPermutations);
	}
	return pvalue;
}


double gwAssociations::calcKbacstP(unsigned sided, unsigned nPermutations, unsigned adaptive)
{
	/*! * the KBAC in score test
	   *Implementation:
	 */
	if (sided == 0 || sided > 2)
		sided = 1;

	//!- trim data by mafs upper-lower bounds
	m_trimXdat();


	unsigned sampleSize = __ydat.size();
	// sample size
	unsigned regionLen = __xdat[0].size();
	// candidate region length
	unsigned nCases = 0;
	// case size

	for (unsigned i = 0; i != __ydat.size(); ++i) {
		if (__ydat[i] == AFFECTED)
			++nCases;
	}
	unsigned nCtrls = sampleSize - nCases;

	vectorF genotypeId(sampleSize);
	//!-Compute unique genotype patterns (string) as ID scores (double)

	for (unsigned i = 0; i != sampleSize; ++i) {

		double vntIdL = 0.0;
		double vntIdR = 0.0;
		const double ixiix = pow(9.0, 10.0);
		unsigned lastCnt = 0;
		unsigned tmpCnt = 0;

		for (unsigned j = 0; j != regionLen; ++j) {

			if (__xdat[i][j] != MISSING_ALLELE && __xdat[i][j] != MAJOR_ALLELE)
				vntIdR += pow(3.0, 1.0 * (j - lastCnt)) * __xdat[i][j];
			else
				continue;
			if (vntIdR >= ixiix) {
				vntIdL = vntIdL + 1.0;
				vntIdR = vntIdR - ixiix;
				lastCnt = lastCnt + tmpCnt + 1;
				tmpCnt = 0;
				continue;
			}else {
				++tmpCnt;
				continue;
			}
		}

		genotypeId[i] = vntIdL + vntIdR * 1e-10;
		// one-to-one "ID number" for a genotype pattern
	}

	// unique genotype patterns
	// use "list" and some STL algorithms

	std::list<double> uniqueId(genotypeId.begin(), genotypeId.end());
	uniqueId.remove(MAJOR_ALLELE);

	if (uniqueId.size() == 0)
		return 1.0;

	uniqueId.sort();
	uniqueId.unique();
	// remove wildtype and get unique genotype patterns
	//!- Unique genotype patterns that occur in the sample


	vectorF uniquePattern(uniqueId.size());
	copy(uniqueId.begin(), uniqueId.end(), uniquePattern.begin());
	uniqueId.clear();

	//  std::cout << uniquePattern << std::endl;

	// count number of sample individuals for each genotype pattern
	unsigned uniquePatternCounts[uniquePattern.size()];
	for (unsigned u = 0; u != uniquePattern.size(); ++u)
		uniquePatternCounts[u] = 0;

	for (unsigned i = 0; i != sampleSize; ++i) {
		// for each sample, identify/count its genotype pattern

		for (unsigned u = 0; u != uniquePattern.size(); ++u) {

			if (genotypeId[i] == uniquePattern[u]) {
				// genotype pattern identified
				++uniquePatternCounts[u];
				// count this genotype pattern
				break;
			}else ;
			// genotype pattern not found -- move on to next pattern
		}
	}


	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;
	while (iPermutation <= nPermutations) {

		// the KBAC statistic. Will be of length 1 or 2
		vectorF kbacStatistics(0);
		// two models
		for (unsigned s = 0; s != sided; ++s) {

			//!- count number of sample cases (for the 1st model, or ctrls for the 2nd model) for each genotype pattern
			unsigned uniquePatternCountsSub[uniquePattern.size()];
			for (unsigned u = 0; u != uniquePattern.size(); ++u)
				uniquePatternCountsSub[u] = 0;
			// genotype pattern counts in cases (for the 1st model, or ctrls for the 2nd model)

			for (unsigned i = 0; i != sampleSize; ++i) {
				if (__ydat[i] == (AFFECTED - 1.0 * s) ) {
					// for each "case (for the 1st model, or ctrls for 2nd model)", identify/count its genotype pattern
					for (unsigned u = 0; u != uniquePattern.size(); ++u) {
						if (genotypeId[i] == uniquePattern[u]) {
							// genotype pattern identified in cases (for the 1st model, or ctrls for 2nd model)
							++uniquePatternCountsSub[u];
							// count this genotype pattern
							break;
						}else ;
						// genotype pattern not found -- move on to next pattern
					}
				}else ;
			}

			//!- KBAC weights
			double uniquePatternWeights[uniquePattern.size()];
			// genotype pattern weights, the hypergeometric distribution cmf
			for (unsigned u = 0; u != uniquePattern.size(); ++u)
				uniquePatternWeights[u] = 0.0;

			for (unsigned u = 0; u != uniquePattern.size(); ++u) {
				if (s == 0)
					uniquePatternWeights[u] = gsl_cdf_hypergeometric_P(uniquePatternCountsSub[u], uniquePatternCounts[u], sampleSize - uniquePatternCounts[u], nCases);
				else
					uniquePatternWeights[u] = gsl_cdf_hypergeometric_P(uniquePatternCountsSub[u], uniquePatternCounts[u], sampleSize - uniquePatternCounts[u], nCtrls);
			}

			//!- regression score test using KBAC weights as regressors
			vectorF regressors(sampleSize, 0.0);
			for (unsigned i = 0; i != sampleSize; ++i) {
				for (unsigned u = 0; u != uniquePattern.size(); ++u) {
					if (genotypeId[i] == uniquePattern[u]) {
						regressors[i] = uniquePatternWeights[u];
						break;
					}else ;
				}
			}
			double xbar = gw_mean(regressors);
			double kbac = m_calcSimpleLogitRegScore(regressors, __ydat, xbar, nCases, 1);
			kbacStatistics.push_back(kbac);
		}

		double statistic = 0.0;
		//!- one model statistic
		if (kbacStatistics.size() == 1) {
			statistic = kbacStatistics[0];
		}
		//!- two model statistic
		else if (kbacStatistics.size() == 2) {
			statistic = fmax(kbacStatistics[0], kbacStatistics[1]);
		}else {
			std::cerr << "ERROR: KBAC statistic error" << std::endl;
			exit(-1);
		}

		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (adaptive != 0)
				pvalue = m_checkAdaptivePvalue(permcount1, permcount2, iPermutation, adaptive, 0);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		random_shuffle(__ydat.begin(), __ydat.end());
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		pvalue = (1.0 + permcount1) / (1.0 + nPermutations);
	}
	return pvalue;
}


double gwAssociations::calcVtP(unsigned sided, unsigned nPermutations, unsigned adaptive)
{
	/*! * Variable threshold method, Price 2010 AJHG <br>
	 * * <em> Unique feature of VT </em>: instead of using fixed threshold + Morris A RV counts, it uses a variable threshold for RV frequency. <br>
	 * * Since it uses Z_max = max(Z_i ...) where Z_i is proportional to standard normal by a coefficient which is not a function of RV counts, it does not matter to include the "known" common variants. <br>
	 * * In brief it is a multiple testing method that analyzes both common and rare variants <br><br>
	 * Implementation:
	 */

	// case size
	unsigned nCases = 0;

	for (unsigned i = 0; i != __ydat.size(); ++i) {
		if (__ydat[i] == AFFECTED)
			++nCases;
	}

	double po = (1.0 * nCases) / (1.0 * __ydat.size());

	// number of variants at each loci
	vectorUI nLocusVariants(__xdat[0].size(), 0);
	for (unsigned j = 0; j != __xdat[0].size(); ++j) {
		for (unsigned i = 0; i != __xdat.size(); ++i) {
			nLocusVariants[j] += (unsigned)__xdat[i][j];
		}
	}


	//!- Compute the unique non-zero variant counts at each loci, the "variable threshold"
	std::list<unsigned> lvts(nLocusVariants.begin(), nLocusVariants.end());
	lvts.remove(0);
	if (lvts.size() == 0)
		return 1.0;
	lvts.sort();
	lvts.unique();
	lvts.push_front(0);
	vectorUI vts(lvts.size());
	copy(lvts.begin(), lvts.end(), vts.begin());
	lvts.clear();


	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= nPermutations) {

		//! - Define <b> 'allZs' </b>, a vector of the Z scores computed under different thresholds
		vectorF allZs(0);
		vectorF zIa(__xdat.size(), 0.0), zIb(__xdat.size(), 0.0);

		//! - Iterate the following for all thresholds:
		for (unsigned t = 1; t != vts.size(); ++t) {

			//! - - Record the index of loci that can be added at the new threshold criteria
			vectorUI idxesAdding(0);
			for (unsigned s = 0; s != nLocusVariants.size(); ++s) {
				if (nLocusVariants[s] > vts[t - 1] && nLocusVariants[s] <= vts[t])
					idxesAdding.push_back(s);
			}

			//!- - For loci passing the threshold criteria, implement Price paper page 3 z(T) formula
			for (unsigned i = 0; i != __xdat.size(); ++i) {

				double zIai = 0.0, zIbi = 0.0;
				for (unsigned j = 0; j != idxesAdding.size(); ++j) {
					unsigned locIdx = idxesAdding[j];
					zIai += ((__ydat[i] - 1.0) - po) * __xdat[i][locIdx];
					zIbi += __xdat[i][locIdx] * __xdat[i][locIdx];
				}

				zIa[i] += zIai;
				zIb[i] += zIbi;
			}

			//! - Now each element in <b> 'allZs' </b> is z(T) for different thresholds
			allZs.push_back(gw_sum(zIa) / sqrt(gw_sum(zIb) ) );
		}

		//! - Compute zmax, the statistic; square it so that it would work for protectives
		if (sided == 0) {
			for (unsigned i = 0; i != allZs.size(); ++i) {
				allZs[i] = allZs[i] * allZs[i];
			}
		}

		double statistic = *max_element(allZs.begin(), allZs.end());

		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (adaptive != 0)
				pvalue = m_checkAdaptivePvalue(permcount1, permcount2, iPermutation, adaptive, sided);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		random_shuffle(__ydat.begin(), __ydat.end());
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		if (sided == 1 || sided == 0)
			pvalue = (1.0 + permcount1) / (1.0 + nPermutations);
		else {
			double permcount = gw_dmin(permcount1, permcount2);
			pvalue = (2.0 * permcount + 2.0) / (1.0 + nPermutations);
		}
	}
	return pvalue;
}


double gwAssociations::calcVtFisherP(unsigned sided, unsigned nPermutations, unsigned adaptive)
{
	/*! * Variable threshold method + Fisher's test <br><br>
	 * The test combines VT and CMC one-sided Fisher's with mid-p corrections. If none of the tests are significant (judged by Fisher's exact p-value) for the original data then no permutation test is pursued <br><br>
	 * Implementation:
	 */

	// number of variants at each loci
	vectorUI nLocusVariants(__xdat[0].size(), 0);

	for (unsigned j = 0; j != __xdat[0].size(); ++j) {
		for (unsigned i = 0; i != __xdat.size(); ++i) {
			nLocusVariants[j] += (unsigned)__xdat[i][j];
		}
	}

	//!- Compute the unique non-zero variant counts at each loci, the "variable threshold"
	std::list<unsigned> lvts(nLocusVariants.begin(), nLocusVariants.end());
	lvts.remove(0);
	if (lvts.size() == 0)
		return 1.0;
	lvts.sort();
	lvts.unique();
	lvts.push_front(0);
	vectorUI vts(lvts.size());
	copy(lvts.begin(), lvts.end(), vts.begin());
	lvts.clear();


	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;
	bool shouldPermutationTest = false;

	while (iPermutation <= nPermutations) {

		//! - Define <b> 'allPs' </b>, a vector of the p-values computed under different thresholds
		vectorF allPs(0);
		vectorF regressorsCurr(__xdat.size(), 0.0);
		//! - Iterate the following for all thresholds:
		for (unsigned t = 1; t != vts.size(); ++t) {

			//! - - Record the index of loci that can be added at the new threshold criteria
			vectorUI idxesAdding(0);
			for (unsigned s = 0; s != nLocusVariants.size(); ++s) {
				if (nLocusVariants[s] > vts[t - 1] && nLocusVariants[s] <= vts[t])
					idxesAdding.push_back(s);
			}

			//!- - For loci passing the threshold, implement CMC one-side Fisher
			vectorF regressors(__xdat.size(), 0.0);

			for (unsigned i = 0; i != regressors.size(); ++i) {
				for (unsigned j = 0; j != idxesAdding.size(); ++j) {
					unsigned locIdx = idxesAdding[j];
					regressors[i] = (__xdat[i][locIdx] == MAJOR_ALLELE || __xdat[i][locIdx] == MISSING_ALLELE) ? regressors[i] : 1.0;
				}
				regressorsCurr[i] = (regressors[i] == 1.0) ? 1.0 : regressorsCurr[i];
			}

			double pfisher = m_calc2X2Fisher(regressorsCurr, __ydat, sided);
			if (pfisher <= __alpha) {
				shouldPermutationTest = true;
			}

			allPs.push_back(-log(pfisher));
		}

		//!- an adaptive approach via Fisher's test
		if (shouldPermutationTest == false && iPermutation == 0) {
			return 1.0;
		}

		//! - Compute the statistic
		double statistic = *max_element(allPs.begin(), allPs.end());

		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (adaptive != 0)
				pvalue = m_checkAdaptivePvalue(permcount1, permcount2, iPermutation, adaptive, sided);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		random_shuffle(__ydat.begin(), __ydat.end());
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		if (sided == 1 || sided == 0)
			pvalue = (1.0 + permcount1) / (1.0 + nPermutations);
		else {
			double permcount = gw_dmin(permcount1, permcount2);
			pvalue = (2.0 * permcount + 2.0) / (1.0 + nPermutations);
		}
	}
	return pvalue;
}


double gwAssociations::calcAsumP(const char moi, unsigned sided, unsigned nPermutations, unsigned adaptive)
{
	/*! * Number of rare variants per site with "protective" site recoded, by Pan and Han (2010) Hum Hered <br>
	 *  * The authors use alpha0 = 0.1 to screen variants that should be recoded. See their paper for details <br>
	 *  * Original paper uses LR score test with approx. distribution, which is equivalent to Armitage's test for trend for this case (Schaid 2002). For RV not many sites are homozygote. In my implementation I use Fisher's 2by2 test instead.
	 * * Here allows user specified lower and upper bounds <br><br>
	 * Implementation:
	 */

	//!- trim data by mafs upper-lower bounds
	m_trimXdat();
	unsigned nCases = 0;
	for (unsigned i = 0; i != __ydat.size(); ++i)
		if (__ydat[i] == AFFECTED)
			++nCases;
	unsigned nCtrls = __ydat.size() - nCases;

	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= nPermutations) {

		//!- Sites that have excess rare variants in controls
		vectorL ctrlExcess(__xdat[0].size(), false);

		for (unsigned j = 0; j != __xdat[0].size(); ++j) {
			double tmpCase = 0.0;
			double tmpCtrl = 0.0;
			for (unsigned i = 0; i != __xdat.size(); ++i) {
				if (__ydat[i] == UNAFFECTED) {
					if (__xdat[i][j] != MISSING_ALLELE && __xdat[i][j] != MAJOR_ALLELE)
						tmpCtrl += __xdat[i][j];
				}else {
					if (__xdat[i][j] != MISSING_ALLELE && __xdat[i][j] != MAJOR_ALLELE)
						tmpCase += __xdat[i][j];
				}
			}
			if (tmpCtrl / (nCtrls * 2.0) > tmpCase / (nCases * 2.0))
				ctrlExcess[j] = true;
			else
				continue;
		}

		//!- Define sites that needs to be recoded
		vectorL recodeSites(__xdat[0].size(), false);

		for (unsigned j = 0; j != __xdat[0].size(); ++j) {
			if (ctrlExcess[j] == false)
				continue;

			vectorF vdat(0);
			for (unsigned i = 0; i != __xdat.size(); ++i)
				vdat.push_back(__xdat[i][j]);

			//!- 2 by 2 Fisher's test
			double pfisher = m_calc2X2Fisher(vdat, __ydat, moi, sided);

			if (pfisher < 0.1)
				recodeSites[j] = true;
			else
				continue;
		}

		vectorF regressors(0);

		//! - Count anrv; recode when necessary
		for (unsigned i = 0; i != __xdat.size(); ++i) {
			//  scan all sample individuals

			double nrv = 0.0;
			for (unsigned j = 0; j != __xdat[0].size(); ++j) {
				double tmpCode = 0.0;
				if (__xdat[i][j] != MISSING_ALLELE)
					tmpCode += __xdat[i][j];
				//!- Recode protective variants as in Pan 2010
				if (recodeSites[j] == true)
					tmpCode = 2.0 - tmpCode;
				// aggregated number of rare variants
				nrv += tmpCode;
			}

			regressors.push_back(nrv);
		}

		//!- Score test implementation for logistic regression model logit(p) = b0 + b1x (derivation of the test see my labnotes vol.2 page 3)

		if (__isDebug)
			std::cout << regressors << std::endl;
		//!- logistic regression score statistic
		double xbar = gw_mean(regressors);
		double statistic = m_calcSimpleLogitRegScore(regressors, __ydat, xbar, nCases, sided);

		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (adaptive != 0)
				pvalue = m_checkAdaptivePvalue(permcount1, permcount2, iPermutation, adaptive, sided);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		random_shuffle(__ydat.begin(), __ydat.end());
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		if (sided == 1 || sided == 0)
			pvalue = (1.0 + permcount1) / (1.0 + nPermutations);
		else {
			double permcount = gw_dmin(permcount1, permcount2);
			pvalue = (2.0 * permcount + 2.0) / (1.0 + nPermutations);
		}
	}
	return pvalue;
}


double gwAssociations::calcCmcqtP(unsigned sided, unsigned nPermutations, unsigned adaptive)
{
	//!- trim data by mafs upper-lower bounds
	m_trimXdat();

	//!- CMC scoring
	vectorF regressors = m_indicateRegionalVariants();

	if (__isDebug)
		std::cout << regressors << std::endl;


	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= nPermutations) {

		//!- two sample t test. Group 1: have cmc score = 1, Group 2: have cmc score = 0
		double statistic = m_calc2sampleT(regressors, __ydat, sided);
		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (adaptive != 0)
				pvalue = m_checkAdaptivePvalue(permcount1, permcount2, iPermutation, adaptive, sided);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		random_shuffle(__ydat.begin(), __ydat.end());
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		if (sided == 1 || sided == 0)
			pvalue = (1.0 + permcount1) / (1.0 + nPermutations);
		else {
			double permcount = gw_dmin(permcount1, permcount2);
			pvalue = (2.0 * permcount + 2.0) / (1.0 + nPermutations);
		}
	}
	return pvalue;
}


double gwAssociations::calcAnrvqtPermP(unsigned sided, unsigned nPermutations, unsigned adaptive)
{

	//!- trim data by mafs upper-lower bounds
	m_trimXdat();
	//!- ANRV scoring
	vectorF regressors = m_countRegionalVariants();
	if (__isDebug)
		std::cout << regressors << std::endl;

	//!- Linear regression test
	double xbar = gw_mean(regressors);
	double ybar = gw_mean(__ydat);

	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= nPermutations) {
		double statistic = m_calcSimpleLinearRegScore(regressors, __ydat, xbar, ybar, 0);
		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (adaptive != 0)
				pvalue = m_checkAdaptivePvalue(permcount1, permcount2, iPermutation, adaptive, sided);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		random_shuffle(__ydat.begin(), __ydat.end());
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		if (sided == 1 || sided == 0)
			pvalue = (1.0 + permcount1) / (1.0 + nPermutations);
		else {
			double permcount = gw_dmin(permcount1, permcount2);
			pvalue = (2.0 * permcount + 2.0) / (1.0 + nPermutations);
		}
	}
	return pvalue;
}


double gwAssociations::calcAnrvqtP(unsigned sided)
{

	//!- trim data by mafs upper-lower bounds
	m_trimXdat();
	//!- ANRV scoring
	vectorF regressors = m_countRegionalVariants();
	if (__isDebug)
		std::cout << regressors << std::endl;

	//!- Linear regression test
	double xbar = gw_mean(regressors);
	double ybar = gw_mean(__ydat);
	double pvalue = m_calcSimpleLinearRegScore(regressors, __ydat, xbar, ybar, sided);
	return pvalue;
}


double gwAssociations::calcAnrvqtP(double yh, double yl, unsigned sided)
{

	//!- trim data by mafs upper-lower bounds
	m_trimXdat();
	//!- ANRV scoring
	vectorF regressors = m_countRegionalVariants();
	if (__isDebug)
		std::cout << regressors << std::endl;

	//!- Linear regression conditional score test
	double pvalue = m_calcConditionalLinearRegScore(regressors, __ydat, yh, yl, sided);
	return pvalue;
}


double gwAssociations::calcTestRareP(unsigned sided, unsigned nPermutations, unsigned adaptive)
{
	/*! * the RVP Statistic: variant counts weighted by poisson kernels. <br>
	 * * It is a permutation based two model test, max(protective model, risk model).  <br>
	 * * See <em> Ionita-Laza and Lange 2011 PLoS Genet. </em> <br><br>
	 * * Implementation (credit to the authors for part of the codes below):
	 */

	//!- trim data by mafs upper-lower bounds
	m_trimXdat();

	unsigned sampleSize = __xdat.size();
	// sample size
	unsigned regionLen = __xdat[0].size();
	// candidate region length

	unsigned nCases = 0;
	for (unsigned i = 0; i != __ydat.size(); ++i) {
		if (__ydat[i] == AFFECTED)
			++nCases;
	}
	unsigned nCtrls = sampleSize - nCases;


	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= nPermutations) {

		double sumR = 0, sumP = 0;
		for (unsigned j = 0; j != regionLen; ++j) {

			//! - Count number of variants in cases/controls at a locus
			unsigned countcs = 0;
			unsigned countcn = 0;
			for (unsigned i = 0; i != sampleSize; ++i) {
				if (__ydat[i] == UNAFFECTED)
					countcn += (unsigned)__xdat[i][j];
				else
					countcs += (unsigned)__xdat[i][j];
			}

			//! - the RVP method. Codes adopted from the author's implementation
			double f = 1.0;
			//float w;
			//k0=(int)(freq*2*nCtrls);
			if (countcs > 0)
				f = gsl_cdf_poisson_P(countcn, nCtrls * (countcs + countcn) / (1.0 * (nCases + nCtrls))) * (1.0 - gsl_cdf_poisson_P(countcs - 1, nCases * (countcs + countcn) / (1.0 * (nCases + nCtrls))));

			if ((1.0 * countcs) / (1.0 * nCases) > (1.0 * countcn) / (1.0 * nCtrls))
				sumR -= log(f);

			f = 1.0;
			//k0=(int)(freq*2*nCases);
			if (countcn > 0)
				f = gsl_cdf_poisson_P(countcs, nCases * (countcs + countcn) / (1.0 * (nCases + nCtrls))) * (1.0 - gsl_cdf_poisson_P(countcn - 1, nCtrls * (countcs + countcn) / (1.0 * (nCases + nCtrls))));
			if ((1.0 * countcn) / (1.0 * nCtrls) > (1.0 * countcs) / (1.0 * nCases))
				sumP -= log(f);
		}

		//!- RVP statistic: The max statistic in the manuscript: R - potentially risk and P - potentially protective
		double statistic = (sided == 1) ? sumR : fmax(sumR, sumP);

		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (adaptive != 0)
				pvalue = m_checkAdaptivePvalue(permcount1, permcount2, iPermutation, adaptive, sided);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		random_shuffle(__ydat.begin(), __ydat.end());
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else
		pvalue = (1.0 + permcount1) / (1.0 + nPermutations);

	return pvalue;
}


double gwAssociations::calcCalphaP(unsigned sided, unsigned nPermutations, unsigned adaptive)
{

	/*! * the c-alpha Statistic: sum of the std. error of variant counts in cases <br>
	 * *  One-sided test, claims to be robust against protective mutations. <br>
	 * * See <em> Ben. Neale et al. 2011 PLoS Genet. </em> <br><br>
	 * * Implementation:
	 */

	//!- trim data by mafs upper-lower bounds
	m_trimXdat();

	unsigned sampleSize = __xdat.size();
	// sample size
	unsigned regionLen = __xdat[0].size();
	// candidate region length

	unsigned nCases = 0;
	for (unsigned i = 0; i != __ydat.size(); ++i) {
		if (__ydat[i] == AFFECTED)
			++nCases;
	}

	double phat = (1.0 * nCases) / (1.0 * sampleSize);


	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= nPermutations) {

		double calpT = 0.0;
		double calpV = 0.0;
		unsigned singletonAll = 0;
		unsigned singletonCases = 0;
		bool isEmptyData = true;

		for (unsigned j = 0; j != regionLen; ++j) {

			//! - Count number of variants in cases/controls at a locus
			unsigned countcs = 0;
			unsigned countcn = 0;
			for (unsigned i = 0; i != sampleSize; ++i) {
				if (__ydat[i] == UNAFFECTED)
					countcn += (unsigned)__xdat[i][j];
				else
					countcs += (unsigned)__xdat[i][j];
			}

			//! - the c-alpha method implementation
			unsigned ni = countcs + countcn;
			if (ni < 2) {
				singletonAll += ni;
				singletonCases += countcs;
				continue;
			}else
				isEmptyData = false;
			//!- * skip singletons

			double niv = ni * phat * (1 - phat);
			calpT += (countcs - ni * phat) * (countcs - ni * phat) - niv;
			for (unsigned u = 0; u <= ni; ++u) {
				double tmess = (u - ni * phat) * (u - ni * phat) - niv;
				calpV += tmess * tmess * gsl_ran_binomial_pdf(u, phat, ni);
			}
		}

		if (singletonAll >= 2) {
			isEmptyData = false;
			double niv = singletonAll * phat * (1 - phat);
			calpT += (singletonCases - singletonAll * phat) * (singletonCases - singletonAll * phat) - niv;
			for (unsigned u = 0; u <= singletonAll; ++u) {
				double tmess = (u - singletonAll * phat) * (u - singletonAll * phat) - niv;
				calpV += tmess * tmess * gsl_ran_binomial_pdf(u, phat, singletonAll);
			}
		} else ;
		//!- * bin singletons

		if (isEmptyData)
			return 1.0;

		//!- c-alpha statistic, two sided.
		double statistic = calpT * calpT / calpV;
		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (adaptive != 0)
				pvalue = m_checkAdaptivePvalue(permcount1, permcount2, iPermutation, adaptive, sided);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		random_shuffle(__ydat.begin(), __ydat.end());
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		if (sided == 1 || sided == 0)
			pvalue = (1.0 + permcount1) / (1.0 + nPermutations);
		else {
			double permcount = gw_dmin(permcount1, permcount2);
			pvalue = (2.0 * permcount + 2.0) / (1.0 + nPermutations);
		}
	}
	return pvalue;
}


double gwAssociations::calcRareCoverP(unsigned sided, unsigned nPermutations, unsigned adaptive)
{
	/*! * RareCover method, 2010 PLoS CompBio <br>
	 * * Implementation:
	 */

	//!- the cut-off to use for the "heuristic greedy algorithm". = 0.5 as suggested by the paper
	const double difQ = 0.5;

	//!- Index of loci that are observed to be polymophic
	vectorUI vntVct(0);

	for (unsigned j = 0; j != __observedMafs.size(); ++j) {
		if (__observedMafs[j] > 0.0)
			vntVct.push_back(j + 1);
	}

	if (vntVct.size() == 0)
		return 1.0;

	unsigned smpSize = __xdat.size();


	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= nPermutations) {

		//!- the current and the next statistics
		double sCurr = 0.0, sNext = 0.0;
		vectorF regressors(__ydat.size(), 0.0);
		//!- the current and the next genotype coding
		vectorF regressorsCurr(__ydat.size(), 0.0);
		vectorUI vntNow = vntVct;


		do {
			sCurr = sNext;
			unsigned rmIdx = 0;
			//!- the "test contributing" variant index, for the vntVct object
			bool rmIdxFlag = false;

			for (unsigned t = 0; t != vntNow.size(); ++t) {

				if (vntNow[t] == 0)
					continue;
				unsigned iIdx = vntNow[t] - 1;
				//!- the index of a variant site

				for (unsigned i = 0; i != smpSize; ++i)
					regressors[i] = (regressorsCurr[i] + __xdat[i][iIdx] > 0) ? MINOR_ALLELE : MAJOR_ALLELE;

				double statistic;

				//! - 2 by 2 one-sided Fisher's test
				if (sided == 1) {
					statistic = -1.0 * log(m_calc2X2Fisher(regressors, __ydat, sided));
				}
				//! - 2 by 2 Chisq test
				else {
					statistic = m_calc2X2Chisq(regressors, __ydat);
				}

				if (statistic > sNext) {
					sNext = statistic;
					rmIdx = t;
					rmIdxFlag = true;
				} else continue;
			}
			//!- Now end up with a properly updated sNext statistic.

			// in this case, sNext is not updated, OK to break the loop.
			if (rmIdxFlag == false)
				break;

			else {
				unsigned rmVnt = vntNow[rmIdx] - 1;
				//!- Update the genotype coding by adding in the contributing locus
				for (unsigned i = 0; i != smpSize; ++i)
					regressorsCurr[i] = regressorsCurr[i] + __xdat[i][rmVnt];

				//!- remove the contributing locus to avoid duplicated visit to it the next time.
				vntNow[rmIdx] = 0;
			}
		} while (sNext - sCurr > difQ);

		//!- Test statistic 'sNext'
		double statistic = sNext;

		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (adaptive != 0)
				pvalue = m_checkAdaptivePvalue(permcount1, permcount2, iPermutation, adaptive, sided);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		random_shuffle(__ydat.begin(), __ydat.end());
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		if (sided == 1 || sided == 0)
			pvalue = (1.0 + permcount1) / (1.0 + nPermutations);
		else {
			double permcount = gw_dmin(permcount1, permcount2);
			pvalue = (2.0 * permcount + 2.0) / (1.0 + nPermutations);
		}
	}
	return pvalue;
}


double gwAssociations::calcWsFisherP(unsigned sided, unsigned nPermutations, unsigned adaptive, bool isMidPvalue, const char moi)
{
	/*! * Weighted sum Fisher test, Wang Shuang and Patrick (2011) <br>
	 * Implementation:
	 */

	//!- trim data by mafs upper-lower bounds
	m_trimXdat();

	unsigned nAllCases = 0;
	for (unsigned i = 0; i != __ydat.size(); ++i)
		if (__ydat[i] == AFFECTED)
			++nAllCases;
	unsigned nAllCtrls = __ydat.size() - nAllCases;


	unsigned regionLen = __xdat[0].size();

	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= nPermutations) {

		//!- Sites that have excess rare variants in cases = 1, in controls = -1, equal = 0
		vectorI excessIndicators(regionLen, 0);
		//!- -log(p)
		vectorF minusLogPvalues(regionLen, -9.0);
		vectorF weights(regionLen, 0.0);

		for (unsigned j = 0; j != regionLen; ++j) {
			// for dealing with missing data
			unsigned nCases = nAllCases;
			unsigned nCtrls = nAllCtrls;
			double allelesCase = 0.0;
			double allelesCtrl = 0.0;
			double multiplier = 1.0;

			for (unsigned i = 0; i != __xdat.size(); ++i) {

				if (__xdat[i][j] == MISSING_ALLELE) {
					if (__ydat[i] == UNAFFECTED) --nCtrls;
					else --nCases;
					continue;
				}

				switch (moi) {
				case 'R':
				{
					if (__ydat[i] == UNAFFECTED)
						allelesCtrl += (__xdat[i][j] == HOMO_ALLELE) ? MINOR_ALLELE : MAJOR_ALLELE;
					else
						allelesCase += (__xdat[i][j] == HOMO_ALLELE) ? MINOR_ALLELE : MAJOR_ALLELE;
				}
				break;
				case 'D':
				{
					if (__ydat[i] == UNAFFECTED)
						allelesCtrl += (__xdat[i][j] != MAJOR_ALLELE) ? MINOR_ALLELE : MAJOR_ALLELE;
					else
						allelesCase += (__xdat[i][j] != MAJOR_ALLELE) ? MINOR_ALLELE : MAJOR_ALLELE;
				}
				break;
				default:
				{
					multiplier = 2.0;
					if (__ydat[i] == UNAFFECTED)
						allelesCtrl += __xdat[i][j];
					else
						allelesCase += __xdat[i][j];
				}
				break;
				}
			}

			if (allelesCtrl / (nCtrls * multiplier) > allelesCase / (nCases * multiplier)) {
				excessIndicators[j] = -1;
				// Weight
				double qi = (allelesCase + 1.0) / (multiplier * nCases + 2.0);
				weights[j] = 1.0 / sqrt((nCases + nCtrls) * multiplier * qi * (1.0 - qi));
				// p-value
				double tmpPvalue =
				    (isMidPvalue)
				    ? (allelesCase > 0.0) * gsl_cdf_hypergeometric_P((unsigned)(allelesCase - 1.0), (unsigned)(allelesCtrl + allelesCase),
						(unsigned)((nCases + nCtrls) * multiplier - allelesCtrl - allelesCase), (unsigned)(nCases * multiplier))
				    + 0.5 * gsl_ran_hypergeometric_pdf((unsigned)(allelesCase), (unsigned)(allelesCtrl + allelesCase),
						(unsigned)((nCases + nCtrls) * multiplier - allelesCtrl - allelesCase), (unsigned)(nCases * multiplier))
					: gsl_cdf_hypergeometric_P((unsigned)(allelesCase), (unsigned)(allelesCtrl + allelesCase),
						(unsigned)((nCases + nCtrls) * multiplier - allelesCtrl - allelesCase), (unsigned)(nCases * multiplier));

				minusLogPvalues[j] = -1.0 * log(tmpPvalue);
			}else if (allelesCtrl / (nCtrls * multiplier) < allelesCase / (nCases * multiplier)) {
				excessIndicators[j] = 1;
				// Weight
				double qi = (allelesCtrl + 1.0) / (multiplier * nCtrls + 2.0);
				weights[j] = 1.0 / sqrt((nCtrls + nCases) * multiplier * qi * (1.0 - qi));
				// p-value
				double tmpPvalue =
				    (isMidPvalue)
				    ? (allelesCtrl > 0.0) * gsl_cdf_hypergeometric_P((unsigned)(allelesCtrl - 1.0), (unsigned)(allelesCtrl + allelesCase),
						(unsigned)((nCases + nCtrls) * multiplier - allelesCtrl - allelesCase), (unsigned)(nCtrls * multiplier))
				    + 0.5 * gsl_ran_hypergeometric_pdf((unsigned)(allelesCtrl), (unsigned)(allelesCtrl + allelesCase),
						(unsigned)((nCases + nCtrls) * multiplier - allelesCtrl - allelesCase), (unsigned)(nCtrls * multiplier))
					: gsl_cdf_hypergeometric_P((unsigned)(allelesCtrl), (unsigned)(allelesCtrl + allelesCase),
						(unsigned)((nCases + nCtrls) * multiplier - allelesCtrl - allelesCase), (unsigned)(nCtrls * multiplier));

				minusLogPvalues[j] = -1.0 * log(tmpPvalue);
			}else
				continue;
		}

		vectorF statistics(2, 0.0);
		for (unsigned j = 0; j != regionLen; ++j) {
			if (excessIndicators[j] == 1) {
				statistics[0] += weights[j] * minusLogPvalues[j];
			}else if (excessIndicators[j] == -1) {
				statistics[1] += weights[j] * minusLogPvalues[j];
			}else
				continue;
		}

		double statistic = 0.0;
		if (sided == 0)
			statistic = fmax(statistics[0], statistics[1]);
		else if (sided == 1)
			statistic = statistics[0];
		else {
			std::cerr << "ERROR: WsFisher alternative error" << std::endl;
			exit(-1);
		}

		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (adaptive != 0)
				pvalue = m_checkAdaptivePvalue(permcount1, permcount2, iPermutation, adaptive, sided);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		random_shuffle(__ydat.begin(), __ydat.end());
		++iPermutation;
	}
	if (pvalue <= 1.0) ;
	else {
		pvalue = (1.0 + permcount1) / (1.0 + nPermutations);
	}
	return pvalue;
}


double gwAssociations::calcSkatP(unsigned nPermutations, unsigned adaptive)
{

#ifdef SKAT_GSL
	gsl_set_error_handler_off();
#endif
	unsigned sided = 1;
	//!- trim data by mafs upper-lower bounds
	m_trimXdat();
	// case size
	unsigned nCases = 0;
	for (unsigned i = 0; i != __ydat.size(); ++i) {
		if (__ydat[i] == AFFECTED)
			++nCases;
	}
	double po = (1.0 * nCases) / (1.0 * __ydat.size());
	//
	// observed maf and weights
	//
	vectorF obsMafs(__xdat[0].size(), 0.0);
	for (unsigned j = 0; j != __xdat[0].size(); ++j) {
		for (unsigned i = 0; i != __xdat.size(); ++i) {
			obsMafs[j] += __xdat[i][j] / (__ydat.size() * 2.0);
		}
	}
	vectorF weights(__xdat[0].size(), 0.0);
	for (unsigned j = 0; j != __xdat[0].size(); ++j) {
		weights[j] = pow(gsl_ran_beta_pdf(obsMafs[j], 1.0, 25.0), 2.0);
	}
	//
	// GWG'
	//
	vector2F gmat = __xdat;
	for (unsigned i = 0; i != gmat.size(); ++i) {
		std::transform(gmat[i].begin(), gmat[i].end(), weights.begin(), gmat[i].begin(), std::multiplies<double>());
	}

#ifdef SKAT_GSL
	gsl_matrix * kmat = gsl_matrix_alloc(__xdat.size(), __xdat.size());

	for (unsigned i = 0; i != __xdat.size(); ++i) {
		for (unsigned k = 0; k != __xdat.size(); ++k) {
			double gmatii = 0.0;
			for (unsigned j = 0; j != gmat[i].size(); ++j) {
				gmatii += gmat[i][j] * __xdat[k][j];
			}
			gsl_matrix_set(kmat, i, k, gmatii);
		}
	}
#else
	vector2F kmat(gmat.size());
	for (unsigned i = 0; i != kmat.size(); ++i) {
		for (unsigned k = 0; k != __xdat.size(); ++k) {
			double gmatii = 0.0;
			for (unsigned j = 0; j != gmat[i].size(); ++j) {
				gmatii += gmat[i][j] * __xdat[k][j];
			}
			kmat[i].push_back(gmatii);
		}
	}
#endif
	vectorF udat = __ydat;
	for (size_t i = 0; i < __ydat.size(); ++i) {
		udat[i] -= (1.0 + po);
	}

	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= nPermutations) {
		//
		// Q
		//
#ifdef SKAT_GSL
		gsl_vector * qvec = gsl_vector_alloc(udat.size());
		gsl_vector_view uvec = gsl_vector_view_array(&udat[0], udat.size());
		int m_err = gsl_blas_dgemv(CblasTrans, 1.0, kmat, &uvec.vector, 0.0, qvec);
		if (m_err != 0) {
			std::cerr << "Error in gsl_blas_dgemv(CblasTrans, 1.0, x, y, 0.0, b)" << std::endl;
			exit(-1);
		}

		gsl_vector_mul(qvec, &uvec.vector);
		double statistic = 0.0;
		for (size_t i = 0; i < __ydat.size(); ++i) {
			statistic += gsl_vector_get(qvec, i);
		}
		gsl_vector_free(qvec);
#else
		vectorF qvec(0);
		for (unsigned j = 0; j != kmat[0].size(); ++j) {
			double tmp = 0.0;
			for (unsigned i = 0; i != udat.size(); ++i) {
				tmp += udat[i] * kmat[i][j];
			}
			qvec.push_back(tmp);
		}

		std::transform(qvec.begin(), qvec.end(), udat.begin(), qvec.begin(), std::multiplies<double>());
		double statistic = std::accumulate(qvec.begin(), qvec.end(), 0.0, std::plus<double>());
#endif
		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (adaptive != 0)
				pvalue = m_checkAdaptivePvalue(permcount1, permcount2, iPermutation, adaptive, sided);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		random_shuffle(udat.begin(), udat.end());
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		if (sided == 1 || sided == 0)
			pvalue = (1.0 + permcount1) / (1.0 + nPermutations);
		else {
			double permcount = gw_dmin(permcount1, permcount2);
			pvalue = (2.0 * permcount + 2.0) / (1.0 + nPermutations);
		}
	}
#ifdef SKAT_GSL
	gsl_matrix_free(kmat);
#endif
	return pvalue;
}


//////////
// Private member functions
//////////

void gwAssociations::m_trimXdat()
{
	vector2F xdat = __xdat;

	__xdat.clear();
	__xdat.resize(xdat.size());

	for (unsigned j = 0; j != __observedMafs.size(); ++j) {
		if (__observedMafs[j] <= __mafLower || __observedMafs[j] > __mafUpper)
			continue;

		else {
			for (unsigned i = 0; i != xdat.size(); ++i)
				__xdat[i].push_back(xdat[i][j]);
		}
	}
	return;
}


void gwAssociations::m_maskWildtypeSibpair()
{
	for (unsigned i = 0; i < (__xdat.size() - 1); ) {
		for (unsigned j = 0; j < __observedMafs.size(); ++j) {
			if (__xdat[i][j] <= __xdat[i + 1][j]) {
				__xdat[i][j] = 0.0;
				__xdat[i + 1][j] = 0.0;
			}
		}
		i += 2;
	}
	return;
}


double gwAssociations::m_calcSimpleLogitRegScore(const vectorF & regressors, const vectorF & responses, double xbar, unsigned nCases, unsigned sided) const
{
	assert(regressors.size() == responses.size());

	//!- Score test implementation for logistic regression model logit(p) = b0 + b1x (derivation of the test see my labnotes vol.2 page 3)

	//double ebo = (1.0 * nCases) / (1.0 * numCtrl);
	//double bo = log(ebo);

	double po = (1.0 * nCases) / (1.0 * responses.size());
	double statistic = 0.0;
	double ss = 0.0, vm1 = 0.0;
	// the score and its variance
	for (unsigned i = 0; i != regressors.size(); ++i) {
		ss += (regressors[i] - xbar) * ((responses[i] - 1.0) - po);
		vm1 += (regressors[i] - xbar) * (regressors[i] - xbar) * po * (1.0 - po);
	}


	if (!fsame(ss, 0.0)) {
		statistic = ss / sqrt(vm1);
		if (sided != 1) {
			statistic = statistic * statistic;
		}
	}

	if (__isDebug) {
		std::cout << std::endl;
		std::cout << xbar << std::endl;
		std::cout << ss << std::endl;
		std::cout << vm1 << std::endl;
		std::cout << statistic << std::endl;
		exit(1);
	}

	gw_round(statistic, 1E-3);
	return statistic;
}


double gwAssociations::m_calc2X2Chisq(const vectorF & regressors, const vectorF & responses) const
{
	assert(regressors.size() == responses.size());

	//! - 2 by 2 Chisq test
	double A0 = 0.0, A1 = 0.0, U0 = 0.0, U1 = 0.0;
	for (unsigned i = 0; i != regressors.size(); ++i) {
		if (regressors[i] == MAJOR_ALLELE && responses[i] == UNAFFECTED)
			U0 += 1.0;
		else if (regressors[i] == MAJOR_ALLELE && responses[i] == AFFECTED)
			A0 += 1.0;
		else if (regressors[i] == MINOR_ALLELE && responses[i] == UNAFFECTED)
			U1 += 1.0;
		else if (regressors[i] == MINOR_ALLELE && responses[i] == AFFECTED)
			A1 += 1.0;
		else {
			std::cerr << "Input data problem in m_calc2X2Chisq(). Now Quit." << std::endl;
			exit(1);
		}
	} // collect the contigency table


	double Aobs = A0 + A1;
	double Uobs = U0 + U1;
	double Tobs = Aobs + Uobs;
	double Obs0 = A0 + U0;
	double Obs1 = A1 + U1;

	double EA0 = (Aobs * Obs0) / Tobs; if (EA0 == 0) EA0 = 0.05;
	double EA1 = (Aobs * Obs1) / Tobs; if (EA1 == 0) EA1 = 0.05;
	double EU0 = (Uobs * Obs0) / Tobs; if (EU0 == 0) EU0 = 0.05;
	double EU1 = (Uobs * Obs1) / Tobs; if (EU1 == 0) EU1 = 0.05;

	double statistic = ( (A0 - EA0) * (A0 - EA0) ) / EA0
	                   + ( (A1 - EA1) * (A1 - EA1) ) / EA1
	                   + ( (U0 - EU0) * (U0 - EU0) ) / EU0
	                   + ( (U1 - EU1) * (U1 - EU1) ) / EU1;

	return statistic;
}


double gwAssociations::m_calc2X2Fisher(const vectorF & regressors, const vectorF & responses, const char moi, unsigned sided) const
{
	assert(regressors.size() == responses.size());

	vectorI twotwoTable(4, 0);

	for (unsigned i = 0; i != regressors.size(); ++i) {

		if (regressors[i] != MAJOR_ALLELE && regressors[i] != MINOR_ALLELE && regressors[i] != HOMO_ALLELE) {
			std::cerr << "Input data problem in m_calc2X2Fisher() (X data have missing entries). Now Quit." << std::endl;
			exit(-1);
		}

		if (responses[i] != AFFECTED && responses[i] != UNAFFECTED) {
			std::cerr << "Input data problem in m_calc2X2Fisher() table (Y data not binary). Now Quit." << std::endl;
			exit(-1);
		}

		switch (moi) {
		case 'R':
		{
			if (responses[i] == AFFECTED) {
				if (regressors[i] != HOMO_ALLELE)
					twotwoTable[0] += 1;
				else
					twotwoTable[1] += 1;
			}else {
				if (regressors[i] != HOMO_ALLELE)
					twotwoTable[2] += 1;
				else
					twotwoTable[3] += 1;
			}
		}
		break;
		case 'D':
		{
			if (responses[i] == AFFECTED) {
				if (regressors[i] == MAJOR_ALLELE)
					twotwoTable[0] += 1;
				else
					twotwoTable[1] += 1;
			}else {
				if (regressors[i] == MAJOR_ALLELE)
					twotwoTable[2] += 1;
				else
					twotwoTable[3] += 1;
			}
		}
		break;
		default:
		{
			if (responses[i] == AFFECTED) {
				if (regressors[i] == MAJOR_ALLELE)
					twotwoTable[0] += 2;
				else
					twotwoTable[1] += (unsigned)regressors[i];
			}else {
				if (regressors[i] == MAJOR_ALLELE)
					twotwoTable[2] += 2;
				else
					twotwoTable[3] += (unsigned)regressors[i];
			}
		}
		break;
		}
	}

	double pvalue2X2 = 1.0;
	if (sided == 1) {
		pvalue2X2 = (twotwoTable[3] > 0) * gsl_cdf_hypergeometric_P((twotwoTable[3] - 1), (twotwoTable[1] + twotwoTable[3]), (twotwoTable[0] + twotwoTable[2]), (twotwoTable[3] + twotwoTable[2]))
		            + 0.5 * gsl_ran_hypergeometric_pdf(twotwoTable[3], (twotwoTable[1] + twotwoTable[3]), (twotwoTable[0] + twotwoTable[2]), (twotwoTable[3] + twotwoTable[2]));
	}else {
		pvalue2X2 = fexact_two_sided_pvalue(twotwoTable);
	}
	return pvalue2X2;
}


double gwAssociations::m_calc2X2Fisher(const vectorF & regressors, const vectorF & responses, unsigned sided) const
{
	assert(regressors.size() == responses.size());

	vectorI twotwoTable(4, 0);

	for (unsigned i = 0; i != regressors.size(); ++i) {

		if (regressors[i] != MAJOR_ALLELE && regressors[i] != MINOR_ALLELE && regressors[i] != HOMO_ALLELE) {
			std::cerr << "Input data problem in m_calc2X2Fisher() (X data have missing entries). Now Quit." << std::endl;
			exit(-1);
		}

		if (responses[i] != AFFECTED && responses[i] != UNAFFECTED) {
			std::cerr << "Input data problem in m_calc2X2Fisher() table (Y data not binary). Now Quit." << std::endl;
			exit(-1);
		}

		if (responses[i] == AFFECTED) {
			if (regressors[i] == MAJOR_ALLELE)
				twotwoTable[0] += 1;
			else
				twotwoTable[1] += 1;
		}else {
			if (regressors[i] == MAJOR_ALLELE)
				twotwoTable[2] += 1;
			else
				twotwoTable[3] += 1;
		}
	}

	double pvalue2X2 = 1.0;
	if (sided == 1) {
		pvalue2X2 = (twotwoTable[3] > 0) * gsl_cdf_hypergeometric_P((twotwoTable[3] - 1), (twotwoTable[1] + twotwoTable[3]), (twotwoTable[0] + twotwoTable[2]), (twotwoTable[3] + twotwoTable[2]))
		            + 0.5 * gsl_ran_hypergeometric_pdf(twotwoTable[3], (twotwoTable[1] + twotwoTable[3]), (twotwoTable[0] + twotwoTable[2]), (twotwoTable[3] + twotwoTable[2]));
	}else {
		pvalue2X2 = fexact_two_sided_pvalue(twotwoTable);
	}

	if (__isDebug) {
		std::cout << twotwoTable << std::endl;
		std::cout << pvalue2X2 << std::endl;
	}

	return pvalue2X2;
}


double gwAssociations::m_calcSimpleLinearRegScore(const vectorF & regressors, const vectorF & responses, double xbar, double ybar, unsigned sided) const
{
	assert(regressors.size() == responses.size());

	//!- Statistic: LSE (MLE) for beta, centered and scaled (bcz E[b] = 0 and sigma = 1 by simulation)
	//!- See page 41 of Kutner's Applied Linear Stat. Model, 5th ed.
	//
	double statistic = 0.0;
	double numerator = 0.0, denominator = 0.0, ysigma = 0.0;
	for (unsigned i = 0; i != regressors.size(); ++i) {
		numerator += (regressors[i] - xbar) * responses[i];
		denominator += pow(regressors[i] - xbar, 2.0);
	}

	if (!fsame(numerator, 0.0)) {
		//!- Compute MSE and V[\hat{beta}]
		//!- V[\hat{beta}] = MSE / denominator
		double b1 = numerator / denominator;
		double b0 = ybar - b1 * xbar;
		//SSE
		for (size_t i = 0; i != regressors.size(); ++i) {
			ysigma += pow(responses[i] - (b0 + b1 * regressors[i]), 2.0);
		}

		double varb = ysigma / ((responses.size() - 2.0) * denominator);
		statistic = b1 / sqrt(varb);
	}

	if (sided == 1) {
		statistic = gsl_cdf_ugaussian_Q(statistic);
	}else if (sided == 2) {
		statistic = statistic * statistic;
		statistic = gsl_cdf_chisq_Q(statistic, 1.0);
	}else ;  // report the score, not the p-value

	return statistic;
}


double gwAssociations::m_calcConditionalLinearRegScore(const vectorF & regressors, const vectorF & responses, double yh, double yl, unsigned sided) const
{
	assert(regressors.size() == responses.size());
	//!- Statistic: score test statistic for beta, given Y~N(0,1) under the null
	//!- a conditional test due to extreme sampling, Lin and Huang 2007
	//!- See my labnotes vol. 1
	//
	double fcdf = gsl_cdf_ugaussian_Q(yl) - gsl_cdf_ugaussian_Q(yh) + 1.0;
	double fpdf_1 = gsl_ran_ugaussian_pdf(yh);
	double fpdf_2 = gsl_ran_ugaussian_pdf(yl);
	double ratio = (fpdf_1 - fpdf_2) / fcdf;
	double numerator = 0.0, denominator = 0.0;
	for (unsigned i = 0; i != regressors.size(); ++i) {
		numerator += regressors[i] * (responses[i] + ratio);
		denominator += pow(regressors[i], 2.0) * (-1.0 + ratio * ratio + (fpdf_1 * yh - fpdf_2 * yl) / fcdf);
	}

	if (fsame(denominator, 0.0)) denominator = 1.0e-10;

	double statistic = numerator / sqrt(fabs(denominator));

	if (sided == 1) {
		statistic = gsl_cdf_ugaussian_Q(statistic);
	}else if (sided == 2) {
		statistic = statistic * statistic;
		std::cout << statistic << std::endl;
		statistic = gsl_cdf_chisq_Q(statistic, 1.0);
	}else ;  // report the score, not the p-value

	return statistic;
}


double gwAssociations::m_calc2sampleT(const vectorF & x1s, const vectorF & x2s, unsigned sided) const
{
	assert(x1s.size() == x2s.size());

	double XA = 0.0, XU = 0.0, SA = 0.0, SU = 0.0;
	vectorF VA(0), VU(0);
	for (unsigned i = 0; i != x1s.size(); ++i) {
		// genotype codings are 0 and 1 while phenotype is continuous
		if (x1s[i] == MAJOR_ALLELE)
			VU.push_back(x2s[i]);
		else if (x1s[i] == MINOR_ALLELE)
			VA.push_back(x2s[i]);
		else {
			std::cerr << "Input data problem in m_calc2sampleT(). Now Quit." << std::endl;
			exit(1);
		}
	}

	XA = gw_mean(VA);
	XU = gw_mean(VU);
	SA = gw_var(VA);
	SU = gw_var(VU);
	if (SA == 0) SA = 1.0e-6;
	if (SU == 0) SU = 1.0e-6;

	double statistic = (XA - XU) / sqrt(SA / VA.size() + SU / VU.size());
	if (sided != 1)
		statistic = statistic * statistic;
	return statistic;
}


vectorF gwAssociations::m_countRegionalVariants() const
{
	/*! * Number of rare variants per site, by Morris A (2009) Genet Epi <br>
	 * * The authors use the MAF cut-off of 5% but here allows user specified lower and upper bounds <br><br>
	 */

	vectorF data(0);

	//!- Counts of variants in a region
	for (unsigned i = 0; i != __xdat.size(); ++i) {

		double nrv = 0.0;
		for (unsigned j = 0; j != __xdat[i].size(); ++j) {
			if (__xdat[i][j] != MISSING_ALLELE)
				nrv += __xdat[i][j];
		}
		data.push_back(nrv);
	}
	return data;
}


vectorF gwAssociations::m_indicateRegionalVariants() const
{

	vectorF data(0);

	//!- Collapsing of variants in a region
	for (unsigned i = 0; i != __xdat.size(); ++i) {

		bool isWild = true;
		for (unsigned j = 0; j != __xdat[i].size(); ++j) {
			// scan all genetic loci of an individual

			//! - Apply indicator function: rare variant found. score it '1' and break the loop
			if (__xdat[i][j] != MAJOR_ALLELE && __xdat[i][j] != MISSING_ALLELE) {
				data.push_back(1.0);
				isWild = false;
				break;
			}else ;
		}

		//! - Otherwise score it '0' after every locus of an indv has been scanned.
		if (isWild == true)
			data.push_back(0.0);
	}

	return data;
}


vectorF gwAssociations::m_indicateRegionalUniqueVariants() const
{

	vectorF data(0);

	//!- Identify variants observed only in cases or in controls. Exclude them otherwise

	vectorL areUnique(__xdat[0].size(), true);

	for (unsigned j = 0; j != __xdat[0].size(); ++j) {

		bool caseFlag = false, ctrlFlag = false;
		for (unsigned i = 0; i != __xdat.size(); ++i) {

			if (__ydat[i] == AFFECTED && caseFlag == false) {
				if (__xdat[i][j] != MAJOR_ALLELE && __xdat[i][j] != MISSING_ALLELE)
					caseFlag = true;
				else ;
			}else if (__ydat[i] == UNAFFECTED && ctrlFlag == false) {
				if (__xdat[i][j] != MAJOR_ALLELE && __xdat[i][j] != MISSING_ALLELE)
					ctrlFlag = true;
				else ;
			}else ;

			if (caseFlag == true && ctrlFlag == true)
				break;
			else
				continue;
		}

		areUnique[j] = ((caseFlag == true && ctrlFlag == false) || (caseFlag == false && ctrlFlag == true));
	}

	//!- Collapsing of variants for unique loci only
	for (unsigned i = 0; i != __xdat.size(); ++i) {

		bool isWild = true;
		for (unsigned j = 0; j != __xdat[i].size(); ++j) {
			// scan all genetic loci of an individual

			// skip non-unique site
			if (areUnique[j] == false)
				continue;
			else ;

			//! - Apply indicator function: rare variant found. score it '1' and break the loop
			if (__xdat[i][j] != MAJOR_ALLELE && __xdat[i][j] != MISSING_ALLELE) {
				data.push_back(1.0);
				isWild = false;
				break;
			}else ;
		}

		//! - Otherwise score it '0' after every locus of an indv has been scanned.
		if (isWild == true)
			data.push_back(0.0);
	}

	return data;
}


double gwAssociations::m_checkAdaptivePvalue(unsigned permcount1, unsigned permcount2, unsigned currentIdx, unsigned checkPoint, unsigned sided) const
{
	if (currentIdx % checkPoint == 0 && checkPoint > 5) {
		//!- adaptive p-value calculation, at an interval of #checkPoint permutations
		// apply the "six-sigma" rule

		double pval = 1.0;
		if (sided == 1 || sided == 0)
			pval = (1.0 + permcount1) / (1.0 + currentIdx);
		else {
			double permcount = gw_dmin(permcount1, permcount2);
			pval = (2.0 * permcount + 2.0) / (1.0 + currentIdx);
		}

		double sd = sqrt(pval * (1.0 - pval) / (1.0 * currentIdx));
		double sixsigma = pval - 6.0 * sd;

		if (sixsigma > __alpha)
			return pval;
		else
			return 9.0;
	}else
		return 9.0;
}


void gwAssociations::m_printPunches(int n) const
{
	for (int i = 0; i != n; ++i)
		std::cout << "#";
	std::cout << "\n" << std::endl;
	return;
}



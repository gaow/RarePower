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


#include <cmath>
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

namespace gpow {

bool gwAssocdata::trimXdat()
{
	for (size_t j = 0; j != __observedMafs.size(); ++j) {
		if (__observedMafs[j] <= __mafLower || __observedMafs[j] > __mafUpper) {
			__observedMafs.erase(__observedMafs.begin() + j);
			for (size_t i = 0; i < __xdat.size(); ++i) {
				__xdat[i].erase(__xdat[i].begin() + j);
			}
			--j;
		}
	}
	return true;
}


bool gwAssocdata::markwildSibpairloci()
{
	for (unsigned i = 0; i < (__xdat.size() - 1); i += 2) {
		for (unsigned j = 0; j < __observedMafs.size(); ++j) {
			//if (__xdat[i][j] <= __xdat[i + 1][j]) {
				// case has less mutant allele than ctrl
			if (__xdat[i][j] == __xdat[i + 1][j]) {
				// case/ctrl are concordent
				__xdat[i][j] = 0.0;
				__xdat[i + 1][j] = 0.0;
			}
		}
	}
	return true;
}


vectorF gwAssocdata::countRegionalVariants() const
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


vectorF gwAssocdata::binariesRegionalVariants() const
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
			}
			;
		}

		//! - Otherwise score it '0' after every locus of an indv has been scanned.
		if (isWild == true)
			data.push_back(0.0);
	}

	return data;
}


vectorF gwAssocdata::binariesRegionalUniqueVariants() const
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
			}
			;

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

			//! - Apply indicator function: rare variant found. score it '1' and break the loop
			if (__xdat[i][j] != MAJOR_ALLELE && __xdat[i][j] != MISSING_ALLELE) {
				data.push_back(1.0);
				isWild = false;
				break;
			}
			;
		}

		//! - Otherwise score it '0' after every locus of an indv has been scanned.
		if (isWild == true)
			data.push_back(0.0);
	}

	return data;
}


//!\brief Collapsing variants simple LR score test statistic

double CmcstP::apply(gwAssocdata & d)
{
	/*! Li&Leal 2008 collapsing. Statistical test is the score test for simple logistic regression model.
	 */

	vectorF & ydat = d.ydat();

	//!- trim data by mafs upper-lower bounds
	d.trimXdat();

	//!- CMC scoring
	vectorF regressors = d.binariesRegionalVariants();

	if (__v)
		std::clog << regressors << std::endl;

	//!- logistic regression score statistic + permutation
	double xbar = gw_mean(regressors);
	unsigned nCases = 0;
	for (unsigned i = 0; i != ydat.size(); ++i)
		if (ydat[i] == AFFECTED)
			++nCases;

	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= __nperm) {
		double statistic = m_gstat.testLogitRegression1(regressors, ydat, xbar, nCases);
		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (__permcheckpnt)
				pvalue = checkP(permcount1, permcount2, iPermutation);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		d.permutate();
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		if (__sided == 1 || __sided == 0)
			pvalue = (1.0 + permcount1) / (1.0 + __nperm);
		else {
			double permcount = gw_dmin(permcount1, permcount2);
			pvalue = (2.0 * permcount + 2.0) / (1.0 + __nperm);
		}
	}
	return pvalue;
}


//!\brief Counts of variants simple LR score test statistic

double AnrvstP::apply(gwAssocdata & d)
{
	/*! Andrew P. Morris 2009 collapsing method. Statistical test is the score test for simple logistic regression model.
	 */
	vectorF & ydat = d.ydat();

	//!- trim data by mafs upper-lower bounds
	d.trimXdat();

	//!- ANRV scoring
	vectorF regressors = d.countRegionalVariants();
	if (__v)
		std::clog << regressors << std::endl;
	//!- logistic regression score statistic + permutation
	double xbar = gw_mean(regressors);
	unsigned nCases = 0;
	for (unsigned i = 0; i != ydat.size(); ++i)
		if (ydat[i] == AFFECTED)
			++nCases;


	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= __nperm) {
		double statistic = m_gstat.testLogitRegression1(regressors, ydat, xbar, nCases);
		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (__permcheckpnt)
				pvalue = checkP(permcount1, permcount2, iPermutation);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		d.permutate();
		++iPermutation;
	}

	if (__v)
		std::clog << permcount1 << " " << permcount2 << std::endl;

	if (pvalue <= 1.0) ;
	else {
		if (__sided == 1 || __sided == 0)
			pvalue = (1.0 + permcount1) / (1.0 + __nperm);
		else {
			double permcount = gw_dmin(permcount1, permcount2);
			pvalue = (2.0 * permcount + 2.0) / (1.0 + __nperm);
		}
	}
	return pvalue;
}


double CmcchiP::apply(gwAssocdata & d)
{
	/*!* By Li&Leal 2008 collapsing <br>
	 * * The original CMC test combine the collapsed rare variants and uncollapsed common variants in a Hotelling's T test. <br>
	 * * Here assume common variants are not analyzed -- the CMC test is now the "Collapsing test" <br>
	 * * Common variants not properly handled here but may be dealt with elsewhere (i.e. regression framework, etc) <br> <br>
	 * Implementation:
	 */
	__sided = 0;
	vectorF & ydat = d.ydat();

	//!- trim data by mafs upper-lower bounds
	d.trimXdat();

	//!- CMC scoring
	vectorF regressors = d.binariesRegionalVariants();
	if (__v)
		std::clog << regressors << std::endl;

	//!- 2 by 2 Chisq test

	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= __nperm) {
		double statistic = m_gstat.chisqtest2X2(regressors, ydat);
		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (__permcheckpnt)
				pvalue = checkP(permcount1, permcount2, iPermutation);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		d.permutate();
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		if (__sided == 1 || __sided == 0)
			pvalue = (1.0 + permcount1) / (1.0 + __nperm);
		else {
			double permcount = gw_dmin(permcount1, permcount2);
			pvalue = (2.0 * permcount + 2.0) / (1.0 + __nperm);
		}
	}
	return pvalue;
}


double CmcfisherP::apply(gwAssocdata & d)
{
	/*!* Similar with the chi-squared test based CMC but uses the Fisher's exact test, no permutation <br> <br>
	 * Implementation:
	 */
	vectorF & ydat = d.ydat();

	//!- trim data by mafs upper-lower bounds
	d.trimXdat();

	//!- CMC scoring
	vectorF regressors = d.binariesRegionalVariants();
	if (__v)
		std::clog << regressors << std::endl;

	//!- 2 by 2 Fisher's exact test, one or two sided
	double pvalue = m_gstat.fishertest2X2(regressors, ydat, __sided);
	return pvalue;
}


double RvefisherP::apply(gwAssocdata & d)
{
	/*!* Carrier frequencies are compared for variants that are only present in cases to those that are only observed controls <br> <br>
	 * Implementation:
	 */
	vectorF & ydat = d.ydat();

	//!- trim data by mafs upper-lower bounds
	d.trimXdat();

	//!- RVE scoring
	vectorF regressors = d.binariesRegionalUniqueVariants();
	if (__v)
		std::clog << regressors << std::endl;

	//!- 2 by 2 Fisher's exact test, one or two sided
	double pvalue = m_gstat.fishertest2X2(regressors, ydat, __sided);
	return pvalue;
}


double WssRankP::apply(gwAssocdata & d)
{
	/*! * Implement Browning 2009 Wss paper <br><br>
	 * Implementation:
	 */
	vectorF & ydat = d.ydat(); vector2F & xdat = d.xdat();

	//!- trim data by mafs upper-lower bounds
	d.trimXdat();
	unsigned nCases = 0;
	// case size

	for (unsigned i = 0; i != ydat.size(); ++i) {
		if (ydat[i] == AFFECTED)
			++nCases;
	}
	unsigned nCtrls = ydat.size() - nCases;

	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= __nperm) {

		//!- Calc MAF in controls
		vectorF weights(0);
		for (unsigned j = 0; j != xdat[0].size(); ++j) {

			unsigned nVariants = 0;
			//! - Count number of mutations in controls
			for (unsigned i = 0; i != xdat.size(); ++i) {
				// Control only
				if (ydat[i] == UNAFFECTED) {
					if (xdat[i][j] != MISSING_ALLELE)
						nVariants += (unsigned)xdat[i][j];
					else {
						std::cerr << "Input data problem in gwAssocdata::calcWssRankP(). Now Quit." << std::endl;
						exit(-1);
					}
				}else
					continue;
			}
			//! - Compute the "q" for the locus
			weights.push_back((nVariants + 1.0) / (2.0 * nCtrls + 2.0));
		}

		//!- Calc the Weights for each locus
		for (unsigned j = 0; j != xdat[0].size(); ++j)
			weights[j] = sqrt(2.0 * xdat.size() * weights[j] * (1.0 - weights[j]));

		vectorF scores(0);

		for (unsigned i = 0; i != xdat.size(); ++i) {
			double score = 0.0;
			//! - Define genetic score of the individual

			for (unsigned j = 0; j != xdat[i].size(); ++j) {
				// scan all loci of an individual
				double tmp = 0.0;
				// Dominant model
				if (__moi == 'D') {
					if (xdat[i][j] == HOMO_ALLELE)
						tmp = 1.0;
				}
				// recessive model
				else if (__moi == 'R') {
					if (xdat[i][j] != MAJOR_ALLELE)
						tmp = 1.0;
				}
				// additive and multiplicative models
				else ;

				//! - Parameterize mode of inheritance I_ij = {0, 1, 2} in Browning paper
				tmp = xdat[i][j] - tmp;
				score += tmp / weights[j];
			}

			//!- Calculate and store the genetic score
			scores.push_back(score);
		}

		if (__v)
			std::clog << scores << std::endl;

		//! - Compute rank test statistic, via <b> Mann_Whitneyu() </b>
		double caseScores[nCases], ctrlScores[nCtrls];
		int tmpa = 0, tmpu = 0;
		for (unsigned i = 0; i != ydat.size(); ++i) {
			if (ydat[i] == AFFECTED) {
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
			if (__permcheckpnt)
				pvalue = checkP(permcount1, permcount2, iPermutation);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		d.permutate();
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		if (__sided == 1 || __sided == 0)
			pvalue = (1.0 + permcount1) / (1.0 + __nperm);
		else {
			double permcount = gw_dmin(permcount1, permcount2);
			pvalue = (2.0 * permcount + 2.0) / (1.0 + __nperm);
		}
	}
	return pvalue;
}


double WssRankPA::apply(gwAssocdata & d)
{
	/*! * Implement Browning 2009 Wss paper normal approximation <br><br>
	 * Implementation:
	 */
	vectorF & ydat = d.ydat(); vector2F & xdat = d.xdat();

	//!- trim data by mafs upper-lower bounds
	d.trimXdat();
	unsigned nCases = 0;
	// case size

	for (unsigned i = 0; i != ydat.size(); ++i) {
		if (ydat[i] == AFFECTED)
			++nCases;
	}
	unsigned nCtrls = ydat.size() - nCases;

	unsigned iPermutation = 0;
	double observedStatistic = 0.0;
	vectorF statistics(0);

	while (iPermutation <= 1000) {

		//!- Calc MAF in controls
		vectorF weights(0);
		for (unsigned j = 0; j != xdat[0].size(); ++j) {

			int nVariants = 0;
			//! - Count number of mutations in controls
			for (unsigned i = 0; i != xdat.size(); ++i) {
				// Control only; consider additive codings only
				if (ydat[i] == UNAFFECTED) {
					if (xdat[i][j] != MISSING_ALLELE)
						nVariants += (int)xdat[i][j];
					else {
						std::cerr << "Input data problem in gwAssocdata::calcWssRankP. Now Quit." << std::endl;
						exit(-1);
					}
				}else
					continue;
			}
			//! - Compute the "q" for the locus
			weights.push_back((nVariants + 1.0) / (2.0 * nCtrls + 2.0));
		}

		//!- Calc the Weights for each locus
		for (unsigned j = 0; j != xdat[0].size(); ++j)
			weights[j] = sqrt(2.0 * xdat.size() * weights[j] * (1.0 - weights[j]));

		vectorF scores(0);

		for (unsigned i = 0; i != xdat.size(); ++i) {
			double score = 0.0;
			//! - Define genetic score of the individual

			for (unsigned j = 0; j != xdat[i].size(); ++j) {
				// scan all loci of an individual
				double tmp = 0.0;
				// Dominant model
				if (__moi == 'D') {
					if (xdat[i][j] == HOMO_ALLELE)
						tmp = 1.0;
				}
				// recessive model
				else if (__moi == 'R') {
					if (xdat[i][j] != MAJOR_ALLELE)
						tmp = 1.0;
				}
				// additive and multiplicative models
				else ;

				//! - Parameterize mode of inheritance I_ij = {0, 1, 2} in Browning paper
				tmp = xdat[i][j] - tmp;
				score += tmp / weights[j];
			}

			//!- Calculate and store the genetic score

			scores.push_back(score);
		}

		if (__v)
			std::clog << scores << std::endl;

		//! - Compute rank test statistic, via <b> Mann_Whitneyu() </b>
		double caseScores[nCases], ctrlScores[nCtrls];
		int tmpa = 0, tmpu = 0;
		for (unsigned i = 0; i != ydat.size(); ++i) {
			if (ydat[i] == AFFECTED) {
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
		d.permutate();
		++iPermutation;
	}

	double mean = gw_mean(statistics);
	double sd = gw_sd(statistics);
	if (sd == 0.0)
		sd = 1.0e-6;

	double statisticstd = (observedStatistic - mean) / sd;
	double pvalue = 9.0;
	if (__sided == 1) {
		pvalue = gsl_cdf_ugaussian_Q(statisticstd);
	}else if (__sided == 2) {
		statisticstd = statisticstd * statisticstd;
		pvalue = gsl_cdf_chisq_Q(statisticstd, 1.0);
	}else {
		std::cerr << "ERROR: \"__sided\" should be either 2 (default two-sided test) or 1 (one-sided test)" << std::endl;
		exit(-1);
	}
	return pvalue;
}


double KbacP::apply(gwAssocdata & d)
{
	/*! * the KBAC Statistic: sum of genotype pattern frequencies differences, weighted by hypergeometric kernel. <br>
	 * * It is a permutation based one/two-sided test.  <br>
	 * * See <em> Liu DJ 2010 PLoS Genet. </em> <br><br>
	 * * Implementation:
	 */
	unsigned model = 1;

	if (__sided == 2) {
		__sided = 0;
		model = 2;
	}

	vectorF & ydat = d.ydat(); vector2F & xdat = d.xdat();

	//!- trim data by mafs upper-lower bounds
	d.trimXdat();


	unsigned sampleSize = ydat.size();
	// sample size
	unsigned regionLen = xdat[0].size();
	// candidate region length
	unsigned nCases = 0;
	// case size

	for (unsigned i = 0; i != ydat.size(); ++i) {
		if (ydat[i] == AFFECTED)
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

			if (xdat[i][j] != MISSING_ALLELE && xdat[i][j] != MAJOR_ALLELE)
				vntIdR += pow(3.0, 1.0 * (j - lastCnt)) * xdat[i][j];
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

	//  std::clog << uniquePattern << std::endl;

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
	while (iPermutation <= __nperm) {

		// the KBAC statistic. Will be of length 1 or 2
		vectorF kbacStatistics(0);
		// two models
		for (unsigned s = 0; s != model; ++s) {

			//!- count number of sample cases (for the 1st model, or ctrls for the 2nd model) for each genotype pattern
			unsigned uniquePatternCountsSub[uniquePattern.size()];
			for (unsigned u = 0; u != uniquePattern.size(); ++u)
				uniquePatternCountsSub[u] = 0;
			// genotype pattern counts in cases (for the 1st model, or ctrls for the 2nd model)

			for (unsigned i = 0; i != sampleSize; ++i) {
				if (ydat[i] == (AFFECTED - 1.0 * s) ) {
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

			if (__v)
				std::clog << kbac << std::endl;

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
			if (__permcheckpnt)
				pvalue = checkP(permcount1, permcount2, iPermutation);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		d.permutate();
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		pvalue = (1.0 + permcount1) / (1.0 + __nperm);
	}
	return pvalue;
}


double KbacstP::apply(gwAssocdata & d)
{
	/*! * the KBAC in score test
	   *Implementation:
	 */

	unsigned model = 1;

	if (__sided == 2) {
		__sided = 0;
		model = 2;
	}

	vectorF & ydat = d.ydat(); vector2F & xdat = d.xdat();

	//!- trim data by mafs upper-lower bounds
	d.trimXdat();


	unsigned sampleSize = ydat.size();
	// sample size
	unsigned regionLen = xdat[0].size();
	// candidate region length
	unsigned nCases = 0;
	// case size

	for (unsigned i = 0; i != ydat.size(); ++i) {
		if (ydat[i] == AFFECTED)
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

			if (xdat[i][j] != MISSING_ALLELE && xdat[i][j] != MAJOR_ALLELE)
				vntIdR += pow(3.0, 1.0 * (j - lastCnt)) * xdat[i][j];
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

	//  std::clog << uniquePattern << std::endl;

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
	while (iPermutation <= __nperm) {

		// the KBAC statistic. Will be of length 1 or 2
		vectorF kbacStatistics(0);
		// two models
		for (unsigned s = 0; s != model; ++s) {

			//!- count number of sample cases (for the 1st model, or ctrls for the 2nd model) for each genotype pattern
			unsigned uniquePatternCountsSub[uniquePattern.size()];
			for (unsigned u = 0; u != uniquePattern.size(); ++u)
				uniquePatternCountsSub[u] = 0;
			// genotype pattern counts in cases (for the 1st model, or ctrls for the 2nd model)

			for (unsigned i = 0; i != sampleSize; ++i) {
				if (ydat[i] == (AFFECTED - 1.0 * s) ) {
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
			double kbac = m_gstat.testLogitRegression1(regressors, ydat, xbar, nCases);
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
			if (__permcheckpnt)
				pvalue = checkP(permcount1, permcount2, iPermutation);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		d.permutate();
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		pvalue = (1.0 + permcount1) / (1.0 + __nperm);
	}
	return pvalue;
}


double VtP::apply(gwAssocdata & d)
{
	/*! * Variable threshold method, Price 2010 AJHG <br>
	 * * <em> Unique feature of VT </em>: instead of using fixed threshold + Morris A RV counts, it uses a variable threshold for RV frequency. <br>
	 * * Since it uses Z_max = max(Z_i ...) where Z_i is proportional to standard normal by a coefficient which is not a function of RV counts, it does not matter to include the "known" common variants. <br>
	 * * In brief it is a multiple testing method that analyzes both common and rare variants <br><br>
	 * Implementation:
	 */
	if (__sided == 2) __sided = 0;
	vectorF & ydat = d.ydat(); vector2F & xdat = d.xdat();

	// case size
	unsigned nCases = 0;

	for (unsigned i = 0; i != ydat.size(); ++i) {
		if (ydat[i] == AFFECTED)
			++nCases;
	}

	double po = (1.0 * nCases) / (1.0 * ydat.size());

	// number of variants at each loci
	vectorUI nLocusVariants(xdat[0].size(), 0);
	for (unsigned j = 0; j != xdat[0].size(); ++j) {
		for (unsigned i = 0; i != xdat.size(); ++i) {
			nLocusVariants[j] += (unsigned)xdat[i][j];
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

	while (iPermutation <= __nperm) {

		//! - Define <b> 'allZs' </b>, a vector of the Z scores computed under different thresholds
		vectorF allZs(0);
		vectorF zIa(xdat.size(), 0.0), zIb(xdat.size(), 0.0);

		//! - Iterate the following for all thresholds:
		for (unsigned t = 1; t != vts.size(); ++t) {

			//! - - Record the index of loci that can be added at the new threshold criteria
			vectorUI idxesAdding(0);
			for (unsigned s = 0; s != nLocusVariants.size(); ++s) {
				if (nLocusVariants[s] > vts[t - 1] && nLocusVariants[s] <= vts[t])
					idxesAdding.push_back(s);
			}

			//!- - For loci passing the threshold criteria, implement Price paper page 3 z(T) formula
			for (unsigned i = 0; i != xdat.size(); ++i) {

				double zIai = 0.0, zIbi = 0.0;
				for (unsigned j = 0; j != idxesAdding.size(); ++j) {
					unsigned locIdx = idxesAdding[j];
					zIai += ((ydat[i] - 1.0) - po) * xdat[i][locIdx];
					zIbi += xdat[i][locIdx] * xdat[i][locIdx];
				}

				zIa[i] += zIai;
				zIb[i] += zIbi;
			}

			//! - Now each element in <b> 'allZs' </b> is z(T) for different thresholds
			allZs.push_back(gw_sum(zIa) / sqrt(gw_sum(zIb) ) );
		}

		//! - Compute zmax, the statistic; square it so that it would work for protectives
		if (__sided == 0) {
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
			if (__permcheckpnt)
				pvalue = checkP(permcount1, permcount2, iPermutation);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		d.permutate();
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		pvalue = (1.0 + permcount1) / (1.0 + __nperm);
	}
	return pvalue;
}


double VtFisherP::apply(gwAssocdata & d)
{
	/*! * Variable threshold method + Fisher's test <br><br>
	 * The test combines VT and CMC one-sided Fisher's with mid-p corrections. If none of the tests are significant (judged by Fisher's exact p-value) for the original data then no permutation test is pursued <br><br>
	 * Implementation:
	 */
	if (__sided == 2) __sided = 0;
	vectorF & ydat = d.ydat(); vector2F & xdat = d.xdat();

	// number of variants at each loci
	vectorUI nLocusVariants(xdat[0].size(), 0);

	for (unsigned j = 0; j != xdat[0].size(); ++j) {
		for (unsigned i = 0; i != xdat.size(); ++i) {
			nLocusVariants[j] += (unsigned)xdat[i][j];
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

	while (iPermutation <= __nperm) {

		//! - Define <b> 'allPs' </b>, a vector of the p-values computed under different thresholds
		vectorF allPs(0);
		vectorF regressorsCurr(xdat.size(), 0.0);
		//! - Iterate the following for all thresholds:
		for (unsigned t = 1; t != vts.size(); ++t) {

			//! - - Record the index of loci that can be added at the new threshold criteria
			vectorUI idxesAdding(0);
			for (unsigned s = 0; s != nLocusVariants.size(); ++s) {
				if (nLocusVariants[s] > vts[t - 1] && nLocusVariants[s] <= vts[t])
					idxesAdding.push_back(s);
			}

			//!- - For loci passing the threshold, implement CMC one-side Fisher
			vectorF regressors(xdat.size(), 0.0);

			for (unsigned i = 0; i != regressors.size(); ++i) {
				for (unsigned j = 0; j != idxesAdding.size(); ++j) {
					unsigned locIdx = idxesAdding[j];
					regressors[i] = (xdat[i][locIdx] == MAJOR_ALLELE || xdat[i][locIdx] == MISSING_ALLELE) ? regressors[i] : 1.0;
				}
				regressorsCurr[i] = (regressors[i] == 1.0) ? 1.0 : regressorsCurr[i];
			}

			double pfisher = m_gstat.fishertest2X2(regressorsCurr, ydat, (__sided) ? 1 : 2);
			if (pfisher <= __alpha) {
				shouldPermutationTest = true;
			}

			allPs.push_back(-log(pfisher));
		}

		//!- an 999999 approach via Fisher's test
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
			if (__permcheckpnt)
				pvalue = checkP(permcount1, permcount2, iPermutation);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		d.permutate();
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		pvalue = (1.0 + permcount1) / (1.0 + __nperm);
	}
	return pvalue;
}


double AsumP::apply(gwAssocdata & d)
{
	/*! * Number of rare variants per site with "protective" site recoded, by Pan and Han (2010) Hum Hered <br>
	 *  * The authors use alpha0 = 0.1 to screen variants that should be recoded. See their paper for details <br>
	 *  * Original paper uses LR score test with approx. distribution, which is equivalent to Armitage's test for trend for this case (Schaid 2002). For RV not many sites are homozygote. In my implementation I use Fisher's 2by2 test instead.
	 * * Here allows user specified lower and upper bounds <br><br>
	 * Implementation:
	 */
	__sided = 0;
	vectorF & ydat = d.ydat(); vector2F & xdat = d.xdat();

	//!- trim data by mafs upper-lower bounds
	d.trimXdat();
	unsigned nCases = 0;
	for (unsigned i = 0; i != ydat.size(); ++i)
		if (ydat[i] == AFFECTED)
			++nCases;
	unsigned nCtrls = ydat.size() - nCases;

	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= __nperm) {

		//!- Sites that have excess rare variants in controls
		vectorL ctrlExcess(xdat[0].size(), false);

		for (unsigned j = 0; j != xdat[0].size(); ++j) {
			double tmpCase = 0.0;
			double tmpCtrl = 0.0;
			for (unsigned i = 0; i != xdat.size(); ++i) {
				if (ydat[i] == UNAFFECTED) {
					if (xdat[i][j] != MISSING_ALLELE && xdat[i][j] != MAJOR_ALLELE)
						tmpCtrl += xdat[i][j];
				}else {
					if (xdat[i][j] != MISSING_ALLELE && xdat[i][j] != MAJOR_ALLELE)
						tmpCase += xdat[i][j];
				}
			}
			if (tmpCtrl / (nCtrls * 2.0) > tmpCase / (nCases * 2.0))
				ctrlExcess[j] = true;
			else
				continue;
		}

		//!- Define sites that needs to be recoded
		vectorL recodeSites(xdat[0].size(), false);

		for (unsigned j = 0; j != xdat[0].size(); ++j) {
			if (ctrlExcess[j] == false)
				continue;

			vectorF vdat(0);
			for (unsigned i = 0; i != xdat.size(); ++i)
				vdat.push_back(xdat[i][j]);

			//!- 2 by 2 Fisher's test
			double pfisher = m_gstat.fishertest2X2(vdat, ydat, 2, __moi);

			if (pfisher < 0.1)
				recodeSites[j] = true;
			else
				continue;
		}

		vectorF regressors(0);

		//! - Count anrv; recode when necessary
		for (unsigned i = 0; i != xdat.size(); ++i) {
			//  scan all sample individuals

			double nrv = 0.0;
			for (unsigned j = 0; j != xdat[0].size(); ++j) {
				double tmpCode = 0.0;
				if (xdat[i][j] != MISSING_ALLELE)
					tmpCode += xdat[i][j];
				//!- Recode protective variants as in Pan 2010
				if (recodeSites[j] == true)
					tmpCode = 2.0 - tmpCode;
				// aggregated number of rare variants
				nrv += tmpCode;
			}

			regressors.push_back(nrv);
		}

		//!- Score test implementation for logistic regression model logit(p) = b0 + b1x (derivation of the test see my labnotes vol.2 page 3)

		if (__v)
			std::clog << regressors << std::endl;
		//!- logistic regression score statistic
		double xbar = gw_mean(regressors);
		double statistic = m_gstat.testLogitRegression1(regressors, ydat, xbar, nCases);

		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (__permcheckpnt)
				pvalue = checkP(permcount1, permcount2, iPermutation);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		d.permutate();
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		pvalue = (1.0 + permcount1) / (1.0 + __nperm);
	}
	return pvalue;
}


double CmcqtP::apply(gwAssocdata & d)
{
	vectorF & ydat = d.ydat();

	//!- trim data by mafs upper-lower bounds
	d.trimXdat();

	//!- CMC scoring
	vectorF regressors = d.binariesRegionalVariants();

	if (__v)
		std::clog << regressors << std::endl;


	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= __nperm) {

		//!- two sample t test. Group 1: have cmc score = 1, Group 2: have cmc score = 0
		double statistic = m_gstat.ttestIndp(regressors, ydat);
		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (__permcheckpnt)
				pvalue = checkP(permcount1, permcount2, iPermutation);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		d.permutate();
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		if (__sided == 1 || __sided == 0)
			pvalue = (1.0 + permcount1) / (1.0 + __nperm);
		else {
			double permcount = gw_dmin(permcount1, permcount2);
			pvalue = (2.0 * permcount + 2.0) / (1.0 + __nperm);
		}
	}
	return pvalue;
}


double AnrvqtPermP::apply(gwAssocdata & d)
{

	vectorF & ydat = d.ydat();

	//!- trim data by mafs upper-lower bounds
	d.trimXdat();
	//!- ANRV scoring
	vectorF regressors = d.countRegionalVariants();
	if (__v)
		std::clog << regressors << std::endl;

	//!- Linear regression test
	double xbar = gw_mean(regressors);
	double ybar = gw_mean(ydat);

	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= __nperm) {
		double statistic = m_gstat.testLnRegression1(regressors, ydat, xbar, ybar);
		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (__permcheckpnt)
				pvalue = checkP(permcount1, permcount2, iPermutation);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		d.permutate();
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		if (__sided == 1 || __sided == 0)
			pvalue = (1.0 + permcount1) / (1.0 + __nperm);
		else {
			double permcount = gw_dmin(permcount1, permcount2);
			pvalue = (2.0 * permcount + 2.0) / (1.0 + __nperm);
		}
	}
	return pvalue;
}


double AnrvqtP::apply(gwAssocdata & d)
{
	vectorF & ydat = d.ydat();

	//!- trim data by mafs upper-lower bounds
	d.trimXdat();
	//!- ANRV scoring
	vectorF regressors = d.countRegionalVariants();
	if (__v)
		std::clog << regressors << std::endl;

	//!- Linear regression test
	double xbar = gw_mean(regressors);
	double ybar = gw_mean(ydat);
	double pvalue = m_gstat.testLnRegression1(regressors, ydat, xbar, ybar);
	if (__sided == 1 || __sided == 0) {
		pvalue = gsl_cdf_ugaussian_Q(pvalue);
	}else {
		pvalue = pvalue * pvalue;
		pvalue = gsl_cdf_chisq_Q(pvalue, 1.0);
	}
	return pvalue;
}


double AnrvqtCondP::apply(gwAssocdata & d)
{
	vectorF & ydat = d.ydat();

	//!- trim data by mafs upper-lower bounds
	d.trimXdat();
	//!- ANRV scoring
	vectorF regressors = d.countRegionalVariants();
	if (__v)
		std::clog << regressors << std::endl;

	double yh = getDoubleVar("yh");
	double yl = getDoubleVar("yl");
	assert(yh >= yl);
	//!- Linear regression conditional score test
	double pvalue = m_gstat.testLnRegressionCond(regressors, ydat, yh, yl);
	if (__sided == 1 || __sided == 0) {
		pvalue = gsl_cdf_ugaussian_Q(pvalue);
	}else {
		pvalue = pvalue * pvalue;
		pvalue = gsl_cdf_chisq_Q(pvalue, 1.0);
	}
	return pvalue;
}


double TestRareP::apply(gwAssocdata & d)
{
	/*! * the RVP Statistic: variant counts weighted by poisson kernels. <br>
	 * * It is a permutation based two model test, max(protective model, risk model).  <br>
	 * * See <em> Ionita-Laza and Lange 2011 PLoS Genet. </em> <br><br>
	 * * Implementation (credit to the authors for part of the codes below):
	 */
	if (__sided == 2) __sided = 0;
	vectorF & ydat = d.ydat(); vector2F & xdat = d.xdat();

	//!- trim data by mafs upper-lower bounds
	d.trimXdat();

	unsigned sampleSize = xdat.size();
	// sample size
	unsigned regionLen = xdat[0].size();
	// candidate region length

	unsigned nCases = 0;
	for (unsigned i = 0; i != ydat.size(); ++i) {
		if (ydat[i] == AFFECTED)
			++nCases;
	}
	unsigned nCtrls = sampleSize - nCases;


	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= __nperm) {

		double sumR = 0, sumP = 0;
		for (unsigned j = 0; j != regionLen; ++j) {

			//! - Count number of variants in cases/controls at a locus
			unsigned countcs = 0;
			unsigned countcn = 0;
			for (unsigned i = 0; i != sampleSize; ++i) {
				if (ydat[i] == UNAFFECTED)
					countcn += (unsigned)xdat[i][j];
				else
					countcs += (unsigned)xdat[i][j];
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
		double statistic = (__sided == 1) ? sumR : fmax(sumR, sumP);

		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (__permcheckpnt)
				pvalue = checkP(permcount1, permcount2, iPermutation);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		d.permutate();
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		pvalue = (1.0 + permcount1) / (1.0 + __nperm);
	}
	return pvalue;
}


double CalphaP::apply(gwAssocdata & d)
{

	/*! * the c-alpha Statistic: sum of the std. error of variant counts in cases <br>
	 * *  Two-sided test, claims to be robust against protective mutations. <br>
	 * * See <em> Ben. Neale et al. 2011 PLoS Genet. </em> <br><br>
	 * * Implementation:
	 */
	vectorF & ydat = d.ydat(); vector2F & xdat = d.xdat();

	//!- trim data by mafs upper-lower bounds
	d.trimXdat();

	unsigned sampleSize = xdat.size();
	// sample size
	unsigned regionLen = xdat[0].size();
	// candidate region length

	unsigned nCases = 0;
	for (unsigned i = 0; i != ydat.size(); ++i) {
		if (ydat[i] == AFFECTED)
			++nCases;
	}

	double phat = (1.0 * nCases) / (1.0 * sampleSize);


	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= __nperm) {

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
				if (ydat[i] == UNAFFECTED)
					countcn += (unsigned)xdat[i][j];
				else
					countcs += (unsigned)xdat[i][j];
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

		//!- c-alpha statistic
		double statistic = calpT / sqrt(calpV);
		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (__permcheckpnt)
				pvalue = checkP(permcount1, permcount2, iPermutation);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		d.permutate();
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		if (__sided == 1 || __sided == 0)
			pvalue = (1.0 + permcount1) / (1.0 + __nperm);
		else {
			double permcount = gw_dmin(permcount1, permcount2);
			pvalue = (2.0 * permcount + 2.0) / (1.0 + __nperm);
		}
	}
	return pvalue;
}


double RareCoverP::apply(gwAssocdata & d)
{
	/*! * RareCover method, 2010 PLoS CompBio <br>
	 * * Implementation:
	 */
	if (__sided == 2) __sided = 0;
	vectorF & ydat = d.ydat(); vector2F & xdat = d.xdat();

	//!- the cut-off to use for the "heuristic greedy algorithm". = 0.5 as suggested by the paper
	const double difQ = 0.5;

	//!- Index of loci that are observed to be polymophic
	vectorUI vntVct(0);
	vectorF & mafs = d.mafs();

	for (unsigned j = 0; j != mafs.size(); ++j) {
		if (mafs[j] > 0.0)
			vntVct.push_back(j + 1);
	}

	if (vntVct.size() == 0)
		return 1.0;

	unsigned smpSize = xdat.size();


	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= __nperm) {

		//!- the current and the next statistics
		double sCurr = 0.0, sNext = 0.0;
		vectorF regressors(ydat.size(), 0.0);
		//!- the current and the next genotype coding
		vectorF regressorsCurr(ydat.size(), 0.0);
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
					regressors[i] = (regressorsCurr[i] + xdat[i][iIdx] > 0) ? MINOR_ALLELE : MAJOR_ALLELE;

				double statistic;

				//! - 2 by 2 one-sided Fisher's test
				if (__sided == 1) {
					statistic = -1.0 * log(m_gstat.fishertest2X2(regressors, ydat, __sided));
				}
				//! - 2 by 2 Chisq test
				else {
					statistic = m_gstat.chisqtest2X2(regressors, ydat);
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
					regressorsCurr[i] = regressorsCurr[i] + xdat[i][rmVnt];

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
			if (__permcheckpnt)
				pvalue = checkP(permcount1, permcount2, iPermutation);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		d.permutate();
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		if (__sided == 1 || __sided == 0)
			pvalue = (1.0 + permcount1) / (1.0 + __nperm);
		else {
			double permcount = gw_dmin(permcount1, permcount2);
			pvalue = (2.0 * permcount + 2.0) / (1.0 + __nperm);
		}
	}
	return pvalue;
}


double WsFisherP::apply(gwAssocdata & d)
{
	/*! * Weighted sum Fisher test, Wang Shuang and Patrick (2011) <br>
	 * Implementation:
	 */
	if (__sided == 2) __sided = 0;
	vectorF & ydat = d.ydat(); vector2F & xdat = d.xdat();

	//!- trim data by mafs upper-lower bounds
	d.trimXdat();

	unsigned nAllCases = 0;
	for (unsigned i = 0; i != ydat.size(); ++i)
		if (ydat[i] == AFFECTED)
			++nAllCases;
	unsigned nAllCtrls = ydat.size() - nAllCases;


	unsigned regionLen = xdat[0].size();

	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= __nperm) {

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

			for (unsigned i = 0; i != xdat.size(); ++i) {

				if (xdat[i][j] == MISSING_ALLELE) {
					if (ydat[i] == UNAFFECTED) --nCtrls;
					else --nCases;
					continue;
				}

				switch (__moi) {
				case 'R':
				{
					if (ydat[i] == UNAFFECTED)
						allelesCtrl += (xdat[i][j] == HOMO_ALLELE) ? MINOR_ALLELE : MAJOR_ALLELE;
					else
						allelesCase += (xdat[i][j] == HOMO_ALLELE) ? MINOR_ALLELE : MAJOR_ALLELE;
				}
				break;
				case 'D':
				{
					if (ydat[i] == UNAFFECTED)
						allelesCtrl += (xdat[i][j] != MAJOR_ALLELE) ? MINOR_ALLELE : MAJOR_ALLELE;
					else
						allelesCase += (xdat[i][j] != MAJOR_ALLELE) ? MINOR_ALLELE : MAJOR_ALLELE;
				}
				break;
				default:
				{
					multiplier = 2.0;
					if (ydat[i] == UNAFFECTED)
						allelesCtrl += xdat[i][j];
					else
						allelesCase += xdat[i][j];
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
				    (__isMidP)
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
				    (__isMidP)
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
		if (__sided == 1)
			statistic = statistics[0];
		else statistic = fmax(statistics[0], statistics[1]);

		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (__permcheckpnt)
				pvalue = checkP(permcount1, permcount2, iPermutation);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		d.permutate();
		++iPermutation;
	}
	if (pvalue <= 1.0) ;
	else {
		pvalue = (1.0 + permcount1) / (1.0 + __nperm);
	}
	return pvalue;
}


double SkatP::apply(gwAssocdata & d)
{

#ifdef SKAT_GSL
	gsl_set_error_handler_off();
#endif
	vectorF & ydat = d.ydat(); vector2F & xdat = d.xdat();
	__sided = 0;
	//!- trim data by mafs upper-lower bounds
	d.trimXdat();
	// case size
	unsigned nCases = 0;
	for (unsigned i = 0; i != ydat.size(); ++i) {
		if (ydat[i] == AFFECTED)
			++nCases;
	}
	double po = (1.0 * nCases) / (1.0 * ydat.size());
	//
	// observed maf and weights
	//
	vectorF obsMafs(xdat[0].size(), 0.0);
	for (unsigned j = 0; j != xdat[0].size(); ++j) {
		for (unsigned i = 0; i != xdat.size(); ++i) {
			obsMafs[j] += xdat[i][j] / (ydat.size() * 2.0);
		}
	}
	vectorF weights(xdat[0].size(), 0.0);
	for (unsigned j = 0; j != xdat[0].size(); ++j) {
		weights[j] = pow(gsl_ran_beta_pdf(obsMafs[j], 1.0, 25.0), 2.0);
	}
	//
	// GWG'
	//
	vector2F gmat = xdat;
	for (unsigned i = 0; i != gmat.size(); ++i) {
		std::transform(gmat[i].begin(), gmat[i].end(), weights.begin(), gmat[i].begin(), std::multiplies<double>());
	}

#ifdef SKAT_GSL
	gsl_matrix * kmat = gsl_matrix_alloc(xdat.size(), xdat.size());

	for (unsigned i = 0; i != xdat.size(); ++i) {
		for (unsigned k = 0; k != xdat.size(); ++k) {
			double gmatii = 0.0;
			for (unsigned j = 0; j != gmat[i].size(); ++j) {
				gmatii += gmat[i][j] * xdat[k][j];
			}
			gsl_matrix_set(kmat, i, k, gmatii);
		}
	}
#else
	vector2F kmat(gmat.size());
	for (unsigned i = 0; i != kmat.size(); ++i) {
		for (unsigned k = 0; k != xdat.size(); ++k) {
			double gmatii = 0.0;
			for (unsigned j = 0; j != gmat[i].size(); ++j) {
				gmatii += gmat[i][j] * xdat[k][j];
			}
			kmat[i].push_back(gmatii);
		}
	}
#endif
	for (size_t i = 0; i < ydat.size(); ++i) {
		ydat[i] -= (1.0 + po);
	}

	unsigned iPermutation = 0;
	unsigned permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	double pvalue = 9.0;

	while (iPermutation <= __nperm) {
		//
		// Q
		//
#ifdef SKAT_GSL
		gsl_vector * qvec = gsl_vector_alloc(ydat.size());
		gsl_vector_view uvec = gsl_vector_view_array(&ydat[0], ydat.size());
		int m_err = gsl_blas_dgemv(CblasTrans, 1.0, kmat, &uvec.vector, 0.0, qvec);
		if (m_err != 0) {
			std::cerr << "Error in gsl_blas_dgemv(CblasTrans, 1.0, x, y, 0.0, b)" << std::endl;
			exit(-1);
		}

		gsl_vector_mul(qvec, &uvec.vector);
		double statistic = 0.0;
		for (size_t i = 0; i < ydat.size(); ++i) {
			statistic += gsl_vector_get(qvec, i);
		}
		gsl_vector_free(qvec);
#else
		vectorF qvec(0);
		for (unsigned j = 0; j != kmat[0].size(); ++j) {
			double tmp = 0.0;
			for (unsigned i = 0; i != ydat.size(); ++i) {
				tmp += ydat[i] * kmat[i][j];
			}
			qvec.push_back(tmp);
		}

		std::transform(qvec.begin(), qvec.end(), ydat.begin(), qvec.begin(), std::multiplies<double>());
		double statistic = std::accumulate(qvec.begin(), qvec.end(), 0.0, std::plus<double>());
#endif
		if (iPermutation == 0)
			observedStatistic = statistic;
		else {
			if (statistic >= observedStatistic)
				++permcount1;
			if (statistic <= observedStatistic)
				++permcount2;
			if (__permcheckpnt)
				pvalue = checkP(permcount1, permcount2, iPermutation);
		}
		if (pvalue <= 1.0)
			break;
		//!- Permutation
		d.permutate();
		++iPermutation;
	}

	if (pvalue <= 1.0) ;
	else {
		if (__sided == 1 || __sided == 0)
			pvalue = (1.0 + permcount1) / (1.0 + __nperm);
		else {
			double permcount = gw_dmin(permcount1, permcount2);
			pvalue = (2.0 * permcount + 2.0) / (1.0 + __nperm);
		}
	}
#ifdef SKAT_GSL
	gsl_matrix_free(kmat);
#endif
	return pvalue;
}


double gwBaseTest::checkP(unsigned pcount1, unsigned pcount2, size_t current) const
{
	//!- adaptive p-value calculation at an interval of #checkPoint permutations
	if (current % __permcheckpnt != 0 || __permcheckpnt < 100) {
		return 9.0;
	}
	double x;
	if (__sided == 1 || __sided == 0) {
		x = 1.0 + pcount1;
	} else {
		x = fmin(pcount1 + 1.0, pcount2 + 1.0);
	}


	double n = current + 1.0;
	double alpha = 0.05;
	double pval = x / n;
	double z = gsl_cdf_gaussian_Pinv(1.0 - alpha / 2.0, 1.0);
	double plw = pval - z * sqrt(pval * (1.0 - pval) / n);
	//OPTION2: six-sigma rule
	//double sd = sqrt(pval * (1.0 - pval) / (1.0 * current));
	//double plw = pval - 6.0 * sd;

	plw = (__sided == 1) ? plw : plw * 2.0;
	if (plw > __alpha) {
		return (__sided == 1) ? pval : pval * 2.0;
	} else {
		return 9.0;
	}
}


}

// =====================================================================================
//
//       Filename:  assoctests.h
//
//    Description:  Statistical tests for association mapping
//
//        Version:  1.0
//        Created:  01/04/2011 10:46:12 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Gao Wang (gw), wangow@gmail.com
//                  Baylor College of Medicine, Texas, USA
//        License:  GNU General Public License <http://www.gnu.org/licenses/>
//                  Copyright (c) 2011, Gao Wang
//
// =====================================================================================


#ifndef ASSOC_TESTS_H
#define ASSOC_TESTS_H

#include <vector>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"
#include "gw_utilities.h"
#include "gw_maths.h"

/*!\brief Rare variants association tests. Notes:
 * <br>
 * The testRare() method was adopted from the author's original code.
 *  Credit goes to Dr. Iuliana Ionita-Laza at Columbia U.
 */

class gwAssociations
{

public:
	/*!\brief Constructor, verbose
	 * \param observedMafs Observed minor allele frequencies in the data-set
	 * \param ydat
	 * \param xdat
	 * \param mafLower Lower bound of MAF that shall be analyzed. Many rare variants methods have it as (0, 0.01]
	 * \param mafUpper Upper bound
	 * \param alpha pre-set (desired) significance level
	 */
	gwAssociations(const vectorF & observedMafs, const vectorF & ydat,
		const vector2F & xdat, double mafLower, double mafUpper, double alpha);


	/*!\brief Deconstrutor, nothing to do
	 */
	~gwAssociations();


	/*!\brief Verbose screen out-put, for debug mainly
	 */
	void setVerbose(int debuglevel);

	void debug(int showWhat) const;


	//!\brief Collapsed variants simple LR model score test statistic
	double calcCmcstP(unsigned sided, unsigned nPermutations, unsigned adaptive);


	//!\brief Counts of variants simple LR model score test statistic
	double calcAnrvstP(unsigned sided, unsigned nPermutations, unsigned adaptive);


	//!\brief Collapsed variants 2 by 2 chisq statistic
	double calcCmcchiP(unsigned sided, unsigned nPermutations, unsigned adaptive);


	//!\brief CMC Fisher's exact test p-value no permutation
	double calcCmcfisherP(unsigned sided);


	//!\brief Cohen's RVE fisher exact test p-value no permutation
	double calcRvefisherP(unsigned sided);


	//!\brief The WSS method rank statistic, page 3 of Browning paper.
	//!- Not using normal approximation
	double calcWssRankP(const char moi, unsigned nCtrls, unsigned sided, unsigned nPermutations, unsigned adaptive);

	//!- Using normal approximation
	double calcWssRankP(const char moi, unsigned nCtrls, unsigned sided);


	//!\brief The KBAC method, Dajiang Liu's paper
	double calcKbacP(unsigned sided, unsigned nPermutations, unsigned adaptive);

	double calcKbacstP(unsigned sided, unsigned nPermutations, unsigned adaptive);


	//!\brief The variable threshold method, Price's paper 2010 AJHG
	double calcVtP(unsigned sided, unsigned nPermutations, unsigned adaptive);

	//!\brief The variable threshold method with Fisher's test
	double calcVtFisherP(unsigned sided, unsigned nPermutations, unsigned adaptive);

	//!\brief The adaptive sum test (in Pan and Han 2010) + permutation
	double calcAsumP(const char moi, unsigned sided, unsigned nPermutations, unsigned adaptive);


	//!\brief CMCQT
	double calcCmcqtP(unsigned sided, unsigned nPermutations, unsigned adaptive);


	//!\brief ANRVQT
	double calcAnrvqtP(unsigned sided);

	//!\brief ANRVQT for Xtreme sampling
	double calcAnrvqtP(double yh, double yl, unsigned sided);

	//!\brief ANRVQT
	double calcAnrvqtPermP(unsigned sided, unsigned nPermutations, unsigned adaptive);


	//!\brief Ionita-Laza and Lange 2011 PLoS Genet.
	double calcTestRareP(unsigned sided, unsigned nPermutations, unsigned adaptive);


	//!\brief c-alpha test, Ben. Neale et al. 2011 PLoS Genet.
	double calcCalphaP(unsigned sided, unsigned nPermutations, unsigned adaptive);


	//!\brief RareCover method, 2010 PLoS CompBio
	double calcRareCoverP(unsigned sided, unsigned nPermutations, unsigned adaptive);

	//!\brief Weighted fisher's test, 2011 Shuang Wang
	double calcWsFisherP(unsigned sided, unsigned nPermutations, unsigned adaptive, bool isMidPvalue, const char moi);

	//!\brief skat test, 2011 xlin
	double calcSkatP(unsigned nPermutations, unsigned adaptive);

private:
	//!\brief 2D object, Sample genotypes
	vector2F __xdat;
	//!\brief 1D object, phenotypes
	vectorF __ydat;
	//!\brief pre-set (desired) significance level
	double __alpha;
	//!\brief Observed minor allele frequency in the data-set
	vectorF __observedMafs;
	//!\brief Upper and lower bounds of MAF that shall be analyzed. Many rare variants methods have it as (0, 0.01]
	double __mafLower;
	//!\brief Upper and lower bounds of MAF that shall be analyzed. Many rare variants methods have it as (0, 0.01]
	double __mafUpper;
	//!\brief Logical variable, flag for debug mode
	bool __isDebug;


	/*!\brief
	 * \param
	 * \return
	 */
	void m_trimXdat();

	void m_maskWildtypeSibpair();

	/*!\brief Simple single variable logistic regression score statistic
	 * \param regressors
	 * \param responses
	 * \param xbar
	 * \param nCases
	 * \return statistic
	 */
	double m_calcSimpleLogitRegScore(const vectorF & regressors, const vectorF & responses, double xbar, unsigned nCases, unsigned sided) const;


	/*!\brief 2 by 2 chisq test statistic
	 * \param regressors
	 * \param responses
	 * \return statistic
	 */
	double m_calc2X2Chisq(const vectorF & regressors, const vectorF & responses) const;


	/*!\brief 2 by 2 Fisher's test statistic
	 * \param regressors
	 * \param responses
	 * \return p-value
	 */
	double m_calc2X2Fisher(const vectorF & regressors, const vectorF & responses, const char moi, unsigned sided) const;

	double m_calc2X2Fisher(const vectorF & regressors, const vectorF & responses, unsigned sided) const;


	/*!\brief Simple linear regression score statisitc
	 * \param regressors
	 * \param responses
	 * \param xbar
	 * \param ybar
	 * \param sided 0 for the score, 1 and 2 for the p-value
	 * \return statistic or p-value
	 */
	double m_calcSimpleLinearRegScore(const vectorF & regressors, const vectorF & responses, double xbar, double ybar, unsigned sided) const;

	double m_calcConditionalLinearRegScore(const vectorF & regressors, const vectorF & responses, double yh, double yl, unsigned sided) const;

	/*!\brief two-sample t statistic
	 * \param x1s
	 * \param x2s
	 * \return statistic
	 */
	double m_calc2sampleT(const vectorF & x1s, const vectorF & x2s, unsigned sided) const;


	/*!\brief ANRV scoring theme
	 * \return
	 */
	vectorF m_countRegionalVariants() const;


	/*!\brief CMC scoring theme
	 * \return
	 */
	vectorF m_indicateRegionalVariants() const;


	/*!\brief RVE scoring theme
	 * \return
	 */
	vectorF m_indicateRegionalUniqueVariants() const;


	/*\brief Adaptive p-value calculation
	   *\return a logical variable indicating whether or not should continue permutations
	 */
	double m_checkAdaptivePvalue(unsigned permcount1, unsigned permcount2, unsigned currentIdx, unsigned checkPoint, unsigned sided) const;


	void m_printPunches(int n) const;

};
#endif ///:~

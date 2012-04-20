//!\file assoctests.h
//!\brief Statistical tests for association mapping
// Copyright 2010 2011 Gao Wang

/*  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *
 *                                                                          *
 *  This program is free software: you can redistribute it and/or modify    *
 *  it under the terms of the GNU General Public License as published by    *
 *  the Free Software Foundation, either version 3 of the License, or       *
 *  (at your option) any later version.                                     *
 *                                                                          *
 *  This program is distributed in the hope that it will be useful,         *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *  GNU General Public License for more details.                            *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License       *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *                                                                          *
 *  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   */  

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

class gwAssociations {

  public:


    /*!\brief Constructor, verbose
     * \param observedMafs Observed minor allele frequencies in the data-set
     * \param ydat
     * \param xdat
     * \param mafLower Lower bound of MAF that shall be analyzed. Many rare variants methods have it as (0, 0.01]
     * \param mafUpper Upper bound
     * \param alpha pre-set (desired) significance level
     */
    gwAssociations(const vectorF& observedMafs, const vectorF& ydat, 
        const vector2F& xdat, double mafLower, double mafUpper, double alpha);


    /*!\brief Deconstrutor, nothing to do
    */
    ~gwAssociations();


    /*!\brief Verbose screen out-put, for debug mainly
     */
    void setVerbose(int verbose);
    void debug(int showWhat) const;


    //!\brief Collapsed variants simple LR model score test statistic
    double calcCmcstP(UINT sided, UINT nPermutations, UINT adaptive);
  

    //!\brief Counts of variants simple LR model score test statistic
    double calcAnrvstP(UINT sided, UINT nPermutations, UINT adaptive);


    //!\brief Collapsed variants 2 by 2 chisq statistic
    double calcCmcchiP(UINT sided, UINT nPermutations, UINT adaptive);


    //!\brief CMC Fisher's exact test p-value no permutation
    double calcCmcfisherP();


    //!\brief Cohen's RVE fisher exact test p-value no permutation
    double calcRvefisherP();
  

    //!\brief The WSS method rank statistic, page 3 of Browning paper.  
    //!- Not using normal approximation
    double calcWssRankP(char* moi, UINT nCtrls, UINT sided, UINT nPermutations, UINT adaptive);
    //!- Using normal approximation
    double calcWssRankP(char* moi, UINT nCtrls);


    //!\brief The KBAC method, Dajiang Liu's paper
    double calcKbacP(UINT sided, UINT nPermutations, UINT adaptive, bool squared);
    

    //!\brief The variable threshold method, Price's paper 2010 AJHG
    double calcVtP(UINT sided, UINT nPermutations, UINT adaptive);


    //!\brief The adaptive sum test (in Pan and Han 2010) + permutation
    double calcAsumP(UINT sided, UINT nPermutations, UINT adaptive);

    
    //!\brief CMCQT 
    double calcCmcqtP(UINT sided, UINT nPermutations, UINT adaptive);
    
    
    //!\brief ANRVQT
    double calcAnrvqtP(UINT sided, UINT nPermutations, UINT adaptive);

      
    //!\brief Ionita-Laza and Lange 2011 PLoS Genet.
    double calcTestRareP(UINT sided, UINT nPermutations, UINT adaptive);

      
    //!\brief c-alpha test, Ben. Neale et al. 2011 PLoS Genet.
    double calcCalphaP(UINT sided, UINT nPermutations, UINT adaptive);


    //!\brief RareCover method, 2010 PLoS CompBio
    double calcRareCoverP(UINT sided, UINT nPermutations, UINT adaptive); 


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

    
    /*!\brief Simple single variable logistic regression score statistic
     * \param regressors
     * \param responses
     * \param xbar
     * \param nCases 
     * \return statistic 
     */ 
    double m_calcSimpleLogitRegScore(const vectorF& regressors, const vectorF& responses, double xbar, UINT nCases) const;


    /*!\brief 2 by 2 chisq test statistic
     * \param regressors
     * \param responses
     * \return statistic
     */ 
    double m_calc2X2Chisq(const vectorF& regressors, const vectorF& responses) const;


    /*!\brief 2 by 2 Fisher's test statistic
     * \param regressors
     * \param responses
     * \return p-value
     */ 
    double m_calc2X2Fisher(const vectorF& regressors, const vectorF& responses) const;


    /*!\brief Simple linear regression score statisitc
     * \param regressors
     * \param responses
     * \param xbar
     * \return statistic
     */ 
    double m_calcSimpleLinearRegScore(const vectorF& regressors, const vectorF& responses, double xbar) const;


    /*!\brief Simple linear regression GOF p-value
     * \param regressors
     * \param responses
     * \param xbar
     * \param ybar
     * \return statistic
     */ 
    double m_calcSimpleLinearRegPvalue(const vectorF& regressors, const vectorF& responses, double xbar, double ybar) const; 


    /*!\brief two-sample t statistic
     * \param x1s
     * \param x2s
     * \return statistic
     */ 
    double m_calc2sampleT(const vectorF& x1s, const vectorF& x2s) const;


    /*!\brief ANRV scoring theme
     * \return 
     */
    vectorF m_countRegionalVariants () const;


    /*!\brief CMC scoring theme
     * \return 
     */
    vectorF m_indicateRegionalVariants () const;
  

    /*!\brief RVE scoring theme
     * \return 
     */
    vectorF m_indicateRegionalUniqueVariants () const;


    /*\brief Adaptive p-value calculation
     *\return a logical variable indicating whether or not should continue permutations
     */
    double m_checkAdaptivePvalue(UINT permcount1, UINT permcount2, UINT currentIdx, UINT checkPoint, UINT sided) const;


    void m_printPunches(int n) const;
};
#endif ///:~

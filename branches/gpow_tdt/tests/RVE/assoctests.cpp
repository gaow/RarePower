//!\file assoctests.cpp
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

#include <iostream>
#include <cmath>
#include <algorithm>
#include <list>

#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"

#include "gw_maths.h"
#include "gw_utilities.h"
#include "assoctests.h"

namespace {
  const double AFFECTED = 2.0, UNAFFECTED = 1.0, HOMO_ALLELE = 2.0, 
        MINOR_ALLELE = 1.0, MAJOR_ALLELE = 0.0, MISSING_ALLELE = -9.0;
}

gwAssociations::gwAssociations( const vectorF& observedMafs, const vectorF& ydat, 
    const vector2F& xdat, double mafLower, double mafUpper, double alpha)
{
  bool isInputOk = (xdat.size() != 0 && ydat.size() == xdat.size() 
      && observedMafs.size() == xdat[0].size() && mafLower >= 0.0 
      && mafUpper <= 1.0 && mafUpper > mafLower && alpha > 0.0 && alpha < 1.0);

  if (isInputOk);
  else {
    std::cerr << "Input data problem in gwAssociations::gwAssociations(). Now Quit." << std::endl;
    exit(-1);
  }

  __xdat = xdat;
  __ydat = ydat;
  __observedMafs = observedMafs;
  __mafLower = mafLower;
  __mafUpper = mafUpper;
  __alpha = alpha;

  __isDebug = false;
  // Debug::
  //std::cout << __observedMafs << std::endl;
}


gwAssociations::~gwAssociations() {}


void gwAssociations::setVerbose(int verbose) 
{
  if (verbose != 0)
    __isDebug = true;
  else __isDebug = false;
  return;
}

void gwAssociations::debug(int showWhat) const 
{
  m_printPunches(60);
  std::cout.precision(9);

  switch (showWhat) {
    case 1 : 
      std::cout << "__observedMafs:\n" << __observedMafs << "\n" << std::endl;
      break;
    case 2 : 
      std::cout << "__ydat:\n" << __ydat << "\n" << std::endl;
      break;
    case 3 : 
      {
        std::cout << "__xdat[0]:\n" << __xdat[0] << "\n" << std::endl; 
        std::cout << "__xdat[__xdat.size()]:\n" << __xdat[__xdat.size()] << "\n" << std::endl;
      }
      break;
     default : 
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

double gwAssociations::calcCmcstP(UINT sided, UINT nPermutations, UINT adaptive) 
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
  UINT nCases = 0;
  for (UINT i = 0; i != __ydat.size(); ++i)
    if (__ydat[i] == AFFECTED)
      ++nCases;

  UINT iPermutation = 0;
  UINT permcount1 = 0, permcount2 = 0;
  double observedStatistic = 0.0;
  double pvalue = 9.0;

  while (iPermutation <= nPermutations) {
    double statistic = m_calcSimpleLogitRegScore(regressors, __ydat, xbar, nCases);
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

  if (pvalue <= 1.0);
  else {
    if (sided == 1) 
      pvalue = (1.0 * permcount1) / (1.0 * nPermutations);
    else {
      double permcount = gw_dmin(permcount1, permcount2);
      pvalue = (2.0 * permcount) / (1.0 * nPermutations);
    }
  }
  return pvalue;
}


//!\brief Counts of variants simple LR score test statistic

double gwAssociations::calcAnrvstP(UINT sided, UINT nPermutations, UINT adaptive) 
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
  UINT nCases = 0;
  for (UINT i = 0; i != __ydat.size(); ++i)
    if (__ydat[i] == AFFECTED)
      ++nCases;


  UINT iPermutation = 0;
  UINT permcount1 = 0, permcount2 = 0;
  double observedStatistic = 0.0;
  double pvalue = 9.0;

  while (iPermutation <= nPermutations) {
    double statistic = m_calcSimpleLogitRegScore(regressors, __ydat, xbar, nCases);
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

  if (pvalue <= 1.0);
  else {
    if (sided == 1) 
      pvalue = (1.0 * permcount1) / (1.0 * nPermutations);
    else {
      double permcount = gw_dmin(permcount1, permcount2);
      pvalue = (2.0 * permcount) / (1.0 * nPermutations);
    }
  }
  return pvalue;
}


double gwAssociations::calcCmcchiP(UINT sided, UINT nPermutations, UINT adaptive)
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

  UINT iPermutation = 0;
  UINT permcount1 = 0, permcount2 = 0;
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

  if (pvalue <= 1.0);
  else {
    if (sided == 1) 
      pvalue = (1.0 * permcount1) / (1.0 * nPermutations);
    else {
      double permcount = gw_dmin(permcount1, permcount2);
      pvalue = (2.0 * permcount) / (1.0 * nPermutations);
    }
  }
  return pvalue;
}


double gwAssociations::calcCmcfisherP() 
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

  //!- 2 by 2 Fisher's exact test, using the <b> fexact_two_sided_pvalue() </b> function
  double pvalue = m_calc2X2Fisher(regressors, __ydat);
  return pvalue;
}


double gwAssociations::calcRvefisherP() 
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

  //!- 2 by 2 Fisher's exact test, using the <b> fexact_two_sided_pvalue() </b> function
  double pvalue = m_calc2X2Fisher(regressors, __ydat);
  return pvalue;
}


double gwAssociations::calcWssRankP(char* moi, UINT nCtrls, UINT sided, UINT nPermutations, UINT adaptive) 
{
  /*! * Implement Browning 2009 Wss paper <br><br>
   * Implementation:
   */
  
  //!- trim data by mafs upper-lower bounds
  m_trimXdat();

  UINT nCases = __ydat.size() - nCtrls;

  
  UINT iPermutation = 0;
  UINT permcount1 = 0, permcount2 = 0;
  double observedStatistic = 0.0;
  double pvalue = 9.0;

  while (iPermutation <= nPermutations) {

    //!- Calc MAF in controls
    vectorF weights(0);
    for (UINT j = 0; j != __xdat[0].size(); ++j) { 

      UINT nVariants = 0;
      //! - Count number of mutations in controls
      for (UINT i = 0; i != __xdat.size(); ++i) {  
        // Control only
        if (__ydat[i] == UNAFFECTED) {  
          if (__xdat[i][j] != MISSING_ALLELE) 
            nVariants += (UINT)__xdat[i][j];  
          else {
            std::cerr << "Input data problem in gwAssociations::calcWssRankP(). Now Quit." << std::endl;
            exit(-1);
          }
        }
        else
          continue;
      }
      //! - Compute the "q" for the locus
      weights.push_back((nVariants + 1.0) / (2.0 * nCtrls + 2.0));
    }

    //!- Calc the Weights for each locus
    for (UINT j = 0; j != __xdat[0].size(); ++j) 
      weights[j] = sqrt(__xdat.size() * weights[j] * (1 - weights[j]));

    vectorF scores(0); 

    for (UINT i = 0; i != __xdat.size(); ++i) {
      double score = 0.0;  
      //! - Define genetic score of the individual

      for (UINT j = 0; j != __xdat[i].size(); ++j) {  
        // scan all loci of an individual
        double tmp = 0.0;  
        // Dominant model
        if (*moi == 'D') {
          if (__xdat[i][j] == HOMO_ALLELE) 
            tmp = 1.0; 
        }
        // recessive model
        else if (*moi == 'R') { 
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


    //! - Compute rank test statistic, via <b> Mann_Whitneyu() </b>
    double caseScores[nCases]; double ctrlScores[nCtrls]; int tmpa = 0; int tmpu = 0;
    for(UINT i = 0; i != __ydat.size(); ++i) { 
      if (__ydat[i] == AFFECTED) {
        caseScores[tmpa] = scores[i];
        ++tmpa;
      }
      else {
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

  if (pvalue <= 1.0);
  else {
    if (sided == 1) 
      pvalue = (1.0 * permcount1) / (1.0 * nPermutations);
    else {
      double permcount = gw_dmin(permcount1, permcount2);
      pvalue = (2.0 * permcount) / (1.0 * nPermutations);
    }
  }
  return pvalue;
}


double gwAssociations::calcWssRankP(char* moi, UINT nCtrls) 
{
  /*! * Implement Browning 2009 Wss paper normal approximation <br><br>
   * Implementation:
   */
  
  //!- trim data by mafs upper-lower bounds
  m_trimXdat();

  int nCases = __ydat.size() - nCtrls;

  UINT iPermutation = 0;
  double observedStatistic = 0.0;
  vectorF statistics(0);

  while (iPermutation <= 1000) {

    //!- Calc MAF in controls
    vectorF weights(0);
    for (UINT j = 0; j != __xdat[0].size(); ++j) { 

      int nVariants = 0;
      //! - Count number of mutations in controls
      for (UINT i = 0; i != __xdat.size(); ++i) {  
        // Control only
        if (__ydat[i] == UNAFFECTED) {  
          if (__xdat[i][j] != MISSING_ALLELE) 
            nVariants += (int)__xdat[i][j];  
          else {
            std::cerr << "Input data problem in gwAssociations::calcWssRankP. Now Quit." << std::endl;
            exit(-1);
          }
        }
        else
          continue;
      }
      //! - Compute the "q" for the locus
      weights.push_back((nVariants + 1.0) / (2.0 * nCtrls + 2.0));
    }

    //!- Calc the Weights for each locus
    for (UINT j = 0; j != __xdat[0].size(); ++j) 
      weights[j] = sqrt(__xdat.size() * weights[j] * (1 - weights[j]));

    vectorF scores(0); 

    for (UINT i = 0; i != __xdat.size(); ++i) {
      double score = 0.0;  
      //! - Define genetic score of the individual

      for (UINT j = 0; j != __xdat[i].size(); ++j) {  
        // scan all loci of an individual
        double tmp = 0.0;  
        // Dominant model
        if (*moi == 'D') {
          if (__xdat[i][j] == HOMO_ALLELE) 
            tmp = 1.0; 
        }
        // recessive model
        else if (*moi == 'R') { 
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


    //! - Compute rank test statistic, via <b> Mann_Whitneyu() </b>
    double caseScores[nCases]; double ctrlScores[nCtrls]; int tmpa = 0; int tmpu = 0;
    for(UINT i = 0; i != __ydat.size(); ++i) { 
      if (__ydat[i] == AFFECTED) {
        caseScores[tmpa] = scores[i];
        ++tmpa;
      }
      else {
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

  double statisticstd = pow((observedStatistic - mean) / sd, 2.0);  
  return gsl_cdf_chisq_Q(statisticstd, 1.0);
}


double gwAssociations::calcKbacP(UINT sided, UINT nPermutations, UINT adaptive, bool squared)
{
  /*! * the KBAC Statistic: sum of genotype pattern frequencies differences, weighted by hypergeometric kernel. <br>
   * * It is a permutation based two-sided test.  <br>
   * * See <em> Liu DJ 2010 PLoS Genet. </em> <br><br>
   *Implementation:
   */

  //!- trim data by mafs upper-lower bounds
  m_trimXdat();

  UINT sampleSize = __ydat.size();
  // sample size
  UINT regionLen = __xdat[0].size();
  // candidate region length
  UINT nCases = 0; 
  // case size

  for(UINT i = 0; i !=  __ydat.size(); ++i) {
    if(__ydat[i] == AFFECTED) 
      ++nCases;
  }

  vectorF genotypeId(sampleSize);   
  // Compute unique genotype patterns (string) as ID scores (double)
  //!-Compute unique genotype patterns (string) as ID scores (double) 

  for (UINT i = 0; i != sampleSize; ++i) {

    double vntIdL = 0.0; 
    double vntIdR = 0.0;
    const double ixiix= pow(9.0, 10.0);
    UINT lastCnt = 0;
    UINT tmpCnt = 0;

    for (UINT j = 0; j != regionLen; ++j) { 

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
      }
      else { 
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
  UINT uniquePatternCounts[uniquePattern.size()];
  for (UINT u = 0; u != uniquePattern.size(); ++u) 
    uniquePatternCounts[u] = 0;

  for (UINT i = 0; i != sampleSize; ++i) {
    // for each sample, identify/count its genotype pattern

    for (UINT u = 0; u != uniquePattern.size(); ++u) {

      if (genotypeId[i] == uniquePattern[u]) {
        // genotype pattern identified
        ++uniquePatternCounts[u];
        // count this genotype pattern
        break;
      }
      else;
      // genotype pattern not found -- move on to next pattern
    }
  }


  UINT iPermutation = 0;
  UINT permcount1 = 0, permcount2 = 0;
  double observedStatistic = 0.0;
  double pvalue = 9.0;
  while (iPermutation <= nPermutations) {

    //!- count number of sample cases for each genotype pattern
    UINT uniquePatternCountsCs[uniquePattern.size()];
    for (UINT u = 0; u != uniquePattern.size(); ++u) 
      uniquePatternCountsCs[u] = 0;
    // genotype pattern counts in cases

    for (UINT i = 0; i != sampleSize; ++i) {
      if (__ydat[i] == AFFECTED) {
        // for each "case", identify/count its genotype pattern
        for (UINT u = 0; u != uniquePattern.size(); ++u) {
          if (genotypeId[i] == uniquePattern[u]) {
            // genotype pattern identified in cases
            ++uniquePatternCountsCs[u];
            // count this genotype pattern
            break;
          }
          else;
          // genotype pattern not found -- move on to next pattern
        }
      }
      else;
    }

    //!- KBAC weights
    double uniquePatternWeights[uniquePattern.size()];
    for (UINT u = 0; u != uniquePattern.size(); ++u) 
      uniquePatternWeights[u] = 0.0;
    // genotype pattern weights
    for (UINT u = 0; u != uniquePattern.size(); ++u)
      // genotype pattern weights, the hypergeometric distribution cmf
      uniquePatternWeights[u] = gsl_cdf_hypergeometric_P(uniquePatternCountsCs[u], uniquePatternCounts[u], sampleSize - uniquePatternCounts[u], nCases);
    //uniquePatternWeights[u] = gw_hypergeometric_cmf(uniquePatternCountsCs[u], uniquePatternCounts[u], sampleSize - uniquePatternCounts[u], nCases);

    //!- KBAC statistic
    double statistic = 0.0;
    for (UINT u = 0; u != uniquePattern.size(); ++u) 
      statistic = statistic + ( (1.0 * uniquePatternCountsCs[u]) / (1.0 * nCases) - (1.0 * (uniquePatternCounts[u] - uniquePatternCountsCs[u])) / (1.0 * (sampleSize - nCases)) ) *  uniquePatternWeights[u];


    if (squared) 
      statistic = statistic * statistic;
    // statistic: sum of genotype pattern frequencies differences in cases vs. controls, weighted by the hypergeometric distribution kernel
    // test is two-sided
    
    if (__isDebug)
      std::cout << statistic << std::endl;

    //FIXME
    //gw_round(statistic, 0.0001);
    
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

  if (pvalue <= 1.0);
  else {
    if (sided == 1) 
      pvalue = (1.0 * permcount1) / (1.0 * nPermutations);
    else {
      double permcount = gw_dmin(permcount1, permcount2);
      pvalue = (2.0 * permcount) / (1.0 * nPermutations);
    }
  }
  return pvalue;
}


double gwAssociations::calcVtP(UINT sided, UINT nPermutations, UINT adaptive) 
{
  /*! * Variable threshold method, Price 2010 AJHG <br>
   * * <em> Unique feature of VT </em>: instead of using fixed threshold + Morris A RV counts, it uses a variable threshold for RV frequency. <br>
   * * Since it uses Z_max = max(Z_i ...) where Z_i is proportional to standard normal by a coefficient which is not a funciton of RV counts, it does not matter to include the "known" common variants. <br> 
   * * In brief it is a multiple testing method that analyzes both common and rare variants <br><br>
   * Implementation:
   */

  UINT nCases = 0; 
  // case size

  for(UINT i = 0; i !=  __ydat.size(); ++i) {
    if(__ydat[i] == AFFECTED) 
        ++nCases;
  }

  double po = (1.0 * nCases) / (1.0 * __ydat.size());


  vectorF vaCnt(0);
  // number of variants at each loci
  for (UINT j = 0; j != __xdat[0].size(); ++j) { 

    int tmpVaN = 0;
    for (UINT i = 0; i != __xdat.size(); ++i)   
      tmpVaN += (int)__xdat[i][j];  

    vaCnt.push_back(tmpVaN);
  }


  std::list<double> lvts(vaCnt.begin(), vaCnt.end());
  lvts.remove(MAJOR_ALLELE); 
  if (lvts.size() == 0) 
    return 1.0;
  lvts.sort(); 
  lvts.unique();
  //!- Compute the unique non-zero variant counts at each loci, the "variable threshold"

  vectorF vts(lvts.size());
  copy(lvts.begin(), lvts.end(), vts.begin());
  lvts.clear();


  UINT iPermutation = 0;
  UINT permcount1 = 0, permcount2 = 0;
  double observedStatistic = 0.0;
  double pvalue = 9.0;

  while (iPermutation <= nPermutations) {

    vectorF allZs(0);
    //! - Define <b> 'allZs' </b>, a vector of the Z scores computed under different thresholds

    for (UINT t = 0; t != vts.size(); ++t) {
      //! - Iterate the following for all thresholds:

      vectorI underVt(0);
      //! - - Record the index of loci that have smaller variants count than the threshold
      for (UINT j = 0; j != vaCnt.size(); ++j) 
        if (vaCnt[j] <= vts[t]) 
          underVt.push_back(j);

      //!- - For loci passing the threshold, implement Price paper page 3 z(T) formula
      vectorF zIa(0), zIb(0);
      for (UINT i = 0; i != __xdat.size(); ++i) {  

        double zIIa = 0, zIIb = 0;
        for (UINT j = 0; j != underVt.size(); ++j) { 
          //  All loci that are under the threshold

          int locIdx = underVt[j];
          if ( __xdat[i][locIdx] != MAJOR_ALLELE ) {

            double cij = __xdat[i][locIdx];
            zIIa = zIIa + ((__ydat[i] - 1.0) - po) * cij;
            zIIb = zIIb + cij * cij;
          }
        }

        zIa.push_back(zIIa); zIb.push_back(zIIb);
      }

      allZs.push_back( gw_sum(zIa) / sqrt( gw_sum(zIb) ) );
      //! - Now each element in <b> 'allZs' </b> is z(T) for different thresholds
    }

    //! - Compute zmax, the statistic 
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

  if (pvalue <= 1.0);
  else {
    if (sided == 1) 
      pvalue = (1.0 * permcount1) / (1.0 * nPermutations);
    else {
      double permcount = gw_dmin(permcount1, permcount2);
      pvalue = (2.0 * permcount) / (1.0 * nPermutations);
    }
  }
  return pvalue;
}


double gwAssociations::calcAsumP(UINT sided, UINT nPermutations, UINT adaptive)
{
  /*! * Number of rare variants per site with "protective" site recoded, by Pan and Han (2010) Hum Hered <br>
   *  * The authors use alpha0 = 0.1 to screen variants that should be recoded. See their paper for details <br>
   *  * Original paper uses LR score test with approx. distribution, which is equivalent to Armitage's test for trend for this case (Schaid 2002). For RV not many sites are homozygote. In my implementation I use Fisher's 2by2 test instead.
   * * Here allows user specified lower and upper bounds <br><br>
   * Implementation:
   */

  //!- trim data by mafs upper-lower bounds
  m_trimXdat();
  UINT nCases = 0;
  for (UINT i = 0; i != __ydat.size(); ++i)
    if (__ydat[i] == AFFECTED)
      ++nCases;

  UINT iPermutation = 0;
  UINT permcount1 = 0, permcount2 = 0;
  double observedStatistic = 0.0;
  double pvalue = 9.0;

  while (iPermutation <= nPermutations) {

    //!- Sites that have excess rare variants in controls
    vectorL ctrlExcess(__xdat[0].size(), false);

    for (UINT j = 0; j != __xdat[0].size(); ++j) {
      double tmpCase = 0.0;
      double tmpCtrl = 0.0;
      for (UINT i = 0; i != __xdat.size(); ++i) {
        if (__ydat[i] == UNAFFECTED) {
          if (__xdat[i][j] != MISSING_ALLELE && __xdat[i][j] != MAJOR_ALLELE) 
            tmpCtrl += __xdat[i][j]; 
        }
        else {
          if (__xdat[i][j] != MISSING_ALLELE && __xdat[i][j] != MAJOR_ALLELE) 
            tmpCase += __xdat[i][j]; 
        }
      }
      if (tmpCtrl > tmpCase) 
        ctrlExcess[j] = true; 
      else 
        continue;
    }

    //!- Define sites that needs to be recoded
    vectorL recodeSites(__xdat[0].size(), false);

    for(UINT j = 0; j != __xdat[0].size(); ++j) { 
      if (ctrlExcess[j] == false) 
        continue;

      vectorF vdat(0);
      for (UINT i = 0; i != __xdat.size(); ++i) 
        vdat.push_back(__xdat[i][j]);  

      //!- 2 by 2 Fisher's test
      double pfisher = m_calc2X2Fisher(vdat, __ydat);       

      if (pfisher < 0.1) 
        recodeSites[j] = true; 
      else 
        continue;
    }

    vectorF regressors(0);

    //! - Count anrv; recode when necessary
    for (UINT i = 0; i != __xdat.size(); ++i) {  
      //  scan all sample individuals

      double nrv = 0.0;
      for(UINT j = 0; j != __xdat[0].size(); ++j) { 
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
    double statistic = m_calcSimpleLogitRegScore(regressors, __ydat, xbar, nCases);

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

  if (pvalue <= 1.0);
  else {
    if (sided == 1) 
      pvalue = (1.0 * permcount1) / (1.0 * nPermutations);
    else {
      double permcount = gw_dmin(permcount1, permcount2);
      pvalue = (2.0 * permcount) / (1.0 * nPermutations);
    }
  }
  return pvalue;
}


double gwAssociations::calcCmcqtP(UINT sided, UINT nPermutations, UINT adaptive) 
{
  //!- trim data by mafs upper-lower bounds
  m_trimXdat();

  //!- CMC scoring
  vectorF regressors = m_indicateRegionalVariants(); 

  if (__isDebug) 
    std::cout << regressors << std::endl;


  UINT iPermutation = 0;
  UINT permcount1 = 0, permcount2 = 0;
  double observedStatistic = 0.0;
  double pvalue = 9.0;

  while (iPermutation <= nPermutations) {

    //!- two sample t test. Group 1: have cmc score = 1, Group 2: have cmc score = 0
    double statistic = m_calc2sampleT(regressors, __ydat);
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

  if (pvalue <= 1.0);
  else {
    if (sided == 1) 
      pvalue = (1.0 * permcount1) / (1.0 * nPermutations);
    else {
      double permcount = gw_dmin(permcount1, permcount2);
      pvalue = (2.0 * permcount) / (1.0 * nPermutations);
    }
  }
  return pvalue;
}


double gwAssociations::calcAnrvqtP(UINT sided, UINT nPermutations, UINT adaptive) 
{

  //!- trim data by mafs upper-lower bounds
  m_trimXdat();
  //!- ANRV scoring
  vectorF regressors = m_countRegionalVariants(); 
  if (__isDebug) 
    std::cout << regressors << std::endl;

  //!- Linear regression test 
  double xbar = gw_mean(regressors);

  UINT iPermutation = 0;
  UINT permcount1 = 0, permcount2 = 0;
  double observedStatistic = 0.0;
  double pvalue = 9.0;

  while (iPermutation <= nPermutations) {
    double statistic = m_calcSimpleLinearRegScore(regressors, __ydat, xbar);
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

  if (pvalue <= 1.0);
  else {
    if (sided == 1) 
      pvalue = (1.0 * permcount1) / (1.0 * nPermutations);
    else {
      double permcount = gw_dmin(permcount1, permcount2);
      pvalue = (2.0 * permcount) / (1.0 * nPermutations);
    }
  }
  return pvalue;
}


double gwAssociations::calcTestRareP(UINT sided, UINT nPermutations, UINT adaptive)
{
  /*! * the RVP Statistic: variant counts weighted by poisson kernels. <br>
   * * It is a permutation based two model test, max(protective model, risk model).  <br>
   * * See <em> Ionita-Laza and Lange 2011 PLoS Genet. </em> <br><br>
   *Implementation (credit to the authors for part of the codes below):
   */

  //!- trim data by mafs upper-lower bounds
  m_trimXdat();

  UINT sampleSize = __xdat.size();
  // sample size
  UINT regionLen = __xdat[0].size();
  // candidate region length

  UINT nCases = 0; 
  for(UINT i = 0; i !=  __ydat.size(); ++i) {
    if(__ydat[i] == AFFECTED) 
      ++nCases;
  }
  UINT nCtrls = sampleSize - nCases;


  UINT iPermutation = 0;
  UINT permcount1 = 0, permcount2 = 0;
  double observedStatistic = 0.0;
  double pvalue = 9.0;

  while (iPermutation <= nPermutations) {

    double sumR = 0, sumP = 0;
    for(UINT j = 0; j != regionLen; ++j) {  

      //! - Count number of variants in cases/controls at a locus
      UINT countcs = 0;
      UINT countcn = 0;
      for(UINT i = 0; i != sampleSize; ++i) { 
        if (__ydat[i] == UNAFFECTED) 
          countcn += (UINT) __xdat[i][j]; 
        else 
          countcs += (UINT)__xdat[i][j]; 
      }

      //! - the RVP method. Codes adopted from the author's implementation
      double f=1.0;
      //float w;  
      //k0=(int)(freq*2*nCtrls);
      if (countcs>0) 
        f=gsl_cdf_poisson_P(countcn,nCtrls*(countcs+countcn)/(1.0*(nCases+nCtrls)))*(1.0-gsl_cdf_poisson_P(countcs-1,nCases*(countcs+countcn)/(1.0*(nCases+nCtrls))));

      if ((1.0*countcs)/(1.0*nCases) > (1.0*countcn)/(1.0*nCtrls)) 
        sumR -= log(f); 

      f=1.0;
      //k0=(int)(freq*2*nCases);
      if (countcn>0) 
        f=gsl_cdf_poisson_P(countcs,nCases*(countcs+countcn)/(1.0*(nCases+nCtrls)))*(1.0-gsl_cdf_poisson_P(countcn-1,nCtrls*(countcs+countcn)/(1.0*(nCases+nCtrls))));
      if ((1.0*countcn)/(1.0*nCtrls) > (1.0*countcs)/(1.0*nCases)) 
        sumP -= log(f);
    }

    //!- RVP statistic: The max statistic in the manuscript: R - potentially risk and P - potentially protective
    double statistic = fmax(sumR,sumP);

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

  if (pvalue <= 1.0);
  else {
    if (sided == 1) 
      pvalue = (1.0 * permcount1) / (1.0 * nPermutations);
    else {
      double permcount = gw_dmin(permcount1, permcount2);
      pvalue = (2.0 * permcount) / (1.0 * nPermutations);
    }
  }
  return pvalue;
}


double gwAssociations::calcCalphaP(UINT sided, UINT nPermutations, UINT adaptive)
{

  /*! * the c-alpha Statistic: sum of the std. error of variant counts in cases <br>
   * *  One-sided test, claims to be robust against protective mutations. <br>
   * * See <em> Ben. Neale et al. 2011 PLoS Genet. </em> <br><br>
   *Implementation:
   */

  //!- trim data by mafs upper-lower bounds
  m_trimXdat();

  UINT sampleSize = __xdat.size();
  // sample size
  UINT regionLen = __xdat[0].size();
  // candidate region length

  UINT nCases = 0; 
  for(UINT i = 0; i !=  __ydat.size(); ++i) {
    if(__ydat[i] == AFFECTED) 
      ++nCases;
  }

  double phat = (1.0 * nCases)/(1.0 * sampleSize);


  UINT iPermutation = 0;
  UINT permcount1 = 0, permcount2 = 0;
  double observedStatistic = 0.0;
  double pvalue = 9.0;

  while (iPermutation <= nPermutations) {

    double calpT = 0.0;
    double calpV = 0.0;
    UINT singletonAll = 0;
    UINT singletonCases = 0;
    bool isEmptyData = true;

    for(UINT j = 0; j != regionLen; ++j) {  

      //! - Count number of variants in cases/controls at a locus
      UINT countcs = 0;
      UINT countcn = 0;
      for(UINT i = 0; i != sampleSize; ++i) { 
        if (__ydat[i] == UNAFFECTED) 
          countcn += (UINT) __xdat[i][j]; 
        else 
          countcs += (UINT) __xdat[i][j]; 
      }

      //! - the c-alpha method implementation
      UINT ni = countcs + countcn;
      if (ni < 2) {
        singletonAll += ni;
        singletonCases += countcs;
        continue;
      } 
      else 
        isEmptyData = false;
      //!- * skip singletons

      double niv =  ni * phat * (1 - phat);
      calpT += (countcs - ni * phat) * (countcs - ni * phat) - niv;
      for (UINT u = 0; u <= ni; ++u) {
        double tmess = (u - ni * phat) * (u - ni * phat) - niv;
        calpV += tmess * tmess * gsl_ran_binomial_pdf(u, phat, ni);
      }
    }

    if (singletonAll >= 2) {
      isEmptyData = false; 
      double niv =  singletonAll * phat * (1 - phat);
      calpT += (singletonCases - singletonAll * phat) * (singletonCases - singletonAll * phat) - niv;
      for (UINT u = 0; u <= singletonAll; ++u) {
        double tmess = (u - singletonAll * phat) * (u - singletonAll * phat) - niv;
        calpV += tmess * tmess * gsl_ran_binomial_pdf(u, phat, singletonAll);
      }
    } else; 
    //!- * bin singletons

    if (isEmptyData) 
      return 1.0; 

    //!- c-alpha statistic, two sided. 
    //!- Note: the original permutation based c-alpha test is one-sided but to play fair here we use the two sided. 
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

  if (pvalue <= 1.0);
  else {
    if (sided == 1) 
      pvalue = (1.0 * permcount1) / (1.0 * nPermutations);
    else {
      double permcount = gw_dmin(permcount1, permcount2);
      pvalue = (2.0 * permcount) / (1.0 * nPermutations);
    }
  }
  return pvalue;
}


double gwAssociations::calcRareCoverP(UINT sided, UINT nPermutations, UINT adaptive) 
{
  /*! * RareCover method, 2010 PLoS CompBio <br>
   * * Implementation:
   */

  //!- the cut-off to use for the "heuristic greedy algorithm". = 0.5 as suggested by the paper
  const double difQ = 0.5;

  //!- Index of loci that are observed to be polymophic
  vectorUI vntVct(0);

  for (UINT j = 0; j != __observedMafs.size(); ++j) {
    if (__observedMafs[j] > 0.0) 
      vntVct.push_back(j + 1);
  }    

  if (vntVct.size() == 0) 
    return 1.0;

  UINT smpSize = __xdat.size();


  UINT iPermutation = 0;
  UINT permcount1 = 0, permcount2 = 0;
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
      UINT rmIdx = 0;
      //!- the "test contributing" variant index, for the vntVct object
      bool rmIdxFlag = false;

      for (UINT t = 0; t != vntNow.size(); ++t) {

        if (vntNow[t] == 0) 
          continue;
        UINT iIdx = vntNow[t] - 1;
        //!- the index of a variant site

        for (UINT i = 0; i != smpSize; ++i) 
          regressors[i] = (regressorsCurr[i] + __xdat[i][iIdx] > 0) ? MINOR_ALLELE : MAJOR_ALLELE;

        //! - 2 by 2 Chisq test
        double statistic = m_calc2X2Chisq(regressors, __ydat);

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
        UINT rmVnt = vntNow[rmIdx] - 1;
        //!- Update the genotype coding by adding in the contributing locus
        for (UINT i = 0; i != smpSize; ++i) 
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

  if (pvalue <= 1.0);
  else {
    if (sided == 1) 
      pvalue = (1.0 * permcount1) / (1.0 * nPermutations);
    else {
      double permcount = gw_dmin(permcount1, permcount2);
      pvalue = (2.0 * permcount) / (1.0 * nPermutations);
    }
  }
  return pvalue;
}

/////////////////////////////////////////////////

void gwAssociations::m_trimXdat()
{
  vector2F xdat = __xdat;
  __xdat.clear();
  __xdat.resize(xdat.size());

  for (UINT j = 0; j != __observedMafs.size(); ++j) {
    if (__observedMafs[j] <= __mafLower || __observedMafs[j] > __mafUpper) 
      continue;

    else {
      for (UINT i = 0; i != xdat.size(); ++i)
        __xdat[i].push_back(xdat[i][j]);
    }
  }  
  return;
}  


double gwAssociations::m_calcSimpleLogitRegScore(const vectorF& regressors, const vectorF& responses, double xbar, UINT nCases) const 
{
  bool isInputOk = (regressors.size() == responses.size());
  if (isInputOk);
  else {
    std::cerr << "Input data problem in m_calcSimpleLogitRegScore(). Now Quit." << std::endl;
    exit(-1);
  }

  //!- Score test implementation for logistic regression model logit(p) = b0 + b1x (derivation of the test see my labnotes vol.2 page 3)

  //double ebo = (1.0 * nCases) / (1.0 * numCtrl);
  //double bo = log(ebo);
  
  double po = (1.0 * nCases) / (1.0 * responses.size());

  double ss = 0.0;
  // the score
  for (UINT i = 0; i != regressors.size(); ++i) 
    ss += (regressors[i] - xbar) * ((responses[i] - 1.0) - po);
  double vm1 = 0.0;
  // variance^-1 of score
  for (UINT i = 0; i != regressors.size(); ++i) 
    vm1 += (regressors[i] - xbar) * (regressors[i] - xbar) * po * (1.0 - po);

  double statistic = ss * sqrt(vm1);

  if (__isDebug) { 
    std::cout << std::endl;
    std::cout << xbar <<  std::endl;
    std::cout << ss <<  std::endl;
    std::cout << vm1 << std::endl;
    exit(1);
  }

  statistic = statistic * statistic;
  
  //!-FIXME: not sure why if I do not do this I get strange type I error problem because I get many number such as 3.72397e-35 etc
  gw_round(statistic, 0.0001);
  
  if (__isDebug)
    std::cout << statistic << std::endl;

  return (statistic);
}


double gwAssociations::m_calc2X2Chisq(const vectorF& regressors, const vectorF& responses) const
{
  bool isInputOk = (regressors.size() == responses.size());
  if (isInputOk);
  else {
    std::cerr << "Input data problem in m_calc2X2Chisq(). Now Quit." << std::endl;
    exit(-1);
  }

  //! - 2 by 2 Chisq test
  double A0 = 0.0, A1 = 0.0, U0 = 0.0, U1 = 0.0;
  for (UINT i = 0; i != regressors.size(); ++i) { 
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

  double statistic = ( ( A0 - EA0 ) * ( A0 - EA0 ) ) / EA0
    + ( ( A1 - EA1 ) * ( A1 - EA1 ) ) / EA1
    + ( ( U0 - EU0 ) * ( U0 - EU0 ) ) / EU0
    + ( ( U1 - EU1 ) * ( U1 - EU1 ) ) / EU1;

  return statistic;
}


double gwAssociations::m_calc2X2Fisher(const vectorF& regressors, const vectorF& responses) const
{
  bool isInputOk = (regressors.size() == responses.size());
  if (isInputOk);
  else {
    std::cerr << "Input data problem in m_calc2X2Fisher() Initialization. Now Quit." << std::endl;
    exit(-1);
  }

  vectorI twotwoTable(4, 0);
  for (UINT i = 0; i != regressors.size(); ++i) { 
    if (responses[i] == AFFECTED) {
      if (regressors[i] == MAJOR_ALLELE) 
        twotwoTable[0] += 1; 
      else if (regressors[i] == MINOR_ALLELE || regressors[i] == HOMO_ALLELE) 
        twotwoTable[1] += 1;
    }
    else if (responses[i] == UNAFFECTED) {
      if (regressors[i] == MAJOR_ALLELE) 
        twotwoTable[2] += 1; 
      else if (regressors[i] == MINOR_ALLELE || regressors[i] == HOMO_ALLELE) 
        twotwoTable[3] += 1; 
    }
    else {
      std::cerr << "Input data problem in m_calc2X2Fisher() table. Now Quit." << std::endl; 
      exit(-1);
    }
  }

  return fexact_two_sided_pvalue(twotwoTable);
}


double gwAssociations::m_calcSimpleLinearRegScore(const vectorF& regressors, const vectorF& responses, double xbar) const 
{
  bool isInputOk = (regressors.size() == responses.size());
  if (isInputOk);
  else {
    std::cerr << "Input data problem in m_calcSimpleLinearRegScore(). Now Quit." << std::endl;
    exit(-1);
  }


  double numerator = 0.0, denominator = 0.0;


  for (UINT i = 0; i != regressors.size(); ++i) {

    numerator += (regressors[i] - xbar) * responses[i];
    denominator += (regressors[i] - xbar) * (regressors[i] - xbar); 	
  }

  if (denominator == 0.0) 
    denominator = 1.0e-6;
  double lseb = numerator / sqrt(denominator);
  //!- Statistic: LSE (MLE) for beta, centered and scaled (bcz E[b] = 0 and sigma = 1 by simulation) See page 41 of Kutner's Applied Linear Stat. Model, 5th ed.
  return lseb * lseb;
}


double gwAssociations::m_calcSimpleLinearRegPvalue(const vectorF& regressors, const vectorF& responses, double xbar, double ybar) const 
{
  bool isInputOk = (regressors.size() == responses.size());
  if (isInputOk);
  else {
    std::cerr << "Input data problem in m_calcSimpleLinearRegPvalue(). Now Quit." << std::endl;
    exit(-1);
  }


  double numerator = 0.0, denominator = 0.0, ysigma = 0.0; 


  for (UINT i = 0; i != regressors.size(); ++i) {

    numerator += (regressors[i] - xbar) * responses[i];
    denominator += (regressors[i] - xbar) * (regressors[i] - xbar); 	
    ysigma += (responses[i] - ybar) * (responses[i] - ybar);
  }
  ysigma = ysigma / responses.size();

  if (denominator == 0.0) 
    denominator = 1.0e-6;        
  //!- Compute MSE
  double varb = ysigma / denominator;
  //!- E[\hat{beta}] = MSE / denominator
  double statistic = (numerator / denominator) * (numerator / denominator) / varb;
  //!- two-sided pvalue
  return gsl_cdf_chisq_Q(statistic, 1.0);
  //!- Statistic: LSE (MLE) for beta, centered and scaled (bcz E[b] = 0 and sigma = 1 by simulation) 
  //!- See page 41 of Kutner's Applied Linear Stat. Model, 5th ed.
}


double gwAssociations::m_calc2sampleT(const vectorF& x1s, const vectorF& x2s) const
{
  bool isInputOk = (x1s.size() == x2s.size());
  if (isInputOk);
  else {
    std::cerr << "Input data problem in m_calc2sampleT(). Now Quit." << std::endl;
    exit(-1);
  }

  double XA = 0.0, XU = 0.0, SA = 0.0, SU = 0.0;
  vectorF  VA(0), VU(0);
  for (UINT i = 0; i != x1s.size(); ++i) {
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

  double statistic = (XA - XU) * (XA - XU) / (SA / VA.size() + SU / VU.size());
  return statistic;

}


vectorF gwAssociations::m_countRegionalVariants () const
{
  /*! * Number of rare variants per site, by Morris A (2009) Genet Epi <br>
   * * The authors use the MAF cut-off of 5% but here allows user specified lower and upper bounds <br><br>
   */

  vectorF data(0);
  //!- Counts of variants in a region
  for (UINT i = 0; i != __xdat.size(); ++i) {  

    double nrv = 0.0; 
    for (UINT j = 0; j != __xdat[i].size(); ++j) {  
      if (__xdat[i][j] != MISSING_ALLELE)   
        nrv += __xdat[i][j];
    }
    data.push_back(nrv);  
  }
  return data;
}


vectorF gwAssociations::m_indicateRegionalVariants () const
{

  vectorF data(0);
  //!- Collapsing of variants in a region
  for (UINT i = 0; i != __xdat.size(); ++i) {  

    bool isWild = true; 
    for (UINT j = 0; j != __xdat[i].size(); ++j) {  
      // scan all genetic loci of an individual

      //! - Apply indicator function: rare variant found. score it '1' and break the loop 
      if (__xdat[i][j] != MAJOR_ALLELE && __xdat[i][j] != MISSING_ALLELE) {  
        data.push_back(1.0);
        isWild = false;
        break;
      }
      else; 
    }

    //! - Otherwise score it '0' after every locus of an indv has been scanned.
    if (isWild == true) 
      data.push_back(0.0);  
  }

  return data;
}


vectorF gwAssociations::m_indicateRegionalUniqueVariants () const
{

  vectorF data(0);

  //!- Identify variants observed only in cases or in controls. Exclude them otherwise

  vectorL areUnique(__xdat[0].size(), true);
  for (UINT j = 0; j != __xdat[0].size(); ++j) {

    bool caseFlag = false, ctrlFlag = false;
    for (UINT i = 0; i != __xdat.size(); ++i) {

      if (__ydat[i] == AFFECTED && caseFlag == false) {
        if (__xdat[i][j] != MAJOR_ALLELE && __xdat[i][j] != MISSING_ALLELE) 
          caseFlag = true; 
        else; 
      } 
      else if (__ydat[i] == UNAFFECTED && ctrlFlag == false) {
        if (__xdat[i][j] != MAJOR_ALLELE && __xdat[i][j] != MISSING_ALLELE) 
          ctrlFlag = true;
        else; 
      }
      else;

      if (caseFlag == true && ctrlFlag == true) 
        break; 
      else 
        continue;
    }

    areUnique[j] = ((caseFlag == true && ctrlFlag == false) || (caseFlag == false && ctrlFlag == true));
  }

  //!- Collapsing of variants for unique loci only  
  for (UINT i = 0; i != __xdat.size(); ++i) {  

    bool isWild = true; 
    for (UINT j = 0; j != __xdat[i].size(); ++j) {  
      // scan all genetic loci of an individual

      // skip non-unique site
      if (areUnique[j] == false) 
        continue; 
      else; 

      //! - Apply indicator function: rare variant found. score it '1' and break the loop 
      if (__xdat[i][j] != MAJOR_ALLELE && __xdat[i][j] != MISSING_ALLELE) {  
        data.push_back(1.0);
        isWild = false;
        break;
      }
      else; 
    }

    //! - Otherwise score it '0' after every locus of an indv has been scanned.
    if (isWild == true) 
      data.push_back(0.0);  
  }

  return data;
}


double gwAssociations::m_checkAdaptivePvalue(UINT permcount1, UINT permcount2, UINT currentIdx, UINT checkPoint, UINT sided) const
{
  if (currentIdx % checkPoint == 0) {
    //!- adaptive p-value calculation, at an interval of 5000 permutations 
    // apply the "six-sigma" rule

    double pval = 1.0;
    if (sided == 1) 
      pval = (1.0 * permcount1) / (1.0 * currentIdx);
    else {
      double permcount = gw_dmin(permcount1, permcount2);
      pval = (2.0 * permcount) / (2.0 * currentIdx);
    }

    double sd = sqrt(pval * (1.0 - pval) / (1.0 * currentIdx));
    double sixsigma = pval - 6.0 * sd;

    if (sixsigma > __alpha) 
      return pval;
    else
      return 9.0;
  }
  else 
    return 9.0;
}


void gwAssociations::m_printPunches(int n) const
{
  for (int i = 0; i != n; ++i)
    std::cout << "#";
  std::cout << "\n" << std::endl;
  return;
}

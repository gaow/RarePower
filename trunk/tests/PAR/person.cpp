//!\file person.cpp
//!\brief "person" object implementation
// Copyright 2010 Gao Wang

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
#include<cfloat>

#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"

#include "gw_maths.h"
#include "gw_utilities.h"
#include "person.h"

namespace {
  const int D_RV = 1, D_CV = 6, P_RV = -1, P_CV = -6, SYNO_RV = 15, SYNO_CV = 65, 
        N_RV = 17, N_CV = 67, MARK_MISSING = 9999, MARK_WILD = 1000;
  const double BOUNDARY = 0.01, AFFECTED = 2.0, UNAFFECTED = 1.0, UNPHENOTYPED = 0.0, 
        MINOR_ALLELE = 1.0, MAJOR_ALLELE = 0.0, MISSING_ALLELE = -9.0;
}

gwPerson::gwPerson() 
{
  __genos.resize(2);
  __phenos.resize(0);
  vectorF zeros(30, 1.0e-5);
  __mafs = zeros;
  __genoFreqs.resize(__mafs.size());
  __pedInfos.resize(6);
  vectorF missings(30, MISSING_ALLELE);
  __genos[0] = missings;
  __genos[1] = missings;
  __pedInfos[0] = 0.0;
  __pedInfos[1] = 0.0;
  __pedInfos[2] = 0.0;
  __pedInfos[3] = 0.0;
  __pedInfos[4] = 0.0;
  __pedInfos[5] = UNPHENOTYPED;
  
  for (UINT i = 0; i != __mafs.size(); ++i) {
    __genoFreqs[i].push_back( (1 - __mafs[i]) * (1 - __mafs[i]) );
    __genoFreqs[i].push_back( __mafs[i] * __mafs[i] );
  }
  
  __selcoefs = zeros;
  vectorI synonymous(30, SYNO_RV);
  __locusAttributes = synonymous;
  vectorUI positions(30, 0);
  __snvPositions = positions;
  __isDebug = false;
  
  // Debug::
  // std::cout << __pedInfos << std::endl;
}


gwPerson::gwPerson( const vectorF& pedInfos, const vectorF& mafs, 
    const vectorF& selcoefs, const vectorUI& snvPositions )
{
  bool isInputOk = (pedInfos.size() == 6 && mafs.size() != 0 && 
      selcoefs.size() == mafs.size() && snvPositions.size() == mafs.size());
  for (UINT i = 0; i != mafs.size(); ++i) {
    if (mafs[i] == 0.0) {
      std::cerr << "MAF" << i << " = 0!" << std::endl;
      isInputOk = false;
      break;
    }
  }

  if (isInputOk);
  else {
    std::cerr << "Error initializing gwPerson object. Please check input parameters. " << std::endl;
    exit(-1);
  }

  __genos.resize(2);
  __phenos.resize(0);
  __genoFreqs.resize(mafs.size());
  __pedInfos.resize(6);
   vectorF missings(mafs.size(), MISSING_ALLELE);
  __genos[0] = missings;
  __genos[1] = missings; 
  __pedInfos[0] = pedInfos[0];
  __pedInfos[1] = pedInfos[1];
  __pedInfos[2] = pedInfos[2];
  __pedInfos[3] = pedInfos[3];
  __pedInfos[4] = pedInfos[4];
  __pedInfos[5] = pedInfos[5];
  __mafs = mafs;

  for (UINT i = 0; i != __mafs.size(); ++i) {
    __genoFreqs[i].push_back( (1 - __mafs[i]) * (1 - __mafs[i]) );
    __genoFreqs[i].push_back( __mafs[i] * __mafs[i] );
  }

  __selcoefs = selcoefs;
  __snvPositions = snvPositions;
  //!- Initialize locus attributes based on selection coeficients
  for (UINT i = 0; i != mafs.size(); ++i) {

    if (selcoefs[i] > 1e-4) {
      if (mafs[i] <= BOUNDARY) 
        __locusAttributes.push_back(D_RV);
      else __locusAttributes.push_back(D_CV);
    } 
    else if (selcoefs[i] < -1e-4) {
      if (mafs[i] <= BOUNDARY) 
        __locusAttributes.push_back(P_RV);
      else __locusAttributes.push_back(P_CV);
    } 
    else {  
      if (mafs[i] <= BOUNDARY)
        __locusAttributes.push_back(SYNO_RV);
      else __locusAttributes.push_back(SYNO_CV);
    }
  }

  __isDebug = false;

  /*
  // Debug::
  std::cout << __pedInfos << std::endl;
  std::cout << __mafs << std::endl;
  std::cout << __selcoefs << std::endl;
  std::cout << __locusAttributes << std::endl;
  */
}


gwPerson::gwPerson( const vectorF& pedInfos, const vectorF& mafs, const vectorF& selcoefs, 
    const vectorUI& snvPositions, const vector2F& confounders, const vector2F& confounderProbs, 
    const vector2F& confounderEffects, const vectorL& areEpistatic )
{
  bool isInputOk = (pedInfos.size() == 6 && mafs.size() != 0 
      && selcoefs.size() == mafs.size() && snvPositions.size() == mafs.size() 
      && confounders.size() == 0 && confounderProbs.size() == confounders.size() 
      && areEpistatic.size() == confounders.size() && confounderEffects.size() == confounders.size());

  for (UINT i = 0; i != mafs.size(); ++i) {
    if (mafs[i] == 0.0) {
      std::cerr << "MAF" << i << " = 0!" << std::endl;
      isInputOk = false;
      break;
    }
  }

  if (isInputOk);
  else {
    std::cerr << "Error initializing gwPerson object. Please check input parameters. " << std::endl;
    exit(-1);
  }

  __genos.resize(2);
  __phenos.resize(0);
  __genoFreqs.resize(mafs.size());
  __pedInfos.resize(6);
  vectorF missings(mafs.size(), MISSING_ALLELE);
  __genos[0] = missings;
  __genos[1] = missings;
  __pedInfos[0] = pedInfos[0];
  __pedInfos[1] = pedInfos[1];
  __pedInfos[2] = pedInfos[2];
  __pedInfos[3] = pedInfos[3];
  __pedInfos[4] = pedInfos[4];
  __pedInfos[5] = pedInfos[5];
  __mafs = mafs;

  for (UINT i = 0; i != __mafs.size(); ++i) {
    __genoFreqs[i].push_back( (1 - __mafs[i]) * (1 - __mafs[i]) );
    __genoFreqs[i].push_back( __mafs[i] * __mafs[i] );
  }

  __selcoefs = selcoefs;
  __snvPositions = snvPositions;
  
  //!- Initialize locus attributes based on selection coeficients
  for (UINT i = 0; i != mafs.size(); ++i) {

    if (selcoefs[i] > 1e-4) {
      if (mafs[i] <= BOUNDARY) 
        __locusAttributes.push_back(D_RV);
      else __locusAttributes.push_back(D_CV);
    } 
    else if (selcoefs[i] < -1e-4) {
      if (mafs[i] <= BOUNDARY) 
        __locusAttributes.push_back(P_RV);
      else __locusAttributes.push_back(P_CV);
    } 
    else {  
      if (mafs[i] <= BOUNDARY)
        __locusAttributes.push_back(SYNO_RV);
      else __locusAttributes.push_back(SYNO_CV);
    }
  }

  __confounders = confounders;
  __confounderProbs = confounderProbs;
  __areEpistatic = areEpistatic;
  __confounderEffects = confounderEffects;

  __isDebug = false;
}


gwPerson::~gwPerson() {}


void gwPerson::setVerbose(int verbose) 
{
  if (verbose != 0)
    __isDebug = true;
  else __isDebug = false;
  return;
}


void gwPerson::generateGenotype(int byGenoFreqs, gsl_rng* gslr)
{
  switch (byGenoFreqs) {

    case 1 :
      {
        for (UINT i = 0; i != __genoFreqs.size(); ++i) {
          //!- Wild-type aa
          double a0 = __genoFreqs[i][0];
          //!- aa + Aa
          double a1 = 1 - __genoFreqs[i][1];

          double runif = gsl_rng_uniform(gslr);
          if (runif < a0) {
            __genos[0][i] = MAJOR_ALLELE; 
            __genos[1][i] = MAJOR_ALLELE;
          }
          else if (runif > a1) {
            __genos[0][i] = MINOR_ALLELE; 
            __genos[1][i] = MINOR_ALLELE;
          }
          else {
            if (runif < (a1 - a0) * 0.5 + a0) {
              __genos[0][i] = MINOR_ALLELE; 
              __genos[1][i] = MAJOR_ALLELE;
            } 
            else {
              __genos[0][i] = MAJOR_ALLELE; 
              __genos[1][i] = MINOR_ALLELE;
            }
          }
        }
      }
      break;

    case 2 :
      {
        for (UINT i = 0; i != __genoFreqs.size(); ++i) {
          //!- Wild-type aa
          double a0 = __genoFreqs[i][0];
          //!- aa + Aa
          double a1 = 1 - __genoFreqs[i][1];

          double runif = gsl_rng_uniform(gslr);
          if (runif < a0) {
            __genos[0][i] = MAJOR_ALLELE; 
            __genos[1][i] = MAJOR_ALLELE;
          }
          else if (runif > a1) {
            __genos[0][i] = MINOR_ALLELE; 
            __genos[1][i] = MINOR_ALLELE;
          }
          else {
            __genos[0][i] = MINOR_ALLELE; 
            __genos[1][i] = MAJOR_ALLELE;
          } 
        }
      }
      break;

    case 3 :
      {
        for (UINT i = 0; i != __genoFreqs.size(); ++i) {
          //!- Wild-type aa
          double a0 = __genoFreqs[i][0];
          //!- aa + Aa
          double a1 = 1 - __genoFreqs[i][1];

          double runif = gsl_rng_uniform(gslr);
          if (runif < a0) {
            __genos[0][i] = MAJOR_ALLELE; 
            __genos[1][i] = MAJOR_ALLELE;
          }
          else if (runif > a1) {
            __genos[0][i] = MINOR_ALLELE; 
            __genos[1][i] = MINOR_ALLELE;
          }
          else {
            __genos[0][i] = MAJOR_ALLELE; 
            __genos[1][i] = MINOR_ALLELE;
          } 
        }
      }
      break;

    default :
      {
        for (UINT i = 0; i != 2; ++i) {
          for (UINT j = 0; j != __mafs.size(); ++j) {
            double runif = gsl_rng_uniform(gslr); 

            // Debug::
            if (__isDebug)
              std::cout << "myunif " << runif << std::endl;

            if (runif < __mafs[j]) 
              __genos[i][j] = MINOR_ALLELE;  
            else 
              __genos[i][j] = MAJOR_ALLELE; 
          }
        }
      }
      break;
  }
  
  if (__isDebug)
    std::cout << __genos << std::endl;

  return;
}


void gwPerson::updateLocusAttributes(const vectorF& proportionsDeleteriousProtective, gsl_rng* gslr) 
{
  bool isInputOk = ( proportionsDeleteriousProtective[0] <= 1.0 
      && proportionsDeleteriousProtective[1] <= 1.0 
      && proportionsDeleteriousProtective[0] >= 0.0 
      && proportionsDeleteriousProtective[1] >= 0.0 ); 
  if (isInputOk) ;
  else {
    std::cerr << "Error input proportion of deleterious/protective variants. Please check input parameters. " << std::endl;
    exit(-1);
  }
  

  for (UINT i = 0; i != __locusAttributes.size(); ++i) {
    double runif = gsl_rng_uniform(gslr); 
    if (__locusAttributes[i] == D_RV) {
      if (runif >= proportionsDeleteriousProtective[0]) 
        __locusAttributes[i] = N_RV; //!- synonymous but neutral to the trait
      else;
    }
    else if (__locusAttributes[i] == P_RV) {
      if (runif >= proportionsDeleteriousProtective[1]) 
        __locusAttributes[i] = N_RV; //!- protective but not beneficial to the trait
      else;
    }
    else;
  }
  
  return;
}


void gwPerson::updateGenotypeFreqs(const vectorF& parsInput, bool isParConstant, char* moi) 
{
  double diseaseStatus = __pedInfos[5];
  //!- MAF of variants of interest (RV's only)
  vectorF rvWeights(0);
  vectorL shouldUpdateLocus(0);
  double par = 0.0;
  
  vectorF marginalOddsRatios(__mafs.size(), 1.0);
  
  //!- is a control, then have to update its MAF by protective variants
  if (diseaseStatus == UNAFFECTED) {
    par = parsInput[1];
    for (UINT i = 0; i != __mafs.size(); ++i) {
      if (__locusAttributes[i] == P_RV ) {
        rvWeights.push_back(1.0 / __mafs[i]);
        shouldUpdateLocus.push_back(true);
      }
      else
        shouldUpdateLocus.push_back(false);
    }
  }
  //!- is a case, then have to update its MAF by deleterious variants
  else if (diseaseStatus == AFFECTED) {
    par = parsInput[0];
    for (UINT i = 0; i != __mafs.size(); ++i) {
      if (__locusAttributes[i] == D_RV ) {
        rvWeights.push_back(1.0 / __mafs[i]);
        shouldUpdateLocus.push_back(true);
      }
      else
        shouldUpdateLocus.push_back(false);
    }
  }
  //!- Return if the phenotype is unknown or uncertain (not coded as 1 or 2) 
  else return;


  //!- Return if there is nothing to update here
  if (rvWeights.size() == 0 )
    return;
  
  double totalWeight = gw_sum(rvWeights);
  for (UINT i = 0; i != rvWeights.size(); ++i)
    rvWeights[i] = rvWeights[i] / totalWeight;

  //!- Marginal PAR for each variant
  vectorF marginalPars(0);
 
  //!- constant effect model
  if (isParConstant) { 
    for (UINT i = 0; i != shouldUpdateLocus.size(); ++i) {
      if (shouldUpdateLocus[i]) 
        marginalPars.push_back(par / (rvWeights.size() * 1.0));
      else
        marginalPars.push_back(0.0);
    }
  }

  //!- variable effects model
  else {
    for (UINT i = 0; i != shouldUpdateLocus.size(); ++i) {
      if (shouldUpdateLocus[i]) 
        marginalPars.push_back( (1.0 /__mafs[i]) / totalWeight * par );
      else
        marginalPars.push_back(0.0);
    }
  }

  if (__isDebug) {

    std::cout << "\nDisease Status and MOI\n" << diseaseStatus << *moi << std::endl;
    std::cout << "\nWeights at each site (only effective for variable PAR model)\n" << rvWeights << std::endl;
    std::cout << "\nLocus has a effective variant (a deleterious RV for cases or protective RV for controls)\n" << shouldUpdateLocus << std::endl;
    std::cout << "\nMarginal PARs at each site\n" << marginalPars << std::endl;
    std::cout << "\nSum of marginal PAR\n" << gw_sum(marginalPars) << std::endl;
    std::cout << "\nInput PAR {p0, p1} (p0 for deleterious RV, p1 for protective RV)\n" << parsInput << std::endl;
  }

  
  for (UINT i = 0; i != __mafs.size(); ++i) {
    // See p4 of Browning's paper
    if (shouldUpdateLocus[i]) { 
      double qu00 = __genoFreqs[i][0];
      double qu01 = 1.0 - __genoFreqs[i][0] - __genoFreqs[i][1];
      double qu11 = __genoFreqs[i][1];
      double qa00 = 0.0, qa01 = 0.0, qa11 = 0.0;
      double r11 = 0.0;
      switch (*moi) {
        case 'R' : 
          // Odds ratio = 1 for genotypes 00, 01. par_11 = mPar
          {
            qa00 = qu00; qa01 = qu01;
            r11 = marginalPars[i] / ((1.0 - marginalPars[i]) * qu11) + 1.0;
            qa11 = r11 * qu11 / (1.0 + (r11 - 1.0) * qu11);

            qa11 = qa11 / (qa00 + qa01 + qa11); 
            qa01 = qa01 / (qa00 + qa01 + qa11); 
            qa00 = qa00 / (qa00 + qa01 + qa11); 
            // Normalization
          }
          break;

        case 'D' :
          // Odds ratio = 1 for genotype 00.
          // OR_01 should be the same as OR_11.
          // Hence, par_01 ~ (qu_01/qu_11)*par_11
          // => par_01 = [qu_01 / (qu_01 + qu_11)] * mPar
          {
            qa00 = qu00;
            double mPar01 = qu01 / (qu01 + qu11) * marginalPars[i];
            double r01 = mPar01 / ((1.0 - mPar01) * qu01) + 1.0;
            r11 = r01;

            qa01 = r01 * qu01 / (1.0 + (r01 - 1.0) * qu01);
            qa11 = r11 * qu11 / (1.0 + (r11 - 1.0) * qu11);

            qa11 = qa11 / (qa00 + qa01 + qa11); 
            qa01 = qa01 / (qa00 + qa01 + qa11); 
            qa00 = qa00 / (qa00 + qa01 + qa11); 
            // Normalization
          }
          break;

        default :
          // additive
          // Odds ratio = 1 for genotype 00.
          // OR_11 should be two times OR_01 (dosage effect).
          // => par_01 = [qu_01 / (qu_01 + 2qu_11)] * mPar
          {
            qa00 = qu00;

            double mPar01 = qu01 / (qu01 + 2.0 * qu11) * marginalPars[i];
            double r01 = mPar01 / ((1.0 - mPar01) * qu01) + 1.0;
            r11 = 2.0 * r01;

            qa01 = r01 * qu01 / (1.0 + (r01 - 1.0) * qu01);
            qa11 = r11 * qu11 / (1.0 + (r11 - 1.0) * qu11);

            qa11 = qa11 / (qa00 + qa01 + qa11); 
            qa01 = qa01 / (qa00 + qa01 + qa11); 
            qa00 = qa00 / (qa00 + qa01 + qa11);
            // Normalization
          }
          break;
      }
      
      marginalOddsRatios[i] = r11;
      __genoFreqs[i][0] = qa00;
      __genoFreqs[i][1] = qa11;
    }
    else continue;
  }
 
  if (__isDebug) {
    std::cout << "\nMarginal ORs at each site\n" << marginalOddsRatios << std::endl;
    std::cout << "\nJoint OR\n" << gw_prod(marginalOddsRatios) << std::endl;
  }

  return ;
}


void gwPerson::updateGenotypeFreqs(const vector2F& genotypeFreqs) 
{
  bool isInputOk = (genotypeFreqs[0].size() == 2 && 
      genotypeFreqs.size() == __genoFreqs.size());
  if (isInputOk);
  else {
    std::cerr << "Error New genotypeFreqs of wrong size. Please check input genotypeFreqs. " << std::endl;
    exit(-1);
  }

  __genoFreqs = genotypeFreqs;
    
  return;
}


double gwPerson::m_calcLocusOddsRatio(double allele1, double allele2, int locusAttribute, 
    double maf, const vectorF& oddsRatios, double baselinef, char* moi) const
{
  double oddsRatio = 1.0;
  if (allele1 == MAJOR_ALLELE && allele2 == MAJOR_ALLELE) 
    //!- return if locus is wild-type
    return oddsRatio;

  if (locusAttribute == P_CV)
    //!- recode for switch statement. Does not make a different in effect
    locusAttribute = D_CV;

  switch (locusAttribute) {
    case D_RV : 
      {
        //!- constant effect for deleterious RV
        if (oddsRatios[0] == 0.0) 
          oddsRatio = oddsRatios[1];
        //!- variable effects for deleterious RV
        else {
          vectorF mafsRv(0);
          for (UINT j = 0; j != __mafs.size(); ++j) {
            if (__locusAttributes[j] == D_RV) 
              mafsRv.push_back(__mafs[j]);
          }
          double mafMin = *min_element(mafsRv.begin(), mafsRv.end());  // Compute highest freq 
          double mafMax = *max_element(mafsRv.begin(), mafsRv.end());  // Compute lowest freq
          oddsRatio = m_calcVariantEffect(maf, oddsRatios[0], oddsRatios[1], mafMin, mafMax, 1.0);
        }
      }
      break;

    case P_RV :
      {
        //!- constant effect for protective RV
        if (oddsRatios[2] == 0.0) 
          oddsRatio = oddsRatios[3];
        //!- variable effects for protective RV
        else {
          vectorF mafsRv(0);
          for (UINT j = 0; j != __mafs.size(); ++j) {
            if (__locusAttributes[j] == P_RV) 
              mafsRv.push_back(__mafs[j]);
          }
          double mafMin = *min_element(mafsRv.begin(), mafsRv.end());  // Compute highest freq 
          double mafMax = *max_element(mafsRv.begin(), mafsRv.end());  // Compute lowest freq 
          oddsRatio = m_calcVariantEffect(maf, oddsRatios[2], oddsRatios[3], mafMin, mafMax, 1.0);
        }
      }
      break;

    case D_CV : 
      {
        //!- constant effect for CV
        oddsRatio = oddsRatios[4];
      }
      break;

    default :
      return oddsRatio;
  }

  //!- Homozygote locus
  if (allele1 == MINOR_ALLELE && allele2 == MINOR_ALLELE) {
    double baselineOdds = m_calcOdds(baselinef);
    double f1 = m_calcPenetrance(oddsRatio * baselineOdds);

    //!- Additive model, f2 = 2*f1 - f0 ==> loc_Odds = f2/(1-f2) = (2f1 - f0) / (1 - 2f1 + f0)  ==>  rOR = (loc_Odds / bOdds)
    if (*moi == 'A') { 
      //!- Notice that for protective variants under additive model 2*f1 may be smaller than f0, have to adjust for this
      if (2 * f1 - baselinef <= 0) 
        oddsRatio = oddsRatio / 2.0;
      else
        oddsRatio = m_calcOdds( (2 * f1 - baselinef) ) / baselineOdds; 
    }
    //!- Multiplicative model, f2 = f1^3 / f0^2 ==> loc_Odds = f2/(1-f2) = (f1^3 / f0^2) / (1 - f1^3 / f0^2) ==> rOR = (loc_Odds / bOdds)
    else if (*moi == 'M') 
      oddsRatio = m_calcOdds( (pow(f1, 3.0) / pow(baselinef, 2.0)) ) / baselineOdds;  
    //!- no update needed for Recessive model and Dominant model for homozygous locus.
    else;
  }
  //!- Heterozygous locus
  else {
    //!- under Recessive model a heterozygous locus has no risk for disease
    if (*moi == 'R')
      oddsRatio = 1.0;
    //!- under non-recessive models do nothing to a heterozygous (keep the rOR as originally being assigned)
    else;
  }

  return oddsRatio;
}


double gwPerson::computeGenotypicEffect(const vectorF& oddsRatios, double baselinef, char* moi) const
{
  bool isInputOk = (oddsRatios.size() == 5 && (oddsRatios[0] >= 1.0 || oddsRatios[0] == 0.0) 
      && oddsRatios[0] < oddsRatios[1] && (oddsRatios[2] < 1.0 && oddsRatios[2] >= 0.0)
      && oddsRatios[2] < oddsRatios[3] && oddsRatios[3] <= 1.0 && oddsRatios[4] != 0.0
      && baselinef > 0.0);
  if (isInputOk) ;
  else {
    std::cerr << "Error input proportion of odds ratios/prevalence. Please check input parameters. " << std::endl;
    exit(-1);
  }

  vectorF locusOddsRatios(__mafs.size());
  for (UINT i = 0; i != __locusAttributes.size(); ++i) {
    double allele1 = __genos[0][i];
    double allele2 = __genos[1][i];
    double tmpOR = m_calcLocusOddsRatio(allele1, allele2, __locusAttributes[i], 
        __mafs[i], oddsRatios, baselinef, moi);
    locusOddsRatios[i] = tmpOR;
  }
  
  if (__isDebug)
    std::cout << locusOddsRatios << std::endl;

  return m_calcOdds(baselinef) * gw_prod(locusOddsRatios);
}


double gwPerson::computeGenotypicEffect(const vectorF& meanShifts) const
{
  bool isInputOk = (meanShifts.size() == 3 && meanShifts[0] >= 0.0 
      && meanShifts[1] >= meanShifts[0] && meanShifts[2] >= 0.0);
  if (isInputOk) ;
  else {
    std::cerr << "Error input QT mean shift parameters. Please check input parameters. " << std::endl;
    exit(-1);
  }
  

  double rvShift = 0.0;
  double cvShift = 0.0;
    
  for (UINT i = 0; i != __locusAttributes.size(); ++i) {
    if (__genos[0][i] == MAJOR_ALLELE && __genos[1][i] == MAJOR_ALLELE) 
      continue;

    if (__locusAttributes[i] == D_RV || __locusAttributes[i] == P_RV) {
      if (meanShifts[0] == 0.0)
        rvShift += meanShifts[1] * ((__locusAttributes[i] < 0) ? -1 : 1) * (__genos[0][i] + __genos[1][i]);
      else {
        vectorF mafsRv(0);
        for (UINT j = 0; j != __mafs.size(); ++j) {
          if (__locusAttributes[j] == D_RV || __locusAttributes[j] == P_RV) 
            mafsRv.push_back(__mafs[j]);
        }
        double mafMin = *min_element(mafsRv.begin(), mafsRv.end());  // Compute highest freq 
        double mafMax = *max_element(mafsRv.begin(), mafsRv.end());  // Compute lowest freq 
        rvShift += m_calcVariantEffect(__mafs[i], meanShifts[0], meanShifts[1], mafMin, mafMax, 0.0) 
          * ((__locusAttributes[i] < 0) ? -1 : 1) * (__genos[0][i] + __genos[1][i]);
      }
    }
    else if (__locusAttributes[i] == D_CV || __locusAttributes[i] == P_CV)
      cvShift += meanShifts[2] * ((__locusAttributes[i] < 0) ? -1 : 1) * (__genos[0][i] + __genos[1][i]);   
  }
    
  double totalShift = rvShift + cvShift; 
  return totalShift;
}
   

vectorUI gwPerson::getMendelianMafIdxes(double percentageCausal) const 
{
  bool isInputOk = (percentageCausal > 0.0 && percentageCausal <= 1.0);
  if (isInputOk) ;
  else {
    std::cerr << "Error. Please check input parameters. " << std::endl;
    exit(-1);
  }

  vectorF sortedMafs = __mafs;
  sort(sortedMafs.begin(), sortedMafs.end());
  UINT causalIdx = (UINT) (__mafs.size() * percentageCausal);
  double causalMaf = sortedMafs[causalIdx];
  // avoid pointing to common variants
  if (causalMaf > BOUNDARY ) {
    std::cerr << "Error. Causal frequencies too high (should not be higher than " << BOUNDARY << "). Please lower the percentageCausal. " << std::endl;
    exit(-1);
  }
  
  vectorUI causalIdxes(0);
  for (UINT i = 0; i != __mafs.size(); ++i) 
    if (__mafs[i] < causalMaf && __locusAttributes[i] == D_RV)
      causalIdxes.push_back(i);
  for (UINT i = 0; i != __mafs.size(); ++i)
    if (__mafs[i] == causalMaf && __locusAttributes[i] == D_RV && causalIdxes.size() <= causalIdx) 
      causalIdxes.push_back(i);
 
  if (__isDebug)
    std::cout << causalIdxes << std::endl;

  return causalIdxes;
}


void gwPerson::updateLocusAttributes(const vectorUI& mendelianMafIdxes, bool isAllelicHeterogeneous)
{
  if (mendelianMafIdxes.size() == 0)
    return;

  if (isAllelicHeterogeneous) {
    for (UINT i = 0; i != __locusAttributes.size(); ++i) {
      bool shouldSkip = false;
      for (UINT u = 0; u != mendelianMafIdxes.size(); ++u) {
        if (i == mendelianMafIdxes[u]) {
          shouldSkip = true;
          break;
        }
      }
      if (shouldSkip == false) {
        if (__mafs[i] <= BOUNDARY)
          __locusAttributes[i] = N_RV;
        else
          __locusAttributes[i] = N_CV;
      }
      else 
        continue;
    }
  }
  else {
    UINT tmp = mendelianMafIdxes[mendelianMafIdxes.size() - 1];
    for (UINT i = 0; i != __locusAttributes.size(); ++i) {
      if (i != tmp) {
        if (__mafs[i] <= BOUNDARY)
          __locusAttributes[i] = N_RV;
        else
          __locusAttributes[i] = N_CV;
      }
      else continue;
    }
  }
  return;
} 

/* too slow to use this strategy
double gwPerson::computeGenotypicEffect(const vectorUI& mendelianMafIdxes, bool isAllelicHeterogeneous, char* moi) const
{
  bool isInputOk = (*moi == 'D' || *moi == 'R');
  if (isInputOk) ;
  else {
    std::cerr << "Input MOI must be 'D' or 'R'. Please check input parameters. Now Quit." << std::endl;
    exit(-1);
  }

  double odds = 0.0;
  if (isAllelicHeterogeneous) {
    switch (*moi) {
      case 'D' :
        {
          for (UINT i = 0; i != mendelianMafIdxes.size(); ++i) {
            UINT tmp = mendelianMafIdxes[i];
            if (__genos[0][tmp] == MINOR_ALLELE || __genos[1][tmp] == MINOR_ALLELE) {
              odds = DBL_MAX;
              break;
            }
            else continue;              
          }
        }
        break;

      default : // recessive case
        {
          bool isMutation1 = false, isMutation2 = false;
          for (UINT i = 0; i != mendelianMafIdxes.size(); ++i) {
            UINT tmp = mendelianMafIdxes[i];
            if (__genos[0][tmp] == MINOR_ALLELE)
              isMutation1 = true;
            if (__genos[1][tmp] == MINOR_ALLELE) 
              isMutation2 = true;
          }
          if (isMutation1 && isMutation2) 
            odds = DBL_MAX;
        }
        break;
    }
  }

  else { // not allelic heterogeneous
    UINT tmp = mendelianMafIdxes[mendelianMafIdxes.size() - 1];
    switch (*moi) {
      case 'D' :
        {
          if (__genos[0][tmp] == MINOR_ALLELE || __genos[1][tmp] == MINOR_ALLELE) 
            odds = DBL_MAX;
        }
        break;

      default : // recessive case
        {
          if (__genos[0][tmp] == MINOR_ALLELE && __genos[1][tmp] == MINOR_ALLELE) 
            odds = DBL_MAX;
        }
        break;
    }
  }
  return odds; 
}
*/

//!\brief make Genotype frequencies in controls; mask out potential mutations for sites with small MAFs, by MOI 
void gwPerson::updateGenotypeFreqs(const vectorUI& mendelianMafIdxes, bool isAllelicHeterogeneous, char* moi) 
{
  if (__pedInfos[5] == AFFECTED)
    return;

  bool isInputOk = (*moi == 'D' || *moi == 'R' || *moi == 'C');
  if (isInputOk) ;
  else {
    std::cerr << "Input MOI must be 'D' or 'R' or 'C' (compound reccesive). Please check input parameters. Now Quit." << std::endl;
    exit(-1);
  }

  if (isAllelicHeterogeneous) {
    switch (*moi) {

      case 'R' : // recessive case
        {
          for (UINT u = 0; u != mendelianMafIdxes.size(); ++u) {
            __genoFreqs[mendelianMafIdxes[u]][0] = pow((1.0 - __mafs[mendelianMafIdxes[u]]), 2.0) / (2.0 * (1.0 - __mafs[mendelianMafIdxes[u]]) * __mafs[mendelianMafIdxes[u]] + pow((1.0 - __mafs[mendelianMafIdxes[u]]), 2.0));
            __genoFreqs[mendelianMafIdxes[u]][1] = 0.0;
          }
        }
        break;
      
      default : // dominant and compound recessive case
        {
          for (UINT u = 0; u != mendelianMafIdxes.size(); ++u) {
            __genoFreqs[mendelianMafIdxes[u]][0] = 1.0;
            __genoFreqs[mendelianMafIdxes[u]][1] = 0.0;
          }
        }
        break;
    }
  }

  else { // not allelic heterogeneous
    UINT mutationIdx = mendelianMafIdxes[mendelianMafIdxes.size() - 1];
    switch (*moi) {

      case 'R' : // recessive case
        {
          __genoFreqs[mutationIdx][0] = pow((1.0 - __mafs[mutationIdx]), 2.0) / (2.0 * (1.0 - __mafs[mutationIdx]) * __mafs[mutationIdx] + pow((1.0 - __mafs[mutationIdx]), 2.0));
          __genoFreqs[mutationIdx][1] = 0.0;
        }
        break;
      
      default : // dominant and compound recessive case 
        {
          __genoFreqs[mutationIdx][0] = 1.0;
          __genoFreqs[mutationIdx][1] = 0.0;
        }
        break;
    }
  }

  return ;
}


void gwPerson::updateGenotypeFreqs(const vectorUI& mendelianMafIdxes, bool isAllelicHeterogeneous, char* moi, gsl_rng* gslr) 
{
  if (__pedInfos[5] != AFFECTED) {
    std::cerr << "Input should have phenotype 'affected'. Now Quit." << std::endl;
    exit(-1);
  }

  bool isInputOk = (*moi == 'D' || *moi == 'R'|| *moi == 'C');
  if (isInputOk) ;
  else {
    std::cerr << "Input MOI must be 'D' or 'R' or 'C' (compound reccesive). Please check input parameters. Now Quit." << std::endl;
    exit(-1);
  }
  
  if (isAllelicHeterogeneous) {
    switch (*moi) {

      case 'R' : // recessive case
        {
          UINT tmpIdx = gsl_rng_uniform_int(gslr, mendelianMafIdxes.size());
          UINT mutationIdx = mendelianMafIdxes[tmpIdx];
          __genoFreqs[mutationIdx][0] = 0.0;
          __genoFreqs[mutationIdx][1] = 1.0;
        }
        break;
      
      default : // dominant and compound recessive case
        {
          UINT tmpIdx = gsl_rng_uniform_int(gslr, mendelianMafIdxes.size());
          UINT mutationIdx = mendelianMafIdxes[tmpIdx];
          __genoFreqs[mutationIdx][0] = 0.0; 
          // Normalize it
          __genoFreqs[mutationIdx][1] = pow(__mafs[mutationIdx], 2.0) / (2.0 *  __mafs[mutationIdx] * (1 -  __mafs[mutationIdx]) + pow(__mafs[mutationIdx], 2.0)); 
        }
        break;
    }
  }

  else { // not allelic heterogeneous
    UINT mutationIdx = mendelianMafIdxes[mendelianMafIdxes.size() - 1];
    switch (*moi) {

      case 'R' : // recessive case
        {
          __genoFreqs[mutationIdx][0] = 0.0;
          __genoFreqs[mutationIdx][1] = 1.0;
        }
        break;

      default : // dominant and compound recessive case
        {
          __genoFreqs[mutationIdx][0] = 0.0;
          __genoFreqs[mutationIdx][1] = pow(__mafs[mutationIdx], 2.0) / (2 *  __mafs[mutationIdx] * (1 -  __mafs[mutationIdx]) + pow(__mafs[mutationIdx], 2.0)); 
        }
        break;
    }
  }
 
  return ;
}

/*
double gwPerson::computeGenotypicEffect(const vectorUI& mendelianMafIdxes, bool isAllelicHeterogeneous, char* moi) const
{
  bool isInputOk = (*moi == 'D' || *moi == 'R');
  if (isInputOk) ;
  else {
    std::cerr << "Input MOI must be 'D' or 'R'. Please check input parameters. Now Quit." << std::endl;
    exit(-1);
  }

  double odds = 0.0;
  if (isAllelicHeterogeneous) {
    switch (*moi) {

      case 'R' : // recessive case
        {
          bool isMutation1 = false, isMutation2 = false;
          for (UINT i = 0; i != mendelianMafIdxes.size(); ++i) {
            UINT tmp = mendelianMafIdxes[i];
            if (__genos[0][tmp] == MINOR_ALLELE)
              isMutation1 = true;
            if (__genos[1][tmp] == MINOR_ALLELE) 
              isMutation2 = true;
          }
          if (isMutation1 && isMutation2) 
            odds = DBL_MAX;
        }
        break;

      default :
        {
          for (UINT i = 0; i != mendelianMafIdxes.size(); ++i) {
            UINT tmp = mendelianMafIdxes[i];
            if (__genos[0][tmp] == MINOR_ALLELE || __genos[1][tmp] == MINOR_ALLELE) {
              odds = DBL_MAX;
              break;
            }
            else continue;              
          }
        }
        break;
    }
  }

  else { // not allelic heterogeneous
    UINT tmp = mendelianMafIdxes[mendelianMafIdxes.size() - 1];
    switch (*moi) {

      case 'R' : // recessive case
        {
          if (__genos[0][tmp] == MINOR_ALLELE && __genos[1][tmp] == MINOR_ALLELE) 
            odds = DBL_MAX;
        }
        break;
        
      default :
        {
          if (__genos[0][tmp] == MINOR_ALLELE || __genos[1][tmp] == MINOR_ALLELE) 
            odds = DBL_MAX;
        }
        break;
    }
  }
  return odds; 
}
*/


void gwPerson::generatePhenotype(double genotypicEffect, bool isQt, gsl_rng* gslr)
{
  if (isQt) {
    double qt = gsl_ran_ugaussian(gslr);
    __pedInfos[5] = qt + genotypicEffect;
  }
  else {
    double penetrance = m_calcPenetrance(genotypicEffect);
    if (gsl_rng_uniform(gslr) <= penetrance) 
      __pedInfos[5] = AFFECTED;
    else 
      __pedInfos[5] = UNAFFECTED;
  }
  return;
}


vectorF gwPerson::getPhenotypes() const 
{
  if (__phenos.size() == 0) {
    std::cerr << "No additional phenotypes. Now quit." << std::endl;
    exit (-1);
  }
  
  if (__isDebug)
    std::cout << __phenos << std::endl;

  return __phenos;
}


double gwPerson::getPhenotype() const 
{
  return __pedInfos[5];
}


void gwPerson::updatePhenotype(double phenotype)
{
  __pedInfos[5] = phenotype;
  return;
}


void gwPerson::updateLocusAttributes(const vectorF& proportionsMissingData, bool shouldMarkMissing, gsl_rng* gslr) 
{

  bool isInputOk = (proportionsMissingData.size() == 4 
      && proportionsMissingData[0] >= 0.0 &&  proportionsMissingData[0] <= 1.0
      && proportionsMissingData[1] >= 0.0 &&  proportionsMissingData[1] <= 1.0
      && proportionsMissingData[2] >= 0.0 &&  proportionsMissingData[2] <= 1.0
      && proportionsMissingData[3] >= 0.0 &&  proportionsMissingData[3] <= 1.0);

  if (isInputOk) ;
  else {
    std::cerr << "Error input proportions of missing data. Please check input parameters. " << std::endl;
    exit(-1);
  }

  //!- If no missing data, skip. otherwise mark it
  if (gw_sum(proportionsMissingData) == 0.0)
    return;

  for (UINT i = 0; i != __locusAttributes.size(); ++i) {
    double runif = gsl_rng_uniform(gslr);

    if (__locusAttributes[i] == D_RV || __locusAttributes[i] == D_CV) {
      if (runif < proportionsMissingData[0]) 
        __locusAttributes[i] *= ((shouldMarkMissing) ? MARK_MISSING : MARK_WILD);
    }
    else if (__locusAttributes[i] == P_RV || __locusAttributes[i] == P_CV) {
      if (runif < proportionsMissingData[1]) 
        __locusAttributes[i] *= ((shouldMarkMissing) ? MARK_MISSING : MARK_WILD);
    }
    else if (__locusAttributes[i] == N_RV || __locusAttributes[i] == N_CV) {
      if (runif < proportionsMissingData[2]) 
        __locusAttributes[i] *= ((shouldMarkMissing) ? MARK_MISSING : MARK_WILD);
    }
    else if (__locusAttributes[i] == SYNO_RV || __locusAttributes[i] == SYNO_CV) {
      if (runif < proportionsMissingData[3]) 
        __locusAttributes[i] *= ((shouldMarkMissing) ? MARK_MISSING : MARK_WILD);
    }
    else
      std::cerr << "Error WARNING: __locusAttributes appears to have been updated. Please make sure you do not call 'updateLocusAttributes(const vectorF& proportionsMissingData, bool shouldMarkMissing, gsl_rng* gslr)' multiple times. " << std::endl;
  }
 
  return;
}


void gwPerson::updateLocusAttributes(const vectorI& locusAttributes)
{
  if (locusAttributes.size() != __locusAttributes.size()) {
    std::cerr << "Error New attributes of wrong length. Please check input locusAttributes. " << std::endl;
    exit(-1);
  }

  __locusAttributes = locusAttributes;
    
  return;
}


void gwPerson::updateLocusAttributes(const vectorL& shouldTrim)
{
  if (shouldTrim.size() != __locusAttributes.size()) {
    std::cerr << "Error Input vector of wrong length. Please check input 'shouldTrim'. " << std::endl;
    exit(-1);
  }

  for (UINT i = 0; i != __locusAttributes.size(); ++i) 
      __locusAttributes[i] *= ((shouldTrim[i]) ? MARK_MISSING : 1);

  return;
}


void gwPerson::summarizeLocusAttributes() const
{
  int nTotal = __locusAttributes.size(), nRare = 0, nProtective = 0, nDeleterious = 0, 
      nNeutral = 0, nSynonymous = 0, nMissing = 0;
  
  vectorI tmpAttributes = __locusAttributes;
  for (UINT i = 0; i != tmpAttributes.size(); ++i) {
    
    bool isMissing = false; 
    if (abs(tmpAttributes[i]) >= MARK_MISSING && tmpAttributes[i] % MARK_MISSING == 0) {
      tmpAttributes[i] = tmpAttributes[i] / MARK_MISSING;
      isMissing = true;
    }
    if (abs(tmpAttributes[i]) >= MARK_WILD && tmpAttributes[i] % MARK_WILD == 0) {
      tmpAttributes[i] = tmpAttributes[i] / MARK_WILD;
      isMissing = true;
    }
    if (isMissing)
      ++nMissing;
  }

  for (UINT i = 0; i != tmpAttributes.size(); ++i) {
    if (tmpAttributes[i] == D_RV || tmpAttributes[i] == P_RV ||
        tmpAttributes[i] == SYNO_RV || tmpAttributes[i] == N_RV )
      ++nRare;
    if (tmpAttributes[i] == P_RV || tmpAttributes[i] == P_CV )
      ++nProtective;
    if (tmpAttributes[i] == D_RV || tmpAttributes[i] == D_CV )
      ++nDeleterious;
    if (tmpAttributes[i] == N_RV || tmpAttributes[i] == N_CV )
      ++nNeutral;
    if (tmpAttributes[i] == SYNO_RV || tmpAttributes[i] == SYNO_CV )
      ++nSynonymous;
  }

  int nCommon = nTotal - nRare;

  std::cout << "total_variant_sites_No\t" << nTotal << std::endl; 
  std::cout << "rare_variant_sites_No\t" << nRare << std::endl; 
  std::cout << "common_variant_sites_No\t" << nCommon << std::endl;
  std::cout << std::endl; 
  std::cout << "synonymous_mutations_No\t" << nSynonymous << std::endl;
  std::cout << "deleterious_mutations_No\t" << nDeleterious << std::endl;
  std::cout << "protective_mutations_No\t" << nProtective << std::endl;
  std::cout << "neutral_mutations_No\t" << nNeutral << std::endl;
  std::cout << std::endl; 
  std::cout << "randomly_labelled_missing\t" << nMissing << std::endl;

  return;
}


void gwPerson::updateGenotype() 
{
  for (UINT i = 0; i != __locusAttributes.size(); ++i) {
    if ( __locusAttributes[i] % (MARK_WILD * MARK_MISSING) == 0 
        ||  __locusAttributes[i] % MARK_MISSING == 0) {
      __genos[0][i] = MISSING_ALLELE;
      __genos[1][i] = MISSING_ALLELE;
    }
    else if ( __locusAttributes[i] % MARK_WILD == 0) {
      __genos[0][i] = MAJOR_ALLELE;
      __genos[1][i] = MAJOR_ALLELE;
    }
    else;
  }
  
  return;
}


void gwPerson::updateGenotype(const vector2F& genotype) 
{  
  bool isInputOk = (genotype.size() == 2 && 
      genotype[0].size() == __genos[0].size());
  if (isInputOk);
  else {
    std::cerr << "Error New genotype of wrong size. Please check input genotype. " << std::endl;
    exit(-1);
  }
  __genos = genotype;
  return;
}


vectorI gwPerson::getLocusAttributes() const
{
  return __locusAttributes;
}


vectorUI gwPerson::getSnvPositions() const
{
  return __snvPositions;
}


vector2F gwPerson::getGenotype(bool isMissingTrimmed, bool isSynoTrimmed, 
    bool isCvTrimmed, vectorUI& untrimmedPositions) const 
{
  vector2F trimmedGenos(2);
  untrimmedPositions.resize(0);
  
  //!- !! Contains tricky codes below ...
  vectorI tmpAttributes = __locusAttributes;
  for (UINT i = 0; i != tmpAttributes.size(); ++i) {

    if (abs(tmpAttributes[i]) >= MARK_MISSING && tmpAttributes[i] % MARK_MISSING == 0)
      tmpAttributes[i] = tmpAttributes[i] / MARK_MISSING;
    if (abs(tmpAttributes[i]) >= MARK_WILD && tmpAttributes[i] % MARK_WILD == 0)
      tmpAttributes[i] = tmpAttributes[i] / MARK_WILD;
  }
  
  if (__isDebug)
    std::cout << "Recoded attributes: \n" << tmpAttributes << std::endl;
  
  for (UINT i = 0; i != tmpAttributes.size(); ++i) {

    //!- check/trim for missing
    if (isMissingTrimmed) {
      if ( __genos[0][i] == MISSING_ALLELE) 
        continue;
      else;
    }
    //!- check/trim for synonymous
    if (isSynoTrimmed) {
      if ( tmpAttributes[i] == SYNO_RV || tmpAttributes[i] == SYNO_CV )
        continue;
      else;
    }
    //!- check/trim for common variants
    if (isCvTrimmed) {
      if ( tmpAttributes[i] == D_CV || tmpAttributes[i] == P_CV 
          || tmpAttributes[i] == N_CV || tmpAttributes[i] == SYNO_CV )
        continue;
      else;
    }

    trimmedGenos[0].push_back(__genos[0][i]);
    trimmedGenos[1].push_back(__genos[1][i]);
    untrimmedPositions.push_back(__snvPositions[i]);
  }
  
  if (__isDebug)
    std::cout << __genos[0].size() - trimmedGenos[0].size() << " sites are trimmed. " << std::endl;
 
  if (trimmedGenos.size() == 0) {
    trimmedGenos.resize(1);
    trimmedGenos[0][0] = MAJOR_ALLELE;
    trimmedGenos[0][1] = MAJOR_ALLELE;
  }

  return trimmedGenos;
}


vector2F gwPerson::getGenotype() const 
{
  return __genos;
}


vector2F gwPerson::getGenotypeFreqs() const 
{
  return __genoFreqs;
}


void gwPerson::debug(int showWhat) const 
{
  m_printPunches(60);
  std::cout.precision(9);

  switch (showWhat) {
    case 1 : 
      std::cout << "__mafs:\n" << __mafs << "\n" << std::endl;
      break;
    case 2 : 
      std::cout << "__genoFreqs:\n" << __genoFreqs << "\n" << std::endl;
      break;
    case 3 : 
      {
        std::cout << "__geno_1:\n" << __genos[0] << "\n" << std::endl; 
        std::cout << "__geno_2:\n" << __genos[1] << "\n" << std::endl;
      }
      break;
    case 4 :
      std::cout << "__LocusAttributes:\n" << __locusAttributes << std::endl;
      break;
    default : 
      {
        std::cout << "__mafs:\n" << __mafs << "\n" << std::endl;
        std::cout << "__genoFreqs:\n" << __genoFreqs << "\n" << std::endl;
        std::cout << "__geno_1:\n" << __genos[0] << "\n" << std::endl; 
        std::cout << "__geno_2:\n" << __genos[1] << "\n" << std::endl;
        std::cout << "__snvPositions:\n" << __snvPositions << std::endl;
        std::cout << "__LocusAttributes:\n" << __locusAttributes << std::endl;
      }
      break;
  }
  
  m_printPunches(60);
  return;
}

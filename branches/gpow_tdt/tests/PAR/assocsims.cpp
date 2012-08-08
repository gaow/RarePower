//!\file assocsims.cpp
//!\brief simulator for disease association studies
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
#include <vector>
#include <string>
#include <cstdlib>
#include <algorithm>

#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"

#include "gw_utilities.h"
#include "person.h"
#include "assocsims.h"

namespace {
  const double AFFECTED = 2.0, UNAFFECTED = 1.0, UNPHENOTYPED = 0.0; 
}
  
//!\brief This program is the core simulator function that simulates complex traits via combined OR model, PAR model, QT model or Extreme QT models.
//use the "gwPerson" class

gwSimulator::gwSimulator(const vectorF& pedInfos, const vectorF& mafs, const vector2F& genoFreqs, 
        const vectorF& selcoefs, const vectorUI& positions)
{
  __pedInfos = pedInfos;
  __mafs = mafs;
  __genoFreqs = genoFreqs;
  __selcoefs = selcoefs;
  __positions = positions;
  __persons.resize(0);
  __mObservedMafs.resize(0);
}
gwSimulator::gwSimulator() {}
gwSimulator::~gwSimulator() {}

//!-# "pop-odds"
void gwSimulator::createGenotypeComplexTraitsAssociations (const vectorF& propFunctionalRv, double baselinef, 
    const vectorF& oddsRatios, char* moi, UINT nPopulation, gsl_rng* gslr) 
{
  //!- Initialize a "template" person object
  gwPerson* personTemplate = new gwPerson (__pedInfos, __mafs, __selcoefs, __positions);
  personTemplate->updateLocusAttributes(propFunctionalRv, gslr);

  bool isInputOk = (nPopulation > 0);
  if (isInputOk);
  else {
    std::cerr << "Population size not valid. Nothing to do. Quit now." << std::endl;
    exit(-1);
  }

  //!- actual prevalence, will be computed based on updated penetrance of each person in the population
  double prevalence = 0.0; 

  //!- Generate population with disease status, accept all generated samples
  for (UINT i = 0; i != nPopulation; ++i) {
    personTemplate->generateGenotype(0, gslr);
    double odds = personTemplate->computeGenotypicEffect(oddsRatios, baselinef, moi);
    personTemplate->generatePhenotype(odds, 0, gslr);
    __persons.push_back(*personTemplate);
    prevalence += odds / (1 + odds);
  }
  std::cerr << "The actual prevalence for cases is [ " << prevalence / (1.0 * nPopulation) << " ]. \n" << std::endl;
  //!- Output the disease prevalence for the binary trait
  //  std::cout << personTemplate->getPhenotype() << std::endl;
  delete personTemplate; 
  return;
}


//!-# "dichot-odds"
void gwSimulator::createGenotypeComplexTraitsAssociations (const vectorF& propFunctionalRv, double baselinef, 
    const vectorF& oddsRatios, char* moi, UINT nCases, UINT nCtrls, UINT nUnphenotyped, gsl_rng* gslr)
{
  //!- Initialize a "template" person object
  gwPerson* personTemplate = new gwPerson (__pedInfos, __mafs, __selcoefs, __positions);
  personTemplate->updateLocusAttributes(propFunctionalRv, gslr);

  bool isInputOk = (nCases > 0 || nCtrls > 0 || nUnphenotyped > 0);
  if (isInputOk);
  else {
    std::cerr << "Case/Control/Cohort size not valid. Nothing to do. Quit now." << std::endl;
    exit(-1);
  }

  //!- Generate case-ctrl samples (resample, reject or accept)
  UINT iCase = 0, iCtrl = 0, iUnphenotyped = 0;

  while (iCase != nCases || iCtrl != nCtrls) {
    personTemplate->generateGenotype(0, gslr);
    double odds = personTemplate->computeGenotypicEffect(oddsRatios, baselinef, moi);
    personTemplate->generatePhenotype(odds, 0, gslr);
    double trait = personTemplate->getPhenotype();

    if (trait == AFFECTED) {
      // get a case
      if (iCase != nCases) {  
        // if case is still in need, collect it.
        __persons.push_back(*personTemplate);
        ++iCase;
      }
      else;  
      // if case is enough, do nothing.
    }
    else {
      // get a control
      if (iCtrl != nCtrls) {  
        // if control is still in need, collect it.
        __persons.push_back(*personTemplate);
        ++iCtrl;
      }
      else;  
      // if control is enough, do nothing.          
    }
  }

  personTemplate->updatePhenotype(UNPHENOTYPED); 

  //!- Generate unphenotyped cohorts
  while (iUnphenotyped != nUnphenotyped) {
    personTemplate->generateGenotype(0, gslr);
    __persons.push_back(*personTemplate);
    ++iUnphenotyped; 
  }
  //  std::cout << personTemplate->getPhenotype() << std::endl;
  delete personTemplate;
  return;
}

 
//!-# "dichot-par"
void gwSimulator::createGenotypeComplexTraitsAssociations (const vectorF& propFunctionalRv, const vectorF& pars, bool isParConst, char* moi,
        UINT nCases, UINT nCtrls, UINT nUnphenotyped, gsl_rng* gslr)
{  
  //!- Initialize a "template" person object
  gwPerson* personTemplate = new gwPerson (__pedInfos, __mafs, __selcoefs, __positions);
  personTemplate->updateLocusAttributes(propFunctionalRv, gslr);
  
  bool isInputOk = (nCases > 0 || nCtrls > 0 || nUnphenotyped > 0);
  if (isInputOk);
  else {
    std::cerr << "Case/Control/Cohort size not valid. Nothing to do. Quit now." << std::endl;
    exit(-1);
  }

  //!- Generate case-ctrl samples. No rejecting because case/ctrl use different underlying MAF
  vector2F genoFreqsRecovery = __genoFreqs;
  personTemplate->updatePhenotype(AFFECTED); 
  for (UINT i = 0; i != nCases; ++i) {
    personTemplate->updateGenotypeFreqs(pars, isParConst, moi); 
    personTemplate->generateGenotype(1, gslr);
    __persons.push_back(*personTemplate);
    personTemplate->updateGenotypeFreqs(genoFreqsRecovery);
  }

  personTemplate->updatePhenotype(UNAFFECTED); 
  for (UINT i = 0; i != nCtrls; ++i) {
    personTemplate->updateGenotypeFreqs(pars, isParConst, moi); 
    personTemplate->generateGenotype(1, gslr);
    __persons.push_back(*personTemplate);
    personTemplate->updateGenotypeFreqs(genoFreqsRecovery);
  }

  //!- Generate unphenotyped cohorts. Use the MAF directly from SRV_batch

  personTemplate->updatePhenotype(UNPHENOTYPED); 
  for (UINT i = 0; i != nUnphenotyped; ++i) {
    personTemplate->generateGenotype(0, gslr);
    __persons.push_back(*personTemplate);
  }
  //  std::cout << personTemplate->getPhenotype() << std::endl;
  delete personTemplate;
  return;
}


//!-# "qt"
void gwSimulator::createGenotypeComplexTraitsAssociations (const vectorF& propFunctionalRv, const vectorF& qtcoefs, UINT nPopulation, gsl_rng* gslr)
{
  //!- Initialize a "template" person object
  gwPerson* personTemplate = new gwPerson (__pedInfos, __mafs, __selcoefs, __positions);
  personTemplate->updateLocusAttributes(propFunctionalRv, gslr);
  bool isInputOk = (nPopulation > 0);
  if (isInputOk);
  else {
    std::cerr << "Population size not valid. Nothing to do. Quit now." << std::endl;
    exit(-1);
  }

  //!- Generate population with QT values, accept all generated samples
  for (UINT i = 0; i != nPopulation; ++i) {
    personTemplate->generateGenotype(0, gslr);
    double meanShift = personTemplate->computeGenotypicEffect(qtcoefs);
    personTemplate->generatePhenotype(meanShift, 1, gslr);
    __persons.push_back(*personTemplate);
  }
  //  std::cout << personTemplate->getPhenotype() << std::endl;
  delete personTemplate;
  return;
}


//!-# "dichot-qt" (finite and infinite)
void gwSimulator::createGenotypeComplexTraitsAssociations (const vectorF& propFunctionalRv, const vectorF& qtcoefs, 
        const vectorF& qtcuts, bool shouldMarkCaseCtrl, UINT nPopulation, UINT nCases, UINT nCtrls, 
        UINT nUnphenotyped, gsl_rng* gslr) 
{
  //!- Initialize a "template" person object
  gwPerson* personTemplate = new gwPerson (__pedInfos, __mafs, __selcoefs, __positions);
  personTemplate->updateLocusAttributes(propFunctionalRv, gslr);

  //!- most tidious sampling scheme to program ... 
  //!- #0 = sample from infinite population, #!0 = sample extremes from finite cohort
  bool isInputOk = (qtcuts.size() == 2 && qtcuts[0] > 0.0 
      && qtcuts[1] < 1.0 && qtcuts[1] >= qtcuts[0]);
  if (isInputOk);
  else {
    std::cerr << "Extreme Qt lower/upper percentiles not valid. Quit now." << std::endl;
    exit(-1);
  }

  isInputOk = ((nCases > 0 || nCtrls > 0) && 
      (nPopulation == 0 || nPopulation >= (nCases + nCtrls)));
  if (isInputOk);
  else {
    std::cerr << "Case/Control/Cohort or population size not valid. Nothing to do. Quit now." << std::endl;
    exit(-1);
  }

  switch (nPopulation) {
    //!- extreme QT from unlimited population, using "qtcuts[0]" as percentile cut-off for ctrls, "qtcuts[1]" for cases
    //!- Will/Wont mark the phoentype as binary, determined by "shouldMarkCaseCtrl"
    case 0 :
      {
        UINT iCase = 0, iCtrl = 0;

        double lower = gsl_cdf_ugaussian_Pinv(qtcuts[0]);
        double upper = gsl_cdf_ugaussian_Pinv(qtcuts[1]);

        while (iCase != nCases || iCtrl != nCtrls) {
          personTemplate->generateGenotype(0, gslr);
          double meanShift = personTemplate->computeGenotypicEffect(qtcoefs);
          personTemplate->generatePhenotype(meanShift, 1, gslr);
          double trait = personTemplate->getPhenotype();

          if (trait >= upper) {
            // get a case
            if (iCase != nCases) {  
              // if case is still in need, collect it.
              if (shouldMarkCaseCtrl) 
                personTemplate->updatePhenotype(AFFECTED);

              __persons.push_back(*personTemplate);
              ++iCase;
            }
            else;  
            // if case is enough, do nothing.
          }

          else if (trait <= lower) {
            // get a control
            if (iCtrl != nCtrls) {  
              // if control is still in need, collect it.
              if (shouldMarkCaseCtrl) 
                personTemplate->updatePhenotype(UNAFFECTED);

              __persons.push_back(*personTemplate);
              ++iCtrl;
            }
            else;  
            // if control is enough, do nothing.          
          }

          else;
        }
      }
      break;

    default :
      {
        //!- extreme QT from given population, using "qtcuts[0]" as percentile cut-off for ctrls, "qtcuts[1]" for cases
        //!- Needs to convert bounds into index (floored by type conversion -- no need to worry about edge problem)
        UINT lower = (UINT) (nPopulation * qtcuts[0]);
        UINT upper = (UINT) (nPopulation * qtcuts[1]);
        if (upper == lower) ++upper;

        vectorF sortedPhenotypes(nPopulation);

        //!- First generate a cohort with QT known, as in "qt" model
        for (UINT i = 0; i != nPopulation; ++i) {
          personTemplate->generateGenotype(0, gslr);
          double meanShift = personTemplate->computeGenotypicEffect(qtcoefs);
          personTemplate->generatePhenotype(meanShift, 1, gslr);
          __persons.push_back(*personTemplate);
          sortedPhenotypes[i] = personTemplate->getPhenotype();
        }

        //!- Need to sort __mPhenotypes to figure out extreme values
        sort (sortedPhenotypes.begin(), sortedPhenotypes.end());

        //!- Collect only extreme samples
        //!- Will/Wont mark the phoentype as binary, determined by "shouldMarkCaseCtrl"
        std::vector<gwPerson>* __personsExtreme = new std::vector<gwPerson>();
        for (UINT i = 0; i != __persons.size(); ++i) {
          if (__persons[i].getPhenotype() >= sortedPhenotypes[upper]) {
            if (shouldMarkCaseCtrl) 
              __persons[i].updatePhenotype(AFFECTED);
            __personsExtreme->push_back(__persons[i]);
          }
          else if (__persons[i].getPhenotype() <= sortedPhenotypes[lower]) {
            if (shouldMarkCaseCtrl) 
              __persons[i].updatePhenotype(UNAFFECTED);
            __personsExtreme->push_back(__persons[i]);
          }
          else continue;
        }

        //!- Return to the extreme samples, a subset of the cohort.           
        __persons = *__personsExtreme;
        delete __personsExtreme;
      }
      break;   
  }
  
  //  std::cout << personTemplate->getPhenotype() << std::endl;
  delete personTemplate;

  return;
}


/*!\brief Generate Mendelian traits genotype and phenotype data*/
/*
void gwSimulator::creatGenotypeMendelianTraitsAssociations (double percentageCausal, bool isAllelicHeterogeneous,
        char* moi, UINT nCases, UINT nHeterCases, UINT nCtrls, gsl_rng* gslr)
{
  //!\brief his program is the core simulator function that simulates Mendelian traits by person's information (not pedigrees)



  //!- Initialize a "template" person object
  gwPerson* personTemplate = new gwPerson (__pedInfos, __mafs, __selcoefs, __positions);
  vectorUI mendelianMafIdxes = personTemplate->getMendelianMafIdxes(percentageCausal);
  personTemplate->updateLocusAttributes(mendelianMafIdxes, isAllelicHeterogeneous);

  //!- Generate case-ctrl samples (resample, reject or accept)
  UINT iCase = 0, iCtrl = 0, iHeterCases = 0;

  while (iCase != nCases || iCtrl != nCtrls || iHeterCases != nHeterCases) {
    personTemplate->generateGenotype(0, gslr);
    double odds = personTemplate->computeGenotypicEffect(mendelianMafIdxes, isAllelicHeterogeneous, moi);
    personTemplate->generatePhenotype(odds, 0, gslr);
    double trait = personTemplate->getPhenotype();

    if (trait == AFFECTED) {
      // get a case
      if (iCase != nCases) {  
        // if case is still in need, collect it.
        __persons.push_back(*personTemplate);
        ++iCase;
      }
      else;  
      // if case is enough, do nothing.
    }
    else {
      // get a control
      if (iCtrl != nCtrls) {  
        // if control is still in need, collect it.
        __persons.push_back(*personTemplate);
        ++iCtrl;
      }
      else if (iHeterCases != nHeterCases) {
        personTemplate->updatePhenotype(AFFECTED);
        __persons.push_back(*personTemplate);
        ++iHeterCases;
      }  
      // if HeterCases is needed collect it.          
      else;
    }
  }
  delete personTemplate;
  return;
}
*/
void gwSimulator::creatGenotypeMendelianTraitsAssociations (double percentageCausal, bool isAllelicHeterogeneous,
        char* moi, UINT nCases, UINT nHeterCases, UINT nCtrls, gsl_rng* gslr)
{  
  //!- Initialize a "template" person object
  gwPerson* personTemplate = new gwPerson (__pedInfos, __mafs, __selcoefs, __positions);
  vectorUI mendelianMafIdxes = personTemplate->getMendelianMafIdxes(percentageCausal);
  personTemplate->updateLocusAttributes(mendelianMafIdxes, isAllelicHeterogeneous);  

  bool isInputOk = (nCtrls > 0 || nCases > 0 || nHeterCases > 0);
  if (isInputOk);
  else {
    std::cerr << "Case/Control/'Heterogenious' Case size not valid. Quit now." << std::endl;
    exit(-1);
  }

  personTemplate->updateGenotypeFreqs(mendelianMafIdxes, isAllelicHeterogeneous, moi);
  vector2F genoFreqsRecovery = personTemplate->getGenotypeFreqs();


  //!- Generate case-ctrl samples. No rejecting because case/ctrl use different underlying MAF

  personTemplate->updatePhenotype(AFFECTED); 
  for (UINT i = 0; i != nCases; ++i) {
    if (*moi != 'C') {
      personTemplate->updateGenotypeFreqs(mendelianMafIdxes, isAllelicHeterogeneous, moi, gslr); 
      personTemplate->generateGenotype(1, gslr);
      __persons.push_back(*personTemplate);
      personTemplate->updateGenotypeFreqs(genoFreqsRecovery);
    }
    else {
      vector2F tmpGenotype(2);
      personTemplate->updateGenotypeFreqs(mendelianMafIdxes, isAllelicHeterogeneous, moi, gslr); 
      personTemplate->generateGenotype(2, gslr);
      tmpGenotype[0] = personTemplate->getGenotype()[0];
      personTemplate->updateGenotypeFreqs(genoFreqsRecovery);
      personTemplate->updateGenotypeFreqs(mendelianMafIdxes, isAllelicHeterogeneous, moi, gslr); 
      personTemplate->generateGenotype(3, gslr);
      tmpGenotype[1] = personTemplate->getGenotype()[1];
      personTemplate->updateGenotype(tmpGenotype);
      __persons.push_back(*personTemplate);
      personTemplate->updateGenotypeFreqs(genoFreqsRecovery);
    }
  }

  personTemplate->updatePhenotype(UNAFFECTED); 
  for (UINT i = 0; i != nCtrls; ++i) {
    personTemplate->generateGenotype(1, gslr);
    __persons.push_back(*personTemplate);
  }

  //!- Generate unphenotyped cohorts. Use the MAF directly from SRV_batch

  personTemplate->updatePhenotype(AFFECTED); 
  for (UINT i = 0; i != nHeterCases; ++i) {
    personTemplate->generateGenotype(1, gslr);
    __persons.push_back(*personTemplate);
  }
  //  std::cout << personTemplate->getPhenotype() << std::endl;
  delete personTemplate;
  return;
}


void gwSimulator::mimicGenotyping ( const vectorF& proportionsMissingData, bool shouldMarkMissing, 
        gsl_rng* gslr )
{
  //!\brief This program further edits the simulated sample by masking out some functional variants, etc.
  
  bool isInputOk = (__persons.size() != 0);
  //!- Check input: make sure the matrix is square
  //!- (otherwise cannot compute MAF because locus information incomplete)
  for (UINT i = 0; i != __persons.size(); ++i) {
    if (__persons[0].getLocusAttributes().size() != __persons[i].getGenotype()[0].size()) {
      isInputOk = false;
      break;
    }
  }

  if (isInputOk);
  else {
    std::cerr << "Input data-set not square or maf cut-offs not valid. Quit now." << std::endl;
    exit(-1);
  }

  //!- 1) Update locus attributes for each "person" object 2) Trim __mGenotypes 3) make data-set
  
  UINT nVariantSites = __persons[0].getGenotype()[0].size();
  vectorL shouldTrimSites(nVariantSites, false);
  
  for (UINT i = 0; i != __persons.size(); ++i) {
    __persons[i].updateLocusAttributes(proportionsMissingData, shouldMarkMissing, gslr);
    __persons[i].updateLocusAttributes(shouldTrimSites);
    __persons[i].updateGenotype();
  }

  return;
}


void gwSimulator::creatPedfileMatrix(bool isSynoTrimmed, bool isCvTrimmed, 
        bool isPedWritten, std::string projectName)
{
  //!\brief This program writes a ped file and make data matrix for assoc. analysis
  
  bool isInputOk = (__persons.size() != 0);
  //!- Check input: make sure the matrix is square. 
  //!- Even with missing data, to write a ped file the __mGenotypes should still be square
  for (UINT i = 0; i != __persons.size(); ++i) {
    if (__persons[0].getLocusAttributes().size() != __persons[i].getGenotype()[0].size()) {
      isInputOk = false;
      break;
    }
  }

  if (isInputOk);
  else {
    std::cerr << "Input data-set not square / require genotype and phenotype vectors be pre-allocated / MAFs size should be 0. Quit now." << std::endl;
    exit(-1);
  }

  __mGenotypes.resize(__persons.size()); 
  __mPhenotypes.resize(__persons.size());

  vectorUI snvPositions(0);
  vector2F tmpGeno = __persons[0].getGenotype(false, isSynoTrimmed, isCvTrimmed, snvPositions);
  UINT  nVariantSites = tmpGeno[0].size();
  vectorUI tmpSnvPositions(0);
  __mObservedMafs.resize( nVariantSites, 0.0);

  //!- Make the data-matrix 
  //!- Calculate observed MAFs while making the matrix
  for (UINT i = 0; i != __persons.size(); ++i) {
    
    __mPhenotypes[i] = __persons[i].getPhenotype();
    if (__mPhenotypes[i] == UNPHENOTYPED)
      __mPhenotypes[i] = UNAFFECTED;

    tmpGeno = __persons[i].getGenotype(false, isSynoTrimmed, isCvTrimmed, tmpSnvPositions);

    if (tmpGeno[0].size() != nVariantSites || tmpSnvPositions.size() != snvPositions.size()) {
      std::cerr << "Data-set should be square after trimmings of CV and synonymous sites. Quit now." << std::endl;
      exit(-1);
    }

    for (UINT j = 0; j != nVariantSites; ++j) {
      if (tmpGeno[0][j] == -9.0 || tmpGeno[1][j] == -9.0)
        __mGenotypes[i].push_back(-9.0);
      else {
        __mGenotypes[i].push_back( tmpGeno[0][j] + tmpGeno[1][j] );
        __mObservedMafs[j] += __mGenotypes[i][j];
      }
    }
  }
  
  for (UINT j = 0; j != nVariantSites; ++j)
    __mObservedMafs[j] = __mObservedMafs[j] / (1.0 * __persons.size());


  //!- Write the data-matrix
  //!- 3 files to write: PED, MAP and LOG
  if (isPedWritten) {
    std::string pedFileName = projectName;
    pedFileName.append(".ped");
    std::string pedFileNameLog = projectName;
    pedFileNameLog.append(".log");    
    std::string mapFileName = projectName;
    mapFileName.append(".map");
    std::ofstream fout(pedFileName.c_str());
    std::ofstream lout(pedFileNameLog.c_str());
    std::ofstream mout(mapFileName.c_str());
    fout.precision(5);
    lout.precision(5);

    for (UINT i = 0; i != __persons.size(); ++i) {
      for (UINT j = 0; j != 5; ++j) fout << 0 << " ";
      fout << __mPhenotypes[i] << " ";
      fout << __mGenotypes[i] << std::endl;
    }
    fout.close();
    std::cout << "\tPED file written [ " << pedFileName <<  " ]." << std::endl;
    
    for (UINT i = 0; i != snvPositions.size(); ++i) 
      mout << "SNV-" << i + 1 << " " << snvPositions[i] << std::endl;
    mout.close();
    std::cout << "\tMAP file written [ " << mapFileName <<  " ]." << std::endl;
    
    lout << pedFileName << " ** log file ** " << std::endl;
    lout << "Sample size:\n" << __persons.size() << std::endl;
    lout << "Length of variant sites:\n" << nVariantSites << std::endl;
    lout << "Are synonymous sites included in the data? (1 = Yes, 0 = No):\n" << (!isSynoTrimmed) << std::endl;
    lout << "Are (underlying) common variant sites included in the data? (1 = Yes, 0 = No):\n" << (!isCvTrimmed) << std::endl;
    lout << "Observed Locus MAF in the sample data:\n" << __mObservedMafs << std::endl;
    lout.close();
    std::cout << "\tLog file written [ " << pedFileNameLog <<  " ]." << std::endl;
  }
  return;
}

vector2F gwSimulator::getGenotypes() const 
{
  return  __mGenotypes;
}
vectorF gwSimulator::getPhenotypes() const
{
  return __mPhenotypes;
}
vectorF gwSimulator::getObservedMafs() const
{
  return __mObservedMafs;
}

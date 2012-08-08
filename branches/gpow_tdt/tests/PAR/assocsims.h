//!\file assocsims.h
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

#ifndef ASSOCSIMS_H
#define ASSOCSIMS_H

#include "gw_utilities.h"
//!\brief simulation of data for association tests

class gwSimulator {

  public:
    gwSimulator();

    gwSimulator(const vectorF& pedInfos, const vectorF& mafs, const vector2F& genoFreqs, 
        const vectorF& selcoefs, const vectorUI& positions);

    //!\fn createGenotypeComplexTraitsAssociations(...)
    //!\fn create_genotype_phenotype_mendelian(...)
    //!\fn mimic_genotyping(...)
    //!\fn create_ped_data_matrix(...)

    /*!\brief Generate complex traits genotype and phenotype data 
     * \param persons i.e. object to create: std::vector<gwPerson>& persons
     * \param pedInfos
     * \param mafs
     * \param genoFreqs
     * \param selcoefs
     * \param positions
     * \param propFunctionalRv 
     * \param simulationTask
     * \param oddsRatios
     * \param baselinef
     * \param qtcoefs
     * \param qtcuts
     * \param shouldMarkCaseCtrl
     * \param pars 
     * \param isParConst 
     * \param moi
     * \param nPopulation
     * \param nCases
     * \param nCtrls
     * \param nUnphenotyped
     * \param gslr gsl_rng* object, initialize GNU/GSL rng
     * */
    ~gwSimulator();
    void createGenotypeComplexTraitsAssociations (const vectorF& propFunctionalRv, double baselinef, const vectorF& oddsRatios,  
        char* moi, UINT nPopulation, gsl_rng* gslr);

    void createGenotypeComplexTraitsAssociations (const vectorF& propFunctionalRv, double baselinef, const vectorF& oddsRatios,  
        char* moi, UINT nCases, UINT nCtrls, UINT nUnphenotyped, gsl_rng* gslr);

    void createGenotypeComplexTraitsAssociations (const vectorF& propFunctionalRv, const vectorF& pars, bool isParConst, char* moi,
        UINT nCases, UINT nCtrls, UINT nUnphenotyped, gsl_rng* gslr);

    void createGenotypeComplexTraitsAssociations (const vectorF& propFunctionalRv, const vectorF& qtcoefs, UINT nPopulation, gsl_rng* gslr);

    void createGenotypeComplexTraitsAssociations (const vectorF& propFunctionalRv, const vectorF& qtcoefs, 
        const vectorF& qtcuts, bool shouldMarkCaseCtrl, UINT nPopulation, UINT nCases, UINT nCtrls, 
        UINT nUnphenotyped, gsl_rng* gslr);

    /*!\brief Generate Mendelian traits genotype and phenotype data*/
    void creatGenotypeMendelianTraitsAssociations (double percentageCausal, bool isAllelicHeterogeneous,
        char* moi, UINT nCases, UINT nHeterCases, UINT nCtrls, gsl_rng* gslr);


    /*!\brief Add noise to genotypes and/or edit genotypes. Involves trimming genotypes, marking missing values, etc. 
     * \param proportionsMissingData	Portion of sites that has missing data, a vector of length 4: {deleterious%, protective%, non-causal%, synonymous%}
     * \param shouldMarkMissing	Logical, if == true then "MARK_MISSING" (multiply the attribute by 1000) for future coding genotype as "MISSING_GENO(-9)" 
     * (otherwise will be coded "MAJOR_ALLELE(0)")
     * \param gslr gsl_rng* object, initialize GNU/GSL rng
     * */
    void mimicGenotyping ( const vectorF& proportionsMissingData, bool shouldMarkMissing, 
        gsl_rng* gslr );


    /*!\brief Create the ped file as output and/or generate data matrix for assoc. analysis. 
     * \param isSynoTrimmed = true
     * \param isCvTrimmed = false
     * \param isPedWritten Whether or not to produce a ped file
     * \param projectName ped filename = projectName.ped
     * */
    void creatPedfileMatrix(bool isSynoTrimmed, bool isCvTrimmed, 
        bool isPedWritten, std::string projectName); 

    vector2F getGenotypes() const;
    vectorF getPhenotypes() const;
    vectorF getObservedMafs() const;


  private:
    //!\brief persons i.e. object to edit. std::vector<gwPerson>& persons
    std::vector<gwPerson> __persons;

    //!\brief 2D object, underlying genotype frequencies    
    vector2F __genoFreqs;
    //!\brief The first 6 columns of a *.ped file 
    vectorF __pedInfos;    
    //!\brief Additional phenotypes, if any
    vectorF __phenos;
    //!\brief Simulated current population MAF
    vectorF __mafs;
    //!\brief Selection coeffient on each locus of the haplotype of simulated current generation
    vectorF __selcoefs;    
    //!\brief positions in the gene region
    vectorUI __positions;
    //!\brief genotypes 2D vector of genotype matrix (to be modified)
    vector2F __mGenotypes;
    //!\brief phenotypes 1D vector of phenotypes (to be modified)
    vectorF __mPhenotypes; 
    //!\brief observedMafs 1D vector of observed MAF for each locus (to be modified)
    vectorF __mObservedMafs;
};
#endif///:~

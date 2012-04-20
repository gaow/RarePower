//: test_person.cpp
// test program
// Gao Wang 2011


#include <iostream>
#include <vector>

#include "gw_utilities.h"
#include "person.h"
using namespace std;

int main()
{
  
  RNG rng;
  gsl_rng* gslr = rng.get();
/*  
  vector<double> pedInfos(6, 1.0);
  pedInfos[5] += 1.0;
  string mafFile("Boyko2008European1p5k.maf");
  string selFile("Boyko2008European1p5k.sel");
  string posFile("Boyko2008European1p5k.pos");
  vector < vector < double > > FreqDat;
  vector < vector < double > > SelDat;
  vector < vector < unsigned int > > PosDat;
  scan_vector2F(selFile, SelDat);
  scan_vector2F(mafFile, FreqDat);
  scan_vector2UI(posFile, PosDat);
  unsigned int Pick = gsl_rng_uniform_int(gslr, FreqDat.size());
  vector<double>& mafs =  FreqDat[Pick];
  vector<double>& selcoefs = SelDat[Pick];
  vector<unsigned int>& positions = PosDat[Pick];
  vector<double> rEitArgDv(2);
  rEitArgDv[0] = 1.0;
  //!- proportion of effective causal variant
  rEitArgDv[1] = 0.8;
  //!- proportion of effective protective variant
  vector<double> rOddsArgDv(5);
  rOddsArgDv[0] = 2.0;
  rOddsArgDv[1] = 4.0;
  //!- odds ratio for causal variants
  rOddsArgDv[2] = 0.2;
  rOddsArgDv[3] = 0.8;
  //!- odds ratio for protective variants
  rOddsArgDv[4] = 1.2;
  //!- odds ratio for common variants
  vector<double> qtbeta(3);
  //!- Locus specific effects w.r.t. standard deviation
  qtbeta[0] = 0.1;
  qtbeta[1] = 0.6;
	//!- effect per rare variant
  qtbeta[2] = 0.05;
	//!- effect per common variant
  vector<double> rParDv(2);
  rParDv[0] = 0.1;
  rParDv[1] = 0.1;
  char sMoh[] = "A"; 
  char* moi = sMoh;
  vector<double> exclPort(4);
  exclPort[0] = 0.2;
  exclPort[1] = 0.2;
  exclPort[2] = 0.2;
  exclPort[3] = 0.2;
  
  vector2F genoFreqs(mafs.size());

  for (UINT i = 0; i != mafs.size(); ++i) {
    genoFreqs[i].push_back( (1 - mafs[i]) * (1 - mafs[i]) );
    genoFreqs[i].push_back( mafs[i] * mafs[i] );
  }
  
  cout << "BEGIN initialization tests" << endl;
  // test 1
  gwPerson person;
 
  // test 2
  person.debug(9); 

  // test 3
  person = gwPerson (pedInfos, mafs, selcoefs, positions);
  person.setVerbose(1);
  
  person.debug(9);
 
 
  // test 4
  cout << "BEGIN updateGenotype tests" << endl;
  person.generateGenotype(0, gslr);
  person.debug(9);
 // person.debug(1);
 // person.debug(2);
  
  //test 5
//  person.generateGenotype(1, gslr);
//  person.debug(1);
//  person.debug(2);
//  person.debug(9);
  
  // test 6
  cout << "BEGIN updateLocusAttributs tests" << endl;
  person.updateLocusAttributes(rEitArgDv, gslr);
  person.debug(4);
 
  // test 7
  cout << "BEGIN PAR models tests" << endl;
  person.updatePhenotype(2.0);
  person.updateGenotypeFreqs(rParDv, false, moi); 
  person.debug(1);
//  person.debug(2);
  person.updateGenotypeFreqs(genoFreqs);
  person.updateGenotypeFreqs(rParDv, true, moi); 
  person.debug(1);
//  person.debug(2);
  person.generateGenotype(1, gslr);
  
  // test 8
  cout << "BEGIN OR/Beta models tests" << endl;
  
  person = gwPerson (pedInfos, mafs, selcoefs, positions);
  //person.setVerbose(1);
  person.generateGenotype(0, gslr);
  person.updateLocusAttributes(rEitArgDv, gslr);
  gwPerson person4 = person;
  double odds = person.computeGenotypicEffect(rOddsArgDv, 0.01, moi);
  double beta = person4.computeGenotypicEffect(qtbeta);
  cout << "Odds: " << odds << endl;
  cout << "beta: " << beta << endl;
  
  // test 9
  cout << "BEGIN Get phenotype tests" << endl;
  person4.generatePhenotype(beta, 1, gslr);
  person.generatePhenotype(odds, 0, gslr);
  vector<double> trait = person.getPhenotypes();
  vector<double> trait4 = person4.getPhenotypes();
  cout << trait << endl;
  cout << trait4 << endl;

  // test 10
  cout << "BEGIN Locus Attributes tests" << endl;
  person.updateLocusAttributes(exclPort, false, gslr);
 

  // test 11
  vector<bool> shouldTrim(mafs.size(), false);
  unsigned int site = gsl_rng_uniform_int(gslr, shouldTrim.size());
  shouldTrim[site] = true;
  person.updateLocusAttributes(shouldTrim);
  person.updateGenotype();

  person.setVerbose(1);
  vector< vector<double> > mygeno;
  vectorUI untrimmedPositions(0);
  mygeno = person.getGenotype(true, true, true, untrimmedPositions);
  cout << mygeno[0] << endl;
  cout << mygeno[1] << endl;
  cout << untrimmedPositions << endl;
  mygeno = person.getGenotype(true, true, false, untrimmedPositions);
  cout << mygeno[0] << endl;
  cout << mygeno[1] << endl;
  cout << untrimmedPositions << endl;
  mygeno = person.getGenotype(true, false, true, untrimmedPositions);
  cout << mygeno[0] << endl;
  cout << mygeno[1] << endl;
  cout << untrimmedPositions << endl;
  mygeno = person.getGenotype(true, false, false, untrimmedPositions);
  cout << mygeno[0] << endl;
  cout << mygeno[1] << endl;
  cout << untrimmedPositions << endl;
  person.summarizeLocusAttributes();
  cout << person.getLocusAttributes() << endl;
  cout << "END of program" << endl;
  //person.debug(4);
  cout << Pick << endl;
  */


  return 0;
}

#include <iostream>
#include <vector>

#include "gw_utilities.h"
#include "person.h"
using namespace std;

int main()
{

  RNG rng;
  gsl_rng* gslr = rng.get(10086);
  string mafFile("Boyko2008European1p5k.maf");
  string selFile("Boyko2008European1p5k.sel");
  string posFile("Boyko2008European1p5k.pos");
  vector2F mafDat;
  vector2F selDat;
  vector2UI posDat;

  scan_vector2F(selFile, selDat);
  scan_vector2F(mafFile, mafDat);
  scan_vector2UI(posFile, posDat);

  UINT dataIdx = gsl_rng_uniform_int(gslr, mafDat.size());
  vector<double>& mafs =  mafDat[dataIdx];
  vector<double>& selcoefs = selDat[dataIdx];
  vector<UINT>& positions = posDat[dataIdx];
  vector2F genoFreqs(mafs.size());

  for (UINT i = 0; i != mafs.size(); ++i) {
    genoFreqs[i].push_back( (1 - mafs[i]) * (1 - mafs[i]) );
    genoFreqs[i].push_back( mafs[i] * mafs[i] );
  }

  //! population attributable risk model
  vector<double> pars(2);
  pars[0] = 0.05;
  pars[1] = 0.0;
  bool isParConst = true;
  char sMoi[] = "A"; 
  char* moi = sMoi;
  vector<double> pedInfos(6, 0.0);
  
  gwPerson person;
  person.debug(9); 
  person = gwPerson (pedInfos, mafs, selcoefs, positions);
  person.setVerbose(1);
  person.debug(9);

  cout << "BEGIN PAR models tests" << endl;
  person.updatePhenotype(2.0);
  person.updateGenotypeFreqs(pars, isParConst, moi); 
  person.debug(1);
  person.debug(2);
  person.updateGenotypeFreqs(genoFreqs);
//  person.updateGenotypeFreqs(pars, false, moi); 
//  person.debug(1);
//  person.debug(2);

    return 0;
}

#include<iostream>
#include<vector>
#include<string>
#include "gw_maths.h"
#include "gw_utilities.h"
#include "assoctests.h"

using namespace std;

int main() 
{
  vector< vector<double> > pedData(0);
  scan_vector2F("NewProject.ped", pedData);
  vector<double> ydat(pedData.size());
  vector< vector<double> > xdat(pedData.size());
  vector< vector<double> > observedMafs(0);
  scan_vector2F("NewProject.log", observedMafs);

  for (unsigned int i = 0; i != pedData.size(); ++i) {
    ydat[i] = pedData[i][5];
    for (unsigned int j = 6; j != pedData[i].size(); ++j) {
      xdat[i].push_back(pedData[i][j]);
    }
  }

//  cout << ydat << endl;
//  cout << xdat << endl;
  cout << ydat.size() << endl;
  cout << xdat.size() << endl;
  cout << xdat[0].size() << endl;
  cout << observedMafs.size() << endl;
  cout << observedMafs[0].size() << endl;
  cout << observedMafs << endl;
  
  double mafLower = 0.0, mafUpper = 0.01;
  double alpha = 0.05;

  gwAssociations associationTest(observedMafs[0], ydat, xdat, mafLower, mafUpper, alpha);

  double pvalue = associationTest.calcRvefisherP();
  cout << pvalue << endl;

  return 0;
}

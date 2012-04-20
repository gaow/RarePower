#include<iostream>
#include<vector>
#include "gw_maths.h"
using namespace std;

int main() 
{
  vector<int> mydat(4);
  mydat[0] = 5;
  mydat[1] = 2;
  mydat[2] = 18;
  mydat[3] = 32;
  cout << fexact_two_sided_pvalue(mydat) << endl;
  return 0;
}

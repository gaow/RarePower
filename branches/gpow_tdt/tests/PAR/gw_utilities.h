//!\file gw_utilities.h
//!\brief useful utilities
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

#ifndef GWUTILITY_H
#define GWUTILITY_H

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include "gsl/gsl_rng.h"

typedef std::vector<double> vectorF;
typedef std::vector< std::vector<double> > vector2F;
typedef std::vector<int> vectorI;
typedef std::vector<unsigned int> vectorUI;
typedef std::vector< std::vector<unsigned int> > vector2UI;
typedef std::vector<bool> vectorL;
typedef unsigned int UINT;

//!\brief setup GNU/RNG
class RNG {

  public:
  //!- the second-generation ranlux generators have the strongest proof of randomness. 
  //!- However turns out gsl_rng_alloc(gsl_rng_ranlxs2) is buggy !!   
  //gsl_rng_mt19937
    RNG() : rng( gsl_rng_alloc(gsl_rng_mt19937) ) {};
    ~RNG() { gsl_rng_free( rng ); }
    gsl_rng* get()
    {  
      // time(NULL): number of seconds since 00:00:00 GMT Jan. 1, 1970
      __seed = static_cast<unsigned long>(time (NULL) + getpid());
      gsl_rng_set(rng, __seed);
      return rng; 
    }    
    gsl_rng* get(const unsigned long seed)
    {  
      __seed = seed;
      gsl_rng_set(rng, __seed);
      return rng; 
    }

  private:
    gsl_rng* rng;
    unsigned long __seed;
};


namespace std {

//!- Dump a vector to screen
template<class T> ostream & operator<<(ostream & out, const vector<T> & vec)
  {
    if (!vec.empty()) {
      typename vector<T>::const_iterator it = vec.begin();
      out << *it;
      for (++it; it != vec.end(); ++it)
        out << " " << *it ;
    }
    return out;
  }
}

//!- scan file into vector2X
void scan_vector2F(std::string filename, vector2F&);
void scan_vector2UI(std::string filename, vector2UI&);
#endif ///:~

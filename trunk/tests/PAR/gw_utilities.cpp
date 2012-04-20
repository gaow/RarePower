//!\file gw_utilities.cpp
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

#include "gw_utilities.h"

void scan_vector2F(std::string filename, vector2F& info) 
{
  //!- scan file into vector2F
  std::cout << "\tSource data: "<< filename << std::endl;
  std::ifstream file;
  file.open(filename.c_str());	
  if (!file) {
    std::cerr << "\tUnable to source [ " << filename << " ]. Quit." << std::endl;
    exit(-1); // terminate with error
  }
  std::string line;

  while ( getline(file, line) ) {
    vectorF data;
    double value;
    std::istringstream iss(line);
    while (iss >> value) data.push_back(value);
    info.push_back(data);
  }
  return;
}
void scan_vector2UI(std::string filename, vector2UI& info) 
{
  //!- scan file into vector2F
  std::cout << "\tSource data: "<< filename << std::endl;
  std::ifstream file;
  file.open(filename.c_str());	
  if (!file) {
    std::cerr << "\tUnable to source [ " << filename << " ]. Quit." << std::endl;
    exit(-1); // terminate with error
  }
  std::string line;

  while ( getline(file, line) ) {
    vectorUI data;
    UINT value;
    std::istringstream iss(line);
    while (iss >> value) data.push_back(value);
    info.push_back(data);
  }
  return;
}

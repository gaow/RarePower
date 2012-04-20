// =====================================================================================
//
//       Filename:  gw_utilities.cpp
//
//    Description:  useful utilities
//
//        Version:  1.0
//        Created:  01/04/2011 10:50:24 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Gao Wang (gw), wangow@gmail.com
//                  Baylor College of Medicine, Texas, USA
//        License:  GNU General Public License <http://www.gnu.org/licenses/>
//                  Copyright (c) 2011, Gao Wang
//
// =====================================================================================


#include "gw_utilities.h"
bool is_file_empty(std::string filename)
{
	std::ifstream file;

	file.open(filename.c_str());
	if (!file)
		return true;

	std::string line;
	UINT length = 0;
	while (getline(file, line) )
		++length;
	file.close();
	if (length > 0)
		return false;
	else
		return true;
}


void scan_vector2F(std::string filename, vector2F & info)
{
	//!- scan file into vector2F
	std::ifstream file;

	file.open(filename.c_str());
	if (!file) {
		std::cerr << "ERROR: Unable to source [ " << filename << " ]. Quit." << std::endl;
		exit(-1); // terminate with error
	}
	std::string line;

	while (getline(file, line) ) {
		vectorF data;
		double value;
		std::istringstream iss(line);
		while (iss >> value) data.push_back(value);
		info.push_back(data);
	}
	file.close();
	return;
}


void scan_vector2UI(std::string filename, vector2UI & info)
{
	//!- scan file into vector2F

	std::ifstream file;

	file.open(filename.c_str());
	if (!file) {
		std::cerr << "ERROR: Unable to source [ " << filename << " ]. Quit." << std::endl;
		exit(-1); // terminate with error
	}
	std::string line;

	while (getline(file, line) ) {
		vectorUI data;
		UINT value;
		std::istringstream iss(line);
		while (iss >> value) data.push_back(value);
		info.push_back(data);
	}
	file.close();
	return;
}


bool fexists(std::string filename)
{
	std::ifstream ifile(filename.c_str());

	return ifile.good();
}


std::vector<std::string> & ssplit(const std::string & s, char delim, std::vector<std::string> & elems)
{
	std::stringstream ss(s);
	std::string item;

	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}


std::vector<std::string> ssplit(const std::string & s, char delim)
{
	std::vector<std::string> elems;
	return ssplit(s, delim, elems);
}


std::string string_replace(std::string src, std::string const & target, std::string const & repl)
{
	// handle error situations/trivial cases

	if (target.length() == 0) {
		// searching for a match to the empty string will result in
		//  an infinite loop
		//  it might make sense to throw an exception for this case
		return src;
	}

	if (src.length() == 0) {
		return src;  // nothing to match against
	}

	size_t idx = 0;

	for (;; ) {
		idx = src.find(target, idx);
		if (idx == std::string::npos) break;

		src.replace(idx, target.length(), repl);
		idx += repl.length();
	}

	return src;
}



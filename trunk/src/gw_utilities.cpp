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
#include <cmath>

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


void progress_bar(unsigned int x, unsigned int N, std::string iterator_name)
{
	// how wide you want the progress meter to be
	int totaldotz = 50;
	double fraction = (1.0 * x) / (1.0 * N);
	// part of the progressmeter that's already "full"
	int dotz = (int)round(fraction * totaldotz);

	// create the "meter"
	int ii = 0;
	char buffer [10];
	int len = sprintf(buffer, "%3.0f%% [", fraction * 100);

	for (int i = 0; i < len; ++i) std::clog << buffer[i];
	// part that's full already
	for ( ; ii < dotz; ii++) {
		std::clog << "\033[0;35m" << ">" << "\033[0m";
	}
	// remaining part (spaces)
	for ( ; ii < totaldotz; ii++) {
		std::clog << "-";
	}
	// and back to line begin via "\r" - do not forget std::flush to avoid output buffering problems!
	std::clog << "] " << iterator_name << " " << x << "\r" << std::flush;
	return;
}



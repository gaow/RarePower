//!\file gw_maths.h
//!\brief useful maths functions
// Copyright 2010 Gao Wang

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

#ifndef GWMATHS_H
#define GWMATHS_H

#include<vector>

double gw_dmax(double , double ); 
int gw_imax(int , int );
double gw_dmin(double , double );
int gw_imin(int , int );
void gw_round(double&myValue, double PRECISION);

//!- Statistics

double gw_sum(const std::vector<double>&); 
double gw_prod(const std::vector<double>&); 
double gw_mean(const std::vector<double>&); 
double gw_var(const std::vector<double>&); 
double gw_sd(const std::vector<double>&);

extern double lnfact_table[]; 
double gw_ln_factorial(double );
double gw_lnchoose(double, double );
double gw_hypergeometric_pmf(unsigned int , unsigned int , unsigned int , unsigned int );
double gw_hypergeometric_cmf(unsigned int , unsigned int , unsigned int , unsigned int );

double fexact_two_sided_pvalue(const std::vector<int>&);
double Mann_Whitneyu(double*, int , double*, int );
#endif
///:~

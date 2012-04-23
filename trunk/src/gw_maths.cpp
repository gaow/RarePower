// =====================================================================================
//
//       Filename:  gw_maths.cpp
//
//    Description:  useful maths functions
//
//        Version:  1.0
//        Created:  01/04/2011 10:49:19 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Gao Wang (gw), wangow@gmail.com
//                  Baylor College of Medicine, Texas, USA
//        License:  GNU General Public License <http://www.gnu.org/licenses/>
//                  Copyright (c) 2011, Gao Wang
//
// =====================================================================================


#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>
#include <cstdlib>
#include <limits>
#include "gw_maths.h"
#include "gw_utilities.h"
#include "fisher2.h"

#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"

//!- Basic

double gw_dmax(double a, double b)
{
	if (a < b) return b;
	else return a;
}


int gw_imax(int a, int b)
{
	if (a < b) return b;
	else return a;
}


double gw_dmin(double a, double b)
{
	if (a < b) return a;
	else return b;
}


int gw_imin(int a, int b)
{
	if (a < b) return a;
	else return b;
}


void gw_round(double & myValue, double PRECISION)
{
	double myRemainder = fmod(myValue, PRECISION);

	if (myRemainder > (0.5 * PRECISION)) {
		myValue = (myValue - myRemainder + PRECISION);
	}else {
		myValue = (myValue - myRemainder);
	}
	return;
}


bool fEqual(double a, double b)
{
	return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
}


//!- Statistics

double gw_sum(const std::vector<double> & mydata)
{
	double dataSum(0.0);

	for (unsigned int i = 0; i != mydata.size(); ++i)
		dataSum = dataSum + mydata[i];
	return dataSum;
}


double gw_prod(const std::vector<double> & mydata)
{
	double dataProd(1.0);

	for (unsigned int i = 0; i != mydata.size(); ++i)
		dataProd = dataProd * mydata[i];
	return dataProd;
}


double gw_mean(const std::vector<double> & mydata)
{
	double dataSum(0.0);

	for (unsigned int i = 0; i != mydata.size(); ++i)
		dataSum = dataSum + mydata[i];
	return dataSum / double(mydata.size());
}


double gw_var(const std::vector<double> & mydata)
{
	double mMean = gw_mean(mydata);
	double dataSum(0.0);

	for (unsigned int i = 0; i != mydata.size(); ++i)
		dataSum = dataSum + pow(mydata[i] - mMean, 2.0);
	return dataSum / double( mydata.size() - 1 );
}


double gw_sd(const std::vector<double> & mydata)
{
	double mVar = gw_var(mydata);

	return sqrt(mVar);
}


//!- Factorials, combinatory and hypergeometric pmf/cmf
/*
   #include "lnfact_table"
   double gw_ln_factorial(double n)
   {

   //    double rlnFact = 0.0;
   //	if (n == 0 || n == 1) return rlnFact;
   //    for (int i = 2; i <= n; ++i) rlnFact += log(double(i));
   //	return rlnFact;


   //  return lnfact_table[int(n)];
   }


   double gw_lnchoose(double n, double k)
   {

    double lnChoose = 0.0;
    if (n == 0 || n == k || k == 0)
    return lnChoose;
    lnChoose = gw_ln_factorial(n) - gw_ln_factorial(n-k) - gw_ln_factorial(k);
    return lnChoose;
   }

   double gw_hypergeometric_pmf(unsigned int k, unsigned int n1, unsigned int n2, unsigned int t)
   {
   // prob(k) =  choose(n1, k) choose(n2, t-k) / choose(n1+n2,t)


     if (k > n1 || k > t || t - k > n2 || t > n1 + n2) {
     cout << "error!" << endl;
     exit(-1);
     }

   double dn1 = double(n1), dk = double(k), dn2 = double(n2), dt = double(t);

   double c1 = gw_lnchoose(dn1,dk);
   double c2 = gw_lnchoose(dn2,dt-dk);
   double c3 = gw_lnchoose(dn1+dn2,dt);
   double hPmf = exp(c1 + c2 - c3);

   return hPmf;
   }


   double gw_hypergeometric_cmf(unsigned int k, unsigned int n1, unsigned int n2, unsigned int t)
   {
   // prob(<= k) = \sum_0^k choose(n1, k) choose(n2, t-k) / choose(n1+n2,t)


   //        if (k > n1 || k > t || t - k > n2 || t > n1 + n2) {
   //        cout << "error!" << endl;
   //        exit(-1);
   //        }
   double hCmf = gw_hypergeometric_pmf(k, n1, n2, t);
   double hPmf = hCmf;
   for (double i = k; i >= 0; --i) {
    double diffactor = i / (n1 - i + 1.0) * (n2 - t + i) / (t - i + 1.0);
    hPmf = hPmf * diffactor;
    hCmf += hPmf;
   }

   return hCmf;
   }

 */


//!- Fisher's test for 2 by 2 tables

double fexact_two_sided_pvalue(const std::vector<int> & twotwoTable)
{
	// This is specific for 2x2 tables. contingency_table = matrix(twotwoTable, 2, 2, byrow = T)

	double contingency_table[4] = { 0, 0, 0, 0 };

	for (int i = 0; i != 4; ++i) contingency_table[i] = twotwoTable[i];
	//stuff for Fisher's test
	int nrow = 2;
	int ncol = 2;
	double expected = -1.0;
	double percnt = 100.0;
	double emin = 0.0;
	double prt = 0.0;
	double pvalue = 0.0;
	int workspace = 300000;

	fexact(&nrow, &ncol, contingency_table, &nrow, &expected, &percnt, &emin, &prt, &pvalue, &workspace);
	return pvalue;
}


//!- Mann-Whitney rank test statistic "U"
// www.alglib.net
// http://en.wikipedia.org/wiki/Mann-Whitney_U
double Mann_Whitneyu(double x[], int n, double y[], int m)
{
	int i;
	int k;
	int t;
	double tmp;
	int tmpi;
	int ns;
	double u;
	int w;

	// Prepare

	if (n <= 5 || m <= 5) {
		std::cout << "Sample size too small" << std::endl;
		exit(1);
	}
	ns = n + m;
	double r[ns - 1];
	int c[ns - 1];
	for (i = 0; i <= n - 1; i++) {
		r[i] = x[i];
		c[i] = 0;
	}
	for (i = 0; i <= m - 1; i++) {
		r[n + i] = y[i];
		c[n + i] = 1;
	}

	// sort {R, C}, QS: smaller scores ranks higher

	if (ns != 1) {
		i = 2;
		do {
			t = i;
			while (t != 1) {
				k = t / 2;
				if (r[k - 1] >= r[t - 1]) {
					t = 1;
				}else {
					tmp = r[k - 1];
					r[k - 1] = r[t - 1];
					r[t - 1] = tmp;
					tmpi = c[k - 1];
					c[k - 1] = c[t - 1];
					c[t - 1] = tmpi;
					t = k;
				}
			}
			i = i + 1;
		} while (i <= ns);
		i = ns - 1;
		do {
			tmp = r[i];
			r[i] = r[0];
			r[0] = tmp;
			tmpi = c[i];
			c[i] = c[0];
			c[0] = tmpi;
			t = 1;
			while (t != 0) {
				k = 2 * t;
				if (k > i) {
					t = 0;
				}else {
					if (k < i) {
						if (r[k] > r[k - 1]) {
							k = k + 1;
						}
					}
					if (r[t - 1] >= r[k - 1]) {
						t = 0;
					}else {
						tmp = r[k - 1];
						r[k - 1] = r[t - 1];
						r[t - 1] = tmp;
						tmpi = c[k - 1];
						c[k - 1] = c[t - 1];
						c[t - 1] = tmpi;
						t = k;
					}
				}
			}
			i = i - 1;
		} while (i >= 1);
	}

	// Compute U

	u = 0;
	w = 1;
	for (i = 0; i <= ns - 1; i++) {
		if (i == 0)
			w = 1;
		else {
			if (r[i] > r[i - 1])
				w = i + 1;
		}

		if (c[i] == 0) {
			//ranks (sum of) for cases
			//std::cout << w << " ";
			//std::cout << r[i] << " ";
			u = u + w;
		}
	}

	//u = n*m+m*(m+1)/2-u;
	//std::cout << u << " ";
	return u;
}


namespace {
const double AFFECTED = 2.0, UNAFFECTED = 1.0, HOMO_ALLELE = 2.0,
             MINOR_ALLELE = 1.0, MAJOR_ALLELE = 0.0, MISSING_ALLELE = -9.0;
}


double gwStats::testLogitRegression1(const std::vector<double> & regressors, const std::vector<double> & responses, double xbar, unsigned nCases) const
{
	assert(regressors.size() == responses.size());

	//!- Score test implementation for logistic regression model logit(p) = b0 + b1x (derivation of the test see my labnotes vol.2 page 3)

	//double ebo = (1.0 * nCases) / (1.0 * numCtrl);
	//double bo = log(ebo);

	double po = (1.0 * nCases) / (1.0 * responses.size());
	double statistic = 0.0;
	double ss = 0.0, vm1 = 0.0;
	// the score and its variance
	for (unsigned i = 0; i != regressors.size(); ++i) {
		ss += (regressors[i] - xbar) * ((responses[i] - 1.0) - po);
		vm1 += (regressors[i] - xbar) * (regressors[i] - xbar) * po * (1.0 - po);
	}


	if (!fEqual(ss, 0.0)) {
		statistic = ss / sqrt(vm1);
	}

	if (m_v) {
		std::clog << std::endl;
		std::clog << xbar << std::endl;
		std::clog << ss << std::endl;
		std::clog << vm1 << std::endl;
		std::clog << statistic << std::endl;
		exit(0);
	}

	gw_round(statistic, 1E-3);
	return statistic;
}


double gwStats::chisqtest2X2(const std::vector<double> & regressors, const std::vector<double> & responses) const
{
	assert(regressors.size() == responses.size());

	//! - 2 by 2 Chisq test
	double A0 = 0.0, A1 = 0.0, U0 = 0.0, U1 = 0.0;
	for (unsigned i = 0; i != regressors.size(); ++i) {
		if (regressors[i] == MAJOR_ALLELE && responses[i] == UNAFFECTED)
			U0 += 1.0;
		else if (regressors[i] == MAJOR_ALLELE && responses[i] == AFFECTED)
			A0 += 1.0;
		else if (regressors[i] == MINOR_ALLELE && responses[i] == UNAFFECTED)
			U1 += 1.0;
		else if (regressors[i] == MINOR_ALLELE && responses[i] == AFFECTED)
			A1 += 1.0;
		else {
			std::cerr << "Input data problem in gstat.chisqtest2X2(). Now Quit." << std::endl;
			exit(1);
		}
	}     // collect the contigency table


	double Aobs = A0 + A1;
	double Uobs = U0 + U1;
	double Tobs = Aobs + Uobs;
	double Obs0 = A0 + U0;
	double Obs1 = A1 + U1;

	double EA0 = (Aobs * Obs0) / Tobs; if (EA0 == 0) EA0 = 0.05;
	double EA1 = (Aobs * Obs1) / Tobs; if (EA1 == 0) EA1 = 0.05;
	double EU0 = (Uobs * Obs0) / Tobs; if (EU0 == 0) EU0 = 0.05;
	double EU1 = (Uobs * Obs1) / Tobs; if (EU1 == 0) EU1 = 0.05;

	double statistic = ( (A0 - EA0) * (A0 - EA0) ) / EA0
	                   + ( (A1 - EA1) * (A1 - EA1) ) / EA1
	                   + ( (U0 - EU0) * (U0 - EU0) ) / EU0
	                   + ( (U1 - EU1) * (U1 - EU1) ) / EU1;

	return statistic;
}


double gwStats::fishertest2X2(const std::vector<double> & regressors, const std::vector<double> & responses, unsigned sided, char moi) const
{
	assert(regressors.size() == responses.size());

	std::vector<int> twotwoTable(4, 0);

	for (unsigned i = 0; i != regressors.size(); ++i) {

		if (regressors[i] != MAJOR_ALLELE && regressors[i] != MINOR_ALLELE && regressors[i] != HOMO_ALLELE) {
			std::cerr << "Input data problem in gstat.fishertest2X2() (X data have missing entries). Now Quit." << std::endl;
			exit(-1);
		}

		if (responses[i] != AFFECTED && responses[i] != UNAFFECTED) {
			std::cerr << "Input data problem in gstat.fishertest2X2() table (Y data not binary). Now Quit." << std::endl;
			exit(-1);
		}

		switch (moi) {
		case 'R':
		{
			if (responses[i] == AFFECTED) {
				if (regressors[i] != HOMO_ALLELE)
					twotwoTable[0] += 1;
				else
					twotwoTable[1] += 1;
			}else {
				if (regressors[i] != HOMO_ALLELE)
					twotwoTable[2] += 1;
				else
					twotwoTable[3] += 1;
			}
		}
		break;
		case 'D':
		{
			if (responses[i] == AFFECTED) {
				if (regressors[i] == MAJOR_ALLELE)
					twotwoTable[0] += 1;
				else
					twotwoTable[1] += 1;
			}else {
				if (regressors[i] == MAJOR_ALLELE)
					twotwoTable[2] += 1;
				else
					twotwoTable[3] += 1;
			}
		}
		break;
		default:
		{
			if (responses[i] == AFFECTED) {
				if (regressors[i] == MAJOR_ALLELE)
					twotwoTable[0] += 2;
				else
					twotwoTable[1] += (unsigned)regressors[i];
			}else {
				if (regressors[i] == MAJOR_ALLELE)
					twotwoTable[2] += 2;
				else
					twotwoTable[3] += (unsigned)regressors[i];
			}
		}
		break;
		}
	}

	double pvalue2X2 = 1.0;
	if (sided == 1) {
		pvalue2X2 = (twotwoTable[3] > 0) * gsl_cdf_hypergeometric_P((twotwoTable[3] - 1), (twotwoTable[1] + twotwoTable[3]), (twotwoTable[0] + twotwoTable[2]), (twotwoTable[3] + twotwoTable[2]))
		            + 0.5 * gsl_ran_hypergeometric_pdf(twotwoTable[3], (twotwoTable[1] + twotwoTable[3]), (twotwoTable[0] + twotwoTable[2]), (twotwoTable[3] + twotwoTable[2]));
	}else {
		pvalue2X2 = fexact_two_sided_pvalue(twotwoTable);
	}

	if (m_v) {
		std::clog << twotwoTable << std::endl;
		std::clog << pvalue2X2 << std::endl;
	}
	return pvalue2X2;
}


double gwStats::testLnRegression1(const std::vector<double> & regressors, const std::vector<double> & responses, double xbar, double ybar) const
{
	assert(regressors.size() == responses.size());

	//!- Statistic: LSE (MLE) for beta, centered and scaled (bcz E[b] = 0 and sigma = 1 by simulation)
	//!- See page 41 of Kutner's Applied Linear Stat. Model, 5th ed.
	//
	double statistic = 0.0;
	double numerator = 0.0, denominator = 0.0, ysigma = 0.0;
	for (unsigned i = 0; i != regressors.size(); ++i) {
		numerator += (regressors[i] - xbar) * responses[i];
		denominator += pow(regressors[i] - xbar, 2.0);
	}

	if (!fEqual(numerator, 0.0)) {
		//!- Compute MSE and V[\hat{beta}]
		//!- V[\hat{beta}] = MSE / denominator
		double b1 = numerator / denominator;
		double b0 = ybar - b1 * xbar;
		//SSE
		for (size_t i = 0; i != regressors.size(); ++i) {
			ysigma += pow(responses[i] - (b0 + b1 * regressors[i]), 2.0);
		}

		double varb = ysigma / ((responses.size() - 2.0) * denominator);
		statistic = b1 / sqrt(varb);
	}

	return statistic;
}


double gwStats::testLnRegressionCond(const std::vector<double> & regressors, const std::vector<double> & responses, double yh, double yl) const
{
	assert(regressors.size() == responses.size());
	//!- Statistic: score test statistic for beta, given Y~N(0,1) under the null
	//!- a conditional test due to extreme sampling, Lin and Huang 2007
	//!- See my labnotes vol. 1
	//
	double fcdf = gsl_cdf_ugaussian_Q(yl) - gsl_cdf_ugaussian_Q(yh) + 1.0;
	double fpdf_1 = gsl_ran_ugaussian_pdf(yh);
	double fpdf_2 = gsl_ran_ugaussian_pdf(yl);
	double ratio = (fpdf_1 - fpdf_2) / fcdf;
	double numerator = 0.0, denominator = 0.0;
	for (unsigned i = 0; i != regressors.size(); ++i) {
		numerator += regressors[i] * (responses[i] + ratio);
		denominator += pow(regressors[i], 2.0) * (-1.0 + ratio * ratio + (fpdf_1 * yh - fpdf_2 * yl) / fcdf);
	}

	if (fEqual(denominator, 0.0)) denominator = 1.0e-10;

	double statistic = numerator / sqrt(fabs(denominator));
	return statistic;
}


double gwStats::ttestIndp(const std::vector<double> & x1s, const std::vector<double> & x2s) const
{
	assert(x1s.size() == x2s.size());

	double XA = 0.0, XU = 0.0, SA = 0.0, SU = 0.0;
	std::vector<double> VA(0), VU(0);
	for (unsigned i = 0; i != x1s.size(); ++i) {
		// genotype codings are 0 and 1 while phenotype is continuous
		if (x1s[i] == MAJOR_ALLELE)
			VU.push_back(x2s[i]);
		else if (x1s[i] == MINOR_ALLELE)
			VA.push_back(x2s[i]);
		else {
			std::cerr << "Input data problem in gstat.ttestIndp(). Now Quit." << std::endl;
			exit(1);
		}
	}

	XA = gw_mean(VA);
	XU = gw_mean(VU);
	SA = gw_var(VA);
	SU = gw_var(VU);
	if (SA == 0) SA = 1.0e-6;
	if (SU == 0) SU = 1.0e-6;

	double statistic = (XA - XU) / sqrt(SA / VA.size() + SU / VU.size());
	return statistic;
}


///:~

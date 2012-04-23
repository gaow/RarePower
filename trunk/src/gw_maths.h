// =====================================================================================
//
//       Filename:  gw_maths.h
//
//    Description:  useful maths functions
//
//        Version:  1.0
//        Created:  01/04/2011 10:49:48 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Gao Wang (gw), wangow@gmail.com
//                  Baylor College of Medicine, Texas, USA
//        License:  GNU General Public License <http://www.gnu.org/licenses/>
//                  Copyright (c) 2011, Gao Wang
//
// =====================================================================================


#ifndef GWMATHS_H
#define GWMATHS_H

#include <vector>
#include <cassert>

double gw_dmax(double, double);

int gw_imax(int, int);

double gw_dmin(double, double);

int gw_imin(int, int);

void gw_round(double & myValue, double PRECISION);

bool fEqual(double a, double b);

//!- Statistics

double gw_sum(const std::vector<double> &);

double gw_prod(const std::vector<double> &);

double gw_mean(const std::vector<double> &);

double gw_var(const std::vector<double> &);

double gw_sd(const std::vector<double> &);

//extern double lnfact_table[];
//double gw_ln_factorial(double );
//double gw_lnchoose(double, double );
//double gw_hypergeometric_pmf(unsigned int , unsigned int , unsigned int , unsigned int );
//double gw_hypergeometric_cmf(unsigned int , unsigned int , unsigned int , unsigned int );

double fexact_two_sided_pvalue(const std::vector<int> &);

double Mann_Whitneyu(double *, int, double *, int);

class gwStats
{
public:
	gwStats() :
		m_v(0)
	{
	}


	~gwStats() {}

	bool setVerbosity(unsigned level)
	{
		m_v = level;
		return true;
	}


	/*!\brief 2 by 2 chisq test statistic
	 * \param regressors
	 * \param responses
	 * \return statistic
	 */
	double chisqtest2X2(const std::vector<double> & regressors, const std::vector<double> & responses) const;

	/*!\brief 2 by 2 Fisher's test statistic
	 * \param regressors
	 * \param responses
	 * \return p-value
	 */
	double fishertest2X2(const std::vector<double> & regressors, const std::vector<double> & responses, unsigned sided, char moi = 'D') const;


	/*!\brief two-sample t statistic
	 * \param x1s
	 * \param x2s
	 * \return statistic
	 */
	double ttestIndp(const std::vector<double> & x1s, const std::vector<double> & x2s) const;


	/*!\brief Simple single variable logistic regression score statistic
	 * \param regressors
	 * \param responses
	 * \param xbar
	 * \param nCases
	 * \return statistic
	 */

	double testLogitRegression1(const std::vector<double> & regressors, const std::vector<double> & responses, double xbar, unsigned nCases) const;


	/*!\brief Simple linear regression score statisitc
	 * \param regressors
	 * \param responses
	 * \param xbar
	 * \param ybar
	 * \return statistic
	 */
	double testLnRegression1(const std::vector<double> & regressors, const std::vector<double> & responses, double xbar, double ybar) const;

	double testLnRegressionCond(const std::vector<double> & regressors, const std::vector<double> & responses, double yh, double yl) const;

private:
	bool m_v;
};
#endif
///:~

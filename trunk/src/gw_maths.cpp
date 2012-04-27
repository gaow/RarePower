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


namespace gpow {

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

}

double LogFactorial(int n)
{
	//http://www.johndcook.com/stand_alone_code.html
	if (n < 0) {
		std::cerr	<< "Invalid input argument ("
		            << n << "); may not be negative" << std::endl;
		exit(-1);

	}else if (n > 254) {
	    double pi = 3.141592653589793;
		double x = n + 1;
		return (x - 0.5) * log(x) - x + 0.5 * log(2 * pi) + 1.0 / (12.0 * x);
	}else {
		double lf[] =
		{
			0.000000000000000,
			0.000000000000000,
			0.693147180559945,
			1.791759469228055,
			3.178053830347946,
			4.787491742782046,
			6.579251212010101,
			8.525161361065415,
			10.604602902745251,
			12.801827480081469,
			15.104412573075516,
			17.502307845873887,
			19.987214495661885,
			22.552163853123421,
			25.191221182738683,
			27.899271383840894,
			30.671860106080675,
			33.505073450136891,
			36.395445208033053,
			39.339884187199495,
			42.335616460753485,
			45.380138898476908,
			48.471181351835227,
			51.606675567764377,
			54.784729398112319,
			58.003605222980518,
			61.261701761002001,
			64.557538627006323,
			67.889743137181526,
			71.257038967168000,
			74.658236348830158,
			78.092223553315307,
			81.557959456115029,
			85.054467017581516,
			88.580827542197682,
			92.136175603687079,
			95.719694542143202,
			99.330612454787428,
			102.968198614513810,
			106.631760260643450,
			110.320639714757390,
			114.034211781461690,
			117.771881399745060,
			121.533081515438640,
			125.317271149356880,
			129.123933639127240,
			132.952575035616290,
			136.802722637326350,
			140.673923648234250,
			144.565743946344900,
			148.477766951773020,
			152.409592584497350,
			156.360836303078800,
			160.331128216630930,
			164.320112263195170,
			168.327445448427650,
			172.352797139162820,
			176.395848406997370,
			180.456291417543780,
			184.533828861449510,
			188.628173423671600,
			192.739047287844900,
			196.866181672889980,
			201.009316399281570,
			205.168199482641200,
			209.342586752536820,
			213.532241494563270,
			217.736934113954250,
			221.956441819130360,
			226.190548323727570,
			230.439043565776930,
			234.701723442818260,
			238.978389561834350,
			243.268849002982730,
			247.572914096186910,
			251.890402209723190,
			256.221135550009480,
			260.564940971863220,
			264.921649798552780,
			269.291097651019810,
			273.673124285693690,
			278.067573440366120,
			282.474292687630400,
			286.893133295426990,
			291.323950094270290,
			295.766601350760600,
			300.220948647014100,
			304.686856765668720,
			309.164193580146900,
			313.652829949878990,
			318.152639620209300,
			322.663499126726210,
			327.185287703775200,
			331.717887196928470,
			336.261181979198450,
			340.815058870798960,
			345.379407062266860,
			349.954118040770250,
			354.539085519440790,
			359.134205369575340,
			363.739375555563470,
			368.354496072404690,
			372.979468885689020,
			377.614197873918670,
			382.258588773060010,
			386.912549123217560,
			391.575988217329610,
			396.248817051791490,
			400.930948278915760,
			405.622296161144900,
			410.322776526937280,
			415.032306728249580,
			419.750805599544780,
			424.478193418257090,
			429.214391866651570,
			433.959323995014870,
			438.712914186121170,
			443.475088120918940,
			448.245772745384610,
			453.024896238496130,
			457.812387981278110,
			462.608178526874890,
			467.412199571608080,
			472.224383926980520,
			477.044665492585580,
			481.872979229887900,
			486.709261136839360,
			491.553448223298010,
			496.405478487217580,
			501.265290891579240,
			506.132825342034830,
			511.008022665236070,
			515.890824587822520,
			520.781173716044240,
			525.679013515995050,
			530.584288294433580,
			535.496943180169520,
			540.416924105997740,
			545.344177791154950,
			550.278651724285620,
			555.220294146894960,
			560.169054037273100,
			565.124881094874350,
			570.087725725134190,
			575.057539024710200,
			580.034272767130800,
			585.017879388839220,
			590.008311975617860,
			595.005524249382010,
			600.009470555327430,
			605.020105849423770,
			610.037385686238740,
			615.061266207084940,
			620.091704128477430,
			625.128656730891070,
			630.172081847810200,
			635.221937855059760,
			640.278183660408100,
			645.340778693435030,
			650.409682895655240,
			655.484856710889060,
			660.566261075873510,
			665.653857411105950,
			670.747607611912710,
			675.847474039736880,
			680.953419513637530,
			686.065407301994010,
			691.183401114410800,
			696.307365093814040,
			701.437263808737160,
			706.573062245787470,
			711.714725802289990,
			716.862220279103440,
			722.015511873601330,
			727.174567172815840,
			732.339353146739310,
			737.509837141777440,
			742.685986874351220,
			747.867770424643370,
			753.055156230484160,
			758.248113081374300,
			763.446610112640200,
			768.650616799717000,
			773.860102952558460,
			779.075038710167410,
			784.295394535245690,
			789.521141208958970,
			794.752249825813460,
			799.988691788643450,
			805.230438803703120,
			810.477462875863580,
			815.729736303910160,
			820.987231675937890,
			826.249921864842800,
			831.517780023906310,
			836.790779582469900,
			842.068894241700490,
			847.352097970438420,
			852.640365001133090,
			857.933669825857460,
			863.231987192405430,
			868.535292100464630,
			873.843559797865740,
			879.156765776907600,
			884.474885770751830,
			889.797895749890240,
			895.125771918679900,
			900.458490711945270,
			905.796028791646340,
			911.138363043611210,
			916.485470574328820,
			921.837328707804890,
			927.193914982476710,
			932.555207148186240,
			937.921183163208070,
			943.291821191335660,
			948.667099599019820,
			954.046996952560450,
			959.431492015349480,
			964.820563745165940,
			970.214191291518320,
			975.612353993036210,
			981.015031374908400,
			986.422203146368590,
			991.833849198223450,
			997.249949600427840,
			1002.670484599700300,
			1008.095434617181700,
			1013.524780246136200,
			1018.958502249690200,
			1024.396581558613400,
			1029.838999269135500,
			1035.285736640801600,
			1040.736775094367400,
			1046.192096209724900,
			1051.651681723869200,
			1057.115513528895000,
			1062.583573670030100,
			1068.055844343701400,
			1073.532307895632800,
			1079.012946818975000,
			1084.497743752465600,
			1089.986681478622400,
			1095.479742921962700,
			1100.976911147256000,
			1106.478169357800900,
			1111.983500893733000,
			1117.492889230361000,
			1123.006317976526100,
			1128.523770872990800,
			1134.045231790853000,
			1139.570684729984800,
			1145.100113817496100,
			1150.633503306223700,
			1156.170837573242400,
		};
		return lf[n];
	}
}
///:~

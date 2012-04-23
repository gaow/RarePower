// =====================================================================================
//
//       Filename:  assoctests.h
//
//    Description:  Statistical tests for association mapping
//
//        Version:  1.0
//        Created:  01/04/2011 10:46:12 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Gao Wang (gw), wangow@gmail.com
//                  Baylor College of Medicine, Texas, USA
//        License:  GNU General Public License <http://www.gnu.org/licenses/>
//                  Copyright (c) 2011, Gao Wang
//
// =====================================================================================


#ifndef ASSOC_TESTS_H
#define ASSOC_TESTS_H

#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <map>
#include <algorithm>
#include <ctime>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"
#include "gw_utilities.h"
#include "gw_maths.h"

namespace gpow {

inline int randbin () { return (rand()%2); }

class gwAssocdata
{
public:
	/*!\brief Constructor
	 * \param observedMafs Observed minor allele frequencies in the data-set
	 * \param ydat
	 * \param xdat
	 * \param mafLower Lower bound of MAF that shall be analyzed. Many rare variants methods have it as (0, 0.01]
	 * \param mafUpper Upper bound
	 * \param alpha pre-set (desired) significance level
	 */
	gwAssocdata() :
		__xdat(0), __ydat(0), __observedMafs(0),
		__mafLower(0.0), __mafUpper(0.0), __permutator(0),
        __v(0), __randbins(0)
	{
	}


	virtual ~gwAssocdata() {}

	virtual gwAssocdata * clone() const
	{
		return new gwAssocdata(*this);
	}


	bool setVerbosity(unsigned level)
	{
		__v = level;
		return true;
	}


	void locals() const
	{
		std::clog.precision(9);
		std::clog << "__observedMafs:\n" << __observedMafs << "\n" << std::endl;
		std::clog << "__ydat:\n" << __ydat << "\n" << std::endl;
		std::clog << "__xdat:\n" << __xdat << "\n" << std::endl;
		exit(0);
	}


	bool setXdat(const vector2F xdat)
	{
		__xdat = xdat;
		if (__ydat.size() != 0) assert(__xdat.size() == __ydat.size());
		if (__v) std::cout << __xdat << std::endl;
		return true;
	}


	bool setYdat(const vectorF ydat)
	{
		__ydat = ydat;
		if (__xdat.size() != 0) assert(__ydat.size() == __xdat.size());
		if (__v) std::cout << __ydat << std::endl;
		return true;
	}


	bool setMafs(const vectorF mafs, double mafLower, double mafUpper)
	{
		__observedMafs = mafs;
		if (__xdat.size() != 0) assert(__observedMafs.size() == __xdat.front().size());
		assert(mafLower >= 0.0 && mafUpper <= 1.0 && mafUpper > mafLower);
		__mafLower = mafLower;
		__mafUpper = mafUpper;
		if (__v) std::cout << __observedMafs << std::endl;
		return true;
	}

    bool setPermutator(unsigned permutator) {
        __permutator = permutator;
        if (__permutator == 1) {
            assert (__ydat.size() % 2 == 0);
            __randbins.resize(__ydat.size() / 2);
            std::generate(__randbins.begin(), __randbins.end(), randbin);
            
        }
        return true;
    }


    bool permutate() {
        switch (__permutator) {
            case 1:
                // sibpare permutation
                // permute the case control status within the pair, not across studies
                {
                    for (size_t i = 0; i < (__ydat.size() - 1); i += 2) {
                        if (__randbins[i/2]) std::iter_swap(__ydat.begin()+i, __ydat.begin()+i+1);
                    }
                    // modify the random binary sequence ... this avoids using gsl rng 
                    random_shuffle(__randbins.begin(), __randbins.end());
                }
                break;
            default:
                random_shuffle(__ydat.begin(), __ydat.end());
                break;
        }
        return true;
    }

	vectorF & ydat() { return __ydat; }
	vector2F & xdat() { return __xdat; }
	vectorF & mafs() { return __observedMafs; }

	/*!\brief
	 * \param
	 * \return
	 */
	bool trimXdat();

	bool markwildSibpairloci();


	/*!\brief ANRV scoring theme
	 * \return
	 */
	vectorF countRegionalVariants() const;


	/*!\brief CMC scoring theme
	 * \return
	 */
	vectorF binariesRegionalVariants() const;


	/*!\brief RVE scoring theme
	 * \return
	 */
	vectorF binariesRegionalUniqueVariants() const;

	/*\brief Adaptive p-value calculation
	   *\return a logical variable indicating whether or not should continue permutations
	 */

private:
	//!\brief 2D object, Sample genotypes
	vector2F __xdat;
	//!\brief 1D object, phenotypes
	vectorF __ydat;
	//!\brief Observed minor allele frequency in the data-set
	vectorF __observedMafs;
	//!\brief Upper and lower bounds of MAF that shall be analyzed. Many rare variants methods have it as (0, 0.01]
	double __mafLower;
	//!\brief Upper and lower bounds of MAF that shall be analyzed. Many rare variants methods have it as (0, 0.01]
	double __mafUpper;
    //!\brief define data permutation theme
    unsigned __permutator;
	//!\brief Verbosity flag
	unsigned __v;
    //!\brief a random vector of binary 0 or 1
    vectorI __randbins;
};


/*!\brief Rare variants association tests. Notes:
 * <br>
 * The testRare() method was adopted from the author's original code.
 *  Credit goes to Dr. Iuliana Ionita-Laza at Columbia U.
 */

class gwBaseTest
{
	typedef std::map<std::string, double> DoubleVars;

public:
	gwBaseTest() : __alpha(0.05), __v(0),
		__sided(0), __moi('A'), m_gstat()
	{
	}


	virtual ~gwBaseTest()
	{
	}


	virtual gwBaseTest * clone() const
	{
		return new gwBaseTest(*this);
	}


	virtual double apply(gwAssocdata & d)
	{
		std::cerr << "gwBaseTest should not be called" << std::endl;
		return 1.0;
	}


	// return a string of the class name
	virtual std::string name()
	{
		return "GWBASETEST";
	}


	bool setVerbosity(const unsigned v)
	{
		__v = v;
		return true;
	}


	bool setA(const double alpha)
	{
		assert(alpha > 0.0 && alpha < 1.0);
		__alpha = alpha;
		return true;
	}


	virtual bool setMOI(const char moi)
	{
		assert(moi == 'A' || moi == 'D' || moi == 'M' || moi == 'C' || moi == 'R');
		__moi = moi;
		return true;
	}


	bool setAlternative(unsigned sided)
	{
		assert(sided == 1 || sided == 2);
		__sided = sided;
		return true;
	}


	bool setPermutations(unsigned nperm, unsigned permcheckpnt)
	{
		assert(nperm > 0);
		__nperm = nperm;
		__permcheckpnt = permcheckpnt;
		return true;
	}


	virtual bool useMidP()
	{
		return true;
	}


	// store/get a double value with name 'name'
	void setVar(const std::string & name, const double value)
	{
		m_doubleVars[name] = value;
	}


	double getDoubleVar(const std::string & name)
	{
		DoubleVars::iterator it = m_doubleVars.find(name);

		if (it == m_doubleVars.end()) {
			std::cerr << "No double variable with name " + name + " can be found" << std::endl;
			exit(1);
		}
		return it->second;
	}


	double checkP(unsigned pcount1, unsigned pcount2, size_t current) const;

protected:
	//!\brief pre-set (desired) significance level
	double __alpha;
	//!\brief Verbosity flag
	unsigned __v;
	//!\brief alternative hypothesis
	unsigned __sided;
	//!\brief number of permutations / check-point
	unsigned __nperm;
	unsigned __permcheckpnt;
	//!\brief moi
	char __moi;
	gwStats m_gstat;
	// arbitrary double type of variables
	DoubleVars m_doubleVars;
};


typedef std::vector<gwBaseTest *> vectora;


//!\brief Collapsed variants simple LR model score test statistic
class CmcstP : public gwBaseTest
{
public:
	CmcstP() : gwBaseTest()
	{
	}


	gwBaseTest * clone() const
	{
		return new CmcstP(*this);
	}


	double apply(gwAssocdata & d);

};


//!\brief Counts of variants simple LR model score test statistic
class AnrvstP : public gwBaseTest
{
public:
	AnrvstP() : gwBaseTest()
	{
	}


	gwBaseTest * clone() const
	{
		return new AnrvstP(*this);
	}


	double apply(gwAssocdata & d);

};


//!\brief Collapsed variants 2 by 2 chisq statistic
class CmcchiP : public gwBaseTest
{
public:
	CmcchiP() : gwBaseTest()
	{
	}


	gwBaseTest * clone() const
	{
		return new CmcchiP(*this);
	}


	double apply(gwAssocdata & d);

};

//!\brief CMC Fisher's exact test p-value no permutation
class CmcfisherP : public gwBaseTest
{
public:
	CmcfisherP() : gwBaseTest()
	{
	}


	gwBaseTest * clone() const
	{
		return new CmcfisherP(*this);
	}


	double apply(gwAssocdata & d);

};

//!\brief Cohen's RVE fisher exact test p-value no permutation
class RvefisherP : public gwBaseTest
{
public:
	RvefisherP() : gwBaseTest()
	{
	}


	gwBaseTest * clone() const
	{
		return new RvefisherP(*this);
	}


	double apply(gwAssocdata & d);

};

//!\brief The WSS method rank statistic, page 3 of Browning paper.
//!- Not using normal approximation
class WssRankP : public gwBaseTest
{
public:
	WssRankP() : gwBaseTest()
	{
	}


	gwBaseTest * clone() const
	{
		return new WssRankP(*this);
	}


	double apply(gwAssocdata & d);

};


//!- Using normal approximation
class WssRankPA : public gwBaseTest
{
public:
	WssRankPA() : gwBaseTest()
	{
	}


	gwBaseTest * clone() const
	{
		return new WssRankPA(*this);
	}


	double apply(gwAssocdata & d);

};


//!\brief The KBAC method, Dajiang Liu's paper
class KbacP : public gwBaseTest
{
public:
	KbacP() : gwBaseTest()
	{
	}


	gwBaseTest * clone() const
	{
		return new KbacP(*this);
	}


	double apply(gwAssocdata & d);

};


class KbacstP : public gwBaseTest
{
public:
	KbacstP() : gwBaseTest()
	{
	}


	gwBaseTest * clone() const
	{
		return new KbacstP(*this);
	}


	double apply(gwAssocdata & d);

};

//!\brief The variable threshold method, Price's paper 2010 AJHG

class VtP : public gwBaseTest
{
public:
	VtP() : gwBaseTest()
	{
	}


	gwBaseTest * clone() const
	{
		return new VtP(*this);
	}


	double apply(gwAssocdata & d);

};


//!\brief The variable threshold method with Fisher's test
class VtFisherP : public gwBaseTest
{
public:
	VtFisherP() : gwBaseTest()
	{
	}


	gwBaseTest * clone() const
	{
		return new VtFisherP(*this);
	}


	double apply(gwAssocdata & d);

};


//!\brief The adaptive sum test (in Pan and Han 2010) + permutation
class AsumP : public gwBaseTest
{
public:
	AsumP() : gwBaseTest()
	{
	}


	gwBaseTest * clone() const
	{
		return new AsumP(*this);
	}


	double apply(gwAssocdata & d);

};


//!\brief CMCQT
class CmcqtP : public gwBaseTest
{
public:
	CmcqtP() : gwBaseTest()
	{
	}


	gwBaseTest * clone() const
	{
		return new CmcqtP(*this);
	}


	double apply(gwAssocdata & d);

};

//!\brief ANRVQT
class AnrvqtPermP : public gwBaseTest
{
public:
	AnrvqtPermP() : gwBaseTest()
	{
	}


	gwBaseTest * clone() const
	{
		return new AnrvqtPermP(*this);
	}


	double apply(gwAssocdata & d);

};


//!\brief ANRVQT
class AnrvqtP : public gwBaseTest
{
public:
	AnrvqtP() : gwBaseTest()
	{
	}


	gwBaseTest * clone() const
	{
		return new AnrvqtP(*this);
	}


	double apply(gwAssocdata & d);

};


//!\brief ANRVQT for Xtreme sampling
class AnrvqtCondP : public gwBaseTest
{
public:
	AnrvqtCondP() : gwBaseTest()
	{
	}


	gwBaseTest * clone() const
	{
		return new AnrvqtCondP(*this);
	}


	double apply(gwAssocdata & d);

};


//!\brief Ionita-Laza and Lange 2011 PLoS Genet.
class TestRareP : public gwBaseTest
{
public:
	TestRareP() : gwBaseTest()
	{
	}


	gwBaseTest * clone() const
	{
		return new TestRareP(*this);
	}


	double apply(gwAssocdata & d);

};


//!\brief c-alpha test, Ben. Neale et al. 2011 PLoS Genet.
class CalphaP : public gwBaseTest
{
public:
	CalphaP() : gwBaseTest()
	{
	}


	gwBaseTest * clone() const
	{
		return new CalphaP(*this);
	}


	double apply(gwAssocdata & d);

};


//!\brief RareCover method, 2010 PLoS CompBio
class RareCoverP : public gwBaseTest
{
public:
	RareCoverP() : gwBaseTest()
	{
	}


	gwBaseTest * clone() const
	{
		return new RareCoverP(*this);
	}


	double apply(gwAssocdata & d);

};


//!\brief Weighted fisher's test, 2011 Shuang Wang
class WsFisherP : public gwBaseTest
{
public:
	WsFisherP() : gwBaseTest(), __isMidP(0)
	{
	}


	gwBaseTest * clone() const
	{
		return new WsFisherP(*this);
	}


	double apply(gwAssocdata & d);

	bool useMidP()
	{
		__isMidP = true;
		return true;
	}


private:
	bool __isMidP;

};

//!\brief skat test, 2011 xlin
class SkatP : public gwBaseTest
{
public:
	SkatP() : gwBaseTest()
	{
	}


	gwBaseTest * clone() const
	{
		return new SkatP(*this);
	}


	double apply(gwAssocdata & d);

};
}
#endif ///:~

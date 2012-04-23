// =====================================================================================
//
//       Filename:  main.h
//
//    Description:  simulations, sampling, and power calculations for gene based association studies.
//
//        Version:  1.0
//        Created:  01/04/2011 10:52:55 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Gao Wang (gw), wangow@gmail.com
//                  Baylor College of Medicine, Texas, USA
//        License:  GNU General Public License <http://www.gnu.org/licenses/>
//                  Copyright (c) 2011, Gao Wang
//
// =====================================================================================

#include "assoctests.h"

std::string check_options(std::string prog_name, std::string & projectName, std::string & gFile, double & boundary, double & neutral_cutoff,
	vectorF & propFunctionalRv, char & moi, std::string & simulationTask, vectorF & oddsRatios, double & baselinef, vectorF & pars,
	bool & isParVariable, vectorF & qtcoefs, vectorF & qtcuts, bool & shouldMarkBin, double & percentageCausal,
	bool & isMendelAlleleFixed, vectorF & propMissingData, double & missingLowMaf, bool & shouldMarkMissing, unsigned & nCases,
	unsigned & nCtrls, unsigned & nPopulation, unsigned & nUnphenotyped, double & propHeterCases, bool & isSynoTrimmed,
	bool & isCvTrimmed, bool & isPedWritten, std::string & test, double & mafLower,
	double & mafUpper, double & alpha, unsigned & nPermutations, unsigned & nReplicates, bool & verbose, bool & quiet, unsigned & seed,
	bool & shouldUseGenPool);

gpow::gwBaseTest * testFactory(std::string const & classname)
{
	if (classname == "CMC" || classname == "CMC-one") {
		return new gpow::CmcfisherP();
	} else if (classname == "WSS" || classname == "WSS-one") {
		return new gpow::WssRankPA();
	} else if (classname == "RVE" || classname == "RVE-one") {
		return new gpow::RvefisherP();
	} else if (classname == "CMCST" || classname == "CMCST-one") {
		return new gpow::CmcstP();
	} else if (classname == "MZ" || classname == "MZ-one") {
		return new gpow::AnrvstP();
	} else if (classname == "CMCPM") {
		return new gpow::CmcchiP();
	} else if (classname == "WSSPM" || classname == "WSSPM-one") {
		return new gpow::WssRankP();
	} else if (classname == "KBAC-one" || classname == "KBAC") {
		return new gpow::KbacP();
	} else if (classname == "KBACST-one" || classname == "KBACST") {
		return new gpow::KbacstP();
	} else if (classname == "VT" || classname == "VT-one") {
		return new gpow::VtP();
	} else if (classname == "VTfisher" || classname == "VTfisher-one") {
		return new gpow::VtFisherP();
	} else if (classname == "aSum") {
		return new gpow::AsumP();
	} else if (classname == "CMCQT" || classname == "CMCQT-one") {
		return new gpow::CmcqtP();
	} else if (classname == "MZQT" || classname == "MZQT-one") {
		return new gpow::AnrvqtP();
	} else if (classname == "MZQTPM" || classname == "MZQTPM-one") {
		return new gpow::AnrvqtPermP();
	} else if (classname == "ExtremeQT" || classname == "ExtremeQT-one") {
		return new gpow::AnrvqtCondP();
	} else if (classname == "RBT" || classname == "RBT-one") {
		return new gpow::TestRareP();
	} else if (classname == "calpha" || classname == "calpha-one") {
		return new gpow::CalphaP();
	} else if (classname == "RareCover" || classname == "RareCover-one") {
		return new gpow::RareCoverP();
	} else if (classname.find("WF") < classname.size()) {
		return new gpow::WsFisherP();
	} else if (classname == "SKAT") {
		return new gpow::SkatP();
	} else {
		return NULL;
	}
}


const char* args_dsc(std::string name, bool empty)
{
	if (empty) return "";
	if (name == "task")
	{
		return "\n\tAnalysis task. \n\t| Type values \"1~7\"\n\t| 1) Case-ctrl samples given odds ratio and prevalence \n\t| 2) Population samples given odds ratio and prevalence \n\t| 3) Case-ctrl samples given population attributable risk \n\t| 4) Quantitative traits samples \n\t| 5) Extreme quantitative traits samples \n\t| 6) Mendelian traits samples \n\t| 7) Affected/unaffected sib-pairs \n\t";
	}
	if (name == "gdata")
	{
		return "\n\tGenetic data files for the simulation to be based on. \n\t| STRING\n\t| Proper gdata.maf, gdata.ann and gdata.pos files need to be provided to the program (gdata.hap file will be needed if --use_haplotype_pool is envoked). \n\t";
	}
	if (name == "pname")
	{
		return "\n\tProject name. \n\t| STRING \n\t| set output file names. \n\t| with \"-v\" option, the program will generate summary files of power calculation \n\t| with \"-z\" option, the program will generate simulated data only\n\t";
	}
	if (name == "f")
	{
		return "\n\tDefinition of rare variants. \n\t| FLOAT\n\t| variants having MAF <= frequency will be considered as \"rare\" variants\n\t";
	}
	if (name == "q")
	{
		return "\n\tProportion of FUNCTIONAL deleterious variants. \n\t| FLOAT\n\t| 0.0 <= fraction <= 1.0\n\t| (1-proportion)x100\% is the proportion of non-causal deleterious variants (noise) \n\t| This is NOT the proportion of deleterious variants, which should have been defined in <gdata.ann> file\n\t";
	}
	if (name == "p")
	{
		return "\n\tProportion of FUNCTIONAL protective variants. \n\t| FLOAT\n\t| 0.0 <= fraction <= 1.0\n\t| (1-proportion)x100\% is the proportion of non-causal protective variants (noise) \n\t| This is NOT the proportion of protective variants, which should have been defined in <gdata.ann> file\n\t";
	}
	if (name == "g")
	{
		return "\n\tMode of inheritance under which the phenotype data is simulated.\n\t| CHAR\n\t| \"A\" (additive), \"D\" (dominant), \"R\" (recessive), \"M\" (multiplicative), \n\t| \"C\" (compound dominant for non-mendelian traits, or compound recessive for mendelian traits)\n\t";
	}
	if (name == "A")
	{
		return "\n\t[task=1,2] Minimum odds ratio for deleterious variants. \n\t| FLOAT\n\t| 0.0 for fixed effect size model; >=1.0 for variable effect sizes model\n\t";
	}
	if (name == "B")
	{
		return "\n\t[task=1,2] Maximum odds ratio for deleterious variants. \n\t| FLOAT\n\t| >=1.0 AND > $OR_deleterious_min\n\t| will be the odds ratio for deleterious variants in fixed effect size model, or the maximum odds ratio in variable effect sizes model\n\t";
	}
	if (name == "C")
	{
		return "\n\t[task=1,2] Minimum odds ratio for protective variants. \n\t| FLOAT\n\t| 0.0 for fixed effect size model; <1.0 for variable effect sizes model\n\t";
	}
	if (name == "D")
	{
		return "\n\t[task=1,2] Maximum odds ratio for protective variants. \n\t| FLOAT\n\t| <=1.0 AND > $OR_protective_min\n\t| will be the odds ratio for protective variants in fixed effect size model, or the maximum odds ratio in variable effect sizes model\n\t";
	}
	if (name == "E")
	{
		return "\n\t[task=1,2] Odds ratio for common variants. \n\t| FLOAT\n\t| =1.0 for neutral, >1.0 for deleterious, <1.0 for protective \n\t| this is the odds ratio for all variants having MAF > $define_rare, i.e., the \"common\" variants \n\t| assuming all common variants have fixed effect size\n\t";
	}
	if (name == "F")
	{
		return "\n\t[task=1,2] Disease prevalence. \n\t| FLOAT \n\t| will be used as baseline penetrance of a gene (baseline penetrance ~= disease prevalence)\n\t";
	}
	if (name == "G")
	{
		return "\n\t[task=3] Total population attributable risk for deleterious variants. \n\t|FLOAT \n\t| 0.0 <= fraction <= 1.0\n\t";
	}
	if (name == "H")
	{
		return "\n\t[task=3] Total population attributable risk for protective variants. \n\t|FLOAT \n\t| 0.0 <= fraction <= 1.0\n\t";
	}
	if (name == "I")
	{
		return "\n\t[task=3] Locus specific population attributable risk is inversely proportional to its MAF (rather than uniformly distributed).\n\t| BOOLEAN \n\t";
	}
	if (name == "J")
	{
		return "\n\t[task=4,5] Minimum mean value shift per variant. \n\t| FLOAT\n\t| 0.0 for fixed effect size model; >0.0 for variable effect sizes model\n\t";
	}
	if (name == "K")
	{
		return "\n\t[task=4,5] Maximum mean value shift per variant. \n\t| FLOAT\n\t| >=0.0 AND > $QT_effect_min \n\t| will be locus effect to quantitative trait for fixed effect size model, or maximum effect for variable effect sizes model\n\t| the mean of the quantitative trait will be shifted by (multiplier*sigma) where sigma is standard deviation of the quantitative trait\n\t| will automatically shift trait values up for deleterious variants, down for protective variants\n\t";
	}
	if (name == "L")
	{
		return "\n\t[task=4,5] Mean value shift for common variants. \n\t| FLOAT\n\t| >=0.0 \n\t| this is the effect for all variants having MAF > $define_rare, i.e., the \"common\" variants \n\t| assuming all common variants have fixed effect size\n\t";
	}
	if (name == "M")
	{
		return "\n\t[task=5] Lower percentile cutoff for quantitative traits in extreme QT sampling. \n\t| FLOAT\n\t| 0.0 <= fraction <= $QT_upper_percentile\n\t";
	}
	if (name == "N")
	{
		return "\n\t[task=5] Upper percentile cutoff for quantitative traits in extreme QT sampling. \n\t| FLOAT\n\t| $QT_lower_percentile <= fraction <= 1.0\n\t";
	}
	if (name == "O")
	{
		return "\n\t[task=5] Re-code extreme quantitative traits using binary codings.\n\t| BOOLEAN\n\t| if envoked, will convert extreme quantitative traits into binary traits\n\t";
	}
	if (name == "P")
	{
		return "\n\t[task=6] Percentage of rare variants being causal in Mendelian traits simulation. \n\t| FLOAT\n\t| 0.0 <= fraction <= 1.0\n\t| any of the variants having the top (fraction)x100\% smallest MAF will be contributable to a Mendelian trait \n\t| this assumes the disease is allelic heterogeneous \n\t";
	}
	if (name == "Q")
	{
		return "\n\t[task=6] Porportion of cases that do not carry the disease allele at the gene region under study. \n\t| FLOAT\n\t| 0.0 <= fraction <= 1.0\n\t";
	}
	if (name == "R")
	{
		return "\n\t[task=6] No allelic heterogeneity for Mendelian traits. \n\t| BOOLEAN\n\t| if envoked, will fix the causal variants of the Mendelian trait to the one that has the (fraction)x100\%-th smallest MAF \n\t| rather than using all variants having the top (fraction)x100\% smallest MAF\n\t";
	}
	if (name == "W")
	{
		return "\n\t[task=2,4,5] Number of total samples. \n\t| INT (>=0)\n\t| for extreme QT simulations (task=5) $num_all_samples>0 will have extreme QT samples collected from this finite cohort.\n\t";
	}
	if (name == "X")
	{
		return "\n\t[task=1,3,5,6] Number of cases. \n\t| INT (>0)\n\t| for extreme QT simulations (task=5) it will be #samples having high QT values from the population when $num_all_samples is set to 0\n\t";
	}
	if (name == "Y")
	{
		return "\n\t[task=1,3,5,6] Number of ctrls. \n\t| INT (>0)\n\t| for extreme QT simulations (task=5) it will be #samples having low QT values from the population when $num_all_samples is set to 0\n\t";
	}
	if (name == "Z")
	{
		return "\n\t[task=1,3,6] Number of unphenotyped cohort ctrls. \n\t| INT (>0)\n\t";
	}
	if (name == "U")
	{
		return "\n\t[task=1,2,4,5] Randomly sample haplotypes from haplotype pool file $gdata.hap, rather than generating haplotypes on the fly. \n\t| BOOLEAN \n\t";
	}
	if (name == "a")
	{
		return "\n\tProportion of deleterious variants missing data. \n\t| FLOAT \n\t| missing genotypes will be coded as wildtype genotype by default\n\t";
	}
	if (name == "b")
	{
		return "\n\tProportion of protective variants missing data. \n\t| FLOAT \n\t| missing genotypes will be coded as wildtype genotype by default\n\t";
	}
	if (name == "c")
	{
		return "\n\tProportion of non-causal variants missing data. \n\t| FLOAT \n\t| missing genotypes will be coded as wildtype genotype by default\n\t";
	}
	if (name == "d")
	{
		return "\n\tProportion of synonymous variants missing data. \n\t| FLOAT \n\t| missing genotypes will be coded as wildtype genotype by default\n\t";
	}
	if (name == "e")
	{
		return "\n\tVariants having MAF < $missing_low_maf will be marked as missing data. \n\t| FLOAT \n\t| note that $missing_low_maf is compared against the haplotype pool, not the sample\n\t";
	}
	if (name == "k")
	{
		return "\n\tRe-code missing data. \n\t| BOOLEAN\n\t| if envoked, will re-code missing data from wildtype genotype to \"-9\", indicating missingness\n\t";
	}
	if (name == "i")
	{
		return "\n\tKeep synonymous variants from analysis. \n\t| BOOLEAN\n\t";
	}
	if (name == "j")
	{
		return "\n\tRemove common variant sites from analysis. \n\t| BOOLEAN\n\t| the \"common\" loci refers to variants in the haplotype pool having MAF > $define_rare \n\t";
	}
	if (name == "l")
	{
		return "\n\tLower bound of observed sample minor allele frequency. \n\t| FLOAT\n\t| loci having observed MAF < $maf_lower will not be analyzed \n\t";
	}
	if (name == "m")
	{
		return "\n\tUpper bound of observed sample minor allele frequency. \n\t| FLOAT\n\t| loci having observed MAF > $maf_upper will not be analyzed \n\t";
	}
	if (name == "n")
	{
		return "\n\tAnnotation value cut-off that defines a variant to be \"neutral\" in evolution (either synonymous or non-coding). \n\t| FLOAT\n\t| loci having annotation value (in gdata.ann) between (-$annotation_cutoff, +$annotation_cutoff) will be regarded \"neutral\" and will not contribute to phenotype \n\t";
	}
	if (name == "t")
	{
		return "\n\tAssociation test method. \n\t| STRING \n\t| CMC, CMC-one, WSS, WSS-one, RVE, RVE-one, CMCST, CMCST-one, WSSPM, WSSPM-one\n\t| MZ, MZ-one, KBAC, KBAC-one, KBACST, KBACST-one, VT, VT-one, VTfisher, VTfisher-one \n\t| aSum, RBT, RBT-one, calpha, RareCover, RareCover-one, WF, WF-one, SKAT\n\t| CMCQT, CMCQT-one, MZQT, MZQT-one, MZQTPM, MZQTPM-one\n\t";
	}
	if (name == "s")
	{
		return "\n\tSignificance level at which power will be evaluated. \n\t| FLOAT\n\t";
	}
	if (name == "r")
	{
		return "\n\tNumber of replicates for power evaluation. \n\t| INT (>0)\n\t";
	}
	if (name == "u")
	{
		return "\n\tNumber of permutations, only applicable to permutation based methods. \n\t| INT (>0)\n\t";
	}
	if (name == "y")
	{
		return "\n\tSeed for random number generator. \n\t| INT (>=0) \n\t| =0 is to use a random seed (seed = system time + process ID)\n\t";
	}
	if (name == "z")
	{
		return "\n\tWrite out simulated genotype-phenotype data to file $pname, rather than calculating power. \n\t| BOOLEAN \n\t";
	}
	if (name == "v")
	{
		return "\n\t Maximal screen output information as well as file output for intermediate statistic such as p-values, etc.\n\t";
	}
	if (name == "x")
	{
		return "\n\t Only output result to screen. No file output.\n\t";
	}
    std::cerr << "Error input" << std::endl;
    return "";
}

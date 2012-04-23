# /home/gw/Toolkit/sfscode/bin/sfs_code 1 20 -L 1 1500 --popSize 7947 -Td 0 0.032968 -Td 0.005285 26.79 -Td 0.291997 7.53768 -TE 0.328237 --selDistType 2 0 1 1 0.206 15400 --sampSize 10000 --outfile out.txt --errfile err.txt --popFreq freq.txt
# http://www.plosgenetics.org/article/info:doi%2F10.1371%2Fjournal.pgen.1000083
import os
import sys
fname=sys.argv[1]

outtxt = [x for x in file(fname+".out")]
freqtxt = [x.split() for x in file(fname+".freq")]

# freq.txt
#Each line will contain:
#the locus carrying the mutation, the site of the mutation within the locus, the generation the
#mutation arose, and the type of mutation (`s'ubstitution, `i'nsertion, `d'eletion, or in`v'ersion)
#followed by the ancestral and derived alleles (respectively). This will be followed by a tab delimited
#list containing the frequency of the mutation in the sample from population 0, the frequency in all
#of population 0, the frequency in the sample from population 1, the frequency in all of population
#1, etc. This information can be combined with the information in the general output le from
#SFS CODE to access most necessary information.
mafs = {}
for item in freqtxt:
    if(float(item[-2])):
        try:
            mafs[item[0]].append(float(item[-1]))
        except:
            mafs[item[0]] = []
            mafs[item[0]].append(float(item[-1]))
    else:
        pass

# out.txt
#The information provided for each mutation are as follows:
#1. locus that the mutation arose on (zero-based)
#2. `A', 'X', 'Y', indicating Autosomal, X-, or Y-linked mutation, respectively
#3. position of mutation in locus (zero based)
#4. generation mutation arose (negative for mutations arising during burn-in)
#5. generation mutation fixed in population (or time of sampling if segregating)
#6. ancestral trinucleotide (middle nucleotide mutated, NOT CODON)
#7. derived nucleotide
#8. 0 or 1 for synonymous or nonsynonymous (respectively; 0 for non-coding)
#9. ancestral amino acid (single character representation; `X' for non-coding)
#10. derived amino acid (single character; `X' for non-coding)
#11. fitness effect (this is s, NOT gamma = PNs)
#12. number of chromosomes (n) that carry the mutation
#13-. . . comma delimited list of the n chromosomes carrying mutation

#info=[outtxt[idx+1] for idx, item in enumerate(outtxt) if item.startswith("MALES") and idx <= len(outtxt) -1]
info = []
tmp = ""
flag = False
for idx, item in enumerate(outtxt):
    if item.startswith("MALES"):
        flag = True
        tmp = ""
        continue
    if not flag:
        continue
    if item.startswith("//"):
        flag = False
        if not len(tmp)==0:
            info.append(tmp)
    elif idx == len(outtxt) - 1:
            tmp += item.split("\n")[0]
            info.append(tmp)
    else:
        tmp += item.split("\n")[0]

maffile = open(fname+"_sfscode"+".maf", 'w')
selfile = open(fname+"_sfscode"+".ann", 'w')
posfile = open(fname+"_sfscode"+".pos", 'w')
synfile = open(fname+"_sfscode"+".syn", 'w')

for idx, iter in enumerate(info):
    lociinfo = [x.split(',') for x in iter.split(";") if x !='\n']
    del lociinfo[len(lociinfo)-1]
    positions = [x[2] for x in lociinfo] #
    generation = [x[3] for x in lociinfo]
    issynonymous = [x[7] for x in lociinfo]
    ancestral_aminoacid = [x[8] for x in lociinfo]
    derived_aminoacid = [x[9] for x in lociinfo]
    selection_coef = [-1.0*(float(x[10])) if float(x[10]) != 0.0 else float(x[10]) for x in lociinfo] #
    maf = mafs[str(idx)] #
    #maf = [maf_dict[x][1] for x in positions if (x in maf_dict && )
    
    if len(maf) != len(positions):
        print "Something wrong with iteration " + str(idx)
        #print len(maf)
        #print maf
        #print len(positions)
        #print positions
    else:
        for item in maf:
            maffile.write(str(item) + "\t")
        maffile.write("\n")
        for item in selection_coef:
            selfile.write(str(item) + "\t")
        selfile.write("\n")
        for item in positions:
            posfile.write(str(item) + "\t")
        posfile.write("\n")
        for item in issynonymous:
            synfile.write(str(item) + "\t")
        synfile.write("\n")
    
maffile.close()
selfile.close()
posfile.close()
synfile.close()

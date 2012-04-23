# command from Deng el al 2011
/home/$USER/Toolkit/sfscode/bin/sfs_code 1 1 -t 0.00031788 -r 0.00031788 -B 10 -L 1 $k --popSize 7947 -Td 0 0.032968 -Td 0.005285 26.79 -Td 0.291997 7.53768 -TE 0.328237 --sampSize 10000 --outfile out.txt --selDistType 2 0 1 1 1.02 800 --popFreq freq.txt
# where k was replaced by either gene length in unit of base pair.

###
# parameters for the model
###

-L 1 1500 # one loci, length is 1500bp
-Td 0 0.032968 # 0, 262/7947
-Td 0.005285 26.79 # 84/7947/2, 7019/262
-Td 0.328237 7.537683 # 5217/7947/2, 52907/7019
-TE 0.03624009 # 576/7947/2 ???
--selDistType 2 0 1 1 0.206 15400 # 0.146*2*52907, it is the dist for gamma, not for s (s has to be scaled)

0.146*2*7947 = 2320.5
0.146*2*262 = 76.504
0.146*2*7019 = 2049.5
0.146*2*52907 = 15400


--sampSize -1 # sample all population
--BURN # default is 5*2*popSize
--theta # 4N*mu, mutation rate 
--rho # 4N*r, recombination rate
--popSize 7947 # not seeing it as necessary if you are simulating a distribution of selective effects where the mean of the distribution is greater than the population size (in absolute value), then the entire population might go extinct.


# my proposed command for Boyko 2009 CEU 
/home/$USER/Toolkit/sfscode/bin/sfs_code 1 20 -L 1 1500 --popSize 7947 -Td 0 0.032968 -Td 0.005285 26.79 -Td 0.328237 7.537683 -TE 0.03624009 -TW 0 2 0 1 1 0.206 76.50 -TW 0.005285 2 0 1 1 0.206 2049.5 -TW 0.328237 2 0 1 1 0.206 15400 --sampSize 100000 --outfile out.txt --errfile err.txt --popFreq freq.txt &


# Ryan's recommendation
./sfs_code 1 20 -t 0.0003773252 -Td 0 0.7218 -Td 0.4878 5.2693 -TE 0.5432 -o sfs_code.txt -W 2 0 0 0 0.184 0.00040244


-W 2 0 0 0 0.184 0.00040244 #  sfs_code uses the version with mean = alpha/beta; 1/(0.16*7778*2) = 0.00040244 where 7778 is number of ancestry population for the African model; Other parameters (see the example on SFS documentation and Boyko2008) Nanc=7895, Nbt=5699, Ncurr=30030, expansion time:874gen ago; bottleneck duration 7703gen

-Td 0 0.7218 # 5699/7895 = 0.7218
-Td 0.4878 5.2693 # 7703/7895/2 = 0.4878404, 30030/5699 = 5.269345
-TE 0.5432 # 874/7895/2 + 0.4878 = 0.5431515
-t 0.0003773252 # 7895*1.2e-8*4 = 0.00037896

dir="/home/$USER/Toolkit/sfscode/bin"
# Using Boyko's model -- my simulation
$dir/sfs_code 1 500 -t 0.0003773252 -L 1 1500 --popSize 7895 -Td 0 0.7218 -Td 0.4878 5.2693 -TE 0.5432 -W 2 0 0 0 0.184 0.00040244 --sampSize 30030 --outfile boyko.out --errfile boyko.err --popFreq boyko.freq &

# MAF calculated from 30030 sample size seems perfectly agree with population MAF


# Using Shamil's model -- my simulation
$dir/sfs_code 1 500 -t 0.000384 -L 1 1500 --popSize 8000 -Td 0 0.9875 -Tg 0.000625 112.5 -TE 0.02375 -W 2 0 0 0 0.184 0.000390625 --sampSize 100000 --outfile kyrukov.out --errfile kyrukov.err --popFreq kyrukov.freq &
# MAF calculated from 200000 sample size would be even smaller than from the population (supposedly the population should be 900000 samples), this is very strange. Using 100000 for now.


-Td 0 0.9875 # 7900/8000
-Tg 0.000625 112.5 # 10/8000/2 = 0.000625, 900000/8000 = 112.5
-TE 0.02375 # 370/8000/2 + 0.000625 = 0.02375
-t 0.000384 # = 1.2e-8*4*8000
#-W 2 0 0 0 0.01 0.0001111425 # 1/(0.562341*8000*2)
-W 2 0 0 0 0.184 0.000390625 # 1/(0.16*8000*2) = 0.000390625


# Using Shamil's model with protective variants ... my simulation

# 25%
$dir/sfs_code 1 500 -t 0.000384 -L 1 1500 --popSize 8000 -Td 0 0.9875 -Tg 0.000625 112.5 -TE 0.02375 -W 2 0.25 0.184 0.000390625 0.184 0.000390625 --sampSize 100000 --outfile kyrukovp25.out --errfile kyrukovp25.err --popFreq kyrukovp25.freq &

# 50%
$dir/sfs_code 1 500 -t 0.000384 -L 1 1500 --popSize 8000 -Td 0 0.9875 -Tg 0.000625 112.5 -TE 0.02375 -W 2 0.5 0.184 0.000390625 0.184 0.000390625 --sampSize 100000 --outfile kyrukovp50.out --errfile kyrukovp50.err --popFreq kyrukovp50.freq &

# 75%
$dir/sfs_code 1 500 -t 0.000384 -L 1 1500 --popSize 8000 -Td 0 0.9875 -Tg 0.000625 112.5 -TE 0.02375 -W 2 0.75 0.184 0.000390625 0.184 0.000390625 --sampSize 100000 --outfile kyrukovp75.out --errfile kyrukovp75.err --popFreq kyrukovp75.freq &

# 100%
$dir/sfs_code 1 500 -t 0.000384 -L 1 1500 --popSize 8000 -Td 0 0.9875 -Tg 0.000625 112.5 -TE 0.02375 -W 2 1 0.184 0.000390625 0 0 --sampSize 100000 --outfile kyrukovp100.out --errfile kyrukovp100.err --popFreq kyrukovp100.freq &

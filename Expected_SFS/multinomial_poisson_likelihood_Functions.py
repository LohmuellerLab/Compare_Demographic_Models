# November 28, 2016

# Author: Beichman, Annabel
# Contact: annabel.beichman@gmail.com
# Purpose: Python code to calculate multinomial and Poisson log-likelihoods
import math

# numHaps is number of haploid individuals in SFS
# model_SFS_freq: SNP frequencies in model SFS
# model_SFS_count: SNP counts in model SFS
# obs_SNP_counts: SNP counts in observed (empirical) SFS

########## Multinomial Log  Likelihood ##########

numHaps = 20 # haploid individuals
def LhoodCalc(model_SFS_freq,obs_SNP_counts,numHaps):
    llSum=0
    for i in range(1,numHaps):
        llSum += obs_SNP_counts[i]*math.log(model_SFS_freq[i])
    return llSum

########## POISSON Log LIKELIHOODS (see Lohmueller 2008) ##########

numHaps = 20 # haploid individuals
def LhoodCalcPoisson(model_SFS_count,obs_SNP_counts,numHaps):
    llSum=0
    for i in range(1,numHaps):
        llSum += obs_SNP_counts[i]*math.log(model_SFS_count[i])
    ll = -np.sum(model_SFS_count) + llSum
    return ll


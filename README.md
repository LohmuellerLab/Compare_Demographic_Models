# Compare_SFS_SMC
Scripts associated with Beichman &amp; Lohmueller (2016), in preparation, comparing SFS and SMC-based demographic inference methods

The following scripts are in this repository:

(1) getExpectedSFS_fromMSMC_inDadi_2_4_8_Haplotypes_Nov28_v1.0.py : Python script to calculate the expected SFS under the demographic histories inferred using MSMC by Schiffels & Durbin (2014)

(2) multinomial_poisson_likelihood_Functions.py : Python script to calculate multinomical and Poisson log likelihood of expected SFSs

(3) simulateGutenkunstExpSFS_OutOfAfrica.py : a modification of Gutenkunst et al.'s (2010) code to calculate the expected SFS under their Out of Africa Model

(4) SimulateGutenkunst_OutOfAfrica_inMaCS.sh : bash script to use MaCS (Chen et al. 2009) to simulate data under the Gutenkunst Out of Africa Model 

(5) SimulateGutenkunst_OutofAfrica_plusExplosiveGrowth_inMaCS_v1.0_Nov28.sh : bash script to use MaCS (Chen et al. 2009) to simulate data under the Gutenkunst Out of Africa Model plus Explosive growth based on parameters from Tennessen et al. (2012). 

(6) Simulate_GutenkunstPlusNeanderthal_inMSMS_v1.0_Nov28.sh : bash script to use MaCS (Chen et al. 2009) to simulate data under the Gutenkunst Out of Africa Model plus Neanderthal admixture based on parameters from Harris & Nielsen (2016).

(7) Simulate_MSMC2HapModels_inMaCS_v1.0_Nov28.sh : bash script to use MaCS (Chen et al. 2009) to simulate data under the demographic models inferred using 2 haplotypes in MSMC


# Compare\_Demographic_Models
### These are the scripts associated with Beichman, Phung & Lohmueller (2017): 
### "Comparison of single genome and allele frequency data reveals discordant demographic histories"
### biorxiv: 


##### _These scripts evaluate the fit of the following published demographic models to empirical data for the following three statistics:_

(1) genome-wide and neutral heterozygosity 

(2) linkage disequilibrium decay

(3) expected site frequency spectrum (SFS)

#### Demographic Models (CEU, CHB, YRI): 

* SFS-based (dadi): Gutenkunst et al. (2009)
* MSMC: Li & Durbin (2011) 
* SMC++: Terhorst, Kamm & Song 2017)

#### Scripts to simulate data in MaCS (Chen et al. 2010) under published demographic models for CEU, CHB, YRI populations

###### **Simulate\_100kb_Blocks/**: Simulate 20,000 x 100kb sequence blocks for comparison with whole genome 1000 Genomes Data (includes file of recombination rates for 100kb windows based on the deCode project)

###### **Simulate\_10kb_Blocks/**: Simulate 6300 x 10kb sequence blocks for comparison with neutral 1000 Genomes Data (includes file of recombination rates based on empirical 10kb neutral windows)

###### **Simulate\_80x30Mb_chromosomes/**: Simulate 80 x 30Mb "chromosomes" to use as input for MSMC (includes file of recombination rates for 100kb windows based on the deCode project)

#### Scripts to get expected site frequency spectrum (SFS) in dadi for each published demographic model 

###### **Expected\_SFS/**: use dadi (Gutenkunst et al. 2009) to get the expected SFS under each demographic model


#### Scripts to calculate linkage disequilibrium decay from simulated data

###### **Calculate\_LD\_Decay\_PhysicalDistance/**: scripts to convert MaCS output to vcf format, calcuate genotype-r2 using vcftools, bin r2 values by physical distance between SNPs, and average. 

#### Scripts to calculate linkage disequilibrium decay from empirical 1000 Genomes data

###### **Tanya Phung's scripts**: https://github.com/tnphung/1000G\_Summary_Stats

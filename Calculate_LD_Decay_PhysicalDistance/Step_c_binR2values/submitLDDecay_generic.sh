#! /bin/bash
#$ -cwd
#$ -l h_rt=5:00:00,h_data=2G,highp
#$ -N smcppLD
#$ -o VCFoutput.txt
#$ -e VCFerror.txt
#$ -m abe
#$ -M ab08028
#$ -t 1-200 


# Feb 11, 2017
# author: Annabel Beichman
# contact: annabel.beichman@gmail.com
# purpose: This script will bin r2 values (output from vcftools geno-r2) by the physical distance
# between the SNPs (you can set number of bins)

source /u/local/Modules/default/init/modules.sh
module load python 

# estimateLDdecay.ABnanEdit.py must be in the folder!

j=$SGE_TASK_ID
echo $j
flag= # model name, signifier
ids= # population ids 

rundate=20170214
for id in $ids
do
for i in {1..100}
do
echo $id
# usage python estimateLDdecay.ABnanEdit.py --input <input vcftools output> --bin <# of bins you want + 1> --outfile <output file name> 
python estimateLDdecay.ABnanEdit.py --input group_${j}_out/group_${j}_block_${i}.${flag}.${id}.vcftools.LD.R2.output.geno.ld --format vcftools --bin 101 --outfile group_${j}_out/group_${j}_block_${i}.${flag}.${id}.vcftools.LD.R2.binned.txt
done
done

sleep 10m


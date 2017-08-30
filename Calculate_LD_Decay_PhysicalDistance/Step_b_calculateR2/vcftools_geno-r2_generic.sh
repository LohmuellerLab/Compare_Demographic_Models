#! /bin/bash
#$ -cwd
#$ -l h_rt=05:00:00,h_data=10G,highp
#$ -N gutvcf
#$ -o /u/scratch/a/ab08028/vcftoolsReports_GUT
#$ -e /u/scratch/a/ab08028/vcftoolsReports_GUT
#$ -m abe
#$ -M ab08028
#$ -t 1-200

# Feb 11, 2017
# author: Annabel Beichman
# contact: annabel.beichman@gmail.com
# purpose: Run vcftools geno-r2 on vcf files of 20,000 x 100kb simulated blocks under a demographic model


source /u/local/Modules/default/init/modules.sh
module load vcftools/0.1.14

## make sure ##info line is gone
#sed -i.bak -e '2d' group_*.output.vcf
### want to run vcf tools
# run from vcfs dir
flag= # if something else is in model  name
ids= # population IDs e.g. "yri ceu chb"
rundate= # date of your simulations
j=$SGE_TASK_ID

mkdir $SCRATCH/${flag}_vcftoolsOutput_blockwise
out=$SCRATCH/${flag}_vcftoolsOutput_blockwise
# infile : run script from inside vcf folder in home directory:
gpin=group_${j}_${flag}_vcfs_${rundate}
mkdir $out/group_${j}_out
gpout=$out/group_${j}_out

for id in $ids
do
for i in {1..100}
do

vcftools --mac 2 --vcf $gpin/group_${j}_block_${i}.$flag.$id.output.vcf --geno-r2 --ld-window-bp 100000 --out $gpout/group_${j}_block_${i}.$flag.$id.vcftools.LD.R2.output

done
done

echo "finished ${SGE_TASK_ID}"

sleep 20m


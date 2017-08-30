#! /bin/bash
#$ -l h_rt=24:00:00,h_data=28G
#$ -N gutsim
#$ -o output.txt
#$ -e error.txt
#$ -m abe
#$ -M ab08028
#$ -t 1-200

# November 28, 2016
# Author: Beichman, Annabel  annabel.beichman@gmail.com
# Contact: annabel.beichman@gmail.com
# This script simulates 20,000 x 100kb windows across 10 individuals under the Gutenkunst et al. (2009) models for CEU, CHB and YRI

j=$SGE_TASK_ID
rundate=`date +%Y%m%d`
model=gutenkunst_inMACS_model_100kb_10ind
recRatesFile=allChrRecRates_Mperbp_100kbWindows.forMACS.20170201.txt
mkdir group_$j.${model}
cd group_$j.${model}
cp ../macs ./
cp ../msformatter ./
mkdir africanSample
mkdir europeanSample
mkdir asianSample
for i in {1..100}
do
window=$((((${j}-1)*100)+$i))
r=`awk -v window=$window '{if(NR==window) print}' $recRatesFile`
ri=`bc <<< "${r}*7310*4"`
# YRI gutenkunst sample
./macs 20 100000 -t 0.000690 -r $ri  -I 3 20 0 0 -n 1 1.682020 -n 2 3.736830 -n 3 7.292050 -eg 0 2 116.010723 -eg 0.0000000001 3 160.246047 -ma x 0.881098 0.561966 0.881098 x 2.797460 0.561966 2.797460 x -ej 0.028985 3 2 -en 0.02898501 2 0.287184 -ema 0.02898502 3 x 7.293140 x 7.293140 x x x x x -ej 0.197963 2 1 -en 0.303501 1 1  > africanSample/group_${j}_block_${i}.gutenkunst.yri.macs.macsFormat.OutputFile.${rundate}.txt 
./msformatter < africanSample/group_${j}_block_${i}.gutenkunst.yri.macs.macsFormat.OutputFile.${rundate}.txt > africanSample/group_${j}_block_${i}.gutenkunst.yri.macs.msFormat.OutputFile.${rundate}.txt
# CEU gutenkunst sample
./macs 20 100000 -t 0.000690 -r $ri  -I 3 0 20 0 -n 1 1.682020 -n 2 3.736830 -n 3 7.292050 -eg 0 2 116.010723 -eg 0.0000000001 3 160.246047 -ma x 0.881098 0.561966 0.881098 x 2.797460 0.561966 2.797460 x -ej 0.028985 3 2 -en 0.02898501 2 0.287184 -ema 0.02898502 3 x 7.293140 x 7.293140 x x x x x -ej 0.197963 2 1 -en 0.303501 1 1  > europeanSample/group_${j}_block_${i}.gutenkunst.ceu.macs.macsFormat.OutputFile.${rundate}.txt 
./msformatter < europeanSample/group_${j}_block_${i}.gutenkunst.ceu.macs.macsFormat.OutputFile.${rundate}.txt > europeanSample/group_${j}_block_${i}.gutenkunst.ceu.macs.msFormat.OutputFile.${rundate}.txt
# CHB gutenkunst sample
./macs 20 100000 -t 0.000690 -r $ri  -I 3 0 0 20 -n 1 1.682020 -n 2 3.736830 -n 3 7.292050 -eg 0 2 116.010723 -eg 0.0000000001 3 160.246047 -ma x 0.881098 0.561966 0.881098 x 2.797460 0.561966 2.797460 x -ej 0.028985 3 2 -en 0.02898501 2 0.287184 -ema 0.02898502 3 x 7.293140 x 7.293140 x x x x x -ej 0.197963 2 1 -en 0.303501 1 1  > asianSample/group_${j}_block_${i}.gutenkunst.chb.macs.macsFormat.OutputFile.${rundate}.txt 
./msformatter < asianSample/group_${j}_block_${i}.gutenkunst.chb.macs.macsFormat.OutputFile.${rundate}.txt > asianSample/group_${j}_block_${i}.gutenkunst.chb.macs.msFormat.OutputFile.${rundate}.txt
done
cd ../

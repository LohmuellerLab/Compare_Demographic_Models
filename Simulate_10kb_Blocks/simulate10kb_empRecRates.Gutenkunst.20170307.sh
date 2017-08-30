# This script can be run in the shell and will simulate 6300 x 10kb blocks (to match the empirical neutral windows identified with the NRE)
# under the Gutenkunst et al. 2009 out of Africa models

rundate=`date +%Y%m%d`
model=gutenkunst_inMACS_model_10kb_10ind
recRatesFile=/u/home/a/ab08028/klohmueldata/annabel_data/psmc_project/hetLDComparisons_2017/10kb_window_simulations_20170307/Mbp_recRates_10kbWindows_6333_20170307.txt
for j in {1..10}
do
mkdir group_$j.${model}
cd group_$j.${model}
cp ../macs ./
cp ../msformatter ./
mkdir africanSample
mkdir europeanSample
mkdir asianSample
for i in {1..630}
do
window=$((((${j}-1)*630)+$i))
r=`awk -v window=$window '{if(NR==window) print}' $recRatesFile`
ri=`bc <<< "${r}*7310*4"`
# YRI gutenkunst sample
./macs 20 10000 -t 0.000690 -r $ri  -I 3 20 0 0 -n 1 1.682020 -n 2 3.736830 -n 3 7.292050 -eg 0 2 116.010723 -eg 0.0000000001 3 160.246047 -ma x 0.881098 0.561966 0.881098 x 2.797460 0.561966 2.797460 x -ej 0.028985 3 2 -en 0.02898501 2 0.287184 -ema 0.02898502 3 x 7.293140 x 7.293140 x x x x x -ej 0.197963 2 1 -en 0.303501 1 1  > africanSample/group_${j}_block_${i}.gutenkunst.yri.macs.macsFormat.OutputFile.${rundate}.txt 
./msformatter < africanSample/group_${j}_block_${i}.gutenkunst.yri.macs.macsFormat.OutputFile.${rundate}.txt > africanSample/group_${j}_block_${i}.gutenkunst.yri.macs.msFormat.OutputFile.${rundate}.txt
# CEU gutenkunst sample
./macs 20 10000 -t 0.000690 -r $ri  -I 3 0 20 0 -n 1 1.682020 -n 2 3.736830 -n 3 7.292050 -eg 0 2 116.010723 -eg 0.0000000001 3 160.246047 -ma x 0.881098 0.561966 0.881098 x 2.797460 0.561966 2.797460 x -ej 0.028985 3 2 -en 0.02898501 2 0.287184 -ema 0.02898502 3 x 7.293140 x 7.293140 x x x x x -ej 0.197963 2 1 -en 0.303501 1 1  > europeanSample/group_${j}_block_${i}.gutenkunst.ceu.macs.macsFormat.OutputFile.${rundate}.txt 
./msformatter < europeanSample/group_${j}_block_${i}.gutenkunst.ceu.macs.macsFormat.OutputFile.${rundate}.txt > europeanSample/group_${j}_block_${i}.gutenkunst.ceu.macs.msFormat.OutputFile.${rundate}.txt
# CHB gutenkunst sample
./macs 20 10000 -t 0.000690 -r $ri  -I 3 0 0 20 -n 1 1.682020 -n 2 3.736830 -n 3 7.292050 -eg 0 2 116.010723 -eg 0.0000000001 3 160.246047 -ma x 0.881098 0.561966 0.881098 x 2.797460 0.561966 2.797460 x -ej 0.028985 3 2 -en 0.02898501 2 0.287184 -ema 0.02898502 3 x 7.293140 x 7.293140 x x x x x -ej 0.197963 2 1 -en 0.303501 1 1  > asianSample/group_${j}_block_${i}.gutenkunst.chb.macs.macsFormat.OutputFile.${rundate}.txt 
./msformatter < asianSample/group_${j}_block_${i}.gutenkunst.chb.macs.macsFormat.OutputFile.${rundate}.txt > asianSample/group_${j}_block_${i}.gutenkunst.chb.macs.msFormat.OutputFile.${rundate}.txt
done
cd ../
done

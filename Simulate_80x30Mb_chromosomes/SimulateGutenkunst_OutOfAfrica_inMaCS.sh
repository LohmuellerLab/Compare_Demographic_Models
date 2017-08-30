#!/bin/bash/
# November 28, 2016
# Author: Beichman, Annabel  annabel.beichman@gmail.com
# Contact: annabel.beichman@gmail.com
# Purpose: simulate 80 x 30Mbp independent segments of genetic sequence 
# in MaCS (Chen et al. 2008) under Gutenkunst et al.Õs (2009) Out of Africa model
# of human demography (Gutenkunst et al. 2009 Supplementary Information)

# make sure python3 is available
# Requires ms2multihetsep.py https://github.com/stschiff/msmc-tools/blob/master/ms2multihetsep.py
# and 
# msformatter https://github.com/gchen98/macs
# to be in your starting folder

# sequence length: (30Mb)
seqLength=30000000

# Path to recombination rate file (relative to average recombination rate)
# (or can just use a single recombination rate)
recRatesFile=[/path/to/recRates/file]

# Theta: 2788.2
# L: 4.04x10^6
# Nref: 7310
# mu: 2.35 x 10 ^ -8
# So theta = 4 N u L 
# but MaCS wants 4*N*u, so
# divide THETA / L -> 2788.2 / 4.04 x 10^6 = 6.90E-04

## AVERAGE (baseline) REC RATE:
# 9.257856e-09

for j in {1..50}
do
mkdir ${j}.gutenkunstMacsSim
cd ${j}.gutenkunstMacsSim
cp ../macs ./
cp ../ms2multihetsep.py ./
cp ../msformatter ./
mkdir africanSample
mkdir europeanSample
mkdir asianSample

	for i in {1..80}
	do
	./macs 2 $seqLength -t 0.000690 -r 0.00027 -R ${recRatesFile}/${i}.300blocks.RATIOS.MACS.txt -I 3 2 0 0 -n 1 1.682020 -n 2 3.736830 -n 3 7.292050 -eg 0 2 116.010723 -eg 0.0000000001 3 160.246047 -ma x 0.881098 0.561966 0.881098 x 2.797460 0.561966 2.797460 x -ej 0.028985 3 2 -en 0.02898501 2 0.287184 -ema 0.02898502 3 x 7.293140 x 7.293140 x x x x x -ej 0.197963 2 1 -en 0.303501 1 1| ./msformatter| python3 ms2multihetsep.py $i $seqLength > africanSample/${i}.af.msmcFormat.OutputFile.txt

	# sampling Europeans
	./macs 2 $seqLength -t 0.000690 -r 0.00027 -R ${recRatesFile}/${i}.300blocks.RATIOS.MACS.txt -I 3 0 2 0 -n 1 1.682020 -n 2 3.736830 -n 3 7.292050 -eg 0 2 116.010723 -eg 0.0000000001 3 160.246047 -ma x 0.881098 0.561966 0.881098 x 2.797460 0.561966 2.797460 x -ej 0.028985 3 2 -en 0.02898501 2 0.287184 -ema 0.02898502 3 x 7.293140 x 7.293140 x x x x x -ej 0.197963 2 1 -en 0.303501 1 1| ./msformatter| python3 ms2multihetsep.py $i $seqLength > europeanSample/${i}.eu.msmcFormat.OutputFile.txt

	# sampling Asians
	./macs 2 $seqLength -t 0.000690 -r 0.00027 -R ${recRatesFile}/${i}.300blocks.RATIOS.MACS.txt -I 3 0 0 2 -n 1 1.682020 -n 2 3.736830 -n 3 7.292050 -eg 0 2 116.010723 -eg 0.0000000001 3 160.246047 -ma x 0.881098 0.561966 0.881098 x 2.797460 0.561966 2.797460 x -ej 0.028985 3 2 -en 0.02898501 2 0.287184 -ema 0.02898502 3 x 7.293140 x 7.293140 x x x x x -ej 0.197963 2 1 -en 0.303501 1 1| ./msformatter| python3 ms2multihetsep.py $i $seqLength > asianSample/${i}.asn.msmcFormat.OutputFile.txt

	done
cd ../
done
#!/bin/bash/
# November 28, 2016
# Author: Beichman, Annabel  annabel.beichman@gmail.com
# Contact: annabel.beichman@gmail.com
# Purpose: Bash code to simulate 80 x 30Mbp independent segments of 
# genetic sequence in MaCS (Chen et al. 2008)
# under Gutenkunst et al.Õs (2009) Out of Africa model of human demography 
# (Gutenkunst et al. 2009, Supplementary Information) 
# modified to include Neanderthal admixture as in Harris & Nielsen (2016)

# recommend running as an array on a cluster to improve speed
# Requires ms2multihetsep.py https://github.com/stschiff/msmc-tools/blob/master/ms2multihetsep.py

scriptDir=[path to script]
msms=[path to MSMS]
outputDir=[path to output directory]

# sequence length: (30Mb)
# macs theta and r are scaled by 4NA, but in msms need to also be scaled by LENGTH 
# L =  30,000,000 (30MB)
# t = 0.000690 * 30,000,000(length)  = 20700
# r = 0.00027 (average from Tanya's files) 
# used r = 1e-08 M per base (slightly higher than rec rate used in MaCS simulations)
# so r = 1e-08 * 4 *7310 = 0.0002924 >> * 30,000,000(Length) = 8,772

seqLength=30000000
model="neanderthal" 
cd $outputDir
for j in {1..50}
do
	mkdir ${j}.${model}.MsPrimeSim
	cd ${j}.${model}.MsPrimeSim
	mkdir africanSample
	mkdir europeanSample
	mkdir asianSample
	for i in {1..80}
		do 
		# sample african
		$msms -ms 2 1 -oFP 0.0000000000E00 -threads 4 -t 20700 -r 8772 -I 3 2 0 0 -n 1 1.682020 -n 2 3.736830 -n 3 7.292050 -eg 0 2 116.010723 -eg 0 3 160.246047 -ma x 0.881098 0.561966 0.881098 x 2.797460 0.561966 2.797460 x -ej 0.028985 3 2 -en 0.028985 2 0.287184 -ema 0.028985 3 x 7.293140 x 7.293140 x x x x x -es 0.0683995 2 0.98 -en 0.0683995 4 0.1367989056 -en 0.0683995 2 0.287184 -ej 0.197963 2 1 -en 0.303501 1 1 -ej 0.6155950752 4 1 -en 0.615596 1 1 | python3 $scriptDir/ms2multihetsep.py $i $seqLength > africanSample/${i}.af.msms.msmcFormat.OutputFile.txt
		# sample asian
		$msms -ms 2 1 -oFP 0.0000000000E00 -threads 4 -t 20700 -r 8772 -I 3 0 0 2 -n 1 1.682020 -n 2 3.736830 -n 3 7.292050 -eg 0 2 116.010723 -eg 0 3 160.246047 -ma x 0.881098 0.561966 0.881098 x 2.797460 0.561966 2.797460 x -ej 0.028985 3 2 -en 0.028985 2 0.287184 -ema 0.028985 3 x 7.293140 x 7.293140 x x x x x -es 0.0683995 2 0.98 -en 0.0683995 4 0.1367989056 -en 0.0683995 2 0.287184 -ej 0.197963 2 1 -en 0.303501 1 1 -ej 0.6155950752 4 1 -en 0.615596 1 1 | python3 $scriptDir/ms2multihetsep.py $i $seqLength > asianSample/${i}.asn.msms.msmcFormat.OutputFile.txt
		# sample european
		$msms -ms 2 1 -oFP 0.0000000000E00 -threads 4 -t 20700 -r 8772 -I 3 0 2 0 -n 1 1.682020 -n 2 3.736830 -n 3 7.292050 -eg 0 2 116.010723 -eg 0 3 160.246047 -ma x 0.881098 0.561966 0.881098 x 2.797460 0.561966 2.797460 x -ej 0.028985 3 2 -en 0.028985 2 0.287184 -ema 0.028985 3 x 7.293140 x 7.293140 x x x x x -es 0.0683995 2 0.98 -en 0.0683995 4 0.1367989056 -en 0.0683995 2 0.287184 -ej 0.197963 2 1 -en 0.303501 1 1 -ej 0.6155950752 4 1 -en 0.615596 1 1 | python3 $scriptDir/ms2multihetsep.py $i $seqLength > europeanSample/${i}.eur.msms.msmcFormat.OutputFile.txt
		done
cd ../
done
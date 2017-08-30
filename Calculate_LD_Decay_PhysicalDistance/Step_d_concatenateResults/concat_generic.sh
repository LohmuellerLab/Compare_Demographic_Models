# Feb 11, 2017
# author: Annabel Beichman
# contact: annabel.beichman@gmail.com
# purpose: This script will concatenate the results of binning your r2 values by physical distance

flag= # model 
ids= # pop ids, e.g. "yri ceu chb"
date=`date +%Y%m%d`
for id in $ids
do
echo $id
> $flag.$id.ALL.LD.BINNED.OUTPUT.${date}.txt
for j in {1..200}
do
echo $j
for i in {1..100}
do
cat group_${j}_out/group_${j}_block_${i}.$flag.${id}.vcftools.LD.R2.binned.txt >> $flag.$id.ALL.LD.BINNED.OUTPUT.${date}.txt
done
done
done

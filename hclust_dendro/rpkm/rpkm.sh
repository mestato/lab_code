#!/bin/bash
./calculate_rpkms.pl
mkdir tmp
mv *.rpkm tmp
cut -f 1 tmp/24_GA-CL.rpkm > tmp/rows
for file in `ls tmp/*.rpkm`
do
cut -f 2 $file > $file.cut
done
paste tmp/rows tmp/*.cut > tmp/all.tsv
cat labels.txt tmp/all.tsv > all_labeled.tsv
cut -f 1-55 all_labeled.tsv > all_final_cut55.tsv
rm tmp/*.rpkm tmp/*.cut tmp/all.tsv tmp/rows
rmdir tmp

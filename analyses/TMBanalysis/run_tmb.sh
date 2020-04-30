#!/bin/bash

set -e
set -o pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

# Step1 changing MAF to BE format 
# Creating BED file format from MAF file so that we can use 
# bedtools to restrict  variants only within the BED file 
input_BED_file=$(echo $1  | sed 's/.maf/.bed/g')
tail -n +3 $1 | cut -f 5- > $input_BED_file

# Step2 intersecting MAFBED with regionsBED 
mafBED_within_target=$(echo $input_BED_file | sed 's/.bed/.withintarget.bed/g')
bedtools intersect -a  $input_BED_file -b  $2 > $mafBED_within_target

#Step3 header to BED file that  is within the target region
cat inputs/header > temp
cat $mafBED_within_target >> temp
mv temp $mafBED_within_target

#Step4 calculate  TMB scores  
outfile=$(echo $mafBED_within_target | rev  | cut -d "/" -f1 | rev| sed 's/.withintarget.bed/.tmbscores.txt/g')
python calculate_TMBscores.py -m $mafBED_within_target -i $3 -o $outfile 



#!/bin/bash

#Download commonly mutated genes

set -e
set -o pipefail

#parse options
while (( "$#" )); do #loop through each option given
  case "$1" in
    -o|--output) #path to output directory
    if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
      outputDir=$2
      shift 2
    else
      echo "Error: Argument for $1 is missing" >&2
      exit 1
    fi
    ;;
    -*|--*=) #unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
  esac
done

#exclude the 50 most commonly mutated genes
excludeUri='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5706417/bin/12920_2017_309_MOESM3_ESM.txt'
excludeFile=./exclude-genes.txt
echo -e "gene\tcount" > $excludeFile
head -n 50 <(curl -s ${excludeUri}) >> $excludeFile

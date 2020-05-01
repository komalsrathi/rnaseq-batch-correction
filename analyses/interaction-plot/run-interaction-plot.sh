#!/bin/bash

#Generate the coocurrence and  disease gene plots

set -e
set -o pipefail

#TODO / Questions
#    change max genes to be option?
#    option to only generate cooc OR disease gene plot?
#    default to not downloading gene exclusion file?
#    option to specify which diseases to plot

#default options
maxGenes=50
vaf=0.05
outputDir=$PWD/outputs
scriptDir=$PWD/scripts

#set up diseases
declare -A disease
disease[all]="All"
disease[Medulloblastoma]="Medulloblastoma"
disease[LGAT]="Low-grade astrocytic tumor"
disease[Ependymoma]="Ependymoma"
disease[HGAT]="High-grade glioma"
disease[DMG]="Diffuse midline glioma"
disease[Ganglioglioma]="Ganglioglioma"
disease[Craniopharyngioma]="Craniopharyngioma"

#parse options
while (( "$#" )); do #loop through each option given
  case "$1" in
    -p|--path) #path to scripts folder, defuault to ./scripts/ dir
    if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
      scriptDir=$2
      shift 2
    else
      echo "Error: Argument for $1 is missing" >&2
      exit 1
    fi
    ;;
    -o|--output) #path to output directory
    if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
      outputDir=$2
      shift 2
    else
      echo "Error: Argument for $1 is missing" >&2
      exit 1
    fi
    ;;
    -s|--samples) #tsv with sample ids and ???
    if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
      samples=$2
      shift 2
    else
      echo "Error: Argument for $1 is missing" >&2
      exit 1
    fi
    ;;
    -d|--metadata) #tsv with sample metadata
    if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
      metadata=$2
      shift 2
    else
      echo "Error: Argument for $1 is missing" >&2
      exit 1
    fi
    ;;
    -m|--maf) #maf file
    if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
      maf=$2
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

#check that samples, metadata, and maf were given
if [ -z "$samples" ]; then
  echo "Error: samples file not given"
fi
if [ -z "$metadata" ]; then
  echo "Error: metadata file not given"
fi
if [ -z "$maf" ]; then
  echo "Error: maf file not given"
fi

#figure out script path, options, inputs, and outputs based on options
filesDir=$outputDir/files
figureDir=$outputDir/figures
mkdir -p $figureDir
mkdir -p $filesDir

#output files basenames
cooccur=${filesDir}/cooccur
geneDisease=${filesDir}/gene_disease
coocPlot=${figureDir}/cooccur
diseasePlot=${figureDir}/gene_disease
combinedPlot=${figureDir}/combined

#exclude the 50 most commonly mutated genes
excludeUri='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5706417/bin/12920_2017_309_MOESM3_ESM.txt'
excludeFile=${filesDir}/FLAGS.tsv
echo -e "gene\tcount" > $excludeFile
head -n 50 <(curl -s ${excludeUri}) >> $excludeFile

#run scripts
for diseaseId in "${!disease[@]}"; do
  echo $diseaseId

  Rscript ${scriptDir}/01-disease-specimen-lists.R \
    --metadata ${metadata} \
    --specimen_list ${samples} \
    --disease "${disease[$diseaseId]}" \
    --outfile ${filesDir}/${diseaseId}.tsv

  Rscript ${scriptDir}/02-process_mutations.R \
    --maf ${maf} \
    --metadata ${metadata} \
    --specimen_list ${filesDir}/${diseaseId}.tsv \
    --exclude_genes $excludeFile \
    --vaf $vaf \
    --min_mutated 5 \
    --max_genes $maxGenes \
    --disease_table ${geneDisease}.${diseaseId}.tsv \
    --out ${cooccur}.${diseaseId}.tsv

  Rscript ${scriptDir}/03-plot_interactions.R \
    --infile ${cooccur}.${diseaseId}.tsv \
    --outfile ${coocPlot}.${diseaseId}.png \
    --disease_table ${geneDisease}.${diseaseId}.tsv \
    --disease_plot ${diseasePlot}.${diseaseId}.png \
    --combined_plot ${combinedPlot}.${diseaseId}.png \
    --plotsize 50
done


import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--maffile', required = True,
                    help = 'path to the histology file')
parser.add_argument('-i', '--histfile', required = True,
                    help = 'path to the histology file')
parser.add_argument('-o', '--outfile', required = True,
                    help = "output notebook")
args = parser.parse_args()

maf_file=pd.read_csv(args.maffile, sep="\t", index_col=False)
maf_file = maf_file.fillna("NA")

histology = pd.read_csv(args.histfile, sep="\t", index_col=False)

sample_tmb_colnames = ["Samplenames", "TMB"]
sample_tmb_df = pd.DataFrame(columns=sample_tmb_colnames)

var_by_sample = maf_file.groupby("Tumor_Sample_Barcode")
for name in var_by_sample:
    sample_name = name[0]
    affected_var = name[1].loc[(name[1]['Variant_Classification'] == 'Missense_Mutation')|( name[1]['Variant_Classification'] ==  'Nonsense_Mutation')]
    count = affected_var.shape[0]
    sample_tmb_df = sample_tmb_df.append({"Samplenames" : sample_name , "TMB" : str(count*1000000/77462866)} , ignore_index=True)
                                          
final_tmb = sample_tmb_df.join(histology.set_index("Kids_First_Biospecimen_ID"), on="Samplenames")[["Samplenames",  "TMB", "short_histology"]]
final_tmb.to_csv(args.outfile, index=False)





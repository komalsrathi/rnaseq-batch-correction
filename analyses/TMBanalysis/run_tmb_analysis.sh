



python3  code/01_setup_mafdb.py \
	-m /Users/kogantit/Documents/TMB/temp.maf  \
	-d /Users/kogantit/Documents/git_repos/d3b-bix-analysis-toolkit/analyses/TMBanalysis/output/sample_var_db.sqlite \
	-t /Users/kogantit/Documents/TMB/xgen-exome-research-panel-targets_hg38_ucsc_liftover.100bp_padded.sort.merged.bed \
	-c /Users/kogantit/Documents/TMB/gencode.v34.annotation.exononly.bed \
	-b mutect2 \
	-i /Users/kogantit/Documents/git_repos/d3b-bix-analysis-toolkit/analyses/TMBanalysis/xgen_intersected_cds.bed


python3 code/02_calculate_tmb_from_mafdb.py \
	-m /Users/kogantit/Documents/TMB/pbta-histologies.tsv \
	-c short_histology \
	-d /Users/kogantit/Documents/git_repos/d3b-bix-analysis-toolkit/analyses/TMBanalysis/output/sample_var_db.sqlite \
	-s Kids_First_Biospecimen_ID \
	-b mutect2 \
	-o /Users/kogantit/Documents/git_repos/d3b-bix-analysis-toolkit/analyses/TMBanalysis/output/sample_out_tmbscores.txt \
	-i /Users/kogantit/Documents/git_repos/d3b-bix-analysis-toolkit/analyses/TMBanalysis/output/xgen_intersected_cds.bed 


python3 code/03_tmbplots.py  \
	-t /Users/kogantit/Documents/git_repos/d3b-bix-analysis-toolkit/analyses/TMBanalysis/output/pbta-snv-mutect2-tmbscores.target.txt \
	-o /Users/kogantit/Documents/git_repos/d3b-bix-analysis-toolkit/analyses/TMBanalysis/output/temp.png




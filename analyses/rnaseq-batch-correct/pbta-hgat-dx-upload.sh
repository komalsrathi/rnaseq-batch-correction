# uncorrected files
aws s3 --profile saml cp output/pbta-hgat-dx-gene-expression-rsem-tpm-uncorrected.rds s3://d3b-bix-dev-data-bucket/hgg-dmg-integration/merged_ngs_files/rna-seq/batch-corrected/
aws s3 --profile saml cp output/pbta-hgat-dx-gene-expression-rsem-fpkm-uncorrected.rds s3://d3b-bix-dev-data-bucket/hgg-dmg-integration/merged_ngs_files/rna-seq/batch-corrected/

# corrected files
aws s3 --profile saml cp output/pbta-hgat-dx-gene-expression-rsem-tpm-corrected.rds s3://d3b-bix-dev-data-bucket/hgg-dmg-integration/merged_ngs_files/rna-seq/batch-corrected/
aws s3 --profile saml cp output/pbta-hgat-dx-gene-expression-rsem-fpkm-corrected.rds s3://d3b-bix-dev-data-bucket/hgg-dmg-integration/merged_ngs_files/rna-seq/batch-corrected/

# combined files
aws s3 --profile saml cp output/pbta-gene-expression-rsem-tpm-collapsed.combined.rds s3://d3b-bix-dev-data-bucket/hgg-dmg-integration/merged_ngs_files/rna-seq
aws s3 --profile saml cp output/pbta-gene-expression-rsem-fpkm-collapsed.combined.rds s3://d3b-bix-dev-data-bucket/hgg-dmg-integration/merged_ngs_files/rna-seq
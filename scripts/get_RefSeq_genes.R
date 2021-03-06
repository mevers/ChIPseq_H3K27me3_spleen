library(biomaRt);

# Load Biomart
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "ensembl.org");
mart <- useDataset("mmusculus_gene_ensembl", mart = mart);

# Attributes
attributes <- listAttributes(mart);
chr <- c(1:19, "X", "Y", "MT");

# Get all genes
genes <- getBM(
    attributes = c("refseq_mrna", "refseq_ncrna", "ensembl_gene_id", "external_gene_name", "entrezgene"),
    filters = "chromosome_name",
    values = chr,
    mart = mart);
genes <- cbind(refseq_gene_id = "", genes);
genes$refseq_gene_id <- ifelse(nchar(genes$refseq_mrna) > 0, genes$refseq_mrna, genes$refseq_ncrna);


# Save genes
write.csv(genes[order(genes$refseq_mrna, genes$ensembl_gene_id), -c(2, 3)], "refseq_gene_IDs.csv", row.names = FALSE);

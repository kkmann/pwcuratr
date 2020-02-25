query_reactome_pathways <- function(genes) {
    valid_ensembl_gene_ids(genes)
    pwcuratr_tbls$ensembl2pathways %>%
        filter(ensembl_id %in% genes) %>%
        pull(reactome_pathway_id) %>%
        unique
}

query_participating_genes <- function(pathways) {
    valid_reactome_pathways(pathways)
    pwcuratr_tbls$ensembl2pathways %>%
        filter(reactome_pathway_id %in% pathways) %>%
        pull(ensembl_id) %>%
        unique
}

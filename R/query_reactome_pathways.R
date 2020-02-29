#' @export
query_reactome_pathways <- function(genes) {
    valid_ensembl_gene_ids(genes)
    tbl_ensembl2reactome %>%
        filter(.data$ensembl_id %in% genes) %>%
        pull(.data$reactome_pathway_id) %>%
        unique
}

#' @export
query_participating_genes <- function(pathways) {
    valid_reactome_pathways(pathways)
    tbl_ensembl2reactome %>%
        filter(.data$reactome_pathway_id %in% pathways) %>%
        pull(.data$ensembl_id) %>%
        unique
}

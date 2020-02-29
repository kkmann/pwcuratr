#' Query reactome for involved biological pathways
#'
#' Get all reactome.org pathways that involve any of the genes provided.
#'
#' @param genes character vector of valid ensembl gene IDs (no version suffix)
#'
#' @return character vector of reactome IDs of the pathways involving at least
#' one of the input genes.
#'
#' @examples
#' query_reactome_pathways(c("ENSG00000130203", "ENSG00000189058"))
#'
#' @include zzz.R
#'
#' @export
query_reactome_pathways <- function(genes) {
    valid_ensembl_gene_ids(genes)
    tbl_ensembl2reactome %>%
        filter(.data$ensembl_id %in% genes) %>%
        pull(.data$reactome_pathway_id) %>%
        unique
}



#' Query reactome for all participating genes in a set of pathways
#'
#' Get all genes participating in the provided list of reactome pathway identifiers.
#'
#' @param pathways character vector of valid reactome IDs
#'
#' @return character vector of ensembl IDs of the genes participating in any of the
#' input pathways.
#'
#' @examples
#' query_participating_genes("R-HSA-425407")
#'
#' @include zzz.R
#'
#' @export
query_participating_genes <- function(pathways) {
    valid_reactome_pathways(pathways)
    tbl_ensembl2reactome %>%
        filter(.data$reactome_pathway_id %in% pathways) %>%
        pull(.data$ensembl_id) %>%
        unique
}

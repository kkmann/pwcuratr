#' Convert gene list to igraph object
#'
#' Convert a character vector of ensembl gene IDs to an igraph object using the
#' pre-stored gene-gene interactions as undirected edges.
#' Interactions can be filtered by their evidence score.
#'
#' @param genes character vector of valid ensembl gene IDs (no version suffix)
#' @param minscore numeric between 0 and 1, minimal score for gene-gene
#' interactions to be used
#'
#' @return igraph object
#'
#' @examples
#' as_igraph(c("ENSG00000130203", "ENSG00000189058"),  0.9)
#'
#' @include zzz.R
#'
#' @export
as_igraph <- function(genes, minscore = 0) {
    genes <- unique(genes)
    tbl_interactions <- tbl_reactome_interactions %>%
        filter(
            .data$a %in% genes,
            .data$b %in% genes,
            .data$score >= minscore
        ) %>%
        select(.data$a, .data$b)
    graph(
        tbl_interactions %>%
            as.matrix() %>%
            t() %>%
            as.character(),
        isolates = genes[
            !(genes %in% tbl_interactions$a) & !(genes %in% tbl_interactions$b)
        ],
        directed = FALSE
    )
}

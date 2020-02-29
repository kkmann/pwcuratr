#' Pruning of a gene list by interaction distance to seed genes
#'
#' Prune a list of genes by maximal edge distance from provided seed genes using
#' gene-gene functional interactions as edges.
#'
#' @param genes character vector of ensembl IDs
#' @param seed_genes character vector of ensembl IDs
#' @param minscore minimal evidence score for edges to include (see [as_igraph()])
#' @param maxedgedistance only genes that are connected to any of the seed genes
#'     via `maxedgedistance` or less edges (interactions) are retained
#'
#' @return a character vector with retained genes, always containes all seed genes
#'
#' @examples
#' prune(c("ENSG00000130203", "ENSG00000189058"), c("ENSG00000130203"), 0.9, 2)
#'
#' @include zzz.R
#'
#' @export
prune <- function(genes, seed_genes, minscore = 0, maxedgedistance = Inf) {
    tbl_interactions <- tbl_reactome_interactions %>%
        filter(
            .data$a %in% genes,
            .data$b %in% genes,
            .data$score >= minscore
        ) %>%
        select(.data$a, .data$b)
    gr <- igraph::graph(
        tbl_interactions %>%
            as.matrix() %>%
            t() %>%
            as.character(),
        isolates = seed_genes[
            !(seed_genes %in% tbl_interactions$a) & !(seed_genes %in% tbl_interactions$b)
        ],
        directed = FALSE
    )
    dists <- igraph::distances(gr, igraph::V(gr), to = seed_genes)
    dists %>%
        as_tibble(rownames = "from") %>%
        pivot_longer(-.data$from, names_to = "to",  values_to = "dist") %>%
        filter(.data$to %in% seed_genes, .data$dist <= maxedgedistance) %>%
        pull(.data$from) %>%
        unique
}

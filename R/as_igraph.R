as_igraph <- function(genes, minscore = 0) {
    genes <- unique(genes)
    tbl_interactions <- pwcuratr_tbls$reactome_interactions %>%
        filter(
            a %in% genes & b %in% genes,
            score >= minscore
        ) %>%
        select(-id, -score)
    igraph::graph(
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

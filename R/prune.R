prune <- function(genes, seed_genes, minscore = 0, maxedgedistance = Inf) {
    tbl_interactions <- pwcuratr_tbls$reactome_interactions %>%
        filter(
            a %in% genes & b %in% genes,
            score >= minscore
        ) %>%
        select(-id, -score)
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
        pivot_longer(-from, names_to = "to",  values_to = "dist") %>%
        filter(to %in% seed_genes, dist <= maxedgedistance) %>%
        pull(from) %>%
        unique
}

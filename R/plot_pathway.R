plot_graph <- function(igraph, seed_genes, title = NULL, titlesize = 10) {
    ggr <- tidygraph::as_tbl_graph(igraph) %>%
        tidygraph::activate(nodes)
    ggr <- left_join(
        ggr,
        ggr %>%
            as_tibble() %>%
            left_join(
                pwcuratr_tbls$ensemble %>%
                    select(external_gene_name, ensembl_gene_id) %>%
                    distinct(),
                by = c(name = "ensembl_gene_id")
            ) %>%
            group_by(name) %>%
            summarise(
                external_name = paste(unique(external_gene_name), collapse = "|")
            ) %>%
            mutate(external_name = ifelse(external_name == "", name, external_name)),
        by = "name"
    )
    ggraph::ggraph(ggr) +
        ggraph::geom_edge_link(
            alpha = .05
        ) +
        ggraph::geom_node_point(
            aes(color = name %in% seed_genes)
        ) +
        ggraph::geom_node_label(
            aes(
                label = external_name,
                color = name %in% seed_genes
            ),
            repel         = TRUE,
            label.padding = unit(0.1, "lines"),
            label.r       = unit(0.1, "lines"),
            show.legend   = FALSE,
            segment.size  = .2,
            box.padding   = unit(0.15, "lines")
        ) +
        scale_color_manual(
            "seed gene",
            values = c(
                "TRUE"  = scales::muted("green", l = 50, c = 70),
                "FALSE" = "black"
            )
        ) +
        {if (!is.null(title)) {ggtitle(label = title)} else {NULL}} +
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text  = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none",
            plot.title = element_text(size = titlesize)
        )
}



plot_pathway <- function(genes, seed_genes, minscore = 0, titlesize = 10) {
    tbl_cmpnts <- as_igraph(genes, minscore = minscore) %>%
        igraph::decompose(min.vertices = 2) %>%
        enframe
    plts <- map2(
        tbl_cmpnts$name, tbl_cmpnts$value,
        ~plot_graph(
            ..2,
            seed_genes,
            title = sprintf("Component %i", ..1),
            titlesize = titlesize
        )
    )
    cowplot::plot_grid(plotlist = plts, ncol = 1)
}

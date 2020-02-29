#' Plot a gene network graph
#'
#' @param igraph an igraph object as returned by [as_igraph()]
#' @param seed_genes character vector of ensembl IDs of the seed genes in the igraph object
#' @param title plot title
#' @param titlesize text size for the title
#'
#' @return a ggplot object
#'
#' @include as_igraph.R
#'
#' @examples
#' gr <- as_igraph(c("ENSG00000130203", "ENSG00000189058"),  0.9)
#' plot_graph(gr, c("ENSG00000130203"))
#'
#' @rdname plot-gene-network
#' @aliases plot_graph
#'
#' @include zzz.R
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_label
#'
#' @export
plot_graph <- function(igraph, seed_genes, title = NULL, titlesize = 10) {

    ggr <- as_tbl_graph(igraph)
    ggr <- left_join(
        ggr,
        ggr %>%
            as_tibble() %>%
            left_join(
                tbl_ensembl %>%
                    select(.data$external_gene_name, .data$ensembl_gene_id) %>%
                    distinct(),
                by = c(name = "ensembl_gene_id")
            ) %>%
            group_by(.data$name) %>%
            summarise(
                external_name = paste(unique(.data$external_gene_name), collapse = "|")
            ) %>%
            mutate(
                external_name = ifelse(.data$external_name == "", .data$name, .data$external_name)
            ),
        by = "name"
    )
    ggraph(ggr, layout = "stress") +
        geom_edge_link(
            alpha = .05
        ) +
        geom_node_point(
            ggplot2::aes(color = .data$name %in% seed_genes)
        ) +
        geom_node_label(
            ggplot2::aes(
                label = .data$external_name,
                color = .data$name %in% seed_genes
            ),
            repel         = TRUE,
            label.padding = grid::unit(0.1, "lines"),
            label.r       = grid::unit(0.1, "lines"),
            show.legend   = FALSE,
            segment.size  = .2,
            box.padding   = grid::unit(0.15, "lines")
        ) +
        ggplot2::scale_color_manual(
            "seed gene",
            values = c(
                "TRUE"  = scales::muted("green", l = 50, c = 70),
                "FALSE" = "black"
            )
        ) +
        {if (!is.null(title)) {ggplot2::ggtitle(label = title)} else {NULL}} +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid      = ggplot2::element_blank(),
            axis.title      = ggplot2::element_blank(),
            axis.text       = ggplot2::element_blank(),
            axis.ticks      = ggplot2::element_blank(),
            legend.position = "none",
            plot.title      = ggplot2::element_text(size = titlesize)
        )
}


#' Plot a gene network graph
#'
#' @param genes character vector of ensembl IDs of the genes in the network
#' @param minscore minimal evidence score for edges to include (see [as_igraph()])
#' @param minnodes minimal number of nodes in component to include in plot
#'
#' @examples
#' plot_pathway(
#'     c("ENSG00000130203", "ENSG00000189058"),
#'     c("ENSG00000130203", "ENSG00000189058"),
#'     0.9
#' )
#'
#' @importFrom cowplot plot_grid
#' @importFrom purrr map2
#'
#' @rdname plot-gene-network
#' @aliases plot_pathway
#'
#' @export
plot_pathway <- function(genes, seed_genes, minscore = 0, titlesize = 10, minnodes = 1) {
    tbl_cmpnts <- as_igraph(genes, minscore = minscore) %>%
        igraph::decompose(min.vertices = minnodes) %>%
        enframe
    plts <- purrr::map2(
        tbl_cmpnts$name, tbl_cmpnts$value,
        ~plot_graph(
            ..2,
            seed_genes,
            title     = sprintf("Component %i", ..1),
            titlesize = titlesize
        )
    )
    plot_grid(plotlist = plts, ncol = 1)
}

#' Plot a gene network graph
#'
#' @param igraph an igraph object as returned by [as_igraph()]
#' @param seed_genes character vector of ENSEMBL IDs of the seed genes in the igraph object
#' @param title plot title
#' @param titlesize text size for the title
#' @param textsize text size for labels
#' @param pointsize point size
#' @param edgesize edge strength
#' @param layout ggraph layout
#' @param layout_args list of argments to pass to ggraph layout function
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
plot_graph <- function(igraph,
                       seed_genes, title = NULL, titlesize = 10,
                       textsize = 3, pointsize = 1, edgesize = .5,
                       layout   = "stress", layout_args = list()
) {
    tbl_components <- tibble::tibble(
            gr = igraph::decompose(igraph)
        ) %>%
        mutate(
            component = as.character(dplyr::row_number())
        ) %>%
        select(.data$component, dplyr::everything()) %>%
        mutate(
            name = purrr::map(.data$gr, ~igraph::get.vertex.attribute(., "name"))
        ) %>%
        select(-.data$gr) %>%
        tidyr::unnest(.data$name)

    ggr <- as_tbl_graph(igraph)
    ggr <- left_join(
            ggr,
            ggr %>%
                as_tibble() %>%
                mutate(
                    external_gene_name = get_external_names(.data$name)
                ) %>%
                group_by(.data$name) %>%
                summarise(
                    external_name = paste(unique(.data$external_gene_name), collapse = "|")
                ) %>%
                mutate(
                    external_name = ifelse(.data$external_name == "unknown", .data$name, .data$external_name)
                ),
            by = "name"
        ) %>%
        left_join(tbl_components, by = "name") %>%
        mutate(
            is_seed = .data$name %in% seed_genes
        )
    do.call(ggraph, args = c(list(ggr), list(layout = layout), layout_args)) +
        geom_edge_link(
            alpha = .05,
            edge_width = edgesize
        ) +
        geom_node_point(
            ggplot2::aes(
                color = .data$component
            ),
            size = pointsize
        ) +
        geom_node_label(
            ggplot2::aes(
                label = ifelse(!.data$is_seed, .data$external_name, ""),
                color = .data$component#.data$name %in% seed_genes
            ),
            repel         = TRUE,
            label.padding = grid::unit(textsize * 0.033, "lines"),
            label.r       = grid::unit(textsize * 0.033, "lines"),
            show.legend   = FALSE,
            segment.size  = textsize * 0.05,
            box.padding   = grid::unit(textsize * 0.033, "lines"),
            size          = textsize
        ) +
        geom_node_label(
            ggplot2::aes(
                label = ifelse(.data$is_seed, .data$external_name, "")
            ),
            repel         = TRUE,
            label.padding = grid::unit(textsize * 0.033, "lines"),
            label.r       = grid::unit(textsize * 0.033, "lines"),
            show.legend   = FALSE,
            segment.size  = textsize * 0.05,
            box.padding   = grid::unit(textsize * 0.033, "lines"),
            size          = textsize
        ) +
        ggplot2::scale_colour_brewer(type = "qual", palette = 1) +
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

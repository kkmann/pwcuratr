globalVariables(c("tbl_reactome_interactions", "tbl_ensembl", "tbl_ensembl2reactome"))

#' A package to facilitate curation of functional gene networks
#'
#' Todo
#'
#' @import tidyverse
#' @import DT
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble enframe
#' @importFrom dplyr filter select left_join group_by pull summarise mutate distinct
#' @importFrom tidyr pivot_longer
#' @importFrom igraph graph
#' @importFrom tidygraph as_tbl_graph activate
#' @importFrom glue glue
#'
#' @name pwcuratr
NULL



#' ENSEMBL ID to gene name
#'
#' This dataset is used to map ENSEMBL IDs to gene names.
#' It is parsed from the biomaRt bioconductor package using ensebl version 97.
#'
#' @format
#' \describe{
#'   \item{external_gene_name}{gene name}
#'   \item{ensembl_gene_id}{ENSEMBL gene ID}
#' }
#'
#' @name tbl_ensembl
#'
#' @source \url{https://bioconductor.org/packages/release/bioc/html/biomaRt.html}
"tbl_ensembl"



#' ENSEMBL ID to reactome.org pathway
#'
#' This dataset is used to map ENSEMBL IDs to reactome.org pathway IDs.
#' It is parsed from a downloaded file of an identfier mapping from
#' https://reactome.org/download.
#'
#' @format
#' \describe{
#'   \item{ensembl_id}{ENSEMBL gene ID}
#'   \item{reactome_id}{reactome.org ID for gene}
#'   \item{reactome_pathway_id}{reactome.org ID for pathway}
#'   \item{description}{reactome.org description of pathway}
#' }
#'
#' @name tbl_ensembl2reactome
#'
#' @source \url{https://reactome.org/download/current/Ensembl2Reactome.txt}
"tbl_ensembl2reactome"



#' Function gene-gene interactions
#'
#' This dataset contains the predicted functional gene-gene interactions.
#' The methodology is described in
#'
#' Wu, G., Feng, X., & Stein, L. (2010). A human functional protein interaction
#' network and its application to cancer data analysis. Genome biology, 11(5), R53.
#'
#' @format
#' \describe{
#'   \item{score}{evidence score}
#'   \item{a}{ENSEMBL gene ID of interactor a}
#'   \item{b}{ENSEMBL gene ID of interactor b}
#' }
#'
#' @name tbl_reactome_interactions
#'
#' @source \url{http://cpws.reactome.org/caBigR3WebApp2018/FIsInGene_122718_with_annotations.txt.zip}
"tbl_reactome_interactions"

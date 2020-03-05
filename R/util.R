valid_ensembl_gene_ids <- function(genes) {
    problematic_genes <- genes[!(genes %in% tbl_ensembl$ensembl_gene_id)]
    if (length(problematic_genes) == 0)
        return(TRUE)
    warning(
        glue("\rensembl gene IDs not found:\n\r{paste(problematic_genes, collapse = '\n\r')}")
    )
}

valid_reactome_pathways <- function(pathways) {
    problematic_pathways <- pathways[!(pathways %in% tbl_ensembl2reactome$reactome_pathway_id)]
    if (length(problematic_pathways) == 0)
        return(TRUE)
    warning(
        glue("\rreactome pathway IDs not found:\n\r{paste(problematic_pathways, collapse = '\n\r')}")
    )
}



find_isolates <- function(genes, minscore = 0) {
    tbl_interactions <- tbl_reactome_interactions %>%
        filter(
            .data$score >= minscore
        ) %>%
        select(.data$a, .data$b)
    unique(genes[
        !(genes %in% tbl_interactions$a) & !(genes %in%  tbl_interactions$b)
    ])
}



#' Convert ENSEMBL IDs to external gene names
#'
#' @param genes character vector of valid ENEMBL gene IDs (no version suffix)
#'
#' @return character vector of matching external gene names
#'
#' @examples
#' get_external_names("ENSG00000278801")
#'
#' @include zzz.R
#'
#' @export
get_external_names <- function(genes) {
    res <- tibble::tibble(
        ensembl_gene_id = genes
    ) %>%
    left_join(
        tbl_ensembl,
        by = "ensembl_gene_id"
    ) %>%
    mutate(
        external_gene_name = ifelse(
            is.na(.data$external_gene_name), "unknown", .data$external_gene_name
        )
    ) %>%
    pull(
        .data$external_gene_name
    )
    if (length(res) != length(genes)) stop("mismatch")
    return(res)
}



# parse_reactome_functional_interactions <- function(
#     file = "inst/data/FIsInGene_122718_with_annotations.txt",
#     outfile = "inst/data/tbl_reactome_interactions.rds"
# ) {
#     vroom::vroom(
#         file,
#         delim = "\t",
#         col_names = c(
#             "Gene1",
#             "Gene2",
#             "Annotation",
#             "Direction",
#             "Score"
#         ),
#         skip = 1
#     ) %>%
#     transmute(
#         a     = Gene1,
#         b     = Gene2,
#         score = Score
#     ) %>%
#     left_join(
#         pwcuratr_tbls$ensemble %>%
#             transmute(
#                 external_gene_name,
#                 a_ensembl_gene_id = ensembl_gene_id
#             ),
#         by = c(a = "external_gene_name")
#     ) %>%
#     left_join(
#         pwcuratr_tbls$ensemble %>%
#             transmute(
#                 external_gene_name,
#                 b_ensembl_gene_id = ensembl_gene_id
#             ),
#         by = c(b = "external_gene_name")
#     ) %>%
#     filter(complete.cases(.)) %>%
#     select(-a, -b) %>%
#     rename(
#         a = a_ensembl_gene_id,
#         b = b_ensembl_gene_id
#     ) %>%
#     distinct(a, b, .keep_all = TRUE) %>%
#     mutate(
#         id = row_number()
#     ) %>%
#     pivot_longer(
#         c(a, b),
#         names_to = "interactor",
#         values_to = "ensembl_id"
#     ) %>%
#     group_by(id) %>%
#     arrange(id, ensembl_id) %>%
#     mutate(
#         interactor = c("a", "b")
#     ) %>%
#     pivot_wider(
#         names_from = "interactor",
#         values_from = "ensembl_id"
#     ) %>%
#     select(-id) %>%
#     ungroup() %>%
#     distinct() %>%
#     readr::write_rds(outfile, compress = "gz")
# }
#
#
#
# parse_ensembl2reactome <- function(
#     file = "inst/data/Ensembl2Reactome_PE_All_Levels.txt.gz",
#     outfile = "inst/data/tbl_ensembl2reactome.rds"
# ) {
#     vroom::vroom(
#         file,
#         delim = "\t",
#         col_names = c(
#             "ensembl_id",
#             "reactome_id",
#             "description",
#             "reactome_pathway_id",
#             "url",
#             "type",
#             "evidence_level",
#             "species"
#         )
#     ) %>%
#     filter(
#         species == "Homo sapiens"
#     ) %>%
#     select(ensembl_id, reactome_id, reactome_pathway_id, evidence_level) %>%
#     distinct() %>%
#     filter(complete.cases(.)) %>%
#     mutate( # remove version suffix (not used consistently)
#         ensembl_id = str_extract(ensembl_id, "[[:alnum:]]+(?=([\\.[0-9]{1,3}]$)*)")
#     ) %>%
#     distinct() %>%
#     write_rds(outfile, compress = 'gz')
# }

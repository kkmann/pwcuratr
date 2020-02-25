.onLoad <- function(libname, pkgname) {
        pwcuratr_tbls <<- list()

        pwcuratr_tbls$ensembl2pathways <<- readr::read_rds(
            glue::glue("{find.package('pwcuratr')}/data/tbl_ensembl2reactome.rds")
        )

        pwcuratr_tbls$ensemble <<- readr::read_rds(
            glue::glue("{find.package('pwcuratr')}/data/tbl_ensemble_97.rds")
        )

        pwcuratr_tbls$reactome_interactions <<- readr::read_rds(
            glue::glue("{find.package('pwcuratr')}/data/tbl_reactome_interactions.rds")
        )
}

# tmp <- pwcuratr_tbls$ensembl2pathways %>%
#         select(reactome_pathway_id) %>%
#         distinct() %>%
#         mutate(
#         description = map_chr(
#                         reactome_pathway_id,
#                         ~httr::GET(glue("https://reactome.org/ContentService/data/query/{.}/displayName")) %>%
#                                 httr::parsed_content()
#                 )
#         )
#
# pwcuratr_tbls$ensembl2pathways %>%
#         left_join(tmp) %>%
#         write_rds("inst/data/tbl_ensembl2reactome.rds", compress = "gz")


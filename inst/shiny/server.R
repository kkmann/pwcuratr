library(shiny)
library(glue)
library(pwcuratr)
library(tidyverse)



server <- function(input, output) {

    tbls <- reactiveValues()

    tbls$interactions_gene_gene <- reactive({
        pwcuratr_tbls$reactome_interactions %>%
            filter(score >= input$minScore)
    })

    tbls$seed_genes <- tibble(
        external_gene_name = character(0L),
        ensembl_gene_id    = character(0L)
    )
    output$tblSeedGenes <- DT::renderDataTable(
        tbls$seed_genes %>%
            mutate(
                ensembl_gene_id = map_chr(
                    ensembl_gene_id,
                    ~glue('<a href=https://www.ensembl.org/id/{.} target="_blank">{.}</a>')
                )
            ) %>%
            rename(
                ENSEMBl = ensembl_gene_id,
                Name    = external_gene_name
            ),
        options = list(
            paging = FALSE
        ),
        escape   = FALSE,
        rownames = FALSE
    )
    tblSeedGenesProxy <- DT::dataTableProxy("tblSeedGenes")
    observeEvent(input$uploadSeedGenesFile, {
        seed_genes_selected_old <- tbls$seed_genes %>%
            filter(
                row_number() %in% input$tblSeedGenes_rows_selected
            ) %>%
            pull(ensembl_gene_id)
        new_data <- read_csv(
                input$uploadSeedGenesFile$datapath[1],
                col_types = cols_only(
                    ensembl_gene_id = col_character()
                )
            ) %>%
            filter(
                !(ensembl_gene_id %in% tbls$seed_genes$ensembl_gene_id)
            ) %>%
            left_join(
                pwcuratr_tbls$ensemble %>%
                    dplyr::select(external_gene_name, ensembl_gene_id) %>%
                    dplyr::filter(complete.cases(.)),
                by = "ensembl_gene_id"
            ) %>%
            dplyr::distinct(external_gene_name, ensembl_gene_id)
        seed_genes_selected_new <- new_data$ensembl_gene_id
        tbls$seed_genes <- bind_rows(
                tbls$seed_genes,
                new_data
            ) %>%
            dplyr::distinct(external_gene_name, ensembl_gene_id) %>%
            dplyr::arrange(external_gene_name)
        seed_genes_selected <- c(
            seed_genes_selected_old,
            seed_genes_selected_new
        )
        tblSeedGenesProxy %>%
            DT::selectRows(
                which(tbls$seed_genes$ensembl_gene_id %in% seed_genes_selected)
            )
    })
    observeEvent(input$addGene, {
        seed_genes_selected_old <- tbls$seed_genes %>%
            filter(
                row_number() %in% input$tblSeedGenes_rows_selected
            ) %>%
            pull(ensembl_gene_id)
        new_data <- tibble(
                external_gene_name = str_extract(input$selectGene, "^.+(?= \\()"),
                ensembl_gene_id    = str_extract(input$selectGene, "(?<=^.{1,99} \\()[[:alnum:]]+(?=\\)$)")
            )
        seed_genes_selected_new <- new_data$ensembl_gene_id
        tbls$seed_genes <- bind_rows(
                tbls$seed_genes,
                new_data
            ) %>%
            dplyr::distinct(external_gene_name, ensembl_gene_id, .keep_all = TRUE) %>%
            dplyr::arrange(external_gene_name)
        seed_genes_selected <- c(
            seed_genes_selected_old,
            seed_genes_selected_new
        )
        tblSeedGenesProxy %>%
            DT::selectRows(
                which(tbls$seed_genes$ensembl_gene_id %in% seed_genes_selected)
            )
    })
    observeEvent(input$uploadPathwayZip, {
        tdir <- tempdir()
        odir <- setwd(tdir)
        unzip(input$uploadPathwayZip$datapath[1])
        tbls$seed_genes <- read_csv(
            "seed_genes.csv",
            col_types = cols(
                external_gene_name = col_character(),
                ensembl_gene_id    = col_character()
            )
        )
        setwd(odir)
        tblSeedGenesProxy %>%
            DT::selectRows(1:nrow(tbls$seed_genes))
    })
    seed_genes <- reactive({
        tbls$seed_genes %>%
        filter(
            row_number() %in% input$tblSeedGenes_rows_selected
        ) %>%
        pull(ensembl_gene_id)
    })
    output$tbl_selectedSeedGenes <- renderTable({
        tbls$seed_genes %>%
            filter(
                ensembl_gene_id %in% seed_genes()
            ) %>%
            arrange(external_gene_name) %>%
            transmute(
                name = sprintf("%s (%s)", external_gene_name, ensembl_gene_id)
            )
        },
        rownames = FALSE,
        colnames = FALSE
    )

    tbls$pathways <- tibble(
        pathway     = character(0L),
        description = character(0L),
        seed_genes  = list(),
        n_genes     = integer(0L)
    )
    output$tbl_candidate_pathways <- DT::renderDataTable({
        tbls$pathways %>%
            mutate(
                pathway = map_chr(
                    pathway,
                    ~glue('<a href=https://reactome.org/content/detail/{.} target="_blank">{.}</a>')
                ),
                seed_genes = map_chr(
                    seed_genes,
                    function(seed_genes) {
                        pwcuratr_tbls$ensemble %>%
                            select(ensembl_gene_id, external_gene_name) %>%
                            filter(ensembl_gene_id %in% seed_genes) %>%
                            distinct() %>%
                            pull(external_gene_name) %>%
                            paste(collapse = ", ")
                    }
                )
            ) %>%
            rename(
                `Reactome ID` = pathway,
                Description   = description,
                `contained seed genes` = seed_genes,
                `total number of genes` = n_genes
            )
        },
        options = list(
            paging = FALSE
        ),
        escape   = FALSE,
        rownames = FALSE
    )
    tblCandidatePathwaysProxy <- DT::dataTableProxy("tbl_candidate_pathways")
    update_pathways <- function() {
        if (is.null(input$uploadPathwayZip$datapath[1])) {
            pathways_selected_old_cache <- tbls$pathways %>%
                filter(
                    row_number() %in% input$tbl_candidate_pathways_rows_selected
                ) %>%
                pull(pathway)
        } else {
            tdir <- tempdir()
            odir <- setwd(tdir)
            unzip(input$uploadPathwayZip$datapath[1])
            pathways_selected_old_cache <- read_csv(
                "pathways.csv",
                col_types = cols(
                    pathway     = col_character(),
                    description = col_character()
                )
            ) %>%
                pull(pathway)
            setwd(odir)
        }
        tbls$pathways <- tibble(
                pathway = query_reactome_pathways(seed_genes())
            ) %>%
            left_join(
                pwcuratr_tbls$ensembl2pathways %>%
                    select(reactome_pathway_id, description) %>%
                    distinct(),
                by = c(pathway = "reactome_pathway_id")
            ) %>%
            mutate(
                seed_genes = map(
                    pathway,
                    ~seed_genes()[seed_genes() %in% query_participating_genes(.)]
                ),
                n_genes = map_int(
                    pathway,
                    ~filter(
                        pwcuratr_tbls$ensembl2pathways,
                        reactome_pathway_id == .,
                    ) %>%
                    select(reactome_pathway_id, ensembl_id) %>%
                    distinct() %>%
                    nrow()
                )
            ) %>%
            arrange(n_genes)

        tblCandidatePathwaysProxy %>%
            DT::selectRows(which(tbls$pathways$pathway %in% pathways_selected_old_cache))
    }
    observeEvent(
        input$btn_updatePathwaySelection, {
        update_pathways()
    })

    output$tbl_candidate_pathway_summary <- renderTable({
        if (any(input$tbl_candidate_pathways_rows_selected)) {
            tbls$pathways %>%
                filter(
                    row_number() %in% input$tbl_candidate_pathways_rows_selected
                ) %>%
                summarise(
                    `#pathways` = as.character(n()),
                    `#genes`    = pwcuratr_tbls$ensembl2pathways %>%
                        filter(reactome_pathway_id %in% pathway) %>%
                        pull(ensembl_id) %>%
                        unique %>%
                        length() %>%
                        as.character(),
                    `seed genes covered` = pwcuratr_tbls$ensemble %>%
                        select(ensembl_gene_id, external_gene_name) %>%
                        filter(ensembl_gene_id %in% unlist(seed_genes)) %>%
                        distinct() %>%
                        pull(external_gene_name) %>%
                        unique() %>%
                        sort() %>%
                        paste(collapse = ", "),
                    `seed genes not covered` = pwcuratr_tbls$ensemble %>%
                        select(ensembl_gene_id, external_gene_name) %>%
                        filter(
                            ensembl_gene_id %in% seed_genes(),
                            !(ensembl_gene_id %in% unlist(seed_genes))
                        ) %>%
                        distinct() %>%
                        pull(external_gene_name) %>%
                        unique() %>%
                        sort() %>%
                        paste(collapse = ", ")
                ) %>%
                pivot_longer(everything())
        } else {
            tibble(
                `#pathways` = "0",
                `#genes`    = "0",
                `seed genes selected` = "",
                `seed genes not selected` = pwcuratr_tbls$ensemble %>%
                    select(ensembl_gene_id, external_gene_name) %>%
                    filter(
                        ensembl_gene_id %in% unlist(tbls$seed_genes$ensembl_gene_id)
                    ) %>%
                    distinct() %>%
                    pull(external_gene_name) %>%
                    unique() %>%
                    sort() %>%
                    paste(collapse = ", ")
            ) %>%
            pivot_longer(everything())
        }
    }, colnames = FALSE)

    candidate_genes <- function(minscore, maxedgedistance) {
        tbls$pathways %>%
            filter(
                row_number() %in% input$tbl_candidate_pathways_rows_selected
            ) %>%
            pull(pathway) %>%
            query_participating_genes %>%
            prune(
                seed_genes      = tbls$seed_genes$ensembl_gene_id,
                maxedgedistance = maxedgedistance,
                minscore        = minscore
            )
    }

    output$plt_select_k <- renderPlot({
        if (length(input$tbl_candidate_pathways_rows_selected) == 0) return(NULL)
        p1 <- tibble(
            maxedgedistance = 0:10,
            n_genes_pruned  = map_int(
                maxedgedistance,
                ~length(candidate_genes(
                    minscore        = input$minScore,
                    maxedgedistance = .
                ))
            )
        ) %>%
        ggplot() +
            aes(maxedgedistance, n_genes_pruned) +
            geom_bar(
                aes(fill = maxedgedistance == input$pruning_distance),
                stat = "identity",
                show.legend = FALSE
            ) +
            scale_y_continuous(
                "",
                limits = c(0L, NA_integer_)
            ) +
            scale_x_continuous("maximal edge distance", breaks = 0:10) +
            scale_color_manual("",
                               values = c(
                                   "TRUE"  = scales::muted("green", l = 60, c = 60),
                                   "FALSE" = "gray"
                               ),
                               aesthetics = "fill"
            ) +
            ggtitle("Number of retained candidate genes") +
            theme_bw() +
            theme(
                panel.grid.minor = element_blank()
            )
        p2 <- tibble(
            maxedgedistance = 0:10,
            n_components_pruned  = map_int(
                maxedgedistance,
                ~candidate_genes(
                    minscore        = input$minScore,
                    maxedgedistance = .
                ) %>%
                as_igraph(minscore = input$minScore) %>%
                igraph::components() %>%
                .[["no"]]
            )
        ) %>%
        ggplot() +
            aes(maxedgedistance, n_components_pruned) +
            geom_bar(
                aes(fill = maxedgedistance == input$pruning_distance),
                stat = "identity",
                show.legend = FALSE
            ) +
            scale_y_continuous(
                "",
                limits = c(0L, NA_integer_)
            ) +
            scale_x_continuous("maximal edge distance", breaks = 0:10) +
            scale_color_manual("",
                values = c(
                   "TRUE"  = scales::muted("green", l = 60, c = 60),
                   "FALSE" = "gray"
                ),
                aesthetics = "fill"
            ) +
            ggtitle("Number of connected components") +
            theme_bw() +
            theme(
                panel.grid.minor = element_blank()
            )
        cowplot::plot_grid(p1, p2, nrow = 1)
    }, height = 150)

    output$component_selector <- renderUI({
        if (length(seed_genes()) == 0) return(NULL)
        if (tbls$pathways %>%
            filter(
                row_number() %in% input$tbl_candidate_pathways_rows_selected
            ) %>%
            pull(pathway) %>%
            length() == 0) {
            cdt_genes <- seed_genes()
        } else {
            cdt_genes <- candidate_genes(
                minscore        = input$minScore,
                maxedgedistance = input$pruning_distance
            )
        }
        n_components <- cdt_genes %>%
            as_igraph() %>%
            igraph::components() %>%
            .[["no"]]
        selectInput("component_selector",
            label    = "graph component",
            choices  = as.character(1:n_components),
            multiple = FALSE,
            width    = "100%",
            selected = "1"
        )
    })
    output$pruning_distance <- renderUI({
        if (is.null(input$uploadPathwayZip$datapath[1])) {
            return(numericInput('pruning_distance',
                label = 'k',
                min = 0, max = NA, value = 1
            ))
        } else {
            tdir <- tempdir()
            odir <- setwd(tdir)
            unzip(input$uploadPathwayZip$datapath[1])
            pruning_distance <- jsonlite::read_json(
                    "parameters.json"
                )$pruning_distance
            setwd(odir)
            return(numericInput('pruning_distance',
                label = 'k',
                min = 0, max = NA, value = pruning_distance
            ))
        }
    })
    output$minscore <- renderUI({
        if (is.null(input$uploadPathwayZip$datapath[1])) {
            return(numericInput('minScore',
                label = 'min. score',
                min = 0, max = 1, step = .01, value = .99
            ))
        } else {
            tdir <- tempdir()
            odir <- setwd(tdir)
            unzip(input$uploadPathwayZip$datapath[1])
            minscore <- jsonlite::read_json(
                "parameters.json"
            )$minscore
            setwd(odir)
            return(numericInput('minScore',
                label = 'min. score',
                min = 0, max = 1, step = .01, value = minscore
            ))
        }
    })
    output$pathwayName <- renderUI({
        if (is.null(input$uploadPathwayZip$datapath[1])) {
            return(textInput("pwcName",
                             "pathway name",
                             value = "my_pathway_cluster",
                             placeholder = "my_pathway_cluster"))
        } else {
            tdir <- tempdir()
            odir <- setwd(tdir)
            unzip(input$uploadPathwayZip$datapath[1])
            name <- jsonlite::read_json(
                "parameters.json"
            )$name
            setwd(odir)
            return(textInput("pwcName", "pathway name", value = name,
                             placeholder = "my_pathway_name"))
        }
    })

    output$plt_pathway <- renderPlot({
        if (is.null(input$component_selector)) return(NULL)
        withProgress({
                candidate_genes(
                    minscore        = input$minScore,
                    maxedgedistance = input$pruning_distance
                ) %>%
                as_igraph() %>%
                igraph::decompose() %>%
                .[[as.integer(input$component_selector)]] %>%
                plot_graph(
                    seed_genes = tbls$seed_genes$ensembl_gene_id
                )
            }, value = 0, message = "plotting...")
        },
        height = function() 72*input$plt_height,
        width  = function() if (is.na(input$plt_width)) "auto" else 72*input$plt_width,
        res    = 72
    )

    output$download <- downloadHandler(
        filename = function() glue("{input$pwcName}.zip"),
        content = function(file) {
            tdir <- dirname(file)
            tbls$seed_genes %>%
                filter(
                    ensembl_gene_id %in% seed_genes()
                ) %>%
                write_csv(path = glue("{tdir}/seed_genes.csv"))
            tbls$pathways %>%
                select(
                    pathway, description
                ) %>%
                filter(
                    row_number() %in% input$tbl_candidate_pathways_rows_selected
                ) %>%
                write_csv(path = glue("{tdir}/pathways.csv"))
            candidate_genes(
                    minscore        = input$minScore,
                    maxedgedistance = input$pruning_distance
                ) %>%
                as_igraph() %>%
                write_rds(path = glue("{tdir}/igraph.rds"), compress = "gz")
            list(
                minscore = input$minScore,
                pruning_distance = input$pruning_distance,
                name = input$pwcName
            ) %>%
            jsonlite::write_json(glue("{tdir}/parameters.json"))
            pth <- setwd(tdir)
            zip(
                "tmp.zip",
                files = c(
                    "seed_genes.csv",
                    "pathways.csv",
                    "igraph.rds",
                    "parameters.json"
                ),
                flags = "-9XjD"
            )
            setwd(pth)
            file.rename(glue("{tdir}/tmp.zip"), file)
        }
    )

}

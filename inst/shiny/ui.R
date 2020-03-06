library(shiny)
library(pwcuratr)
library(tidyverse)

ui <- navbarPage("Pathway Cluster Curation Tool", id = "navbarpage",
    tabPanel("Info", value = "info_panel",
        withMathJax(includeMarkdown('markdown/info.md'))
    ),
    tabPanel("Define Seed Genes", value = "seed_gene_panel",
        sidebarLayout(
            sidebarPanel(
                fileInput('uploadPathwayZip', 'restore Pathway Cluster',
                    accept = c('.zip')
                ),
                hr(),
                fileInput('uploadSeedGenesFile', 'add seed genes from file',
                    accept = c('.csv')
                ),
                hr(),
                selectizeInput("selectGene",
                    label  = "add gene directly",
                    choices = tbl_ensembl %>%
                        transmute(
                            gene = sprintf("%s (%s)", external_gene_name, ensembl_gene_id)
                        ) %>%
                        pull(gene) %>%
                        unique()
                ),
                actionButton("addGene", label = "add gene", width = "100%"),
                hr(),
                tags$b("selected seed genes"),
                tableOutput("tbl_selectedSeedGenes")
            ),
            mainPanel(
                includeMarkdown('markdown/seed_genes_header.md'),
                br(), br(),
                DT::dataTableOutput('tblSeedGenes')
            )
        )
    ),
    tabPanel("Select Pathways", value = "select_pathway_panel",
        sidebarLayout(
            sidebarPanel(
                actionButton("btn_updatePathwaySelection",
                    label = "update from seed gene selection",
                    width = "100%"
                ),
                hr(),
                h4('current selection'),
                tableOutput('tbl_candidate_pathway_summary')
            ),
            mainPanel(
                includeMarkdown('markdown/select_pathways_header.md'),
                br(), br(),
                DT::dataTableOutput('tbl_candidate_pathways')
            )
        )
    ),
    tabPanel("Prune & Plot",
        sidebarLayout(
            sidebarPanel(
                uiOutput("pathwayName"),
                hr(),
                fluidRow(
                    column(6, uiOutput('pruning_distance')),
                    column(6, uiOutput('minscore'))
                ),
                hr(),
                uiOutput("component_selector"),
                fluidRow(
                    column(6, numericInput('plt_height', label = 'height (in)',
                                           min = 4, max = 50, value = 7
                    )),
                    column(6, numericInput('plt_width', label = 'width (in)',
                                           min = 4, max = 50, value = NA_integer_
                    )),
                ),
                fluidRow(
                    column(6, numericInput('plt_textsize', label = 'text-size',
                                           min = .5, max = 5, value = 4, step = .1
                    )),
                    column(6, numericInput('plt_pointsize', label = 'point-size',
                                           min = .5, max = 2, value = 2.5, step = .1
                    ))
                ),
                fluidRow(
                    column(6, numericInput('plt_edgesize',
                        label = 'edge-size',
                        min = .1,
                        max = 1,
                        value = .5,
                        step = .1
                    )),
                    column(6, numericInput('plt_bbox',
                        label = 'bounding box',
                        min = 1,
                        max = 50,
                        value = 3,
                        step = 1
                    ))
                ),
                fluidRow(
                    column(6, actionButton("plot", "plot", width = "100%")),
                    column(6, downloadButton("downloadPlot", label = "download"))
                ),
                hr(),
                downloadButton("download",
                    label = "download pathway cluster"
                )
            ),
            mainPanel(
                includeMarkdown('markdown/prune_header.md'),
                verticalLayout(
                    plotOutput("plt_select_k", height = 150),
                    plotOutput('plt_pathway')
                )
            )
        )
    )
)

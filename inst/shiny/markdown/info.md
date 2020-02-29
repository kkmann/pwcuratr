# Candidate Gene Pathway Clusters

**tl;dr:** elicit candidate gene networks based on seed-genes.

This tool allows the defintion of candidate gene networks based on a set of 
seed genes, a combination of [reactome](https://reactome.org/) pathway
information and predicted functional gene-gene interactions following

> Wu, G., Feng, X., & Stein, L. (2010). A human functional protein interaction network and its application to cancer data analysis. Genome biology, 11(5), R53.

This version uses the 2018 release of the function interaction dataset directly
downloaded from [https://reactome.org/download-data](https://reactome.org/download-data).
Internally all mappings are done via an [ENSEMBL](https://www.ensembl.org/index.html)
version 97 snapshot downloaded using the [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) package.

The app is part of the R package [pwcuratr](https://github.com/kkmann/pwcuratr).


## Define Seed-Genes

The first tab allows three different options for specifying the set of 
seed genes.

1. If a curated pathway cluster already exists and has been downloaded from the
app before, the .zip file can be uploaded and all settings are restored to 
continue working.
2. A set of genes can be uploaded directly as .csv file. The file must contain a
column 'ensembl_gene_id' with ENSEMBL gene (or transcript) identifiers.
Version suffixes are not supported since they are not consitently used in the
underlying databases.
All genes will be added to the current pool of potential seed genes and 
automatically be selected.
3. Individual genes can be added to the current pool of seed genes via the 
dropdown selector. 
The selector support fuzzy search over gene names and ENSEMBL gene identifiers when text is entered directly.


## Select Pathways

This tab allows to match the selected set of seed genes to reactome.org
biological pathway identifiers. 
The process must be started manually by clicking on the button 
'update from seed gene selection'.
All reactome pathways in which at least one of the seed genes participates 
will be listed with clickable links to the respective reactome.org page.
After revising the set of seed genes, the button must be clicked again but 
all currently selected pathways are retained if they contain at least one of
the revised seed genes.
If pathway cluster was restored by uploading a saved .zip file, 
the previously selected pathways will be restored when clicking 
'update from seed gene selection'.

The pathways are sorted by number of participating genes (ascending) to
prioritize specific (small) pathways.


## Pruning & Plotting

To further increase the specificity of the functional interaction neighborhood
defined by the pathway cluster, 
the resulting set of candidate genes (all genes contained in the selected
reactome pathways plus the original seed genes)
can be pruned.
Top this end, we use a recent (2018) release of a set of predicted functional 
gene-gene interactions using methodology for Wu et al, 2010.
Every interaction is scored between 0 (low confidence) and 1 (high confidence).
The set of interactions used to define a functional neighborhood is controlled by
specifying the 'minimal score' input.
Finally, the parameter 
$k$ 
controls the maximal number of edges between any
of the seed genes and a candidate gene.
I.e., all candidate genes whose closest connect to any of the seed genes is greater
than 
$k$ 
will be pruned from the set.

To guide the selection of the two pruning parameters, 
the number of retained genes and the number of connected components in the
final graph are plotted against 
$k$ 
with the current selection highlighted in 
green.

The individual connected components of the graph can be selected for plotting
via a dropdown menu on the left and the dimensions of the plot can be tailored 
manually to the number of genes. 

To save the current curated pathway cluster, simply click the download button.
The .zip file can then be used to restore the same pathway cluster by uploading
the file on the tab 'Define Seed Genes'.

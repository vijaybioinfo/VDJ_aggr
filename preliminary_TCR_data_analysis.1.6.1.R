############    -------------   TCR Data    -------------    ############
############    --------   Preliminary analyses   -------    ############
# By Vicente Fajardo

cat('############    -------------   TCR Data    -------------    ############\n')
cat('############    --------   Preliminary analyses   -------    ############\n')
cat('By Vicente Fajardo\n\n')
cat('Version: 1.6.0\n')

# Version: 1
# Subversion: 6
# Sub-subversion: 1
# Updates.
# ---> Version updates:
#   Same stable version.
# ---> Subversion updates:
#   * Removed various plots depicting clonotype sharing non-appropriatetly.
#   * Various Venn diagrams included to show clonotype sharing.
# ---> Sub-subversion updates:
#   * One further section was added to help correct the cell barcode suffixes in the table listing the relationships between the cell barcodes and the clonotypes. This process is solely based on the new correspondences that should be provided as a csv file (argument 'sffxs.corrs.file').


cat('\n\n')
### --------------------------- Libraries --------------------------- ###
cat('### --------------------------- Libraries --------------------------- ###\n')
cat('Importing libraries...\n')
library(Seurat)
library(optparse)
library(ggplot2)
# library(reshape2)
library(data.table)
library(english)
library(VennDiagram)
library(stringr)
library(gtools)
library(tidyr)
source('/home/vfajardo/scripts/functions/R_handy_functions.R')
source('/home/vfajardo/scripts/functions/R_visualization_functions.1.5.R')
cat('Libraries imported!\n')


cat('\n\n')
### --------------------------- Arguments --------------------------- ###
cat('### --------------------------- Arguments --------------------------- ###\n')
# Declaring arguments to parse from command line ----------------------->
option.list <- list(
  make_option(opt_str="--TCRContigs", type="character", default=NULL, dest="cells.clons.info.file", help="Absolute path to 'filtered_contig annotations.csv' file output from 10X."),
  make_option(opt_str="--TCRClonotypes", type="character", default=NULL, dest="clons.info.file", help="Absolute path to file describing clonotypes info and output from 10X (clonotypes.csv file)."),
  make_option(opt_str="--ReportsPath", type="character", default=NULL, dest="reports.path", help="Absolute path to output file."),
  make_option(opt_str="--SeuratObj", type="character", default=NULL, dest="seurat.obj.file", help="Absolute path to RDS seurat object file."),
  make_option(opt_str="--SffxsCorrs", type="character", default=NULL, dest="sffxs.corrs.file", help="Character, absolute path to csv file listing the cell barcode suffix correspondences between the seurat object and the table listing the relationships between the cell barcode and the clonotypes. All suffixes must exist in either element of the list."),
  make_option(opt_str="--ExpThold", type="integer", default=2, dest="expansion.thold", help="Threshold to account for clonal expansion."),
  make_option(opt_str="--AmountClons", type="integer", default=5, dest="amount.clons", help="Clonotype size to consider a clonotype whether as highly expanded or not\n."),
  make_option(opt_str="--Tags", type="character", default=NULL, dest="tags.of.int", help="Possible tags of interest. If input, they should already be part of the seurat object's meta data."),
  make_option(opt_str="--OutputObj", type="logical", default=TRUE, dest="output.obj", help="Logical, indicates whether the modified seurat object should be output or not."),
  make_option(opt_str="--PrjSffx", type="character", default=NULL, dest="prj.sffx", help="Tag suffix to add to the modified seurat object to be output. If not provided, the date the program is run will be added.")
)

# Getting arguments from command line and setting their values to their respective variable names.
opt.parser <- OptionParser(option_list=option.list)
opt <- parse_args(opt.parser)

# Moving options to their own variables
cells.clons.info.file <- opt$cells.clons.info.file
clons.info.file <- opt$clons.info.file
reports.path <- opt$reports.path
seurat.obj.file <- opt$seurat.obj.file
expansion.thold <- opt$expansion.thold
amount.clons <- opt$amount.clons
tags.of.int <- opt$tags.of.int
output.obj <- opt$output.obj
prj.sffx <- opt$prj.sffx
if(!is.null(tags.of.int)) tags.of.int <- eval(expr=parse(text=tags.of.int))
sffxs.corrs.file <- opt$sffxs.corrs.file
this.color.scale <- c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000') # provided by Ciro.
# Presenting parameters ------------------------------------------------>
cat('Absolute path to TCR contigs file (describing relations between TCR clonotypes and barcodes):', cells.clons.info.file, '\n')
cat('Absolute path to TCR file describing clonotypes enrichment:', clons.info.file, '\n')
cat('Absolute path to seurat object RDS file:', seurat.obj.file, '\n')
cat('Absolute path to output directory:', reports.path, '\n')
cat('Threshold to call clonally expanded cells:', expansion.thold, '\n')
cat(paste0('Tags to investigate further:\n\t', paste0(tags.of.int, collapse='\n\t'), '\n'))
cat('\n\n')


### --------------------------- Functions --------------------------- ###
# 1 -------------------------------------------------------------------->
# Name: Output clonotypes types values barplots across tag values.
# Description:
# Given a set of tag values defined in the seurat object within the environment, this function outputs to an output directory the frequency for expanded and non-expanded clonotypes from the cells across those tag values.
# Arguments ------------------------->
# output.path
# tag.values - Tag values to consider.
# cells.input - Logical indicating wether the output should be considered at the clonotypes (unqiue clonotype values per tag value) or cell (non-unique clonotypes) level.
# Function:

output.barplot <- function(output.path, tag.values, cells.input=TRUE){
  tag.clons.nos <- sapply(X=tag.values, FUN=function(tmp.tag.val){
    val.clons.info <- seurat.obj@meta.data[seurat.obj@meta.data[, tmp.tag]==tmp.tag.val, c('clonotype.tag', 'TCR.tag')]
    val.clons <- val.clons.info$clonotype.tag
    # If this is required for cell level, we don't need to restrict the analysis to unique clonotypes.
    if(!cells.input) val.clons <- unique(val.clons)
    # No matter what, we need to filter NA values.
    val.clons <- val.clons[!is.na(val.clons)]
    val.clons.no <- length(val.clons)
    # Take clonotypes with expanded tags.
    val.exp.clons <- val.clons.info$clonotype.tag[val.clons.info$TCR.tag=='Expanded']
    val.exp.clons <- val.exp.clons[!is.na(val.exp.clons)]
    val.exp.clons.no <- ifelse(test=cells.input, yes=length(val.exp.clons), no=length(unique(val.exp.clons)))
    # And same for non-expanded.
    val.non.exp.clons <- val.clons.info$clonotype.tag[val.clons.info$TCR.tag=='Non-expanded']
    val.non.exp.clons <- val.non.exp.clons[!is.na(val.non.exp.clons)]
    val.non.exp.clons.no <- ifelse(test=cells.input, yes=length(val.non.exp.clons), no=length(unique(val.non.exp.clons)))
    # Return results
    to.output <- c(val.clons.no, val.exp.clons.no, val.non.exp.clons.no)
    return(to.output)
  })
  tag.clons.nos <- as.data.frame(tag.clons.nos, stringsAsFactors=FALSE)
  tag.clons.nos$type <- c('all', 'expanded', 'non.expanded')
  # Get a backup of the table to output.
  to.output <- tag.clons.nos
  # Get tidy data.
  tag.clons.nos <- gather(data=tag.clons.nos, key=tag.value, value=clonotypes, tag.values)
  # Output barplots for this data.
  # @ Raw clonotypes numbers (across tag values) barplots.
  tmp.ggplot <- ggplot(data=tag.clons.nos, aes(x=tag.value, y=clonotypes, fill=type)) + geom_bar(stat='identity', position=position_dodge(width=1), width=0.6) + labs(x=tmp.tag.lab, y=ifelse(test=cells.input, yes='Cell count', no='Clonotype count'), fill='Type')
  tmp.file.name <- paste0(output.path, '/TagByClonotypeExpansionDodgedBarplot.pdf')
  pdf(file=tmp.file.name)
  print(tmp.ggplot + theme_bw())
  dev.off()
  # @ Clonotypes percentages (across tag values) barplots.
  # Filter certain clonotype type if required.
  tag.clons.nos <- tag.clons.nos[tag.clons.nos$type!='all', ]
  tmp.ggplot <- ggplot(data=tag.clons.nos, aes(x=tag.value, y=clonotypes, fill=type)) + geom_bar(stat='identity', position='fill', width=0.6) + scale_y_continuous(labels = scales::percent_format()) + coord_flip() + labs(x=tmp.tag.lab, y=ifelse(test=cells.input, yes='Cell percent', no='Clonotype percent'), fill='Type')
  tmp.file.name <- paste0(output.path, '/TagByClonotypeExpansionFilledBarplot.pdf')
  pdf(file=tmp.file.name)
  print(tmp.ggplot + theme_bw())
  dev.off()
  # ---> Output.
  # Return main table.
  return(to.output)
}

# 2 -------------------------------------------------------------------->
# Name: Heatmap for clonotype sharing among tag values.
# Description:
# Given a set of tag values defined in the seurat object within the environment, this function outputs to an output directory the heatmap depicting clonpotype sharing among tag values.
# Arguments ------------------------->
# output.preffix - Output file name preffix.
# tag.values - Tag values to consider.
# do.norm - Row-wise normalization.
# Function:

output.heatmap <- function(tag.values, output.preffix, do.norm=FALSE){
  # Get sharing table
  sharing.table <- sapply(X=tag.values, FUN=function(main.tag.val){
    main.tag.val.clons <- unique(seurat.obj@meta.data[, 'clonotype.tag'][seurat.obj@meta.data[, tmp.tag]==main.tag.val])
    val.sharing <- sapply(X=tag.values, FUN=function(second.tag.val){
      second.tag.val.clons <- unique(seurat.obj@meta.data[, 'clonotype.tag'][seurat.obj@meta.data[, tmp.tag]==second.tag.val])
      return(length(base::intersect(x=main.tag.val.clons, y=second.tag.val.clons)))
    })
    # Do normalization if required.
    if(do.norm) val.sharing <- val.sharing*100/length(main.tag.val.clons)
    # Then, return.
    return(val.sharing)
  })
  # Change diagonal values.
  diag(sharing.table) <- NA
  # Output table.
  tmp.file.name <- paste0(output.preffix, 'SharingTable.csv')
  write.csv(file=tmp.file.name, x=sharing.table)
  # Heatmap.
  melted.table <- melt(sharing.table)
  colnames(melted.table) <- c('Clusters', 'clusters', 'value')
  tmp.ggplot <- ggplot(data=melted.table, aes(x=Clusters, y=clusters, fill=value)) + geom_tile()
  tmp.file.name <- paste0(output.preffix, 'SharingHeatMap.pdf')
  pdf(file=tmp.file.name)
  print(tmp.ggplot)
  dev.off()
}


cat('\n\n')
### ------------------------- Data Loading ------------------------- ###
cat('### ------------------------- Data Loading ------------------------- ###\n')
# Seurat object ------------------------------------------------------->
seurat.obj <- readRDS(file=seurat.obj.file)
# Check that for the seurat object, we have dim. reductions methods applied.
if(!typeof(seurat.obj@reductions$umap)=="S4") stop('There is not either UMAP or tSNE applied to the input seurat object.')
# Check tags of interest are defined in seurat object's metadata.
if(!all(tags.of.int %in% colnames(seurat.obj@meta.data))){
  tags.of.int <- tags.of.int[!tags.of.int %in% colnames(seurat.obj@meta.data)]
  tmp.error <- paste0('Not all tags of interest defined in seurat object\'s metadata. Undefined tags follow:\n\t', paste0(tags.of.int, collapse='\n\t'), '\n')
  stop(tmp.error)
}
cat('Seurat object has been loaded.\n')
# TCR data files ------------------------------------------------------>
# ---> Clonotypes info
clons.info <- read.csv(file=clons.info.file, stringsAsFactors=FALSE)
cat('Both, TCR data and seurat object read to R objects. Check for warnings or errors if any.\n\n')
# ---> Cells-clonotypes relationships info.
cells.clons.info <- read.csv(file=cells.clons.info.file, stringsAsFactors=FALSE)
# Other files --------------------------------------------------------->
# ---> File listing the suffixes correspondences between the seurat object (corresponding to the gene expression data aggregation table) and the table listing the relationships between the cell barcodes and the clonotypes (corresponding to the TCR data aggregation table).
sffxs.corrs <- if(!is.null(sffxs.corrs.file)) read.csv(file=sffxs.corrs.file, stringsAsFactors=FALSE) else NULL


cat('\n\n')
### ---------------------- Data preprocessing ----------------------- ###

# ---> Aggregation suffixes.
# Change suffix if necessary.
if(!is.null(sffxs.corrs)){
  # Check suffixes correspondance table does meet the feature requirements.
  to.check <- c('aggr.table', 'seurat.obj')
  if(!all(to.check %chin% colnames(sffxs.corrs))){
    to.check <- to.check[!to.check %chin% colnames(sffxs.corrs)]
    tmp.err <- paste0('Suffixes correspondance table is not appropriately formatted. It lacks essential column(s):\n', paste0(to.check, collapse='\n'), '\n'); stop(tmp.err)
  }
  to.check <- length(unique(sffxs.corrs$aggr.table))==nrow(sffxs.corrs) & length(unique(sffxs.corrs$seurat.obj))==nrow(sffxs.corrs)
  if(!to.check){ tmp.err <- paste0('Values listed in the suffixes correspondance table should be unique for both elements, the suerat object and the aggregation table.\n'); stop(tmp.err) }
  # Check suffixes listed in the suffixes correspondance table are defined in both elements, the seurat object and the table listing the relationships between cell barcodes and clonotypes.
  to.check <- sffxs.corrs$seurat.obj %in% unique(str_extract(string=Cells(seurat.obj), pattern='\\d+$'))
  if(!all(to.check)){
    to.check <- sffxs.corrs$seurat.obj[!to.check]
    tmp.err <- paste0('Not all suffixes defined in the suffixes correspondance table for the seurat object element are appropriately defined in the actual input seurat object. Specifically, next suffixes are missing in the actual seurat object:\n', paste0(to.check, collapse='\n'), '\n'); stop(tmp.err)
  }
  to.check <- sffxs.corrs$aggr.table %in% unique(str_extract(string=cells.clons.info$barcode, pattern='\\d+$'))
  if(!all(to.check)){
    to.check <- sffxs.corrs$aggr.table[!to.check]
    tmp.err <- paste0('Not all suffixes defined in the suffixes correspondance table for the aggregation table element are appropriately defined in the table listing the relationships between cell barcodes and clonotypes. Specifically, next suffixes are missing in the latter:\n', paste0(to.check, collapse='\n'), '\n'); stop(tmp.err)
  }
  # Proceed to transform suffixes found in the table listing the relationships between cells and clonotypes.
  cells.clons.info$ref.barcode <- cells.clons.info$barcode
  for(idx in 1:nrow(sffxs.corrs)){
    tmp.pttn <- paste0('-', sffxs.corrs[idx, 'aggr.table'], '$'); tmp.replacement <- paste0('-', sffxs.corrs[idx, 'seurat.obj'])
    idxs.to.modify <- str_detect(string=cells.clons.info$ref.barcode, pattern=tmp.pttn)
    cells.clons.info[idxs.to.modify, 'barcode'] <- str_replace_all(string=cells.clons.info[idxs.to.modify, 'ref.barcode'], pattern=tmp.pttn, replacement=tmp.replacement)
  }
  cells.clons.info$ref.barcode <- NULL
}

# ---> Cells-clonotypes relationships info (filtering).
# Keep track of relationships deleted due to a chain different than alpha or beta.
chains.info <- as.data.frame(table(cells.clons.info$chain))
colnames(chains.info) <- c('chain', 'freq')
# Also, keep track of how many relationships we had previous to any filtering.
pre.filtering.cells <- sum(Cells(seurat.obj) %in% cells.clons.info$barcode)
pre.filtering.rels <- nrow(cells.clons.info)
# Filter out non-productive or not-of-interest contigs from the cells-clonotypes info.
# Conditions:
# * For barcodes called as cells.
# * For clonotype contigs marked with high confidence.
# * For chains ultimately defined as TRA or TRB.
# * For productive clonotypes (see 10X support website for a broader definition of 'productive').
# Also, we'll take out this info since it has been taken into consideration already.
cells.to.keep <- cells.clons.info$is_cell=='True' & cells.clons.info$high_confidence=='True' & (cells.clons.info$chain=='TRA' | cells.clons.info$chain=='TRB') & cells.clons.info$productive=='True'
feats.to.keep <- c('barcode', 'contig_id', 'chain', 'v_gene', 'd_gene', 'j_gene', 'cdr3', 'cdr3_nt', 'reads', 'umis', 'raw_clonotype_id')
cells.clons.info <- cells.clons.info[cells.to.keep, feats.to.keep]


cat('\n\n')
### ------------------------- Main program ------------------------- ###
cat('### ------------------------- Main program ------------------------- ###\n\n')

cat('\n\n')
### ------------------ General reports and cleaning ----------------- ###
cat('### ------------------ General reports and cleaning ----------------- ###\n')

gen.reports.path <- paste0(reports.path, '/general_reports_and_cleaning')
create.dir(dir.path=gen.reports.path, path.desc='General reports and cleaning')

# ---> General reports.
cat('# ---> General reports.\n')
# First of all, report chains info previous to filtering.
tmp.file.name <- paste0(gen.reports.path, '/ChainTypesFreqBe4Filtering.csv')
write.csv(file=tmp.file.name, x=chains.info)
# Then, more general reports.
# @ Total number of cells.
total.cells <- length(Cells(seurat.obj))
tmp.int <- total.cells
tmp.text <- 'Cell barcodes in this dataset'
to.output <- data.table(concept=tmp.text, count=tmp.int, percentage=100)
# @ Cells with TCR info.
tmp.int <- sum(Cells(seurat.obj) %in% cells.clons.info$barcode)
tmp.pct <- round(x=tmp.int*100/total.cells, digits=3)
tmp.text <- 'Cell barcodes annotated at least with one TCR (alpha; beta or both) chain'
to.output <- rbind(to.output, data.table(concept=tmp.text, count=tmp.int, percentage=tmp.pct))
# @ Cells with total TCR info.
tmp.int <- sum(colnames(seurat.obj) %in% cells.clons.info[cells.clons.info$chain=='TRA', 'barcode'] & colnames(seurat.obj) %in% cells.clons.info[cells.clons.info$chain=='TRB', 'barcode'])
tmp.pct <- round(x=tmp.int*100/total.cells, digits=3)
tmp.text <- 'Cell barcodes annotated with both both kinds of TCR chains (alpha and beta)'
to.output <- rbind(to.output, data.table(concept=tmp.text, count=tmp.int, percentage=tmp.pct))
# @ TRB chains.
tmp.int <- sum(colnames(seurat.obj) %in% cells.clons.info[cells.clons.info$chain=='TRB', 'barcode'])
tmp.pct <- round(x=tmp.int*100/total.cells, digits=3)
tmp.text <- 'Cell barcodes annotated at least with a TCR beta chain sequence'
to.output <- rbind(to.output, data.table(concept=tmp.text, count=tmp.int, percentage=tmp.pct))
# @ TRA chains.
tmp.int <- sum(colnames(seurat.obj) %in% cells.clons.info[cells.clons.info$chain=='TRA', 'barcode'])
tmp.pct <- round(x=tmp.int*100/total.cells, digits=3)
tmp.text <- 'Cell barcodes annotated at least with a TCR alpha chain sequence'
to.output <- rbind(to.output, data.table(concept=tmp.text, count=tmp.int, percentage=tmp.pct))
# Output general report.
# Prepare file.
tmp.file.name <- paste0(gen.reports.path, '/BarcodesTCRsAssociationsSummary.csv')
fwrite(file=tmp.file.name, x=to.output)
cat('Summary of the relations between cell barcodes and TCR sequences have been output.\n\n')

# ---> Filtering of Cells-clonotypes associations.
cat('# ---> Filtering of Cells-clonotypes associations.\n')
# For downstream analyses, we'll keep just the clonotypes with a valid associated barcode (i.e., those called as cells with good data quality.)
cell.clons.info <- cells.clons.info[cells.clons.info$barcode %in% Cells(seurat.obj), ]
cat('For this to be done, we will just take into consideration the clonotypes for which a valid barcode (i.e., one called a cell) was assigned.\n')
# Then, based on that data structure, keep just the clonotype data for the cells of interest (i.e., those captured in the seurat object).
clons.info <- clons.info[clons.info$clonotype_id %in% cells.clons.info$raw_clonotype_id, ]
# ---> Update clonotype proportion.
# Also, update clonotype proportions considering just the ones found for our cells in the seurat object.
total.freq <- sum(clons.info$frequency)
clons.info$updated_proportion <- clons.info$frequency/total.freq
# ---> Clonally expanded and non-clonally expanded clonotypes.
# Add column to tell either the clonotype is or not clonally expanded.
clons.info$expanded <- clons.info$frequency >= expansion.thold
# The total amount of clonally expanded clonotypes.
total.expanded <- sum(clons.info$expanded)
# ---> Cells with multiple clonotypes assigned to them.
# Idenitify cells with multiple clonotypes assigned to them. If any, we get rid of their clonotypes for downstream analyses and report them to an output file.
no.clons.per.cell <- sapply(X=Cells(seurat.obj), FUN=function(cell){
  return(length(unique(cells.clons.info$raw_clonotype_id[cells.clons.info$barcode==cell])))
})
cells.with.mult.clons <- names(no.clons.per.cell)[no.clons.per.cell>1]
# Retrieve data for cells with multiple assigned clonotypes and filter them out since they cause trouble for downstream analyses (mainly, we're interested in assigning a single clonotype for each cell).
if(length(cells.with.mult.clons)>0) cells.clons.info <- cells.clons.info[!cells.clons.info$barcode %in% cells.with.mult.clons, ]
# ---> Add new info to general reports.
to.output <- data.table(concept=c('Total frequency', 'Unique clonotypes count', 'Clonotypes showing clonal expansion', 'Cell barcodes with more than one clonotype associated to them'), count=c(total.freq, nrow(clons.info), total.expanded, length(cells.with.mult.clons)))
tmp.file.name <- paste0(gen.reports.path, '/BarcodesAndTCRAssociationsSummary.csv')
fwrite(x=to.output, file=tmp.file.name)
cat('Preprocessing step complete.\n')


cat('\n\n')
### --------------------- Clonotypes proportions -------------------- ###
cat('### --------------------- Clonotypes proportions -------------------- ###\n')

clons.props.path <- paste0(reports.path, '/clonotypes_proportions')
create.dir(dir.path=clons.props.path, path.desc='Clonotypes proportons')

# ---> Clonally expanded vs not clonally expanded proportion
cat('# ---> Clonally expanded vs not clonally expanded proportion\n')
# Create data frame in the appropriate format to output a piechart.
tmp.df <- data.frame(Group=c('Clonally expanded', 'Non clonally expanded'), value=c(total.expanded, nrow(clons.info)-total.expanded))
# Output to a file
tmp.file.name <- paste0(clons.props.path, '/PieChartExpandedvsNonExpanded.pdf')
output.piechart(tmp.df=tmp.df, file.name=tmp.file.name, piechart.title='Proportions of clonally expanded and non-clonally expanded TCRs.')
cat('Piechart depicting the proportions of clonally expanded and not clonally expanded TCRs output.\n')

# ---> Single frequencies of clonally expanded TCRs
cat('# ---> Single frequencies of clonally expanded TCRs\n')
# Same, but for clonotypes clonally expanded, the one with the highes frequency, the second one and so on according to the amount.of.clons input variable.
# For this purpose, order the rows of the clonotypes info according to clonotypes' frequency.
expanded.clons.freq <- clons.info[clons.info$expanded, c('frequency', 'updated_proportion')] %>% arrange(desc(frequency))
tmp.sum <- sum(expanded.clons.freq[(amount.clons+1):nrow(expanded.clons.freq), 'frequency'])
tmp.df <- data.frame(Group=c(ordinal(1:amount.clons), 'Rest'), value=c(expanded.clons.freq[1:amount.clons, 'frequency'], tmp.sum))
tmp.file.name <- paste0(clons.props.path, '/PieChartFreqsForExpandedClons.pdf')
output.piechart(tmp.df=tmp.df, file.name=tmp.file.name, text.col='black', piechart.title='Frequencies of expanded clonotypes')
cat('Piechart depicting the proportions of the', amount.clons, 'TCRs with the higest frequencies and the rest output.\n\n')

# ---> Highly clonally expanded TCRs.
cat('# ---> Highly clonally expanded TCRs.\n')
# Proportions of highly clonally expanded TCRs (freq. > 4) vs barely clonally expanded TCRs (freq. < 4).
no.highly.exp.clons <- sum(expanded.clons.freq$frequency > 4)
no.barely.exp.clons <- nrow(expanded.clons.freq) - no.highly.exp.clons
tmp.df <- data.frame(Group=c('freq.>4', 'freq.<=4'), value=c(no.highly.exp.clons, no.barely.exp.clons))
tmp.file.name <- paste0(clons.props.path, '/PieChartPropsHighlyvsBarelyExpClons.pdf')
output.piechart(tmp.df=tmp.df, file.name=tmp.file.name, text.col='white', piechart.title='Highly expanded vs other clonotypes')


### ------------------ Mark clonotypes on UMAP/tSNE ----------------- ###
cat('### ------------------ Mark clonotypes on UMAP/tSNE ----------------- ###\n')

dim.reduction.path <- paste0(reports.path, '/dim_reduction')
create.dir(dir.path=dim.reduction.path, path.desc='Dimensionality reduction')

# ---> Add tags to seurat object.
cat('# ---> Add tags to seurat object.\n')
# Add TCR info to seurat object meta data.
# For each cell, we'll get:
#   TRA.tag, TRB.tag, TCR.tag (combination of TRA and TRB) each with three possible values: Expanded, Non-expanded or NA
#   CMC (chain-related multiplet clonotype) tag.
#   Chains sequences tags, for both, aa and nt, levels.

clons.info.for.cells <- as.data.frame(t(sapply(X=Cells(seurat.obj), FUN=function(cell){
  # Chains info.
  what.chain <- cells.clons.info$chain[cells.clons.info$barcode==cell]
  # Clonotype ID.
  clon.id <- unique(cells.clons.info$raw_clonotype_id[cells.clons.info$barcode==cell])
  if(length(clon.id)>0){
    # Expansion info.
    if(clons.info$frequency[clons.info$clonotype_id==clon.id]>=expansion.thold) expanded <- 'Expanded' else expanded <- 'Non-expanded'
    # Chains sequences.
    TRA.nt.chains <- get.cdr3s.from.table(cdr3s.table=clons.info[clons.info$clonotype_id==clon.id, ], cdr3.col='cdr3s_nt', chain.ptn='TRA')
    TRB.nt.chains <- get.cdr3s.from.table(cdr3s.table=clons.info[clons.info$clonotype_id==clon.id, ], cdr3.col='cdr3s_nt', chain.ptn='TRB')
    TRA.aa.chains <- get.cdr3s.from.table(cdr3s.table=clons.info[clons.info$clonotype_id==clon.id, ], cdr3.col='cdr3s_aa', chain.ptn='TRA')
    TRB.aa.chains <- get.cdr3s.from.table(cdr3s.table=clons.info[clons.info$clonotype_id==clon.id, ], cdr3.col='cdr3s_aa', chain.ptn='TRB')
    # Clonotype size.
    clon.size <- unique(clons.info[clons.info$clonotype_id==clon.id, 'frequency'])
    clon.prop <- unique(clons.info[clons.info$clonotype_id==clon.id, 'proportion'])
    # CMC-related info.
    if(length(TRA.nt.chains) > 1 | length(TRB.nt.chains) > 1) cmc.tag <- TRUE else cmc.tag <- FALSE
    if(length(unique(what.chain))==2){
      to.return <- c(clon.id, rep(expanded, times=3), cmc.tag, paste0(TRA.nt.chains, collapse=';'), paste0(TRB.nt.chains, collapse=';'), paste0(TRA.aa.chains, collapse=';'), paste0(TRB.aa.chains, collapse=';'), clon.size, clon.prop)
    }else{
      # The TCR tag will indicate if there's any chain marked as expanded.
      if('TRA' %in% what.chain) to.return <- c(clon.id, expanded, NA, expanded, cmc.tag, paste0(TRA.nt.chains, collapse=';'), NA, paste0(TRA.aa.chains, collapse=';'), NA, clon.size, clon.prop) else to.return <- c(clon.id, NA, expanded, expanded, cmc.tag, NA, paste0(TRB.nt.chains, collapse=';'), NA, paste0(TRB.aa.chains, collapse=';'), clon.size, clon.prop)
    }
  }else{
    to.return <- rep(NA, times=11)
  }
  names(to.return) <- c('clonotype.tag', 'TRA.tag', 'TRB.tag', 'TCR.tag', 'CMC.tag', 'TRA.nt.chains.tag', 'TRB.nt.chains.tag', 'TRA.aa.chains.tag', 'TRB.aa.chains.tag', 'clon.size.tag', 'clon.proportion.tag')
  return(to.return)
})))
clons.info.for.cells[, 'CMC.tag'] <- as.logical(clons.info.for.cells[, 'CMC.tag'])
clons.info.for.cells[, 'clonotype.tag'] <- as.character(clons.info.for.cells[, 'clonotype.tag'])
clons.info.for.cells[, 'clon.size.tag'] <- as.numeric(as.character(clons.info.for.cells[, 'clon.size.tag']))
clons.info.for.cells[, 'clon.proportion.tag'] <- as.numeric(as.character(clons.info.for.cells[, 'clon.proportion.tag']))
seurat.obj@meta.data <- cbind(seurat.obj@meta.data, clons.info.for.cells)

# ---> Explore distribution for clone size.
# Density plot.
tmp.ggplot <- ggplot(data=seurat.obj@meta.data, aes(x=clon.size.tag)) + geom_density(alpha=0.6, fill='#ff4d4d') + scale_x_continuous(trans='log10')
tmp.file.name <- paste0(dim.reduction.path, '/CloneSizeDistribution.pdf')
pdf(file=tmp.file.name)
print(tmp.ggplot + theme_bw())
dev.off()
# Summary stats.
up.thold <- quantile(x=seurat.obj@meta.data[, 'clon.size.tag'], probs=0.95, na.rm=TRUE)

# ---> Depict new tags on Dimensionality reduction maps.
cat('# ---> Depict new tags in a DimPlot.\n')
# -@ Clone size.
tmp.file.name <- paste0(dim.reduction.path, '/CloneSizeOnFeaturePlot.pdf')
tmp.caption <- paste0('Maximum scale value is 95 percentile, ', up.thold)
cells.to.depict <- Cells(seurat.obj)[!is.na(seurat.obj@meta.data[, 'clon.size.tag']) & seurat.obj@meta.data[, 'clon.size.tag'] >= expansion.thold]
tmp.col.gradient <- if(up.thold > expansion.thold) scale_colour_gradientn(trans='log2', colours=this.color.scale, limits=c(expansion.thold, up.thold)) else scale_colour_gradientn(trans='log2', colours=this.color.scale)
tmp.ggplot <- FeaturePlot(seurat.obj, reduction='umap', feature='clon.size.tag', cells=cells.to.depict, min.cutoff=3, max.cutoff=up.thold) + scale_alpha(0.7) + tmp.col.gradient + labs(title='Clone size', x='UMAP 1', y='UMAP 2', col='Size', caption=tmp.caption)
pdf(file=tmp.file.name)
print(tmp.ggplot)
dev.off()

# -@ Categorical variables.
# All cells with TCR info.
tmp.file.name <- paste0(dim.reduction.path, '/TCRTagsOnDimPlot.pdf')
pdf(file=tmp.file.name)
for(tmp.tag in paste0(c('TRA', 'TRB', 'TCR', 'CMC'), '.tag')) print(DimPlot(seurat.obj, reduction='umap', group.by=tmp.tag, cells=Cells(seurat.obj)[!is.na(seurat.obj@meta.data[,tmp.tag])])+labs(title=tmp.tag)+scale_alpha(0.5))
dev.off()
# Just expanded cells.
tmp.file.name <- paste0(dim.reduction.path, '/TCRTagsOnDimPlotJustForExpandedCells.pdf')
pdf(file=tmp.file.name)
for(tmp.tag in paste0(c('TRA', 'TRB', 'TCR', 'CMC'), '.tag')){
  cells.to.depict <- Cells(seurat.obj)[seurat.obj@meta.data[,tmp.tag]=='Expanded' | seurat.obj@meta.data[,tmp.tag]==TRUE]
  cells.to.depict <- cells.to.depict[!is.na(cells.to.depict)]
  print(DimPlot(seurat.obj, reduction='umap', group.by=tmp.tag, cells=cells.to.depict)+labs(title=tmp.tag)+scale_alpha(0.7))
}
dev.off()
# -@ Numerical variables.
tmp.file.name <- paste0(dim.reduction.path, '/TCRTagsOnFeaturePlot.pdf')
pdf(file=tmp.file.name)
for(tmp.tag in paste0(c('clon.size', 'clon.proportion'), '.tag')){
  cells.to.depict <- Cells(seurat.obj)[!is.na(seurat.obj@meta.data[ ,tmp.tag])]
  # Normal scale.
  print(FeaturePlot(seurat.obj, reduction='umap', feature=tmp.tag, cells=cells.to.depict)+labs(title=tmp.tag)+scale_alpha(0.7))
  # Log scale.
  print(FeaturePlot(seurat.obj, reduction='umap', feature=tmp.tag, cells=cells.to.depict)+labs(title=tmp.tag)+scale_alpha(0.7)+scale_colour_gradient(trans='log2')+labs(subtitle='log2 scale'))
}
dev.off()


cat('\n\n')
if(!is.null(tags.of.int)){
### --------------------- Tag-specific analyses --------------------- ###
cat('### --------------------- Tag-specific analyses --------------------- ###\n')

tags.path <- paste0(reports.path, '/tag_specific_analysis')
create.dir(dir.path=tags.path, path.desc='Tags-specific analyses')

# Regarding TCR info.
for(tmp.tag in tags.of.int){
  cat('\n---> Process for tag', tmp.tag, '\n')
  # ------ TCR expansion accross tag values ------ #
  # This tag path
  tmp.tag.path <- paste0(tags.path, '/', tmp.tag)
  tmp.tag.lab <- str_to_title(str_replace(string=str_replace_all(string=tmp.tag, pattern='\\.', replacement=' '), pattern=' tag', replacement=''))
  # Get all tag values (different than NA).
  tag.values <- mixedsort(unique(as.character(seurat.obj@meta.data[, tmp.tag])), decreasing=FALSE)
  tag.values <- tag.values[!is.na(tag.values)]
  # Width length to increase to output plots for a tag based on its count of values.
  plus.width <- if(length(tag.values) < 6) 0 else length(tag.values)*0.2

  # if(is.factor(seurat.obj@meta.data[, tmp.tag])) seurat.obj@meta.data[, tmp.tag] <- factor(x=seurat.obj@meta.data[, tmp.tag], levels=tag.values)
  seurat.obj@meta.data[, tmp.tag] <- factor(x=seurat.obj@meta.data[, tmp.tag], levels=tag.values)

  # Skip in case there aren't at least two tag values to evaluate.
  if(length(tag.values)<2) next else create.dir(dir.path=tmp.tag.path, path.desc=tmp.tag)

  # ---> Clone size across cells in the same group.
  clones.tag.path <- paste0(tmp.tag.path, '/clone_size')
  create.dir(dir.path=clones.tag.path, path.desc=paste0(tmp.tag, ' for clone size'))
  # Box plot.
  tmp.ggplot <- ggplot(data=seurat.obj@meta.data[!is.na(seurat.obj@meta.data[, tmp.tag]), ], aes_string(y='clon.size.tag', x=tmp.tag, fill=tmp.tag)) + geom_boxplot(outlier.shape=NA, alpha=0.7, width=0.3) + scale_y_continuous(trans='log10')  + labs(title=paste0('Clone size across groups'), x=tmp.tag.lab, y='Clone size')
  tmp.file.name <- paste0(clones.tag.path, '/CloneSizeAcrossTagValuesBoxplots.pdf')
  pdf(tmp.file.name, width=7+plus.width)
  print(tmp.ggplot + theme_bw() + theme(legend.position='none'))
  dev.off()
  # Violin plot.
  tmp.ggplot <- vln.plot(seurat.obj=seurat.obj, feature='clon.size.tag', groups.tag=tmp.tag, na.rm=TRUE, groups.of.int=NULL, filter.tags=NULL, color='median', trim.val=TRUE)
  tmp.ggplot <- tmp.ggplot + labs(x=tmp.tag, y='Clone size', fill='Size median') + theme_bw()
  tmp.file.name <- paste0(clones.tag.path, '/CloneSizeAcrossTagValuesViolinPlots.pdf')
  pdf(tmp.file.name, width=7+plus.width)
  print(tmp.ggplot + theme_bw())
  dev.off()

  # ---> Expansion for cells across tag values.
  tmp.cells.tag.path <- paste0(tmp.tag.path, '/cells_expansion')
  create.dir(dir.path=tmp.cells.tag.path, path.desc=paste0(tmp.tag, ' for cells expansion'))
  # Tag by clonotypes info for cells.
  # i.e., we don't care about unique clonotypes.
  tmp <- output.barplot(output.path=tmp.cells.tag.path, tag.values=tag.values, cells.input=TRUE)
  # This should be changed in order to use tmp rather than the next matrix.
  cells.clons.by.tag <- as.matrix(table(seurat.obj@meta.data[, c('TCR.tag', tmp.tag)]))
  # Normalization per tag.
  tag.totals <- colSums(cells.clons.by.tag)
  norm.factors <- tag.totals/min(tag.totals)
  normt.cells.clons.by.tag <- t(apply(X=cells.clons.by.tag, MARGIN=1, FUN=function(tag.row) tag.row/norm.factors))
  # Output both tables.
  tmp.file.name <- paste0(tmp.cells.tag.path, '/TagByClonotypeCellsExpansionTable.csv')
  write.csv(file=tmp.file.name, x=cells.clons.by.tag)
  tmp.file.name <- paste0(tmp.cells.tag.path, '/TagByClonotypeCellsExpansionTableAfterNormalization.csv')
  write.csv(file=tmp.file.name, x=normt.cells.clons.by.tag)

  # ---> Expansion for clonotypes across tag values.
  tmp.clons.tag.path <- paste0(tmp.tag.path, '/clons_expansion')
  create.dir(dir.path=tmp.clons.tag.path, path.desc=paste0(tmp.tag, ' for clonotypes expansion'))
  tmp <- output.barplot(output.path=tmp.clons.tag.path, tag.values=tag.values, cells.input=FALSE)

  # --- Clonotypes sharing across tag values --- #
  tmp.share.tag.path <- paste0(tmp.tag.path, '/clonotypes_sharing')
  create.dir(dir.path=tmp.share.tag.path, path.desc=paste0(tmp.tag, ' for sharing among values'))

  # ---> Venn diagram.
  # Clonotypes sharing among tag values' cells.
  # Process to be applied for different counts of sets to take into account. When there are too many values for a given tag, we keep only the first 8 values to make the process more managable.
  set.counts <- if(length(tag.values)>=5) 2:5 else 2:length(tag.values)
  if(length(tag.values)>8) tag.values <- tag.values[1:8]
  for(set.count in set.counts){
    set.combns <- combn(x=tag.values, m=set.count, simplify=FALSE)
    tmp.file.name <- paste0(tmp.share.tag.path, '/ClonotypesSharingByVennDiagram_ForSetCount-', set.count, '.pdf')
    pdf(file=tmp.file.name)
    for(set.combn in set.combns){
      # Set data.
      combn.lab <- paste0('Combination ', paste0(set.combn, collapse='-'))
      tmp.data.1 <- lapply(X=set.combn, FUN=function(tag.val){
        tag.val.cells <- seurat.obj@meta.data[, tmp.tag]==tag.val & !is.na(seurat.obj@meta.data[, 'clonotype.tag'])
        to.output <- unique(seurat.obj@meta.data[tag.val.cells, 'clonotype.tag'])
        return(to.output)
      })
      names(tmp.data.1) <- set.combn
      # Calculate P-value based on Fisher's exact test.
      #     Process done based on the instructions from: https://rdrr.io/bioc/GeneOverlap/man/GeneOverlap.html
      if(set.count==2){
        tmp.data.2 <- matrix(nrow=2, data=c(
          nrow(clons.info) - length(union(x=tmp.data.1[[1]], y=tmp.data.1[[2]])),
          length(setdiff(x=tmp.data.1[[1]], y=tmp.data.1[[2]])),
          length(setdiff(x=tmp.data.1[[2]], y=tmp.data.1[[1]])),
          length(intersect(x=tmp.data.1[[1]], y=tmp.data.1[[2]]))
        ))
        test.results <- fisher.test(x=tmp.data.2, alternative='greater')
        tmp.sub <- paste0('Fisher\'s exact test: ', round(x=test.results$p.value, digits=2))
      }else{
        tmp.sub <- ''
      }
      # Get Venn diagram.
      tmp.plot <- venn.diagram(x=tmp.data.1, file=NULL, main=combn.lab, sub=tmp.sub, force.unique=TRUE, na='remove')
      # tmp.plot <- venn.diagram(x=tmp.data.1, file=NULL, main=combn.lab, force.unique=TRUE, na='remove', hyper.test=TRUE, total.population=nrow(clons.info))
      grid.newpage(); grid.draw(tmp.plot)
    }
    dev.off()
  }

  # ---> Cells' sharing according to cells' tag values.
  # Retrieve data of interest as a data table.
  tmp.data.1 <- as.data.table(seurat.obj@meta.data[, c('clonotype.tag', tmp.tag)]); colnames(tmp.data.1) <- c('clonotype.tag', 'tmp.tag')
  tmp.data.2 <- tmp.data.1[, .(sharing.degree.tag=uniqueN(tmp.tag)), by=.(clonotype.tag)]
  tmp.data.2[, tcr.sharing.tag:=ifelse(test=sharing.degree.tag>1, yes='Shared', no='Unique')]
  tmp.data <- merge(tmp.data.1, tmp.data.2, by='clonotype.tag', sort=FALSE, all=TRUE)
  tmp.data <- tmp.data[!is.na(clonotype.tag)]
  # Depicted as a barplot.
  tmp.col.scale <- binary.col.scale; names(tmp.col.scale) <- c('Shared','Unique')
  tmp.ggplot <- ggplot(data=tmp.data) + scale_y_continuous(expand=c(0, 0))
  tmp.ggplot.1 <- tmp.ggplot + geom_bar(aes(x=tmp.tag, fill=as.character(sharing.degree.tag)), position='fill', width=0.7) + scale_fill_brewer(palette='Set1') + labs(x=tmp.tag.lab, y='Proportion', fill='Sharing degree')
  tmp.ggplot.2 <- tmp.ggplot + geom_bar(aes(x=tmp.tag, fill=tcr.sharing.tag), position='fill', width=0.7) + scale_fill_manual(values=tmp.col.scale) + labs(x=tmp.tag.lab, y='Proportion', fill='Sharing status')
  # Output.
  tmp.file.name <- paste0(tmp.share.tag.path, '/CellsSharing.pdf')
  pdf(file=tmp.file.name)
  print(tmp.ggplot.1 + theme_bw())
  print(tmp.ggplot.2 + theme_bw())
  dev.off()
}
# ---> Remove trash produced to create venn diagrams.
files.to.rm <- list.files(path=getwd(), pattern='VennDiagram.+\\.log', all.files=FALSE, full.names=TRUE, recursive=FALSE, include.dirs=FALSE)
file.remove(files.to.rm)
cat('Tag-specific analyses finished!\n')
}

if(output.obj){
  cat('\n\n')
  ### ------------------------- Seurat object ------------------------- ###
  cat('### ------------------------- Seurat object ------------------------- ###\n')
  # Add a suffix if no value is provided.
  if(is.null(prj.sffx)) prj.sffx <- paste0('ChangedOn_', Sys.Date())
  # Save seurat object changing the seurat object name.
  tmp.file.name <- sub(x=seurat.obj.file, pattern='\\.RDS', ignore.case=TRUE, replacement=paste0('_WithTCRTags_', prj.sffx, '.RDS'))
  saveRDS(object=seurat.obj, file=tmp.file.name)
  cat('Seurat object saved!\n')
}

cat('\n\nSession info:')
sessionInfo()

cat('\n\nProgram finished.\n')

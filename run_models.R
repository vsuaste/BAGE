#!/usr/bin/env Rscript

user_lib <- "~/Rlibs"
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))

install_if_missing <- function(pkgs) {
  for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
      library(pkg, character.only = TRUE)
    }
  }
}

pkgs <- c("tidyverse",
          "stats",
          "BGLR",
          "optparse",
          "caret",
          "BiocManager")

install_if_missing(pkgs)
BiocManager::install(version='3.20')
BiocManager::install("DESeq2")
BiocManager::install("org.Hs.eg.db")
source('helpers.R')

option_list = list(
  make_option(c("-c", "--counts"), type="character", default=NULL, 
              help="counts file name", metavar="character"),
  make_option(c("-s", "--samples"), type="character", default="healthy_samples_filt_v5.tsv", 
              help="samples file name", metavar="character"),
  make_option(c("-m", "--out_metrics"), type="character", default="metrics", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-o", "--out_dir"), type="character", default="./report", 
              help="output path directory [default= %default]", metavar=NULL),
  make_option(c("-b", "--best_index"), type="character", default="best_index", 
              help="best index file name [default= %default]", metavar=NULL),
  make_option(c("-t", "--age_transformation"), type="character", default="none", 
              help="transformation performed to age", metavar="character"),
  make_option(c("-T", "--counts_transformation"), type="character", default="none", 
              help="transformation performed to counts", metavar="character"),
  make_option(c("-g", "--out_genes"), type="character", default="sorted_genes", 
              help="output genes file", metavar="character"),
  make_option(c("-n", "--num_genes"), type="integer", default=2000, 
              help="number of genes in the model", metavar="integer"),
  make_option(c("-N", "--num_runs"), type="integer", default=5, 
              help="k-fold number", metavar="integer"),
  make_option(c("-i", "--interaction"), type = "logical", default = FALSE,
              help = "Include gene-gender interaction"),
  make_option(c("-p", "--step"), type = "integer", default = 1,
              help = "increment to the number of genes considered in the model"),
  make_option(c("-f", "--filter"), type = "logical", default = TRUE,
              help = "filter low counts genes"),
  make_option(c("-l", "--loo"), type = "logical", default = FALSE,
              help = "leave one out prediction"),
  make_option(c("-d", "--deseq"), type = "logical", default = FALSE,
              help = "deseq selection variable case"),
  make_option(c("-v", "--variables"), type="character", default="raw_age_deseq_selection_loo.tsv", 
              help="deseq variables selection file", metavar="character"),
  make_option(c("-D", "--dataset"), type = "logical", default = FALSE,
              help = "include dataset in the model")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

counts <- opt$counts
samples <- opt$samples
out_metrics <- opt$out 

all_counts <- read_tsv(opt$counts)
all_counts <- column_to_rownames(all_counts, var="GENEID")

samples_healthy <- read_tsv(opt$samples)
samples_healthy <- column_to_rownames(samples_healthy, var="id")

counts_healthy <- all_counts[  , rownames(samples_healthy) ]


## create directories
if (!dir.exists(opt$out_dir)) {
  dir.create(opt$out_dir, recursive = TRUE)
}

metrics_dir <- paste0(opt$out_dir,"/metrics")
genes_dir <- paste0(opt$out_dir,"/genes")
best_index_dir <- paste0(opt$out_dir,"/best")

if (!dir.exists(metrics_dir)) {
  dir.create(metrics_dir, recursive = TRUE)
}

if (!dir.exists(genes_dir)) {
  dir.create(genes_dir, recursive = TRUE)
}

if (!dir.exists(best_index_dir)) {
  dir.create(best_index_dir, recursive = TRUE)
}

set.seed(37)

folds <- createFolds(samples_healthy$age, k=opt$num_runs, list = TRUE, returnTrain = FALSE)

# transform age
samples_healthy$age_t <- samples_healthy$age

if(opt$age_transformation == "log"){
  samples_healthy$age_t <- log(samples_healthy$age)
}

if(opt$age_transformation == "sqrt"){
  samples_healthy$age_t <- sqrt(samples_healthy$age)
}

# transform counts
all_counts_t <- counts_healthy

#filter zeros
if(opt$filter ==TRUE){
  all_counts_t <- prefilter_counts(counts_healthy, zero_percentage = 0.4, entrez = FALSE )
}else{
  print("No filtered counts")
}


if(opt$counts_transformation == "log"){
  all_counts_t  <- log10( all_counts_t  +1)
}


if(opt$deseq == TRUE){
  selected_genes <- read_tsv(opt$variables)
}

if(opt$dataset){
  form <- age_t ~ gender+ dataset +gene_expression
}else{ 
  form <- age_t ~ gender+gene_expression
}

#number of times the model will run
ntimes <- 1 + (opt$num_genes - 50)/opt$step

## run model 
#opt$num_runs
for( k in 1:opt$num_runs){
  metrics_path <- paste0( metrics_dir,"/metrics_", k, ".tsv" )
  best_index_path <- paste0(best_index_dir, "/best_index_", k, ".tsv" )
  genes_path <- paste0( genes_dir, "/sorted_genes_", k, ".tsv" )
  
  test_samples <- folds[[k]]
  train_samples <- samples_healthy[-test_samples,]
  
  ## select variable only from training set
  counts_filt <- all_counts_t[ , rownames(train_samples) ]
  counts_filt_raw <- counts_filt
  
  
  if(opt$deseq){
    sorted_genes <- data.frame(
      id = selected_genes[,k, drop=FALSE]
    )
    colnames(sorted_genes)<- c("id")
  }else{

    results <- do.call(rbind, lapply(rownames(counts_filt_raw), process_gene,
                                 counts = counts_filt_raw,
                                 samples = train_samples,
                                 form = form))

    sorted_genes <-  results[order(results$pv), ]                             
    
    write_tsv(sorted_genes, genes_path)
  }

  ngenes_seq <- 50 + (0:(ntimes-1)) * opt$step

  results_list <- lapply(ngenes_seq, function(ng) {
    genes <- sorted_genes[1:ng, ,drop=FALSE]
    info_genes <- data.frame(genes = genes$id)
    predicted_values(all_counts_t, info_genes, samples_healthy, test_samples,
                    t = opt$age_transformation,
                    dataset = opt$dataset,
                    interact = opt$interaction)
  })

  # bind results into a matrix
  predicted_matrix <- do.call(cbind, results_list)

  metrics <- compute_metrics(predicted_matrix, samples_healthy$age ,test_samples, opt$loo)
  
  write_tsv(metrics, metrics_path)
  
  if(opt$loo == FALSE){
      
    max_index <- apply(metrics, 2, which.max)
    min_index <- apply(metrics, 2, which.min)
    
    max_val <- apply(metrics, 2, max)
    min_val <- apply(metrics, 2, min)
    
    best_index <- data.frame(best_index =c(max_index[ 1:4 ], min_index[5:10]),
                             best_value =c(max_val[ 1:4 ], min_val[5:10] ) )
    
    best_index <- rownames_to_column(best_index, var="metric")
    
    write_tsv(best_index, best_index_path) 
  }
}






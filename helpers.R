get_entrez <- function(ids_to){
  mapIds(org.Hs.eg.db,
         keys = ids_to,
         column="ENTREZID",
         keytype = "ENSEMBL",
         multiVals = "first")
}


prefilter_counts <- function(counts_all, zero_percentage=0.7, entrez=TRUE ){
  
  counts_f <- counts_all
  
  if(entrez){
    genes <- data.frame(ensembl = rownames(counts_all) )
    genes$entrez <- get_entrez(genes$ensembl)
    counts_f1 <- filter(genes, !is.na(entrez))
    counts_f <- counts_all[counts_f1$ensembl,]
  }

  # Set the threshold (70%)
  threshold <- zero_percentage
  
  # Calculate the proportion of less than 5 in each row
  row_zero_proportion <- rowSums(counts_f < 5) / ncol(counts_f)
  
  # Filter the rows where the proportion of lowcounts is less than or equal to 70%
  counts_f <- counts_f[row_zero_proportion <= threshold, ]
  
  #return
  counts_f
}

# Calculate R-squared
r_squared <- function(data, predicted_data) {
  ss_res <- sum((data - predicted_data)^2)  # Residual sum of squares
  ss_tot <- sum((data - mean(data))^2)  # Total sum of squares
  return(1 - (ss_res / ss_tot))
}

# Calculate RMSE
rmse <- function(data, predicted_data) {
  sqrt(mean((data - predicted_data)^2))  # Root of Mean Squared Error
}

predicted_values <- function( all_counts, info_genes, samples_healthy, test_samples, t = "none", dataset=FALSE, interact = FALSE ){
  selected_counts <- info_genes
  counts_filt <- all_counts[ selected_counts$genes , rownames(samples_healthy) ]
  
  X <-as.matrix(counts_filt)
  Xt <- t(X)
  Xt <- scale(Xt)
  G <- tcrossprod(Xt ) /ncol(Xt)
  X_G=model.matrix(~0+gender, data=samples_healthy)
  X_I= model.matrix(~0+X_G:Xt)
  

  
  if(interact){
    ETA <- list(Env1= list(X=X_G, model="FIXED"),
                GenInt= list(X=X_I, model="BayesB"),
                Gen=list(K = G, model = 'RKHS')
    )
  }else if(dataset){
    X_d= model.matrix(~0+dataset, data=samples_healthy)
    K_env <- tcrossprod(scale(X_d))/ncol(X_d)
    ETA <- list(Env1= list(X=X_G, model="FIXED"),
                 Env5= list(K=K_env , model="RKHS"),
                Gen=list(K = G, model = 'RKHS')
    )
  }else{
    ETA <- list(Env1= list(X=X_G, model="FIXED"),
                Gen=list(K = G, model = 'RKHS')
    )
  }

  yNA <- samples_healthy$age
  
  if(t == "sqrt"){
    yNA <- sqrt(samples_healthy$age )
  }
  
  if(t == "log"){
    yNA <- log(samples_healthy$age )
  }
  
  
  yNA[test_samples] <- NA
  
  gblup_bayes <- BGLR(y = yNA, ETA = ETA, nIter = 12000, burnIn = 2000, verbose = FALSE)
  predicted_val <- (gblup_bayes$yHat)
  
  if(t == "sqrt"){
    predicted_val <- (gblup_bayes$yHat)^2  
  }
  
  if(t == "log"){
    predicted_val <- exp(gblup_bayes$yHat) 
  }
  predicted_val
}

compute_metrics <- function(predicted_matrix, true_age ,test_samples, loo=FALSE){
  
  n_models <- ncol(predicted_matrix)
  
  if(loo == FALSE){
    metrics_matrix <- matrix( nrow=n_models , ncol=10)
    col_names <- c("test_r2",
                  "train_r2",
                  "test_cor",
                  "train_cor",
                  "test_mae",
                  "train_mae",
                  "test_med",
                  "train_med",
                  "tst_rmse",
                  "train_rmse")
    
    tst <- test_samples
    for(i in 1:n_models){
      predicted_val <- predicted_matrix[,i]
      metrics <- c(  
        r_squared(samples_healthy$age[tst], predicted_val[tst] ),
        r_squared(samples_healthy$age[-tst], predicted_val[-tst] ),
        cor(samples_healthy$age[tst], predicted_val[tst]),
        cor(samples_healthy$age[-tst], predicted_val[-tst] ),
        mean(abs(samples_healthy$age[tst]-predicted_val[tst] )),
        mean(abs(samples_healthy$age[-tst]-predicted_val[-tst] )),
        median(abs(samples_healthy$age[tst]-predicted_val[tst] )),
        median(abs(samples_healthy$age[-tst]-predicted_val[-tst] )),
        rmse(samples_healthy$age[tst], predicted_val[tst] ),
        rmse(samples_healthy$age[-tst],predicted_val[-tst] )
      )
      metrics_matrix[i,] <- metrics 
      
    }
    colnames(metrics_matrix) <- col_names
    as.data.frame(metrics_matrix)
  }else{
    metrics_matrix <- matrix( nrow=n_models , ncol=10)
    col_names <- c("error",
                  "test_index",
                  "train_r2",
                  "train_cor",
                  "test_mae",
                  "train_mae",
                  "test_med",
                  "train_med",
                  "tst_rmse",
                  "train_rmse")
    
    tst <- test_samples
    
    for(i in 1:n_models){
      predicted_val <- predicted_matrix[,i]
      metrics <- c(  
        (true_age[tst[1] ] - predicted_val[tst[1] ]  ),
        tst[1],
        r_squared(samples_healthy$age[-tst], predicted_val[-tst] ),
        cor(samples_healthy$age[-tst], predicted_val[-tst] ),
        mean(abs(samples_healthy$age[tst]-predicted_val[tst] )),
        mean(abs(samples_healthy$age[-tst]-predicted_val[-tst] )),
        median(abs(samples_healthy$age[tst]-predicted_val[tst] )),
        median(abs(samples_healthy$age[-tst]-predicted_val[-tst] )),
        rmse(samples_healthy$age[tst], predicted_val[tst] ),
        rmse(samples_healthy$age[-tst],predicted_val[-tst] )
      )
      metrics_matrix[i,] <- metrics 
    }
    colnames(metrics_matrix) <- col_names
    as.data.frame(metrics_matrix)
  }
}

process_gene <- function(gene, counts, samples, form) {
  gene_data <- as.numeric(counts[gene, ])
  samples$gene_expression <- gene_data
  
  fit <- lm(form, data = samples)
  coef_summary <- summary(fit)$coefficients
  
  if (!"gene_expression" %in% rownames(coef_summary)) {
    return(NULL)  
  }
  
  pv <- coef_summary["gene_expression", "Pr(>|t|)"]
  preds <- fit$fitted.values
  
  data.frame(
    id = gene,
    pv = pv,
    r2 = r_squared(samples$age, preds),
    rmse = rmse(samples$age, preds),
    cor = cor(samples$age, preds),
    stringsAsFactors = FALSE
  )
}
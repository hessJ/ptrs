# -- masthead
masthead = paste("**********************************
Polytranscript risk scoring (PTRS)
Jonathan L. Hess, PhD
SUNY Upstate Medical University
Syracuse, NY, USA
e-mail: hessjo@upstate.edu
**********************************")

# print masthead
cat(masthead)

# --- load packages
suppressMessages(require(plyr))
suppressMessages(require(data.table))

# -- Load gene weights for a phentype and format columns
load_weights = function(file, gene_col, weight_col, p_col, sep="\t", header = T, calc_weight = FALSE, beta_col = '', se_col = ''){
  # -- check parameters exist
  if(is.null(file) == TRUE){stop("Weight file required!")}
  if(is.null(gene_col) == TRUE){stop("Column denoting gene name required!")}
  if(is.null(p_col) == TRUE){stop("Column denoting p-value required!")}
  
  # -- load file
  load_file =  fread(file, header = header, sep=sep)
  
  weight_file <<- load_file
  
  if(calc_weight == FALSE){
    if(is.null(weight_col) == TRUE){stop("Column denoting weight required!")}
  # -- keep columns required for scoring
  munge_file <<- load_file[,colnames(load_file) %in% c(gene_col, weight_col, p_col),with=FALSE]
  
  # -- rename columns
  names(munge_file)[names(munge_file) %in% gene_col] = "gene_id"
  names(munge_file)[names(munge_file) %in% weight_col] = "weight"
  names(munge_file)[names(munge_file) %in% p_col] = "pvalue"
  
  } 
  
  if(calc_weight == TRUE){
    
    if(beta_col == ""){stop("Warning! Column with beta coefficient required to calculate weight.")}
    if(se_col == ""){stop("Warning! Column with standard error required to calculate weight.")}
    
    load_file$weight = load_file[,colnames(load_file) %in% beta_col,with=FALSE]/load_file[,colnames(load_file) %in% se_col,with=FALSE]
    
    munge_file <<- load_file[,colnames(load_file) %in% c(gene_col, "weight", p_col),with=FALSE]
    
    # -- rename columns
    names(munge_file)[names(munge_file) %in% gene_col] = "gene_id"
    names(munge_file)[names(munge_file) %in% p_col] = "pvalue"
    
  }
  
  # -- object available to environment 
  return(munge_file[order(munge_file$pvalue, decreasing = FALSE), ])
}

# Function to calculate risk scores using PTRS algorithm
ptrs = function(dat = NULL, weight_table = NULL, p_thres = c(0.001, 0.01, 0.05, 0.1, 0.5, 1.0), scale_scores = TRUE){
  # -- check existance of parameters
  if(is.null(dat) == TRUE){stop("Data frame with gene expression data required!")}
  if(is.null(weight_table) == FALSE) {
    ptrs_weights = weight_table
  } else {stop("Warning! A data frame with gene weights is required.")}
  
  # -- format as data frame
  dat = data.frame(dat)
  
  # -- check that rownames match gene_id column in ptrs_weights
  common_gene_id = intersect(rownames(dat), ptrs_weights$gene_id)
  if(length(common_gene_id) < 1){stop("No common gene IDs were found! Please check that gene IDs are in same style between target data set and weights file.")}
  
  # -- reorder dat according to munge_file
  dat = dat[rownames(dat) %in% common_gene_id, ]
  munge_file = ptrs_weights[ptrs_weights$gene_id %in% common_gene_id, ]
  
  if(nrow(munge_file) != nrow(dat)){stop("Unexpected number of rows. Exiting.")}
  dat = dat[match(munge_file$gene_id, rownames(dat)), ]
  
  # -- score across thresholds
  scores = list()
  for(x in 1:length(p_thres)){
    message("\rScoring threshold: ", p_thres[[x]])
    weight_sub = munge_file[munge_file$pvalue <= p_thres[[x]], ]
    dat_sub = dat[rownames(dat) %in% weight_sub$gene_id, ]
    dat_sub = data.frame(t(dat_sub))
    scoreMatrix = sweep(dat_sub, MARGIN = 2, STATS = weight_sub$weight, `*`)
    if(scale_scores == TRUE){
    scores[[x]] = scale(rowSums(scoreMatrix, na.rm=TRUE))} else {
      scores[[x]] = rowSums(scoreMatrix, na.rm=TRUE)
    }
  }
  score_df = data.frame(do.call(cbind, scores))
  colnames(score_df) = paste("ptrs_",1:length(p_thres), sep="")
  return(score_df)
}


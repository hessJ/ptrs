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

# -- functions
load_weights = function(file, gene_col, weight_col, p_col, sep="\t", header = T){
  # -- check parameters exist
  if(is.null(file) == TRUE){stop("Weight file required!")}
  if(is.null(gene_col) == TRUE){stop("Column denoting gene name required!")}
  if(is.null(weight_col) == TRUE){stop("Column denoting weight required!")}
  if(is.null(p_col) == TRUE){stop("Column denoting p-value required!")}
  
  # -- load file
  load_file =  fread(file, header = header, sep=sep)
  
  # -- keep columns required for scoring
  munge_file <<- load_file[,colnames(load_file) %in% c(gene_col, weight_col, p_col),with=FALSE]
  
  # -- rename columns
  names(munge_file)[names(munge_file) %in% gene_col] = "gene_id"
  names(munge_file)[names(munge_file) %in% weight_col] = "weight"
  names(munge_file)[names(munge_file) %in% p_col] = "pvalue"
  
  # -- object available to environment 
  munge_file <<- munge_file[order(munge_file$pvalue, decreasing = FALSE), ]
}

ptrs = function(dat, p_thres = c(0.001, 0.01, 0.05, 0.1, 0.5, 1.0), scale_scores = TRUE){
  # -- check existance of parameters
  if(is.null(dat) == TRUE){stop("Data frame with gene expression data required!")}
  
  # -- format as data frame
  dat = data.frame(dat)
  
  # -- check that rownames match gene_id column in munge_file
  common_gene_id = intersect(rownames(dat), munge_file$gene_id)
  if(length(common_gene_id) < 1){stop("No common gene IDs were found! Please check that gene IDs are in same style between target data set and weights file.")}
  
  # -- reorder dat according to munge_file
  dat = dat[rownames(dat) %in% common_gene_id, ]
  munge_file = munge_file[munge_file$gene_id %in% common_gene_id, ]
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
  colnames(score_df) = paste("ptrs.",p_thres, sep="")
  return(score_df)
}


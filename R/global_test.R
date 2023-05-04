#' Fetches global test statistic
#'
#' Tabulates the global test statistic according to Hoetelling's T2
#'
#' @param samples_used An array of CLR samples
#' @param conditions A vector of conditions
#' @return A single test statistic and p-value.
#'
#'
#' @export
global_test <- function(matrix_samples, conditions){
  mc.samples <- dim(matrix_samples)[3]
  if(is.matrix(conditions)){
    mmX <- conditions
    if(sum(mmX[,1]==1) == nrow(mmX)){
      mmX <- mmX[,-1]
    }
    mmX <- as.data.frame(mmX)
    q <- ncol(mmX)
    f_stats <- matrix(NA, nrow = q, ncol = mc.samples)
    p_vals <- matrix(NA, nrow = q, ncol = mc.samples)
  } else{
    q <- length(unique(conditions))
    mmX <- as.data.frame(model.matrix(~conditions -1))
    f_stats <- matrix(NA, nrow = q-1, ncol = mc.samples)
    p_vals <- matrix(NA, nrow = q-1, ncol = mc.samples)
  }
  
  ## adding entity information
  
  form <- as.formula(paste0("mc.s~", paste(paste(names(mmX), "entity", sep="*"), collapse = "+")))
  
  
  ##Expanding the design matrix
  mmX_expand <- matrix(NA, nrow = dim(matrix_samples)[1]*dim(matrix_samples)[2], ncol = q+2)
  for(j in 1:q){
    mmX_expand[,j] <- rep(mmX[,j], each = dim(matrix_samples)[1])
  }
  mmX_expand[,(q+1)] <- rep(paste0("entity", 1:dim(matrix_samples)[1]), dim(matrix_samples)[2])
  
  message("Caution: this may take awhile if the number of Monte Carlo samples is large.")
  for(i in 1:mc.samples){
    mc.s <- vector()
    for (j in 1:dim(matrix_samples)[2]){
      mc.s <- c(mc.s, matrix_samples[, j,i])
    }    
    mmX_expand[,(q+2)] <- mc.s
    colnames(mmX_expand) <- c(colnames(mmX), "entity", "mc.s")
    mmX_expand <- data.frame(mmX_expand)
    mod <- summary(aov(formula = form,  data = mmX_expand))[[1]]
    
    ## selecting rows and averaging
    inds <- which(str_detect(rownames(mod),":entity|entity:"))
    f_stats[,i] <- mod[inds,4]
    p_vals[,i] <- mod[inds,5]
  }
  
  cond_names <- str_trim(noquote(sub(":entity|entity:", "", rownames(mod)[inds])))
  return(list(conds = cond_names, f_stat = rowMeans(f_stats), p_val = rowMeans(p_vals)))
}

samples_mat <- function(samples){
  D <- nrow(samples[[1]])
  n_samples <- ncol(samples[[1]])
  N <- length(samples)
  lambda_mat <- array(as.numeric(unlist(samples)), dim=c(D, n_samples, N))
  lambda_mat <- aperm(lambda_mat, c(1,3,2))
  return(lambda_mat)
}

gm <- function(x, na.rm = TRUE){
  exp(mean(log(x[x > 0]), na.rm=na.rm))
}

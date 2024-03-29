#' Calculate glm test statistics using a \code{model.matrix}
#'
#' \code{aldex.glm} calculates the expected values for each coefficient of a
#'  glm model on the data returned by \code{aldex.clr}. This function
#'  requires the user to define a model with \code{model.matrix}.
#'
#' @param clr An \code{ALDEx2} object. The output of \code{aldex.clr}.
#' @inheritParams fetch
#' @param ... Arguments passed to \code{glm}.
#'
#' @return Returns a data.frame of the average
#'  coefficients and their p-values for each feature,
#'  with FDR appended as a \code{BH} column.
#'
#' @author Thom Quinn
#'
#' @seealso
#'  \code{\link{fetch}},
#'  \code{\link{aldex.clr}},
#'  \code{\link{aldex.ttest}},
#'  \code{\link{aldex.kw}},
#'  \code{\link{aldex.glm}},
#'  \code{\link{aldex.effect}},
#'  \code{\link{aldex.corr}},
#'  \code{\link{selex}}
#'
#'
#' @examples
#' data(selex)
#' #subset for efficiency
#' selex <- selex[1201:1600,]
#' covariates <- data.frame("A" = sample(0:1, 14, replace = TRUE),
#'                          "B" = c(rep(0, 7), rep(1, 7)))
#' mm <- model.matrix(~ A + B, covariates)
#' x <- aldex.clr(selex, mm, mc.samples=1, denom="all")
#' glm.test <- aldex.glm(x)
aldex.glm <- function(clr, verbose=FALSE, ...){

  # Use clr conditions slot instead of input
  conditions <- clr@conds

  lr2glm <- function(lr, conditions, ...){

    if( !is(conditions, "matrix") &
       !("assign" %in% names(attributes(conditions)))){

      stop("Please define the aldex.clr object for a model.matrix 'conditions'.")
     }

    if(nrow(lr) != nrow(conditions)){

      stop("Input data and 'model.matrix' should have same number of rows.")
    }

    # Build the glm models
    model. <- conditions
    glms <- apply(lr, 2, function(x){
      glm(x ~ model., ...)
    })

    # Extract coefficients and p-values
    extract <- function(model){
      x <- coef(summary(model))
      coefs <- lapply(1:nrow(x), function(i){
        y <- x[i,,drop=FALSE]
        colnames(y) <- paste(rownames(y), colnames(y))
        y})
      do.call("cbind", coefs)
    }

    # Combine to make data.frame
    extracts <- lapply(glms, extract)
    df <- do.call("rbind", extracts)
    rownames(df) <- colnames(lr)
    df <- as.data.frame(df)

    # Create new data.frame for FDR
    pvals <- colnames(df)[grepl("Pr\\(>", colnames(df))]
    df.bh <- df[,pvals]
    colnames(df.bh) <- paste0(colnames(df.bh), ".BH")
    for(j in 1:ncol(df.bh)){
      df.bh[,j] <- p.adjust(df.bh[,j])
    }

    # Merge results with FDR
    cbind(df, df.bh)
  }

  # Keep a running sum of lr2glm instances
  # verbose was throwing an error 'the condition has length > 1'
  if(verbose[1] == TRUE) message("running tests for each MC instance:")
  mc <- fetch::getMonteCarloInstances(clr)
  k <- fetch::numMCInstances(clr)
  r <- 0
  for(i in 1:k){

    if(verbose[1] == TRUE ){ numTicks <- progress(i, k, numTicks) }
    mci_lr <- t(sapply(mc, function(x) x[, i]))
    r <- r + lr2glm(mci_lr, conditions, ...)
  }

  r / k # return expected
}

#' Calculate Wilcoxon Rank Sum test and Welch's t-test statistics
#'
#' \code{aldex.ttest} calculates the expected values of the Wilcoxon Rank Sum
#'  test and Welch's t-test on the data returned by \code{aldex.clr}.
#'
#' @param clr An \code{ALDEx2} object. The output of \code{aldex.clr}.
#' @inheritParams fetch
#' @param paired.test Toggles whether to calculate paired tests.
#' @param hist.plot Toggles whether to plot a histogram of p-values for the
#' first Dirichlet Monte Carlo instance.
#' @param verbose Boolean. If TRUE, outputs extra notes during computation.
#'
#' @return Returns a \code{data.frame} with the following information:
#' \item{we.ep}{ a vector containing the expected p-value of Welch's t-test
#'  for each feature }
#' \item{we.eBH}{ a vector containing the corresponding expected value of the
#'  Benjamini-Hochberg corrected p-value for each feature }
#' \item{wi.ep}{ a vector containing the expected p-value of the Wilcoxon Rank Sum test
#'  for each feature }
#' \item{wi.eBH}{ a vector containing the corresponding expected value of the
#'  Benjamini-Hochberg corrected p-value for each feature }
#'
#' @author Greg Gloor
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
#' @examples
#' data(selex)
#' #subset for efficiency
#' selex <- selex[1201:1600,]
#' conds <- c(rep("NS", 7), rep("S", 7))
#' x <- aldex.clr(selex, conds, mc.samples=2, denom="all")
#' ttest.test <- aldex.ttest(x)
aldex.ttest <- function(clr, paired.test=FALSE, hist.plot=FALSE, verbose=FALSE) {
  
  # Use clr conditions slot instead of input
  conditions <- clr@conds
  if(length(unique(conditions)) != 2){
    stop("Please define the aldex.clr object for a vector of two unique 'conditions'.")
  }
  
  # get dimensions, names, etc from the input data
  smpl.ids <- getSampleIDs(clr)
  feature.names <- getFeatureNames(clr)
  feature.number <- numFeatures(clr)
  mc.instances <- numMCInstances(clr)
  
  conditions <- as.factor( conditions )
  levels     <- levels( conditions )
  
  if ( length( conditions ) !=  numConditions(clr) ){
    stop(paste("mismatch btw 'length(conditions)' and 'length(names(clr))'. len(condtitions):",
               length(conditions),"len(names(clr)):",numConditions(clr)))}
  if ( length( levels ) != 2 ) stop("only two condition levels are currently supported")
  
  # generate the comparison sets from the condition levels
  levels <- vector( "list", length( levels ) )
  names( levels ) <- levels( conditions )
  sets <- names(levels)
  setAsBinary <- as.numeric(conditions == sets[1])
  setA <- which(conditions == sets[1])
  setB <- which(conditions == sets[2])
  
  # set up the t-test result containers
  wi.p.matrix <- as.data.frame(matrix(1, nrow = feature.number, ncol = mc.instances))
  wi.BH.matrix <- wi.p.matrix # duplicate result container
  we.p.matrix <- wi.p.matrix # duplicate result container
  we.BH.matrix <- wi.p.matrix # duplicate result container
  
  # mc.i is the i-th Monte-Carlo instance
  if(verbose) message("running tests for each MC instance:")
  mc.all <- getMonteCarloInstances(clr)
  for(mc.i in 1:mc.instances){
    
    if(verbose) numTicks <- progress(mc.i, mc.instances, numTicks)
    
    # generate a matrix of i-th Monte-Carlo instance, columns are samples, rows are features
    t.input <- sapply(mc.all, function(y){y[, mc.i]})
    
    wi.p.matrix[, mc.i] <- wilcox.fast(t.input, setAsBinary, paired.test)
    wi.BH.matrix[, mc.i] <- p.adjust(wi.p.matrix[, mc.i], method = "BH")
    
    we.p.matrix[, mc.i] <- t.fast(t.input, setAsBinary, paired.test)
    we.BH.matrix[, mc.i] <- p.adjust(we.p.matrix[, mc.i], method = "BH")
  }
  
  if(hist.plot == TRUE){
    
    par(mfrow=c(2,2))
    hist(we.p.matrix[,1], breaks=99, main="Welch's P values Instance 1")
    hist(wi.p.matrix[,1], breaks=99, main="Wilcoxon P values Instance 1")
    hist(we.BH.matrix[,1], breaks=99, main="Welch's BH values Instance 1")
    hist(wi.BH.matrix[,1], breaks=99, main="Wilcoxon BH values Instance 1")
    par(mfrow=c(1,1))
  }
  
  # get the Expected values of p, q and lfdr
  we.ep <- rowMeans(we.p.matrix) # rowMeans is faster than apply()!!
  we.eBH <- rowMeans(we.BH.matrix)
  wi.ep <- rowMeans(wi.p.matrix)
  wi.eBH <- rowMeans(wi.BH.matrix)
  
  z <- data.frame(we.ep, we.eBH, wi.ep, wi.eBH)
  rownames(z) <- getFeatureNames(clr)
  return(z)
}

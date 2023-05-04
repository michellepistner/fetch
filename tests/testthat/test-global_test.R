test_that("globalFetch", {
  library(plyr)
  library(dplyr)
  library(fetch)
  metadata_filtered = mouse_metadata %>%
    as.data.frame() %>%
    filter(!is.na(P_C_F))

  samples_to_keep = metadata_filtered$X.SampleID

  mouse_counts = as.data.frame(mouse_counts)
  counts_subset = mouse_counts[,names(mouse_counts) %in% samples_to_keep]
  counts_subset = counts_subset[,colSums(counts_subset)>=5000]
  rows_to_keep = rowSums(mouse_counts >= 10) >= 10
  other_sum = colSums(counts_subset[!(rows_to_keep),])
  otu_filtered = rbind(counts_subset[rows_to_keep,], other_sum)

  samples_to_keep = names(otu_filtered)
  metadata_filtered = metadata_filtered %>%
    filter(X.SampleID %in% samples_to_keep)

  ###Prepping the data
  conds = matrix(metadata_filtered$P_C_F, nrow = 1)

  set.seed(1)
  out = fetch(otu_filtered, c(conds), denom = "all", test = "global")

  expect_equal(round(c(out$f_stat),3), round(c(0.99592),3))
})


test_that("localFetch", {
  set.seed(1)
  Y = rmultinom(10, 10000, prob = seq(1:100))
  X = matrix(sample(1:2,10,replace=TRUE),ncol = 10)
  out =fetch(Y, c(X), denom = "iqlr", test = "t")

  expect_equal(out[1,1], -5.000539551)
})

test_that("localFetchvsAldex", {
  set.seed(1)
  Y = rmultinom(10, 10000, prob = seq(1:15))
  X = matrix(sample(1:2,10,replace=TRUE),ncol = 10)
  out =fetch(Y, c(X), denom = "all", test = "t")

  out2 = ALDEx2::aldex(Y,c(X),bayesEst = FALSE)
  sig2 = ifelse(out2[,11] > 0.05,0,1)
  sig = ifelse(out[,11] > 0.05,0,1)

  expect_equal(sig, sig2)
})


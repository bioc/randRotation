test_that("Test df.estimation", {
  set.seed(0)
  # Dataframe of phenotype data (sample information)
  pdata = data.frame(batch = rep(1:3, c(10,10,10)),
                     phenotype = rep(c("Control", "Cancer"), c(5,5)))
  features = 100

  # Matrix with random gene expression data
  edata = matrix(rnorm(features * nrow(pdata)), features)
  rownames(edata) = paste("feature", 1:nrow(edata))
  colnames(edata) = paste("sample", 1:ncol(edata))

  mod1 = model.matrix(~phenotype, pdata)



  mapping <- function(Y, batch, mod){
    capture.output(suppressMessages(Y.tmp <- sva::ComBat(dat = Y, batch = batch,
                                                         mod = mod)
    ))
    Y.tmp
  }

  dfs <- df_estimate(
    data = edata,
    features = 1,
    mapping = mapping,
    batch = pdata$batch, mod = mod1)

  expect_equal(median(dfs), 27.6, tolerance = 1)



  mapping <- function(Y, batch, mod){
    limma::removeBatchEffect(x = Y, batch = batch, design = mod)
  }

  dfs <- df_estimate(
    data = edata,
    features = 1,
    mapping = mapping,
    batch = pdata$batch, mod = mod1)

  expect_equal(median(dfs), 28)



  mapping <- function(Y){
    abs(Y)
  }

  dfs <- df_estimate(
    data = edata,
    features = 1,
    mapping = mapping)

  expect_equal(median(dfs), 30)



  mapping <- function(Y){
    abs(Y-100)
  }

  dfs <- df_estimate(
    data = edata,
    features = 1,
    mapping = mapping)

  expect_equal(median(dfs), 30)



  mapping <- function(Y){
    floor(Y)
  }

  dfs <- df_estimate(
    data = edata,
    features = 1,
    mapping = mapping)

  expect_equal(median(dfs), 0)



  mapping <- function(Y, A){
    Y %*% A
  }

  dfs <- df_estimate(
    data = edata,
    features = 1,
    mapping = mapping,
    A = randpermut(ncol(edata)))

  expect_equal(median(dfs), 30)



  mapping <- function(Y){
    -Y
  }

  dfs <- df_estimate(
    data = edata,
    features = 1,
    mapping = mapping)

  expect_equal(median(dfs), 30)



  mapping <- function(Y){
    matrix(0, nrow(Y), ncol(Y))
  }

  dfs <- df_estimate(
    data = edata,
    features = 1,
    mapping = mapping)

  expect_equal(median(dfs), 0)


  })

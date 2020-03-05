test_that("Test execution of rotate.stat (including parallel computation mode), p.fdr with \"fdr.q\" and \"fdr.qu\"", {
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

  capture.output(init1 <- init.batch.randrot(edata, mod1, -1, pdata$batch))


  statistic <- function(Y, batch, mod, coef){
    capture.output(suppressMessages(Y.tmp <- sva::ComBat(Y, batch = batch, mod)))

    fit1 <- lm.fit(mod, t(Y.tmp))
    t(abs(coef(fit1)[coef,, drop = FALSE]))
  }

  res1 <- rotate.stat(initialised.obj = init1,
                      R = 100,
                      statistic = statistic,
                      batch = pdata$batch, mod = mod1, coef = 1:2,
                      parallel = FALSE)


  expect_equal(sort(p.fdr(res1))[1:10], c(0.0047,0.0096,0.0143,0.0198,0.0286,0.0293,0.0324,0.0366,0.0429,0.0474), tolerance = 1e-7)
  expect_equal(sort(p.fdr(res1, "fdr.q"))[1:10], c(0.4700000,0.9600000,0.9715805,0.9936457,0.9996935,0.9999904,0.9999931,0.9999982,0.9999998,0.9999999), tolerance = 1e-7)
  expect_equal(sort(p.fdr(res1, "fdr.qu"))[1:10], c(0.4700000,0.9600000,0.9715805,0.9936457,0.9996935,0.9999904,0.9999931,0.9999982,0.9999998,0.9999999), tolerance = 1e-7)


  ### Test for exact dimensions and for parallel processing
  res1 <- rotate.stat(initialised.obj = init1, R = 100, statistic = statistic,
                      batch = pdata$batch, mod = mod1, coef = 1:2,
                      parallel = TRUE)
  expect_equal(dim(p.fdr(res1)), c(features,2))


  #### If "statistic" returns only 1 number, not one for each feature --> this still has to work properly !

  statistic <- function(Y, batch, mod, coef){
    capture.output(suppressMessages(Y.tmp <- sva::ComBat(Y, batch = batch, mod)))

    fit1 <- lm.fit(mod, t(Y.tmp))
    t(abs(coef(fit1)[coef,1, drop = FALSE]))
  }

  res1 <- rotate.stat(initialised.obj = init1, R = 10, statistic = statistic,
                      batch = pdata$batch, mod = mod1, coef = 1:2,
                      parallel = FALSE)
  expect_equal(dim(p.fdr(res1)), c(1,2))

  res1 <- rotate.stat(initialised.obj = init1, R = 10, statistic = statistic,
                      batch = pdata$batch, mod = mod1, coef = 2,
                      parallel = FALSE)
  expect_equal(dim(p.fdr(res1)), c(1,1))


  #### If "statistic" returns only 1 number, not one for each feature (this still has to work properly !)
  #### Additionally --> R = 1

  statistic <- function(Y, batch, mod, coef){
    capture.output(suppressMessages(Y.tmp <- sva::ComBat(Y, batch = batch, mod)))

    fit1 <- lm.fit(mod, t(Y.tmp))
    t(abs(coef(fit1)[coef,1, drop = FALSE]))
  }

  res1 <- rotate.stat(initialised.obj = init1, R = 10, statistic = statistic,
                      batch = pdata$batch, mod = mod1, coef = 2,
                      parallel = FALSE)

  expect_equal(dim(p.fdr(res1)), c(1,1))

  })

test_that("Test execution of rotateStat (including parallel computation mode), pFdr with \"fdr.q\" and \"fdr.qu\"", {
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

  capture.output(init1 <- initBatchRandrot(edata, mod1, -1, pdata$batch))


  statistic <- function(Y, batch, mod, coef){
    capture.output(suppressMessages(Y.tmp <- sva::ComBat(Y, batch = batch, mod)))

    fit1 <- lm.fit(mod, t(Y.tmp))
    t(abs(coef(fit1)[coef,]))
  }

  res1 <- rotateStat(initialised.obj = init1,
                      R = 100,
                      statistic = statistic,
                      batch = pdata$batch, mod = mod1, coef = 1:2,
                      parallel = FALSE)


  expect_equal(sort(pFdr(res1))[1:10], c(0.00529947, 0.00919908, 0.01529847, 0.01999800, 0.02839716, 0.02939706, 0.03169683, 0.03689631, 0.04279572, 0.04799520), tolerance = 1e-7)
  expect_equal(sort(pFdr(res1, "fdr.q"))[1:10], c(0.5200000, 0.9100000, 0.9742616, 0.9942451, 0.9997224, 0.9999845, 0.9999913, 0.9999984, 0.9999995, 0.9999999), tolerance = 1e-7)
  expect_equal(sort(pFdr(res1, "fdr.qu"))[1:10], c(0.5200000, 0.9100000, 0.9742616, 0.9942451, 0.9997224, 0.9999845, 0.9999913, 0.9999984, 0.9999995, 0.9999999), tolerance = 1e-7)


  ### Test for exact dimensions and for parallel processing
  res1 <- rotateStat(initialised.obj = init1, R = 10, statistic = statistic,
                      batch = pdata$batch, mod = mod1, coef = 1:2,
                      parallel = TRUE)
  expect_equal(dim(pFdr(res1)), c(features,2))


  #### If "statistic" returns only 1 number, not one for each feature --> this still has to work properly !

  statistic <- function(Y, batch, mod, coef){
    capture.output(suppressMessages(Y.tmp <- sva::ComBat(Y, batch = batch, mod)))

    fit1 <- lm.fit(mod, t(Y.tmp))
    t(abs(coef(fit1)[coef,1]))
  }

  res1 <- rotateStat(initialised.obj = init1, R = 10, statistic = statistic,
                      batch = pdata$batch, mod = mod1, coef = 1:2,
                      parallel = FALSE)
  expect_equal(dim(pFdr(res1)), c(1,2))

  res1 <- rotateStat(initialised.obj = init1, R = 10, statistic = statistic,
                      batch = pdata$batch, mod = mod1, coef = 2,
                      parallel = FALSE)
  expect_equal(dim(pFdr(res1)), c(1,1))


  #### If "statistic" returns only 1 number, not one for each feature (this still has to work properly !)
  #### Additionally --> R = 1

  statistic <- function(Y, batch, mod, coef){
    capture.output(suppressMessages(Y.tmp <- sva::ComBat(Y, batch = batch, mod)))

    fit1 <- lm.fit(mod, t(Y.tmp))
    t(abs(coef(fit1)[coef,1]))
  }

  res1 <- rotateStat(initialised.obj = init1, R = 1, statistic = statistic,
                      batch = pdata$batch, mod = mod1, coef = 2,
                      parallel = FALSE)

  expect_equal(dim(pFdr(res1)), c(1,1))

  })

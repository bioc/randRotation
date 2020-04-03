test_that("Test df.estimation and idempot", {
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

  capture.output(init1 <- init.batch.randrot(Y = edata, X = mod1, coef.h = 2,
                                             batch = pdata$batch))

  mapping <- function(Y, batch, mod){
    capture.output(suppressMessages(Y.tmp <- sva::ComBat(dat = Y, batch = batch,
                                                         mod = mod)
    ))

    Y.tmp
  }

  idempot(edata, mapping = mapping, batch = pdata$batch, mod = mod1)

  dfs <- df.estimate(initialised.obj = init1,
              mapping = mapping,
              batch = pdata$batch, mod = mod1,
              R = 100)

  expect_equal(median(dfs), 29, tolerance = 1)

  })

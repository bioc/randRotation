test_that("Random rotation gives same coefficients if all coefs are set as coef.d ?", {
  set.seed(0)
  # Dataframe of phenotype data (sample information)
  pdata = data.frame(phenotype = rep(c("Control", "Cancer"), c(5,5)))
  features = 10

  # Matrix with random gene expression data
  edata = matrix(rnorm(features * nrow(pdata)), features)
  rownames(edata) = paste("feature", 1:nrow(edata))

  mod1 = model.matrix(~phenotype, pdata)

  rr1 = init.randrot(edata, mod1, -1)


  fit1 = lm.fit(mod1, t(edata))

  edata.rotated = randrot(rr1)
  fit.rotated = lm.fit(mod1, t(edata.rotated))

  expect_equal(coef(fit1), coef(fit.rotated))



  })


test_that("Model transformation by function contrastModel must result in equal contrast estimates.", {
  set.seed(0)
  # Dataframe of group data (sample information)
  pdata = data.frame(group = rep(c("1", "2", "3"), c(5,5,5)))
  features = 1

  # Matrix with random gene expression data
  edata <- matrix(rnorm(features * nrow(pdata)), nrow = features)
  rownames(edata) <- paste0("feature", 1:nrow(edata))

  mod1 <- model.matrix(~0+group, pdata)

  # weights
  w <- runif(nrow(pdata))


  # ordinary contrast estimates
  fit1 <- lm.wfit(x = mod1, y = t(edata), w = w)
  cont1 <- cbind(c1 = c(1,0,-1), c2 = c(0,-1,1))
  gamma1 <- drop(t(cont1) %*% coef(fit1))


  # transformed model
  mod.cont <- contrastModel(mod1, cont1)
  fit2 <- lm.wfit(x = mod.cont, y = t(edata), w = w)
  gamma2 <- coef(fit2)[attr(mod.cont, "coef.h")]

  expect_equal(gamma1, gamma2)

  })


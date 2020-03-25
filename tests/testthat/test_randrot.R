test_that("Random rotation gives same coefficients if all coefs are set as coef.d ?", {
  set.seed(0)
  # Dataframe of phenotype data (sample information)
  pdata = data.frame(phenotype = rep(c("Control", "Cancer"), c(5,5)))
  features = 10

  # Matrix with random gene expression data
  edata = matrix(rnorm(features * nrow(pdata)), nrow = features)
  rownames(edata) = paste0("feature", 1:nrow(edata))

  mod1 = model.matrix(~phenotype, pdata)

  rr1 = init.randrot(Y = edata, X = mod1, coef.h = -1)

  fit1 = lm.fit(x = mod1, y = t(edata))

  edata.rotated = randrot(rr1)
  fit.rotated = lm.fit(x = mod1, y = t(edata.rotated))

  expect_equal(coef(fit1), coef(fit.rotated))



  #### Test with weights specified

  weights = matrix(runif(length(edata))+0.2, nrow = nrow(edata))

  rr1 = init.randrot(Y = edata, X = mod1, coef.h = -1, weights = weights)

  coefs <- vapply(seq_len(features),
                  function(i) coef(lm.wfit(x = mod1, y = edata[i,], w = weights[i,])),
                  numeric(ncol(mod1)))

  edata.rot = randrot(rr1)

  coefs.rot <- vapply(seq_len(features),
                  function(i) coef(lm.wfit(x = mod1, y = edata.rot[i,], w = weights[i,])),
                  numeric(ncol(mod1)))

  expect_equal(coefs, coefs.rot)



  #### Test with correlation matrix and weights specified

  library(lme4)

  colnames(pdata) <- "individual"
  pdata$individual <- as.factor(rep(c("2", "1"), c(5,5)))
  df1 <- cbind(pdata, t(edata[1,,drop = FALSE]), w = weights[1,],
               phenotype = as.factor(c("Cancer.tissue", "Control.tissue")))

  ### add "individual" and "phenotype" effect to avoid numerical problems
  df1$feature1 = df1$feature1 + as.numeric(df1$individual)
  df1$feature1 = df1$feature1 + as.numeric(df1$phenotype)

  ### We take "phenotype" as hypothesis coefficient, hence it is omitted in the
  ### "coef.d only" model

  me0 <- lmer(feature1 ~ 1 + (1|individual), data = df1, weights = w)

  ### Note: The correct ("identical") parametrisation of the models is very
  ### important for data rotation !
  me1 <- lmer(feature1 ~ 1 + phenotype + (1|individual), data = df1, weights = w)
  mod1 <- as.matrix(cbind(getME(me1, "X"), getME(me1, "Z")))

  # Resort columns to avoid the following warning in "R CMD check":
  # warning("coef.h does not correspond to the last columns in the design matrix.")
  mod1 <- mod1[,c(1,3,4,2)]

  # Remove one column of the random effect to avoid the following warning in
  # "R CMD check":
  # warning("Partitioned design matrix X[,coef.d] does not have full rank.")
  mod1 <- mod1[,c(1,2,4)]

  rr1 <- init.randrot(Y = t(df1$feature1), X = mod1, coef.h = ncol(mod1), weights = t(df1$w))

  edata.rot <- randrot(rr1)

  df2 <- df1
  df2$feature1 <- edata.rot[1,]

  me0.rot <- lmer(feature1 ~ 1 + (1|individual), data = df2, weights = w)

  expect_equal(ranef(me0), ranef(me0.rot), tolerance = 1e-5)
  expect_equal(fixef(me0), fixef(me0.rot), tolerance = 1e-5)


  })


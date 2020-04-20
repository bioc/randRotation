test_that("Matrix decomposition is orthogonal ?", {
  des1 = cbind(1, rep(0:1, 5))
  d1 = X_decomp(des1)
  x1 = cbind(d1$Xd, d1$Xhe)
  expect_equal(x1 %*% t(x1), diag(ncol(x1)))
})


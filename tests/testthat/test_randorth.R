test_that("Random orthogonal matrices are orthogonal ?", {
  sapply(1:3, function(i){
    r1 = randorth(i)
    expect_equal(r1 %*% t(r1), diag(ncol(r1)))
    r1 = randorth(i, "unitary")
    expect_equal(r1 %*% Conj(t(r1)), diag(ncol(r1))+0i)
  })
})

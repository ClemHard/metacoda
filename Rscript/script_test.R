library(testthat)

source("Rscript/script1.R");

## Test norm data
test_that("norm_data detects nonpositive values", {
  expect_warning(norm_data(c(1, -1)), "Non positive values are present in the data.")
  expect_warning(norm_data(c(1, 0)), "Non positive values are present in the data.")
})

test_that("normdata returns matrices", {
  expect_equal(norm_data(c(1, 1)), matrix(c(1, 1), nrow = 1))
  expect_equal(norm_data(matrix(c(1, 1, 1, 4), nrow = 2)),
               matrix(c(1, 1, 1, 4), nrow = 2))
})

## Test closure

test_that("closure works properly", {
  expect_equal(closure(c(1, 1)), matrix(c(0.5, 0.5), nrow = 1))
  expect_equal(closure(matrix(c(1, 1, 1, 4), nrow = 2)),
               matrix(c(0.5, 0.2, 0.5, 0.8), nrow = 2))
  expect_equal(closure(c(1, 1), k = 4), matrix(c(2, 2), nrow = 1))
})

## Test Inner product

test_that("Inner product works properly", {
  expect_equal(Inner_product(matrix(c(1.2,7.1,1.9,6.2),nrow=2), matrix(c(8.1,0.3,9.2,7.5),nrow=2)), matrix(c(0.029258, -0.00863, 0.7395, -0.2181), nrow=2), tolerance=1e-3)
  expect_equal(Inner_product(c(7.1, 6.2),c(0.3, 7.5)), as.matrix(-0.21815),tolerance=1e-5)
  expect_equal(Inner_product(c(12.8, 3.2),c(18.3, 73.5)), as.matrix(-0.96374),tolerance=1e-5)

  expect_error(Inner_product(c(1.2, 1.9),c(8.1, 9.2,2.3)), "x and y should have the same dimension")
  
})


## Test pertubation

test_that("closure works properly", {
  expect_equal(perturbation(matrix(c(0.5, 0.8, 0.5, 0.2), nrow = 2), c(2, 3)), matrix(c(0.4, 0.727, 0.6, 0.272), nrow=2), tolerance=1e-2)
  expect_equal(perturbation(matrix(c(0.5, 0.8, 0.5, 0.2), nrow = 2), 3 ), matrix(c(0.5, 0.8, 0.5, 0.2), nrow = 2), tolerance=1e-2)

})

## Test power

test_that("power works properly", {
  expect_equal(power(matrix(c(0.84, 0.89, 0.53, 0.245), nrow = 2), 2.6), matrix(c(0.768, 0.9662, 0.2319, 0.03376), nrow = 2), tolerance=1e-2)
  expect_equal(power(matrix(c(0.5, 0.8, 0.5, 0.2), nrow = 2),  2), matrix(c(0.5, 0.9411765, 0.5, 0.05882353), nrow = 2), tolerance=1e-2)
  
})

## Test clr

test_that("clr works properly", {
  expect_equal(clr(matrix(c(0.84, 0.89, 0.53, 0.245 ,0.665, 0.452), nrow = 2)), matrix(c(0.23137, 0.65583, -0.2291, -0.63412, -0.00223, -0.02170), nrow = 2), tolerance=1e-2)
  expect_equal(clr(matrix(c(0.5, 0.8, 0.5, 0.2), nrow = 2)), matrix(c(0, 0.69314, 0, -0.6931), nrow = 2), tolerance=1e-2)
  
})

## Test clr inverse

test_that("clr inverse works properly", {
  expect_equal(clr_inverse(matrix(c(0.84, 0.89, 0.53, 0.245 ,0.665, 0.452), nrow = 2)), matrix(c(0.38866, 0.46083, 0.2850, 0.24178, 0.32626, 0.29738), nrow = 2), tolerance=1e-2)
  expect_equal(clr_inverse(matrix(c(0.5, 0.8, 0.5, 0.2), nrow = 2)), matrix(c(0.5, 0.6456, 0.5, 0.35434), nrow = 2), tolerance=1e-2)
  
})


## Test Base_SIGMA_matrix

test_that("Base_SIGMA_matrix inverse works properly", {
  expect_equal(Base_SIGMA_matrix(3), matrix(c(-0.70710, -0.40824, 0.7071, -0.40824, 0, 0.8164), nrow = 2), tolerance=1e-2)
  expect_equal(Base_SIGMA_matrix(5), matrix(c(-0.70710, -0.40824, -0.28867, -0.22360, 0.70710, -0.40824, -0.28867, -0.2236, 0,  0.81649, -0.28867, -0.22360, 0, 0, 0.86602, -0.2236, 0, 0, 0, 0.8944), nrow = 4), tolerance=1e-2)
  
})


## Test Base_binary_matrix

test_that("Base_SIGMA_matrix inverse works properly", {
  expect_equal(Base_binary_matrix(3), matrix(c(1, 1, -1, 1, 0, -1), nrow = 2), tolerance=1e-2)
  expect_equal(Base_binary_matrix(5), matrix(c(1, 1, 1, 1, -1, 1, 1, 1, 0,  -1, 1, 1, 0, 0, -1, 1, 0, 0, 0, -1), nrow = 4), tolerance=1e-2)
  
})


## Test SIGMA_matrix

test_that("SIGMA_matrix inverse works properly", {
  expect_equal(SIGMA_matrix(matrix(c(1, 1, 1, -1, -1, 0), nrow=2)), matrix(c(0.40824, 0.7071, 0.40824, -0.70710, -0.81649, 0), nrow = 2), tolerance=1e-2)
  expect_equal(SIGMA_matrix(matrix(c(1, 1, 1, -1, -1 ,0 ,1, -1, -1 ),nrow=3)), matrix(c(0.40824, 0.81649, 0.70710, -0.81649, -0.40824, 0, 0.40824, -0.40824, -0.70710), nrow = 3), tolerance=1e-2)
  
})

## Test ilr

test_that("ilr works properly", {
  expect_equal(ilr(matrix(c(1,2,3,4,2,1,1,6,1, 5,4,10), nrow=4)), matrix(c(0.49012, -0.49012, -0.77683, 0.28670, -0.28297, 1.03112, 0.68339, 0.58261), nrow = 4), tolerance=1e-2)
  expect_equal(ilr(matrix(c(1.2, 5.6, 4.6, 3.2, 8.2 ,1.2, 4.3 ,1.56,1, 8.45, 96, 4), nrow=4)), matrix(c(1.35892, -1.0892, -0.04768, -0.50803, -0.9334, 0.96479, 2.50828, 0.47550), nrow = 4), tolerance=1e-2)
  
})

## Test center_data

test_that("SIGMA_matrix inverse works properly", {
  expect_equal(center_data(matrix(c(1,2,3,4,2,1,1,6,1, 5,4,10), nrow=4)), matrix(c(0.28249, 0.23754, 0.47996), nrow = 1), tolerance=1e-2)
  expect_equal(center_data(matrix(c(1.2, 5.6, 4.6, 3.2, 8.2 ,1.2, 4.3 ,1.56,1, 8.45, 96, 4), nrow=4)), matrix(c(0.2327, 0.21033, 0.55694), nrow = 1), tolerance=1e-2)
  
})


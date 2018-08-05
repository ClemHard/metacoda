library(testthat)

source("Rscript/coda.R");
source("Rscript/test_groupe.R")
source("Rscript/bootstrap.R")

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
  expect_equal(Base_binary_matrix(3), matrix(c(-1, -1, 1, -1, 0, 1), nrow = 2), tolerance=1e-2)
  expect_equal(Base_binary_matrix(5), matrix(c(-1, -1, -1, -1, 1, -1, -1, -1, 0,  1, -1, -1, 0, 0, 1, -1, 0, 0, 0, 1), nrow = 4), tolerance=1e-2)
  
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


## Test ilr_inverse

test_that("ilr_inverse works properly", {
  expect_equal(ilr_inverse(matrix(c(1,2,3,4,2,1,1,6,1, 5,4,10), nrow=4)), matrix(c(2.3167e-02, 4.9454e-04, 7.2938e-04, 4.9256e-08, 9.5293e-02, 8.3671e-03, 5.0759e-02, 1.4099e-05, 0.544, 0.00692, 0.02070, 0.00129, 0.33732, 0.98421, 0.92780, 0.99869), nrow = 4), tolerance=1e-2)
  expect_equal(ilr_inverse(matrix(c(1.2, 1.6, 0.6, 1.2, 2.12 ,1.2, 1.3 ,1.56,1, 8.45, 2.6, 4), nrow=4)), matrix(c(1.8045e-02, 1.1435e-05, 1.5833e-02, 2.127e-03, 0.09849, 0.000109, 0.03698, 0.01161, 0.56560, 0.00015, 0.11893, 0.03358, 0.3178, 0.9997, 0.8282, 0.9526), nrow = 4), tolerance=1e-2)
  
})


## Test center_data

test_that("SIGMA_matrix inverse works properly", {
  expect_equal(center_data(matrix(c(1,2,3,4,2,1,1,6,1, 5,4,10), nrow=4)), matrix(c(0.28249, 0.23754, 0.47996), nrow = 1), tolerance=1e-2)
  expect_equal(center_data(matrix(c(1.2, 5.6, 4.6, 3.2, 8.2 ,1.2, 4.3 ,1.56,1, 8.45, 96, 4), nrow=4)), matrix(c(0.2327, 0.21033, 0.55694), nrow = 1), tolerance=1e-2)
  
})

x <- matrix(c(48.29,48.83,45.61,45.5,49.27,46.53,48.12,47.93,46.96,49.16,48.41,47.9,48.45,48.98,48.74,49.61,49.2,
              2.33,2.47,1.7,1.54,3.3,1.99,2.34,2.32,2.01,2.73,2.47,2.24,2.35,2.48,2.44,3.03,2.5, 11.48,12.38,8.33,8.17,12.10,9.49,11.43,11.18,9.9,12.54,11.8,11.17,11.64,12.05,11.6,12.91,12.32,
              1.59,2.15,2.12,1.6,1.77,2.16,2.26,2.46,2.13,1.83,2.81,2.41,1.04,1.39,1.38,1.6,1.26,
              10.03,9.41,10.02,10.44,9.89,9.79,9.46,9.36,9.72,10.02,8.91,9.36,10.37,10.17,10.18,9.68,10.13, 0.18,0.17,0.17,0.17,0.17,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.17,0.18,
              13.58,11.08,23.06,23.87,10.46,19.28,13.65,14.33,18.31,10.05,12.52,14.64,13.23,11.18,12.35,8.84,10.51,
              9.85,10.64,6.98,6.79,9.65,8.18,9.87,9.64,8.58,10.55,10.18,9.58,10.13,10.83,10.45,10.96,11.05,
              1.9,2.02,1.33,1.28,2.25,1.54,1.89,1.86,1.58,2.09,1.93,1.82,1.89,1.73,1.67,2.24,2.02,
              0.44,0.47,0.32,0.31,0.65,0.38,0.46,0.45,0.37,0.56,0.48,0.41,0.45,0.8,0.79,0.55,0.48,
              0.23,0.24,0.16,0.15,0.3,0.18,0.22,0.21,0.19,0.26,0.23,0.21,0.23,0.24,0.23,0.27,0.23),
            nrow=17,
            dimnames = list(1:17,
                            paste0("v", 1:11)))

## Test variation_matrix

test_that("center_data works properly", {
  

  var_x <- matrix(c(0.000, 0.013, 0.006, 0.039, 0.001, 0.001, 0.049, 0.007, 0.009, 0.031, 0.012,
   0.013, 0.000, 0.003, 0.061, 0.020, 0.017, 0.109, 0.005, 0.002, 0.016, 0.000,
   0.006, 0.003, 0.000, 0.053, 0.012, 0.009, 0.089, 0.000, 0.002, 0.018, 0.002,
   0.039, 0.061, 0.053, 0.000, 0.047, 0.037, 0.057, 0.057, 0.053, 0.099, 0.062,
   0.001, 0.020, 0.012, 0.047, 0.000, 0.001, 0.040, 0.013, 0.016, 0.036, 0.018,
   0.001, 0.017, 0.009, 0.037, 0.001, 0.000, 0.042, 0.010, 0.013, 0.035, 0.016,
   0.049, 0.109, 0.089, 0.057, 0.040, 0.042, 0.000, 0.092, 0.097, 0.138, 0.106,
   0.007, 0.005, 0.000, 0.057, 0.013, 0.010, 0.092, 0.000, 0.003, 0.017, 0.004,
   0.009, 0.002, 0.002, 0.053, 0.016, 0.013, 0.097, 0.003, 0.000, 0.025, 0.002,
   0.031, 0.016, 0.018, 0.099, 0.036, 0.035, 0.138, 0.017, 0.025, 0.000, 0.015,
   0.012, 0.000, 0.002, 0.062, 0.018, 0.016, 0.106, 0.004, 0.002, 0.015, 0.000),nrow=11, byrow=TRUE)
  expect_equal(normalised_variation_matrix(x), var_x, tolerance=1e-2)
 
})



## Test totvar

test_that("totvar inverse works properly", {
  expect_equal(totvar(x), 0.316, tolerance=1e-2)
})


## Test marginal_univariate_distributions

test_that("marginal_univariate_distributions works properly", {
  
  set.seed(1)
  x <- mvrnorm(500, c(0, 2.3), matrix(c(1.1, 0,0, 5.1), nrow = 2))
  test <- marginal_univariate_distributions(ilr_inverse(x))
  
  expect_equal(test$Anderson_Darling, c(0.2368, 0.33341), tolerance=1e-2)
  expect_equal(test$Cramer_von_Mises, c(0.03092, 0.06432), tolerance=1e-2)
  expect_equal(test$Watson, c(0.02953, 0.06398), tolerance=1e-2)
  expect_equal(test$marginale_normale, c(1, 2), tolerance=1e-2)
})

## Test bivariate_angle_distribution

test_that("bivariate_angle_distribution works properly", {
  
  set.seed(1)
  x <- mvrnorm(400, c(0, 2.3, -6.1), matrix(c(1.1, 0,0, 0.1, 5, 0.3, 0, 0.2, 5), nrow = 3))
  test <- Bivariate_angle_distribution(ilr_inverse(x))
  
  expect_equal(test$Anderson_Darling, matrix(c(NA, NA, NA, 0.68517, NA, NA ,0.5184, 1.3479, NA), nrow=3), tolerance=1e-2)
  expect_equal(test$Cramer_von_Mises, matrix(c(NA, NA, NA, -4.9167, NA, NA ,4.3443, 8.311, NA), nrow=3), tolerance=1e-2)
  expect_equal(test$Watson, matrix(c(NA, NA, NA, -4.9012, NA, NA ,4.3321, 8.2877, NA), nrow=3), tolerance=1e-2)

})


## Test Raduis_test

test_that("Raduis_test works properly", {
  
  set.seed(1)
  x <- mvrnorm(400, c(0, 2.3, -6.1), matrix(c(1.1, 0,0, 0.1, 5, 0.3, 0, 0.2, 5), nrow = 3))
  test <- Raduis_test(ilr_inverse(x))
  
  expect_equal(test$Anderson_Darling, 0.3567, tolerance=1e-2)
  expect_equal(test$Cramer_von_Mises, 0.0401, tolerance=1e-2)
  expect_equal(test$Watson, 0.03998, tolerance=1e-2)
  
})


## Test testing

test_that("testing works properly", {
  
  set.seed(1)
  
  u1 <- c(0, 2.3, -6.1)
  u2 <- c(1, 5, 2)
  Sigma1 <- matrix(c(1.1, 0,0, 0.1, 5, 0.3, 0, 0.2, 5), nrow = 3)
  Sigma2 <- matrix(c(1.6, 0.8, 1, 0, 8, 2, 0, 2, 5), nrow = 3)
  
  data1 <- mvrnorm(500, u1, Sigma1)
  data2 <- mvrnorm(500, u1, Sigma1)
  test <- testing(ilr_inverse(data1), ilr_inverse(data2), case=1)
  
  expect_equal(test$statistic, 7.0388, tolerance=1e-2)
  expect_equal(test$quantile, 12.591, tolerance=1e-2)
  expect_equal(test$result, TRUE)
  
  
  data2 <- mvrnorm(500, u2, Sigma1)
  test <- testing(ilr_inverse(data1), ilr_inverse(data2), case=1)
  
  expect_equal(test$statistic, 1491.17, tolerance=1e-2)
  expect_equal(test$quantile, 12.591, tolerance=1e-2)
  expect_equal(test$result, FALSE)
  
  
  
  
  data1 <- mvrnorm(500, u2, Sigma2)
  test <- testing(ilr_inverse(data1), ilr_inverse(data2), case=1)
  
  expect_equal(test$statistic, 80.944, tolerance=1e-2)
  expect_equal(test$quantile, 12.591, tolerance=1e-2)
  expect_equal(test$result, FALSE)
  
  
  
  

  data1 <- mvrnorm(500, u2, Sigma1)
  data2 <- mvrnorm(500, u1, Sigma1)
  
  test <- testing(ilr_inverse(data1), ilr_inverse(data2), case=2)
  
  expect_equal(test$statistic, 7.592, tolerance=1e-2)
  expect_equal(test$quantile, 7.814, tolerance=1e-2)
  expect_equal(test$result, TRUE)
  
  
  data2 <- mvrnorm(500, u1, Sigma2)
  test <- testing(ilr_inverse(data1), ilr_inverse(data2), case=2)
  
  expect_equal(test$statistic, 85.264, tolerance=1e-2)
  expect_equal(test$quantile, 7.814, tolerance=1e-2)
  expect_equal(test$result, FALSE)
  
  
  
  
  data1 <- mvrnorm(500, u1, Sigma2)
  data2 <- mvrnorm(500, u1, Sigma1)
  
  test <- testing(ilr_inverse(data1), ilr_inverse(data2), case=3)
  
  expect_equal(test$statistic, 0.8905, tolerance=1e-2)
  expect_equal(test$quantile, 7.814, tolerance=1e-2)
  expect_equal(test$result, TRUE)
  
  
  data2 <- mvrnorm(500, u2, Sigma1)
  test <- testing(ilr_inverse(data1), ilr_inverse(data2), case=2)
  
  expect_equal(test$statistic, 78.162, tolerance=1e-2)
  expect_equal(test$quantile, 7.814, tolerance=1e-2)
  expect_equal(test$result, FALSE)
  
})



## Test count_to_proportion

test_that("count_to_proportion works properly", {
  expect_equal(count_to_proportion(matrix(c(1, 2, 3, 1.2, 1.3, 6.2, 4.2, 6.2, 0.25),nrow=3)), matrix(c(0.1562, 0.2105, 0.3174, 0.18750, 0.13684, 0.6560, 0.65625, 0.652631, 0.02645), nrow=3), tolerance=1e-2)
})



## Test MAP

test_that("MAP works properly", {
  expect_equal(MAP(matrix(c(1, 2, 3, 1.2, 1.3, 6.2, 4.2, 6.2, 0.25),nrow=3)), matrix(c(0.2127, 0.2400, 0.3212, 0.2340, 0.1840, 0.57831, 0.5531, 0.5760, 0.1004), nrow=3), tolerance=1e-2)
})



## Test max_abundance_value_OTU

test_that("max_abundance_value_OTU works properly", {
  set.seed(20283)
  
  data <- matrix(sample(0:5, 100, replace=TRUE), nrow=10)
  data_result <- data.frame(zero_inflated=c(0.3, 0.4, 0.4, 0.4, 0.6, 0.2, 0.3, 0.3, 0.3, 0.3), value=c(5, 3, 5, 1, 4, 2, 5, 1, 0, 2))
  expect_identical(max_abundance_value_OTU(data), data_result)
  
  
  data <- matrix(sample(0:5, 100, replace=TRUE), nrow=10)
  data_result <- data.frame(zero_inflated=c(0.3, 0.3, 0.2, 0.3, 0.3, 0.2, 0.3, 0.4, 0.4, 0.4), value=c(3, 2, 3, 4, 1, 1, 4, 4, 5, 1))
  expect_identical(max_abundance_value_OTU(data), data_result)
  
  unif <- runif(5,0,10)
  data <- matrix(sample(unif, 100, replace=TRUE), nrow=10)
  data_result <- data.frame(zero_inflated=c(0.3, 0.4, 0.4, 0.5, 0.3, 0.3, 0.4, 0.4, 0.5, 0.3), value=c(unif[2], unif[1], unif[2], unif[4], unif[5], unif[4], unif[2], unif[4], unif[1], unif[3]))
  expect_identical(max_abundance_value_OTU(data), data_result)
  
})




## Test zero_inflated

test_that("zero_inflated works properly", {
  set.seed(20283)
  
  data <- matrix(sample(0:5, 200, replace=TRUE), nrow=20)
  classification <- sample(1:4, 20, replace=TRUE)
  data_result <- tibble(OTU=c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8,
                              8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10) ,
                        cluster=c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4),
                        value=c(3, 3, 3, 3, 1, 1, 1, 1, 4, 4, 4, 4, 1, 1, 1, 1, 0, 0, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1, 1, 4, 4, 4, 4, 5, 5, 5, 5),
                        zero=c(0.488, 0.164, 0.423, 0, 0, 0.328, 0.282, 0.968, 0, 0.820, 0.282, 0.323, 0.488, 0,
                                0.423, 0.323, 0, 0.492, 0.282, 0, 0.488, 0.164, 0.141, 0.323, 0.732, 0, 0.141, 0.323,
                                0.244, 0.328, 0.141, 0.323, 0.244, 0.328, 0.282, 0.645, 0.488, 0.328, 0, 0.323),
                        nb_sample_cluster=c(4, 6, 7, 3, 4, 6, 7, 3, 4, 6, 7, 3, 4, 6, 7, 3, 4, 6, 7, 3, 4, 6, 7, 3, 4, 6, 7, 3, 4, 6, 7, 3, 4, 6, 7, 3, 4, 6, 7, 3),
                        zero_inflated_coeff=c(0.30, 0.30, 0.30, 0.30, 0.35, 0.35, 0.35, 0.35, 0.40, 0.40, 0.40, 0.40, 0.30, 0.30, 0.30, 0.30, 0.25,
                                               0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.35, 0.35,
                                               0.35, 0.35, 0.25, 0.25, 0.25, 0.25)) %>% as.matrix()
  
  
  expect_equal(zero_inflated(data, classification) %>% as.matrix(), data_result, tolerance=1e-2)
  
})


## Test check_table

test_that("check_table properly", {
  
  table1 <- table(c(1,2,3,1,2,2,3,2,1,2,2), c("e","r",rep("e",5), rep("t",4)))
  table2 <- table(c(1,2,3,1,2,2,3,2,1,2,2), c("e",rep("z",3),rep("r",5), rep("t",2)))
  
  table_result1 <- matrix(c(2, 2, 2, 0, 1, 0, 1, 3, 0, 0, 0, 0), nrow=3, dimnames = list(c("1", "2", "3"), c("e", "r", "t", "z"))) %>% as.table()
  table_result2 <- matrix(c(1, 0, 0, 1, 3, 1, 0, 2, 0, 1, 1, 1), nrow=3, dimnames = list(c("1", "2", "3"), c("e", "r", "t", "z"))) %>% as.table()
  
  expect_equal(check_table(table1, table2)[[1]], table_result1)
  expect_equal(check_table(table1, table2)[[2]], table_result2, check.attributes=FALSE)
  expect_equal(check_table(table1, table2) %>% length(), 2)
})



test_that("table_to_percentage_table", {
  
  table1 <- table(c(rep(0, 100), rep(1, 120), rep(3, 238)), c(rep("e", 144), rep("z", 314)))
  table_result <- matrix(c(69.44, 30.56, 0, 0, 24.2, 75.8), nrow=3, dimnames = list(c(0, 1, 3), c("e", "z")))
  expect_equal(table_to_percentage_table(table1), as.table(table_result), tolerance=1e-1, check.attributes=FALSE)
  
})


test_that("list_table_sum", {
  
  table1 <- table(c(rep(0, 100), rep(1, 120), rep(3, 238)), c(rep("e", 144), rep("z", 314)))
  table3 <- table(c(rep(1, 100), rep(0, 120), rep(4, 238)), c(rep("i", 144), rep("z", 314)))
  table2 <- table(c(1,2,3,1,2,2,3,2,1,2,2), c("e",rep("z",3),rep("r",5), rep("t",2)))
  table4 <- table(c(1,2,3,1,2,2,3,2,1,2,2), c("e",rep("p",3),rep("q",5), rep("t",2)))
  l <- list(list(table1, table2), list(table3, table4))
  
  table_result <- matrix(c(100, 44, 0, 0, 44, 100, 0, 0, 76, 76, 238, 238), nrow=4, dimnames = list(c(0, 1, 3, 4), c("e", "i", "z"))) %>% as.table()
  table_result1 <- matrix(c(2, 0, 0, 1, 1, 1, 1, 3, 1, 1, 3, 1, 0, 4, 0, 1, 1, 1), nrow=3, dimnames = list(c(1, 2, 3), c("e", "p", "q", "r", "t", "z"))) %>% as.table()
  
  expect_error(list_table_sum(table1), "l isn't a list")
  expect_warning(list_table_sum(list()), "list is empty")
  expect_equal(list(table1) %>% list_table_sum(), table1)
  expect_equal(list(list(table1, table2), list(table3, table4)) %>% list_table_sum(), list(table_result, table_result1))
  
})


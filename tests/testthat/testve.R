context("ve")

# Data
data(data_temp)
tolerance <- 0.001
set.seed(1)

# ve
logisticFit <- glm(disease_any ~ nAb1, data = data_temp, family = binomial())
logisticVE <- ve(logisticFit, data_temp, nboot = 500)

coxFit <- coxph(Surv(time_event, disease_any) ~ nAb1, data = data_temp)
coxVE <- ve(coxFit, data_temp, nboot = 500)

test_that("ve", 
          {expect_equal(logisticVE$VE, 0.2025612, tolerance = tolerance)
           expect_equal(logisticVE$CI$LB, -1.789451, tolerance = tolerance)
           expect_equal(logisticVE$CI$UB, 2.667624, tolerance = tolerance)
           expect_equal(coxVE$VE, 0.1198574, tolerance = tolerance)
           expect_equal(coxVE$CI$LB, -1.764782, tolerance = tolerance)
           expect_equal(coxVE$CI$UB, 2.070932, tolerance = tolerance)
           })

# glmParametricSampling
set.seed(1)
Data.vaccinated <- filter(data_temp, vaccine == 1)
Data.control <- filter(data_temp, vaccine == 0)
logisticFit <- glm(disease_any ~ nAb1, data = data_temp, family = binomial())
efficacySet <- glmParametricSampling(logisticFit, nboot = 500, Data.vaccinated, Data.control)
CI <- lapply(EfficacyCI(efficacySet),"*", 100)

test_that("glmParametricSampling", 
          {expect_equal(CI$mean, 0.2468702, tolerance = tolerance)
           expect_equal(CI$median, 0.08851155, tolerance = tolerance)
           expect_equal(CI$CILow, -1.789451, tolerance = tolerance)
           expect_equal(CI$CIHigh, 2.667624, tolerance = tolerance)
           expect_equal(length(efficacySet), 500)
           })

# coxphParametricSampling
set.seed(1)
Data.vaccinated <- filter(data_temp, vaccine == 1)
Data.control <- filter(data_temp, vaccine == 0)
coxFit <- coxph(Surv(time_event, disease_any) ~ nAb1, data = data_temp)
efficacySet <- coxphParametricSampling(coxFit, nboot = 500, Data.vaccinated, Data.control)
CI <- lapply(EfficacyCI(efficacySet),"*", 100)

test_that("coxphParametricSampling", 
          {expect_equal(CI$mean, 0.07762033, tolerance = tolerance)
            expect_equal(CI$median, 0.01538626, tolerance = tolerance)
            expect_equal(CI$CILow, -1.934508, tolerance = tolerance)
            expect_equal(CI$CIHigh, 2.330776, tolerance = tolerance)
            expect_equal(length(efficacySet), 500)
          })
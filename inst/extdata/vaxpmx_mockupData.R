set.seed(1)
ID <- seq(1,600)
nAb1 <- rnorm(600, mean = 5, sd = 2)
nAb2 <- rnorm(600, mean = 8, sd = 1)
group <- rbinom(n = 600, size = 1, prob = 0.5)
vaccine <- rbinom(n = 600, size = 1, prob = 0.5)
type_disease <- rbinom(n = 600, size = 2, prob = 0.1)
disease_any <- as.numeric(type_disease>0)
disease_type1 <- as.numeric(type_disease==1)
disease_type2 <- as.numeric(type_disease==2)

data_temp <- data.frame(ID, nAb1, nAb2, group, vaccine, type_disease, disease_any)

# save(data_temp, file = "data/data_temp.rda")
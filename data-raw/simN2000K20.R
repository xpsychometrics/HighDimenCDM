## code to prepare `simN2000K20` dataset goes here
# 1. Set a seed for reproducibility
set.seed(123)

# 2. Simulate the data
Q20 = readRDS('data-raw/Q20.rds')
K = ncol(Q20)
gs.parms = matrix(runif(nrow(Q20) * 2, 0.05, 0.25),ncol = 2)
cutoffs <- m <- rep(0,K)
vcov <- matrix(0.5,K,K)
diag(vcov) <- 1

df = GDINA::simGDINA(N = 2000, Q = Q20,gs.parm = gs.parms, att.dist = "mvnorm",
                 mvnorm.parm=list(mean = m, sigma = vcov,cutoffs = cutoffs))
simN2000K20 = list(data = df$dat, Q = Q20, true.params = df$catprob.parm)
usethis::use_data(simN2000K20, overwrite = TRUE)

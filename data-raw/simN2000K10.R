## code to prepare `simN2000K10` dataset goes here
# 1. Set a seed for reproducibility
set.seed(123)

# 2. Simulate the data
Q10 = readRDS('data-raw/Q10.rds')
gs.parms = matrix(runif(nrow(Q10) * 2, 0.05, 0.25),ncol = 2)
K = ncol(Q10)
cutoffs <- m <- rep(0,K)
vcov <- matrix(0.5,K,K)
diag(vcov) <- 1

df = GDINA::simGDINA(N = 2000, Q = Q10,gs.parm = gs.parms,, att.dist = "mvnorm",
                 mvnorm.parm=list(mean = m, sigma = vcov,cutoffs = cutoffs))

simN2000K10 = list(data = df$dat, Q = Q10, true.params = df$catprob.parm)
usethis::use_data(simN2000K10, overwrite = TRUE)

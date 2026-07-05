set.seed(123)

n <- 2000
rho_eu <- 0.1
rho_vw <- 0.4
Sigma <- matrix(c(1, rho_eu, 0,
                  rho_eu, 1, rho_vw,
                  0, rho_vw, 1), ncol = 3, byrow = TRUE)
latent <- mvtnorm::rmvnorm(n, sigma = Sigma)
colnames(latent) <- c("struct_error", "latent_endog", "instrument")

x_exog <- rnorm(n)
endog <- 0.2 * latent[, "instrument"] + latent[, "latent_endog"]
struct_error <- latent[, "struct_error"]

response <- 1 + 2 * endog + 1.5 * x_exog + struct_error

sim_endog <- data.frame(
  y = response,
  z_endog = endog,
  x_exog = x_exog,
  w_instr = latent[, "instrument"]
)

if (!dir.exists("data")) dir.create("data")

save(sim_endog, file = "data/sim_endog.rda", compress = "xz")

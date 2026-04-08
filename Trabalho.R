library(stats)
library(ggplot2)

#### Pergunta 2 ####
# Alinea a)


# Função para gerar trajetórias do processo de Ornstein-Uhlenbeck (comum a todas as alíneas)
generate_trajectories <- function(alpha, beta, sigma, X0, TT, dt, N, n_trajectories) {
  t <- seq(0, TT, by = dt)
  
  # Matrizes para armazenar as trajetórias
  trajectories_transition <- matrix(0, nrow = n_trajectories, ncol = N + 1)
  trajectories_milstein <- matrix(0, nrow = n_trajectories, ncol = N + 1)
  
  # Condições iniciais
  trajectories_transition[, 1] <- X0
  trajectories_milstein[, 1] <- X0
  
  # Criação das trajetórias
  set.seed(42)
  for (i in 1:n_trajectories) {
    for (j in 1:N) {
      # Método das densidades de transição
      mean <- trajectories_transition[i, j] * exp(-alpha * dt) + beta * (1 - exp(-alpha * dt))
      var <- (sigma^2 * (1 - exp(-2 * alpha * dt))) / (2 * alpha)
      trajectories_transition[i, j + 1] <- rnorm(1, mean, sqrt(var))
      
      # Esquema de Milstein
      X_t <- trajectories_milstein[i, j]
      dW <- rnorm(1, 0, sqrt(dt))
      trajectories_milstein[i, j + 1] <- X_t + alpha * (beta - X_t) * dt + 
        sigma * dW + 0.5 * sigma^2 * (dW^2 - dt)
    }
  }
  
  # Devolve as trajetórias e o tempo
  list(transition = trajectories_transition, milstein = trajectories_milstein, time = t)
}

# Alínea b)

n_trajectories <- 1

params_a <- generate_trajectories(
  alpha = 0.1,
  beta = 20,
  sigma = 0.3,
  X0 = 20,
  TT = 100,
  dt = 0.01,
  N = 10000, # T /dt = 10 000
  n_trajectories = n_trajectories
)

# Cores definidas para transição e milstein
cores_transicao <- c("red", "green", "blue")
cores_milstein <- c("purple", "orange", "deeppink")

# Ajustamos as margens para o gráfico caber 
par(mar = c(4, 4, 2, 1))

# Gráfico inicial com primeira trajetória
plot(
  params_a$time, params_a$transition[1, ],
  type = "l", col = cores_transicao[1],
  xlab = "Tempo", ylab = expression(X[t]),
  main = "Trajetórias do Processo Ornstein-Uhlenbeck",
  ylim = range(params_a$transition, params_a$milstein)
)
lines(params_a$time, params_a$milstein[1, ], col = cores_milstein[1], lty = 2)

# Adicionar restantes trajetórias, se existirem
if (n_trajectories_a > 1) {
  for (i in 2:n_trajectories_a) {
    cor_t <- cores_transicao[(i - 1) %% length(cores_transicao) + 1]
    cor_m <- cores_milstein[(i - 1) %% length(cores_milstein) + 1]
    lines(params_a$time, params_a$transition[i, ], col = cor_t)
    lines(params_a$time, params_a$milstein[i, ], col = cor_m, lty = 2)
  }
}
# Legendas automáticas
leg_labels <- c(
  sprintf("Transição %d", 1:n_trajectories_a),
  sprintf("Milstein %d", 1:n_trajectories_a)
)
leg_colors <- c(
  sapply(1:n_trajectories_a, function(i) cores_transicao[(i - 1) %% length(cores_transicao) + 1]),
  sapply(1:n_trajectories_a, function(i) cores_milstein[(i - 1) %% length(cores_milstein) + 1])
)
leg_lty <- c(rep(1, n_trajectories_a), rep(2, n_trajectories_a))

legend("topright", legend = leg_labels, col = leg_colors, lty = leg_lty, cex = 0.8)


# Alínea c)

# Função para estimativa inicial via mínimos quadrados
initial_estimates <- function(trajectory, dt) {
  n <- length(trajectory) - 1
  alpha_init <- 0.1  # Valor inicial aproximado baseado na alínea a
  beta_init <- mean(trajectory[2:(n+1)] - trajectory[1:n] * exp(-alpha_init * dt)) / (1 - exp(-alpha_init * dt))
  residuals <- trajectory[2:(n+1)] - (trajectory[1:n] * exp(-alpha_init * dt) + beta_init * (1 - exp(-alpha_init * dt)))
  sigma_init <- sqrt(mean(residuals^2) * (2 * alpha_init) / (1 - exp(-2 * alpha_init * dt)))
  c(alpha = alpha_init, beta = beta_init, sigma = sigma_init)
}

# Função para calcular as derivadas da log-verosimilhança
compute_derivatives <- function(trajectory, alpha, beta, sigma, dt) {
  n <- length(trajectory) - 1
  Delta <- dt
  
  # Variância da transição
  variance_term <- (sigma^2 / (2 * alpha)) * (1 - exp(-2 * alpha * Delta))
  if (variance_term <= 0 || is.na(variance_term)) return(c(dL_da = NA, d2L_da2 = NA))
  
  # Resíduos
  expected <- trajectory[1:n] * exp(-alpha * Delta) + beta * (1 - exp(-alpha * Delta))
  residuals <- trajectory[2:(n+1)] - expected
  
  # Termo de estabilização para evitar divisões por zero
  exp_term <- exp(-2 * alpha * Delta)
  denom <- 1 - exp_term
  if (abs(denom) < 1e-10) denom <- 1e-10  # Limite inferior para evitar NA
  
  # Primeira derivada em relação a alpha
  dL_da <- n / (2 * alpha) - (n * Delta * exp_term) / denom +
    sum(residuals^2 / (2 * variance_term) * (Delta * exp_term) / denom)
  
  # Segunda derivada em relação a alpha
  d2L_da2 <- -n / (2 * alpha^2) + (2 * n * Delta^2 * exp_term) / (denom^2) +
    sum(residuals^2 / (2 * variance_term) * ((Delta^2 * exp_term) / (denom^2) -
                                               (Delta * exp_term) / denom * (1 / variance_term) * (sigma^2 / alpha) * (1 - exp_term)))
  
  if (is.na(dL_da) || is.na(d2L_da2) || is.infinite(dL_da) || is.infinite(d2L_da2)) {
    return(c(dL_da = NA, d2L_da2 = NA))
  }
  
  c(dL_da = dL_da, d2L_da2 = d2L_da2)
}

# Função para estimadores de máxima verosimilhança usando Newton-Raphson para alpha
mle_estimators <- function(trajectory, dt, initial_params, max_iter = 100, tol = 1e-6) {
  n <- length(trajectory) - 1
  alpha <- initial_params["alpha"]
  beta <- initial_params["beta"]
  sigma <- initial_params["sigma"]
  
  # Iteração de Newton-Raphson para alpha
  for (iter in 1:max_iter) {
    
    # Calcular derivadas
    derivs <- compute_derivatives(trajectory, alpha, beta, sigma, dt)
    dL_da <- derivs["dL_da"]
    d2L_da2 <- derivs["d2L_da2"]
    
    # Verificar se as derivadas são válidas
    if (is.na(dL_da) || is.na(d2L_da2) || is.infinite(dL_da) || is.infinite(d2L_da2)) {
      warning("Derivadas inválidas na iteração ", iter)
      break
    }
    
    # Atualizar alpha
    alpha_new <- alpha - dL_da / d2L_da2
    
    # Verificar convergência
    if (!is.na(alpha_new) && abs(alpha_new - alpha) < tol) {
      # Após convergência de alpha, otimizar beta e sigma
      expected <- trajectory[1:n] * exp(-alpha * dt) + beta * (1 - exp(-alpha * dt))
      residuals <- trajectory[2:(n+1)] - expected
      beta <- mean(residuals + trajectory[1:n] * exp(-alpha * dt)) / (1 - exp(-alpha * dt))
      variance_term <- (sigma^2 / (2 * alpha)) * (1 - exp(-2 * alpha * dt))
      sigma <- sqrt(sum(residuals^2) / (n * variance_term) * (sigma^2 / 2))
      break
    }
    
    alpha <- alpha_new
    
    # Atualizar beta e sigma
    expected <- trajectory[1:n] * exp(-alpha * dt) + beta * (1 - exp(-alpha * dt))
    residuals <- trajectory[2:(n+1)] - expected
    beta <- mean(residuals + trajectory[1:n] * exp(-alpha * dt)) / (1 - exp(-alpha * dt))
    variance_term <- (sigma^2 / (2 * alpha)) * (1 - exp(-2 * alpha * dt))
    sigma <- sqrt(sum(residuals^2) / n * (2 * alpha) / (1 - exp(-2 * alpha * dt)))
  }
  
  if (iter == max_iter) {
    warning("Newton-Raphson não convergiu após ", max_iter, " iterações")
  }
  
  list(alpha = alpha, beta = beta, sigma = sigma, iterations = iter)
}

# Usar a trajetória gerada na alínea a
trajectory <- params_a$transition[1, ]
dt <- 0.01

# Estimativa inicial
initial_params <- initial_estimates(trajectory, dt)

# Calcular estimadores de máxima verosimilhança
estimates <- mle_estimators(trajectory, dt, initial_params)

# Resultados
cat("Estimativas de Máxima Verosimilhança:\n")
cat(sprintf("alpha = %.4f (real = 0.1)\n", estimates$alpha)) # 0.1
cat(sprintf("beta = %.4f (real = 20)\n", estimates$beta)) # 19.6483
cat(sprintf("sigma = %.4f (real = 0.3)\n", estimates$sigma)) # 0.3032
cat(sprintf("Iterações: %d\n", estimates$iterations)) # 1
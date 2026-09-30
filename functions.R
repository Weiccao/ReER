#---- basic functions ----#
weightMatrix <- function(tau, X, Y, beta) {
  epsilon <- Y - X %*% beta
  wDiag <- ifelse(epsilon < 0, 1 - tau, tau) %>% as.numeric()
  return(wDiag) # 保存对角线元素，简化计算
}

generateData <- function(N, simType = c('sim1', 'sim2'), uDist = "norm", doScale = FALSE, tau, 
                        seeds) {
  set.seed(seeds)

  beta <- c(1,2,2)
  gamma <- if (simType == 'sim1') c(1, 0, 0) else c(1, 0, 1)
  
  x1 <- runif(N); x2 <- runif(N)
  Xmat <- cbind(x1, x2)
  #X <- cbind(1, if (doScale) scale(Xmat, scale = FALSE) else Xmat)
  X <- cbind(1, Xmat)
  #Xc <- cbind(1, scale(Xmat, scale = FALSE))
  
  u <- if (uDist == 'norm') rnorm(N) else rt(N, df = 3)
  Y <- X %*% beta + X %*% gamma * u
  #Yc <- Xc %*% beta + Xc %*% gamma * u
  betaTrue <- beta + gamma * (if (uDist == 'norm') enorm(tau) else et(tau,df =3))
  
  return(list(X = X, Y = as.vector(Y), 
              betaTrue = betaTrue))
}

#---- matrix functions ----#
uMatrix <- function(tau, X, Y, beta) {
  W_vec <- weightMatrix(tau, X, Y, beta)
  #return(crossprod(X, W %*% X) / length(Y))
  return(crossprod(X * sqrt(W_vec)/length(Y)))
}

vMatrix <- function(tau, X, Y, beta, covType = c('Heter', 'Homo')) {
  W_vec <- weightMatrix(tau, X, Y, beta)
  epsilon <- Y - X %*% beta
  if (covType == 'Heter') {
    W <- diag(W_vec)
    V <- t(X) %*% W %*% epsilon %*% t(epsilon) %*% W %*% X / length(Y)
  } else {
    #w <- diag(W)
    V <- (sum((W_vec * epsilon)^2) *  crossprod(X)) / length(Y)
  }
  return(V)
}

qMatrix <- function(X, Y, beta, tau, covType = c('Heter', 'Homo')) {
  U <- uMatrix(tau, X, Y, beta)
  V <- vMatrix(tau, X, Y, beta, covType)
  
  if (rcond(U) > 1e-10) {
    U_inv <- solve(U)
  } else {
    U_inv <- ginv(U)
  }
  return(U_inv %*% V %*% U_inv)
}

#---- algorithm ----#
# 自定义的重加权expectile regression
er <- function(TAU, X, Y, max_iter = 50, tol = 1e-6) {
  # 初始化beta
  start_time <- Sys.time()
  beta <- solve(crossprod(X, X),crossprod(X, Y))
  
  for (iter in 1:max_iter) {
    w <- weightMatrix(TAU, X, Y, beta) 
    beta_new <- solve(
      crossprod(X, w * X),
      crossprod(X, w * Y)
    )
    if (max(abs(beta_new - beta)) < tol) break
    beta <- beta_new
  }
  
  end_time <- Sys.time()
  elapsed <-  as.numeric(difftime(end_time, start_time, units = "secs"))
  return(list(beta = beta, elapsed = elapsed))

}

## Song (2021) method
DCER <- function(X, Y, beta, Q_cum, Q_beta_cum, tau, covType = c('Heter', 'Homo')) {
  start_time <- Sys.time()
  Q_k <- qMatrix(X, Y, beta, tau, covType)
  if (rcond(Q_k) > 1e-10) {
    Q_k_inv <- solve(Q_k)
  } else {
    Q_k_inv <- ginv(Q_k)
  }
  Q_cum <- Q_cum + Q_k_inv
  Q_beta_cum <- Q_beta_cum + Q_k_inv%*% beta
  
  try({
    beta_new <- solve(Q_cum, Q_beta_cum)
  }, silent = TRUE)
  
  if (!exists("beta_new") || any(is.na(beta_new))) {
    beta_new <- ginv(Q_cum) %*% Q_beta_cum
  }
  
  # if (rcond(Q_cum) > 1e-10) {
  #   beta_new <- solve(Q_cum, Q_beta_cum)
  #   # beta_new <- solve(Q_cum) %*% Q_beta_cum
  # } else {
  #   beta_new <- ginv(Q_cum) %*% Q_beta_cum
  # }
  
  end_time <- Sys.time()
  elapsed <-  as.numeric(difftime(end_time, start_time, units = "secs"))
  return(list(beta_new = beta_new, Q_cum = Q_cum, Q_beta_cum = Q_beta_cum,
              elapsed = elapsed))
}

## Pan (2025) method
PAER <- function(X, Y, beta, Sigma_cum, Theta_cum) {
  start_time <- Sys.time()
  Sigma_cum <- Sigma_cum + crossprod(X)
  Theta_cum <- Theta_cum + crossprod(X, X %*% beta)
  
  try({
    beta_new <- solve(Sigma_cum, Theta_cum)
  }, silent = TRUE)
  
  if (!exists("beta_new") || any(is.na(beta_new))) {
    beta_new <- ginv(Sigma_cum) %*% Theta_cum
  }
  
  end_time <- Sys.time()
  elapsed <-  as.numeric(difftime(end_time, start_time, units = "secs"))
  
  return(list(beta_new = beta_new, Sigma_cum = Sigma_cum, Theta_cum = Theta_cum, 
              elapsed = elapsed))
}


# 只更新一次的ReER
ReER <- function(X, Y, beta, Z_cum, tau) {
  beta_iter <- beta
  start_time <- Sys.time()
  W_elapsed <- 0
  
  W_vec <- weightMatrix(tau, X, Y, beta_iter)
  W_mat <- crossprod(X, X * W_vec)
  #W_end <- Sys.time()
  #W_elapsed <- W_elapsed + as.numeric(difftime(W_end, W_time, units = "secs"))
  Hess <- Z_cum + W_mat
  Grad <- Z_cum %*% beta_iter + crossprod(X, Y * W_vec)
  try({
    beta_new <- solve(Hess, Grad)
  }, silent = TRUE)
  
  if (!exists("beta_new") || any(is.na(beta_new))) {
    beta_new <- ginv(Hess) %*% Grad
  }
  
  #elapsed <- max(0, as.numeric(difftime(end_time, start_time, units = "secs")) - W_elapsed)
  W_vec <- weightMatrix(tau, X, Y, beta_new)
  W_mat <- crossprod(X, X * W_vec)
  Z_cum <- Z_cum + W_mat
  end_time <- Sys.time()
  elapsed <-  as.numeric(difftime(end_time, start_time, units = "secs"))
  return(list(beta_new = beta_new, Z_cum = Z_cum,
              elapsed = elapsed, time_W = W_elapsed))
}

#---- assessment criteria ----#
summarize_table <- function(beta_true, beta_est) {
  # beta_true: 模拟中每轮的真值矩阵 (sim_time × p)
  # beta_est: 对应算法的估计矩阵 (sim_time × p)
  
  sim_time <- nrow(beta_true)
  p <- ncol(beta_true)
  
  Bias_vec <- numeric(p)
  Var_vec <- numeric(p)
  MSE_vec <- numeric(p)
  
  for (j in 1:p) {
    error <- beta_est[, j] - beta_true[, j]
    Bias_vec[j] <- mean(error)
    Var_vec[j] <- var(beta_est[, j])
    MSE_vec[j] <- mean(error^2)
  }
  
  summary_mat <- rbind(Bias = Bias_vec, Var = Var_vec, MSE = MSE_vec)
  return(summary_mat)
}


plot_summary_results <- function(summary_dir, save_plot = TRUE) {
  mse_df <- read_csv(file.path(summary_dir, "mse_summary.csv"), show_col_types = FALSE)
  time_df <- read_csv(file.path(summary_dir, "time_summary.csv"), show_col_types = FALSE)
  
  # 绘制时间
  p_time <- ggplot(time_df, aes(x = K, y = Time, color = Method)) +
    geom_line(aes(linetype = Method), size = 1) +
    geom_point(aes(shape = Method), size = 2) +
    labs(title = "Computation Time vs K", x = "K", y = "Average Time (s)") +
    theme_bw(base_size = 14)
  
  # 绘制beta的MSE
  mse_df$beta <- factor(mse_df$beta, labels = c("beta0", "beta1", "beta2"))
  
  p_mse <- ggplot(mse_df, aes(x = K, y = mean_MSE, color = Method, group = Method)) +
    geom_line(aes(linetype = Method), linewidth = 1) +
    geom_point(aes(shape = Method), size = 2) +
    facet_wrap(~ beta, labeller = label_parsed) +
    labs(title = "Mean MSE vs K by Method",
         x = "K", y = "Mean MSE") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom")
  
  
  if (save_plot) {
    plot_dir <- file.path(summary_dir, "plots")
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(filename = file.path(plot_dir, "MSE_vs_K_by_beta.png"), plot = p_mse, width = 8, height = 4)
    ggsave(filename = file.path(plot_dir, "MSE_vs_K_by_beta.pdf"), plot = p_mse, width = 8, height = 4)
    
    ggsave(filename = file.path(plot_dir, "Time.png"), plot = p_time, width = 7, height = 5)
    ggsave(filename = file.path(plot_dir, "Time.pdf"), plot = p_time, width = 7, height = 5)
    message("Saved plots to ", plot_dir)
  }
  
  return(list(mse_plot = p_mse, time_plot = p_time))
}


#---- online algorithm ----#
data_split <- function(X, Y, random.Split = FALSE, K){ # random.Split = FALSE 确保stream场景是包含关系
  #----------- 数据拆分 -----------
  #set.seed(123) 
  N <- length(Y)
  if(random.Split){
    shuffle_idx <- sample(N)  # 随机打乱索引
    X <- X[shuffle_idx, , drop = FALSE]
    Y <- Y[shuffle_idx]
  }
  
  ## 按照等量分批
  indices <- split(1:N, cut(1:N, breaks = K, labels = FALSE))
  
  X_split <- lapply(indices, function(idx) X[idx, , drop = FALSE])
  Y_split <- lapply(indices, function(idx) Y[idx])
 
  
  return(list(X_k = X_split, Y_k = Y_split))
  
}

online_alg <- function(X, Y, random.Split = FALSE,
                       K, tau, 
                       regFormula, covType = c('Heter', 'Homo'),
                       tol = 1e-10, max_iter = 100) {
  
  #----------- 初始化 -----------
  ## 初始化 β（第1批）
  X_1 <- X[[1]] 
  Y_1 <- Y[[1]]
  
  ## 初始化时间
  start_time <- Sys.time()
  
  ## 第一批次估计
  m <- er(tau, X_1, Y_1)
  beta_1 <- m$beta
  
  ## 三种算法的累积统计量
  
  ### ReER
  Z_cum <- 0
  W <- weightMatrix(tau, X_1, Y_1, beta_1)
  Z_cum <- Z_cum + t(X_1) %*% W %*% X_1
  
  ### PAER
  Sigma_cum <- Theta_cum <- 0
  Sigma_cum <- Sigma_cum + crossprod(X_1)
  Theta_cum <- Theta_cum + crossprod(X_1, X_1 %*% beta_1)
  
  ### DCER
  Q_cum <- Q_beta_cum <- 0
  Q_1 <- qMatrix(X_1, Y_1, beta_1, tau, covType = 'Heter')
  Q_cum <- Q_cum + Q_1
  Q_beta_cum <- Q_beta_cum + Q_1 %*% beta_1
  
  ### DCER2
  Q2_cum <- Q2_beta_cum <- 0
  Q2_1 <- qMatrix(X_1, Y_1, beta_1, tau, covType = 'Homo')
  Q2_cum <- Q2_cum + Q2_1
  Q2_beta_cum <- Q2_beta_cum + Q2_1 %*% beta_1
  
  #total_time_ReER <- total_time_DCER <- total_time_DCER2 <- total_time_PAER <- Sys.time() - start_time
  
  total_time_ReER <- total_time_DCER <- total_time_DCER2 <- total_time_PAER <- 0
  
  #----------- 更新估计 -----------
  # ReER估计路径（假设为矩阵，每行为一个beta）
  ReER_path <- matrix(NA, nrow = K, ncol = length(beta_1))
  ReER_path[1, ] <- beta_1
  
  for (k in 2:K) {
    X_k <- X[[k]]
    Y_k <- Y[[k]]
    Data_k <- as.data.frame(cbind(Y = Y_k, X_k[,-1]))
    
    ## static estimate
    m <- er(tau, X_k, Y_k)
    beta_k <- m$beta
    
    ### --- DCER ---
    DCER_k <- DCER(X_k, Y_k, beta_k, Q_cum, Q_beta_cum, tau, covType = 'Heter')
    DCER_est <- DCER_k$beta_new %>% t()
    Q_cum <- DCER_k$Q_cum
    Q_beta_cum <- DCER_k$Q_beta_cum
    total_time_DCER <- total_time_DCER + DCER_k$total_time_online + m$elapsed
    
    ### --- DCER2 ---
    # m <- expectreg.ls(formula = regFormula, data = Data_k, expectiles = tau)
    # beta_k <- c(m$intercepts, unlist(m$coefficients))
    
    DCER2_k <- DCER(X_k, Y_k, beta_k, Q_cum, Q_beta_cum, tau, covType = 'Homo')
    DCER2_est <- DCER2_k$beta_new %>% t()
    Q2_cum <- DCER2_k$Q_cum
    Q2_beta_cum <- DCER2_k$Q_beta_cum
    total_time_DCER2 <- total_time_DCER2 + DCER2_k$total_time_online + m$elapsed
    
    ### --- PAER ---
    PAER_k <- PAER(X_k, Y_k, beta_k, Sigma_cum, Theta_cum)
    PAER_est <- PAER_k$beta_new %>% t()
    Sigma_cum <- PAER_k$Sigma_cum
    Theta_cum <- PAER_k$Theta_cum
    total_time_PAER <- total_time_PAER + PAER_k$total_time_online + m$elapsed
    
    ### --- ReER ---
    beta_old <- ReER_path[(k - 1), ] # 区别于前两个模型 ReER 需要用到上一步骤的估计
    ReER_k <- ReER(X_k, Y_k, beta_old, Z_cum, tau = tau, tol, max_iter)
    ReER_est <- ReER_k$beta_new %>% t()
    Z_cum <- ReER_k$Z_cum
    total_time_ReER <- total_time_ReER + ReER_k$total_time_online
    ReER_path[k, ] <- ReER_est
  }
  
  beta_est <- list(ReER = ReER_est, DCER = DCER_est, PAER = PAER_est, DCER2 = DCER2_est)
  
  return(list(
    beta_est = beta_est,
    total_time_ReER = total_time_ReER,
    total_time_DCER = total_time_DCER,
    total_time_PAER = total_time_PAER,
    total_time_DCER2 = total_time_DCER2
  ))
}

plot_denisty <- function(all_estimates, save_dir) {
  cat("Generating Distribution vs Theory plots...\n")
  
  plot_dir <- file.path(save_dir, "plots_distribution")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 遍历每个 K
  for (k_str in names(all_estimates)) {
    # 1. 提取所有模拟的 ReER 估计值 (200 x dim 矩阵)
    beta_list <- all_estimates[[k_str]]$beta_hat
    reer_estimates <- do.call(rbind, lapply(beta_list, function(df) {
      as.numeric(df[df$Method == "ReER", -1])
    }))

    # Extract the true coefficient vector.  The vertical reference line in
    # each panel makes the finite-sample bias visible as the displacement
    # between the true value and the center of the estimated distribution.
    true_list <- all_estimates[[k_str]]$beta_true
    true_mat <- do.call(rbind, lapply(true_list, function(v) as.numeric(v)))
    true_betas <- colMeans(true_mat)
    
    # 2. 提取所有模拟的 ReER 理论 SE (200 x dim 矩阵)
    se_list <- all_estimates[[k_str]]$se_hat
    reer_ses <- do.call(rbind, lapply(se_list, function(v) as.numeric(v)))
    
    # 3. 计算 200 次模拟的均值 (作为正态分布的均值 mu)
    mean_betas <- colMeans(reer_estimates)
    
    # 4. 计算 200 次模拟的理论 SE 的均值 (作为正态分布的标准差 sigma)
    # 注意：这里使用 Mean_SE 更好的反映理论值，也可以尝试使用 Emp_SD 做对比
    mean_ses <- colMeans(reer_ses)
    
    # 整理数据用于 ggplot
    dim <- ncol(reer_estimates)
    plot_data_list <- list()
    curve_data_list <- list()
    true_line_data_list <- list()
    
    for (d in 1:dim) {
      beta_label <- paste0("beta[", d-1, "]")
      
      # 直方图数据
      plot_data_list[[d]] <- data.frame(
        Value = reer_estimates[, d],
        Parameter = beta_label
      )
      
      # 正态曲线数据
      mu <- mean_betas[d]
      sigma <- mean_ses[d] # 使用理论 SE 的均值
      
      # 生成曲线的 X 轴范围 (均值上下 4 个标准差)
      x_range <- seq(mu - 4*sigma, mu + 4*sigma, length.out = 200)
      curve_data_list[[d]] <- data.frame(
        x = x_range,
        y = dnorm(x_range, mean = mu, sd = sigma),
        Parameter = beta_label
      )

      true_line_data_list[[d]] <- data.frame(
        TrueBeta = true_betas[d],
        Parameter = beta_label
      )
    }
    
    df_plot <- do.call(rbind, plot_data_list)
    df_curve <- do.call(rbind, curve_data_list)
    df_true <- do.call(rbind, true_line_data_list)
    
    # 绘制图形
    p <- ggplot(df_plot, aes(x = Value)) +
      # 1. 绘制直方图 (y轴使用 density)
      geom_histogram(aes(y = after_stat(density)), 
                     bins = 30, fill = "#69b3a2", color = "#e9ecef", alpha = 0.7) +
      # 2. 叠加理论正态分布曲线
      geom_line(data = df_curve, aes(x = x, y = y), 
                color = "#404040", size = 1, linetype = "dashed") +
      # 3. 真实参数值；与直方图/正态曲线中心的距离直观反映 bias
      geom_vline(data = df_true, aes(xintercept = TrueBeta),
                 color = "#C44E52", linewidth = 0.8, linetype = "solid") +
      # 分面展示每个 beta
      facet_wrap(~Parameter, scales = "free", labeller = label_parsed) +
      theme_minimal() +
      labs(
        x = "Estimated Value",
        y = "Density"
      ) +
      theme_bw(base_size = 14) + 
      theme(
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold")
      )
    
    # 保存图片
    file_name <- paste0("denisty", k_str, ".pdf")
    ggsave(file.path(plot_dir, file_name), p, width = 8, height = 3, dpi = 300)
  }
  
  cat("Plots saved to:", plot_dir, "\n")
}

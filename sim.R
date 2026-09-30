#---- load package ----#
rm(list = ls()); gc()

library(expectreg)
library(dplyr)
library(ggplot2)
library(readr)
library(MASS)
library(openxlsx)
library(readr)
library(writexl)
library(purrr)
library(tidyr)
source('functions.R')


#---- 多次仿真模拟 ----#
# 模拟逻辑如下
## 生成数据200次，根据不同的K进行分割 --> 每次模拟数据不同K的数据相同
run_experiment <- function(sim_time = 200, dgpType = c('fixed','stream'), random.Split = FALSE,
                           N = 100000, nk = 100,
                           simType = c('sim1', 'sim2'), uDist = "norm", doScale = FALSE, regFormula, tau = 0.25, 
                           covType = c('Heter', 'Homo'), tol = 1e-8, max_iter = 100, 
                           save_dir = "./results/") {
  
  
  simType <- simType
  dgpType <- dgpType
  sim_name <- paste0(simType, '-', uDist)
  save_dir_simulate <- file.path(save_dir, sim_name, 'simulate')
  save_dir_sum <- file.path(save_dir, sim_name, 'summarize')
  dir.create(save_dir_simulate, showWarnings = FALSE, recursive = TRUE)
  dir.create(save_dir_sum, showWarnings = FALSE, recursive = TRUE)
  form <- as.formula(regFormula)
  
  if (dgpType == 'fixed') {
    K_all <- c(100, 200, 500, 1000)
    N_all <- rep(N, length(K_all))
  } else {
    K_all <- c(50, 100, 200, 500)
    N_all <- nk * K_all
  }
  
  se_records <- data.frame()
  mse_records <- data.frame()
  time_records <- data.frame()
  all_estimates <- list()  # 存储每次模拟的估计值
  
  # ------ 模拟开始 ------#
  for (sim in 1:sim_time) {
    cat("Running simulation #", sim, "\n")
    
    max_N <- max(N_all)
    data <- generateData(N = max_N, simType = simType, uDist = uDist,
                         doScale = doScale, tau = tau, seeds = 123 + sim)
    X_full <- data$X
    Y_full <- data$Y
    Beta_true <- data$betaTrue
    dim <- length(Beta_true)
    
    for (i in seq_along(K_all)) {
      k_val <- K_all[i]
      N_cur <- N_all[i]
      
      cat("\tEvaluating K =", k_val, "(N =", N_cur, ")\n")
      X <- X_full[1:N_cur, ]
      Y <- Y_full[1:N_cur]
      
      # Oracle
      Data_o <- as.data.frame(cbind(Y, X))
      names(Data_o) <- c('Y','x0','x1','x2')
      
      start_time_Oracle <- Sys.time()
      m_oracle <- er(tau, X, Y)
      Oracle_beta <- m_oracle$beta
      end_time_Oracle <- Sys.time()
      Time_Oracle <- as.numeric(difftime(end_time_Oracle, start_time_Oracle, units = "secs"))
      
      # Online algorithms
      ## split data
      Data_split <- data_split(X, Y,random.Split, K = k_val)
      
      X_split <- Data_split$X_k; Y_split <- Data_split$Y_k
      
      #----------- 初始化 -----------
      ## 初始化时间 
      total_time_ReER <- total_time_DCER <- total_time_PAER <- 0
      #total_time_ReER2 <- 0
      
      ## 初始化 β（第1批）
      X_1 <- X_split[[1]]
      Y_1 <- Y_split[[1]]
      
      ## 第一批次估计
      m_1<- er(tau, X_1, Y_1)
      beta_1 <- m_1$beta
      elapse <- m_1$elapsed
      total_time_ReER <- total_time_ReER + elapse
      total_time_DCER <- total_time_DCER + elapse
      total_time_PAER <- total_time_PAER + elapse
      #total_time_ReER2 <- total_time_ReER2 + elapse
      
      ## 三种算法的累积统计量
      
      ### ReER 只迭代一次
      Z_cum <- V_cum <- 0
      start_time <- Sys.time()
      w_vec <- weightMatrix(tau, X_1, Y_1, beta_1)
      W_mat <- crossprod(X_1, X_1 * w_vec)
      Z_cum <- Z_cum + W_mat
      end_time <- Sys.time()
      elapse <- as.numeric(difftime(end_time, start_time, units = "secs"))
      total_time_ReER <- total_time_ReER + elapse
      eps_vec <- as.numeric(Y_1 - X_1 %*% beta_1)
      V_cum <- V_cum + crossprod(X_1, (w_vec^2 * eps_vec^2) * X_1)
      
      ### PAER
      Sigma_cum <- Theta_cum <- 0
      start_time <- Sys.time()
      Sigma_cum <- Sigma_cum + crossprod(X_1)
      Theta_cum <- Theta_cum + crossprod(X_1, X_1 %*% beta_1)
      end_time <- Sys.time()
      elapse <- as.numeric(difftime(end_time, start_time, units = "secs"))
      total_time_PAER <- total_time_PAER + elapse
      
      ### DCER
      Q_cum <- Q_beta_cum <- 0
      start_time <- Sys.time()
      Q_1 <- qMatrix(X_1, Y_1, beta_1, tau, covType = 'Homo') # 默认求和在外
      Q_1_inv <- ginv(Q_1)
   
      Q_cum <- Q_cum + Q_1_inv
      Q_beta_cum <- Q_beta_cum + Q_1_inv%*% beta_1
      end_time <- Sys.time()
      elapse <- as.numeric(difftime(end_time, start_time, units = "secs"))
      total_time_DCER <- total_time_DCER + elapse
      
      
      #----------- Online learning -----------
      ### ReER 需要用到上一步估计需要存储该值
      ReER_path <- matrix(NA, nrow = k_val, ncol = dim)
      ReER_se_path <- matrix(NA, nrow = k_val, ncol = dim)
      
      ReER_path[1, ] <- beta_1
      H_inv_1 <- ginv(Z_cum) 
      cov_1 <- H_inv_1 %*% V_cum %*% H_inv_1
      ReER_se_path[1, ] <- sqrt(pmax(0, diag(cov_1)))
      
      for(k in 2:k_val){
        X_k <- X_split[[k]]
        Y_k <- Y_split[[k]]
        
        beta_old <- ReER_path[(k - 1), ] # 区别于前两个模型 ReER 需要用到上一步骤的估计
        ReER_k <- ReER(X_k, Y_k, beta_old, Z_cum, tau = tau)
        ReER_est <- ReER_k$beta_new %>% t()
        Z_cum <- ReER_k$Z_cum
        total_time_ReER <- total_time_ReER + ReER_k$elapsed
        ReER_path[k, ] <- ReER_est
        w_vec_k <- weightMatrix(tau, X_k, Y_k, as.numeric(ReER_k$beta_new))
        eps_vec_k <- as.numeric(Y_k - X_k %*% as.numeric(ReER_k$beta_new))
        V_cum <- V_cum + crossprod(X_k, (w_vec_k^2 * eps_vec_k^2) * X_k)
        
        # 计算当前步的 SE (Sandwich: H^-1 * V * H^-1)
        H_inv_k <- ginv(Z_cum)
        cov_k <- H_inv_k %*% V_cum %*% H_inv_k
        ReER_se_path[k, ] <- sqrt(pmax(0, diag(cov_k)))
        
        m_k<- er(tau, X_k, Y_k)
        beta_k <- m_k$beta
        elapse <- m_k$elapsed
        # if(dgpType == 'fixed'){
        #   elapse <- elapse/k_val
        # } # 固定样本估计时间只记录一次
        # total_time_ReER <- total_time_ReER + elapse
        total_time_DCER <- total_time_DCER + elapse
        total_time_PAER <- total_time_PAER + elapse
        
        ## PAER
        PAER_k <- PAER(X_k, Y_k, beta_k, Sigma_cum, Theta_cum)
        PAER_est <- PAER_k$beta_new %>% t()
        Sigma_cum <- PAER_k$Sigma_cum
        Theta_cum <- PAER_k$Theta_cum
        total_time_PAER <- total_time_PAER + PAER_k$elapsed
        
        ## DCER
        DCER_k <- DCER(X_k, Y_k, beta_k, Q_cum, Q_beta_cum, tau, covType = 'Homo')
        DCER_est <- DCER_k$beta_new %>% t()
        Q_cum <- DCER_k$Q_cum
        Q_beta_cum <- DCER_k$Q_beta_cum
        total_time_DCER <- total_time_DCER + DCER_k$elapsed
        
      }
      
      beta_online_est <- list(DCER = DCER_est, PAER = PAER_est,
                              ReER = ReER_est)
      ReER_final_se   <- ReER_se_path[k_val, ]
      est_methods <- c("Oracle", "DCER", "PAER", "ReER")

      time_vec <- c(Time_Oracle,total_time_DCER, 
                    total_time_PAER, total_time_ReER)
      
       # if(dgpType == 'stream'){
       #   time_vec[-1] <- time_vec[-1]/k_val
       # } # 流数据场景用一批次的时间
      
      # 记录估计值
      if (is.null(all_estimates[[as.character(k_val)]])) {
        all_estimates[[as.character(k_val)]] <- list(beta_hat = list(), se_hat = list(), 
                                                     beta_true = list(), time = list())
      }
      
      est_df <- data.frame(Method = est_methods,
                           t(sapply(est_methods, function(m) {
                             if (m == "Oracle") Oracle_beta else beta_online_est[[m]]
                           })))
      all_estimates[[as.character(k_val)]]$beta_hat[[sim]] <- est_df
      all_estimates[[as.character(k_val)]]$se_hat[[sim]] <- ReER_final_se
      all_estimates[[as.character(k_val)]]$beta_true[[sim]] <- Beta_true
      time_df <- data.frame(Method = est_methods, Time = time_vec)
      all_estimates[[as.character(k_val)]]$time[[sim]] <- time_df
    }
  }
  
  # ------ 模拟结束保存结果 ------#
  # 保存每个 K 的估计值为 Excel 多 sheet
  for (k_str in names(all_estimates)) {
    beta_list <- all_estimates[[k_str]]$beta_hat
    time_list <- all_estimates[[k_str]]$time
    true_list <- all_estimates[[k_str]]$beta_true
    
    # 生成 sheet 数据
    sheet_list <- list()
    
    # True 值
    true_mat <- do.call(rbind, true_list)
    colnames(true_mat) <- paste0("beta", 0:(ncol(true_mat) - 1))
    sheet_list[["True"]] <- as.data.frame(true_mat)
    
    # 估计值
    methods <- c("Oracle", "DCER", "PAER","ReER")
    for (method in methods) {
      est_mat <- do.call(rbind, lapply(beta_list, function(df) {
        as.numeric(df[df$Method == method, -1])
      }))
      colnames(est_mat) <- paste0("beta", 0:(ncol(est_mat) - 1))
      sheet_list[[method]] <- as.data.frame(est_mat)
    }
    
    # 时间记录
    time_mat <- do.call(rbind, lapply(time_list, function(df) {
      df$Time[match(methods, df$Method)]
    }))
    colnames(time_mat) <- methods
    sheet_list[["Time"]] <- as.data.frame(time_mat)
    
    # 写入 Excel
    out_path <- file.path(save_dir_simulate, paste0("estimates_K", k_str, ".xlsx"))
    openxlsx::write.xlsx(sheet_list, file = out_path)
    
    
    # ----基于输出的excel汇总结果-----
    method_summaries <- list()
    time_means <- colMeans(time_mat)
    se_mat_all <- do.call(rbind, all_estimates[[k_str]]$se_hat)
    for (d in 1:ncol(true_mat)) {
      # 创建空的 data frame 用于保存每个 beta 的结果
      summary_df <- data.frame(Method = character(),
                               Betatrue = numeric(),
                               Betahat = numeric(),
                               Bias = numeric(),
                               Var = numeric(),
                               MSE = numeric(),
                               Mean_SE = numeric(), 
                               Time = numeric(),
                               stringsAsFactors = FALSE)
      
      # 计算每个方法的指标
      for (i in seq_along(methods)) {
        m <- methods[i]
        est_mat <- sheet_list[[m]]
        
        # 获取当前维度的估计值
        est_d <- est_mat[, d]; est_ <-mean(est_d)
        # 获取真实值
        beta_d <- true_mat[, d]; beta_ <- mean(beta_d)
        
        # 计算 Bias, Var, MSE
        beta_diff <- est_d - beta_d
        bias  <- mean(beta_diff)
        var_  <- var(est_d)
        mse   <- mean(beta_diff^2)
        time_ <- time_means[i]
        
        m_se <- if(m == "ReER") mean(se_mat_all[, d]) else NA
        
        summary_df <- rbind(summary_df, 
                            data.frame(Method = m,
                                       BetaTrue = mean(beta_d),
                                       BetaHat = mean(est_d),
                                       Bias = bias,
                                       Var = var_,
                                       MSE = mse,
                                       Mean_SE = m_se,
                                       Time = time_means[i]))
      }
      
      # 保存每个 beta 的总结
      method_summaries[[paste0("beta_", d - 1)]] <- summary_df
      mse_records <- rbind(mse_records, data.frame(K = k_str, beta = d - 1, Method = summary_df$Method, mean_MSE = summary_df$MSE))
      
      se_records <- rbind(se_records, data.frame(
        K = k_str, 
        beta = d - 1, 
        Method = "ReER", 
        Mean_Beta = mean(sheet_list[["ReER"]][, d]),
        Mean_SE = mean(se_mat_all[, d])))
    }
    time_records <- rbind(time_records, data.frame(K = k_str, Method = methods, Time = time_means))
    
    out_path <- file.path(save_dir_sum, paste0("estimates_K", k_str, ".xlsx"))
    openxlsx::write.xlsx(method_summaries, file = out_path)
    cat("Saved summarized result to:", out_path, "\n\n")
    
  }
  write_csv(mse_records, file.path(save_dir_sum, "mse_summary.csv"))
  write_csv(time_records, file.path(save_dir_sum, "time_summary.csv"))
  write_csv(se_records, file.path(save_dir_sum, "reer_se_summary.csv"))
  cat("Saved MSE & Time summary CSVs\n")
  
  
  # ---- 绘制折线图 ----#
  plot_summary_results(summary_dir = save_dir_sum)
  plot_denisty(all_estimates = all_estimates, 
                            save_dir = save_dir_sum)
}

# ==============================================================================
# 测试部分
dgpType = 'fixed'
simType = "sim1"
uDist = "norm"
doScale = FALSE
regFormula = 'Y ~ x1 + x2'
Target_tau = 0.25
save_dir = "./test/"
N = 100000
random.Split = FALSE
tol = 1e-8
max_iter = 100
sim_time = 20

run_experiment(
  sim_time = 20, dgpType = 'fixed',
  simType = "sim1",
  uDist = "norm",
  doScale = FALSE,
  regFormula = 'Y ~ x1 + x2',
  tau = Target_tau,
  save_dir = "/Users/weiccao/Desktop/test/"
)

run_experiment(
  sim_time = 30, dgpType = 'stream',
  simType = "sim2",
  uDist = "t",
  doScale = FALSE,
  regFormula = 'Y ~ x1 + x2',
  tau = Target_tau,
  save_dir = "/Users/weiccao/Desktop/test/"
)

# ==============================================================================
# 模拟 N = 10000，固定场景
run_experiment(
  sim_time = 500, dgpType = 'fixed',
  simType = "sim1",
  uDist = "norm",
  doScale = FALSE,
  regFormula = 'Y ~ x1 + x2',
  tau = Target_tau,
  save_dir = "/Users/weiccao/Desktop/fix/"
)

# summary(save_dir = "./S2/fix/", simType = "sim1",
#         uDist = "norm")
# plot_all_beta_metrics(save_dir = "./S2/fix/", simType = "sim1",
#                       uDist = "norm")

run_experiment(
  sim_time = 500, dgpType = 'fixed',
  simType = "sim1",
  uDist = "t",
  doScale = FALSE,
  regFormula = 'Y ~ x1 + x2',
  tau = Target_tau,
  save_dir = "/Users/weiccao/Desktop/fix/"
)

# summary(save_dir = "./S2/fix/", simType = "sim1",
#         uDist = "t")
# plot_all_beta_metrics(save_dir = "./S2/fix/", simType = "sim1",
#                       uDist = "t")

run_experiment(
  sim_time = 500, dgpType = 'fixed',
  simType = "sim2",
  uDist = "norm",
  doScale = FALSE,
  regFormula = 'Y ~ x1 + x2',
  tau = Target_tau,
  save_dir = "/Users/weiccao/Desktop/fix/"
)

# summary(save_dir = "./S2/fix/", simType = "sim2",
#         uDist = "norm")
# plot_all_beta_metrics(save_dir = "./S2/fix/", simType = "sim2",
#                       uDist = "norm")

run_experiment(
  sim_time = 500, dgpType = 'fixed',
  simType = "sim2",
  uDist = "t",
  doScale = FALSE,
  regFormula = 'Y ~ x1 + x2',
  tau = Target_tau,
  save_dir = "/Users/weiccao/Desktop/fix/"
)

# summary(save_dir = "./S2/fix/", simType = "sim2",
#         uDist = "t")
# plot_all_beta_metrics(save_dir = "./S2/fix/", simType = "sim2",
#                       uDist = "t")

# ==============================================================================

# 模拟 nk = 300，stream场景
run_experiment(
  sim_time = 500, dgpType = 'stream',
  simType = "sim1",
  uDist = "norm",
  doScale = FALSE,
  regFormula = 'Y ~ x1 + x2',
  tau = Target_tau,
  save_dir = "/Users/weiccao/Desktop/stream/"
)

# summary(save_dir = "./S2/stream/", simType = "sim1",
#         uDist = "norm")
# plot_all_beta_metrics(save_dir = "./S2/stream/", simType = "sim1",
#                       uDist = "norm")

run_experiment(
  sim_time = 500, dgpType = 'stream',
  simType = "sim1",
  uDist = "t",
  doScale = FALSE,
  regFormula = 'Y ~ x1 + x2',
  tau = Target_tau,
  save_dir = "/Users/weiccao/Desktop/stream/"
)
# summary(save_dir = "./S2/stream/", simType = "sim1",
#         uDist = "t")
# plot_all_beta_metrics(save_dir = "./S2/stream/", simType = "sim1",
#                       uDist = "t")

run_experiment(
  sim_time = 500, dgpType = 'stream',
  simType = "sim2",
  uDist = "norm",
  doScale = FALSE,
  regFormula = 'Y ~ x1 + x2',
  tau = Target_tau,
  save_dir = "/Users/weiccao/Desktop/stream/"
)

# summary(save_dir ="./S2/stream/", simType = "sim2",
#         uDist = "norm")
# plot_all_beta_metrics(save_dir = "./S2/stream/", simType = "sim2",
#                       uDist = "norm")

run_experiment(
  sim_time = 500, dgpType = 'stream',
  simType = "sim2",
  uDist = "t",
  doScale = FALSE,
  regFormula = 'Y ~ x1 + x2',
  tau = Target_tau,
  save_dir = "/Users/weiccao/Desktop/stream/"
)

# summary(save_dir = "./S2/stream/", simType = "sim2",
#         uDist = "t")
# plot_all_beta_metrics(save_dir = "./S2/stream/", simType = "sim2",
#                       uDist = "t")

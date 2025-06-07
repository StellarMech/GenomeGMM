#' Estimate Genome Size using Gaussian Mixture Model
#'
#' This function reads a k-mer frequency histogram file, fits a Gaussian mixture model (GMM), estimates genome size, and generates a PDF visualization.
#'
#' @param data_file Path to the input k-mer frequency histogram file. The file should have two columns: the first column is k-mer frequency, the second column is the count of k-mers with that frequency.
#' @param n Number of Gaussian components (integer), usually set according to the expected ploidy or number of main peaks.
#'
#' @return No return value. The function outputs the estimated genome size in the console and saves a PDF plot in the working directory.
#'
#' @examples
#' # Example usage:
#' data_file <-"examples/4merged.histo.txt"
#' GenomeGMM (data_file, n = 4)
#'
#' @export
GenomeGMM <- function(data_file, n = 6, ...) {
# 检查数据文件是否存在
if (!file.exists(data_file)) {
  stop(paste("数据文件'", data_file, "'不存在，请检查文件路径。"))
}
# 读取完整数据
full_data <- read.table(data_file, header = FALSE)
# 将 y 值除以 1e6 以防止溢出
full_data[, 2] <- full_data[, 2] / 1e6 
# 查找第一个满足条件的点作为起始点
start_index <- 1
for (i in seq_len(nrow(full_data) - 1)) {
    current_y <- full_data[i, 2]
    next_y <- full_data[i + 1, 2]
    if (next_y >= current_y * 4/5) {
        start_index <- i
        break
    }
}
# 过滤数据，保留从起始点开始及之后的所有点
filtered_data <- full_data[start_index:nrow(full_data), ]
x_filtered <- filtered_data[, 1]
y_filtered <- filtered_data[, 2]
# 添加频率上限（核心修改点）
if (n == 6) {
  max_freq <- 100000  # n=6 专用高频截断阈值
} else {
  max_freq <- 100000  # 其他 n 值（如 n=2）的默认阈值
}  # 高频截断阈值（可调整，类似findGSE的maxright4fit）
effective_max_freq <- min(max_freq, max(x_filtered))  # 不超过实际数据中的最大频率
# 计算基于完整过滤数据的概率密度（用于后续转换）
total_samples_filtered <- sum(y_filtered)
bin_width <- ifelse(length(unique(diff(x_filtered))) == 1, diff(x_filtered)[1], 1)
histogram_density <- y_filtered / total_samples_filtered / bin_width
# 找到概率密度首次小于0.001的x值
right_index <- which(histogram_density < 0.00035)[1]
if (!is.na(right_index)) {
  right_value <- x_filtered[right_index]
  fit_range <- start_index:which(full_data[, 1] == right_value)
  #cat("自动确定的拟合右值为：", right_value, "\n")
} else {
  # 如果没有小于0.001的点，则用全部数据
  fit_range <- start_index:nrow(full_data)
  #cat("未找到概率密度小于0.001的点，使用全部数据\n")
}
# 截取用于拟合的数据
fitted_data <- full_data[fit_range, ]
x_fit <- fitted_data[, 1]
y_fit <- fitted_data[, 2]
total_x_times_y <- sum(filtered_data[, 1] * filtered_data[, 2])

# x_fit, y_fit 已经是用于拟合的数据
y_smooth <- stats::filter(y_fit, rep(1, 5)/5, sides=2)
y_smooth[is.na(y_smooth)] <- y_fit[is.na(y_smooth)]
# 找局部极大值
is_peak <- diff(sign(diff(y_smooth))) < 0
peak_indices <- which(is_peak) + 1  # +1是因为diff后长度变短
# 只保留明显的主峰（如高度大于全局最大值的20%）
peak_heights <- y_fit[peak_indices]
main_peak_indices <- peak_indices[peak_heights > max(y_fit) * 0.2]
first_peak_x <- x_fit[main_peak_indices[1]]
last_peak_x <- x_fit[main_peak_indices[length(main_peak_indices)]]
theoretical_miu <- seq(first_peak_x, last_peak_x, length.out = n)
# 初始化参数
  quantiles <- quantile(x_fit, seq(1/n, 1, length.out = n))
  alpha <- rep(1/n, n)
  miu <- as.numeric(quantiles)
  main_region <- full_data[full_data[, 1] > 10, ]
  hist_peak_index <- which.max(main_region[, 2])
  first_peak <- main_region[hist_peak_index, 1]
  sigma <- rep(sd(x_fit), n)
  #cat("\n【初始参数】\n")
#cat("初始混合系数 alpha:", alpha, "\n")
#cat("初始均值 miu:", miu, "\n")
#cat("初始标准差 sigma:", sigma, "\n")

# EM 算法参数设置（双重收敛条件）
if (n == 2) {
  max_iter <- 50000     # 增加迭代次数上限
  threshold <- 1e-8    # 参数收敛阈值
  rel_error <- 0.1      # 峰值倍数允许的相对误差（10%）
} else {
  max_iter <- 100000    # 增加迭代次数上限
  threshold <- 1e-8   # 参数收敛阈值
  rel_error <- 0.1      # 峰值倍数允许的相对误差（10%）
}
converged <- FALSE
both_conditions_met <- FALSE  # 同时满足两个条件的标记

# 异常值处理
y_safe <- y_fit
y_safe[y_safe < 0] <- 0

# EM 算法
prob <- matrix(0, nrow = length(x_fit), ncol = n)
weight <- matrix(0, nrow = length(x_fit), ncol = n)

for (step in 1:max_iter) {
  # E-step
  for (j in 1:n) {
    prob_j <- dnorm(x_fit, miu[j], sigma[j])
    prob_j[is.na(prob_j)] <- 1e-9
    weight[, j] <- alpha[j] * prob_j
  }
  row_weight <- rowSums(weight)
  row_weight[row_weight <= 0] <- 1e-9
  prob <- weight / row_weight
  
  old_alpha <- alpha
  old_miu <- miu
  old_sigma <- sigma
  
  # M-step
  for (j in 1:n) {
    resp_weights <- y_safe * prob[, j]
    total_weight <- sum(resp_weights)
    
    if (total_weight < 1e-9) {
      #cat("重置第", j, "个分量参数\n")
      sigma[j] <- sd(x_fit)
    }
     # 先更新均值和方差
    miu[j] <- sum(x_fit * resp_weights) / total_weight
    #miu[j] <- min(max(miu[j], theoretical_miu[j] * 0.95), theoretical_miu[j] * 1.05)
    target <- first_peak * j
    miu[j] <- min(max(miu[j], target * 0.8), target * 1.2)
    sigma[j] <- sqrt(sum((x_fit - miu[j])^2 * resp_weights) / total_weight)
    sigma[j] <- max(sigma[j], 1)
    sigma[j] <- min(sigma[j], 14)
}
# 再统一更新权重
  for (j in 1:n) {
  resp_weights <- y_safe * prob[, j]
  total_weight <- sum(resp_weights)
  alpha[j] <- total_weight / sum(y_safe)
  }
  # 1. 检查参数变化是否小于收敛阈值
  delta <- max(c(
    max(abs(alpha - old_alpha)),
    max(abs(miu - old_miu)),
    max(abs(sigma - old_sigma))
  ))
  params_converged <- (delta < threshold)
  peak_ratio_converged <- TRUE  # 单峰时无需检查
  # 3. 仅当满足时终止迭代
  if (params_converged) {
    #cat("✅ 算法在第", step, "步终止：参数变化小于阈值（", formatC(threshold, digits=2, format="e"), 
        #"）\n", sep="")
    converged <- TRUE
    break
  }
}
if (!converged) {
  #cat("⚠️ 警告：EM算法在最大迭代步数", max_iter, "内未收敛！\n")
 # cat("当前参数变化最大值为：", delta, "\n")
}

# 输出 EM 迭代后的参数
#cat("\n【EM 迭代后的参数】\n")
#cat("混合系数 alpha:", round(alpha, 4), "\n")
#cat("均值 miu:", round(miu, 4), "\n")
#cat("标准差 sigma:", round(sigma, 4), "\n")

# 绘图函数（使用截取的数据）
draw_plot <- function(x, y, alpha, miu, sigma, genome_size_bp) {
  k <- length(alpha)
  bin_width <- ifelse(length(unique(diff(x))) == 1, diff(x)[1], 0.5)
  #x_vals <- seq(min(x) - bin_width, max(x) + bin_width, length.out = 1000)
  x_max <- max(x)+50  # 或 max(x) + 50
  x_vals <- seq(min(x) - bin_width, x_max, length.out = 1000)

  total_samples <- sum(y)
  histogram_density <- y / total_samples / bin_width

  components <- lapply(1:k, function(j) alpha[j] * dnorm(x_vals, miu[j], sigma[j]))
  total_density <- Reduce(`+`, components)
  
  y_max <- max(c(histogram_density, total_density, sapply(components, max)), na.rm = TRUE) * 1.2
  
  # 创建图形
  plot(NULL, 
       xlim = c(min(x), x_max),
       ylim = c(0, y_max), 
       main = "Gaussian mixture model fitting results", 
       xlab = "kmer frequency", 
       ylab = "probability density"
  )

  rect(
    xleft = x - bin_width/2, 
    ybottom = 0, 
    xright = x + bin_width/2, 
    ytop = histogram_density, 
    col = "lightblue", 
    border = "white"
  )

  lines(x_vals, total_density, col = "red", lwd = 2)
  #colors <- c("blue", "green3", "orange", "purple", "brown", "pink", "cyan", "black")
  #colors <- rainbow(k)
  colors <- c(
  "#1B1B7A", "#15616D", "#78290F", "#3E065F", "#1B263B",
  "#2D6A4F", "#4B1D3F", "#2C3639", "#2B2D42", "#22223B",
  "#1A535C", "#3D405B", "#264653", "#283618", "#3A0CA3"
)[1:k]
  for (j in 1:k) {
    lines(x_vals, components[[j]], col = colors[j], lty = 2, lwd = 1.5)
  }
  
  # 将基因组大小转换为 Mb 和 Gb
  genome_size_mb <- genome_size_bp / 1e6
  genome_size_gb <- genome_size_bp / 1e9
  
  # 创建图例并添加基因组大小文本
  legend_size <- 1.2  # 图例文字大小
  legend("topright",
         legend = c("Total Density", paste("Component", 1:k)),
         col = c("red", colors)[1:(k + 1)],
         lty = c(1, rep(2, k)),
         lwd = 2,
         bty = "n",
         cex = legend_size
  )
  
  # 添加基因组大小文本到图例下方
  genome_text <- paste(
    "Genome size estimation:\n",
    signif(genome_size_bp, 5), "bp (\n",
    signif(genome_size_mb, 5), "Mb (\n",
    signif(genome_size_gb, 5), "Gb)"
  )
  text(x = max(x_vals) * 0.98, y = max(histogram_density) * 0.3,
       labels = genome_text,
       col = "black", cex = legend_size, font = 2, adj = c(1, 0))
}

# 计算拟合曲线密度（使用截取的数据）
k <- length(alpha)
x_vals <- seq(min(x_fit) - bin_width, max(x_fit) + bin_width, length.out = 1000)
components <- lapply(1:k, function(j) alpha[j] * dnorm(x_vals, miu[j], sigma[j]))
total_density <- Reduce(`+`, components)

# 计算指定自然数对应拟合曲线上的密度值
natural_numbers <- c(10, 20, 30)  # 要计算的自然数数组
density_values <- sapply(natural_numbers, function(x) {
  sum(sapply(1:k, function(j) alpha[j] * dnorm(x, miu[j], sigma[j])) )
})

# 计算kmer深度（使用截取的数据）
if (n == 2) {
  peak_index <- which.max(total_density)
  kmer_depth <- x_vals[peak_index]
} else {
  is_peak <- diff(sign(diff(total_density))) < 0
  peak_indices <- which(is_peak) + 1  # 调整索引
  
  if (length(peak_indices) > 0) {
    first_peak_index <- peak_indices[1]
    kmer_depth <- x_vals[first_peak_index]
    #cat("找到第一个峰，位置在索引", first_peak_index, "，频率为", kmer_depth, "\n")
  } else {
    #warning("未找到局部最大值，使用全局最大值作为替代")
    peak_index <- which.max(total_density)  # 定义peak_index用于兼容输出
    kmer_depth <- x_vals[peak_index]
  }
}

# 转换回频数（使用完整过滤数据的总数）
fitted_frequency <- total_density * total_samples_filtered * bin_width
# 计算拟合数据在左阈值到右阈值范围内的每个自然数的密度值
fit_range_upper <- max(fitted_data[, 1])  # 使用拟合数据的最大x值作为右边界
natural_x <- seq(start_index, fit_range_upper, by = 1)  # 生成自然数序列
# 计算每个自然数对应的密度值
density_natural_x <- approx(x_vals, total_density, xout = natural_x)$y
# 转换为频数
frequency_natural_x <- density_natural_x * total_samples_filtered * bin_width
# 计算加权频数之和
fit_weighted_frequencies <- natural_x * frequency_natural_x
fit_total_kmers <- sum(fit_weighted_frequencies)
# 输出结果
#cat("过滤后拟合部分的kmer总数 (N, 包含f<=拟合右值):", fit_total_kmers * 1e6, "\n")

# 计算拟合曲线以后的数据（即fit_range_upper之后）的加权频数之和，仅使用自然数
remaining_data <- full_data[full_data[, 1] > fit_range_upper, ]
remaining_x <- remaining_data[, 1]
remaining_y <- remaining_data[, 2]
# 计算加权频数之和
remaining_weighted_frequencies <- remaining_x * remaining_y
remaining_total_kmers <- sum(remaining_weighted_frequencies)
# 计算总数（拟合区域和拟合曲线以后数据的加权频数之和）
total_kmers <- fit_total_kmers + remaining_total_kmers
# 输出结果
cat("\n【最终计算结果】\n")
#cat("低频过滤阈值 (f >=):", start_index, "\n")
#cat("拟合曲线右值界限:", fit_range_upper, "\n")
#cat("拟合峰对应的kmer频率:", kmer_depth, "\n")  # 直接使用kmer_depth
#cat("调整后的平均测序深度 (m):", kmer_depth, "\n")
#cat("过滤后拟合部分的kmer总数 (N, 包含f<=拟合右值):", fit_total_kmers * 1e6, "\n")
#cat("原始数据加权频数之和 (f>拟合右值):", remaining_total_kmers * 1e6, "\n")
cat("基因组大小预估 (G = N/m):", total_kmers * 1e6 / kmer_depth, "bp\n")
# 保存并绘制图形
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
filename <- paste0("em_fitted_gmm_n", n, "_", timestamp, ".pdf")
pdf(filename, width = 8, height = 6)
draw_plot(x_fit, y_fit, alpha, miu, sigma, genome_size_bp = total_kmers * 1e6 / kmer_depth)
dev.off()
}
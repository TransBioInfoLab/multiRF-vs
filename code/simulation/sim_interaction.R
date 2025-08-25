# Load related function
dir.base <- "."
script <- list.files( 
  path = file.path(dir.base,"function"),
  pattern = "[.]R$", 
  full.names = T,
  recursive = T
)
for (f in script) source(f)

set.seed(0)
n_reps     <- 50
num_signal <- 7

doParallel::registerDoParallel(10)
all_metrics <- plyr::llply(
  seq_len(n_reps),
  .fun = function(i) {
    sim   <- simulate_multivar(seed = i, binary = F)
    imps  <- get_importances(sim, type = "regression")
    true_idx <- 1:num_signal
    met_i <- t(sapply(imps, evaluate_global, true_idx = true_idx))
    df_i  <- as.data.frame(met_i)
    df_i$Method    <- rownames(df_i)
    df_i$Replicate <- i
    df_i
  }, .parallel = T
)

df_metrics <- do.call(rbind, all_metrics)
rownames(df_metrics) <- NULL

write_csv(
  data.frame(df_metrics),
  file.path(dir.base, "results/simulation/results_reg_with_noise.csv")
)


library(ggplot2)
library(patchwork)

df_metrics$Method <- factor(df_metrics$Method, levels = c(paste0("IMD_a0", 0:9), "IMD", "MD", "RFuni", "GBM", "XGBoost", "PMDCCA", "RGCCA", "SPLS"))

# 1. Make each plot without its own legend position
p1 <- ggplot(df_metrics, aes(x=Method, y=AUC)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "grey") + theme_light() +
  labs(title="Ranking AUC") +
  theme(legend.position="bottom") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(df_metrics, aes(x=Method, y=PRAUC)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "grey") + theme_light() +
  labs(title="PR-AUC") +
  theme(legend.position="bottom") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p3 <- ggplot(df_metrics, aes(x=Method, y=Recall)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "grey") + theme_light() +
  labs(title="Recall@k") +
  theme(legend.position="bottom") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p4 <- ggplot(df_metrics, aes(x=Method, y=Precision)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "grey") + theme_light() +
  labs(title="Precision@k") +
  theme(legend.position="bottom") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p5 <- ggplot(df_metrics, aes(x=Method, y=F1)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "grey") + theme_light() +
  labs(title="F1@k") +
  theme(legend.position="bottom") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine them in a 2×2, collect one legend
grid <- (p1 + p2) / (p3 + p4 + p5) +
  plot_layout(guides = "collect")  # gather legends

# Combine with shared legend and title
final_plot <- grid +
  plot_annotation(
    title = "Performance Comparison of Interaction Mix Reg Scenario",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  ) &
  theme(legend.position = "bottom")

# Draw
print(final_plot)

ggsave(
  filename = file.path(dir.base, "results/simulation/plot/results_with_noise_reg.pdf"),
  width = 12,
  height = 8
)


# Load related function
dir.base <- "."
script <- list.files( 
  path = file.path(dir.base,"function"),
  pattern = "[.]R$", 
  full.names = T,
  recursive = T
)
for (f in script) source(f)


n <- c(100, 200, 200)
p <- c(200, 500, 1000)
psel = c(20, 30, 50)
n_reps <- 50

parameter <- list(
  s1 = list(rho = 0)
)

mod = "SimpleNonLinear"
table <- data.frame()

for(k in c(2,3)) {
  for(i in c(1:3)){
    sim.para <- list(
      n = n[i],
      p = p[i],
      rho = 0,
      psel = psel[i],
      j = k
    )
    
    doParallel::registerDoParallel(10)
    all_metrics <- plyr::llply(
      seq_len(n_reps),
      .fun = function(j) {
        dat <- sim.nonlinear2(
          n = sim.para$n,
          p = sim.para$p,
          rho = sim.para$rho,
          psel = sim.para$psel,
          j = sim.para$j,
          seed = j
        )
        keep.list <-  rep(list(sim.para$psel),sim.para$j)
        names(keep.list) <- names(dat$X)
        if(k == 2) connect_list <- list(c("X1", "X2")) else connect_list <- NULL
        imps  <- get_all_imp(dat, keep.list = keep.list, connect_list = connect_list)
        true_idx <- rep(c(rep(1,psel[i]), rep(0,p[i] - psel[i])), times = k)
        imps2 <- do.call(c,imps[1:3])
        imps <- c(imps2, imps[4:6])
        met_i <- t(sapply(imps, evaluate_varsel, true_idx = true_idx))
        df_i  <- as.data.frame(met_i)
        df_i$Method    <- rownames(df_i)
        df_i$Replicate <- j
        df_i
      }, .parallel = T
    )
    
    df_metrics <- do.call(rbind, all_metrics)
    rownames(df_metrics) <- NULL
    df_metrics$Scenario <- paste0(mod, "_n", n[i], "_p", p[i], "_psel", psel[i], "_j", k)
    write_csv(
      data.frame(df_metrics),
      file.path(dir.base, paste0("results/simulation/latent/results_latent_model", paste0(mod, "_n", n[i], "_p", p[i], "_psel", psel[i], "_j", k), ".csv"))
    )
    
    message("Finished evaluation ", paste0(mod, "_n", n[i], "_p", p[i], "_psel", psel[i], "_j", k))
    table <- rbind(table, df_metrics)
    
  }
  
}

write_csv(
  data.frame(table),
  file.path(dir.base, "results/simulation/results_latent_model.csv")
)

## Plot
table <- read_csv(
  file.path(dir.base, "results/simulation/results_latent_model.csv")
)

library(gtable)
library(gt)
library(ggpubr)
library(ggh4x)

tb1 <- table %>% 
  mutate(Setting = ifelse(str_extract(Scenario, "j2|j3") == "j2", "Two Omics", "Three Omics"),
         Scenario = ifelse(grepl("p200", Scenario), "S1", ifelse(grepl("p500", Scenario), "S2", "S3")),
         alpha = case_when(
           grepl("IMD_a01", Method) ~ "0",
           grepl("IMD_a05", Method) ~ "0.5",
           grepl("IMD.",  Method) ~ "1",
           TRUE                      ~ "Other"
         ) %>% factor(levels = c("0","0.5","1","Other")),
         Method = case_when(
           grepl("filter",   Method) ~ "Filter",
           grepl("mixture",  Method) ~ "Mixture",
           grepl("test", Method) ~ "Trans",
           TRUE                       ~ Method
         ),
         Method = factor(Method, levels = c("Filter", "Mixture", "Trans", "PMDCCA", "RGCCA", "SPLS")),
         Scenario = factor(Scenario),
         Setting  = factor(Setting, levels = c("Two Omics", "Three Omics"))) %>% 
  group_by(alpha, Scenario, Method, Setting) %>% 
  dplyr::summarise(Precision_se = sd(Precision),
                   Precision = mean(Precision),
                   Recall_se = sd(Recall),
                   Recall = mean(Recall),
                   "PR-AUC_se" = sd(PRAUC),
                   "PR-AUC" = mean(PRAUC),
                   Size_se = sd(Size),
                   Size = mean(Size)) %>%
  ungroup() %>% 
  mutate_if(is.numeric, ~formatC(.,digits = 2, format = "f")) %>% 
  dplyr::mutate(Precision = paste0(Precision," (",Precision_se,")"),
                Recall = paste0(Recall," (",Recall_se,")"),
                "PR-AUC" = paste0(`PR-AUC`," (",`PR-AUC_se`,")"),
                "Model Size" = paste0(Size," (",Size_se,")"),
                .keep = "unused") %>% 
  pivot_wider(names_from = c("Setting"), values_from = c("PR-AUC", "Precision", "Recall",  "Model Size"), names_vary = "slowest") %>%
  # mutate(Method = fct_relevel(Method, c("MRF-IMD-filter", "MRF-IMD-mixture",  "PMDCCA", "RGCCA", "SPLS"))) %>%
  group_by(alpha, Scenario) %>%
  gt(rowname_col = "Method") |>
  tab_spanner_delim(delim = "_", reverse = T) |> 
  cols_align(
    align = "center"
  )

write_csv(
  tb1$`_data`,
  file.path(dir.base, "results/simulation/results_latent_model_summ.csv")
)


table_long <- pivot_longer(table, cols = c("PRAUC", "AUC", "Recall", "Precision", "F1", "Size"),
                           names_to = "Evaluation Metrics", 
                           values_to = "Values") %>%
  mutate(Setting = ifelse(str_extract(Scenario, "j2|j3") == "j2", "Two Omics", "Three Omics"),
         Scenario = ifelse(grepl("p100", Scenario), "S1", ifelse(grepl("p500", Scenario), "S2", "S3"))) %>% 
  filter(`Evaluation Metrics` %in% c("PRAUC", "Precision", "Recall"))



tbl2 <- table_long %>%
  mutate(
    alpha = case_when(
      grepl("IMD_a01", Method) ~ "0.1",
      grepl("IMD_a05", Method) ~ "0.5",
      grepl("IMD.",  Method) ~ "1",
      TRUE                      ~ "Other"
    ) %>% factor(levels = c("0.1","0.5","1","Other")),
    Method = case_when(
      grepl("filter",   Method) ~ "Filter",
      grepl("mixture",  Method) ~ "Mixture",
      grepl("test", Method) ~ "Trans",
      TRUE                       ~ Method
    ),
    alpha  = factor(alpha, levels = c("0.1","0.5","1","Other")),
    Method = factor(Method, levels = c("Filter", "Mixture", "Trans", "PMDCCA", "RGCCA", "SPLS")),
    Scenario = factor(Scenario),
    Setting  = factor(Setting, levels = c("Two Omics", "Three Omics"))
  )

alpha_labels <- c(
  # `0.1`     = "MRF-IMD^{alpha==0.1}",
  # `0.5`   = "MRF-IMD^{alpha==0.5}",
  `1`     = "MRF-IMD",
  Other   = "\"Other Benchmark Methods\""
)

ggboxplot(
  data          = tbl2 %>% filter(alpha %in% c("1", "Other")),
  x             = "Setting",
  y             = "Values",
  xlab  = "Setitngs",
  ylab = "Scenarios",
  fill          = "Method",
  outlier.size  = .5,
  palette       = "npg"
) +geom_vline(xintercept = c(1.5), 
              linetype = "dashed", color = "grey30") + 
  ggh4x::facet_nested(
    Scenario ~ alpha + `Evaluation Metrics` ,     
    scales = "free_x", 
    space  = "free_x",
    strip  = ggh4x::strip_nested() ,                     # ggh4x helper for clean nested strips
    labeller = labeller(alpha = as_labeller(alpha_labels, label_parsed),
                        .default = label_value)
  ) +
  scale_y_continuous(limits = c(0,1)) +
  labs(
    title = "Performance Comparison of Latent Model Case",
    fill  = NULL
  ) +
  theme_bw() +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(
    plot.title        = element_text(hjust = 0.5, face = "bold", size = 14),
    # axis.text.x       = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y       = element_text(size = 10),
    # axis.title        = element_blank(),
    legend.position   = "top",
    legend.text = element_text(size = 12),
    legend.direction = "horizontal",
    legend.box       = "horizontal",
    strip.background  = element_rect(fill = "white", color = "grey20"),
    panel.border      = element_rect(color = "gray70", fill = NA),
    panel.spacing.x   = unit(0, "cm"),
    panel.spacing.y   = unit(0, "cm")
  ) 


ggsave(
  file.path(dir.base, "results/simulation/plot/results_latent_model.pdf"),
  width = 12,
  height = 8
)


alpha_labels <- c(
  `0.1`     = "MRF-IMD^{alpha==0.1}",
  `0.5`   = "MRF-IMD^{alpha==0.5}"
  # `1`     = "MRF-IMD",
  # Other   = "\"Other Benchmark Methods\""
)

ggboxplot(
  data          = tbl2 %>% filter(alpha %in% c("0.1", "0.5")),
  x             = "Setting",
  y             = "Values",
  xlab  = "Setitngs",
  ylab = "Scenarios",
  fill          = "Method",
  outlier.size  = .5,
  palette       = "npg"
) +geom_vline(xintercept = c(1.5), 
              linetype = "dashed", color = "grey30") + 
  ggh4x::facet_nested(
    Scenario ~ alpha + `Evaluation Metrics` ,     
    scales = "free_x", 
    space  = "free_x",
    strip  = ggh4x::strip_nested() ,                     # ggh4x helper for clean nested strips
    labeller = labeller(alpha = as_labeller(alpha_labels, label_parsed),
                        .default = label_value)
  ) +
  scale_y_continuous(limits = c(0,1)) +
  labs(
    title = "Performance Comparison of Latent Model Case in Different Alpha",
    fill  = NULL
  ) +
  theme_bw() +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(
    plot.title        = element_text(hjust = 0.5, face = "bold", size = 14),
    # axis.text.x       = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y       = element_text(size = 10),
    # axis.title        = element_blank(),
    legend.position   = "top",
    legend.text = element_text(size = 12),
    legend.direction = "horizontal",
    legend.box       = "horizontal",
    strip.background  = element_rect(fill = "white", color = "grey20"),
    panel.border      = element_rect(color = "gray70", fill = NA),
    panel.spacing.x   = unit(0, "cm"),
    panel.spacing.y   = unit(0, "cm")
  ) 


ggsave(
  file.path(dir.base, "results/simulation/plot/results_latent_model_imd_alpha.pdf"),
  width = 12,
  height = 8
)







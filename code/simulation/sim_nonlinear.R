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

parameter <- list(
  s1 = list(p.group = 2, q.group = 5, p.b = 0),
  s2 = list(p.group = 2, q.group = 10, p.b = 0)
)
n_reps <- 50

mod = "NonLinearReg"
table <- data.frame()

for(i in c(1:3)){
  k <- 1
  for(para in parameter){
    sim.para <- list(
      n = n[i],
      p = p[i],
      p.group = para$p.group,
      q.group = para$q.group,
      p.b = para$p.b
    )
    
    doParallel::registerDoParallel(10)
    all_metrics <- plyr::llply(
      seq_len(n_reps),
      .fun = function(j) {
        dat <-  sim.nonlinear3(
          n = sim.para$n,
          p = sim.para$p,
          p.group = sim.para$p.group,
          q.group = sim.para$q.group,
          p.b = sim.para$p.b,
          seed = j
        )
        keep.list <-  list(X1 = 2 * sim.para$p.group, X2 = sim.para$q.group)
        names(keep.list) <- names(dat$X)
        imps  <- get_all_imp(dat, keep.list = keep.list, connect_list = list(c("X1", "X2")))
        true_idx <- c(
          c(rep(1, (2 * sim.para$p.group * sim.para$q.group)), rep(0, sim.para$p - (2 * sim.para$p.group * sim.para$q.group))),
          c(rep(1, sim.para$q.group), rep(0, sim.para$p - sim.para$q.group))
        )
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
    df_metrics$Scenario <- paste0(mod, "_n", n[i], "_p", p[i], "_", names(parameter)[k])
    write_csv(
      data.frame(df_metrics),
      file.path(dir.base, paste0("results/simulation/nonlinear/results_nonlinear_model", paste0(mod, "_n", n[i], "_p", p[i], "_", names(parameter)[k]), ".csv"))
    )
    
    table <- rbind(table, df_metrics)
    k <- k + 1
   }
}
  

write_csv(
  data.frame(table),
  file.path(dir.base, "results/simulation/results_nonlinear_model.csv")
)

table <- read_csv(
  file.path(dir.base, "results/simulation/results_nonlinear_model.csv")
)

library(gtable)
library(gt)
library(ggpubr)
library(ggh4x)

tb1 <- table %>% 
  mutate(Setting = ifelse(str_extract(Scenario, "s1|s2") == "s1", "Setting 1", "Setting2"),
         Scenario = ifelse(grepl("p100", Scenario), "S1", ifelse(grepl("p500", Scenario), "S2", "S3")),
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
         Method = factor(Method, levels = c("Filter", "Mixture", "Trans", "PMDCCA", "RGCCA", "SPLS"))) %>% 
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
  file.path(dir.base, "results/simulation/results_nonlinear_model_summ.csv")
)

table_long <- pivot_longer(table, cols = c("PRAUC", "AUC", "Recall", "Precision", "F1", "Size"),
                           names_to = "Evaluation Metrics", 
                           values_to = "Values") %>%
  mutate(Setting = ifelse(str_extract(Scenario, "s1|s2") == "s1", "Setting 1", "Setting2"),
         Scenario = ifelse(grepl("p100", Scenario), "S1", ifelse(grepl("p500", Scenario), "S2", "S3"))) %>% 
  filter(`Evaluation Metrics` %in% c("PRAUC", "Precision", "Recall"))

tbl2 <- table_long %>%
       mutate(
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
             alpha  = factor(alpha, levels = c("0","0.5","1","Other")),
             Method = factor(Method, levels = c("Filter", "Mixture", "Trans", "PMDCCA", "RGCCA", "SPLS")),
             Scenario = factor(Scenario),
             Setting  = factor(Setting)
         )

alpha_labels <- c(
  # `0`     = "MRF-IMD^{alpha==0.1}",
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
    title = "Performance Comparison of Nonlinear Regression Case",
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
  file.path(dir.base, "results/simulation/plot/results_nonlinear_model.pdf"),
  width = 12,
  height = 8
)

alpha_labels <- c(
  `0`     = "MRF-IMD^{alpha==0.1}",
  `0.5`   = "MRF-IMD^{alpha==0.5}"
  # `1`     = "MRF-IMD",
  # Other   = "\"Other Benchmark Methods\""
)

ggboxplot(
  data          = tbl2 %>% filter(alpha %in% c("0", "0.5")),
  x             = "Setting",
  y             = "Values",
  xlab  = "Setitngs",
  ylab = "Scenarios",
  fill          = "Method",
  outlier.size  = .5,
  palette       = "npg"
) +
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
    title = "Performance Comparison of Nonlinear Regression Case in Different Alpha",
    fill  = NULL
  ) +
  theme_bw() +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(
    plot.title        = element_text(hjust = 0.5, face = "bold", size = 14),
    # axis.text.x       = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y       = element_text(size = 10),
    legend.position   = "top",
    legend.direction = "horizontal",
    legend.box       = "horizontal",
    strip.background  = element_rect(fill = "white", color = "grey20"),
    panel.border      = element_rect(color = "gray70", fill = NA),
    panel.spacing.x   = unit(0, "cm"),
    panel.spacing.y   = unit(0, "cm")
  ) 

ggsave(
  file.path(dir.base, "results/simulation/plot/results_nonlinear_imd_alpha.pdf"),
  width = 12,
  height = 8
)

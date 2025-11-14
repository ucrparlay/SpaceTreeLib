library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(latex2exp)

wash_name <- function(t) {
  # NOTE: Wash data
  t$solver[t$solver == "PTree-H"] <- "SPaC-H"
  t$solver[t$solver == "PTree-Z"] <- "SPaC-Z"
  t$solver[t$solver == "OrthTree"] <- "P-Orth"
  t$solver[t$solver == "KdTree"] <- "PkdTree"

  return(t)
}

plot <- function(t, bench, plot_name, rt, dis_type, titleSize, axisSize, legendSize) {
  # t <- wash_name(t)
  # colnames(t)[5] <- "k=1"
  # colnames(t)[8] <- "k=10"
  # colnames(t)[11] <- "k=100"
  knn_columns <- c(paste0(dis_type, "1"), paste0(dis_type, "10"), paste0(dis_type, "100"))
  colnames(t)[which(colnames(t) == knn_columns[1])] <- "k=1"
  colnames(t)[which(colnames(t) == knn_columns[2])] <- "k=10"
  colnames(t)[which(colnames(t) == knn_columns[3])] <- "k=100"
  knn_column_idx <- c(
    which(colnames(t) == "k=1"),
    which(colnames(t) == "k=10"),
    which(colnames(t) == "k=100")
  )
  # View(t)
  # nameLevel <- c("OrthTree", "PTree-H", "KdTree", "ZdTree", "PTree-Z", "CPAM-H", "CPAM-Z", "Boost")
  t$solver[t$solver == "P-Orth"] <- "P-Orth (Ours)"
  t$solver[t$solver == "SPaC-H"] <- "SPaC-H (Ours)"
  nameLevel <- c("P-Orth (Ours)", "SPaC-H (Ours)", "PkdTree", "ZdTree", "CPAM-H", "Boost", "SPaC-Z", "CPAM-Z")
  plot_solvers <- c("P-Orth (Ours)", "PkdTree", "SPaC-H (Ours)", "ZdTree", "CPAM-H", "Boost")
  # nameLevel <- c("P-Orth", "SPaC-H", "PkdTree", "ZdTree", "CPAM-H", "Boost", "SPaC-Z", "CPAM-Z")
  shapeLevel <- c(1, 0, 2, 5, 7, 14)
  # colors <- c(
  # as.character(lineColorMap[nameLevel[[1]]]),
  # as.character(lineColorMap[nameLevel[[2]]]),
  # as.character(lineColorMap[nameLevel[[3]]]),
  # as.character(lineColorMap[nameLevel[[4]]])
  # )
  # shapeLevel <- c(1, 2, 0, 3, 7, 8, 5, 14)
  # print(colors)
  t$solver <- factor(t$solver, levels = nameLevel)

  t <- t %>%
    filter(benchmark == bench) %>%
    filter(solver %in% plot_solvers) %>%
    filter(ratio == rt) %>%
    select(solver, knn_column_idx[1], knn_column_idx[2], knn_column_idx[3]) %>%
    pivot_longer(!solver, names_to = "K", values_to = "time") %>%
    arrange(desc(solver)) # NOTE: make Ours in front

  print(t)
  # t <- t %>% mutate(across(everything(), function(x) {
  #   replace(x, which(x < 0), NA)
  # }))

  t$K <- factor(t$K, levels = unique(t$K))
  names(t)[names(t) == "solver"] <- "Solvers" # NOTE: change column title

  # View(t)

  plt <- ggplot(t, aes(
    x = K, y = time, color = Solvers,
    shape = Solvers,
    group = Solvers
  )) +
    # geom_point(size = 3, alpha = 0.8) + # point property
    geom_point(size = 3, stroke = 1) + # point property
    geom_line(size = 0.8) +
    # scale_y_log10(breaks = c(0.1, 1, 10), label = c("0.1", "1.0", "10")) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x, n = 2),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
    ) +
    scale_shape_manual(values = shapeLevel) +
    ggtitle(plot_name) +
    theme(
      plot.margin = unit(c(.1, .1, .1, .1), "cm"),
      plot.title = element_text(hjust = 0.5, size = titleSize, color = "black"), # title size
      axis.text.x = element_text(size = axisSize, color = "black"), # axe size
      axis.text.y = element_text(size = axisSize, color = "black"),
      axis.title = element_blank(),
      legend.text = element_text(size = legendSize),
      legend.key.width = unit(0.3, "cm"), # Width of legend keys,
      legend.title = element_blank(),
      # axis.title = element_text(margin = margin(0.5, 0.1, 0.1, 0.1), size = 8),
    ) +
    guides(color = guide_legend(nrow = 1), line = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
    annotation_logticks(sides = "l") +
    # scale_color_viridis_d()
    scale_color_locuszoom()
  # coord_fixed(ratio = 1) +
  # scale_color_brewer(palette = "Set1")
  # scale_color_manual(values = colors)

  return(plt)
}

# show_col(pal_locuszoom("default")(7))

titleSize <- 12
axisSize <- 13
legendSize <- 13
faxisSize <- 14

t <- read.csv("data/incre_update/insert_2.csv")
t <- wash_name(t)
# sapply(t, typeof) // check column type
gs_ids <- plot(t, "Uniform_by_x", "Sweepline-IND", 0.0001, "IDS", titleSize, axisSize, legendSize)

gv_ids <- plot(t, "Varden", "Varden-IND", 0.0001, "IDS", titleSize, axisSize, legendSize)

gu_ids <- plot(t, "Uniform", "Uniform-IND", 0.0001, "IDS", titleSize, axisSize, legendSize)

gs_ods <- plot(t, "Uniform_by_x", "Sweepline-OOD", 0.0001, "ODS", titleSize, axisSize, legendSize)

gv_ods <- plot(t, "Varden", "Varden-OOD", 0.0001, "ODS", titleSize, axisSize, legendSize)

gu_ods <- plot(t, "Uniform", "Uniform-OOD", 0.0001, "ODS", titleSize, axisSize, legendSize)
glist <- list(gs_ids, gv_ids, gu_ids, gs_ods, gv_ods, gu_ods)
# glist <- list(gs_ods, gv_ods, gu_ods)
g <- ggpubr::ggarrange(
  plotlist = glist,
  nrow = 2,
  ncol = 3,
  common.legend = TRUE,
  legend = "top"
) + theme(
  plot.margin = grid::unit(c(0.1, 0.1, 0.1, 0.1), "mm")
)


annotate_figure(g,
  bottom = text_grob(TeX("$k$-NN search"), size = faxisSize),
  left = text_grob("Time in seconds (log-scale)", size = faxisSize, rot = 90)
)

ggsave("plots/fig4_knn.pdf", width = 7, height = 4.8, dpi = 1800)

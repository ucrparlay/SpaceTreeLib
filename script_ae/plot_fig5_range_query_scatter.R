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

plot <- function(t, bench, name, X, Y, titleSize, axisSize, legendSize, title) {
  # t <- wash_name(t)
  # t$solver[t$solver == "test_count"] <- "Ours (range count)"
  # t$solver[t$solver == "Ours"] <- "Ours"

  # nameLevel <- c("Ours", "Ours (range count)", "Log-tree", "BHL-tree", "CGAL")
  # nameLevel <- c("P-Orth", "SPaC-H", "PkdTree", "ZdTree", "CPAM-H", "Boost", "SPaC-Z", "CPAM-Z")
  t$solver[t$solver == "P-Orth"] <- "P-Orth (Ours)"
  t$solver[t$solver == "SPaC-H"] <- "SPaC-H (Ours)"
  nameLevel <- c("P-Orth (Ours)", "SPaC-H (Ours)", "PkdTree", "ZdTree", "CPAM-H", "Boost", "SPaC-Z", "CPAM-Z")
  plot_solvers <- c("P-Orth (Ours)", "PkdTree", "SPaC-H (Ours)", "ZdTree", "CPAM-H", "Boost")
  # shapeLevel <- c(1, 0, 2, 5, 7, 14)
  shapeLevel <- c(16, 15, 17, 5, 7, 14)
  # colors <- c(
  # as.character(lineColorMap[nameLevel[[1]]]),
  # as.character(lineColorMap[nameLevel[[2]]]),
  # as.character(lineColorMap[nameLevel[[3]]]),
  # as.character(lineColorMap[nameLevel[[4]]]),
  # as.character(lineColorMap[nameLevel[[5]]])
  # )
  # print(colors)
  t$solver <- factor(t$solver, levels = nameLevel)
  # t <- t %>% slice_sample(prop = 0.1)
  t <- t %>% filter(solver %in% plot_solvers)

  ls <- 10
  ms <- 20
  rs <- 10
  l <- 100
  r <- 10000
  tleft <- t %>%
    filter(benchmark == bench) %>%
    select(solver, type, size, time) %>%
    filter(size < l) %>%
    group_by(solver) %>%
    slice_sample(n = ls) %>%
    arrange(desc(solver))

  tBridge <- t %>%
    filter(benchmark == bench) %>%
    select(solver, type, size, time) %>%
    filter(size >= l & size <= r) %>%
    group_by(solver) %>%
    arrange(size) %>%
    slice_sample(n = ms) %>%
    # slice(round(seq(1, n(), length.out = 3))) %>%
    arrange(desc(solver))

  tright <- t %>%
    filter(benchmark == bench) %>%
    select(solver, type, size, time) %>%
    filter(size > r) %>%
    group_by(solver) %>%
    slice_sample(n = rs) %>%
    arrange(desc(solver))
  # t <- t %>%
  #   filter(benchmark == bench) %>%
  #   select(solver, type, size, time) %>%
  #   filter(size < l) %>%
  #   group_by(solver) %>%
  #   arrange(size) %>%
  #   # slice(ifelse(type == 1, seq(mid_proportion * n()), seq(proportion * n()))) %>%
  #   slice_sample(n=2) %>%
  #   arrange(desc(solver))

  # print(t)
  t <- bind_rows(tleft, tBridge, tright) %>%
    arrange(desc(solver))
  # names(t)[names(t) == "solver"] <- "Solvers" # NOTE: change column title

  g <- ggplot(t, aes(
    x = size, y = time, color = solver,
    shape = solver, group = solver
  )) +
    # geom_point(size = 3, alpha = 0.8) + # point property
    geom_point(size = 2.5, stroke = 0.5) + # point property
    # scale_y_log10() +
    # scale_x_log10() +
    # scale_y_log10(labels = function(x) format(x, scientific = TRUE)) +
    scale_y_log10(
      limits = c(1e-6, 1e-3),
      breaks = scales::trans_breaks("log10", function(x) 10^x, n = 3),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    # scale_x_log10(breaks = c(1e1, 1e3, 1e5)) +
    scale_x_log10(
      limits = c(10, 1e5),
      breaks = scales::trans_breaks("log10", function(x) 10^x, n = 3),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_shape_manual(values = shapeLevel) +
    ggtitle(name) +
    theme(
      plot.margin = unit(c(.1, .1, .1, .1), "cm"),
      plot.title = element_text(hjust = 0.5, size = titleSize, color = "black"), # title size
      axis.text.x = element_text(size = axisSize, color = "black"), # axe size
      axis.text.y = element_text(size = axisSize, color = "black"),
      axis.title = element_blank(),
      legend.text = element_text(size = legendSize),
      legend.title = element_blank(),
      legend.position = "top",
      legend.key.width = unit(0.3, "cm"),
      # axis.title = element_text(margin = margin(0.5, 0.1, 0.1, 0.1), size = 8),
    ) +
    guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1, override.aes = list(size = 3, stroke = 1))) +
    annotation_logticks(sides = "lb") +
    scale_color_locuszoom()
  # scale_color_viridis_d("rocket")

  # coord_fixed(ratio = 1) +
  # scale_color_brewer(palette = "Set1")
  # scale_color_manual(values = colors)
  return(g)
}

titleSize <- 12.5
axisSize <- 13
legendSize <- 13
faxisSize <- 13
t <- read.csv("data/range_query_log/range_query_log.csv")
t <- wash_name(t)
gux <- plot(t, "Uniform_by_x", "Sweepline", 3.32, 0.9, titleSize, axisSize, legendSize)
gv <- plot(t, "Varden", "Varden", 3.32, 3.5, titleSize, axisSize, legendSize)
gu <- plot(t, "Uniform", "Uniform", 3.32, 3.5, titleSize, axisSize, legendSize)

glist <- list(gux, gv, gu)
g <- ggpubr::ggarrange(
  plotlist = glist,
  nrow = 1,
  ncol = 3,
  common.legend = TRUE,
  legend = "top"
) + theme(
  plot.margin = grid::unit(c(0.5, 0.1, 0.1, 0.1), "mm")
)

# gv
annotate_figure(g,
  bottom = text_grob(TeX("Output Size (log-scale)"), size = faxisSize),
  left = text_grob("Time in seconds (log-scale)", size = faxisSize, rot = 90)
)

ggsave("plots/fig5_range_query.pdf", width = 7, height = 3, dpi = 1800)

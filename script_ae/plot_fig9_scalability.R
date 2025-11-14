library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(latex2exp)
library(scales)

wash_name <- function(t) {
  # NOTE: Wash data
  t$solver[t$solver == "PTree-H"] <- "SPaC-H"
  t$solver[t$solver == "PTree-Z"] <- "SPaC-Z"
  t$solver[t$solver == "OrthTree"] <- "P-Orth"
  t$solver[t$solver == "KdTree"] <- "PkdTree"

  return(t)
}

plot <- function(name, type, t, titleSize, axisSize, legendSize) {
  # t <- wash_name(t)

  # NOTE: level
  cs <- c("1", "4", "16", "32", "64", "112", "224")
  t$threads <- as.character(t$threads)
  t$threads <- factor(t$threads, levels = cs)
  # nameLevel <- c("OrthTree", "PTree-H", "KdTree", "ZdTree", "PTree-Z", "CPAM-H", "CPAM-Z", "Boost")
  # nameLevel <- c("P-Orth", "SPaC-H", "PkdTree", "ZdTree", "CPAM-H", "Boost", "SPaC-Z", "CPAM-Z")
  t$solver[t$solver == "P-Orth"] <- "P-Orth (Ours)"
  t$solver[t$solver == "SPaC-H"] <- "SPaC-H (Ours)"
  nameLevel <- c("P-Orth (Ours)", "SPaC-H (Ours)", "PkdTree", "ZdTree", "CPAM-H", "Boost", "SPaC-Z", "CPAM-Z")
  plot_solvers <- c("P-Orth (Ours)", "PkdTree", "SPaC-H (Ours)", "ZdTree", "CPAM-H", "Boost")
  shapeLevel <- c(1, 0, 2, 5, 7, 14)
  # colors <- c(
  # as.character(lineColorMap[nameLevel[[1]]]),
  # as.character(lineColorMap[nameLevel[[2]]]),
  # as.character(lineColorMap[nameLevel[[3]]]),
  # as.character(lineColorMap[nameLevel[[4]]])
  # )
  # shapeLevel <- c(1, 2, 9, 8)
  t$solver <- factor(t$solver, levels = nameLevel)

  ourValues <- t %>%
    filter(solver == "SPaC-H (Ours)", benchmark == name, threads == 1) %>%
    select(build, insert, delete)

  t <- t %>%
    filter(benchmark == name) %>%
    select(solver, build, insert, delete, threads) %>%
    group_by(solver) %>%
    arrange(threads, .by_group = TRUE) %>%
    mutate(build = ourValues$build[1] / build) %>%
    mutate(insert = ourValues$insert[1] / insert) %>%
    mutate(delete = ourValues$delete[1] / delete) %>%
    arrange(desc(solver))

  title <- ""
  if (type == "build") {
    title <- paste("Build tree, ", ifelse(name == "Uniform_by_x", "Sweepline", name))
  } else if (type == "insert") {
    title <- paste("Batch insert, ", ifelse(name == "Uniform_by_x", "Sweepline", name))
  } else if (type == "delete") {
    title <- paste("Batch delete, ", ifelse(name == "Uniform_by_x", "Sweepline", name))
  }

  cs_new <- c("1", "4", "16", "32", "64", "112", "112h")
  g <- ggplot(t, aes(
    x = threads, y = .data[[type]], color = solver,
    shape = solver, group = solver
  )) +
    # geom_point(size = 3, alpha = 0.8) + # point property
    geom_point(size = 2, stroke = 0.7) + # point property
    geom_line() +
    scale_x_discrete(labels = cs_new) +
    # scale_y_continuous(trans = "log2", breaks = c(1, 2, 4, 8, 16, 24, 48, 96)) +
    # scale_y_continuous(trans = "log2", breaks = c(0.1, 0.8, 1)) +
    scale_y_log10(
      # breaks = scales::trans_breaks("log10", function(x) 10^x, n = 3),
      # labels = scales::trans_format("log10", scales::math_format(10^.x))
      breaks = c(0.01, 0.1, 1.0, 10.0, 50.0),
      labels = label_number(accuracy = 1)
    ) +
    # scale_y_log10() +
    scale_shape_manual(values = shapeLevel) +
    ggtitle(title) +
    theme(
      plot.margin = unit(c(0, 0.15, -0.05, 0), "cm"),
      plot.title = element_text(hjust = 0.5, size = titleSize, margin = margin(, , 2, ), color = "black"), # title size
      axis.text.x = element_text(size = axisSize, angle = 30, color = "black"), # axe size
      axis.text.y = element_text(size = axisSize, color = "black"),
      axis.title = element_blank(),
      legend.text = element_text(size = legendSize),
      legend.title = element_blank(),
      legend.key.width = unit(0.3, "cm"),
      legend.box.margin = margin(0, 0, 0, 0)

      # legend.background = element_blank(),
      # legend.box.background = element_rect(colour = "grey")
      # legend.title = element_text(size = legendSize),
      # axis.title = element_text(margin = margin(0.5, 0.1, 0.1, 0.1), size = 8),
    ) +
    annotation_logticks(sides = "l") +
    scale_color_locuszoom()
  # coord_fixed(ratio = 1) +
  # scale_color_brewer(palette = "Set1")
  # scale_color_manual(values = colors)
  return(g)
}

titleSize <- 8.5
axisSize <- 8
legendSize <- 9
faxisSize <- 10
t <- read.csv("data/scalability/scalability.csv")
t <- wash_name(t)
# t_pargeo <- read.csv("../data/scalability_pargeo.csv")
# t <- t %>% bind_rows(t_pargeo)
# View(t)

gsb <- plot("Uniform_by_x", "build", t, titleSize, axisSize, legendSize)
gsi <- plot("Uniform_by_x", "insert", t, titleSize, axisSize, legendSize)
gsd <- plot("Uniform_by_x", "delete", t, titleSize, axisSize, legendSize)
gub <- plot("Varden", "build", t, titleSize, axisSize, legendSize)
gui <- plot("Varden", "insert", t, titleSize, axisSize, legendSize)
gud <- plot("Varden", "delete", t, titleSize, axisSize, legendSize)
gvb <- plot("Uniform", "build", t, titleSize, axisSize, legendSize)
gvi <- plot("Uniform", "insert", t, titleSize, axisSize, legendSize)
gvd <- plot("Uniform", "delete", t, titleSize, axisSize, legendSize)

glist <- list(gsb, gub, gvb, gsi, gui, gvi, gsd, gud, gvd)
# glist <- list(gsb, gub, gvb, gsi, gui, gvi)
g <- ggpubr::ggarrange(
  plotlist = glist,
  nrow = 3,
  ncol = 3,
  common.legend = TRUE,
  legend = "top"
) + theme(
  plot.margin = grid::unit(c(0, 0, 0, 0), "mm")
)

annotate_figure(g,
  bottom = text_grob(TeX("Number of processors"), size = faxisSize),
  left = text_grob("Normalized speedup (log-scale)", size = faxisSize, rot = 90)
)
ggsave("plots/fig9_scalability.pdf", width = 4.5, height = 4.8, dpi = 1800)


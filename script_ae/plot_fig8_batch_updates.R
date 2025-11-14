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

plot <- function(t, bench, title, titleSize, axisSize, legendSize, insert_tag) {
  t <- wash_name(t)

  t$solver[t$solver == "P-Orth"] <- "P-Orth (Ours)"
  t$solver[t$solver == "SPaC-H"] <- "SPaC-H (Ours)"
  nameLevel <- c("P-Orth (Ours)", "SPaC-H (Ours)", "PkdTree", "ZdTree", "CPAM-H", "Boost", "SPaC-Z", "CPAM-Z")
  plot_solvers <- c("P-Orth (Ours)", "PkdTree", "SPaC-H (Ours)", "ZdTree", "CPAM-H")
  # nameLevel <- c("P-Orth", "SPaC-H", "PkdTree", "ZdTree", "CPAM-H", "Boost", "SPaC-Z", "CPAM-Z")
  shapeLevel <- c(1, 0, 2, 5, 7, 14)

  t$solver <- factor(t$solver, levels = nameLevel)
  # colnames(t)[c(seq(from = 5, to = 17))] <- c(
  #   0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000
  # )
  # print(t)
  cinsert <- c("solver", "i0.1", "i0.2", "i0.5", "i1", "i2", "i5", "i10", "i20", "i50", "i100", "i200", "i500", "i1000")
  cdelete <- c("solver", "d0.1", "d0.2", "d0.5", "d1", "d2", "d5", "d10", "d20", "d50", "d100", "d200", "d500", "d1000")
  # insert_range <- ifelse(insert_tag == "1", cinsert, cdelete)
  # print(insert_range)
  # print("helloworld")
  if (insert_tag == "1") {
    insert_range <- cinsert
  } else {
    insert_range <- cdelete
  }
  t <- t %>%
    filter(benchmark == bench) %>%
    filter(solver %in% plot_solvers) %>%
    # select(-benchmark) %>%
    select(insert_range) %>%
    pivot_longer(!solver, names_to = "BatchSize", values_to = "time") %>%
    arrange(desc(solver))
  # print(t)
  # print(length(plot_solvers))

  t$BatchSize <- rep(c(
    0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000
  ), length(plot_solvers)) # WARN:should keep the repetation same as # solvers
  # names(t)[names(t) == "solver"] <- "Solvers" # NOTE: change column title


  g <- ggplot(t, aes(
    x = BatchSize, y = time, color = solver,
    shape = solver, group = solver
  )) +
    # geom_point(size = 3, alpha = 0.8) + # point property
    geom_point(size = 2, stroke = 1) + # point property
    geom_line() +
    # scale_x_discrete(labels = c(
    #   "0.1", "0.2", "0.5", "1",
    #   "2", "5", "10",
    #   "20", "50", "100", "200", "500", "1000"
    # )) +
    scale_y_log10(
      # breaks = c(0.1, 10)
      breaks = scales::trans_breaks("log10", function(x) 10^x, n = 2),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x, n = 5),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
      # breaks = c(0.1, 1, 10, 100, 1000),
      # labels = c("0.1", "1", "10", "100", "1000")
    ) +
    scale_shape_manual(values = shapeLevel) +
    ggtitle(title) +
    theme(
      plot.margin = unit(c(.05, .15, .05, .05), "cm"),
      plot.title = element_text(hjust = 0.5, size = titleSize, margin = margin(0.1, 0.1, 4, 0.1)), # title size
      axis.text.x = element_text(size = axisSize, color = "black"), # axe size
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title = element_blank(),
      legend.text = element_text(size = legendSize),
      # legend.title = element_text(margin = margin(0.1, 0.1, 0.1, 0.1)),
      legend.title = element_blank(),
      legend.key.width = unit(0.3, "cm"),
      legend.box.margin = margin(0, 0, 0, 0)
      # legend.background = element_blank(),
      # legend.box.background = element_rect(colour = "grey")
      # legend.title = element_text(size = legendSize),
      # axis.title = element_text(margin = margin(0.5, 0.1, 0.1, 0.1), size = 8),
    ) +
    annotation_logticks(sides = "lb") +
    scale_color_locuszoom()
  # coord_fixed(ratio = 1) +
  # scale_color_brewer(palette = "Set1")
  # scale_color_manual(values = colors)

  # scale_color_npg()
  # scale_color_tableau()
  # scale_color_tableau("Classic 10 Medium")
  return(g)
}

titleSize <- 9.9
axisSize <- 10
legendSize <- 11
faxisSize <- 11
# PERF: READ updates
ti <- read.csv("data/batch_updates/batch_updates.csv")
# td <- read.csv("../data/delete.csv")
giv <- plot(ti, "Uniform", "Batch insert, Uniform", titleSize, axisSize, legendSize, "1")
gis <- plot(ti, "Uniform_by_x", "Batch insert, Sweepline", titleSize, axisSize, legendSize, "1")
giu <- plot(ti, "Varden", "Batch insert, Varden", titleSize, axisSize, legendSize, "1")
gdv <- plot(ti, "Uniform", "Batch delete, Uniform", titleSize, axisSize, legendSize, "0")
gds <- plot(ti, "Uniform_by_x", "Batch delete, Sweepline", titleSize, axisSize, legendSize, "0")
gdu <- plot(ti, "Varden", "Batch delete, Varden", titleSize, axisSize, legendSize, "0")
glist <- list(gis, giv, giu, gds, gdv, gdu)
g <- ggpubr::ggarrange(
  plotlist = glist,
  nrow = 2,
  ncol = 3,
  common.legend = TRUE,
  legend = "top"
) + theme(
  plot.margin = grid::unit(c(0, 0, 0, 0), "mm")
)

annotate_figure(g,
  # bottom = text_grob(TeX("Batch size ($\\times 10$M)"), size = faxisSize),
  bottom = text_grob(TeX("Batch size ($\\times 1 M$)"), size = faxisSize),
  left = text_grob("Time in seconds (log-scale)", size = faxisSize, rot = 90)
)
ggsave("plots/fig8_batch_update.pdf", width = 5.5, height = 4, dpi = 3600)

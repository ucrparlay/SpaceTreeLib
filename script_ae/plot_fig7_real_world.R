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

nameLevel <- c("P-Orth", "ZdTree", "SPaC-H", "SPaC-Z", "CPAM-H", "CPAM-Z", "Boost", "PkdTree")

gen_build <- function(t, bench) {
  build <- t %>%
    filter(benchmark == bench) %>%
    filter(op == "Insert") %>%
    filter(ratio == 1.0) %>%
    select(solver, tot) %>%
    arrange(solver)

  return(build)
}

gen_update <- function(t, bench, op_name) {
  # update <- t %>%
  #   filter(benchmark == bench, op == op_name) %>%
  #   filter(ratio == 0.0001) %>%
  #   select(solver, tot, ratio) %>%
  #   pivot_wider(names_from = ratio, values_from = tot) %>%
  #   arrange(solver) %>%
  #   select(!solver)

  if (op_name == "Insert") {
    update_query <- t %>%
      filter(benchmark == bench, op == op_name) %>%
      filter(ratio == 0.0001) %>%
      select(solver, tot, IDS10, RCL) %>%
      arrange(solver) %>%
      select(!solver)
  } else {
    update_query <- t %>%
      filter(benchmark == bench, op == op_name) %>%
      filter(ratio == 0.0001) %>%
      select(solver, tot) %>%
      arrange(solver) %>%
      select(!solver)
  }

  # nt <- update %>% bind_cols(update_query)
  return(update_query)
}

gen_one_row <- function(t, bench) {
  t$solver <- factor(t$solver, levels = nameLevel)
  bt <- gen_build(t, bench)
  bt <- cbind(benchmark = bench, bt)
  it <- gen_update(t, bench, "Insert")
  dt <- gen_update(t, bench, "Delete")
  # Reorder: build, insert, delete, then the rest from insert (10NN, RG)
  return(bt %>% bind_cols(it[, 1, drop = FALSE]) %>% bind_cols(dt) %>% bind_cols(it[, 2:3]))
}

t <- read.csv("data/real_world/real_world.csv")

t <- wash_name(t)

# print(t)
tv <- gen_one_row(t, "Cosmo50")
tuni <- gen_one_row(t, "OSM")
t <- bind_rows(tv, tuni)

# Rename columns
colnames(t) <- c("Benchmark", "Solver", "Build", "Insert", "Delete", "10NN", "RG")

write.csv(t, "plots/fig7_real_world.csv", row.names = FALSE)
# print(t)
# rownames(df) <- NULL
# df
# write.csv(df, "../data/merged_summary.csv", row.names = FALSE)
# t

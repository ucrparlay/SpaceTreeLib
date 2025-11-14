library(dplyr)
library(tidyr)

wash_name <- function(t) {
  # NOTE: Wash data
  t$solver[t$solver == "PTree-H"] <- "SPaC-H"
  t$solver[t$solver == "PTree-Z"] <- "SPaC-Z"
  t$solver[t$solver == "OrthTree"] <- "P-Orth"
  t$solver[t$solver == "KdTree"] <- "PkdTree"

  t$benchmark[t$benchmark == "Uniform_by_x"] <- "Sweepline"

  return(t)
}

# nameLevel <- c("OrthTree", "PTree-H", "PTree-Z", "KdTree", "ZdTree", "CPAM-H", "CPAM-Z", "Boost")
nameLevel <- c("P-Orth", "ZdTree", "SPaC-H", "SPaC-Z", "CPAM-H", "CPAM-Z", "Boost", "PkdTree")

gen_build <- function(t, bench) {
  build <- t %>%
    filter(benchmark == bench) %>%
    filter(ratio == 1.0) %>%
    select(solver, tot, IDS10, ODS10, RCL, RRL) %>%
    arrange(solver)

  return(build)
}

gen_update <- function(t, bench) {
  update <- t %>%
    filter(benchmark == bench) %>%
    filter(ratio %in% c(0.1, 0.01, 0.001, 0.0001)) %>%
    select(solver, tot, ratio) %>%
    pivot_wider(names_from = ratio, values_from = tot) %>%
    arrange(solver) %>%
    select(!solver)
  # print(update)

  update_query <- t %>%
    filter(benchmark == bench) %>%
    filter(ratio == 0.0001) %>%
    select(solver, IDS10, ODS10, RCL, RRL) %>%
    arrange(solver) %>%
    select(!solver)
  # print(update_query)

  nt <- update %>% bind_cols(update_query)
  return(nt)
}

gen_one_row <- function(insert_t, delete_t, bench) {
  insert_t$solver <- factor(insert_t$solver, levels = nameLevel)
  delete_t$solver <- factor(delete_t$solver, levels = nameLevel)
  bt <- gen_build(insert_t, bench)
  bt <- cbind(benchmark = bench, bt)
  it <- gen_update(insert_t, bench)
  dt <- gen_update(delete_t, bench)
  return(bt %>% bind_cols(it) %>% bind_cols(dt))
  # return(bt)
}

insert_t <- read.csv("data/incre_update/insert_2.csv")
delete_t <- read.csv("data/incre_update/delete_2.csv")
# insert_t <- read.csv("data/insert_3.csv")
# delete_t <- read.csv("data/delete_3.csv")
insert_t <- wash_name(insert_t)
delete_t <- wash_name(delete_t)

tu <- gen_one_row(insert_t, delete_t, "Sweepline")
tv <- gen_one_row(insert_t, delete_t, "Varden")
tuni <- gen_one_row(insert_t, delete_t, "Uniform")
t <- bind_rows(tuni, tu, tv)
colnames(t) <- c("Benchmark", "Solver", "Build", "IND10NN", "OOD10NN", "Count", "List", "Ins10%","Ins1%", "Ins0.1%", "Ins0.01%", "IND10NN", "OOD10NN", "Count", "List","Del10%","Del1%", "Del0.1%", "Del0.01%",  "IND10NN", "OOD10NN", "Count", "List")
write.csv(t, "plots/fig3_summary.csv", row.names = FALSE)

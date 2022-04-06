library(dplyr)
library(reshape2)
library(ggplot2)

res = readLines("benchmark_res.txt")
header = c("Size", "Dataset", "F77time", "F77err", "F77nops", "Cpptime", "Cpperr", "Cppnops")

## Benchmark result for symmetric solver
sym_start = grep("eigs_sym", res) + 4
sym_end   = grep("========", res)[2] - 1
sym_dat   = strsplit(res[sym_start:sym_end], " ")
sym_dat   = lapply(sym_dat, function(x) as.numeric(x[x != ""]))
sym_dat   = as.data.frame(do.call(rbind, sym_dat))
colnames(sym_dat) = header

sym_medtime = sym_dat %>%
    select(Size, Dataset, F77time, Cpptime) %>%
    melt(id.vars = c("Size", "Dataset"),
         variable.name = "Package", value.name = "Time") %>%
    group_by(Size, Dataset, Package) %>%
    summarize(Medtime = median(Time)) %>%
    mutate(Package = ifelse(grepl("F77", Package), "ARPACK", "Spectra"))

sym_medtime$Size = paste("Matrix size:", sym_medtime$Size)

ggplot(sym_medtime, aes(x = factor(Dataset), y = Medtime)) +
    geom_bar(aes(fill = Package), width = 0.75, position = "dodge", stat = "identity") +
    facet_wrap(~ Size, scales = "free", ncol = 2) +
    xlab("Matrix ID") + ylab("Median Elapsed Time (ms)") +
    ggtitle("Symmetric Eigen Solver") +
    theme_bw(base_size = 20) +
    theme(plot.title = element_text(hjust = 0.5))

## Benchmark result for general solver
gen_start = grep("eigs_gen", res) + 4
gen_end   = grep("========", res)[4] - 1
gen_dat   = strsplit(res[gen_start:gen_end], " ")
gen_dat   = lapply(gen_dat, function(x) as.numeric(x[x != ""]))
gen_dat   = as.data.frame(do.call(rbind, gen_dat))
colnames(gen_dat) = header

gen_medtime = gen_dat %>%
    select(Size, Dataset, F77time, Cpptime) %>%
    melt(id.vars = c("Size", "Dataset"),
         variable.name = "Package", value.name = "Time") %>%
    group_by(Size, Dataset, Package) %>%
    summarize(Medtime = median(Time)) %>%
    mutate(Package = ifelse(grepl("F77", Package), "ARPACK", "Spectra"))

gen_medtime$Size = paste("Matrix size:", gen_medtime$Size)

ggplot(gen_medtime, aes(x = factor(Dataset), y = Medtime)) +
    geom_bar(aes(fill = Package), width = 0.75, position = "dodge", stat = "identity") +
    facet_wrap(~ Size, scales = "free", ncol = 2) +
    xlab("Matrix ID") + ylab("Median Elapsed Time (ms)") +
    ggtitle("General Eigen Solver") +
    theme_bw(base_size = 20) +
    theme(plot.title = element_text(hjust = 0.5))

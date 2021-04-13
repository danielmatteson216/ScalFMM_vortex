source("common.R")
source("gen_times_taskdep.R")
source("gen_efficiencies_taskdep.R")
source("gen_speedup_plots.R")
source("gen_parallel_efficiency_plots.R")
source("gen_normalized_time_plots.R")
source("gen_comm.R")
source("gen_gantt.R")
source("gen_speed_plots.R")
source("gen_group_size_plots.R")


data <- get_data_subset("loutre.db", 0L, 0L, "False", get_bsize_reference())
data <- subset(data, global_time >= 0)


# plot normalized time : simple mpi, implicit optimaux, explicit optimaux
output <- "output/normalized-time-final.pdf"
all_algo <- list("implicit-variable-bsize", "explicit-variable-bsize", "simple-mpi", "starpu")
d1 <- subset(data, algo %in% all_algo)
gen_normalized_time_plot(d1, all_algo, 0L)

# plot speed : implicit, implicit limited, explicit
output <- "output/normalized-time-limited.pdf"
all_algo <- list("implicit", "explicit", "implicit limited", "starpu")
d1 <- subset(data, algo %in% all_algo)
gen_normalized_time_plot(d1, all_algo, 0L)


# plot normalized time : implicit, explicit, implicit optimaux, explicite optimaux, simple mpi
output <- "output/normalized-time-bsize.pdf"
all_algo <- list("implicit-variable-bsize", "explicit-variable-bsize", "simple-mpi", "starpu", "implicit", "explicit")
# TODO récupérer l'ensemble des valeurs où il y a les 5
d1 <- subset(data, algo %in% all_algo)
d1 <- subset(d1, model == "cube" & npart == "70000000")
gen_normalized_time_plot(d1, all_algo, 0L)


#Plot for group size
data_all <- get_data_subset("loutre.db", 0L, 0L, "False", 0L)

# plot speed group size : implicit, explicit
output <- "output/bsize-speed.pdf"
all_algo <- list("implicit", "explicit")
d1 <- subset(data_all, algo %in% all_algo)
all_nnode <- unique(subset(d1, bsize != get_bsize_reference())$nnode)
tmp <- subset(d1, bsize != get_bsize_reference())
tmp <- unique(tmp[c("nnode", "npart")])
d1 <- subset(d1, d1[c("nnode","npart")] %IN% tmp)
gen_group_size_plot_speed(d1, 0L)

# plot comm group size : implicit
output <- "output/bsize-volume.pdf"
all_algo <- list("implicit")
d1 <- subset(data_all, algo %in% all_algo)
all_nnode <- unique(subset(d1, bsize != get_bsize_reference())$nnode)
tmp <- subset(d1, bsize != get_bsize_reference())
tmp <- unique(tmp[c("nnode", "npart")])
d1 <- subset(d1, d1[c("nnode","npart")] %IN% tmp)
gen_group_size_plot_comm(d1, 0L)


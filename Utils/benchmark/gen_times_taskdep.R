library(plyr)
library(reshape)
library(ggplot2)


gen_times_taskdep_plot <- function(data, algo_wanted, model_wanted)
{
    # Sort data to have task, runtime and idle.
	subdata <- subset(data, model == model_wanted & algo == algo_wanted)
	subdata$rmem <- NULL
	subdata$communication_vol <- NULL
	subdata$global_time <- NULL
	subdata <- melt(subdata, id=c("model", "algo", "nnode", "nthreads", "npart","height","bsize"))
	subdata <- rename(subdata, c("variable"="event", "value"="duration"))

    g <- ggplot(data=subdata, aes(x=nnode, y=duration, fill=event))
	g <- g + geom_bar(stat="identity", position="fill")

	#Pour le titre
	g <- g + facet_wrap(npart ~ height ~ bsize, scales="free",
						labeller = labeller(npart = as_labeller(npart_labeller),
											height = as_labeller(height_labeller),
											bsize = as_labeller(group_size_labeller),
											.default=label_both,
											.multi_line=FALSE))
    # Set colors.
	breaks <- c('idle_time', 'communication_time', 'runtime_time', 'scheduling_time', 'task_time')
	labels <- c('Idle', 'Communication', 'Runtime', 'Scheduling', 'Task')
	colors <- c(
		'task_time'      = "#619dff",
		'runtime_time'   = "#01ba38",
		'idle_time'      = "#f9766e",
		'scheduling_time'   = "#fdbc43",
		'communication_time'   = "#9D2FEA"
	)
	g <- g + scale_fill_manual(name="Time", breaks=breaks,
							   labels=labels, values=colors)

    # Set title and X/Y labels.
    g <- g + xlab("Number of nodes")
    g <- g + ylab("% of time")
	g <- g + get_theme()

	output <- paste(get_output_directory(), "/", model_wanted, "-", algo_wanted, "-times.pdf", sep="")

    # Save generated plot.
    ggsave(output, g, width=29.7, height=21, units=c("cm"), device=cairo_pdf)
}

#Use this function to normalize
compute_timings <- function(data, n, h, m)
{
    # Select data
    sdata <- subset(data, npart == n & height == h & model == m)

    # Select Tt(1)
    tt_1 <- subset(sdata, event == "task_time" & nthreads == 1)

    # Compute task efficiency
    tt_p <- subset(sdata, event == "task_time")
    tt_p$event <- "task"
    tt_p$efficiency <- tt_p$duration / tt_1$duration

    # Compute scheduling efficiency
    ts_p <- subset(sdata, event == "scheduling_time")
    ts_p$event <- "scheduling"
    ts_p$efficiency <- ts_p$duration / tt_1$duration

    # Compute runtime efficiency
    tr_p <- subset(sdata, event == "runtime_time")
    tr_p$event <- "runtime"
    tr_p$efficiency <- tr_p$duration / tt_1$duration

    # Compute pipeline efficiency
    ti_p <- subset(sdata, event == "idle_time")
    ti_p$event <- "idle"
    ti_p$efficiency <- ti_p$duration / tt_1$duration

    # Merge all efficiencies
    sdata <- rbind(tt_p, ts_p, tr_p, ti_p)

    return (sdata)
}

gen_times_taskdep <- function(dbfile)
{
    # Cube (volume)
    data <- get_data_subset(dbfile, 0L, 0L, "False", get_bsize_reference())

	all_algorithm <- unique(data$algo)
	all_model <- unique(data$model)
	for (i in 1:length(all_algorithm))
	{
		if(all_algorithm[i] != "simple-mpi")
		{
			for (j in 1:length(all_model))
			{
				gen_times_taskdep_plot(data, all_algorithm[i], all_model[j])
			}
		}
	}
}

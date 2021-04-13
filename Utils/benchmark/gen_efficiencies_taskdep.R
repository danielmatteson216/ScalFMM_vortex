library(reshape)
library(ggplot2)
library(plyr)

gen_efficiencies_plot <- function(output, data)
{
    g <- ggplot(data=data, aes_string(x="nnode", y="efficiency",
                color="event", group="event"))
    g <- g + geom_line()
    g <- g + geom_point(aes_string(color="event"), size=2)
    g <- g + facet_wrap(npart ~ height ~ bsize, scales="free",
                        labeller = labeller(npart = as_labeller(npart_labeller),
                                            height = as_labeller(height_labeller),
											bsize = as_labeller(group_size_labeller),
                                            .default=label_both,
                                            .multi_line=FALSE))
    # Set colors.
    g <- g + scale_color_manual(name="Efficiencies",
                                breaks=get_breaks_timings(),
                                labels=get_labels_timings(),
                                values=get_colors_timings())

    # Set title and X/Y labels.
    g <- g + xlab("Number of node")
    g <- g + ylab("Efficiency")
	g <- g + get_theme()

    # Set y-axis range
    g <- g + ylim(c(0.0, 1.10))

    # Save generated plot.
    ggsave(output, g, width=29.7, height=21, units=c("cm"), device=cairo_pdf)
}

compute_efficiency <- function(data, n)
{
    # Select data
    #sdata <- subset(data, npart == n & height == h & model == m)
    sdata <- subset(data, npart == n)


    # Compute task efficiency
	tt_1 <- subset(sdata, event == "task_time" & nnode == 1)
    et <- subset(sdata, event == "task_time")
	if(nrow(tt_1))
	{
		print("There is no one node reference for this algorithm")
	}
    et$efficiency <- tt_1$duration / et$duration

    # Compute scheduling efficiency
    es <- subset(sdata, event == "scheduling_time")
    es$efficiency <- et$duration / (et$duration + es$duration)

    # Compute runtime efficiency
    er <- subset(sdata, event == "runtime_time")
    er$efficiency <- (et$duration + es$duration) / (et$duration + er$duration + es$duration)

    # Compute communication efficiency
    ec <- subset(sdata, event == "communication_time")
    ec$efficiency <- (et$duration + es$duration + er$duration) / (et$duration + er$duration + es$duration + ec$duration)

    # Compute pipeline efficiency
    ep <- subset(sdata, event == "idle_time")
    ep$event <- "pipeline_time" # idle is weird.
    ep$efficiency <- (et$duration + er$duration + es$duration + ec$duration) / (et$duration + er$duration + ep$duration + es$duration + ec$duration)

    # Add new rows for the parallel efficiency
    ndata <- subset(sdata, event == "task_time")
    ndata$event <- "parallel_time"
    sdata <- rbind(sdata, ndata)

    # Compute parallel efficiency
    e <- subset(sdata, event == "parallel_time")
    e$efficiency <- et$efficiency * er$efficiency * ep$efficiency * es$efficiency * ec$efficiency

    # Merge all efficiencies
    sdata <- rbind(et, er, ep, es, ec, e)

    return (sdata)
}

gen_efficiency <- function(data_init, algo_wanted, model_wanted)
{
	data <- subset(data_init, algo == algo_wanted & model == model_wanted)
	if(nrow(data))
	{
		sdata <- NULL
		all_nparts <- unique(data$npart)
		for (i in 1:length(all_nparts))
		{
			if(i == 1)
			{
				sdata <- compute_efficiency(data, all_nparts[i])
			}
			else
			{
				sdata <- rbind(sdata, compute_efficiency(data, all_nparts[i]))
			}
		}
		output <- paste(get_output_directory(), "/", model_wanted, "-", algo_wanted, "-efficiencies.pdf", sep="")
		gen_efficiencies_plot(output, sdata)
	}
}

gen_efficiencies <- function(dbfile)
{
    data_init <- get_data_subset(dbfile, 0L, 0L, "False", get_bsize_reference())
	data_init <- subset(data_init, algo != get_one_node_reference_algorithm())
	data <- melt(data_init, id=c("model", "algo", "nnode", "nthreads", "npart","height","bsize"))

	data <- rename(data, c("variable"="event", "value"="duration"))

	all_algorithm <- unique(data$algo)
	all_model <- unique(data$model)
	for (i in 1:length(all_algorithm))
	{
		for (j in 1:length(all_model))
		{
			gen_efficiency(data, all_algorithm[i], all_model[j])
		}
	}
}

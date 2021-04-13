library(plyr)
library(ggplot2)
library(scales)

gen_group_size_plot_comm <- function(d, model_wanted)
{
	if (is.character(model_wanted))
		d <- subset(d, model == model_wanted)
	d$global_time <- d$global_time/1000
	d <- d[d$communication_vol >= 0,]
	if(nrow(d) > 0)
	{
	g <- ggplot(data=d,aes_string(x="bsize", y="communication_vol", color="algo", group="algo"))
    g <- g + geom_line()
	g <- g + facet_wrap(model ~ npart ~ height ~ nnode, scales="free",
						labeller = labeller(npart = as_labeller(npart_labeller),
											height = as_labeller(height_labeller),
											nnode = as_labeller(nnode_labeller),
											model = as_labeller(data_distribution_labeller),
											.default=label_both,
											.multi_line=FALSE))

    # Set our own colors, linetypes and point shapes.
	g <- g + scale_color_manual(name="Algorithm",
								breaks=get_breaks_runtime(),
								labels=get_labels_runtime(),
								values=get_colors_runtime())

    # Set X/Y labels.
	g <- g + xlab("bsize")
	g <- g + ylab("Communication volume (MB)")
	g <- g + get_theme()

    # Save generated plot.
	if (!exists("output"))
		output <- paste(get_output_directory(), "/", model_wanted, "-bsize-volume.pdf", sep="")
	print(output)
	ggsave(output, g, width=29.7, height=21, units=c("cm"), device=cairo_pdf)
	}
}

gen_group_size_plot_speed <- function(d, model_wanted)
{
	if (is.character(model_wanted))
		d <- subset(d, model == model_wanted)
	d$global_time <- d$global_time/1000
	d <- d[d$global_time > 0,]
	if(nrow(d) > 0)
	{
		g <- ggplot(data=d,aes_string(x="bsize", y="global_time", color="algo", group="algo"))
		g <- g + geom_line()
		g <- g + facet_wrap(model ~ npart ~ height ~ nnode, scales="free",
							labeller = labeller(npart = as_labeller(npart_labeller),
												height = as_labeller(height_labeller),
												nnode = as_labeller(nnode_labeller),
												model = as_labeller(data_distribution_labeller),
												.default=label_both,
												.multi_line=FALSE))

		# Set our own colors, linetypes and point shapes.
		g <- g + scale_color_manual(name="Algorithm",
									breaks=get_breaks_runtime(),
									labels=get_labels_runtime(),
									values=get_colors_runtime())

		# Set X/Y labels.
		g <- g + xlab("bsize")
		g <- g + ylab("Time (s)")
		g <- g + get_theme()

		# Save generated plot.
		if (!exists("output"))
			output <- paste(get_output_directory(), "/", model_wanted, "-bsize-speed.pdf", sep="")
		print(output)
		ggsave(output, g, width=29.7, height=21, units=c("cm"), device=cairo_pdf)
	}
}

gen_group_size <- function(dbfile)
{
    data <- get_data_subset(dbfile, 0L, 0L, "False", 0L)

	data <- subset(data, algo == "explicit" | algo == "implicit")
	tmp <- subset(data, bsize != get_bsize_reference())
	tmp <- unique(tmp[c("nnode", "npart")])
	data <- subset(data, data[c("nnode","npart")] %IN% tmp)

	all_model <- unique(data$model)
	for (i in 1:length(all_model))
	{
		gen_group_size_plot_speed(data, all_model[i])
		gen_group_size_plot_comm(data, all_model[i])
	}
}

library(plyr)
library(ggplot2)
library(scales)

gen_speed_plot <- function(d, model_wanted)
{
	d <- subset(d, model == model_wanted)
	d$global_time <- d$global_time/1000
	g <- ggplot(data=d,aes_string(x="nnode", y="global_time", color="algo"))
    g <- g + geom_line()
	g <- g + facet_wrap(npart ~ height, scales="free",
						labeller = labeller(npart = as_labeller(npart_labeller),
											height = as_labeller(height_labeller),
											.default=label_both,
											.multi_line=FALSE))

    # Set our own colors, linetypes and point shapes.
	g <- g + scale_color_manual(name="Algorithm",
								breaks=get_breaks_runtime(),
								labels=get_labels_runtime(),
								values=get_colors_runtime())

    # Set X/Y labels.
	g <- g + xlab("Number of nodes")
	g <- g + ylab("Time (s)")
	g <- g + get_theme()

    # Save generated plot.
	output <- paste(get_output_directory(), "/", model_wanted, "-speed.pdf", sep="")
	print(output)
	ggsave(output, g, width=29.7, height=21, units=c("cm"), device=cairo_pdf)
}

gen_speed <- function(dbfile)
{
    data <- get_data_subset(dbfile, 0L, 0L, "False", get_bsize_reference())
	data <- subset(data, algo != get_one_node_reference_algorithm())
	data <- subset(data, global_time >= 0)

	all_model <- unique(data$model)
	for (i in 1:length(all_model))
	{
		gen_speed_plot(data, all_model[i])
	}
}


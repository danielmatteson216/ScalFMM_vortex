library(plyr)
library(ggplot2)

gen_comm_plot <- function(db, d_breaks, model_wanted)
{
	db <- subset(db, model == model_wanted)

	db <- db[db$communication_vol > 0,]
	if (nrow(db) > 0)
	{
		g <- ggplot(data=db,aes_string(x="nnode", y="communication_vol", color="algo", group="algo"))
		g <- g + geom_line()
		g <- g + facet_wrap(npart ~ height ~ bsize, scales="free",
							labeller = labeller(npart = as_labeller(npart_labeller),
												height = as_labeller(height_labeller),
												bsize = as_labeller(group_size_labeller),
												.default=label_both,
												.multi_line=FALSE))

		# Set our own colors, linetypes and point shapes.
		g <- g + scale_color_manual(name="Algorithm",
									breaks=get_breaks_runtime(),
									labels=get_labels_runtime(),
									values=get_colors_runtime())

		# Set X/Y labels.
		g <- g + xlab("Number of nodes")
		g <- g + ylab("Volume of Communication (MB)")
		g <- g + get_theme()

		# Save generated plot.
		output <- paste(get_output_directory(), "/", model_wanted, "-comm.pdf", sep="")
		ggsave(output, g, width=29.7, height=21, units=c("cm"), device=cairo_pdf)
	}
}

gen_comm <- function(dbfile)
{
    data <- get_data_subset(dbfile, 0L, 0L, "False", get_bsize_reference())
	data <- subset(data, algo != get_one_node_reference_algorithm() & algo != "simple-mpi" & bsize == get_bsize_reference())
	all_model <- unique(data$model)
	for (i in 1:length(all_model))
	{
		gen_comm_plot(data, unique(data$algo), all_model[i])
	}
}


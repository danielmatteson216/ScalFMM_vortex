library(plyr)
library(ggplot2)

calc_normalized_time <- function(data, ref_name)
{
	# TODO: put starpu algorithm instead
    dataref <- subset(data, algo == get_one_node_reference_algorithm())

    # XXX: Most likely suboptimal but it works as expected!
	data$efficiency <- 0
    for (i in 1:length(ref_name)) {
        for (j in 1:nrow(data)) {
            if (data$algo[j] == ref_name[i]) {
                tmp_ref <- subset(dataref, npart == data$npart[j] &
                                           height == data$height[j])

                seq_time <- subset(tmp_ref, nnode == 1)
                tid = as.integer(as.vector(data$nnode[j]))
				if(nrow(tmp_ref) == 0)
				{
					data$efficiency[j] <- NA
				}
				else if(data$global_time[j] == 0)
				{
					data$efficiency[j] <- 0
				}
				else
				{
					data$efficiency[j] <- seq_time$global_time / (data$global_time[j] * tid)
				}
            }
        }
    }
    return (data)
}

gen_normalized_time_plot <- function(db, d_breaks, model_wanted)
{
	if (is.character(model_wanted))
		db <- subset(db, model == model_wanted)
	#Compute normalized time with one node reference
    db <- calc_normalized_time(db, d_breaks) 
	#Then remove one node reference because it's only available on one node
	db <- subset(db, algo != get_one_node_reference_algorithm())
	db <- db[!(is.na(db$efficiency)),]

	if(nrow(db) > 0)
	{
		g <- ggplot(data=db,aes_string(x="nnode", y="efficiency", color="algo"))
		g <- g + geom_line()
		g <- g + facet_wrap(model ~ npart ~ height, scales="free",
							labeller = labeller(npart = as_labeller(npart_labeller),
												height = as_labeller(height_labeller),
												model = as_labeller(data_distribution_labeller),
												.default=label_both,
												.multi_line=FALSE))

		# Set our own colors, linetypes and point shapes.
		g <- g + scale_color_manual(name="Algorithm",
									breaks=get_breaks_runtime(),
									labels=get_labels_runtime(),
									values=get_colors_runtime())

		# Set X/Y labels.
		g <- g + xlab("Number of nodes")
		g <- g + ylab("Normalized efficiency")
		g <- g + get_theme()


		# Save generated plot.
		if (!exists("output"))
			output <- paste(get_output_directory(), "/", model_wanted, "-normalized-time.pdf", sep="")
		print(output)
		ggsave(output, g, width=29.7, height=21, units=c("cm"), device=cairo_pdf)
	}
}

gen_normalized_time <- function(dbfile)
{
    data <- get_data_subset(dbfile, 0L, 0L, "False", get_bsize_reference())
	data <- subset(data, global_time >= 0)
	all_model <- unique(data$model)
	#Get all algorithm without the reference algorithm
	all_algo <- unique(subset(data, algo != get_one_node_reference_algorithm())$algo)
	for (i in 1:length(all_model))
	{
		gen_normalized_time_plot(data, all_algo, all_model[i])
	}
}


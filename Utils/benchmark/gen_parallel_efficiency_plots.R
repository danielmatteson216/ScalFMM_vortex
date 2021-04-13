library(plyr)
library(ggplot2)

calc_parallel_efficiency <- function(data, ref_name)
{
    # XXX: Most likely suboptimal but it works as expected!
	# NOTE does it really make sense to compare with an mpi version on only one node ?
    for (i in 1:length(ref_name)) {
        data_ref <- subset(data, algo == ref_name[i])

        for (j in 1:nrow(data)) {
            if (data$algo[j] == ref_name[i]) {
                tmp_ref <- subset(data_ref, npart == data$npart[j])
                tmp_ref <- subset(tmp_ref, height == data$height[j])

                seq_time <- subset(tmp_ref, nnode == 1)
                tid = as.integer(as.vector(data$nnode[j]))
				if(nrow(subset(tmp_ref, nnode == 1)) == 0)
				{
					data$efficiency[j] <- NA
				}
				else if(data$global_time[j] == 0)
				{
					data$efficiency[j] <- NA
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

gen_pareff_plot <- function(db, d_breaks, model_wanted)
{
	db <- subset(db, model == model_wanted)
    db <- calc_parallel_efficiency(db, d_breaks)

	db <- db[!(is.na(db$efficiency)),]
	if(nrow(db) > 0)
	{
		g <- ggplot(data=db,aes_string(x="nnode", y="efficiency", color="algo"))
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
		g <- g + ylab("Parallel efficiency")
		g <- g + get_theme()

		# Save generated plot.
		output <- paste(get_output_directory(), "/", model_wanted, "-parallel-efficiency.pdf", sep="")
		ggsave(output, g, width=29.7, height=21, units=c("cm"), device=cairo_pdf)
	}
}

gen_pareff <- function(dbfile)
{
    file  <- paste(dbfile, sep="")
    data <- get_data_subset(dbfile, 0L, 0L, "False", get_bsize_reference())
	data <- subset(data, algo != get_one_node_reference_algorithm())
	data <- subset(data, global_time >= 0)
	all_model <- unique(data$model)
	for (i in 1:length(all_model))
	{
		gen_pareff_plot(data, unique(data$algo), all_model[i])
	}
}

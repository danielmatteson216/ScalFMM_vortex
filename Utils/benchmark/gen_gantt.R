library(plyr)
library(ggplot2)

starts_with <- function(str, pattern)
{
	if(nchar(pattern) == 0)
	{
		return (FALSE)
	}
	return (pattern == substr(str, 1, nchar(pattern)))
}
create_ressource_order <- function(df, nnode, nthreads)
{
	ret <- c()
	if(nnode == 1)
	{
		for(nth in 0:(nthreads-1))
		{
			ret <- c(ret, paste("CPU", nth, sep=""))
		}
	}
	else
	{
		for(nno in 0:(nnode-1))
		{
			for(nth in 0:(nthreads-1))
			{
				ret <- c(ret, paste(nno, "_CPU", nth, sep=""))
			}
		}
	}
	return (ret)
}
set_node_name <- function(df, nnode)
{
	for(i in 1:nrow(df))
	{
		node_number <- substr(df$ResourceId[i], 1, grep("_", df$ResourceId[i]))
		name_origin <- paste("Node ", node_number, sep="")
		df$Origin[i] <- name_origin
	}
}
read_trace <- function(file, nnode, nthreads, start_profiling, stop_profiling)
{
	#Read file
	df<- read.table(file, header=TRUE, sep=",", strip.white=TRUE, fill=TRUE)
	#Remove column I don't need
	df = df[!(names(df) %in% c("Nature","Type", "Depth", "Footprint", "JobId", "Tag", "Params"))]
	#Init origin
	df$Origin <- as.factor("Node 0")

	#Remove uninteresting state (at least, uninteresting for me)
	def_states<-c("Initializing","Deinitializing", "Idle", "Overhead","Nothing","Freeing","Allocating","WritingBack","FetchingInput","PushingOutput","Callback","Progressing","Unpartitioning","AllocatingReuse","Reclaiming","DriverCopy","DriverCopyAsync","Scheduling","Executing", "dplgsy", "execute_on_all_wrapper", "Su")
	df<-df[!(df$Value %in% def_states),]

	#Get only the CPUs event
	df$ResourceId <- as.factor(df$ResourceId)
	resource_order <- create_ressource_order(df, nnode, nthreads)
	df$ResourceId <- factor(df$ResourceId, levels = rev(resource_order))
	#Remove invalide ressources if there is some
	df <- subset(df, !is.na(ResourceId))

	#Get only event between start_profiling and stop_profiling
	df <- subset(df, df$Start > start_profiling)
	df <- subset(df, df$Start < stop_profiling)

	#Set value as a string
	df$Value <- as.character(df$Value)

	#When there is more than one node put the good node on origin
	#Dont do it on one node, because CPU name doesn't start with 0_CPU1, but with CPU1, so all node will be on Node 0
	if(nnode > 1)
	{
		#Reset origin with ResourceId, so the replace will work
		#Other way to do it may exist, but I didn't manage to do it in a less weird way, sorry for that
		df$Origin <- df$ResourceId
		df$Origin <- as.character(df$Origin)
		#For the number of node, replace origin
		for(i in 0:(nnode-1))
		{
			name_origin <- paste("Node ", i, sep="")
			cpu_origin <- paste(i, "_C", sep="")
			df$Origin[starts_with(df$Origin, cpu_origin)] <- name_origin
		}
	}

	#This may change given what you wanna see
	df$Value[starts_with(df$Value, "L2L")] <- "L2L"
	df$Value[starts_with(df$Value, "M2L")] <- "M2L"
	df$Value[starts_with(df$Value, "M2M")] <- "M2M"
	df$Value[starts_with(df$Value, "P2P")] <- "P2P"
	df$Value[df$Value == ""] <- NA

	df$SimpleValue <- df$Value
	df$SimpleValue[df$SimpleValue == "P2M"] <- "Far Field"
	df$SimpleValue[df$SimpleValue == "M2M"] <- "Far Field"
	df$SimpleValue[df$SimpleValue == "M2L"] <- "Far Field"
	df$SimpleValue[df$SimpleValue == "L2L"] <- "Far Field"
	df$SimpleValue[df$SimpleValue == "L2P"] <- "Far Field"
	df$SimpleValue[df$SimpleValue == "P2P"] <- "Near Field"

	#Shift the start, if some event has been removed
	m <- min(df$Start)
	df$Start <- df$Start - m
	df$End <- df$Start+df$Duration

	#Return
	df
} 
gen_simple_gantt_plot <- function(data, model, algo, nnode, npart, bsize, height)
{
	output <- paste(get_output_directory(), "/", model, "-", algo, "-", nnode, "N-", npart/1000000, "M-h", height, "-bs", bsize, "-simple-gantt.pdf", sep="")
	title <- paste(model, " ", algo, " ", nnode, "N ", npart/1000000, "M", " h = ", height, " bs = ", bsize, sep="")
	breaks <- c('Sleeping', 'Far Field', 'Near Field')
	labels <- c('Sleeping', 'Far Field', 'Near Field')
	colors <- c(
		'Sleeping'      = "#f9766e",
		'Far Field'   = "#61f1ff",
		'Near Field'   = "#619dff"
	)
	g <- ggplot(data,aes(x=Start,xend=End, y=factor(ResourceId), yend=factor(ResourceId),color=SimpleValue)) 
	g <- g + theme_bw()
	g <- g + geom_segment(size=8)
	g <- g + ggtitle(title)
	g <- g + ylab("Resource") + xlab("Time [ms]")
	g <- g + scale_color_manual(name="Action", breaks=breaks, labels=labels, values=colors)
	g <- g + facet_wrap(~Origin, ncol=1, scales="free_y") 
	g <- g + scale_y_discrete(breaks=NULL)
    ggsave(output, g, width=29.7, height=21, units=c("cm"), device=cairo_pdf)
}
gen_gantt_plot <- function(data, model, algo, nnode, npart, bsize, height)
{
	output <- paste(get_output_directory(), "/", model, "-", algo, "-", nnode, "N-", npart/1000000, "M-h", height, "-bs", bsize, "-gantt.pdf", sep="")
	title <- paste(model, " ", algo, " ", nnode, "N ", npart/1000000, "M", " h = ", height, " bs = ", bsize, sep="")
	g <- ggplot(data,aes(x=Start,xend=End, y=factor(ResourceId), yend=factor(ResourceId),color=Value)) 
	g <- g + theme_bw()
	g <- g + geom_segment(size=8)
	g <- g + ggtitle(title)
	g <- g + ylab("Resource") + xlab("Time [ms]")
	g <- g + scale_color_brewer(palette="Set1")
	g <- g + facet_wrap(~Origin, ncol=1, scales="free_y") 
	g <- g + scale_y_discrete(breaks=NULL)
    ggsave(output, g, width=29.7, height=42, units=c("cm"), device=cairo_pdf)
}
gen_gantt_grab_data <- function(df)
{
	file <- paste(df$filename[1], sep="")
	data <- read_trace(file, df$nnode[1], df$nthreads[1], df$start_profiling[1], df$stop_profiling[1])
	return (data)
}
gen_gantt <- function(dbfile)
{
	df<- read.table(dbfile, header=TRUE, sep=",", strip.white=TRUE, fill=TRUE,
				  dec=".", colClasses=
		c("factor",	# model
		"factor",	# algorithm
		"integer",	# nnode
		"integer",	# nthread
		"integer",	# npart
		"integer",	# height
		"integer",	# bsize
		"factor",	# filename
		"numeric",	# start_profiling
		"numeric"	# stop_profiling
		))
	all_algorithm <- unique(df$algo)
	all_model <- unique(df$model)
	all_nnode <- unique(df$nnode)
	all_npart <- unique(df$npart)
	all_bsize <- unique(df$bsize)
	all_height <- unique(df$height)
	for (bsi in 1:length(all_bsize))
	{
		for (hei in 1:length(all_height))
		{
			for (mod in 1:length(all_model))
			{
				for (alg in 1:length(all_algorithm))
				{
					for (nno in 1:length(all_nnode))
					{
						for (npa in 1:length(all_npart))
						{
							#Get the corresponding line
							tmp_df <- subset(df, model == all_model[mod] & algo == all_algorithm[alg] & nnode == all_nnode[nno] & npart == all_npart[npa] & bsize == all_bsize[bsi] & height == all_height[hei])
							if(nrow(tmp_df) > 0)
							{
								look <- paste("Grab data ", all_model[mod],all_algorithm[alg], all_nnode[nno], "N", all_npart[npa], all_height[hei], all_bsize[bsi], sep="-")
								#Print something so I know were is the script
								print(look) 
								data <- gen_gantt_grab_data(tmp_df)
								gen_gantt_plot(data, all_model[mod], all_algorithm[alg], all_nnode[nno], all_npart[npa], all_bsize[bsi], all_height[hei])
								gen_simple_gantt_plot(data, all_model[mod], all_algorithm[alg], all_nnode[nno], all_npart[npa], all_bsize[bsi], all_height[hei])
								rm(data)
								gc()
							}
						}
					}

				}
			}
		}
	}
}

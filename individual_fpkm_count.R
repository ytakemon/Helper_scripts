#get FPKM counts of each animal for a list of genes from cufflinks stdout

individual_fpkm_count <- function(cufflinks_dir, sample_list, gene_list){

	library(dplyr) #for using arrange(), but could be used for more in the future when I learn dplyr

	final_df <- data.frame( animal_id = character(), gene_id = character(), FPKM = numeric(), stringsAsFactors = FALSE)

	for (i in 1 : length(sample_list)){
		
		data <- read.delim(file = paste0(cufflinks_dir, sample_list[i], "/genes.fpkm_tracking"), sep = "", header = TRUE, stringsAsFactors = FALSE)
		data <- subset(data, !duplicated(gene_id)) #tracking ID contains the gene names	
		data <- data[data$gene_id %in% gene_list,]

		if (dim(data)[1] != length(gene_list)){

			warning(paste( "Some genes in gene_list were not found in", sample_list[1], sep = " "))

		}

		data$animal_id <- sample_list[i]
		data <- data[,c("animal_id", "gene_id", "FPKM")] # reorder to match final_df arrangement
		final_df <- rbind(final_df, data)

	}

	final_df <- arrange(final_df, gene_id)
	rownames(final_df) <- NULL #reset rownames numbering
	return(final_df)

}



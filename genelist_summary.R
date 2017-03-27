genelist_summary <- function(gene_list, list_dir){
	source("Helper_scripts/STDEV_SE_NORM_functions.R") 

	UCSC_head <- c("tracking_id", "class_code", "nearest_ref_id", 
				"gene_id", "gene_short_name", "tss_id", "locus",
				"length", "coverage", "FPKM", "FPKM_conf_lo",
				"FPKM_conf_hi", "FPKM_status")

	#establish base data frame with first gene
	file_WT <- read.delim(paste(list_dir, "/WT_",gene_list[1],"_fpkm.txt", sep = ""), sep = "\t", header = FALSE)
	file_KO <- read.delim(paste(list_dir, "KO_",gene_list[1],"_fpkm.txt", sep = ""), sep = "\t", header = FALSE)
	colnames(file_WT) <- UCSC_head
	colnames(file_KO) <- UCSC_head
	file_comp <- file_WT[,c("tracking_id", "FPKM")]
	file_comp <- rbind(file_comp, file_KO[,c("tracking_id", "FPKM")])
	file_comp$geno <- "WT"
	file_comp$geno[(dim(file_WT)[1] + 1) : (dim(file_comp)[1])] <- "KO"
	file_comp$gene <- gene_list[1]
	summary_data <- file_comp #will be appending to summary_data in loop starting from the second gene on the list

	for (i in 2:length(gene_list)){ #append to summary_data created by the second gene
		
		file_WT <- read.delim(paste(list_dir, "WT_",gene_list[i],"_fpkm.txt", sep = ""), sep = "\t", header = FALSE)
		file_KO <- read.delim(paste(list_dir, "KO_",gene_list[i],"_fpkm.txt", sep = ""), sep = "\t", header = FALSE)
		colnames(file_WT) <- UCSC_head
		colnames(file_KO) <- UCSC_head
		file_comp <- file_WT[,c("tracking_id", "FPKM")]
		file_comp <- rbind(file_comp, file_KO[,c("tracking_id", "FPKM")])
		file_comp$geno <- "WT"
		file_comp$geno[(dim(file_WT)[1] + 1) : (dim(file_comp)[1])] <- "KO"
		file_comp$gene <- gene_list[i]
		summary_data <- rbind(summary_data, file_comp)

	}

	summary_data <- summarySEwithin(summary_data, measurevar = "FPKM", withinvars = c("geno", "gene"))
	summary_data$FPKMnorm_sub_SD <- (summary_data$FPKM_norm - summary_data$sd)
	return(summary_data)
}

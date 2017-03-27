#currenlty can only handle upto 6 samples

cufflinks_fpkm_df <- function(cufflinks_dir, sample_list){

	rm( list = ls(pattern = "data_"))

	for (i in 1 : length(sample_list)){
		
		data <- read.delim(file = paste0(cufflinks_dir, sample_list[i], "/genes.fpkm_tracking"), sep = "", header = TRUE, stringsAsFactors = FALSE)
		data <- subset(data, !duplicated(tracking_id)) #tracking_id contains gene name

		if (any (grep("CUFF.", data$tracking_id) == TRUE )){

			data <- data[-grep("CUFF.", data$tracking_id), ]
			rownames(data) <- make.names(data$tracking_id)
			data <- data[order(rownames(data)),]		
			assign(paste0("data_", i ), data)

		} else {

			rownames(data) <- make.names(data$tracking_id)
			data <- data[order(rownames(data)),]		
			assign(paste0("data_", i ), data)
		}
	}
	
	list <- ls(pattern = "data_")

	if (length(sample_list) == 2){
		
		common <- intersect(data_1$tracking_id, data_2$tracking_id)
		data_1 <- data_1[common, ]
		data_2 <- data_2[common, ]

		FPKM <- array(0, c(dim(data_1)[1], length(sample_list)), dimnames = list( common, sample_list))
		FPKM[,1] <- data_1$FPKM
		FPKM[,2] <- data_2$FPKM
		FPKM <- FPKM[complete.cases(FPKM),]
		return(FPKM)

	} else if (length(sample_list) == 3){
		common <- intersect(intersect(data_1$tracking_id, data_2$tracking_id), data_3$tracking_id)
		data_1 <- data_1[common, ]
		data_2 <- data_2[common, ]
		data_3 <- data_3[common, ]

		FPKM <- array(0, c(dim(data_1)[1], length(sample_list)), dimnames = list( common, sample_list))
		FPKM[,1] <- data_1$FPKM
		FPKM[,2] <- data_2$FPKM
		FPKM[,3] <- data_3$FPKM
		FPKM <- FPKM[complete.cases(FPKM),]
		return(FPKM)
		
	} else if (length(sample_list) == 4){
		common <- intersect(intersect(intersect(data_1$tracking_id, data_2$tracking_id), data_3$tracking_id), data_4$tracking_id)
		data_1 <- data_1[common, ]
		data_2 <- data_2[common, ]
		data_3 <- data_3[common, ]
		data_4 <- data_4[common, ]

		FPKM <- array(0, c(dim(data_1)[1], length(sample_list)), dimnames = list( common, sample_list))
		FPKM[,1] <- data_1$FPKM
		FPKM[,2] <- data_2$FPKM
		FPKM[,3] <- data_3$FPKM
		FPKM[,4] <- data_4$FPKM
		FPKM <- FPKM[complete.cases(FPKM),]
		return(FPKM)

	}else if (length(sample_list) == 5){
		common <- intersect(intersect(intersect(intersect(data_1$tracking_id, data_2$tracking_id), data_3$tracking_id), data_4$tracking_id), data_5$tracking_id)
		data_1 <- data_1[common, ]
		data_2 <- data_2[common, ]
		data_3 <- data_3[common, ]
		data_4 <- data_4[common, ]
		data_5 <- data_5[common, ]

		FPKM <- array(0, c(dim(data_1)[1], length(sample_list)), dimnames = list( common, sample_list))
		FPKM[,1] <- data_1$FPKM
		FPKM[,2] <- data_2$FPKM
		FPKM[,3] <- data_3$FPKM
		FPKM[,4] <- data_4$FPKM
		FPKM[,5] <- data_5$FPKM
		FPKM <- FPKM[complete.cases(FPKM),]
		return(FPKM)

	}else if (length(sample_list) == 6){
		common <- intersect(intersect(intersect(intersect(intersect(data_1$tracking_id, data_2$tracking_id), data_3$tracking_id), data_4$tracking_id), data_5$tracking_id), data_6$tracking_id)
		data_1 <- data_1[common, ]
		data_2 <- data_2[common, ]
		data_3 <- data_3[common, ]
		data_4 <- data_4[common, ]
		data_5 <- data_5[common, ]
		data_6 <- data_6[common, ]

		FPKM <- array(0, c(dim(data_1)[1], length(sample_list)), dimnames = list( common, sample_list))
		FPKM[,1] <- data_1$FPKM
		FPKM[,2] <- data_2$FPKM
		FPKM[,3] <- data_3$FPKM
		FPKM[,4] <- data_4$FPKM
		FPKM[,5] <- data_5$FPKM
		FPKM[,6] <- data_6$FPKM
		FPKM <- FPKM[complete.cases(FPKM),]
		return(FPKM)
	}
}

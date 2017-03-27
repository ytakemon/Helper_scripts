cor_cufflink_fpkm <- function(sample_list, cufflinks_dir){

        original_name <- deparse(substitute(sample_list))
        len <- length(sample_list)
        rm( list = ls(pattern = "data_"))

        for (i in 1 : len){

                data <- read.delim(file = paste0(cufflinks_dir, sample_list[i], "/genes.fpkm_tracking"), sep = "", header = TRUE, stringsAsFactors = FALSE)
                data <- data[-grep("CUFF.", data$tracking_id), ]

                if (any( data$tss_id == "-")){

                        data <- data[-grep("-", data$tss_id), ]
                        rownames(data) <- make.names(data$tss_id)
                        data <- data[order(rownames(data)),]
                        assign(paste0("data_", i ), data)

                } else {

                        rownames(data) <- make.names(data$tss_id)
                        data <- data[order(rownames(data)),]
                        assign(paste0("data_", i ), data)
                }
        }

        list <- ls(pattern = "data_")

        if (len == 2){

                common <- intersect(data_1$tss_id, data_2$tss_id)
                data_1 <- data_1[common, ]
                data_2 <- data_2[common, ]

                FPKM <- array(0, c(dim(data_1)[1], len), dimnames = list( common, sample_list))
                FPKM[,1] <- data_1$FPKM
                FPKM[,2] <- data_2$FPKM
                FPKM <- FPKM[complete.cases(FPKM),]
                return(cor(FPKM))

        } else if (len == 3){
                common <- intersect(intersect(data_1$tss_id, data_2$tss_id), data_3$tss_id)
                data_1 <- data_1[common, ]
                data_2 <- data_2[common, ]
                data_3 <- data_3[common, ]

                FPKM <- array(0, c(dim(data_1)[1], len), dimnames = list( common, sample_list))
                FPKM[,1] <- data_1$FPKM
                FPKM[,2] <- data_2$FPKM
                FPKM[,3] <- data_3$FPKM
                FPKM <- FPKM[complete.cases(FPKM),]
                return(cor(FPKM))

        } else if (len == 4){
                common <- intersect(intersect(intersect(data_1$tss_id, data_2$tss_id), data_3$tss_id), data_4$tss_id)
                data_1 <- data_1[common, ]
                data_2 <- data_2[common, ]
                data_3 <- data_3[common, ]
                data_4 <- data_4[common, ]

                FPKM <- array(0, c(dim(data_1)[1], len), dimnames = list( common, sample_list))
                FPKM[,1] <- data_1$FPKM
                FPKM[,2] <- data_2$FPKM
                FPKM[,3] <- data_3$FPKM
                FPKM[,4] <- data_4$FPKM
                FPKM <- FPKM[complete.cases(FPKM),]
                return(cor(FPKM))

        }
}

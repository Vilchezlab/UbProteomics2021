
options(java.parameters = "-Xmx56000m")
install.packages("gplots")
install.packages("psych")
install.packages("colorRamps")
install.packages("matrixStats")
install.packages("ggplot2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pcaMethods")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute")

install.packages("tmvtnorm",method="curl")
install.packages("norm",method="curl")
install.packages("devtools")
install.packages("stringr")
library("stringr")
library("gplots")
library("psych")
library('colorRamps')
library(matrixStats)
library(ggplot2)
library(dplyr)
library(stringr)
install.packages("usethis")
library("devtools")
install_github("vqv/ggbiplot")
library(ggbiplot)
library(stringr)
library(ggfortify)

install.packages("mvtnorm")
install.packages("Matrix")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("stats4")

install.packages("stats4")
install.packages("gmm")
install.packages("sandwich")
install.packages("norm")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("preprocessCore")

library("preprocessCore")
install.packages("xlsx")

Sys.setenv(JAVA_HOME='C:/Program Files/Java/jre1.8.0_221')
library("xlsx")
library("tibble")
library(reshape2) 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute")
library("imputeLCMD")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("genefilter")

library("genefilter")

options(download.file.method = "curl")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pcaMethods")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("preprocessCore")
install.packages("survMisc")
library(survMisc)
Impute_Matrix <- function(subset,group_name) {
  #generates imputed matrix of given matrix
  
  subset <- apply(subset,c(1,2), function(x){as.numeric(x)})
  
  colnames <- colnames(subset)
  
  model <- model.Selector(subset)
  #	model_count <- table(model[1])
  one_index <- which(unlist(model[1]) %in% 1)
  zero_index <- which(unlist(model[1]) %in% 0)
  
  
  imputed_matrix <- impute.MAR.MNAR(subset,model[1], method.MAR = "MLE", method.MNAR = "QRILC")
  colnames(imputed_matrix) <- colnames 
  
  #part of subset that contains NAs
  #Missing at random
  MAR_imputed_subset <- subset[one_index,]
  a <- MAR_imputed_subset[rowSums(is.na(MAR_imputed_subset)) > 0,]
  a_imput <- imputed_matrix[rownames(a),]
  a_sd_imput <- rowSds(a_imput)
  a_means_imput <- rowMeans(a_imput)
  
  #Missing NOT at random
  MNAs_imputed_subset <- subset[zero_index,]
  zero <- MNAs_imputed_subset[rowSums(is.na(MNAs_imputed_subset)) > 0,]
  zero_imput <- imputed_matrix[rownames(zero),]
  zero_sd_imput <- rowSds(zero_imput)
  zero_means_imput <- rowMeans(zero_imput)
  
  #part of subset with no NAs
  b <- subset[rowSums(is.na(subset)) == 0,]
  b_mean <- rowMeans(b)
  b_sd <- rowSds(b)
  
  return(list(imputed_matrix=imputed_matrix))
  #	return(list(imputed_matrix=imputed_matrix,number_MNAS=number_MNAS_matrix,number_MAR=number_MAR_matrix,NA_count=interval_matrix,model=model,plot=p))
}





Z_normalization <- function(log2_intensity_matrix_use) {
  # Function to do a global Z-normalization
  
  row_names <- colnames(log2_intensity_matrix_use)
  
  matrix_mean <- mean(as.matrix(log2_intensity_matrix_use),na.rm=TRUE)
  matrix_standard_dev <- sd(as.matrix(log2_intensity_matrix_use),na.rm=TRUE)
  
  Z_standard_matrix <- vector()
  for (col_id in 1:ncol(log2_intensity_matrix_use)) {
    
    col_vector <- log2_intensity_matrix_use[,col_id]
    col_mean <- mean(col_vector,na.rm=TRUE)
    col_stdev <- sd(col_vector,na.rm=TRUE)
    
    ##Z-normalize every column separately
    new_vector <- (col_vector-col_mean)/col_stdev
    Z_standard_matrix <-cbind(Z_standard_matrix,new_vector)
  }
  
  ## Back transformation of normalized data into value range of original data
  Z_standard_matrix_back <- Z_standard_matrix * matrix_standard_dev + matrix_mean
  colnames(Z_standard_matrix_back) <- row_names
  return(Z_standard_matrix_back)
}

orig_data <- read.csv("GlyGly (K)Sites.txt",sep = "\t",stringsAsFactors=FALSE,header=TRUE,blank.lines.skip = TRUE)

intensity_data <- orig_data[grep("Intensity",colnames(orig_data))]
intensity_data <- intensity_data[,-1]
filter_index <- grep('___',colnames(intensity_data))
intensity_data_take <- intensity_data[,-filter_index]
colnames <- str_replace(colnames(intensity_data_take),"Intensity.","")
colnames(intensity_data_take) <- colnames

intensity_data_annotate <- cbind(orig_data["Reverse"],orig_data["Potential.contaminant"],orig_data["Proteins"],orig_data["Leading.proteins"],orig_data["Protein.names"],orig_data["Gene.names"],orig_data["Positions.within.proteins"],orig_data["Localization.prob"],intensity_data_take)

intensity_data_annotate_filtered <- intensity_data_annotate[rowSums(cbind(orig_data["Reverse"],orig_data["Potential.contaminant"])=="+")<1,]

intensity_data_use <- intensity_data_annotate_filtered[,9:ncol(intensity_data_annotate_filtered)]

intensity_data_annotation <- intensity_data_annotate_filtered[,3:8]

Non_zero_counts_perROW <- rowSums(intensity_data_use != 0)

Non_Zero_counts_perCOLUMN <- colSums(intensity_data_use != 0)
write.table(as.matrix(Non_Zero_counts_perCOLUMN),file="Identified_GlyGly_sites.txt",quote=FALSE,col.names=FALSE,sep="\t")
intensity_data_use_counts <- rbind(intensity_data_use,Non_Zero_counts_perCOLUMN)

intensity_data_use_counts_COLUMN_filtered <- intensity_data_use_counts[,intensity_data_use_counts[nrow(intensity_data_use_counts),] >= 100][1:nrow(intensity_data_use),]

log2_intensity_data_use_counts_COLUMN_filtered <- log2(intensity_data_use_counts_COLUMN_filtered)
log2_intensity_data_use_counts_COLUMN_filtered[log2_intensity_data_use_counts_COLUMN_filtered == -Inf] <- 0


Zero_counts_perROW <- rowSums(intensity_data_use_counts_COLUMN_filtered == 0)
intensity_data_use_counts_COLUMN_filtered_n <- cbind(intensity_data_annotation,intensity_data_use_counts_COLUMN_filtered,Zero_counts_perROW)
intensity_data_use_counts_COLUMN_filtered_n <- data.frame(intensity_data_use_counts_COLUMN_filtered_n)
intensity_data_use_counts_COLUMNROW_filtered_annot <- intensity_data_use_counts_COLUMN_filtered_n[intensity_data_use_counts_COLUMN_filtered_n$Zero_counts_perROW < ncol(intensity_data_use_counts_COLUMN_filtered)-2,]

intensity_data_use_counts_COLUMNROW_filtered <- intensity_data_use_counts_COLUMNROW_filtered_annot[,7:(ncol(intensity_data_use_counts_COLUMN_filtered_n)-1)]
intensity_data_annotation_filtered <- intensity_data_use_counts_COLUMNROW_filtered_annot[,1:6]

log2_intensity_automaticFiltered <- log2(intensity_data_use_counts_COLUMNROW_filtered)
log2_intensity_automaticFiltered <- as.matrix(replace(log2_intensity_automaticFiltered, log2_intensity_automaticFiltered==-Inf, 0))


NA_log2_intensity_automaticFiltered <- as.matrix(replace(log2_intensity_automaticFiltered,log2_intensity_automaticFiltered==0,NA))
Z_standard_matrix_back <- Z_normalization(NA_log2_intensity_automaticFiltered)

Z_standard_matrix_back_annotate <- cbind(intensity_data_annotation_filtered,Z_standard_matrix_back)
write.table(Z_standard_matrix_back_annotate,file="Z-Normalized_Intensities_all.txt",quote=FALSE,row.names=FALSE,sep="\t")

ind <- grep("G8JY83",as.vector(as.matrix(Z_standard_matrix_back_annotate["Proteins"])))
Ubi_rows <- Z_standard_matrix_back_annotate[ind,]

Ubi_rows[1,c(6:ncol(Ubi_rows))]


Ubi_rows_col <- Ubi_rows[,6:ncol(Ubi_rows)]
Ubi_rows_col[is.na(Ubi_rows_col)] <- 0

Groups <- str_replace(colnames(Z_standard_matrix_back),"_[0-9][0-9]$","")
groups <- as.matrix(cbind(colnames(Z_standard_matrix_back),Groups))
rownames(groups) <- groups[,1]
groups <- as.matrix(groups[,-1])
colnames(groups) <- "Groups"
write.table(groups,"Groups.txt",quote=FALSE,sep="\t",col.names=FALSE)
group_header <- str_replace(colnames(Z_standard_matrix_back),"_[0-9][0-9]$","")
anova_input <- Z_standard_matrix_back
new_matrix <- vector()
filtered_matrix_notUse <- vector()

row_names_new <- vector()
row_names_new_NOTuse <- vector()

intensity_data_annotation_filtered_tukey <- vector()

index_vector <- vector()
Z_standard_matrix_filtered <- vector()
p_value_matrix_merged <- vector()
start <- 1
for (row_index in 1:nrow(anova_input)) {
  #for (row_index in 1:10) {
  
  
  
  #	print(row_index)
  row_values <- as.vector(anova_input[row_index,])
  #	print(anova_input[row_index,])
  
  k <- row_values
  names(k) <- group_header
  
  #test if anova makes sense based on the data
  a <- aggregate(row_values ~ group_header, data=as.matrix(k), mean, na.action = na.omit)
  
  #	aggregate(row_values ~ group_header, data=as.matrix(k), function(x) c(mean = mean(x), sd = sd(x)))
  
  if (nrow(a) > 1) {
    
    anova_output <- aov(row_values ~ group_header,na.action=na.exclude)
    #		print(summary(anova_output))
    
    annotation_row <- intensity_data_annotation_filtered[row_index,]
    
    if (length(summary(anova_output)[[1]]) == 5) {
      
      #			print(length(summary(anova_output)[[1]]))
      p_value <- summary(anova_output)[[1]][5][[1]][1]
      
      if (p_value < 0.5) {
        
        index_vector <- c(index_vector,row_index)
        
        tukey_output <- TukeyHSD(anova_output)
        #				print(tukey_output)
        
        Z_standard_matrix_filtered <- rbind(Z_standard_matrix_filtered,row_values)
        
        p_value_matrix <- as.matrix(tukey_output$group_header[,4])
        rownames(p_value_matrix) <- rownames(tukey_output$group_header)
        
        if (start != 0) {
          
          p_value_matrix_merged <- p_value_matrix
        } else {
          #					print("_______________")
          #					print(p_value_matrix_merged)
          #					print(p_value_matrix)
          m <- merge(p_value_matrix_merged,p_value_matrix,by="row.names",all=TRUE)
          
          p_value_matrix_merged <- m
          rownames(p_value_matrix_merged) <- m[,1]
          p_value_matrix_merged <- p_value_matrix_merged[,-1]
          #					print(p_value_matrix_merged)
          
          #					print("....................")
        }
        
        
        intensity_data_annotation_filtered_tukey <- rbind(intensity_data_annotation_filtered_tukey,annotation_row)
        
        #		row_names_new <- c(row_names_new,row_names[row_index])
        #		new_matrix <- rbind(new_matrix,row_values)
        #		tukey_use <- cbind(tukey_use,tukey_output$group_header[,4])
        
        #	} else {
        #		row_names_new_NOTuse <- c(row_names_new_NOTuse,row_names[row_index])
        #		filtered_matrix_notUse <- rbind(filtered_matrix_notUse,row_values)
        #		tukey_NOTuse <- cbind(tukey_NOTuse,tukey_output$group_header[,4])
        
        
        start <- 0
      }
    }
  }
}


colnames(Z_standard_matrix_filtered) <- colnames(anova_input)

tukey_filtered <- cbind(intensity_data_annotation_filtered_tukey,t(p_value_matrix_merged))
Z_standard_matrix_filtered_annotated <- cbind(intensity_data_annotation_filtered_tukey,Z_standard_matrix_filtered)

write.table(tukey_filtered,"Tukey_results.txt",quote=FALSE,row.names=FALSE,sep="\t")
write.table(Z_standard_matrix_filtered_annotated,"ZnormalizedData_filtered.txt",quote=FALSE,row.names=FALSE,sep="\t")


groups_NOTuse <- read.csv("GlyGly_Samples_NOTuse.txt",sep = "\t",header=FALSE,check.names = FALSE)
NOTuse_sample_list <- list()
for (a in 1:nrow(groups_NOTuse)) {
  
  line <- groups_NOTuse[a,]
  if (is.null(NOTuse_sample_list[[as.vector(line[,2])]])) {
    
    NOTuse_sample_list[[as.vector(line[,2])]] <- as.vector(line[,1])
  } else {
    
    sample_vector <- NOTuse_sample_list[[as.vector(line[,2])]]
    sample_vector <- c(sample_vector,as.vector(line[,1]))
    NOTuse_sample_list[[as.vector(line[,2])]] <- sample_vector
  }
}
all_matrix_norm <- t(Z_standard_matrix_back)
matrix_names_norm <- unique(groups[,1])
s<- merge(groups,all_matrix_norm,by="row.names")
matrix_data_list_norm <- list()
imputed_matrix_data_list <- list()
imputed_matrix <- vector()
matrix_sample_filtered <- vector()
NA_count_vector_list <- list()


for (groupname in matrix_names_norm) {
  
  if (groupname %in% s$Groups) {
    matrix_data_norm <- subset(s,s$Groups==groupname)
    
    if (! is.null(NOTuse_sample_list[[groupname]])) {
      
      sample_vector_NOTuse <- as.vector(NOTuse_sample_list[[groupname]])
      
      n <- 1
      while (n <= length(sample_vector_NOTuse)) {
        sample_name <- sample_vector_NOTuse[n]
        sample_index <- grep(sample_name,t(matrix_data_norm)[1,])
        matrix_data_norm <- t(t(matrix_data_norm)[,-(sample_index)])
        n <- n + 1
      }
      
    }
    #		print(ncol(t(matrix_data_norm)))
    
    if ( ncol(t(matrix_data_norm)) > 0) {
      
      matr_for_imput <- t(matrix_data_norm)
      print(head(matr_for_imput))
      colnames <- matr_for_imput[1,]
      matr_for_imput <- matr_for_imput[-c(1,2),]
      
      t_matrix_data_NA <- apply(matr_for_imput, MARGIN = 1, FUN = function(x) length(x[is.na(x)]) )
      NA_count_vector_list[[groupname]] <- t_matrix_data_NA
      
      colnames(matr_for_imput) <- colnames
      print(head(matr_for_imput))
      
      Impute_output <- Impute_Matrix(matr_for_imput,groupname)
      imputed_matrix_data <- t(Impute_output[[1]])
      
      #		print(head(Impute_output[[1]]))
      
      imputed_matrix <- cbind(imputed_matrix,Impute_output[[1]])
      
      
      imputed_matrix_data_list[[groupname]] <- imputed_matrix_data
      
      matrix_sample_filtered <- cbind(matrix_sample_filtered,matr_for_imput)
      matrix_data_list_norm[[groupname]] <- matrix_data_norm
    }
    
  }
}


imputed_matrix <- apply(imputed_matrix,c(1,2), function(x){as.numeric(x)})
write.table(cbind(intensity_data_annotation_filtered,imputed_matrix),"Z-Normalized_Data_matrix_imputed.txt",quote=FALSE,sep="\t",row.names=FALSE)

matrix_sample_filtered <- apply(matrix_sample_filtered,c(1,2), function(x){as.numeric(x)})
write.table(cbind(intensity_data_annotation_filtered,matrix_sample_filtered),"Z-Normalized_Data_matrix.txt",quote=FALSE,sep="\t",row.names=FALSE)


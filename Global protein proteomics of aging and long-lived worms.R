#! /usr/bin/Rscript --vanilla 

options(java.parameters = "-Xmx66000m")


library("gplots")
library("psych")
library('colorRamps')
library('stringr')
library("xlsx")
library('matrixStats')
library("ggfortify")
library("plyr")
library("imputeLCMD")

library("tidyverse")
library("tibble")



#Impute_subsets <- function(matrix,subset_index,group_name) {

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




#Function filters the FC-matrix according to the given number of allowed NAs in each row

Filter_matrix_for_NA_values <- function(FC_matrix,annotation_matrix,NA_threshold) {

	NA_counts_FC_vector <- apply(FC_matrix, MARGIN = 1, FUN = function(x) length(x[is.na(x)]) )

	NA_FC_matrix <- cbind(FC_matrix,NA_counts_FC_vector,annotation_matrix)
	NA_Count_index <- grep("NA_counts_FC_vector",colnames(NA_FC_matrix))

	#filter proteins which have NA_threshold missing values at maximum
	Matrix_filtered <- NA_FC_matrix[as.numeric(NA_FC_matrix[,NA_Count_index]) <= NA_threshold,]

	FC_matrix_filtered <- Matrix_filtered[,1:ncol(FC_matrix)]
	FC_matrix_filtered <- apply(FC_matrix_filtered,c(1,2), function(x){as.numeric(x)})

	NA_counts_filtered <- as.numeric(Matrix_filtered[,ncol(FC_matrix)+1])
	filtered_annotation_matrix <- Matrix_filtered[,(ncol(FC_matrix)+2):ncol(Matrix_filtered)]


	colnames(FC_matrix_filtered) <- colnames(FC_matrix)

	FC_matrix_filtered_heatmap <- FC_matrix_filtered
	FC_matrix_filtered_heatmap[is.na(FC_matrix_filtered_heatmap)] <- 0


	return(list(FC_matrix_filtered=FC_matrix_filtered,NA_counts_filtered=NA_counts_filtered,annotation=filtered_annotation_matrix))
}









#setwd('/nas/nas9_prot/MS-DataAnalysis/analyzed_data/2018/773_1316')



#R_script <- "ProteomeComparison_V4_Seda_773_1316.R"

#info <- paste("The project was analysed with script: \n /home/corinnaklein/Arbeitsfläche/Corinna/Proteomics/Projects/R_scripts/",R_script,sep="")
#write.table(info,'R_whatIdid.log',sep = "\t",quote=FALSE)


report <- read_csv('report_againUSE.csv')
report[report == "Filtered"] <- NA


annotation <- cbind(report$PG.Genes,report$PG.ProteinAccessions,report$PG.ProteinDescriptions)
colnames(annotation) <- c("Genes","ProteinAccessions","ProteinDescriptions")

PG.Quantity_index <- grep('PG.Quantity',colnames(report))
PG_Quantity_data <- as.matrix(report[,PG.Quantity_index])
PG_Quantity_data_log2 <- log2(apply(PG_Quantity_data,c(1,2), function(x){as.numeric(x)}))



new_names <- str_replace_all(colnames(report[,PG.Quantity_index]),".PG.Quantity","")
#new_names <- str_replace_all(new_names,"\\[[0-9]+\\]","log2")
new_names <- str_replace_all(new_names,"\\[[0-9]+\\] ","")


colnames(PG_Quantity_data_log2) <- new_names
report[,PG.Quantity_index] <- PG_Quantity_data_log2


lfq <- PG_Quantity_data_log2



# groups_NOTuse <- read.csv('Group_NOTuse.txt',sep = "\t",header=FALSE,check.names = FALSE)
# groups_NOTuse <- list()
# NOTuse_sample_list <- list()
# for (a in 1:nrow(groups_NOTuse)) {
# 
# 	line <- groups_NOTuse[a,]
# 	if (is.null(NOTuse_sample_list[[as.vector(line[,2])]])) {
# 
# 		NOTuse_sample_list[[as.vector(line[,2])]] <- as.vector(line[,1])
# 	} else {
# 
# 		sample_vector <- NOTuse_sample_list[[as.vector(line[,2])]]
# 		sample_vector <- c(sample_vector,as.vector(line[,1]))
# 		NOTuse_sample_list[[as.vector(line[,2])]] <- sample_vector
# 	}
# }
NOTuse_sample_list <- c()

Groups <- str_replace(colnames(lfq),"_[0-9][0-9]$","")
Groups<-str_replace(Groups,"DMSO_[0-9][0-9]_[0-9]$","DMSO")
groups <- as.matrix(cbind(colnames(lfq),Groups))


rownames(groups) <- groups[,1]
groups <- as.matrix(groups[,-1])
colnames(groups) <- "Groups"
write.table(groups,"Groups.txt",quote=FALSE,sep="\t",col.names=FALSE)
group_header <- str_replace(colnames(lfq),"_[0-9][0-9]$","")
group_header <- str_replace(group_header,"DMSO_[0-9][0-9]_[0-9]$","DMSO")



groups <- read.csv('Groups.txt',sep = "\t",header=FALSE,check.names = FALSE,row.names=1)
colnames(groups) <- c("Groups")
matrix_names <- unique(groups[,1])


all_matrix <- t(lfq)
s<- merge(groups,all_matrix,by="row.names")

imputed_matrix <- vector()
matrix_data_list <- list()
imputed_matrix_data_list <- list()
matrix_sample_filtered <- vector()

for (groupname in matrix_names) {
	matrix_data <- subset(s,s$Groups==groupname)


	if (! is.null(NOTuse_sample_list[[groupname]])) {

		sample_vector_NOTuse <- as.vector(NOTuse_sample_list[[groupname]])

		n <- 1
		while (n <= length(sample_vector_NOTuse)) {
			sample_name <- sample_vector_NOTuse[n]
			sample_index <- grep(sample_name,t(matrix_data)[1,])
			matrix_data <- t(t(matrix_data)[,-(sample_index)])
			n <- n + 1
		}
	}

	matr_for_imput <- t(matrix_data)

	colnames <- matr_for_imput[1,]
	matr_for_imput <- matr_for_imput[-c(1,2),]

	colnames(matr_for_imput) <- colnames

	Impute_output <- Impute_Matrix(matr_for_imput,groupname)
	imputed_matrix_data <- t(Impute_output[[1]])

	print(head(Impute_output[[1]]))

	imputed_matrix <- cbind(imputed_matrix,Impute_output[[1]])


	imputed_matrix_data_list[[groupname]] <- imputed_matrix_data

	matrix_data_list[[groupname]] <- matrix_data

	matrix_sample_filtered <- cbind(matrix_sample_filtered,matr_for_imput)
}


groups_NOTuse <- read.csv('Group_NOTuse.txt',sep = "\t",header=FALSE,check.names = FALSE)

groupsToExclude<- unique(groups_NOTuse[,2])
groupsToExcludelength<-c(1:length(groupsToExclude))
excludeIndex <- c()
for (x in groupsToExcludelength){
  excludeIndex=cbind(excludeIndex,grep(groupsToExclude[x],colnames(matrix_sample_filtered)))
}

matrix_sample_filtered <- matrix_sample_filtered[,-excludeIndex]


matrix_sample_filtered <- apply(matrix_sample_filtered,c(1,2), function(x){as.numeric(x)})

imputed_matrix <- imputed_matrix[,-excludeIndex]
imputed_matrix <- apply(imputed_matrix,c(1,2), function(x){as.numeric(x)})


write.table(cbind(annotation,imputed_matrix),"Data_matrix_imputed.txt",quote=FALSE,sep="\t",row.names=FALSE)


write.table(cbind(annotation,matrix_sample_filtered),"Data_matrix.txt",quote=FALSE,sep="\t",row.names=FALSE)



#matrix_sample_filtered_zero <- matrix_sample_filtered

#matrix_sample_filtered_zero[is.na(matrix_sample_filtered_zero)] <- 0


#pdf("Heatmap_filteredData.pdf")
#heatmap.2(as.matrix(matrix_sample_filtered_zero),Rowv = TRUE,Colv=TRUE,distfun = dist,hclustfun = hclust,dendrogram="col",col=colorRampPalette(c("white","red","grey","blue")),keysize =1,trace="none",cexCol = 0.8,margins = c(8, 6),main="Heatmap imputed log2 intensities\n at least 3 valid values per protein",labRow = "")
#dev.off()

#PCA_generation_new1(matrix_sample_filtered,groups,"PCA_filtered.pdf","PCA")






pdf("Heatmap_imputed_n.pdf")
heatmap.2(as.matrix(imputed_matrix),Rowv = TRUE,Colv=TRUE,distfun = dist,hclustfun = hclust,dendrogram="col",col=colorRampPalette(c("white","red","grey","blue")),keysize =1,trace="none",cexCol = 0.45,margins = c(8, 2),main="Heatmap imputed log2 intensities\n at least 3 valid values per protein",labRow = "")
dev.off()

PCA_generation_new1(imputed_matrix,groups,"PCA_imputed_n.pdf","PCA")


compare <- as.matrix(compare)


diff_threshold <- 1
FDR_threshold <- 0.05

pdf("Vulcano_plot_n.pdf")
result_matrix <- matrix(, nrow = ncol(matrix_data_list[[1]])-2, ncol = 0)
for (line in 1:nrow(compare)) {
#	print(line)
#	print(compare)

	group_name1 <- as.vector(compare[line,1])
	group_name2 <- as.vector(compare[line,2])

	group1 <- t(matrix_data_list[[group_name1]])[-(1:2),]
	group1_use <- apply(group1,c(1,2), function(x){as.numeric(x)})
	sample_names1 <- as.vector(t(matrix_data_list[[group_name1]])[1,])
	sample_names1 <- paste("log2",sample_names1,sep=" ")

	group2 <- t(matrix_data_list[[group_name2]])[-(1:2),]
	group2_use <- apply(group2,c(1,2), function(x){as.numeric(x)})
	sample_names2 <- as.vector(t(matrix_data_list[[group_name2]])[1,])
	sample_names2 <- paste("log2",sample_names2,sep=" ")

	plus_three_vector <- vector()
	plus_two_vector <- vector()
	t_diff_vector <- vector()
	p_val_vector <- vector()
	adjust_pval_vector <- vector()

	for (i in 1:nrow(group1)){

		group1_vector <- group1_use[i,]
		group2_vector <- group2_use[i,]
		if ( (sum(length(which(is.na(group1_vector)))) < length(group1_vector)-1 && length(group1_vector) >= 2) && (sum(length(which(is.na(group2_vector)))) < length(group2_vector)-1 && length(group2_vector) >= 2)) {

			output <- t.test(group1_vector,group2_vector)
			p_val <- as.vector(unlist(output['p.value']))
			p_val_vector[i] <- p_val

			mean1 <- unlist(output["estimate"])[[1]]
			mean2 <- unlist(output["estimate"])[[2]]
			t_diff <- round((mean1 - mean2),digits = 5)

#			t <- as.vector(unlist(output["statistic"]))
			t_diff_vector[i] <- t_diff
			plus_three_vector[i] <- ""

		} else {
			p_val_vector[i] <- NA
			t_diff_vector[i] <- NA

			NAs_group1 <- sum(length(which(is.na(group1_vector))))
			NAs_group2 <- sum(length(which(is.na(group2_vector))))
			if ( ( NAs_group1 == length(group1_vector)  && (NAs_group2 == 0 || (NAs_group2 == 1 && length(group2_vector) > 2)) )  || ( (NAs_group1 == 0 || (NAs_group1 == 1 && length(group1_vector) > 2) ) && NAs_group2 == length(group2_vector)) ) {
				plus_three_vector[i] <- "++"
			} else {
				plus_three_vector[i] <- ""
			}

		}

	} 

	adjust_pval_vector <- p.adjust(p_val_vector,"BH")
	
	sign_p_vector <- vector()
	sign_diff_vector <- vector()
	sign_vector <- vector()

	for (index in 1:length(t_diff_vector)) {

#		p_val <- adjust_pval_vector[index]
		p_val <- p_val_vector[index]
		diff <- t_diff_vector[index]

		if (! is.na(p_val)) {

			if ( (p_val <= FDR_threshold) && (abs(diff) >= diff_threshold)) {

				sign_p_vector <- c(sign_p_vector,-log10(p_val))
				sign_diff_vector <- c(sign_diff_vector,diff)
				sign_vector[index] <- "+"

			} else {

				sign_vector[index] <- ""
			}
		} else {

			sign_vector[index] <- ""
		}
	}


	result <- cbind(group1_use,group2_use,round(p_val_vector,digits= 6),round(adjust_pval_vector,digits= 6),t_diff_vector,sign_vector,plus_three_vector)

	plot(t_diff_vector,-log10(p_val_vector),main=paste(group_name1," vs ",group_name2,sep=""),xlab=paste("log2 Fold change\n > 1 enriched in ",group_name1,"; < -1 enriched in",group_name2),ylab="-log10 p-value",cex=0.5)
	points(sign_diff_vector,sign_p_vector,col="red",cex=0.5)

	colnames(result) <- c(sample_names1,sample_names2,paste("p-val ",group_name1,"_vs_",group_name2,sep=""),"Adjust p-val (BH)",paste("Difference_",group_name1,"_vs_",group_name2,sep=""),paste("Sign p-val:",FDR_threshold,"Diff:",diff_threshold,sep=" "),"Sign")

	result_matrix <- cbind(result_matrix,result)
}
dev.off()





annotated_result_matrix <- cbind(annotation,result_matrix)
#if (name_test != 0) {

#	if (ibaq_test != 0) {

#		colnames(annotated_result_matrix) <- c("Majority protein IDs","Protein names","Gene names","Score",colnames(iBAC_matrix_use_log10),colnames(result_matrix))
#	} else {
#		colnames(annotated_result_matrix) <- c("Majority protein IDs","Protein names","Gene names","Score",colnames(result_matrix))
#	}
#} else {
#	if (ibaq_test != 0) {
#		colnames(annotated_result_matrix) <- c("Majority protein IDs","Score",colnames(iBAC_matrix_use_log10),colnames(result_matrix))
#	} else {
#		colnames(annotated_result_matrix) <- c("Majority protein IDs","Score",colnames(result_matrix))
#	}
#}




xlsx.addTitle<-function(sheet, rowIndex, title, titleStyle){
  rows <-createRow(sheet,rowIndex=rowIndex)
  sheetTitle <-createCell(rows, colIndex=1)
  setCellValue(sheetTitle[[1,1]], title)
  setCellStyle(sheetTitle[[1,1]], titleStyle)
}




excel_file_name <- "Comparison_Results_n.xlsx"
wb<-createWorkbook(type="xlsx")

sheet <- createSheet(wb, sheetName = "Results")
workbook_start_position <- 1
TABLE_ROWNAMES_STYLE <- CellStyle(wb) + Font(wb, heightInPoints=10, isBold=FALSE)
TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, heightInPoints=10, isBold=TRUE, color ="9", name="Arial") + Border(color="black", position=c("TOP", "BOTTOM"), pen=c("BORDER_THIN", "BORDER_THIN")) + Fill(foregroundColor="#0069AA") + Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER",vertical="VERTICAL_CENTER")
TITLE_STYLE <- CellStyle(wb)+ Font(wb,  heightInPoints=8,isBold=FALSE)

annotated_result_matrix <- annotated_result_matrix[order(annotated_result_matrix[,(ncol(annotated_result_matrix))],annotated_result_matrix[,(ncol(annotated_result_matrix)-1)],annotated_result_matrix[,4],decreasing=TRUE),]
addDataFrame(as.data.frame(annotated_result_matrix),sheet, startRow=workbook_start_position, startColumn=1,colnamesStyle = TABLE_COLNAMES_STYLE,rownamesStyle = TABLE_ROWNAMES_STYLE,row.names = FALSE)
setColumnWidth(sheet, colIndex=c(1:ncol(annotated_result_matrix)), colWidth=16)

saveWorkbook(wb,excel_file_name)




#################results for imputed data

pdf("Vulcano_plot_imputed_n.pdf")

result_matrix <- matrix(, nrow = ncol(imputed_matrix_data_list[[1]]), ncol = 0)
for (line in 1:nrow(compare)) {
#	print(line)
	group_name1 <- as.vector(compare[line,1])
	group_name2 <- as.vector(compare[line,2])

	group1 <- t(imputed_matrix_data_list[[group_name1]])
	group1_use <- apply(group1,c(1,2), function(x){as.numeric(x)})
	sample_names1 <- paste("log2",colnames(group1_use),sep=" ")

	group2 <- t(imputed_matrix_data_list[[group_name2]])
	group2_use <- apply(group2,c(1,2), function(x){as.numeric(x)})
	sample_names2 <- paste("log2",colnames(group2_use),sep=" ")

	plus_three_vector <- vector()
	plus_two_vector <- vector()
	t_diff_vector <- vector()
	p_val_vector <- vector()
	adjust_pval_vector <- vector()

	for (i in 1:nrow(group1)){

		group1_vector <- group1_use[i,]
		group2_vector <- group2_use[i,]
		if ( (sum(length(which(is.na(group1_vector)))) < length(group1_vector)-1 && length(group1_vector) >= 2) && (sum(length(which(is.na(group2_vector)))) < length(group2_vector)-1 && length(group2_vector) >= 2)) {

			output <- t.test(group1_vector,group2_vector)
			p_val <- as.vector(unlist(output['p.value']))
			p_val_vector[i] <- p_val

			mean1 <- unlist(output["estimate"])[[1]]
			mean2 <- unlist(output["estimate"])[[2]]
			t_diff <- round((mean1 - mean2),digits = 5)

#			t <- as.vector(unlist(output["statistic"]))
			t_diff_vector[i] <- t_diff
			plus_three_vector[i] <- ""

		} else {
			p_val_vector[i] <- NA
			t_diff_vector[i] <- NA

			NAs_group1 <- sum(length(which(is.na(group1_vector))))
			NAs_group2 <- sum(length(which(is.na(group2_vector))))
			if ( ( NAs_group1 == length(group1_vector)  && (NAs_group2 == 0 || (NAs_group2 == 1 && length(group2_vector) > 2)) )  || ( (NAs_group1 == 0 || (NAs_group1 == 1 && length(group1_vector) > 2) ) && NAs_group2 == length(group2_vector)) ) {
				plus_three_vector[i] <- "++"
			} else {
				plus_three_vector[i] <- ""
			}

		}

	} 

	adjust_pval_vector <- p.adjust(p_val_vector,"BH")

	
	sign_p_vector <- vector()
	sign_diff_vector <- vector()
	sign_vector <- vector()

	for (index in 1:length(t_diff_vector)) {

		p_val <- p_val_vector[index]
		diff <- t_diff_vector[index]

		if (! is.na(p_val)) {

			if ( (p_val <= FDR_threshold) && (abs(diff) >= diff_threshold)) {

				sign_p_vector <- c(sign_p_vector,-log10(p_val))
				sign_diff_vector <- c(sign_diff_vector,diff)
				sign_vector[index] <- "+"

			} else {

				sign_vector[index] <- ""
			}
		}
	}


	result <- cbind(group1_use,group2_use,round(p_val_vector,digits= 6),round(adjust_pval_vector,digits= 6),t_diff_vector,sign_vector,plus_three_vector)

	plot(t_diff_vector,-log10(p_val_vector),main=paste(group_name1," vs ",group_name2,sep=""),xlab=paste("log2 Fold change\n > 1 enriched in ",group_name1,"; < -1 enriched in",group_name2),ylab="-log10 p-value",cex=0.5)
	points(sign_diff_vector,sign_p_vector,col="red",cex=0.5)

	colnames(result) <- c(sample_names1,sample_names2,paste("p-val ",group_name1,"_vs_",group_name2,sep=""),"Adjust p-val (BH)",paste("Difference_",group_name1,"_vs_",group_name2,sep=""),paste("Sign p-val:",FDR_threshold,"Diff:",diff_threshold,sep=" "),"Sign")

	result_matrix <- cbind(result_matrix,result)
}
dev.off()


annotated_result_matrix <- cbind(annotation,result_matrix)
annotated_result_matrix <- annotated_result_matrix[order(annotated_result_matrix[,(ncol(annotated_result_matrix))],annotated_result_matrix[,(ncol(annotated_result_matrix)-1)],annotated_result_matrix[,4],decreasing=TRUE),]
write.table(annotated_result_matrix,"Comparison_Results_imputed_n.txt",quote=FALSE,row.names=FALSE,sep="\t")

#col_names <- colnames(annotated_result_matrix)
#col_names <- str_replace_all(col_names,"d\\$","")
#col_names <- str_replace_all(col_names,"\\."," ")
#colnames(annotated_result_matrix) <- col_names


excel_file_name <- "Comparison_Results_imputed_n.xlsx"
wb<-createWorkbook(type="xlsx")

sheet <- createSheet(wb, sheetName = "Results")
workbook_start_position <- 1
TABLE_ROWNAMES_STYLE <- CellStyle(wb) + Font(wb, heightInPoints=10, isBold=FALSE)
TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, heightInPoints=10, isBold=TRUE, color ="9", name="Arial") + Border(color="black", position=c("TOP", "BOTTOM"), pen=c("BORDER_THIN", "BORDER_THIN")) + Fill(foregroundColor="#0069AA") + Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER",vertical="VERTICAL_CENTER")
TITLE_STYLE <- CellStyle(wb)+ Font(wb,  heightInPoints=8,isBold=FALSE)


addDataFrame(as.data.frame(annotated_result_matrix),sheet, startRow=workbook_start_position, startColumn=1,colnamesStyle = TABLE_COLNAMES_STYLE,rownamesStyle = TABLE_ROWNAMES_STYLE,row.names = FALSE)
setColumnWidth(sheet, colIndex=c(1:ncol(annotated_result_matrix)), colWidth=16)

saveWorkbook(wb,excel_file_name)





#################results for imputed data

pdf("Vulcano_plot_imputed.pdf")

result_matrix <- matrix(, nrow = ncol(imputed_matrix_data_list[[1]]), ncol = 0)
for (line in 1:nrow(compare)) {
#	print(line)
	group_name1 <- as.vector(compare[line,][,1])
	group_name2 <- as.vector(compare[line,][,2])

	group1 <- t(imputed_matrix_data_list[[group_name1]])
	group1_use <- apply(group1,c(1,2), function(x){as.numeric(x)})
	sample_names1 <- paste("log2 LFQ",colnames(group1_use),sep=" ")

	group2 <- t(imputed_matrix_data_list[[group_name2]])
	group2_use <- apply(group2,c(1,2), function(x){as.numeric(x)})
	sample_names2 <- paste("log2 LFQ",colnames(group2_use),sep=" ")

	NA_count_vector_group1 <- NA_count_vector_list[[group_name1]]
	NA_count_vector_group2 <- NA_count_vector_list[[group_name2]]


	plus_three_vector <- vector()
	plus_two_vector <- vector()
	t_diff_vector <- vector()
	p_val_vector <- vector()
	adjust_pval_vector <- vector()

	for (i in 1:nrow(group1)){

		NAs_group1_nonimputed <- NA_count_vector_group1[i]
		NAs_group2_nonimputed <- NA_count_vector_group2[i]

		group1_vector <- group1_use[i,]
		group2_vector <- group2_use[i,]

		if ( (NAs_group1_nonimputed < ceiling(length(group1_vector)/2) && length(group1_vector) >= 2 ) || ( NAs_group2_nonimputed < ceiling(length(group2_vector)/2) && length(group2_vector) >= 2 ) ) {

#		if ( (sum(length(which(is.na(group1_vector)))) < length(group1_vector)-1 && length(group1_vector) >= 2) && (sum(length(which(is.na(group2_vector)))) < length(group2_vector)-1 && length(group2_vector) >= 2)) {

			output <- t.test(group1_vector,group2_vector)
			p_val <- as.vector(unlist(output['p.value']))
			p_val_vector[i] <- p_val

			mean1 <- unlist(output["estimate"])[[1]]
			mean2 <- unlist(output["estimate"])[[2]]
			t_diff <- round((mean1 - mean2),digits = 5)

#			t <- as.vector(unlist(output["statistic"]))
			t_diff_vector[i] <- t_diff
#			plus_three_vector[i] <- ""

		} else {
			p_val_vector[i] <- NA
			t_diff_vector[i] <- NA
#			plus_three_vector[i] <- ""

#			NAs_group1 <- sum(length(which(is.na(group1_vector))))
#			NAs_group2 <- sum(length(which(is.na(group2_vector))))

#			if ( ( NAs_group1 == length(group1_vector)  && (NAs_group2 == 0 || (NAs_group2 == 1 && length(group2_vector) > 2)) )  || ( (NAs_group1 == 0 || (NAs_group1 == 1 && length(group1_vector) > 2) ) && NAs_group2 == length(group2_vector)) ) {

#				plus_three_vector[i] <- "++"
#			} else {
#				plus_three_vector[i] <- ""
#			}

		}

	} 

	adjust_pval_vector <- p.adjust(p_val_vector,"BH")

	
	sign_p_vector <- vector()
	sign_diff_vector <- vector()
	sign_vector <- vector()

	for (index in 1:length(t_diff_vector)) {

		p_val <- p_val_vector[index]
		diff <- t_diff_vector[index]

		if (! is.na(p_val)) {

			if ( (p_val <= FDR_threshold) && (abs(diff) >= diff_threshold)) {

				sign_p_vector <- c(sign_p_vector,-log10(p_val))
				sign_diff_vector <- c(sign_diff_vector,diff)
				sign_vector[index] <- "+"

			} else {

				sign_vector[index] <- ""
			}
		} else {
			sign_vector[index] <- ""
		}
	}


	result <- cbind(group1_use,group2_use,round(p_val_vector,digits= 6),round(adjust_pval_vector,digits= 6),t_diff_vector,sign_vector)

	plot(t_diff_vector,-log10(p_val_vector),main=paste(group_name1," vs ",group_name2,sep=""),xlab=paste("log2 Fold change\n > 1 enriched in ",group_name1,"; < -1 enriched in",group_name2),ylab="-log10 p-value",cex=0.5)
	points(sign_diff_vector,sign_p_vector,col="red",cex=0.5)

	colnames(result) <- c(sample_names1,sample_names2,paste("p-val ",group_name1,"_vs_",group_name2,sep=""),"Adjust p-val (BH)",paste("Difference_",group_name1,"_vs_",group_name2,sep=""),paste("Sign p-val:",FDR_threshold,"Diff:",diff_threshold,sep=" "))

	result_matrix <- cbind(result_matrix,result)
}
dev.off()


annotated_result_matrix <- cbind(annotation_matrix[,1:3],result_matrix)

#col_names <- colnames(annotated_result_matrix)
#col_names <- str_replace_all(col_names,"d\\$","")
#col_names <- str_replace_all(col_names,"\\."," ")
#colnames(annotated_result_matrix) <- col_names


excel_file_name <- "Comparison_Results_imputed.xlsx"
wb<-createWorkbook(type="xlsx")

sheet <- createSheet(wb, sheetName = "Results")
workbook_start_position <- 1
TABLE_ROWNAMES_STYLE <- CellStyle(wb) + Font(wb, heightInPoints=10, isBold=FALSE)
TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, heightInPoints=10, isBold=TRUE, color ="9", name="Arial") + Border(color="black", position=c("TOP", "BOTTOM"), pen=c("BORDER_THIN", "BORDER_THIN")) + Fill(foregroundColor="#0069AA") + Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER",vertical="VERTICAL_CENTER")
TITLE_STYLE <- CellStyle(wb)+ Font(wb,  heightInPoints=8,isBold=FALSE)

annotated_result_matrix <- annotated_result_matrix[order(annotated_result_matrix[,(ncol(annotated_result_matrix))],annotated_result_matrix[,(ncol(annotated_result_matrix)-1)],annotated_result_matrix[,4],decreasing=TRUE),]
addDataFrame(as.data.frame(annotated_result_matrix),sheet, startRow=workbook_start_position, startColumn=1,colnamesStyle = TABLE_COLNAMES_STYLE,rownamesStyle = TABLE_ROWNAMES_STYLE,row.names = FALSE)
setColumnWidth(sheet, colIndex=c(1:ncol(annotated_result_matrix)), colWidth=16)

sheet <- createSheet(wb, sheetName = "Annotation")
workbook_start_position <- 1
TABLE_ROWNAMES_STYLE <- CellStyle(wb) + Font(wb, heightInPoints=10, isBold=FALSE)
TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, heightInPoints=10, isBold=TRUE, color ="9", name="Arial") + Border(color="black", position=c("TOP", "BOTTOM"), pen=c("BORDER_THIN", "BORDER_THIN")) + Fill(foregroundColor="#0069AA") + Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER",vertical="VERTICAL_CENTER")
TITLE_STYLE <- CellStyle(wb)+ Font(wb,  heightInPoints=8,isBold=FALSE)
addDataFrame(as.data.frame(annotation_matrix),sheet, startRow=workbook_start_position, startColumn=1,colnamesStyle = TABLE_COLNAMES_STYLE,rownamesStyle = TABLE_ROWNAMES_STYLE,row.names = FALSE)
setColumnWidth(sheet, colIndex=c(1:ncol(annotation_matrix)), colWidth=16)

sheet <- createSheet(wb, sheetName = "Parameters")
workbook_start_position <- 1
TABLE_ROWNAMES_STYLE <- CellStyle(wb) + Font(wb, heightInPoints=10, isBold=FALSE)
TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, heightInPoints=10, isBold=TRUE, color ="9", name="Arial") + Border(color="black", position=c("TOP", "BOTTOM"), pen=c("BORDER_THIN", "BORDER_THIN")) + Fill(foregroundColor="#0069AA") + Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER",vertical="VERTICAL_CENTER")
TITLE_STYLE <- CellStyle(wb)+ Font(wb,  heightInPoints=8,isBold=FALSE)


addDataFrame(as.data.frame(parameters),sheet, startRow=workbook_start_position, startColumn=1,colnamesStyle = TABLE_COLNAMES_STYLE,rownamesStyle = TABLE_ROWNAMES_STYLE,row.names = FALSE)
setColumnWidth(sheet, colIndex=c(1:ncol(parameters)), colWidth=16)

summary_print <- cbind(as.vector(summary$Raw.file),as.vector(summary$Experiment),as.vector(summary$MS),as.vector(summary$MS.MS),as.vector(summary$MS.MS.Identified),as.vector(summary$Peptide.Sequences.Identified))

#colnames(summllary_print) <- c("Raw files","Experiment","MS","MS/MS","Identified MS/MS","Identified Peptide Sequences")

print_rows_length <- length(summary$Experiment[summary$Experiment != ""])
summary_print <- summary_print[1:print_rows_length,]

sheet <- createSheet(wb, sheetName = "Summary")
workbook_start_position <- 1
TABLE_ROWNAMES_STYLE <- CellStyle(wb) + Font(wb, heightInPoints=10, isBold=FALSE)
TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, heightInPoints=10, isBold=TRUE, color ="9", name="Arial") + Border(color="black", position=c("TOP", "BOTTOM"), pen=c("BORDER_THIN", "BORDER_THIN")) + Fill(foregroundColor="#0069AA") + Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER",vertical="VERTICAL_CENTER")
TITLE_STYLE <- CellStyle(wb)+ Font(wb,  heightInPoints=8,isBold=FALSE)

addDataFrame(as.data.frame(summary_print),sheet, startRow=workbook_start_position, startColumn=1,colnamesStyle = TABLE_COLNAMES_STYLE,rownamesStyle = TABLE_ROWNAMES_STYLE,row.names = FALSE)
setColumnWidth(sheet, colIndex=c(1:ncol(summary_print)), colWidth=16)

saveWorkbook(wb,excel_file_name)








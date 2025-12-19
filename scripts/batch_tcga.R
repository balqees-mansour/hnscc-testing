library("sva") #Note this exercise requires sva (>= v3.36.0) which is only available for R (>= v4.x)
library("ggplot2")
library("gridExtra")
library("edgeR")
library("UpSetR")


# age , gender, tissue of origin, race 
# Add the condition column based on the prefix of the id
sample_info_all <- colData(hnsc_no_metastatic)

uncorrected_data <- counts_data
sample_names = colnames(uncorrected_data)

#review data structure
head(uncorrected_data)
dim(uncorrected_data)



#define conditions, library methods, and replicates
conditions = sample_info_all$tissue_type
gender = as.character(sample_info_all$gender)
origin = as.character(sample_info_all$tissue_or_organ_of_origin)

replicates = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4)

#calculate principal components for the uncorrected data
pca_uncorrected_obj = prcomp(uncorrected_data[,sample_names])

#pull PCA values out of the PCA object
pca_uncorrected = as.data.frame(pca_uncorrected_obj[2]$rotation)

#assign labels to the data frame
pca_uncorrected[,"condition"] = conditions


#plot the PCA
#create a classic 2-dimension PCA plot (first two principal components) with conditions and library methods indicated
cols <- c("Tumor" = "#481567FF", "Normal" = "#1F968BFF")
p1 = ggplot(data = pca_uncorrected, aes(x = PC1, y = PC2, color = conditions, shape = gender))
p1 = p1 + geom_point(size = 3)
p1 = p1 + stat_ellipse(type = "norm", linetype = 2)
p1 = p1 + labs(title = "PCA, RNA-seq counts for TCGA samples", color = "conditions", shape="gender")
p1 = p1 + scale_colour_manual(values = cols)


p2 = ggplot(data = pca_uncorrected, aes(x = PC1, y = PC2, color = conditions, shape = origin))
p2 = p2 + geom_point(size = 3)
p2 = p2 + stat_ellipse(type = "norm", linetype = 2)
p2 = p2 + labs(title = "PCA, RNA-seq counts for TCGA samples", color = "conditions", shape="origin")
p2 = p2 + scale_colour_manual(values = cols)



























#now run ComBat_seq gender
corrected_data = ComBat_seq(counts = as.matrix(uncorrected_data[,sample_names]), batch = platform, group = conditions)

corrected <- ComBat_seq(counts=as.matrix(uncorrected_data[, sample_names]), batch=platform, group=NULL)

#compare dimensions of corrected and uncorrected data sets
dim(uncorrected_data)
dim(corrected)


corrected == uncorrected_data
corrected[1:10,1:3]
uncorrected_data[1:10,1:3]

#calculate principal components for the uncorrected data
pca_corrected_obj = prcomp(corrected[, sample_names])

#pull PCA values out of the PCA object
pca_corrected = as.data.frame(pca_corrected_obj[2]$rotation)

#assign labels to the data frame
pca_corrected[,"condition"] = conditions
pca_corrected[,"ages"] = ages
pca_corrected[,"gender"] = gender
pca_corrected[,"ethinicity"] = ethinicity

#as above, create a PCA plot for comparison to the uncorrected data
cols <- c("HNSCC" = "#481567FF", "normal" = "#1F968BFF")
p2 = ggplot(data = pca_corrected, aes(x = PC1, y = PC2, color = condition, shape = conditions))
p2 = p2 + geom_point(size = 3)
p2 = p2 + stat_ellipse(type = "norm", linetype = 2)
p2 = p2 + labs(title = "PCA, RNA-seq counts for  (batch corrected data)", color = "Condition", shape = "conditions")
p2 = p2 + scale_colour_manual(values = cols)









#now run ComBat_seq ethinicity
corrected_data = ComBat_seq(counts = as.matrix(uncorrected_data[,sample_names]), batch = ethinicity, group = conditions)

#join the gene and chromosome names onto the now corrected counts from ComBat_seq
corrected_data = cbind(uncorrected_data[, c("Gene", "Chr")], corrected_data)

#compare dimensions of corrected and uncorrected data sets
dim(uncorrected_data)
dim(corrected_data)










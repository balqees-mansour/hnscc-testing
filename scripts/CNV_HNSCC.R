


# Load required library
library(TCGAbiolinks)


# Step 2: Query segment-level CNV data (available in GDC)
query.cnv <- GDCquery(
  project = "TCGA-HNSC",
  data.category = "Copy Number Variation",
  data.type = "Masked Copy Number Segment",  # Primary available CNV type [web:20][web:12]
  workflow.type = "DNAcopy",  # Or "ASCAT", check results from Step 1
  sample.type = c("Primary Tumor"),
  access = "open"
)

print(query.cnv)

# Step 3: Download (smaller files than gene-level)
GDCdownload(query.cnv, method = "api", directory = "GDCdata")

# Step 4: Prepare data
cnv.seg <- GDCprepare(query.cnv, directory = "GDCdata")

# After running GDCprepare
# The returned 'cnv.seg' for segment data is a data frame or GRanges object with rows as segments

# Check the class and structure
class(cnv.seg)
head(cnv.seg)

# Convert to data.frame if it is GRanges
if ("GRanges" %in% class(cnv.seg)) {
  cnv.df <- as.data.frame(cnv.seg)
} else {
  cnv.df <- cnv.seg
}

# 'cnv.df' should have columns like:
# Chromosome, Start, End, Sample, Segment_Mean, etc.

# Filter amplifications
amps <- cnv.df[cnv.df$Segment_Mean > 0.2, ]
# Filter deletions
dels <- cnv.df[cnv.df$Segment_Mean < -0.2, ]

# Example overview
table(cnv.df$Chromosome)
summary(cnv.df$Segment_Mean)






library(dplyr)
library("DESeq2")

# Open files 
GPL15615 <- read.delim("GPL15615.txt")
GSE38246 <- read.delim("GSE38246.txt")
metadata <- read.delim("metadata.txt")

# Merge by ID to obtain the Gene names in the expression data file, and clean it with the grep function. Convert the ORF columm into the Raw names

file1 <- merge(GPL15615,GSE38246,by="ID")
file1<-file1[,-grep("ID|ORGANISM|SEQUENCE|Reporter.Group|Control.Type|Composite.Element.Comment|Composite.Element.Comment.1|SPOT_ID",colnames(file1))]
rownames(file1) = make.names(file1$ORF, unique=TRUE)
file1<-file1[,-grep("ORF",colnames(file1))]

dim(file1)
[1] 45696   122

#Convert null to 0
file1[file1=='NULL'] <- 0

#Clean the metadata file
metadata<-metadata[,-grep("X.Sample_title|X.Sample_geo_accesion|X.Sample_characteristics_ch2|X.Sample_geo_accesion",colnames(metadata))]
rownames(metadata) = metadata$ID_REF

#Normalization Start, the data would be normalizaed according to the condition, in this case the disease stage
count1 <- DESeqDataSetFromMatrix(countData=file1, colData=metadata, design= ~ condition) 



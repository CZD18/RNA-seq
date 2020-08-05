library("dplyr")
library("DESeq2")
library("ggplot2")

# Open files 
GSE87505 <- read.delim("GSE87505_All_Count.tsv")
metadata <- read.delim("metadata.txt")

# Merge by ID to obtain the Gene names in the expression data file, and clean it with the grep function. Convert the X columm into the Raw names (Special consideration
# to avoid same rownames

file1 <- GSE87505[apply(GSE87505[,-1], 1, function(x) !all(x==0)),]

 dim(GSE87505)
#[1] 21742   109
> dim(file1)
#[1] 21733   109

rownames(file1) = make.names(file1$X, unique=TRUE)
file1<-file1[,-grep("X",colnames(file1))]

dim(file1)
#[1] 21733   108

#Convert null to 0
file1[file1=='NULL'] <- 0

#Clean the metadata file
rownames(metadata) =  make.names(metadata$ID, unique=TRUE)

#Normalization Start, the data would be normalizaed according to the condition, in this case the disease stage
count1 <- DESeqDataSetFromMatrix(countData=file1, colData=metadata, design= ~ condition) 
count.object <- DESeq(count1)
#Then we normalize using VST method
Vn <- vst(count.object)
data = assay(Vn)

#Save the normalized data as a .TXT file
write.table(data, sep="\t",file="Normalized.txt", row.names=TRUE,col.names=NA,quote=FALSE)
                        
#Generating a PCA plot from the normalized data
plotPCA(Vn, intgroup=c("condition")) 

 #Generating the VST transformation plot
df <- as_data_frame(assay(Vn)[, 1:2]) %>% mutate(transformation = "vst")
colnames(df)[1:2] <- c("x", "y")
lvls <- c("vst")
df$transformation <- factor(df$transformation, levels=lvls)
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) + coord_fixed() + facet_grid( . ~ transformation) 
                        
#The plot result is under the name "VST.pdf"

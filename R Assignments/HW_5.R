###### Task 1
### a) Load the data in “Gene_Expression_Data.xlsx”, “Gene_Information.csv”, and “Sample_Information.tsv” into R.
gene_exp <- read.csv("C:/Users/icwin/OneDrive - Indiana University/B573 Prog SciInformatics/R/Gene_Expres_Data.csv")
gene_info <- read.csv("C:/Users/icwin/OneDrive - Indiana University/B573 Prog SciInformatics/Python/Gene_Information.csv")
samp_info <- read.delim("C:/Users/icwin/OneDrive - Indiana University/B573 Prog SciInformatics/Python/Sample_Information.tsv", sep = "\t")
head(gene_exp)
head(gene_info)
head(samp_info)

### b) Change the sample names from the “Gene_Expression_Data.xlsx”, based upon the phenotype presented in “Sample_Information.tsv”
# create a mapping between sample names in the gene_exp & sample names in the samp_info
sample_name_mapping <- c(
  'GSM820516' = 'tumor_1',
  'GSM820517' = 'normal_1',
  'GSM820518' = 'tumor_2',
  'GSM820519' = 'normal_2',
  'GSM820520' = 'tumor_3',
  'GSM820521' = 'normal_3',
  'GSM820522' = 'tumor_4',
  'GSM820523' = 'normal_4',
  'GSM820524' = 'tumor_5',
  'GSM820525' = 'normal_5',
  'GSM820526' = 'tumor_6',
  'GSM820527' = 'normal_6',
  'GSM820528' = 'tumor_7',
  'GSM820529' = 'normal_7',
  'GSM820530' = 'tumor_8',
  'GSM820531' = 'normal_8',
  'GSM820532' = 'tumor_9',
  'GSM820533' = 'normal_9'
)
# replace column names in the data frame using the mapping
colnames(gene_exp) <- sapply(colnames(gene_exp), function(x) {
  if (x %in% names(sample_name_mapping)) {
    sample_name_mapping[x]
  } else {
    x
  }
})
print(colnames(gene_exp))

### c) Split the merged data from part b, into to 2 parts, based upon their labeled phenotype (ie. tumor or normal)
# get first column name (probe ID)
probe_id <- colnames(gene_exp)[1]
# filter columns with 'tumor' or  probe ID
gene_exp_tumor <- gene_exp[, grepl("tumor", colnames(gene_exp)) | colnames(gene_exp) == probe_id]
# filter columns with 'normal' name or  probe ID
gene_exp_normal <- gene_exp[, grepl("normal", colnames(gene_exp)) | colnames(gene_exp) == probe_id]
# check data frames
print(head(gene_exp_tumor))
print(head(gene_exp_normal))

### d) find the avg expression for all genes from the 2 data sets from part d
# calculate mean of columns
gene_exp_normal$probe_id_avg_norm <- rowMeans(gene_exp_normal[, -1])
# calculate mean of columns
gene_exp_tumor$probe_id_avg_tumor <- rowMeans(gene_exp_tumor[, -1])
# check data frames
head(gene_exp_normal)
head(gene_exp_tumor)

### e) Determine the log2 fold change for each Probe between the two groups log2((Tumour – Control) / Control)
# calculate log2 fold change for normal dataset
gene_exp_normal$log2_fold_change <- log2((gene_exp_tumor$probe_id_avg_tumor - gene_exp_normal$probe_id_avg_norm) / gene_exp_normal$probe_id_avg_norm)
# calculate log2 fold change for tumor dataset
gene_exp_tumor$log2_fold_change <- log2((gene_exp_tumor$probe_id_avg_tumor - gene_exp_normal$probe_id_avg_norm) / gene_exp_normal$probe_id_avg_norm)
# Check data frames
head(gene_exp_normal)
head(gene_exp_tumor)


### f) Use the data from part e and “Gene_Information.csv” to identify all genes fold change magnitude (absolute value) was greater than 5
# add Gene_Name column to gene_exp_normal and gene_exp_tumor
gene_exp_normal$Gene_Name <- gene_info$Symbol
gene_exp_tumor$Gene_Name <- gene_info$Symbol
# calculate abs value of fold change
gene_exp_normal$abs_fold_change <- abs(gene_exp_normal$fold_change)
gene_exp_tumor$abs_fold_change <- abs(gene_exp_tumor$fold_change)
# mark genes as significant if the abs fold change > 5
gene_exp_normal$significant <- ifelse(gene_exp_normal$abs_fold_change > 5, "Yes", "No")
gene_exp_tumor$significant <- ifelse(gene_exp_tumor$abs_fold_change > 5, "Yes", "No")
# check data frames
head(gene_exp_normal)
head(gene_exp_tumor)

### g) Add a column to the result of part f to include if the gene was higher expressed in “Normal” or “Tumor” samples
# Check row counts to ensure alignment -- i was getting an error here
if (nrow(gene_exp) == nrow(gene_exp_tumor)) {
  # Assign the fold_change column from gene_exp_tumor to gene_exp
  gene_exp$fold_change <- gene_exp_tumor$fold_change
} else {
  stop("Row counts of gene_exp and gene_exp_tumor do not match!")
}
# validate fold_change is numeric and calculate abs_fold_change -- further validation for error
gene_exp$fold_change <- as.numeric(gene_exp$fold_change)  # Ensure numeric type
gene_exp$abs_fold_change <- abs(gene_exp$fold_change)  # Calculate absolute fold change
# create subset of genes with abs_fold_change > 5
sub_gene_exp <- gene_exp[gene_exp$abs_fold_change > 5, ]
# merge other columns from gene_exp_tumor into sub_gene_exp
if ("probe_id" %in% colnames(gene_exp) && "probe_id" %in% colnames(gene_exp_tumor)) {
  sub_gene_exp <- merge(
    sub_gene_exp,
    gene_info[, c("Probe_ID", "Symbol", "Chromosome")],
    y.x = colnames(gene_exp)[1],  # Column to merge on from sub_gene_exp (usually 'Probe_ID' or similar)
    by.y = "Probe_ID",  # Column to merge on from gene_info
    all.x = TRUE
  )
} else {
  stop("probe_id column is missing in one or both data frames!")
}
# add column for higher expression based on fold_change
sub_gene_exp$higher_expressed <- ifelse(sub_gene_exp$fold_change > 0, "tumor", "normal")
# check data frame
str(sub_gene_exp)
head(sub_gene_exp)

###### Task 2
### a) Perform exploratory data analysis on the genes from part 1g
str(gene_exp_normal)
summary(gene_exp_normal)
str(gene_exp_tumor)
summary(gene_exp_tumor)
# distribution of expression levels for normal genes
hist(
  gene_exp_normal$probe_id_avg_norm,
  breaks = 30,
  main = "Distribution of Gene Expression in Normal Samples",
  xlab = "Expression Level",
  col = "blue"
)
# distribution of expression levels for tumor genes
hist(
  gene_exp_tumor$probe_id_avg_tumor,
  breaks = 30,
  main = "Distribution of Gene Expression in Tumor Samples",
  xlab = "Expression Level",
  col = "red"
)

# b. Create a histogram showing the distribution of the number of differentially expressed genes (DEGs) by chromosome
# get Probe_ID from gene_exp
probe_id <- gene_exp$Probe_ID
# filter gene_info to include where Probe_ID is in list of probe_id
genes_inf_2 <- gene_info[gene_info$Probe_ID %in% probe_id, ]
# count number of each chromosome
chromosome_counts <- table(genes_inf_2$Chromosome)
# conver counts to data frame 
chromosome_counts_df <- as.data.frame(chromosome_counts)
colnames(chromosome_counts_df) <- c("Chromosome", "Count")
# plot histogram
library(ggplot2)
ggplot(chromosome_counts_df, aes(x = Chromosome, y = Count)) +
  geom_histogram(stat = "identity", binwidth = 1, fill = "purple", color = "blue") +
  theme_minimal() +
  labs(
    title = "Chromosome Distribution of Filtered Genes",
    x = "Chromosome",
    y = "Number of Genes"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability

### c) Make another histogram showing the distribution of DEGs by chromosome segregated by sample type (Normal or Tumor)
# subset tumor gene expression
tumor_gene_exp <- gene_exp[gene_exp$higher_expressed == "tumor", ]
# filter gene info for tumor genes
get <- gene_info$Probe_ID %in% probe_id
genes_inf_2 <- gene_info[get, ]
# plot chromosome distribution for tumor genes
barplot(table(genes_inf_2$Chromosome), main = "Tumor Chromosome Distribution", xlab = "Chromosome", ylab = "Frequency")
# subset normal gene expression
normal_gene_exp <- gene_exp[gene_exp$higher_expressed == "normal", ]
# filter gene info for normal genes
get <- gene_info$Probe_ID %in% probe_id
genes_inf_2 <- gene_info[get, ]
# plot chromosome distribution for normal genes
barplot(table(genes_inf_2$Chromosome), main = "Normal Chromosome Distribution", xlab = "Chromosome", ylab = "Frequency")

### d) Create a bar chart showing the percentages of the DEGs that are upregulated (higher) in Tumor samples and down regulated (lower) in Tumor samples
# count number of higher_expressed genes
higher_expressed_counts <- table(sub_gene_exp$higher_expressed)
# create a bar plot for the counts
barplot(
  higher_expressed_counts,
  main = "Distribution of Higher Expressed Genes",
  xlab = "Higher Expressed",
  ylab = "Count",
  col = "blue",
  las = 1 # Rotate axis labels horizontally for better readability
)

### e) Use the raw data from part 1b to create a heatmap visualizing gene expression by sample
# Load required libraries
library(ggplot2)
library(reshape2)
# prep data for heatmap
heatmap_data <- gene_exp[, sapply(gene_exp, is.numeric)]  # get only numeric columns
rownames(heatmap_data) <- gene_exp$Probe_ID  # set Probe_ID as row names
# change data to a long format for ggplot2
heatmap_data_long <- melt(as.matrix(heatmap_data), varnames = c("Gene", "Sample"), value.name = "Expression")
# create the heatmap
ggplot(heatmap_data_long, aes(x = Sample, y = Gene, fill = Expression)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Expression Level") +  # Use viridis color scale
  labs(
    title = "Heatmap of Gene Expression by Sample",
    x = "Sample",
    y = "Gene"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    axis.text.y = element_blank()  # Optionally hide y-axis labels if there are many genes
  )

### f) Use the same data from the previous part to create a clustermap visualizing gene expression by sample
# Load required library
install.packages("pheatmap")
library(pheatmap)
# Subset the data for the first 100 genes and 10 samples
reduced_heatmap_data <- heatmap_data[1:100, 1:10]
# Create a clustered heatmap
pheatmap(
  reduced_heatmap_data,
  cluster_rows = TRUE,           # Cluster rows (genes)
  cluster_cols = TRUE,           # Cluster columns (samples)
  main = "Clustermap of Gene Expression by Sample",
  legend = TRUE,
  legend_labels = "Expression Level",
  fontsize_row = 6,              # Adjust row font size for readability
  fontsize_col = 8,              # Adjust column font size
  angle_col = 45                 # Rotate column labels
)
write.csv(sub_gene_exp, "significant_genes_output.csv", row.names = FALSE)
### g) Write a few sentence explaining the findings of your analysis, feel free to reference any of visualizations
# My answer in in the readme file

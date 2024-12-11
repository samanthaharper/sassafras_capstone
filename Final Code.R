#Install packages and load library

pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
rm(list = ls(all = TRUE))

pacman::p_load(tidyverse,
               vcfR,
               ggplot2, 
               kableExtra,
               BiocManager,
               DESeq2,
               pheatmap,
               RColorBrewer,
               vsn,
               gridExtra,
               ggrepel,
               e1071,
               caret,
               randomForest,
               ranger,
               multtest,
               janitor
)

#Load VCF data
vcf <- read.vcfR( "Data/SNP.vcf", verbose = FALSE )

#examine VCF object
vcf

#Convert data to tibble
data <- vcfR2tidy(vcf)

#Isolate Counts Data
gdata <- data$gt

#Load Metadata
sass_metadata <- read_csv("Data/sass_metadata.csv", col_types = cols(Longitude = col_number(), Latitude = col_number(), `Number of Samples` = col_integer(), Temperature = col_number(), Precipitation = col_integer(), `Altitude (m)` = col_integer()))

#Clean Variables
#Change variables so they don't start with a number
gdata$chrom_info <- paste("S", gdata$ChromKey, "_", gdata$POS)
#Remove white space
gdata$chrom_info <- gsub(" ", "", gdata$chrom_info)
#Remove other variables
gdata <- gdata %>% select(chrom_info, Indiv, gt_AD)
#Pivot wider to fit 
gdata <- gdata %>% pivot_wider(names_from = Indiv, values_from = gt_AD)

#set up matrix
gdata <- as.data.frame(gdata)
rownames(gdata) <- gdata$chrom_info
gdata <- gdata[ ,-1]

#create proper metadata
indivs <- unique(data$gt$Indiv)
indivs <- as.data.frame(indivs)
indivs$id <- substr(indivs$indivs, 1, 2)
sass_metadata$id <- substr(sass_metadata$Code, 1, 2)
joined_meta <- inner_join(sass_metadata, indivs, by = 'id')
joined_meta <- subset(joined_meta, select = -c(Code))

#Join metadata with sample names
joined_meta <- as.data.frame(joined_meta)
rownames(joined_meta) <- joined_meta$indivs
joined_meta <- subset(joined_meta, select = -c(indivs))

#rename column to exclude the space
joined_meta$Altitude <- joined_meta$`Altitude (m)`

#Create Design Object
dsgnObject <- DESeqDataSetFromMatrix(countData = gdata, 
                                     colData = joined_meta,
                                     design = ~ Longitude + Latitude + Temperature + Altitude)
dim(dsgnObject)

#Set up VST
runVST <- function(dsgnObject, blind, fitType, makePlot = TRUE, writeTable = FALSE, writeRData = FALSE) {
  ## Perform the VST
  
  # Check if the fitType is the regularized log:
  if(fitType == "rlog") {
    vsData <- rlog(dsgnObject, blind = blind)
  }
  ## Otherwise:
  else {
    vsData <- varianceStabilizingTransformation(dsgnObject, 
                                                blind = blind, 
                                                fitType = fitType)
  }
  
  if(makePlot == TRUE) {
    # Plot the effect of the VS transform:
    p1 <- meanSdPlot(assay(dsgnObject), plot = F)
    p1 <- p1$gg + ggtitle("Before Variance Stabilization") + 
      scale_fill_gradient(low = "cadetblue", high = "purple") + 
      theme_bw() + theme(legend.position = "bottom")
    p2 <- meanSdPlot(assay(vsData), plot = F)
    p2 <- p2$gg + ggtitle("After Variance Stabilization") + 
      scale_fill_gradient(low = "cadetblue", high = "purple") + 
      theme_bw() + theme(legend.position = "bottom")
    grid.arrange(p1, p2, nrow=1)
  }
  
  if(writeTable == TRUE) {
    # Write the data for future use, if needed:
    write.table(assay(vsData),
                file = "vst.txt",
                sep="\t", 
                quote=F, 
                row.names=T)
  }
  if(writeRData == TRUE) {
    save(vsData, file="vst_all_timepoints.Rdata")
  }
  return(vsData)
}

#Examine Log-Log Plot
meanCounts <- rowMeans(assay(dsgnObject))      ## Per locus, what is the average expression
varCounts <- apply(assay(dsgnObject), 1, var)  ## Apply the variance function by margin = 1, which is rows

plot(log(varCounts) ~ log(meanCounts), 
     ylab = "Natural-log Variance in Gene Expression", 
     xlab = "Natural-log Mean Expression", 
     main = "\nLog-Log plot of variance by mean for each gene\n should be approximately linear.\n", 
     pch = 16, 
     cex = 0.75)
abline(lm(log(varCounts+0.0001) ~ log((meanCounts+0.0001))), 
       col = "#a8325e", 
       lty = 2, 
       lwd = 2)
#VST 1
#runVST(dsgnObject, blind = FALSE, fitType = "local", makePlot = TRUE, writeTable = FALSE, writeRData = TRUE)

#VST 2
runVST(dsgnObject, blind = FALSE, fitType = "local", makePlot = TRUE, writeTable = FALSE, writeRData = TRUE)

#VST 3
#runVST(dsgnObject, blind = FALSE, fitType = "mean", makePlot = TRUE, writeTable = FALSE, writeRData = FALSE)

#VST 4
#runVST(dsgnObject, blind = FALSE, fitType = "rlog", makePlot = TRUE, writeTable = FALSE, writeRData = FALSE)

#Split data

#Do Splits from Project 1
#So we don't have to reinvent the wheel here
doSplits <- function(vst, algorithm, splitRatio, filterCutoff) {
  ### @vst = vst dataset as extracted from DESeq2
  ### @algorithm = ML algorithm used; currently set up for rf and svm
  ### @splitRatio = the split ratio to employ (training size)
  ### @filterCutoff = the filter cutoff for median number of VST gene counts
  
  ## According to the Valabas et al. (2019) paper, make sure that we are filtering in TRAINING set only! 
  
  # Extract the VST data and transpose
  tVST <- t(assay(vst))
  
  # We do have gene names, e.g., TRNAV-UAC that are malformed for ranger and randomForest. We will fix that before proceeding:
  for (c in 1:ncol(tVST)) {
    colName <- colnames(tVST)[c]
    colName <- gsub("-", "_", colName)
    colName -> colnames(tVST)[c]
  }
  
  ## Add the metadata as columns & merge
  df1 <- cbind(colData(dsgnObject)[1], colData(dsgnObject)[2], colData(dsgnObject)[4], colData(dsgnObject)[5], colData(dsgnObject)[7])       ## We don't need the size factors
  tVST <- merge(tVST, df1, by = "row.names")
  
  ## The merge turns back into a dataframe and removes the sample names from the rows; let's put them back:
  rownames(tVST) <- tVST[,1]
  tVST <- tVST[,-1]
  
  if(algorithm == "svm") {
    ## Make the factors unordered
    tVST <- tVST %>% 
      mutate_if(is.ordered, factor, ordered = FALSE)
  }
  
  ## Create the data partitions
  ind <- createDataPartition(y = tVST[, c("id")],     ## Treatment is evenly distributed
                             p = splitRatio,                    ## % into training
                             list = FALSE)                      ## don't return a list
  train <- tVST[ind, ]
  test <- tVST[-ind,]
  
  ## Now apply the filtering:
  # Calculate row medians of VST gene counts
  medians <- rowMedians(assay(vst))
  
  # Filter the features out of train:
  train <- train[, medians > filterCutoff]  
  print(paste0("After filtering, the number of genes remaining in the dataset are: ", ncol(train)))
  
  splits <- list(train, test)
  return(splits)
}

# Load in data
load("vst_all_timepoints.Rdata")

#Split data
#Do before processing to avoid data leakage - adjusted splitRatio to allow for adequate data in the test set
splits <- doSplits(vst = vsData, algorithm = "rf", splitRatio = 0.8, filterCutoff = 5)
train <- splits[[1]]
test <- splits[[2]]

#Create a genlight object
x <- vcfR2genlight(vcf)
x

#

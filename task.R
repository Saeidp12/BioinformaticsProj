if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("BiocParallel")

library(BiocParallel)

# On windows:
SnowParam()

library(magrittr)
library(dplyr)
library(data.table)
library(ggplot2)

library(Biobase)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(enrichplot)
library(DESeq2)
require(DOSE)
library(pathview)


############## reading the counts data
counts = data.frame(read.table('E:/2022/Freelance/Jobs/002-200/count.txt', header=T))
#Then we omit additional characters in the column names of the counts data frame
colnames(counts) <- sapply(colnames(counts), function(x){sub(".counts","",x)})
colnames(counts) <- sapply(colnames(counts), function(x){sub("X","",x)})
counts[,1:3] %>% head
# Then we set the first column (gene) as the rownames of counts
rownames(counts) = counts$gene
counts = counts[,-1]


# reading the pheno data
pheno = data.frame(read.table('E:/2022/Freelance/Jobs/002-200/pheno.txt', header=T))
pheno %>% head
pheno_2 = pheno
# Then we set the first column of phenotypic data as the rownames
rownames(pheno_2) = pheno_2$sample
pheno_2 = pheno_2[,-1]

# here you can see the levels of lithium and diagnosis
levels(factor(pheno$lithium)) # 0: non-lithium user, 1: lithium user
levels(factor(pheno$diagnosis)) # BP1: patients treated with lithum, BP2: patients not treated with lithium 


# We create a function to group BP patients together. 
bp_test = function(x){
  if(x=='Control'){
    return(0)
  }
  if(x == 'BP1' || x == 'BP2'){
    return(1)
  }
  else{
    return(NA)
  }
}
# Then we create another variable that indicates if a sample observation
# is a patient or control observation
# bipolar disorder => 1: bipolar patient, 0: Control
pheno_2[,'Bipolar'] = factor(unlist(lapply(pheno_2[,2], bp_test)))


########### Preprocessing

# first we numeralize the age column
pheno_2$age = factor(as.numeric(pheno_2$age))

# then we factorize diagnosis, sex, lithium and bipolar
# sex
pheno_2$sex = factor(pheno_2$sex, levels=c('F','M'))
# diagnosis
pheno_2$diagnosis = factor(pheno_2$diagnosis, levels=c('BP1', 'BP2', 'Control'))
# lithium
pheno_2$lithium = factor(pheno_2$lithium, levels=c(0,1))

# the order of rows in pheno should be equal to the order of columns in counts
all(rownames(pheno_2) == colnames(counts)) # Errors 
counts_2 = counts[,rownames(pheno_2)]
all(rownames(pheno_2) == colnames(counts_2))
# Differential Expression Analysis while correcting for effects of age, sex, and lithium
dds <- DESeqDataSetFromMatrix(countData = counts_2,
                              colData = pheno_2,
                              design = ~ age + sex + lithium + Bipolar)
dds

# Pre-filtering
keep <- rowSums(counts(dds)) >= 10000
dds <- dds[keep,]
dim(dds)

########### DEG
# here we should use parallel computing to get the results faster
dds <- DESeq(dds, parallel=T)
res <- results(dds, alpha=0.05)
res

########### Plot dispersion vs normalized counts
plotDispEsts(dds)

########### Histogram of P-values
# with ggplot
ggplot(as(res, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)


# with barplot
use <- res$baseMean > metadata(res)$filterThreshold
h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")

barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))


############# MA Plot
plotMA(res, ylim=c(-2,2))


# We can order our results table by the smallest p value:


plot(-log10(res$padj), res$log2FoldChange)

############ List top 25
top25_genes = rownames(res[order(res$padj)[1:25], ]) # in Ensembl Format


############ Gene Ontology
uniprotKeys <- head(keys(org.Hs.eg.db, keytype="UNIPROT"))
cols <- c("ENSEMBL", "PATH")
AnnotationDbi::select(org.Hs.eg.db, keys=uniprotKeys, columns=cols, keytype="UNIPROT")


humanGeneUniverse <- as.character(unique(AnnotationDbi::select(org.Hs.eg.db,
                                                               keys=keys(org.Hs.eg.db),
                                                               column="ENSEMBL")$ENSEMBL))
head(humanGeneUniverse)


# Define genes (symbols) of interest that we want to analyze
my_genes = unlist(lapply(top25_genes, substr, 1, 15))


# put the genes that are in my_genes and humanGeneUniverse
geneList <- factor(as.integer(humanGeneUniverse %in% my_genes))
names(geneList) <- humanGeneUniverse

geneList %>% head

# Create a topGO Object
GOdata <- new("topGOdata", 
              ontology="BP", 
              allGenes=geneList, 
              nodeSize=5, # the minimum number of genes for each GO term.
              annotationFun=annFUN.org, 
              mapping="org.Hs.eg.db", 
              ID="ensembl")


GOresult <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultTable <- GenTable(GOdata, GOresult, topNodes = 50, numChar=200)
resultTable

#pValue.classic <- score(GOresult)
#head(pValue.classic)

############# GO Hierarchy
showSigOfNodes(GOdata, score(GOresult), firstSigNodes = 9, useInfo = 'all')

# Gene Set Enrichment Analysis
# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"

#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)


######### KEGG Pathways
# we want the log2 fold change 
original_gene_list = res$log2FoldChange
names(original_gene_list) = unlist(lapply(rownames(res), substr, 1, 15))

# omit any NA values 
gene_list = na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)


# It's normal to produce warnings
gse <- gseGO(geneList=gene_list, 
             #ont ="BP", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             #minGSSize = 3, 
             #maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             #verbose = TRUE, 
             OrgDb = 'org.Hs.eg.db', 
             pAdjustMethod = "BH")


# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb='org.Hs.eg.db')

# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

# Create a new dataframe having only the genes that were successfully mapped
df2 = res[unlist(lapply(rownames(res), substr, 1, 15)) %in% dedup_ids$ENSEMBL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID


# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)


kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               keyType       = "kegg")

head(kk2, 10)


enrichplot::dotplot(kk2, showCategory = 2, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)


############ Top 100
res_ordered = res[order(res$padj, decreasing=TRUE),]
rownames(res_ordered) = unlist(lapply(rownames(res_ordered), substr, 1, 15))

top100 = res_ordered[1:100,]

counts_100 = counts_2[rownames(top100),]
rownames(dds) = unlist(lapply(rownames(dds), substr, 1, 15))
dds100 = dds[rownames(top100,)]
############ Classification

#BiocManager::install("MLSeq")
library(MLSeq)



# Train Test Split
n <- ncol(dds100) # number of samples
p <- nrow(dds100) # number of features=
class = data.frame(Bipolar=as.numeric(dds100$Bipolar), row.names = colnames(dds100))

# number of samples for test set (30% test, 70% train).
nTest <- ceiling(n*0.3)
set.seed(2)
ind <- sample(n, nTest, FALSE)

# train and test split
data.train <- as.matrix(counts_100[,-ind]+1)
data.test <- as.matrix(as.matrix(counts_100[,ind]+1))
classtr <- factor(class[-ind,])
classtr <- data.frame(condition=classtr, row.names=colnames(data.train))
classts <- factor(class[ind,])
classts <- DataFrame(condition=classts, row.names=colnames(data.test))

# Design = ~ 1 indicates that there is no need to create a differential testing design
# like before
data.trainS4 <- DESeqDataSetFromMatrix(countData=data.train, colData=classtr,
                                       design= ~ 1)
data.testS4 <- DESeqDataSetFromMatrix(countData=data.test, colData=classts,
                                      design= ~ 1)
######## svm radial
# The control variable sets the method for cross validation, number of partitions and how
# many times to repeat the process. classProbs gives us the probabilities that an observation
# belongs to different classes
ctrl <- trainControl(method = "repeatedcv",
                     number=2,
                     repeats=2,
                     classProbs = FALSE)
fit.svm = classify(data=data.trainS4, method='svmRadial', preProcessing = 'deseq-rlog',
                    control = ctrl)
show(fit.svm)
confusionMat(fit.svm)

availableMethods()

######### lda
fit.lda = classify(data=data.trainS4, method='lda', preProcessing = 'deseq-rlog',
                   control = ctrl)
show(fit.lda)
confusionMat(fit.lda)

######### glm
fit.glm = classify(data=data.trainS4, method='glm', preProcessing = 'deseq-rlog',
                   control = ctrl)
show(fit.glm)
confusionMat(fit.glm)

######### neural networks
fit.nnet = classify(data=data.trainS4, method='nnet', preProcessing = 'deseq-rlog',
                   control = ctrl)
show(fit.nnet)
confusionMat(fit.nnet)

######### Decision trees
fit.rpart = classify(data=data.trainS4, method='rpart', preProcessing = 'deseq-rlog',
                   control = ctrl)
show(fit.rpart)
confusionMat(fit.rpart)

######### Gradient Boost Machine
fit.gbm = classify(data=data.trainS4, method='gbm', preProcessing = 'deseq-rlog',
                   control = ctrl)
show(fit.gbm)
confusionMat(fit.gbm)


################# Parameter Tuning
set.seed(2)
# Here tuneLength sets the number of levels for parameter tuning
# ref="1" means that our reference class in the target variable (Bipolar) is 1 which
# equivalent for being a patient of bipolar disorder
######### Support vector machines with radial basis function kernel
fit.svm <- classify(data = data.trainS4, method = "svmRadial",
                    preProcessing = "deseq-vst", ref = "1", tuneLength = 10,
                    control = ctrl)
show(fit.svm)
######### LDA 
fit.lda <- classify(data = data.trainS4, method = "lda",
                    preProcessing = "deseq-vst", ref = "1", tuneLength = 10,
                    control = ctrl)
show(fit.lda)
######### Generalized Linear Model
fit.glm <- classify(data = data.trainS4, method = "glm",
                    preProcessing = "deseq-vst", ref = "1", tuneLength = 10,
                    control = ctrl)
show(fit.glm)
######### Neural Network
fit.nnet <- classify(data = data.trainS4, method = "nnet",
                    preProcessing = "deseq-vst", ref = "1", tuneLength = 10,
                    control = ctrl)
show(fit.nnet)
######### Decision Tree with caret
fit.rpart <- classify(data = data.trainS4, method = "rpart",
                    preProcessing = "deseq-vst", ref = "1", tuneLength = 10,
                    control = ctrl)
show(fit.rpart)
######### Gradient Boost Machines
fit.gbm <- classify(data = data.trainS4, method = "gbm",
                    preProcessing = "deseq-vst", ref = "1", tuneLength = 10,
                    control = ctrl)
show(fit.gbm)













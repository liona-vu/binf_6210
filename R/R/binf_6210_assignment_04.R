##### Assignment 04 

####  Course: BINF6210

####  By: Liona Vu

####  gihub link https://github.com/liona-vu/binf_6210

####  Last Updated: 05 Dec 2025

#### BACKGROUND INFO  ####
## The Muridae group is a taxonomic class containing well-known research model organisms such as both Mus musculus (house mouse) and the Rattus norvegicus (brown rat). In fact, several genes such as Cytb a mitochondrial gene, Irbp, a gene known for its function in vision, and Rag1, a gene involved in the immune function, are some key genes that have been well-studied and thus their roles have been robustly established in several biological processes. 

### RATIONALE
## Due to their well-known role, it can be said that these genes should have many deposits on NCBI. Therefore, we can use these genes to generate robust classifiers, either for both genes or species.

#### PROJECT OBJECTIVES 
## The objective is to compare different kinds of classifiers using k-mers features on gene sequences to determine if each can accurately train and predict on unseen data.

#### HYPOTHESIS
## I predict genes classifiers, with kmer dinucleotide features, can accurately train data and be used as gene and species classifiers.

#### Step 1: LOADING IN OUR REQUIRED PACKAGES ####
library(tidyverse)
library(rentrez)
library(Biostrings)
library(randomForest)
## Uncomment and run the following lines first if caret or rPOC are not installed
## install.packages("caret")
## install.packages("pROC")
library(caret)
library(pROC)
library(reshape2)
library(vegan)
library(ggplot2)
theme_minimal()

#### STEP 2: LOADING IN OUR DATA WITH RENTREZ PACKAGE ####
## Creating a function requiring class and gene input to download data from NCBI and with all IDs
function_search <- function(taxon, gene_name) {
  ## Input checks to ensure that taxon and gene name inputs must be in string format
  if(!is.character(taxon)) {
    stop("Taxon input must be in string format!")
  }
  if(!is.character(gene_name)) {
    stop("Gene_name input must be in string format!")
}
## Searching for sequences. Added a bp limit of 1000-2000 bp or else downloading will take forever, some sequences on NCBI are over 100,000,000 bp long! (I actually went and checked, some sequences were whole genome shotgun sequenced)
  search_term <- paste((taxon), "[ORGN] AND ", gene_name, 
                       " AND 1000:1500[SLEN]", sep = "") 
 ## search_temp <- entrez_search(db = "nuccore", term = search_term)
 ## max_hits <- search_temp$count
  search_hits <- entrez_search(db = "nuccore", term = search_term, retmax = 9999,
                               use_history = TRUE)
  return(search_hits)
}

## Quick sanity check to ensure that the function works as intended with a small example taxon and gene
fn_check <- function_search(taxon = "Mus musculus", gene_name = "Cytb[gene]") 
fn_check$retmax #1624 instead of having the default retmax and ID value of 20
class(fn_check) # esearch and list
## Also performed a sanity check on NCBI and it also returned 1624 search results

## No longer needed so remove
rm(fn_check)

## Enabling a vector of gene names to be downloaded from NCBI 
genes_terms <- c("Cytb[gene] and complete CDS", #complete CDS for stringency, else it will give over 20,000 results!
                 "Irbp[gene]",
                 "Rag1[gene]")

## Creating a loop to download all the NCBI datasets with genes of interest
ncbi_muridae <- list()
for (i in genes_terms) {
  cat("***** Downloading", i, "nucleotide database from NCBI ***** \n")
  print(function_search(taxon = "Muridae", gene_name = i))
  ncbi_muridae[[i]] <- function_search(taxon = "Muridae", gene_name = i)
}

## Checking to ensure that the looping was successful in grabbing all the information needed
## Performing a quick preview check for the Muridae taxonomic class by taking a look at the NCBI summary
class(ncbi_muridae) #is a list
length(ncbi_muridae) #3, which is expected since we have 3 genes
ncbi_muridae$`Cytb[gene] and complete CDS`$count #740
ncbi_muridae$`Irbp[gene]`$count #1804
ncbi_muridae$`Rag1[gene]`$count #854

## Since the console throws errors if we download too many IDs from NCBI and results are in the high hundreds, will randomly sample 200 from all gene ID list. Create a function to generate 200 random ids for each gene
function_set_seed <- function(ids, seed) {
  set.seed(seed)
  sampling_seed <- sample(ids, size = 200)
}

## Running the function_set_seed for each gene
cytb_seed <- function_set_seed(ncbi_muridae$`Cytb[gene] and complete CDS`$ids, seed = 314)
irbp_seed <- function_set_seed(ncbi_muridae$`Irbp[gene]`$ids, seed = 314)
rag1_seed <- function_set_seed(ncbi_muridae$`Rag1[gene]`$ids, seed = 314)

## Fetching data for Cytb gene
muridae_fetch_cytb <- entrez_fetch(db ="nuccore", id = cytb_seed, rettype = "fasta")

## Repeating the same for the Irpb and rag1 genes
muridae_fetch_irbp <- entrez_fetch(db ="nuccore", id = irbp_seed, rettype = "fasta")
muridae_fetch_rag1 <- entrez_fetch(db ="nuccore", id = rag1_seed, rettype = "fasta")

## Checking if the class is correct 
class(muridae_fetch_cytb) # long character class
class(muridae_fetch_irbp)
class(muridae_fetch_rag1)

## Writing the fetch variables to actual files on the computer to keep a copy of the data. Commented out to prevent overwriting.
## write(x = muridae_fetch_cytb, file = "../data/muridae_cytb.fasta", sep = "\n")
## write(x = muridae_fetch_irbp, file = "../data/muridae_irbp.fasta", sep = "\n")
## write(x = muridae_fetch_rag1, file = "../data/muridae_rag1.fasta", sep = "\n")

## Checking directory to ensure that the files are properly written to the data file
list.files(path = "../data")

#No longer need so remove, declutter our enviroment
rm(muridae_fetch_cytb,muridae_fetch_irbp, muridae_fetch_rag1, cytb_seed, irbp_seed, rag1_seed, genes_terms)

## Read it back in as DNA StringSet using the readDNAStringSet() function from Biostrings package
muridae_cytb <- readDNAStringSet("../data/muridae_cytb.fasta")
muridae_irbp <- readDNAStringSet("../data/muridae_irbp.fasta")
muridae_rag1 <- readDNAStringSet("../data/muridae_rag1.fasta")

## Checking if the new files are now Biostrings objects
class(muridae_cytb) 
class(muridae_irbp)
class(muridae_rag1) # all are Biostrings, DNAstringset class

## Checking the top few names of the new file
head(names(muridae_cytb)) 
length(muridae_cytb) # 200 which matches the number sampled earlier
length(muridae_irbp) # 200
length(muridae_rag1) # 200

## Created a dataframe with the names of the samples and the associated sequences, and its associate gene named to be used later
df_cytb <- data.frame(samples = c(names(muridae_cytb)), gene = "Cytb", sequences = paste(muridae_cytb))
df_irbp <- data.frame(samples = c(names(muridae_irbp)), gene = "Irbp", sequences = paste(muridae_irbp))
df_rag1 <- data.frame(samples = c(names(muridae_rag1)), gene = "Rag1", sequences = paste(muridae_rag1))

## Combining the above 3 dataframes to one singular dataframe
df_combined <- bind_rows(df_cytb, df_irbp, df_rag1)

## Check for our data class, dimensions, and sample names
dim(df_combined) # dimensions has 600 rows with 3 columns. The 600 matches with the sum of all combined IDs numbers 200+200+200
class(df_combined) #is a dataframe
table(df_combined$gene) #Checking the genes and whether each column is 200, which it is
names(df_combined) # Checking out the names

## Remove to declutter environment
rm(df_cytb, df_irbp, df_rag1)

#### STEP 3: FILTERING THE DATA ####
## Filter in case there are any Na values in names, species_name, genes, and sequence columns
conflicted::conflicts_prefer(dplyr::filter)
df_combined_filtered <- df_combined %>%
         filter(!is.na(samples), 
                !is.na(gene),
                !is.na(sequences))
        
## Sanity check for Nas
sum(is.na(df_combined_filtered$samples)) #Returns 0
sum(is.na(df_combined_filtered$gene)) #Returns 0
sum(is.na(df_combined_filtered$sequences)) #Returns 0

## Check to see if ambiguous nucleotides are present in the sequences column which are marked by "N" and how many are present 
sum(grepl(pattern = "N", df_combined_filtered$sequences)) # Some sequences have Ns
str_count(string = df_combined_filtered$sequences, pattern = "N") #some sequences have really high N counts

## Since there are some sequences that contain "Ns" and some sequences have a high count of Ns, I need to determine how many Ns are present and to filter out these sequences that have more than 5% Ns of the entire sequence. Also, trimming off potential Ns that are leading and trailing the sequences.
df_seq <- df_combined_filtered %>%
         mutate(sequences_2 = str_remove(sequences, "^[-N]+")) %>%
         mutate(sequences_2 = str_remove(sequences_2, "[-N]+$")) %>%
         mutate(sequences_2 = str_remove_all(sequences_2, "-+")) %>%
         filter(str_count(sequences_2, "N") <= 0.05 * str_count(sequences))

## Comparing the summary statistics of pre vs post filtering.
summary(nchar(df_combined_filtered$sequences))
summary(nchar(df_seq$sequences_2))
## The stats have stayed the same. A too steep of a change indicates that the filtering was too stringent.

## STEP 4: EXPLORATORY ANALYSIS ####
##Looking at the distribution of how long each gene is
par(mfrow = c(1,3)) ## Allows graph to be graphed side by side

hist(nchar(df_seq$sequences_2[df_seq$gene == "Cytb"]), xlab = "Sequence Length (bp)", main = "Frequencies of Cytb Sequence Lengths") 

hist(nchar(df_seq$sequences_2[df_seq$gene == "Irbp"]), xlab = "Sequence Length (bp)", main = "Frequencies of Irbp Sequence Lengths") 

hist(nchar(df_seq$sequences_2[df_seq$gene == "Rag1"]), xlab = "Sequence Length (bp)", main = "Frequencies of Rag1 Sequence Lengths") 

#Confirming by also plotting qq plot in base R and plotting a line through
qqnorm(nchar(df_seq$sequences_2[df_seq$gene == "Cytb"]), main = "QQ plot of Cytb")
qqline(nchar(df_seq$sequences_2[df_seq$gene == "Cytb"]), col = "red", lwd = 2)

qqnorm(nchar(df_seq$sequences_2[df_seq$gene == "Irbp"]), main = "QQ plot of Irbp")
qqline(nchar(df_seq$sequences_2[df_seq$gene == "Irbp"]), col = "red", lwd = 2)

qqnorm(nchar(df_seq$sequences_2[df_seq$gene == "Rag1"]), main = "QQ plot of Rag1")
qqline(nchar(df_seq$sequences_2[df_seq$gene == "Irbp"]), col = "red", lwd = 2)

## From some histograms, particularly Cytb and Rag1, it can be said that the gene sequences are not normally distributed. This is to be expected since genes sequences fall within a certain window of sequence lengths. Some sequences were on the larger end, due to potentially having different strains of certain Muridae species. This was also seen in the qqplots, where the data is not normally distributed.

#### STEP 5: CALCULATING SEQUENCE FEATURES AND TRAINING CLASSIFICATION MODEL WITH RANDOM FORESTS ####
## Converting sequences into the Biostrings class
df_seq$sequences_2 <- DNAStringSet(df_seq$sequences_2)
class(df_seq$sequences_2) #is a Biostrings class

## Calculating the kmer frequency (dinucleotide) for each column
df_seq <- cbind(df_seq, as.data.frame(dinucleotideFrequency(df_seq$sequences_2, as.prob = TRUE)))
View(df_seq)

## Next, I will generate a training and validation set to train our model.
## Converting Biostring format back to character data for tidyverse functions.
df_seq$sequences_2 <- as.character(df_seq$sequences_2)
class(df_seq$sequences_2)

## Setting seed to produce randomization, grouping by genes, and want to make this analysis reproducible in creating a validation set which will be 25% of our dataset
set.seed(123)
df_validation <- df_seq %>%
  group_by(gene) %>%
  sample_frac(0.25)

## Ensuring that we have 25% of our dataset, should be around 50 since the filtered dataset had roughly a little less than 200 counts per gene
table(df_validation$gene)

## Doing simple math to check if it is really 25% of our dataset for Cytb, Irbp, and Rag1 genes
sum(df_validation$gene == "Cytb")/sum(df_seq$gene == "Cytb")*100 # 25.1%
sum(df_validation$gene == "Irbp")/sum(df_seq$gene == "Irbp")*100 # 25 %
sum(df_validation$gene == "Rag1")/sum(df_seq$gene == "Rag1")*100 # 25 %

## Need to ensure that the each samples downloaded from NCBI are unique in order to filter out rows from the validation set for our training set, the following should indicate TRUE if the dim and length of unique samples are the same
dim(df_seq)[1] == length(unique(df_seq$samples)) 
#each row has a unique "identity"

## Setting seed to make this analysis reproducible and creating a training set
set.seed(321)
df_training <- df_seq %>%
  filter(!samples %in% df_validation$samples) %>%
  group_by(gene)
  
## Ensuring that we have actually have 75% of our dataset
table(df_training$gene)

## Doing simple math to check if it is really 75% of our dataset for Cytb, Irbp, and Rag1 genes
sum(df_training$gene == "Cytb")/sum(df_seq$gene == "Cytb")*100 # 74.87%
sum(df_training$gene == "Irbp")/sum(df_seq$gene == "Irbp")*100 # 75 %
sum(df_training$gene == "Rag1")/sum(df_seq$gene == "Rag1")*100 # 75 %

## Build a classifier for each gene using with dinucleotide kmers. The response variable is gene, we are trying to predict which gene a sequence belongs to
conflicted::conflicts_prefer(base::as.factor)
gene_classifier_kmer <- randomForest::randomForest(x = df_training[,5:20], y = as.factor(df_training$gene), ntree = 50, importance = TRUE)

## Looking at results, relative importance, left out of bag, error rate, and confusion matrix
gene_classifier_kmer
gene_classifier_kmer$oob.times
gene_classifier_kmer$err.rate
gene_classifier_kmer$confusion

## Now I can run the classifier to the unseen data
prediction_validation_kmers <- predict(gene_classifier_kmer, df_validation[, c(2, 5:20)])
table(observed = df_validation$gene, predicted = prediction_validation_kmers)
## From the random forest model, it seems that the dinucleotide kmer model classifier did a great job at correctly identifying all 3 genes.

## Now, I want to see which feature kmer feature contributed the most to the random forest model. Plotting a kmer frequency heatmap. Need to first convert to kmer counts
kmer_counts_df <- as.data.frame(dinucleotideFrequency(DNAStringSet(df_training$sequences_2), as.prob = TRUE))
head(kmer_counts_df)

## Adding gene names to dataframe
kmer_counts_df$gene <- df_training$gene 

## Converting to longer dataframe format to plot graph in ggplot
kmer_long <- melt(kmer_counts_df)

## Plotting heatmap with ggplot
ggplot(kmer_long, aes(x = variable, y = gene, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  #scale_fill_gradient(low = "black", high = "red") +
  labs(title = "Heatmap of K-mer Dinucleotide Contributions", x = "Dinucleotides", y = "Genes") +
  theme_minimal()
## From looking at the heatmap, it seems that different kmers proportions contributed for all 3 genes classifications.

## Saving plot, uncomment to save
## today_date = Sys.Date()
## ggsave(filename = paste0("../figs/", as.character(today), "_heatmap_kmer_random_forest.png"))

## Plotting ROC curve
## Create function to calculate roc probabilities
calculate_roc_probs <- function(negative, positive) {
  probs <- predict(gene_classifier_kmer, df_training[, c(2, 5:20)], type = "prob", levels =c(as.character(negative), as.character(positive)))[,2]
  return(probs)
}
#Calculating the rocs probabilities
predict_prob_cytb <- calculate_roc_probs(Irbp, Cytb)
predict_prob_irbp <- calculate_roc_probs(Cytb, Irbp)
predict_prob_rag1 <- calculate_roc_probs(Irbp, Rag1)

## Build a roc curve
roc_obj_cytb <- roc(df_training$gene, predict_prob_cytb, levels =c("Irbp", "Cytb"))
roc_obj_irbp <- roc(df_training$gene, predict_prob_irbp, levels =c("Cytb", "Irbp"))
roc_obj_Rag1 <- roc(df_training$gene, predict_prob_rag1, levels =c("Irbp", "Rag1"))

#Plot the ROC curve side by side
par(mfrow = c(1,3))
plot(roc_obj_cytb, main = "ROC Curve of Cytb", col = "pink", xlab = "False positive rate (FPR)", ylab = "True Positive Rate (TPR)", lwd = 3)
plot(roc_obj_irbp, main = "ROC Curve of Irbp", col = "turquoise", xlab = "False positive rate (FPR)", ylab = "True Positive Rate (TPR)", lwd = 3)
plot(roc_obj_Rag1, main = "ROC Curve of Rag1", col = "green", xlab = "False positive rate (FPR)", ylab = "True Positive Rate (TPR)", lwd = 3)
#ROC graph identifies all genes as perfect classifier

## Other plots that I was interested in, for example the out of bag error rate as the number of trees increased.
par(mfrow = c(1,1))
plot(gene_classifier_kmer, main = "Random Forest Out of Bag Error Rate of each Gene")

## Add a legend 
legend("topright",
       legend = c(colnames(gene_classifier_kmer$err.rate)),
       col = 1:(ncol(gene_classifier_kmer$err.rate)),
       lty = 1)
## Seems that error rate decreased drastically after about 5 trees.

#### STEP 6: TRANING AND CLASSIFICATION OF OTHER METHODS USING CARET PACKAGE ####
## After using random forest to train and classify our data, I now want to see if this holds true for other types of classifications. One example is called partial least squares, (PLS).

## First, must get rid of columns that contain string characters that will not be used in the model else, the function train in caret will not work, only works with numerical data.
numeric_cols <- sapply(df_training, is.numeric)

## From the above, I can see that the first 5 columns are not numeric as indicated by FALSE while the numeric columns are true. I also want to keep the gene column. Creating a new dataframe to keep gene column while removing string columns.
df_training_pls <- df_training[, c("gene", names(df_training)[numeric_cols])]
class(df_training_pls)

## Training the model with a partial least squares discriminant analysis (PLSDA)
plsFit <- caret::train(gene ~ ., data = df_training_pls, method = "pls",
  ## Center and scale the predictors for the training
  ## set and all future samples.
  preProc = c("center", "scale"))

## Checking the output of the PLS model
plsFit
plsFit$results #Accuracy increases as the number of components increases up to a maximum of 2
plsFit$resample # Default boostrapping is 25, so it shows the accuracy value of each resampling

## Now, will see how well our model does on unseen (validation data)
plsClasses <- predict(plsFit, newdata = df_validation)
str(plsClasses) #Checking the structure of our validation 

## Computes the confusion matrix and the statistics of the PLS model fit
confusionMatrix(data = plsClasses, as.factor(df_validation$gene))
## Accuracy was also at 100%, similar to the random forest model

## Simple ggplot to showcase the relationship between the performance values and number of PLS components
ggplot(plsFit) +
  theme_minimal() +
  labs(title = "Relationship between Number of Components and Accuracy",
       x = "Number of Components",
       y = "Accuracy")

## Warning message is due to the caret package using the old aes_string() instead of the new aes(), should not affect the analysis (and should probably report this issue to the developers...)

## Saving plot, uncomment to save
## ggsave(filename = paste0("../figs/", as.character(today_date), "_PLS_accuracy_component_plot.png"))

## See which k-mer variable / features contributed to each machine learning model for each gene by using the varImp function from the caret package
varImp(plsFit) %>%
  plot(main = "K-mer Variable Importance Plot for PLS Classifier")
## Again, similar to the RF, different kmers distributuons contribute differently to the model.

#### STEP 7: SPECIES LEVEL CLASSIFICATION (ADDITIONAL ANALYSIS) ####
#Since PLS and random forests all have been at 100% classification accuracy on genes at training on data and predicting on unseen data, I want to see if the classifiers can accurately classify Muridae species. Assignment instructions indicate to attempt a more difficult problem.
#Reading in the full set of cytb set with all 749 genes from web ncbi data, and convert to dataframe
df_cytb_2 <- readDNAStringSet("../data/muridae_cytb_740_entries.fasta")
df_cytb_2 <- data.frame(samples = c(names(df_cytb_2)), gene = "Cytb", sequences = paste(df_cytb_2))

#checking out the dataframe
dim(df_cytb_2) #740 rows with 3 columns
summary(df_cytb_2)

## Since our samples column contains the species name but it is pretty messy, will recreate a new column to their actual species name by using regex expressions. Species names always start with a capital letter, followed by a space, then lowercase letters
species_name <- str_extract(string = df_cytb_2$samples, pattern = "[A-Z][a-z]+ [a-z]+")
table(species_name) # Checking the distribution of the species names

## Combining our species name to our dataframe
df_cytb_2 <- cbind(df_cytb_2, species_name)
names(df_cytb_2) #checking if the species name got added
dim(df_cytb_2) #4 columns which is good

## Some entries do not have actual full species name, such as "Aethomys sp". This is not specific enough for the purpose of the species classifier therefore, will filter them out. Also, filter in case there are any Na values in names, species_name, genes, and sequence columns. Also, filtering out ambiguous nucleotide sequences that are greater than 5% and any trailing and leading Ns.
df_cytb_2_filtered <- df_cytb_2 %>%
  filter(!is.na(species_name),
         !is.na(sequences),
         !is.na(gene)) %>%
  filter(!str_detect(species_name, pattern = "[A-Z][a-z]+ sp")) %>%
  mutate(sequences_2 = str_remove(sequences, "^[-N]+")) %>%
  mutate(sequences_2 = str_remove(sequences_2, "[-N]+$")) %>%
  mutate(sequences_2 = str_remove_all(sequences_2, "-+")) %>%
  filter(str_count(sequences_2, "N") <= 0.05 * str_count(sequences))

## Sanity check for Nas and vague species names
sum(is.na(df_cytb_2_filtered$species_name)) #Returns 0
sum(is.na(df_cytb_2_filtered$sequences)) #Returns 0
sum(is.na(df_cytb_2_filtered$gene)) #Returns 0
table(df_cytb_2_filtered$species_name) # No more vague species names
dim(df_cytb_2_filtered) #number of rows decreased, some vague species names were filtered out, also 5th column with species name showed up

## Convert to biostrings class
df_cytb_2_filtered$sequences_2 <- DNAStringSet(df_cytb_2_filtered$sequences_2)
class(df_cytb_2_filtered$sequences_2) #is a Biostrings class

## Calculating the kmer frequency (dinucleotide) for each column
df_seq_cytb <- cbind(df_cytb_2_filtered, as.data.frame(dinucleotideFrequency(df_cytb_2_filtered$sequences_2, as.prob = TRUE)))
View(df_seq_cytb)

## Convert back to character class for tidyverse
df_seq_cytb$sequences_2 <- as.character(df_seq_cytb$sequences_2)
class(df_seq_cytb$sequences_2)

## Setting seed to produce randomization, grouping by species, and want to make this analysis reproducible in creating a validation set.
set.seed(647)
df_validation_cytb <- df_seq_cytb %>%
  group_by(species_name) %>%
  sample_frac(0.25)

## Ensuring that we have roughly 25% of our training dataset
length(df_validation_cytb$species_name)/length(df_cytb_2_filtered$species_name) #23.5 %

## Setting seed to make this analysis reproducible and creating a training set
set.seed(647)
df_training_cytb <- df_seq_cytb %>%
  dplyr::filter(!samples %in% df_validation_cytb$samples) %>%
  dplyr::group_by(gene)

## Ensuring that we have actually have 75% of our dataset
length(df_training_cytb$species_name)/length(df_seq_cytb$species_name) #76.4 %

## Build a classifier for each gene using the dinucleotide kmers. The response variable is species, we are trying to classify species based on Cytb gene sequence
conflicted::conflicts_prefer(base::as.factor)
species_classifier_cytb <- randomForest::randomForest(x = df_training_cytb[, 6:21], y = as.factor(df_training_cytb$species_name), ntree = 50, importance = TRUE)

## Looking at the results, confusion matrix, and the error rate
species_classifier_cytb
species_classifier_cytb$confusion
species_classifier_cytb$err.rate

## Prediction on unseen data
prediction_cytb <- predict(species_classifier_cytb, df_validation_cytb[, c(4, 6:21)])
table(observed = df_validation_cytb$species_name, predicted = prediction_cytb)
## By looking at the table above, it seems that it the model did a good job at predicting species, there were a few misclassified species though such as Apodemus draco misclassified as Apodemus ilex.

## Relevant diagnostic plots for classification such as checking out error rates
plot(species_classifier_cytb, 
     main = "Error rates vs number of trees in species classifier")

## Add a caption to the bottom of the graph 
mtext(text = "Each line represents a different Muridae species", side = 1, adj = 0, line = 4) 
## This is much busier, as compared to the gene classifier graph

## Now I am plotting a heatmap to see which dinucleotide frequency contributed the most for each species classification during training
kmer_cytb_df <- as.data.frame(dinucleotideFrequency(DNAStringSet(df_training_cytb$sequences_2), 
                                                    as.prob = TRUE))
## Checking for the top of dataframe
head(kmer_cytb_df)

## Adding species names to dataframe
kmer_cytb_df$species_name <- df_training_cytb$species_name 

## Convert to longer dataframe format to plot graph
kmer_cytb_long <- melt(kmer_cytb_df)
head(kmer_cytb_long)

## Plot heatmap with ggplot
ggplot(kmer_cytb_long, aes(x = variable, y = species_name, fill = value)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdPu") +
  labs(title = "Heatmap of K-mer Dinucleotide Contributions", x = "Dinucleotides k-mer", y = "Muridae species") 

## Save to figs folder, uncomment to save
## ggsave(filename = paste0("../figs/", as.character(today), "_kmer_heatmap_Muridae_classifier_training.png"))

## Plotting a multidimensional scaling plot to see clustering of each species classification
mds_cytb <- dist(df_training_cytb[,6:21], method = "euclidean")
cmds_cytb <- cmdscale(mds_cytb)

## Plotting with base R plot
plot(cmds_cytb, col = as.factor(df_training_cytb$species_name), pch = 19, xlab = "MDS Dimension 1", ylab = "MDS Dimension 2")

## Performing a permanova test to see the how well the classifier did and whether the kmer classifier is significant. permanova does not assume normality
permanova_test <- adonis2(mds_cytb ~ species_name, data = df_training_cytb)
permanova_test

## Add the statistic to the plot
mtext(paste0("R² value = ", round(permanova_test$R2[1], 3), 
                       "\n p = ", permanova_test$`Pr(>F)`[1],
                        "\n F statistic = ", round(permanova_test$F[1], 3)), side = 3, line = -3, adj = 1)

## Calculate the mean values for each species with aggregate which computes the mean of each species and returns a dataframe
species_mean <- aggregate(cmds_cytb, by = list(df_training_cytb$species_name), FUN = mean)
is.data.frame(species_mean)

## Commented out due to cluttering the plot with all species names. Uncomment it to see the results if you are curious
## text(species_mean$V1, species_mean$V2, labels = species_mean$Group.1, cex = 0.8, font = 2)

## Checking if I can pull the mean values for Mus musculus
species_mean$V1[species_mean$Group.1 == "Mus musculus"]

## Checking which species have high counts
high_counts <- sort(table(df_training_cytb$species_name), decreasing = TRUE)
high_counts

## Since there are too many species to be labeled, only label the top 6.
top_species <- head(names(high_counts), 6)
top_species # Check the top 6 

## Use a loop to label the MDS plot
for (i in top_species) {
  text(species_mean$V1[species_mean$Group.1 == i],
       species_mean$V2[species_mean$Group.1 == i],
       labels = i,
       cex = 0.8, font = 2)
  }
## From the plot, a yellow cluster can be seen on the right corresponding to Mus musculus. Mus triton makes a nice cluster in green and the same results for Rattus norvegicus in magenta.

#### EXTRA CLASSIFIER FOR CURIOSITY ####
#Trying a different classifier called the Regularized Discriminant Analysis (RDA). Need to make our own custom grid to pass through the tuneGrid parameter, with gamma and lambda values since this is what RDA requires to classify. Avoiding 0 for gamma since there will be many warning messages (about 57) that show up.
rdaGrid = data.frame(gamma = seq(0.1, 1.0, length = 5), lambda = 3/4)

## Running our RDA classifier
rda_fit <- caret::train(gene ~ ., 
                        data = df_training_pls, 
                        method = "rda",
                        tuneGrid = rdaGrid,
                        preProc = c("center", "scale"))

## See the classifier results and plots the results
rda_fit
plot(rda_fit)
## Seems like the classifier is 100% accurate at training our data

## Now testing on unseen data to get the confusion matrix
pls_classes_rda <- predict(rda_fit, newdata = df_validation)

## Checking internal structure of model
str(pls_classes_rda)

## Computes the confusion matrix and the statistics of the RDA model fit
confusionMatrix(data = pls_classes_rda, as.factor(df_validation$gene))
## Accuracy was also at 100%, similar to the random forest and PLS

## Now to plot the data, need to run the following to plot ROC curve
rda_probs <- predict(rda_fit, newdata = df_validation, type = "prob")

## Generating labels for our ROC curves calculations
true_labels <- df_validation$gene

## Generating ROC curves for each gene class from the RDA classifier where 'response' is the true labels, 'predictor' is the predicted probabilities for each class, the levels' specifies the negative and positive classes for binary ROC calculation
roc_cytb <- roc(response = true_labels, predictor = rda_probs$Cytb, levels =c("Irbp", "Cytb"))
roc_irbp <- roc(response = true_labels, predictor = rda_probs$Irbp, levels = c("Cytb", "Irbp"))
roc_rag1 <- roc(response = true_labels, predictor = rda_probs$Rag1, levels = c("Cytb", "Rag1"))

## Creating dataframe for ggplot input
df_cytb <- data.frame(TPR = roc_cytb$sensitivities, FPR = 1 - roc_cytb$specificities, Class = "Cytb")
df_irbp <- data.frame(TPR = roc_irbp$sensitivities, FPR = 1 - roc_irbp$specificities, Class = "Irbp")
df_rag1 <- data.frame(TPR = roc_rag1$sensitivities, FPR = 1 - roc_rag1$specificities, Class = "Rag1")

## Plotting each gene in their own graph due to all of them overlapping each other, will be hard to see if graphed at once
ggplot(df_cytb, aes(x = FPR, y = TPR)) +
  geom_line(linewidth = 1, colour = "pink") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(title = "ROC Curve for RDA Classifier (Cytb)",
       x = "False Positive Rate",
       y = "True Positive Rate")

## ggsave(filename = paste0("../figs/", as.character(today_date), "_roc_curve_RDA_for_cytb.png"))

ggplot(df_irbp, aes(x = FPR, y = TPR)) +
  geom_line(linewidth = 1, colour = "turquoise") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(title = "ROC Curve for RDA Classifier (Irbp)",
       x = "False Positive Rate",
       y = "True Positive Rate")
## ggsave(filename = paste0("../figs/", as.character(today_date), "_roc_curve_RDA_for_Irbp.png"))

ggplot(df_rag1, aes(x = FPR, y = TPR)) +
  geom_line(linewidth = 1, colour = "lightgreen") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(title = "ROC Curve for RDA Classifier (Rag1)",
       x = "False Positive Rate",
       y = "True Positive Rate")

## Saving plot, uncomment to save
## ggsave(filename = paste0("../figs/", as.character(today_date), "_roc_curve_RDA_for_rag1.png"))


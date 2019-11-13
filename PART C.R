#Using the genus Daphnia, identify whether the two genes (COI and 16S) yield the same or different phylogenetic hypotheses
#LOAD PACKAGES ----
library(rentrez) #entrez_search, entrez_fetch
library(seqinr) 
library(Biostrings)  #readDNAStringSet
library(stringr) #word, str_extract
library(tidyverse)
library(stringi)
library(ape) #as.DNAbin, dist.dna
library(RSQLite)
library(DECIPHER) #BrowseSeqs, IdClusters
library(phangorn) #NJ
library(dendextend) #tanglegram

#Use functions as part of the code as opposed to repeat the    same block of code in the single script, because functions can reduce duplication of codes, improve clarity of the code, and also able to reuse the code.  a single script, because it reduces duplication of codes, improve clarity of the code, and also able to reuse the code. Here creates a function called clusters to build dendrogram through aligning and performing a distance matrix
clusters <- function(x){
  #covnert the class of the sequences into Biostrings
  x <- DNAStringSet(x)
  #perform alignment using muscle
  Alignment <- DNAStringSet(muscle::muscle(x, log = "log.tx", verbose = T), use.names = TRUE)
  #display aligned sequences through web browser 
  visual <- BrowseSeqs(Alignment)
  #convert the class from the package Biostrings in BioConductor for consistency within ape package. Converting this data type is essential for further downstream analysis.
  dnaBin <- as.DNAbin(Alignment)
  #computes a matrix of pairwise distances from DNA sequences using a model of DNA evolution. 
  distanceMatrix <- dist.dna(dnaBin, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)
  #cluster according to single linkage method
  cluster <- IdClusters(distanceMatrix, method = "single", cutoff= 0.02, showPlot = TRUE, type = "both", verbose = TRUE)
  #combine the results
  result <- c(cluster,visual)
  #Returns results variable to the global environment
  return(result)
}
get_info <-function(trait){
  #how many observations
  trait_num<- length(unique(trait))
  #names of the observations
  trait_names <- unique(trait)
  #combine the results
  result <- c(trait_num,trait_names)
  #Returns results variable to the global environment
  return(result)
}
#Getting Data from NCBI:
search_result <- entrez_search(db = "pubmed", term = "Daphnia")
#search.result #this is a mistake done in previous code the correct variable is search_result
search_result
#search for sequences of Daphnia within two genes COI and 16S
COI.search <- entrez_search(db = "nuccore", term = "(Daphnia[ORGN] AND COI[TITL]) NOT (genome[TITL])", retmax = 310)
gene16S.search <- entrez_search(db = "nuccore", term = "(Daphnia[ORGN] AND 16S[TITL]) NOT (genome[TITL])", retmax = 400)

#dowloded at 5:16 on wednesday 10/30/2019
# to get fasta format
COI.fetch <- entrez_fetch(db = "nuccore", id = COI.search$ids, rettype = "fasta")
gene16S.fetch <- entrez_fetch(db = "nuccore", id = gene16S.search$ids, rettype = "fasta")
#write the original data to hard disk. \n means to insert a newline in the text
write(COI.fetch, "COI.fetch.fasta", sep = "\n")
write(gene16S.fetch, "gene16S.fetch.fasta", sep = "\n")
#Function to read DNAStringSet from the recieved file which was in fasta format
string.SetCOI <- readDNAStringSet("COI.fetch.fasta")
string.Set16S <- readDNAStringSet("gene16S.fetch.fasta")
#2-Build dataframes:
#create  dataframes for COI gene and 16S gene
dfCOI <- data.frame(COI.Title = names(string.SetCOI), COI.Sequence = paste(string.SetCOI))
df16S <- data.frame(gene16S.Title = names(string.Set16S), gene16S.Sequence = paste(string.Set16S))
#word is in stringr package to get 2nd and 3rd words in the Title details
dfCOI$Species.Name <- word(dfCOI$COI.Title, 2L, 3L) 
df16S$Species.Name <- word(df16S$gene16S.Title, 2L, 3L)
#The accession number is in the first word of the details in Title
dfCOI$Uniqe.Identifier <- word(dfCOI$COI.Title, 1L)
df16S$Uniqe.Identifier <- word(df16S$gene16S.Title, 1L)
#dfCOI$Gene.Name <- str_extract(dfCOI$COI.Title, "COI.*") # this gave a wrong result
str_detect(dfCOI$COI.Title, "COI")#check whether all the sequences have "COI" gene name in the column Title , according to the results all of them have COI gene in the COI.Title    column, so it is quite easy to extract the word COI from the COI_Title column. 
str_detect(df16S$gene16S.Title, "16S")

dfCOI$Gene.Name <- str_extract(dfCOI$COI.Title, "COI")
df16S$Gene.Name <- str_extract(df16S$gene16S.Title, "16S")
#Perform preliminary data checking.
hist(str_length(dfCOI$COI.Sequence))#check the histogram of sequence lengths of COI
hist(str_length(df16S$gene16S.Sequence))#check the histogram of sequence lengths of 16S
#remove the samples with extreamly large sequences
df16S_filtered <- df16S %>% 
  filter(str_length(gene16S.Sequence)<600)
hist(str_length(df16S_filtered$gene16S.Sequence))#check the histogram of after removing the samples with extremely large sequences of 16S
#check the unique species names
get_info(dfCOI$Species.Name)
get_info(df16S_filtered$Species.Name)
#By looking at this result, I found that there is a species name other than Daphnia.So when I check the dataset there were 10 samples from Anaompoda sp., therefore, rows from 98:107 which are Anompoda sp. were removed for downstream analysis
#remove outliars
dfCOI_filtered <- dfCOI[-c(98:107), ]

#perform an alignment to the whole dataset, to check whether it contains outliers
#I have already use the function  clusters which I have built at the top of the code. Now it is quite easy as no need to repeat the whole thing to get the cluster. Just need to replace x,  
clusters(dfCOI_filtered$COI.Sequence)
title("Phylogeny of COI gene")

clusters(df16S_filtered$gene16S.Sequence)
title("Phylogeny of 16S gene")



#Using the function get_info, which I have built at the top of the code.Now it is easy as no need to repeat the whole code. Just need to replace the trait, 
get_info(dfCOI_filtered$Gene.Name)
get_info(df16S_filtered$Gene.Name)
get_info(dfCOI_filtered$Species.Name)
get_info(df16S_filtered$Species.Name)

#perform an alignment and buid a dendrogram for the subset
#use the function clusters which I have created at the top of the code
#clusters(dfCOISubset$COI.Sequence)
#title("Phylogeny of COI gene")

#clusters(df16sSubset$gene16S.Sequence)
#title("Phylogeny of 16S gene")
#randomly select sequences per species for COI and 16S genes to identify whether the two genes (COI and 16S) yield the same or different phylogenetic hypotheses
DaphniaCOI <- dfCOI_filtered %>%
  group_by(Species.Name) %>%
  sample_n(1)
Daphnia16S <- df16S_filtered %>%
  group_by(Species.Name) %>%
  sample_n(1)
#In the previous code the species that overlap within the two gens were selected manually, but it is convenient to use the function merge to get all data in one datafram for further analysis
dfOverlap <- merge(DaphniaCOI, Daphnia16S, by = "Species.Name", all = F)
#check the number of species and species names that can be seen in both datasets (COI and 16S)
get_info(dfOverlap$Species.Name)

#convert the data type
dfOverlap$COI.Sequence <- DNAStringSet(dfOverlap$COI.Sequence)
dfOverlap$gene16S.Sequence <- DNAStringSet(dfOverlap$gene16S.Sequence)
#names function uses to set the names of nucleotides with the species names, which gives species name as the tip labels in phylogenetic tree for downstream analysis
names(dfOverlap$COI.Sequence) <- dfOverlap$Species.Name
names(dfOverlap$gene16S.Sequence) <- dfOverlap$Species.Name
#perform alignment and build dendrogram for the over lapped species to built tanglegram
#Here I'm not using the function clusters as it won't give a dendrogram with assigned variable as I need the assigned cluster variable to perform tanglegram 
DaphniaCOIOverlap.Alignment <- DNAStringSet(muscle::muscle(dfOverlap$COI.Sequence, log = "log.tx", verbose = T), use.names = TRUE)
Daphnia16SOverlap.Alignment <- DNAStringSet(muscle::muscle(dfOverlap$gene16S.Sequence, log = "log.tx", verbose = T), use.names = TRUE)
#check the alignment
lapply(DaphniaCOIOverlap.Alignment, str_count, "-")
lapply(Daphnia16SOverlap.Alignment, str_count, "-")
BrowseSeqs(DaphniaCOIOverlap.Alignment)
BrowseSeqs(Daphnia16SOverlap.Alignment)
#replace the name "Daphnia" with "D." for a clear visualization
names(DaphniaCOIOverlap.Alignment) <- gsub("Daphnia", "D. ", names(DaphniaCOIOverlap.Alignment))
names(Daphnia16SOverlap.Alignment) <- gsub("Daphnia", "D. ", names(Daphnia16SOverlap.Alignment))

dnaBin.DaphniaCOIOverlap <- as.DNAbin(DaphniaCOIOverlap.Alignment)
dnaBin.Daphnia16SOverlap <- as.DNAbin(Daphnia16SOverlap.Alignment)
distanceMatrixDaphniaCOIOverlap <- dist.dna(dnaBin.DaphniaCOIOverlap, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)
distanceMatrixDaphnia16SOverlap <- dist.dna(dnaBin.Daphnia16SOverlap, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)

#COIcluster_NJ1 <- NJ(distanceMatrixDaphniaCOIOverlap)
#gene16Scluster_NJ1 <- NJ(distanceMatrixDaphnia16SOverlap)
#plot(COIcluster_NJ1, main = "Neighbor Joining Phylogenetic tree for COI")
#plot(gene16Scluster_NJ1, main = "Neighbor Joining Phylogenetic tree for 16S")

#tanglegram(COIcluster_NJ1, gene16Scluster_NJ1)
#tanglegram(cluster_NJ1, cluster16s) # Error in ape::as.hclust.phylo(object) : the tree is not ultrametric
#So I use  IdClusters to draw dendrograms
clusters.DaphniaCOIOverlap <- IdClusters(distanceMatrixDaphniaCOIOverlap,
                                         method = "single",
                                         cutoff= 0.02,
                                         showPlot = TRUE,
                                         type = "dendrogram",
                                         verbose = TRUE)


clusters.Daphnia16SOverlap <- IdClusters(distanceMatrixDaphnia16SOverlap,
                                         method = "single",
                                         cutoff= 0.02,
                                         showPlot = TRUE,
                                         type = "dendrogram",
                                         verbose = TRUE)
tanglegram(clusters.DaphniaCOIOverlap, clusters.Daphnia16SOverlap) 
#A dendlist is a function in dendextend package, which produces the dendlist class. This function uses to compare two dendrograms (clusters.DaphniaCOIOverlap, clusters.Daphnia16SOverlap) by chaining them together.
dl <- dendlist(clusters.DaphniaCOIOverlap, clusters.Daphnia16SOverlap)
tanglegram(dl, common_subtrees_color_lines = TRUE, highlight_distinct_edges  = TRUE, highlight_branches_lwd = FALSE, main = "Tanglegram of COI and 16S marker",main_left = "COI gene", main_right = "16S gene", margin_inner=9, cex_main = 1.3)
#all.equal function makes a global comparison of two dendrograms. This will result the difference in branch heights.
all.equal(clusters.DaphniaCOIOverlap, clusters.Daphnia16SOverlap)


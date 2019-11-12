#pART C
#1-Getting the Data from NCBI:
search_result <- entrez_search(db = "pubmed", term = "Daphnia")
# There are 1688 hits but the maximmum I can get are 310!
search.result
COI.search <- entrez_search(db = "nuccore", term = "(Daphnia[ORGN] AND COI[Gene]) NOT (genome[TITL])", retmax = 310)
# When I tried to get retmax = 6066 equel to 16S retmax I got this error message
##Error in entrez_check(response) : 
#HTTP failure 414, the request is too large. For large requests, try using web history as described in the rentrez tutorial. 
#so the maximmum data availble are 310.
#dowloded at 5:16 on wednesday 10/30/2019
# to get fasta format
COI.fetch <- entrez_fetch(db = "nuccore", id = COI.search$ids, rettype = "fasta")
write(COI.fetch, "COI.fetch.fasta", sep = "\n")

string.SetCOI <- readDNAStringSet("COI.fetch.fasta")

#2-Building a dataframe:

dfCOI <- data.frame(COI.Title = names(string.SetCOI), COI.Sequence = paste(string.SetCOI))
dfCOI$Species.Name <- word(dfCOI$COI.Title, 2L, 3L)
dfCOI$Uniqe.Identifier <- word(dfCOI$COI.Title, 1L)
dfCOI$Gene.Name <- str_extract(dfCOI$COI.Title, "COI.*")

# 3- Creating Multiple Sequence Alignment:
#Converting TO DNAStringSet
dfCOI$COI.Sequence <- DNAStringSet(dfCOI$COI.Sequence)
COI.Alignment <- DNAStringSet(muscle::muscle(dfCOI$COI.Sequence, log = "log.tx", verbose = T), use.names = TRUE)
#Check the alignment 
BrowseSeqs(COI.Alignment)
length(COI.Alignment[[1]])  #790 

lapply(COI.Alignment, str_count, "-")
length(unique(dfCOI$Gene.Name))
unique(dfCOI$Gene.Name)
#There are 12 gene after the alignment checking and blastin the outliers, i found that I have 10 samle from Anaompoda sp., I need to clean the data and remove rows from 98:107 which are Anompoda sp.
dfCOISubset <- dfCOI[-c(98:107), ]
dfCOISubset #ROWS were removed from the original data
unique(dfCOISubset$Species.Name) 
length(unique(dfCOISubset$Species.Name)) #37
dfCOISubset$COI.Sequence <- DNAStringSet(dfCOISubset$COI.Sequence)
#MSA FOR the data subset
COI.Alignment1 <- DNAStringSet(muscle::muscle(dfCOISubset$COI.Sequence, gapopen = -400, use.names = TRUE))
BrowseSeqs(COI.Alignment1)
# I will consider COI.Alignment1 for the clustering
#3- cluster by COI marker
dnaBin.COI <- as.DNAbin(COI.Alignment1)

distanceMatrixCOI <- dist.dna(dnaBin.COI, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)


clusters.COI <- IdClusters(distanceMatrixCOI,
                           method = "single",
                           cutoff= 0.02,
                           showPlot = TRUE,
                           type = "clusters",
                           verbose = TRUE)
class(clusters.COI)
clusters.COI
unique(clusters.COI[[1]])
##Data visualization of my clusters  
clusters.COIDend <- IdClusters(distanceMatrixCOI,
                               method = "single",
                               cutoff= 0.02,
                               showPlot = TRUE,
                               type = "dendrogram",
                               verbose = TRUE)
# cluster by OTUs
Data.COI <- cbind(clusters.COI, dfCOISubset)

Data.COISeQ <- Data.COI %>% 
  group_by(cluster) %>% 
  sample_n(1)
COI.Alignment2 <- DNAStringSet(muscle::muscle(Data.COISeQ$COI.Sequence, log = "log.tx", verbose = T), use.names = TRUE)

dnaBin.COISubset <- as.DNAbin(COI.Alignment2)

distanceMatrixCOI2 <- dist.dna(dnaBin.COISubset, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)

cluster_NJ1 <- NJ(distanceMatrixCOI2)

plot(cluster_NJ1, "f", FALSE, cex = 1, main = "Neighbor Joining Phylogenetic tree for COI" )


# To compare between COI & 16S
#COI

COIGene <- Data.COISeQ %>%
  group_by(Species.Name) %>%
  sample_n(1)
COIGene1 <- COIGene[c(2, 7, 12, 13, 14, 15, 22, 23, 26), ]
names(COIGene1$COI.Sequence) <- COIGene1$Species.Name
COI.Alignment7 <- DNAStringSet(muscle::muscle(COIGene1$COI.Sequence, log = "log.tx", verbose = T), use.names = TRUE)

names(COI.Alignment7) <- gsub("Daphnia", "D. ", names(COI.Alignment7))

dnaBin.COIGene <- as.DNAbin(COI.Alignment7)

distanceMatrixCOI7 <- dist.dna(dnaBin.COIGene, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)

COIcluster_NJ1 <- NJ(distanceMatrixCOI7)

plot(COIcluster_NJ1, main = "Neighbor Joining Phylogenetic tree for COI")



# to compare 16S & COI I subset my data by species names and build another phylogenetic tree for 16s
Daphnia16S <- Data.allSeQ %>%
  group_by(Species.Name) %>%
  sample_n(1)
Daphnia.16S <- Daphnia16S[c(1, 7, 12, 14, 16, 17, 24, 29, 32), ]
names(Daphnia.16S$SixteenS.Sequence) <- Daphnia.16S$Species.Name

SixteenS.Alignment9 <- DNAStringSet(muscle::muscle(Daphnia.16S$SixteenS.Sequence, log = "log.tx", verbose = T), use.names = TRUE)
names(SixteenS.Alignment9) <- gsub("Daphnia", "D. ", names(SixteenS.Alignment9))

dnaBin.16Gene <- as.DNAbin(SixteenS.Alignment9)

distanceMatrix16S9 <- dist.dna(dnaBin.16Gene, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)
SixteenScluster_NJ1 <- NJ(distanceMatrix16S9)


plot(SixteenScluster_NJ1, main = "Neighbor Joining Phylogenetic tree for 16S")


tanglegram(SixteenScluster_NJ1, COIcluster_NJ1)

tanglegram(cluster_NJ1, cluster16s)

all.equal.dendrogram(SixteenScluster_NJ1, COIcluster_NJ1)
un fort
tanglegram(cluster_NJ1, cluster16s)


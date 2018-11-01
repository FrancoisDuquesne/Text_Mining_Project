if (!require("easyPubMed")) install.packages("easyPubMed")
library(XML)
library(easyPubMed)
library(ggplot2)

############### PART 1: extraction des textes

##### Option 1: 698 documents via une requete
Querry_String <- "aids"
Ids <- get_pubmed_ids(Querry_String)  
# get_pubmed_ids() decompose les characteristiques de Querry_String en nombre de char et mot clef (ex: AND) et compose une querry comprehensible pour fetch_pubmed_data().
papers <- fetch_pubmed_data(Ids)
# papers est une structure type html qui contien tt les info pour tt les documents (ici +-700 doc) qui ont etes recueillis par fetch_pubmed_data()
# print(papers)

##### Option 2: 52349 documents via importation du fichier xml.
# papers <- xmlParse(file = "/home/francois/Desktop/pubmed18n0924.xml")
# papers est une structure type html qui contien tt les info pour tt les documents (ici +-52300 doc) qui ont etes recueillis par fetch_pubmed_data()


# A partir de papers, on va chercher l'Abstract. 
# /!\ Probleme: les abstracts sont composes de plusieurs parties:
# <AbstractText Label=\"BACKGROUND AND OBJECTIVE\" NlmCategory=\"OBJECTIVE\">
# <AbstractText Label=\"MATERIALS AND METHODS\" NlmCategory=\"METHODS\">
# <AbstractText Label=\"RESULTS\" NlmCategory=\"RESULTS\">
# <AbstractText Label=\"CONCLUSION\" NlmCategory=\"CONCLUSIONS\">
# le code suivant ne fait pas une bonne partition du coup.


Abstract <- unlist(xpathApply(papers, "//AbstractText", saveXML))
# head(Abstract)
# position des "<AbstractText>"
Abstract_pos <- regexpr('<AbstractText>.*<\\/AbstractText>', Abstract)

# on enleve le <AbstractText> au debut de l' Abstract.
Abstract <- substr(Abstract, Abstract_pos + 14, Abstract_pos + attributes(Abstract_pos)$match.length - 16) 
# head(Abstract)

Abstract <- Abstract[Abstract!=""]
# head(Abstract)


############### PART 2: text mining

# Approche Bag of words:
NbrDoc<-100
raw <- Abstract[1:NbrDoc]
# print(raw)

library(quanteda)

# Tokenize
tokens <- tokens(raw, what = "word", 
                       remove_numbers = TRUE, remove_punct = TRUE,
                       remove_symbols = TRUE, remove_hyphens = TRUE)

# for bigrams.
# test.tokens <- tokens_ngrams(test.tokens, n = 1:2)

# minimize capital letters
tokens <- tokens_tolower(tokens)

# stopwords
stop<-stopwords()
new_stopwords<-append(stop,c("fig.","eq.","abstracttext"))
tokens <- tokens_select(tokens, new_stopwords, selection = "remove")

# stem
# tokens <- tokens_wordstem(tokens, language = "english")
# print(tokens)

# Create our first bag-of-words model.
tokens.dfm <- dfm(tokens, tolower = FALSE)

# Transform to a matrix and inspect.
tokens.matrix <- as.matrix(tokens.dfm)
# View(tokens.matrix[1:NbrDoc, 1:100])
# dim(tokens.matrix)

# Tokenfrequence
# In corpus
freq <- sort(colSums(tokens.matrix), decreasing=TRUE)
wf <- data.frame(word=names(freq), freq=freq)

# In specific document
Doc<-5
freqInDoc <- sort(tokens.matrix[Doc,], decreasing=TRUE)
wfindoc <- data.frame(word=names(freqInDoc), freq=freqInDoc)

# plot word frequence
pl <- ggplot(subset(wf, freq > 1) ,aes(word, freq))
# pl <- ggplot(subset(wfindoc, freq > 1) ,aes(word, freq))
pl <- pl + geom_bar(stat="identity", fill="darkred", colour="white")
pl + theme(axis.text.x=element_text(angle=90, hjust=1)) + ggtitle("Uni-Gram Frequency")

# Word Cloud
# library(wordcloud)
# set.seed(100)
# wordcloud(names(freq), freq, min.freq=2, colors=brewer.pal(6, "Dark2"))

# Our function for calculating relative term frequency (TF)
term.frequency <- function(row) {
  row / sum(row)
}

# Our function for calculating inverse document frequency (IDF)
inverse.doc.freq <- function(col) {
  corpus.size <- length(col)
  doc.count <- length(which(col > 0))
  
  log10(corpus.size / doc.count)
}

# Our function for calculating TF-IDF.
tf.idf <- function(x, idf) {
  x * idf
}

# First step, normalize all documents via TF.
tokens.df <- apply(tokens.matrix, 1, term.frequency)
# dim(tokens.df)
# View(tokens.df[1:100, 1:NbrDoc])

# Second step, calculate the IDF vector that we will use - both
tokens.idf <- apply(tokens.matrix, 2, inverse.doc.freq)
str(tokens.idf)

# Lastly, calculate TF-IDF for our training corpus.
tokens.tfidf <-  apply(tokens.df, 2, tf.idf, idf = tokens.idf)
# dim(tokens.tfidf)
# View(tokens.tfidf[1:25, 1:NbrDoc])

# Transpose the matrix
tokens.tfidf <- t(tokens.tfidf)
# dim(tokens.tfidf)
# View(tokens.tfidf[1:NbrDoc, 1:25])

# Check for incopmlete cases.
incomplete.cases <- which(!complete.cases(tokens.tfidf))
raw[incomplete.cases]

# Fix incomplete cases
tokens.tfidf[incomplete.cases,] <- rep(0.0, ncol(tokens.tfidf))
# dim(tokens.tfidf)
sum(which(!complete.cases(tokens.tfidf)))

# Make a clean data frame.
tokens.tfidf.df <- data.frame(tokens.tfidf)
names(tokens.tfidf.df) <- make.names(names(tokens.tfidf.df))

library(irlba)

# Perform SVD. Specifically, reduce dimensionality down to 300 columns
# for our latent semantic analysis (LSA).
irlba <- irlba(t(tokens.tfidf), nv = 30, maxit = 1000)

# Take a look at the new feature data up close.
# View(irlba$v)


# SVD
#-----
# results_SVD<- svd(t(tokens.tfidf))
# eig1<-results_SVD$u[,1]
# eig2<-results_SVD$u[,2]
# eig3<-results_SVD$u[,3]

# Make 3D plot
#----------------
#atribution de nom de lignes
rownames(irlba$v)<-row.names(tokens.tfidf)
eig1<-irlba$v[,1]
eig2<-irlba$v[,2]
eig3<-irlba$v[,3]

# 2D Plot:
# plot(eig1,eig2,col="blue")
# text(eig1,eig2,row.names(irlba$v), cex=0.6, pos=4, col="red")

# 3D Plot:
library("scatterplot3d")

s3d<-scatterplot3d(eig1,eig2,eig3, pch = 16, color="steelblue")
text(s3d$xyz.convert(eig1,eig2,eig3), labels = rownames(irlba$v),
     cex= 0.7, col = "red")

# library(car) # faut aussi installer lib("rgl")
# # 3D plot with the regression plane
# scatter3d(x = eig1, y = eig2, z = eig3)
Querry_String <- 'drug'
index <- match(Querry_String, rownames(tokens.df))
eig_querry <- irlba$u[index,]


s3d<-scatterplot3d(eig1,eig2,eig3, pch = 16, color="steelblue")
text(s3d$xyz.convert(eig1,eig2,eig3), labels = rownames(irlba$v),
     cex= 0.7, col = "green")
s3d$points3d(eig_querry[1],eig_querry[2],eig_querry[3],pch=16,color="red")

# This function computes the euclidean distance between the querry and each document
euc.dist <- function(docs,querry){ 
  dimDocs <- dim(docs)
  squareSum <- 0
  euc.dist <- vector(length=dimDocs[1])
  for (i in (1:dimDocs[1])) {
    for (j in (1:dimDocs[2])){ # TODO: remove the loop and calculate with the whole vectors
      squareSum <- squareSum + (docs[i,j] - querry[j])^2
    }  
    euc.dist[i] <- squareSum ^ 0.5
    squareSum <- 0
  }
  return(euc.dist)
}
# Calculate distance, order and name the rows
distMatrix <- euc.dist(irlba$v, eig_querry)
distDF <- as.data.frame(distMatrix)
rownames(distDF)<-row.names(irlba$v)
distDF <- distDF[order(-distDF$distMatrix), ,drop=FALSE]

# Compute 10 most relevant (tf-idf) words in each documents
bestWords <- function(tokens.tfidf,docId){
  docTfidf <- tokens.tfidf[docId,]
  colnames(docTfidf) <- names(tokens.tfidf)
  docTfidf <- docTfidf[order(-docTfidf),drop=FALSE]
  return(names(docTfidf[1:10]))
}

# matrix with 10 most relevant keywords of the 10 nearests documents 
names <- row.names(distDF)
Result <- matrix(nrow=10,ncol=10)
rownames(Result) <- names[1:10]
for (i in (1:10)) {
  index <- match(names[i],rownames(irlba$v))
  Result[i,] <- bestWords(tokens.tfidf,index)
}
View(Result)
# Clustering
#--------------
res <- kmeans(irlba$v,centers = 2)
plot(irlba$v,col = res$cluster , pch = res$cluster)
points(res$centers, col = 1:5, pch = 8)

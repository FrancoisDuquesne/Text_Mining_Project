library(XML)
library(easyPubMed)
library(ggplot2)

############### PART 1: extraction des textes

##### Option 1: 698 documents via une requete
Querry_String <- "mouse"
Ids <- get_pubmed_ids(Querry_String)  
papers <- fetch_pubmed_data(Ids)

##### Option 2: 52349 documents via importation du fichier xml.
# papers <- xmlParse(file = "/home/francois/Desktop/pubmed18n0924.xml")


# Preprocessing des abstracts
#---------------------------
library(stringr)

# Remove XML text

# Abstract <- unlist(xpathApply(papers, "//Abstract", saveXML))
Article <- unlist(xpathApply(papers, "//PubmedArticle", saveXML))
Abstract <- vector()
for (i in 1:length(Article)) {
  Abstract[i] <- sub(".*<Abstract>\n        ", "", Article[i])
  Abstract[i] <- sub("\n      </Abstract>\n.*", "", Abstract[i])
  # Abstract[i] <- substr(Abstract[i], 14, nchar(Abstract[i]) - 12)
  Abstract[i] <- sub("<CopyrightInformation>.*</CopyrightInformation>", "", Abstract[i])
  if (startsWith(Abstract[i], "<AbstractText>")) {
    Abstract[i] <- substr(Abstract[i], 15, nchar(Abstract[i]) - 18)
  } else {
    Abstract[i] <- str_replace_all(Abstract[i], "<AbstractText.*\">", "")
    Abstract[i] <- gsub("</AbstractText>\n       ","",Abstract[i])
  }
}

# Remove too long or too short Abstracts
for (i in 1:length(Abstract)) {
  if (nchar(Abstract[i])<100 || nchar(Abstract[i])>4000) {
    Abstract[i]<- ""
  }
}
Abstract<-Abstract[Abstract!=""]

#Just to visualyse abstract lengths
a<- vector()
for (i in 1:length(Abstract)) {a[i]<-nchar(Abstract[i])}
plot(a)


############### PART 2: text mining

# Approche Bag of words:
NbrDoc<-100
raw <- Abstract
# print(raw)

library(quanteda)

# Tokenize
tokens <- tokens(raw, what = "word", 
                       remove_numbers = TRUE, remove_punct = TRUE,
                       remove_symbols = TRUE, remove_hyphens = FALSE)

# for bigrams.
# test.tokens <- tokens_ngrams(test.tokens, n = 1:2)

# minimize capital letters
tokens <- tokens_tolower(tokens)

# stopwords
stop<-stopwords()
new_stopwords<-append(stop,c("fig.","eq.","abstracttext","e.g"))
tokens <- tokens_select(tokens, new_stopwords, selection = "remove")
tokens <- tokens_select(tokens,min_nchar = 2, selection ="keep")

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
# pl <- ggplot(subset(wf, freq > 1) ,aes(word, freq))
# # pl <- ggplot(subset(wfindoc, freq > 1) ,aes(word, freq))
# pl <- pl + geom_bar(stat="identity", fill="darkred", colour="white")
# pl + theme(axis.text.x=element_text(angle=90, hjust=1)) + ggtitle("Uni-Gram Frequency")

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
# raw[incomplete.cases]

# Fix incomplete cases
tokens.tfidf[incomplete.cases,] <- rep(0.0, ncol(tokens.tfidf))
# dim(tokens.tfidf)
# sum(which(!complete.cases(tokens.tfidf)))

# Make a clean data frame.
tokens.tfidf.df <- data.frame(tokens.tfidf)
names(tokens.tfidf.df) <- make.names(names(tokens.tfidf.df))

library(irlba)

# Perform SVD. Specifically, reduce dimensionality down to 300 columns
# for our latent semantic analysis (LSA).
irlba <- irlba(t(tokens.tfidf), nv = 10, maxit = 1000)

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
# library("scatterplot3d")

# s3d<-scatterplot3d(eig1,eig2,eig3, pch = 16, color="steelblue")
# text(s3d$xyz.convert(eig1,eig2,eig3), labels = rownames(irlba$v),
#      cex= 0.7, col = "red")

# library(car) # faut aussi installer lib("rgl")
# # 3D plot with the regression plane
# scatter3d(x = eig1, y = eig2, z = eig3)

################ Querries (based on SVD) ###################
#-----------------------------------------------------------

# give a positive querry: as a vector of strings ('querry','querry',...)
posQuerry_String <- c('vaccine','test')
# give a negative querry:
negQuerry_String <- c('eom','cardiac')

posIndex <- vector(length = length(posQuerry_String))
for (i in (1:length(posQuerry_String))) {
  posIndex[i] <- match(posQuerry_String[i], rownames(tokens.df))
}
eig_posQuerry <- irlba$u[posIndex,]

if (negQuerry_String[1] != ""){
  negIndex <- vector(length = length(negQuerry_String))
  for (i in (1:length(negQuerry_String))) {
    negIndex[i] <- match(negQuerry_String[i], rownames(tokens.df))
  }
  eig_negQuerry <- irlba$u[negIndex,]
}

# 
# s3d<-scatterplot3d(eig1,eig2,eig3, pch = 16, color="steelblue")
# text(s3d$xyz.convert(eig1,eig2,eig3), labels = rownames(irlba$v),
#      cex= 0.7, col = "green")
# s3d$points3d(eig_posQuerry[1],eig_posQuerry[2],eig_posQuerry[3],pch=16,color="red")

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
posdistMatrix <- matrix(nrow = length(posQuerry_String),ncol=dim(irlba$v)[1])
posdist <- rep(1,length = dim(irlba$v)[1])
for (i in (1:length(posQuerry_String))) {
  posdistMatrix[i,] <- euc.dist(irlba$v, eig_posQuerry[i,])
  posdist <- posdist * posdistMatrix[i,]
}


if (negQuerry_String[1] != ""){
  negdistMatrix <- matrix(1L,nrow = length(negQuerry_String),ncol=dim(irlba$v)[1])
  negdist <- rep(1,length = dim(irlba$v)[1])
  for (i in (1:length(negQuerry_String))) {
    negdistMatrix[i,] <- euc.dist(irlba$v, eig_negQuerry[i,])
    negdist <- negdist * negdistMatrix[i,]
  }
  distMatrix <- posdist / negdist
}else distMatrix <- posdist


distDF <- as.data.frame(distMatrix)
rownames(distDF)<-row.names(irlba$v)
distDF <- distDF[order(distDF$distMatrix), ,drop=FALSE]



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
# res <- kmeans(irlba$v,centers = 2)
# plot(irlba$v,col = res$cluster , pch = res$cluster)
# points(res$centers, col = 1:5, pch = 8)

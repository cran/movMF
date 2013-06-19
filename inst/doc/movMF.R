### R code from vignette source 'movMF.Rnw'

###################################################
### code chunk number 1: movMF.Rnw:76-89
###################################################
options(width=65, prompt = "R> ", continue = "+  ", useFancyQuotes = FALSE)
cache <- TRUE
library("slam")
library("lattice")
library("movMF")
ltheme <- canonical.theme("pdf", FALSE)
for (i in grep("padding", names(ltheme$layout.heights))) {
  ltheme$layout.heights[[i]] <- 0.2
}
for (i in grep("padding", names(ltheme$layout.widths))) {
  ltheme$layout.widths[[i]] <- 0
}
lattice.options(default.theme = ltheme)


###################################################
### code chunk number 2: corpus-package
###################################################
if (!require("corpus.useR.2008.abstracts", quietly = TRUE)) {
  install.packages("corpus.useR.2008.abstracts", 
                   repos = "http://datacube.wu.ac.at/", type = "source")
}
data("useR_2008_abstracts", package = "corpus.useR.2008.abstracts")


###################################################
### code chunk number 3: movMF.Rnw:559-560
###################################################
cat(paste("R> ", prompt(movMF, filename = NA)$usage[[2]]))


###################################################
### code chunk number 4: movMF.Rnw:1146-1147 (eval = FALSE)
###################################################
## if (!require("corpus.useR.2008.abstracts", quietly = TRUE)) {
##   install.packages("corpus.useR.2008.abstracts", 
##                    repos = "http://datacube.wu.ac.at/", type = "source")
## }
## data("useR_2008_abstracts", package = "corpus.useR.2008.abstracts")


###################################################
### code chunk number 5: movMF.Rnw:1171-1183
###################################################
library("tm")
abstracts_titles <- 
  apply(useR_2008_abstracts[,c("Title", "Abstract")],
        1, paste, collapse = " ")
useR_2008_abstracts_corpus <- Corpus(VectorSource(abstracts_titles))
useR_2008_abstracts_DTM <- 
  DocumentTermMatrix(useR_2008_abstracts_corpus,
                     control = list(
                       tokenize = "MC",
                       stopwords = TRUE,
                       stemming = TRUE,
                       wordLengths = c(3, Inf)))


###################################################
### code chunk number 6: movMF.Rnw:1199-1202
###################################################
library("slam")
ColSums <- col_sums(useR_2008_abstracts_DTM > 0)
sort(ColSums, decreasing = TRUE)[1:10]


###################################################
### code chunk number 7: movMF.Rnw:1210-1213
###################################################
useR_2008_abstracts_DTM <- 
  useR_2008_abstracts_DTM[, ColSums >= 5 & ColSums <= 90]
useR_2008_abstracts_DTM


###################################################
### code chunk number 8: movMF.Rnw:1218-1219
###################################################
useR_2008_abstracts_DTM <- weightTfIdf(useR_2008_abstracts_DTM)


###################################################
### code chunk number 9: fit-movMF (eval = FALSE)
###################################################
## set.seed(2008)
## library("movMF")
## Ks <- c(1:5, 10, 20)
## splits <- sample(rep(1:10, length.out = nrow(useR_2008_abstracts_DTM)))
## useR_2008_movMF <- 
##   lapply(Ks, function(k) 
##          sapply(1:10, function(s) {
##            m <- movMF(useR_2008_abstracts_DTM[splits != s,], 
##                       k = k, nruns = 20)
##            logLik(m, useR_2008_abstracts_DTM[splits == s,])
##          }))
## useR_2008_movMF_common <- 
##   lapply(Ks, function(k) 
##          sapply(1:10, function(s) {
##            m <- movMF(useR_2008_abstracts_DTM[splits != s,], 
##                       k = k, nruns = 20,
##                       kappa = list(common = TRUE))
##            logLik(m, useR_2008_abstracts_DTM[splits == s,])
##          }))


###################################################
### code chunk number 10: movMF.Rnw:1254-1268
###################################################
if(cache & file.exists("movMF.rda")) {
  load("movMF.rda")
  set.seed(2008)
  library("movMF")
  Ks <- c(1:5, 10, 20)
} else {
set.seed(2008)
library("movMF")
Ks <- c(1:5, 10, 20)
splits <- sample(rep(1:10, length.out = nrow(useR_2008_abstracts_DTM)))
useR_2008_movMF <- 
  lapply(Ks, function(k) 
         sapply(1:10, function(s) {
           m <- movMF(useR_2008_abstracts_DTM[splits != s,], 
                      k = k, nruns = 20)
           logLik(m, useR_2008_abstracts_DTM[splits == s,])
         }))
useR_2008_movMF_common <- 
  lapply(Ks, function(k) 
         sapply(1:10, function(s) {
           m <- movMF(useR_2008_abstracts_DTM[splits != s,], 
                      k = k, nruns = 20,
                      kappa = list(common = TRUE))
           logLik(m, useR_2008_abstracts_DTM[splits == s,])
         }))
if(cache) {
  save(useR_2008_movMF, useR_2008_movMF_common,
  file = "movMF.rda")
} else {
  if(file.exists("movMF.rda")) file.remove("movMF.rda")
}
}


###################################################
### code chunk number 11: movMF.Rnw:1281-1293
###################################################
logLiks <- data.frame(logLik = c(unlist(useR_2008_movMF),
                        unlist(useR_2008_movMF_common)),
                      K = c(rep(Ks, sapply(useR_2008_movMF, length)),
                        rep(Ks, sapply(useR_2008_movMF_common, length))),
                      Dataset = seq_len(length(useR_2008_movMF[[1]])),
                      Method = factor(rep(1:2, each = length(unlist(useR_2008_movMF))),
                        1:2, c("free", "common")))
logLiks$logLik <- logLiks$logLik - rep(rep(with(logLiks, tapply(logLik, Dataset, mean)), length(Ks)), 2)
print(xyplot(logLik ~ K | Method, data = logLiks, groups = Dataset, type = "l", lty = 1, 
             xlab = "Number of components", ylab = "Predictive log-likelihood",           
             strip = strip.custom(factor.levels  = 
               expression(paste("Free ", kappa), paste("Common ", kappa)))))


###################################################
### code chunk number 12: movMF.Rnw:1306-1308
###################################################
best_model <- movMF(useR_2008_abstracts_DTM, k = 2, nruns = 20,
                    kappa = list(common = TRUE))


###################################################
### code chunk number 13: movMF.Rnw:1314-1316
###################################################
apply(coef(best_model)$theta, 1, function(x) 
      colnames(coef(best_model)$theta)[order(x, decreasing = TRUE)[1:10]])


###################################################
### code chunk number 14: movMF.Rnw:1329-1336
###################################################
clustering <- predict(best_model)
keywords <- useR_2008_abstracts[,"Keywords"]
keywords <- sapply(keywords, 
                   function(x) sapply(strsplit(x, ", ")[[1]], function(y) 
                                      strsplit(y, "-")[[1]][1]))
tab <- table(Keyword = unlist(keywords), 
             Component = rep(clustering, sapply(keywords, length)))


###################################################
### code chunk number 15: movMF.Rnw:1342-1344
###################################################
tab <- tab[rowSums(tab) > 8,]
tab


###################################################
### code chunk number 16: movMF.Rnw:1357-1361
###################################################
library("vcd")
mosaic(tab, shade = TRUE, 
       labeling_args = list(rot_labels = 0, just_labels = c("center", "right"), 
         pos_varnames = c("left", "right"), rot_varnames = 0))



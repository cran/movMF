\name{useR_2008_abstracts}
\alias{useR_2008_abstracts}
\encoding{UTF-8}
\docType{data}
\title{useR! 2008 Abstracts}
\description{
  Abstracts and some metadata for presentations at the useR! 2008 conference
  held August 12-14, 2008 at Technische Universität Dortmund, Germany.
}
\usage{data("useR_2008_abstracts")}
\format{
  A data frame with 177 observations on the following 5 character
  variables.
  \describe{
    \item{\code{Abstract}:}{the text of the abstract.}
    \item{\code{Title}:}{the title of the presentation.}
    \item{\code{Author}:}{the authors of the presentation, collapsed with
      \samp{ and }.}
    \item{\code{Session}:}{an identifier indicating the session the
      presentation was slotted into (leading \samp{foc} and \samp{kal}
      indicate useR! Focus and useR! Kaleidoscope sections, respectively).}
    \item{\code{Keywords}:}{keywords for the presentation, collapsed with
      \samp{, }.}
  }
}
\details{
  Abstracts were obtained in PDF format from
  \url{https://www.r-project.org/conferences/useR-2008/abstracts/},
  converted to text using \command{pdftotext}, and hand-edited.
  Metadata were provided by the conference organizers.
}
\examples{
data("useR_2008_abstracts")

words <-
    with(useR_2008_abstracts,
         strsplit(Abstract, "[[:space:]]+"))
## (A poor word tokenizer ...)
n_of_words <- sapply(words, length)

authors <-
    with(useR_2008_abstracts,
         strsplit(Author, " and ", fixed = TRUE))
n_of_authors <- sapply(authors, length)

## Do more authors write longer abstracts?
boxplot(n_of_words ~ n_of_authors)

## Session structure:
sessions <-
    with(useR_2008_abstracts,
         sub("-[[:digit:]].*", "", Session))
sort(unique(sessions))    
## Numbers of focus, invited and kaleidoscope presentations:
table(sub("-.*", "", sessions))
}
\keyword{datasets}

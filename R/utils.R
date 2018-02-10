##' @title formatPhrase
##'
##' @param phrase a phrase with taxon name(s)
##' @param taxon a taxon name to be italisized
##' @param ... parameters passed to strsplit
##' @export
##' @author Chenhao Li, Guangchuang Yu
##' @description generate an expression for the phrase with the given taxon italisized
formatPhrase <- function(sentence, taxon, ...){
    sentence.chunks <- strsplit(x=sentence, split=taxon, ...)[[1]]
    return(bquote(.(sentence.chunks[1])~italic(.(taxon) )~.(sentence.chunks[2])))
}



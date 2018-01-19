# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

##' @title parseMetaphlanTSV
##'
##' @param file a single metaphlan table
##' @param delimiter delimiter to separate taxonomic anotations
##' @param index the column number of taxonomic annotation
##' @return a phylo4d object
##' @importClassesFrom  phylobase phylo4d
##' @importFrom magrittr "%>%"
##' @import dplyr
##' @author Chenhao Li
##' @export
##' @description parse a MetaPhlan table to a tree object
parseMetaphlanTSV <- function(file, index=1, delimiter='\\|'){
    taxtab <- read.table(file, sep='\t', stringsAsFactors=FALSE) %>%
        slice(-grep('unclassified', .[,index]))
    tax_chars <- c('k', 'p', 'c', 'o', 'f', 'g', 's', 't')
    tax_split <- strsplit(taxtab$V1, delimiter)    ## split into different taxonomy levels
    child <- sapply(tax_split, tail, n=1)
    tax_class <- do.call(rbind, strsplit(child, '__'))[,1]
    parent <- sapply(tax_split, function(x) ifelse(length(x)>1, x[length(x)-1], 'root'))
    isTip <- !child %in% parent
    index <- c()
    index[isTip] <- 1:sum(isTip)
    index[!isTip] <- (sum(isTip)+1):length(isTip)
    ## tips comes first
    mapping <- data.frame(id=index, row.names=child, isTip, taxaAbun=taxtab$V2)
    edges <- cbind(mapping[parent,]$id, mapping$id)
    edges <- edges[!is.na(edges[,1]),]

    a <- 1
    b <- 1
    nodeSize <- a*log(mapping$taxaAbun[order(mapping$id)]) + b
    nodeClass <- factor(tax_class[order(mapping$id)], levels = rev(tax_chars))

    tr <- phylo4d(x = edges,
                  edge.length = rep(1, nrow(edges)),
                  node.label = rownames(mapping)[!mapping$isTip],
                  tip.label = rownames(mapping[mapping$isTip,]),
                  all.data = data.frame(nodeSize, nodeClass)
    )
    return(tr)
}

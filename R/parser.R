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
##' @param tax.profile a single metaphlan table (either a file path or a a data.frame)
##' @param delimiter delimiter to separate taxonomic anotations
##' @param index the column number of taxonomic annotation
##' @param header whether tax.profile provided as a file path contains a header
##' @param node.size.scale the parameter 'a' controlling node size: nodeSize=a*log(relative_abundance)+b
##' @param node.size.offset the parameter 'b' controlling node size: nodeSize=a*log(relative_abundance)+b
##' @return a treeio::treedata object
##' @importFrom treeio treedata
##' @importFrom magrittr "%>%"
##' @import ape
##' @import dplyr
##' @author Chenhao Li, Guangchuang Yu
##' @export
##' @description parse a MetaPhlan table to a tree object
parseMetaphlanTSV <- function(tax.profile, index=1, header=FALSE, delimiter='\\|', node.size.scale=1, node.size.offset=1){
    if (is.character(tax.profile)) {
        taxtab <- read.table(tax.profile, sep='\t', stringsAsFactors=FALSE, header=header)
    }else{
        taxtab <- tax.profile
    }
    names(taxtab)[index] <- 'tax'
    names(taxtab)[-index] <- 'rel_abun'
    taxtab <- taxtab %>% slice(-grep('unclassified', .[,index]))
    tax_chars <- c('k', 'p', 'c', 'o', 'f', 'g', 's', 't')
    tax_split <- strsplit(taxtab$tax, delimiter)    ## split into different taxonomy levels
    child <- vapply(tax_split, tail, n=1, '')
    tax_class <- do.call(rbind, strsplit(child, '__'))[,1]
    parent <- vapply(tax_split, function(x) ifelse(length(x)>1, x[length(x)-1], 'root'), '')
    isTip <- !child %in% parent
    index <- c()
    index[isTip] <- 1:sum(isTip)
    index[!isTip] <- (sum(isTip)+1):length(isTip)
    ## tips comes first
    mapping <- data.frame(node=index, row.names=child, isTip, taxaAbun=taxtab$rel_abun)
    edges <- cbind(mapping[parent,]$node, mapping$node)
    edges <- edges[!is.na(edges[,1]),]

    a <- node.size.scale
    b <- node.size.offset
    mapping$nodeSize <- a*log(mapping$taxaAbun) + b
    mapping$nodeClass <- factor(tax_class, levels = rev(tax_chars))

    mapping <- mapping[order(mapping$node),]

    node.label <- rownames(mapping)[!mapping$isTip]
    phylo <- structure(list(edge = edges,
                            node.label = node.label,
                            tip.label = rownames(mapping[mapping$isTip,]),
                            edge.length=rep(1, nrow(edges)),
                            Nnode = length(node.label)
                            ),
                       class = "phylo")

    d <- mapping %>% select_(~-isTip)
    treedata(phylo = phylo, data = as_data_frame(d))
}

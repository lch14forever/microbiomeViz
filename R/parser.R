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
    taxtab$tax <- as.character(taxtab$tax)
    taxtab <- taxtab %>% dplyr::slice(grep('unclassified', .[,index], invert=TRUE)) # remove unclassified taxa
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

    d <- mapping %>% dplyr::select_(~-isTip)
    treeio::treedata(phylo = phylo, data = as_data_frame(d))
}

################################################################################

##' @title parsePhyloseq
##'
##' @description This function takes a \code{\link[phyloseq]{phyloseq-class}}
##' and parse it to a \code{\link[treeio]{treedata}} object, that can be further
##' visualized by the \code{\link{tree.backbone}} function. The tax_table slot
##' must not be empty in the phyloseq object.
##'
##' @param physeq A phyloseq object. The tax_table slot must not be empty.
##' @param use_abundance Boolean variable whether to set node sizes according to
##' abundance. If FALSE, all nodes will have the same size.
##' @param node.size.scale the parameter 'a' controlling node size: nodeSize=a*log(relative_abundance)+b
##' @param node.size.offset the parameter 'b' controlling node size: nodeSize=a*log(relative_abundance)+b
##' @return a treeio::treedata object
##' @importFrom treeio treedata
##' @importFrom magrittr "%>%"
##' @import dplyr
##' @author Chenghao Zhu, Chenhao Li, Guangchuang Yu
##' @export
##'
##' @seealso \code{\link{tree.backbone}}, \code{\link{parseMetaphlanTSV}},
##' \code{\link[phyloseq]{phyloseq-class}}
##'
##' @examples
##' data("GlobalPatterns")
##' GP = GlobalPatterns
##' otu_table = otu_table(GP)
##' tax_table = tax_table(GP)
##'
##' # Use the OTUs that make up 99% of the total
##' lib_size = colSums(otu_table)
##' mat = sapply(1:ncol(otu_table), function(i)
##'     otu_table[,i]/lib_size[i] >= 0.01)
##' ind = rowSums(mat)>=1
##'
##' otu_table = otu_table[ind,]
##' tax_table = tax_table[ind,]
##'
##' physeq = phyloseq(otu_table, tax_table) %>%
##'     fix_duplicate_tax()
##'
##' tr = parsePhyloseq(physeq, use_abundance = F)
##' p = tree.backbone(tr, size=1)

parsePhyloseq <- function(physeq,
                          use_abundance = TRUE,
                          node.size.scale = 1,
                          node.size.offset = 1){
    if (!requireNamespace("phyloseq", quietly = TRUE)) {
        stop("Package \"phyloseq\" is needed for this function. Please install it.",
             call. = FALSE)
    }

    # convert taxa_are_rows to TRUE if not
    if(!phyloseq::taxa_are_rows(physeq)){
        otu = phyloseq::otu_table(physeq)
        otu = phyloseq::otu_table(t(otu@.Data), taxa_are_rows = TRUE)
        phyloseq::otu_table(physeq) = otu
    }

    taxtab <- tryCatch(
        phyloseq::tax_table(physeq),
        error = function(e){
            stop("The tax_table is required to draw the cladogram")
        }
    )
    if (use_abundance) {
      otutab <- phyloseq::otu_table(physeq)
    } else {
      row.names = rownames(phyloseq::tax_table(physeq))
      otutab = matrix(rep(0, length(row.names)), ncol = 1)
      rownames(otutab) <- row.names
    }
    # Sometimes the upper level taxonomy is NA, for example:
    # r__Root|k__Bacteria|p__Proteobacteria|c__Alphaproteobacteria|o__Rickettsiales|NA|g__CandidatusPelagibacter
    # Remove all the labels after NA, because they not trustable.
    for(i in 2:ncol(taxtab)){
      na.idx.upper.taxa <- is.na(taxtab[,i])
      na.idx.current.taxa <- is.na(taxtab[, i - 1])
      idx.to.update <- (!na.idx.upper.taxa) & (na.idx.current.taxa)
      if (any(idx.to.update)) {
        taxtab[idx.to.update, i] = NA
      }
    }

    # add taxonomy level if not already have
    if(!grepl("^k__", taxtab[1,1])){
        for(i in 1:ncol(taxtab)){
            tax_level = tolower(strsplit(colnames(taxtab)[i],'')[[1]][1])
            taxtab[,i] = paste0(tax_level, "__", taxtab[,i])
        }
    }

    # summarize taxa
    if(use_abundance){
        otutab_2 = otutab %>%
            rowSums %>%
            data.frame()
        names(otutab_2) = "value"
    }else{
        otutab_2 = data.frame(
            row.names = rownames(otutab),
            value = rep(0, nrow(otutab))
        )
    }
    physeq_2 = phyloseq::phyloseq(
        otu_table(otutab_2, taxa_are_rows = TRUE),
        taxtab
    )

    treetable = data.frame(taxonomy = "r__Root", value = 0)
    if(use_abundance){
        treetable$value = 100
    }

    for(i in 1:ncol(taxtab)){
        summarized = summarize_taxa(physeq_2, level=i)
        summarized = tidytree::filter(summarized, !grepl("NA$", summarized$taxonomy))
        if(use_abundance){
            summarized = mutate(
                summarized,
                value = value/sum(value) * 100
            )
        }
        treetable = rbind(treetable, summarized)
    }
    if(!use_abundance) treetable$value = treetable$value + 5

    parseMetaphlanTSV(treetable,
                      index = 1,
                      header = FALSE,
                      delimiter = "\\|",
                      node.size.scale = node.size.scale,
                      node.size.offset = node.size.offset)
}


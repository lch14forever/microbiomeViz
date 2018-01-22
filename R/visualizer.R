##' @title clade.anno
##'
##' @param gtree a ggtree object
##' @param anno.data a 2 column data.frame of annotation information. It has columns of clade name and color used for highlighting.
##' @return a ggtree object
##' @importFrom treeio treedata
##' @importFrom magrittr "%>%"
##' @import ggtree
##' @import dplyr
##' @author Chenhao Li, Guangchuang Yu
##' @export
##' @description annotate a ggtree plot to highlight certain clades
clade.anno <- function(gtree, anno.data){
    hilight.color <- anno.data$color
    node_list <- anno.data$node
    node_ids <- (gtree$data %>% filter(label %in% node_list ))$node
    for(i in 1:length(node_ids)){
        n <- node_ids[i]
        color <- hilight.color[i]
        mapping <- gtree$data %>% filter(node == n)
        lab <- mapping$label
        angle <- mapping$angle
        angle <- ifelse(angle>180, angle+90, angle-90)
        gtree <-
            gtree + geom_hilight(node=n, fill=color, alpha=0.3, extend=1) +
            geom_cladelabel(angle=angle, node=n, label=lab,
                            offset=1, barsize=0, hjust=0.5)
    }
    gtree
}

##' @title tree.backbone
##'
##' @param tree a treeio::treedata object
##' @param size branch width
##' @param layout tree layout
##' @param shape clade node shape
##' @param fill clade node fill
##' @param color clade node fill
##' @return a ggtree object
##' @importFrom treeio treedata
##' @import ggtree
##' @author Chenhao Li, Guangchuang Yu
##' @export
##' @description basic tree (backbone) plotting utility
tree.backbone <- function(tree, size=2, layout='circular', shape=21, fill='white', color='black'){
    ggtree(tree, size=size, layout = layout)  +
        geom_point(aes(size=I(nodeSize)), shape=shape, fill=fill, color=color)
}

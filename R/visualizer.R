##' @title clade.anno
##'
##' @param gtree a ggtree object
##' @param anno.data a 2 column data.frame of annotation information. It has columns of clade name and color used for highlighting.
##' @param alpha alpha parameter for shading
##' @param anno.depth more specific clades will be shown on the side
##' @param anno.x x position of annotations
##' @param anno.y y position of annotations
##' @return a ggtree object
##' @importFrom magrittr "%>%"
##' @importFrom treeio treedata
##' @import ggtree
##' @import ggplot2
##' @import dplyr
##' @author Chenhao Li, Guangchuang Yu, Chenghao Zhu
##' @export
##' @description annotate a ggtree plot to highlight certain clades
clade.anno <- function(gtree, anno.data, alpha=0.2, anno.depth=3, anno.x=10, anno.y=40){

    short.labs <- letters
    get_offset <- function(x) {(x*0.2+0.2)^2}
    get_angle <- function(node){
        data <- gtree$data
        sp <- tidytree::offspring(data, node)$node
        sp2 <- c(sp, node)
        sp.df <- data[match(sp2, data$node),]
        mean(range(sp.df$angle))
    }
    anno.data <- arrange(anno.data, node)
    hilight.color <- anno.data$color
    node_list <- anno.data$node
    node_ids <- (gtree$data %>% filter(label %in% node_list ) %>% arrange(label))$node
    anno <- rep('white', nrow(gtree$data))
    ## add hilight ... duplicated code
    ## FIXME: duplicated code...
    for(i in 1:length(node_ids)){
        n <- node_ids[i]
        color <- hilight.color[i]
        anno[n] <- color
        mapping <- gtree$data %>% filter(node == n)
        nodeClass <- as.numeric(mapping$nodeClass)
        offset <- get_offset(nodeClass)
        gtree <-
            gtree + geom_hilight(node=n, fill=color, alpha=alpha,
                                 extend=offset)
    }
    gtree$layers <- rev(gtree$layers)
    gtree <- gtree + geom_point2(aes(size=I(nodeSize)), fill=anno, shape=21)
    ## add labels
    short.labs.anno <- NULL
    for(i in 1:length(node_ids)){
        n <- node_ids[i]
        mapping <- gtree$data %>% filter(node == n)
        nodeClass <- as.numeric(mapping$nodeClass)
        if(nodeClass <= anno.depth){## species and strains
            lab <- short.labs[1]
            short.labs <- short.labs[-1]
            # short.labs.anno <- paste0(short.labs.anno, sep='\n', paste0(lab, ': ', mapping$label))
            if(is.null(short.labs.anno)){
                short.labs.anno = data.frame(lab=lab, annot = mapping$label, stringsAsFactors = F)
            }else{
                short.labs.anno = rbind(short.labs.anno,
                                        c(lab, mapping$label))
            }
        }
        else{
            lab <- mapping$label
        }
        offset <- get_offset(nodeClass) - 0.4
        angle <- get_angle(n) + 90
        gtree <- gtree +
            geom_cladelabel(node=n, label=lab, angle=angle,
                            fontsize=1.5+sqrt(nodeClass),
                            offset=offset, barsize=0, hjust=0.5)
    }
    if(is.null(short.labs.anno)){return(gtree)}
    ## add short labels
    anno_shapes = sapply(short.labs.anno$lab, utf8ToInt)
    gtree + geom_point(data = short.labs.anno,
                       aes(x=0, y=0, shape = factor(annot)),
                       size=0, stroke=0) +
        guides(
            shape = guide_legend(override.aes = list(size=3, shape=anno_shapes))
        ) +
        theme(legend.position = c(1.2,0.5),
              legend.title = element_blank())
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

##' @title lefse.lda.plot
##'
##' @param lefse.result the output file from lefse
##' @param negate.class the class to be shown on the negative side (default: everything show on the positive side)
##' @param lda.threshold features with LDA score less than this will be removed
##' @param group a data.frame with two columns for the subgroups of the classes (first column: the names of the features, second column: the subgroup these features belonging to)
##' @import ggplot2
##' @import dplyr
##' @importFrom readr read_tsv
##' @export
##' @author Chenhao Li, Guangchuang Yu
##' @description create a lefse LDA plot
lefse.lda.plot <- function(lefse.result, negate.class=NULL, lda.threshold=NULL, group=NULL){
    input <- read_tsv(lefse.result, col_names = c('feature','dummy','class','lda','pvalue'))

    tmp <- filter(input, !is.na(class)) %>% mutate(lda.scale=lda)
    if(!is.null(lda.threshold)){
        tmp <- filter(tmp, lda>lda.threshold)
    }
    if(!is.null(negate.class) && intersect(negate.class, tmp$class)==negate.class){
        tmp <- mutate(tmp, lda=ifelse(class %in% c(negate.class), lda*-1, lda)) %>%
            mutate(class=factor(class, levels=c(negate.class, setdiff(class, negate.class))))
    }
    tmp$group.origin <- ''
    if(!is.null(group)){
        colnames(group) <- c("X1","group")
        group$X1 <- as.character(group$X1) ## in case imported as factor
        tmp <- merge(tmp, group, by=1, all.x=TRUE)
        tmp$group[is.na(tmp$group)] <- 'Others'
        tmp$group.origin <- tmp$group
    }
    tmp <- mutate(tmp, group=paste0(class, group))

    plot.dat <- dplyr::arrange(tmp, class, lda) %>%
        mutate(group=factor(group, levels=rev(unique(group)), ordered = TRUE)) %>%
        mutate(feature=factor(feature, levels=feature, ordered = TRUE))

    tmp <- select(plot.dat, group, group.origin) %>% unique()
    anno <- tmp$group.origin
    names(anno) <- tmp$group

    p <- ggplot(plot.dat, aes(x=feature, y=lda, fill=class, label=feature)) +
        geom_bar(stat='identity') +
        geom_text(aes(y=0, x=feature), hjust=ifelse(plot.dat$lda<0, -0, 1), nudge_y = -sign(plot.dat$lda)*0.1) +
        coord_flip() +
        labs(x=NULL) +
        theme_classic() +
        facet_grid(group~., scale='free_y', space = "free_y",
                   labeller=labeller(group=anno)) +
        theme(legend.title = element_blank(), legend.key.size = unit(1, 'cm'), legend.position = 'left',
              axis.text.y = element_blank(),
              axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              strip.background = element_blank(),
              strip.text.y = element_blank(),
              panel.spacing = unit(0, "lines"))

    if(!is.null(group)){
        p <- p + theme(strip.text.y = element_text(angle=0, margin = margin(0,3,0,3, "cm")))
    }
    p
}

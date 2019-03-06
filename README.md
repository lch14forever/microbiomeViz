# microbiomeViz -- An R package for visualizing microbiome data

## Authors: 

Chenhao Li, Guangchuang Yu, Chenghao Zhu

## Description:

- parse taxonomic profiling table
- GraPhlAn like tree visualization and annotation
- support cladograme visualization from phyloseq objects and .biom files

## Showcase:

![](http://lchblogs.netlify.com/post/2018-04-20-r-microbiomeviz_example_files/figure-html/unnamed-chunk-5-1.png)

## Example:

### Install:
```{r}
## ## run the following command first if you don't have ggtree installed.
## setRepositories(ind=1:2)
devtools::install_github("lch14forever/microbiomeViz")
```

### Parse a MetaPhlAn table

```{r}
library(microbiomeViz)
data("SRS014459_Stool_profile")
tr <- parseMetaphlanTSV(SRS014459_Stool_profile)
```
### Create a backbone

```{r}
p <- tree.backbone(tr)
p
```
![](https://lichenhao.netlify.com/post/2018-01-18-r-metagenomeViz_files/figure-html/unnamed-chunk-4-1.png)

### Add annotation

```{r}
anno.data <- data.frame(node=c("g__Roseburia", "c__Clostridia", "s__Bacteroides_ovatus"),
                       color='red', stringsAsFactors = FALSE)
p <- clade.anno(p, anno.data)
p
```
![](http://lichenhao.netlify.com/post/2018-01-18-r-metagenomeViz_files/figure-html/unnamed-chunk-5-1.png)

### From phyloseq objects

```{r}
library(phyloseq)
data("GlobalPatterns")
GP = GlobalPatterns

GP = transform_sample_counts(GlobalPatterns, function(otu) otu/sum(otu))
GP = filter_taxa(GP, function(x) max(x)>=0.01,TRUE)
GP = fix_duplicate_tax(GP)

tr = parsePhyloseq(GP)
p = tree.backbone(tr, size=1)
```

### From .biom files

The phyloseq package can import the otu table from .biom files.

```{r}
rich_dense_biom  = system.file("extdata", "rich_dense_otu_table.biom",  package="phyloseq")
rich_dense = import_biom(rich_dense_biom, parseFunction=parse_taxonomy_greengenes)
tr = parsePhyloseq(rich_dense)
p = tree.backbone(tr, size=1)
```

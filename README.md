# microbiomeViz -- An R package for visualizing microbiome data

## Authors: 

Chenhao Li, Guangchuang Yu

## Functions:

- parse taxonomic profiling table
- GraPhlAn like tree visualization and annotation

## TODO list:

- [ ] Add support for biom format

## Example:

### Install:
```{r}
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
![](http://lchblogs.netlify.com/post/2018-01-18-r-metagenomeViz_files/figure-html/unnamed-chunk-4-1.png)

### Add annotation

```{r}
anno.data <- data.frame(node=c("g__Roseburia", "c__Clostridia", "s__Bacteroides_ovatus"),
                       color='red', stringsAsFactors = FALSE)
p <- clade.anno(p, anno.data)
p
```
![](http://lchblogs.netlify.com/post/2018-01-18-r-metagenomeViz_files/figure-html/unnamed-chunk-5-1.png)

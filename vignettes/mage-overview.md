---
title: "Mage Overview"
author: "Charles Bernard"
date: "2018"
output: 
  rmarkdown::html_vignette:
    fig_caption: yes
    keep_md: yes
    number_sections: yes
---



# Prerequisites

The `mage` package has been designed to capture but foremost to classify any significant association 
between the expression of 2 genes in a single cell rna sequencing dataset.

## What does it mean when 2 genes are said to be associated?

The best way to figure yourself out what an association between 2 genes really means would be probably to 
visualize it as a scatterplot with the x and y axis corresponding to their respective expression,
and in which each point would represent a cell. 

As a matter of example, let us produce three scatterplots:

- a linear correlation
- a non-coexistence relationship
- a random relationship


```r
plot(1:10, xlab = "Gene A", ylab = "Gene B", main = "Linear Correlation");
plot(c(0:4, rep(0,5)), c(rep(0,5), 0:4), xlab = "Gene C", ylab = "Gene D", main = "Non-Coexistence");
plot(runif(10, 0, 10), runif(10, 0, 10), xlab = "Gene E", ylab = "Gene F", main = "Random");
```

![](mage-overview_files/figure-html/three_scatterplots-1.png)![](mage-overview_files/figure-html/three_scatterplots-2.png)![](mage-overview_files/figure-html/three_scatterplots-3.png)

Intuitevely, we can quickly see that the two first scatterplots result from a meaningful relationship 
between 2 genes while the third one depict a weaker association. With this respect, the two first gene
pairs must be considered as "significantly associated". 

When 2 genes are said to be "significantly associated", this somehow means that their relationship 
is "remarkable" and clearly distinguishable from noise.

The MINE statistics have been developped first, to estimate the strength of an association between 2 variables
no matter what type of relationship may link them together and second, to characterize any given association
by testing some propreties of the shape of the relationship *(non-linearity, complexity, non-monotonicity ...)*.

The `mage` package rely on these statistics to cluster the associations by type of relationships.

## What is the output of mage?

... to be written

# Tutorial Presentation

This tutorial provides a quick overview of the different features this package can offer. 
It shows how to compute the association strength for each unique pair of genes that 
can be composed from a gene expression matrix, how to characterize them and more essentially 
how to make the most of the variety of the `mage` functions to cluster the significant associations 
based on the type of their relationships *(i.e linear correlation, non-coexistence, sinusoids ...)*.

In a near future, this tutorial will include a series of post-hoc analyses to explain how the standard
output from mage can be used to extract meaningful informations from a biological point of view ...

# Setup

The first step is to load the mage R package as well as the toy dataset that will be used for this tutorial.
The dataset is a gene expression matrix (raw count) used by the Seurat R package as canonical input for
the code examples. It is composed of 80 cells for 230 genes. 


```r
# Load the mage library
library(mage);

# Load the gene expression matrix (raw count) that will be
# used for this tutorial. 
data(pbmc_small_raw_data);
gene_exp_mat <- pbmc_small_raw_data;

# Show all the informations about this dataset:
help(pbmc_small_raw_data);
```

# Compute the set of MINE scores for all unique gene pairs

Our toy dataset contains 230 genes, which result exactly in $230(230 - 1)/2 = 26335$ unique gene pairs


```r
# Notice that the scores are computed in parallel. 
scores <- compute_scores(gene_exp_mat,
                         n.cores = parallel::detectCores() - 1);
#>  * Compute the MINE scores for all pairs of variables ...
#>      this may take quite some time ...
#>  * Compute the Pearson Correlation for all pairs of variables ...
#>  * Compute the Spearman Correlation for all pairs of variables ...

# MEV happens to be redundant with MIC so we can actually get rid of it
scores <- scores[, -"MEV (functionality)"];
```

This is what the table of scores look like

```r
knitr::kable(scores[1:5, ]);
```



X gene    Y gene      MIC (strength)   MIC-p^2 (nonlinearity)   MAS (non-monotonicity)   MCN (complexity)   PEARSON p (linear correlation)   SPEARMAN rho (rank correlation)
--------  ---------  ---------------  -----------------------  -----------------------  -----------------  -------------------------------  --------------------------------
TYROBP    FCER1G           0.8447045                0.2892461                0.0369912           2.000000                        0.7452908                         0.8779817
CST3      LST1             0.8429379                0.6403864                0.1279759           2.584963                        0.4500573                         0.8237421
HLA-DRA   HLA-DRB1         0.8386994               -0.0659080                0.1139512           2.584963                        0.9511085                         0.8718017
LYZ       CST3             0.8344578                0.3241772                0.0896524           2.000000                        0.7143393                         0.8303932
CST3      AIF1             0.8277897                0.7260431                0.1427831           2.584963                        0.3189773                         0.7956390

# input for the Seurat R package examples.


Vignettes are long form documentation commonly included in packages. Because they are part of the distribution of the package, they need to be as compact as possible. The `html_vignette` output type provides a custom style sheet (and tweaks some options) to ensure that the resulting html is as small as possible. The `html_vignette` format:

- Never uses retina figures
- Has a smaller default figure size
- Uses a custom CSS stylesheet instead of the default Twitter Bootstrap style

## Vignette Info

Note the various macros within the `vignette` section of the metadata block above. These are required in order to instruct R how to build the vignette. Note that you should change the `title` field and the `\VignetteIndexEntry` to match the title of your vignette.

## Styles

The `html_vignette` template includes a basic CSS theme. To override this theme you can specify your own CSS in the document metadata as follows:

    output: 
      rmarkdown::html_vignette:
        css: mystyles.css

## Figures

The figure sizes have been customised so that you can easily put two images side-by-side. 


```r
plot(1:10)
plot(10:1)
```

![](mage-overview_files/figure-html/unnamed-chunk-1-1.png)![](mage-overview_files/figure-html/unnamed-chunk-1-2.png)

You can enable figure captions by `fig_caption: yes` in YAML:

    output:
      rmarkdown::html_vignette:
        fig_caption: yes

Then you can use the chunk option `fig.cap = "Your figure caption."` in **knitr**.

## More Examples

You can write math expressions, e.g. $Y = X\beta + \epsilon$, footnotes^[A footnote here.], and tables, e.g. using `knitr::kable()`.


                      mpg   cyl    disp    hp   drat      wt    qsec   vs   am   gear   carb
------------------  -----  ----  ------  ----  -----  ------  ------  ---  ---  -----  -----
Mazda RX4            21.0     6   160.0   110   3.90   2.620   16.46    0    1      4      4
Mazda RX4 Wag        21.0     6   160.0   110   3.90   2.875   17.02    0    1      4      4
Datsun 710           22.8     4   108.0    93   3.85   2.320   18.61    1    1      4      1
Hornet 4 Drive       21.4     6   258.0   110   3.08   3.215   19.44    1    0      3      1
Hornet Sportabout    18.7     8   360.0   175   3.15   3.440   17.02    0    0      3      2
Valiant              18.1     6   225.0   105   2.76   3.460   20.22    1    0      3      1
Duster 360           14.3     8   360.0   245   3.21   3.570   15.84    0    0      3      4
Merc 240D            24.4     4   146.7    62   3.69   3.190   20.00    1    0      4      2
Merc 230             22.8     4   140.8    95   3.92   3.150   22.90    1    0      4      2
Merc 280             19.2     6   167.6   123   3.92   3.440   18.30    1    0      4      4

Also a quote using `>`:

> "He who gives up [code] safety for [code] speed deserves neither."
([via](https://twitter.com/hadleywickham/status/504368538874703872))

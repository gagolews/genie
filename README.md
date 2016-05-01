# Genie (R Package)
## A New, Fast, and Outlier Resistant Hierarchical Clustering Algorithm

[![Build Status](https://travis-ci.org/gagolews/genie.png?branch=master)](https://travis-ci.org/gagolews/genie)

The time needed to apply a hierarchical clustering algorithm
is most often dominated by the number of computations of a pairwise
dissimilarity measure. Such a constraint, for larger data sets,
puts  at a disadvantage the use of all the classical linkage
criteria but the single linkage one. However, it is known that the single
linkage clustering algorithm is very sensitive to outliers, produces highly
skewed dendrograms, and therefore usually does not reflect the true
underlying data structure - unless the clusters are well-separated.
To overcome its limitations, we proposed a new hierarchical clustering linkage
criterion called *Genie*. Namely, our algorithm links two clusters in such 
a way that a chosen economic inequity measure (e.g., the Gini or Bonferroni 
index) of the cluster sizes does not increase drastically above a given 
threshold. Benchmarks indicate a high practical usefulness of the introduced 
method: it most often outperforms the Ward or average linkage in terms of
the clustering quality while retaining the single linkage speed.
The algorithm is easily parallelizable and thus may be run
on multiple threads to speed up its execution further on.
Its memory overhead is small: there is no need to precompute the complete
distance matrix to perform the computations in order to obtain a desired
clustering.

A detailed description of the algorithm will be available in a forthcoming
paper of ours.

**Authors**: [Marek Gagolewski](http://www.gagolewski.com/), 
[Maciej Bartoszuk](http://bartoszuk.rexamine.com), and 
[Anna Cena](http://cena.rexamine.com)

**Homepage**: http://www.gagolewski.com/software/genie/

**CRAN entry**: http://cran.r-project.org/web/packages/genie/

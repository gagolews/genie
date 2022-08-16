# Genie (R Package)

> This project has been superseded by
[genieclust](https://genieclust.gagolewski.com)
(see also: [GitHub](https://github.com/gagolews/genieclust/)),
which features a faster and more feature-rich implementation
of Genie (available for both R and Python).


## A Fast and Robust Hierarchical Clustering Algorithm

The time needed to apply a hierarchical clustering algorithm
is most often dominated by the number of computations of a pairwise
dissimilarity measure. Such a constraint, for larger data sets,
puts the use of all the classical linkage criteria at a disadvantage,
with the exception of the single linkage one. However, it is known that the single
linkage clustering algorithm is very sensitive to outliers, produces highly
skewed dendrograms and therefore usually does not reflect the true
underlying structure of analysed data - unless the clusters are well-separated.
To overcome its limitations, we proposed a hierarchical clustering linkage
criterion called *Genie*. Namely, our algorithm links two clusters in such
a way that the Gini measure of inequity of the cluster sizes
does not exceed a given threshold.
This method most often outperforms the Ward or average linkage in terms of
the clustering quality on benchmark data. At the same time,
Genie retains the high speed of the single linkage approach,
therefore it is also suitable for analysing larger data sets.
The algorithm is easily parallelisable and thus may be run
on multiple threads to speed up its execution further on.
Its memory overhead is small: there is no need to precompute the complete
distance matrix to perform the computations in order to obtain a desired
clustering.

A detailed description of the algorithm can be found in:

Gagolewski M., Bartoszuk M., Cena A., Genie: A new, fast, and outlier-resistant
hierarchical clustering algorithm, *Information Sciences* **363**, 2016, 8–23.
[doi:10.1016/j.ins.2016.05.003](https://dx.doi.org/10.1016/j.ins.2016.05.003).

See also:

Gagolewski M., genieclust: Fast and robust hierarchical clustering,
*SoftwareX* **15**, 2021, 100722.
[doi:10.1016/j.softx.2021.100722](https://dx.doi.org/10.1016/j.softx.2021.100722).


**Authors**: [Marek Gagolewski](https://www.gagolewski.com/),
[Maciej Bartoszuk](http://bartoszuk.rexamine.com), and
[Anna Cena](http://cena.rexamine.com)

**CRAN entry**: <https://cran.r-project.org/web/packages/genie/>

**genieclust**: <https://genieclust.gagolewski.com/>,
<https://github.com/gagolews/genieclust/>,
<https://cran.r-project.org/web/packages/genieclust/>

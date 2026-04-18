# Assess Consistency of Metrics for Simulated TLS Plots

Plots empirical linear charts of density (*N*, trees/ha) and basal area
(*G*, \\{m}^{2}/ha\\) estimates (derived from simulated TLS plots) as a
function of plot size (estimation-size charts) for different plot
designs (circular fixed area, k-tree and angle-count), through
continuous size increments (radius, k and BAF respectively). Size
increments are set at 0.1 m, 1 tree and 0.1 \\{m}^{2}/ha\\ for fixed
area, k-tree and angle-count plot designs, respectively. These
size-estimation line charts represent the consistency in predicting the
stand variables across different values of radius, k and BAF.
Size-estimation charts can be drawn for individual sample plots
(including all plots together in the same charts) or for mean values
(global mean computed for all the sample plots, or for group means if
different strata are considered). Finally, different plot designs can be
compared if specified in the arguments, producing one size-estimation
chart per variable (*N* and *G*).

## Usage

``` r
estimation.plot.size(tree.tls,
                     plot.parameters = data.frame(radius.max = 25,
                                                  k.max = 50,
                                                  BAF.max = 4),
                     dbh.min = 4,
                     average = FALSE, all.plot.designs = FALSE)
```

## Arguments

- tree.tls:

  Data frame with information of trees detected from TLS point cloud
  data in the same format as
  [`tree.detection.single.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.single.scan.md)
  and
  [`tree.detection.multi.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.multi.scan.md)
  ‘Value’.

- plot.parameters:

  Optional data frame containing parameters for circular fixed area,
  k-tree and angle-count plot designs. The parameters are as follows:

&nbsp;

- `radius.max`: maximum plot radius (m) considered for circular fixed
  area plots. If the `radius.max` specified is larger than the farthest
  tree from the plot centre, the horizontal distance from the farthest
  tree will be considered the maximum `radius`. By default, the
  `radius.max` will be 25 m.

- `k.tree.max`: maximum number of trees considered for k-tree plots. If
  `k.tree.max` specified is larger than the maximum number of trees of
  the densest plot, this number of trees will be considered the maximum
  `k.tree.max`. By default, `k.tree.max` is 50.

- `BAF.max`: maximum basal area factor (\\{m}^{2}/ha\\) considered for
  angle-count plots. By default, `BAF.max` is 4.

&nbsp;

- dbh.min:

  Optional minimum *dbh* (cm) considered for detecting trees. By default
  it will be set at 4 cm.

- average:

  Logical; if `TRUE`, plot means values and standard deviation of
  estimations will be represented. By default, it will be set as
  `FALSE`.

- all.plot.designs:

  Logical; if `TRUE`, charts for each plot design are drawn together. By
  default, it will be set as `FALSE`.

## Details

If there are strata in the `tree.tls` argument, they will be
differentiated in charts with different colours. Strata must be
specified in a numeric column named `stratum`.

The `all.plot.designs` argument only works for single strata, and
therefore if there are additional strata in the `tree.tls` argument,
they will be considered equal.

The outputs of this function are inspired by Fig. 3 of Brunner and
Gizachew (2014).

## Value

`Invisible NULL`

## References

Brunner, A., & Gizachew, B. (2014). Rapid detection of stand density,
tree positions, and tree diameter with a 2D terrestrial laser scanner.
*European Journal of Forest Research*, **133(5)**, 819-831.

## Author

Juan Alberto Molina-Valero and Adela Martínez-Calvo.

## Note

Mean values are relevant when plots are representing homogenous strata.

Note that this is an option for choosing the best plot design when field
data are not available. Otherwise, using
[`correlations`](https://molina-valero.github.io/FORTLS/reference/correlations.md),
[`relative.bias`](https://molina-valero.github.io/FORTLS/reference/relative.bias.md)
and
[`optimize.plot.design`](https://molina-valero.github.io/FORTLS/reference/optimize.plot.design.md)
will be more desirable for obtaining the best possible plot design.

## See also

[`tree.detection.single.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.single.scan.md),
[`tree.detection.multi.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.multi.scan.md),
[`tree.detection.several.plots`](https://molina-valero.github.io/FORTLS/reference/tree.detection.several.plots.md)

## Examples

``` r
# \donttest{

# Loading dataset with trees detected from TLS single-scans

data("Rioja.data")

tree.tls <- Rioja.data$tree.tls


# Without strata and plot parameters by default

estimation.plot.size(tree.tls)
#> Computing simulations for plot: '1'
#> Warning: For TLS plot '1', maximum 'k' was reduced to 35, since this is the number of trees in the plot.
#>  (0.38 secs)
#> Computing simulations for plot: '2'
#> Warning: For TLS plot '2', maximum 'k' was reduced to 47, since this is the number of trees in the plot.
#>  (0.4 secs)
#> Computing simulations for plot: '3'
#> Warning: For TLS plot '3', maximum 'k' was reduced to 44, since this is the number of trees in the plot.
#>  (0.38 secs)
#> Computing simulations for plot: '4'
#> Warning: For TLS plot '4', maximum 'k' was reduced to 38, since this is the number of trees in the plot.
#>  (0.39 secs)
#> Computing simulations for plot: '5'
#> Warning: For TLS plot '5', maximum 'k' was reduced to 44, since this is the number of trees in the plot.
#>  (0.4 secs)
#> Computing simulations for plot: '6'
#> Warning: For TLS plot '6', maximum 'k' was reduced to 36, since this is the number of trees in the plot.
#>  (0.37 secs)
#> Computing simulations for plot: '7'
#> Warning: For TLS plot '7', maximum 'k' was reduced to 33, since this is the number of trees in the plot.
#>  (0.35 secs)
#> Computing simulations for plot: '8'
#> Warning: For TLS plot '8', maximum 'k' was reduced to 46, since this is the number of trees in the plot.
#>  (0.39 secs)
#> Computing simulations for plot: '9'
#> Warning: For TLS plot '9', maximum 'k' was reduced to 38, since this is the number of trees in the plot.
#>  (0.35 secs)
#> Computing simulations for plot: '10'
#> Warning: For TLS plot '10', maximum 'k' was reduced to 26, since this is the number of trees in the plot.
#>  (0.32 secs)
#> Computing simulations for plot: '11'
#> Warning: For TLS plot '11', maximum 'k' was reduced to 31, since this is the number of trees in the plot.
#>  (0.37 secs)
#> Computing simulations for plot: '12'
#> Warning: For TLS plot '12', maximum 'k' was reduced to 36, since this is the number of trees in the plot.
#>  (0.38 secs)
#> Computing simulations for plot: '13'
#> Warning: For TLS plot '13', maximum 'k' was reduced to 40, since this is the number of trees in the plot.
#>  (0.38 secs)
#> Computing simulations for plot: '14'
#> Warning: For TLS plot '14', maximum 'k' was reduced to 38, since this is the number of trees in the plot.
#>  (0.39 secs)
#> Computing simulations for plot: '15'
#> Warning: For TLS plot '15', maximum 'k' was reduced to 36, since this is the number of trees in the plot.
#>  (0.39 secs)
#> Computing simulations for plot: '16'
#> Warning: For TLS plot '16', maximum 'k' was reduced to 36, since this is the number of trees in the plot.
#>  (0.39 secs)


estimation.plot.size(tree.tls, average = TRUE)
#> Computing simulations for plot: '1'
#> Warning: For TLS plot '1', maximum 'k' was reduced to 35, since this is the number of trees in the plot.
#>  (0.4 secs)
#> Computing simulations for plot: '2'
#> Warning: For TLS plot '2', maximum 'k' was reduced to 47, since this is the number of trees in the plot.
#>  (0.41 secs)
#> Computing simulations for plot: '3'
#> Warning: For TLS plot '3', maximum 'k' was reduced to 44, since this is the number of trees in the plot.
#>  (0.39 secs)
#> Computing simulations for plot: '4'
#> Warning: For TLS plot '4', maximum 'k' was reduced to 38, since this is the number of trees in the plot.
#>  (0.4 secs)
#> Computing simulations for plot: '5'
#> Warning: For TLS plot '5', maximum 'k' was reduced to 44, since this is the number of trees in the plot.
#>  (0.4 secs)
#> Computing simulations for plot: '6'
#> Warning: For TLS plot '6', maximum 'k' was reduced to 36, since this is the number of trees in the plot.
#>  (0.37 secs)
#> Computing simulations for plot: '7'
#> Warning: For TLS plot '7', maximum 'k' was reduced to 33, since this is the number of trees in the plot.
#>  (0.36 secs)
#> Computing simulations for plot: '8'
#> Warning: For TLS plot '8', maximum 'k' was reduced to 46, since this is the number of trees in the plot.
#>  (0.4 secs)
#> Computing simulations for plot: '9'
#> Warning: For TLS plot '9', maximum 'k' was reduced to 38, since this is the number of trees in the plot.
#>  (0.36 secs)
#> Computing simulations for plot: '10'
#> Warning: For TLS plot '10', maximum 'k' was reduced to 26, since this is the number of trees in the plot.
#>  (0.33 secs)
#> Computing simulations for plot: '11'
#> Warning: For TLS plot '11', maximum 'k' was reduced to 31, since this is the number of trees in the plot.
#>  (0.38 secs)
#> Computing simulations for plot: '12'
#> Warning: For TLS plot '12', maximum 'k' was reduced to 36, since this is the number of trees in the plot.
#>  (0.4 secs)
#> Computing simulations for plot: '13'
#> Warning: For TLS plot '13', maximum 'k' was reduced to 40, since this is the number of trees in the plot.
#>  (0.39 secs)
#> Computing simulations for plot: '14'
#> Warning: For TLS plot '14', maximum 'k' was reduced to 38, since this is the number of trees in the plot.
#>  (0.41 secs)
#> Computing simulations for plot: '15'
#> Warning: For TLS plot '15', maximum 'k' was reduced to 36, since this is the number of trees in the plot.
#>  (0.4 secs)
#> Computing simulations for plot: '16'
#> Warning: For TLS plot '16', maximum 'k' was reduced to 36, since this is the number of trees in the plot.
#>  (0.41 secs)


estimation.plot.size(tree.tls, all.plot.designs = TRUE)
#> Computing simulations for plot: '1'
#> Warning: For TLS plot '1', maximum 'k' was reduced to 35, since this is the number of trees in the plot.
#>  (0.39 secs)
#> Computing simulations for plot: '2'
#> Warning: For TLS plot '2', maximum 'k' was reduced to 47, since this is the number of trees in the plot.
#>  (0.42 secs)
#> Computing simulations for plot: '3'
#> Warning: For TLS plot '3', maximum 'k' was reduced to 44, since this is the number of trees in the plot.
#>  (0.39 secs)
#> Computing simulations for plot: '4'
#> Warning: For TLS plot '4', maximum 'k' was reduced to 38, since this is the number of trees in the plot.
#>  (0.41 secs)
#> Computing simulations for plot: '5'
#> Warning: For TLS plot '5', maximum 'k' was reduced to 44, since this is the number of trees in the plot.
#>  (0.41 secs)
#> Computing simulations for plot: '6'
#> Warning: For TLS plot '6', maximum 'k' was reduced to 36, since this is the number of trees in the plot.
#>  (0.38 secs)
#> Computing simulations for plot: '7'
#> Warning: For TLS plot '7', maximum 'k' was reduced to 33, since this is the number of trees in the plot.
#>  (0.52 secs)
#> Computing simulations for plot: '8'
#> Warning: For TLS plot '8', maximum 'k' was reduced to 46, since this is the number of trees in the plot.
#>  (0.38 secs)
#> Computing simulations for plot: '9'
#> Warning: For TLS plot '9', maximum 'k' was reduced to 38, since this is the number of trees in the plot.
#>  (0.34 secs)
#> Computing simulations for plot: '10'
#> Warning: For TLS plot '10', maximum 'k' was reduced to 26, since this is the number of trees in the plot.
#>  (0.33 secs)
#> Computing simulations for plot: '11'
#> Warning: For TLS plot '11', maximum 'k' was reduced to 31, since this is the number of trees in the plot.
#>  (0.37 secs)
#> Computing simulations for plot: '12'
#> Warning: For TLS plot '12', maximum 'k' was reduced to 36, since this is the number of trees in the plot.
#>  (0.37 secs)
#> Computing simulations for plot: '13'
#> Warning: For TLS plot '13', maximum 'k' was reduced to 40, since this is the number of trees in the plot.
#>  (0.37 secs)
#> Computing simulations for plot: '14'
#> Warning: For TLS plot '14', maximum 'k' was reduced to 38, since this is the number of trees in the plot.
#>  (0.39 secs)
#> Computing simulations for plot: '15'
#> Warning: For TLS plot '15', maximum 'k' was reduced to 36, since this is the number of trees in the plot.
#>  (0.39 secs)
#> Computing simulations for plot: '16'
#> Warning: For TLS plot '16', maximum 'k' was reduced to 36, since this is the number of trees in the plot.
#>  (0.39 secs)



  # }
```

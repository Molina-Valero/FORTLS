# Relative Bias Between Field Estimations and TLS metrics

Computes relative bias between variables estimated from field data and
their TLS counterparts derived from TLS data. Field estimates and TLS
metrics for a common set of plots are required in order to compute
relative bias. These data must come from any of the three different plot
designs currently available (circular fixed area, k-tree and
angle-count) and correspond to plots with incremental values for the
plot design parameter (radius, k and BAF, respectively). In addition to
computing relative bias, interactive line charts graphically
representing the values obtained between each field estimate and its
related TLS metrics are also generated.

## Usage

``` r
relative.bias(simulations,
              variables = c("N", "G", "d", "dg", "d.0", "h", "h.0"),
              save.result = TRUE, dir.result = NULL)
```

## Arguments

- simulations:

  List containing variables estimated from field data and metrics
  derived from TLS data. The structure and format must be analogous to
  output returned by the
  [`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md)
  function. In particular, teh list must include at least one of the
  following named elements:

  - `fixed.area.plot`: same description and format as indicated for same
    named element in `simulations` argument of
    [`correlations`](https://molina-valero.github.io/FORTLS/reference/correlations.md)
    function. The only difference is the columns required when it is
    included in the argument: in adittion to `id`, `radius` and at least
    one of the field estimate columns, it must include at least one TLS
    counterpart for each field estimate condidered.

  - `k.tree.plot`: same description and format as indicated for same
    named element in `simulations` argument of
    [`correlations`](https://molina-valero.github.io/FORTLS/reference/correlations.md)
    function. The only difference is the columns required when it is
    included in the argument: in adittion to `id`, `k` and at least one
    of the field estimate columns, it must include at least one TLS
    counterpart for each field estimate considered.

  - `angle.count.plot`: same description and format as indicated for the
    same named element in the `simulations` argument of
    [`correlations`](https://molina-valero.github.io/FORTLS/reference/correlations.md)
    function. The only difference is the columns that are required when
    this element is included in the argument: in adittion to `id`, `BAF`
    and at least one of the field estimate columns, it must contain at
    least one TLS counterpart for each field estimate considered.

- variables:

  Optional character vector naming field estimates for which the
  relative bias between them and all their available TLS counterparts
  will be computed. If this argument is specified by the user, it must
  contain at least one of the following character strings: “`N`”, “`G`”,
  “`V`”, “`d`”, “`dg`”, “`dgeom`”, “`dharm`”, “`d.0`”, “`dg.0`”,
  “`dgeom.0`”, “`dharm.0`”, “`h`”, “`hg`”, “`hgeom`”, “`hharm`”,
  “`h.0`”, “`hg.0`”,  
  “`hgeom.0`”, or “`hharm.0`”. If this argument is not specified by the
  user, it will be set to
  `c("N", "G", "V", "d", "dg", "d.0", "h", "h.0")` by default. In both
  cases, all the elements in `simulations` argument must include at
  least the columns corresponding to the field estimations specified in
  the `variables` argument.

- save.result:

  Optional logical which indicates whether or not the output files
  described in ‘Output Files’ section must be saved in `dir.result`. If
  this argument is not specified by the user, it will be set to `TRUE`
  by default and, as a consequence, the output files will be saved.

- dir.result:

  Optional character string naming the absolute path of an existing
  directory where files described in ‘Output Files’ section will be
  saved. `.Platform$file.sep` must be used as the path separator in
  `dir.result`. If this argument is not specified by the user, and
  `save.result` is `TRUE`, it will be set to `NULL` by default and, as a
  consequence, the current working directory of the R process will be
  assigned to `dir.result` during the execution.

## Details

For each radius, k or BAF value (according to the currently available
plot designs: circular fixed area, k-tree and angle-count), this
function computes the relative bias between each variable estimated from
field data, and specified in the `variables` argument, and their
counterparts derived from TLS data, and existing in the data frames
included in the `simulations` argument. TLS metrics considered
*counterparts* for each field estimate are detailed below (see
[`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md)
‘Value’ function for details about used notation):

- TLS counterparts for `N` are `N.tls`, `N.hn`, `N.hr`, `N.hn.cov`,
  `N.hr.cov` and `N.sh` in the fixed area and k-tree plot designs; and
  `N.tls` and `N.pam` in the angle-count plot design.

- TLS counterparts for `G` are `G.tls`, `G.hn`, `G.hr`, `G.hn.cov`,
  `G.hr.cov` and `G.sh` in the fixed area and k-tree plot designs; and
  `G.tls` and `G.pam` in the angle-count plot design.

- TLS counterparts for `V` are `V.tls`, `V.hn`, `V.hr`, `V.hn.cov`,
  `V.hr.cov` and `V.sh` in the fixed area and k-tree plot designs; and
  `V.tls` and `V.pam` in the angle-count plot design.

- TLS counterparts for `d`, `dg`, `dgeom`, `dharm`, `d.0`, `dg.0`,
  `dgeom.0`, and `dharm.0` are, respectively: `d.tls`, `dg.tls`,
  `dgeom.tls`, `dharm.tls`, `d.0.tls`, `dg.0.tls`, `dgeom.0.tls`,  
  and `dharm.0.tls` in any of the three available plot designs.

- TLS counterparts for `h`, `hg`, `hgeom`, `hharm`, `h.0`, `hg.0`,
  `hgeom.0`, and `hharm.0` are, respectively `h.tls`, `hg.tls`,
  `hgeom.tls`, `hharm.tls`, `h.0.tls`, `hg.0.tls`, `hgeom.0.tls`, and  
  `hharm.0.tls` in any of the three available plot designs. In adittion,
  `P99` is also taken into account as a counterpart for all these field
  estimates.

The relative bias between a field estimation and any of its TLS
counterparts is estimated as follows \$\$\frac{\frac{1}{n}\sum\_{i =
1}^{n}{y\_{i}} - \frac{1}{n}\sum\_{i =
1}^{n}{x\_{i}}}{\frac{1}{n}\sum\_{i = 1}^{n}{x\_{i}}} \* 100\$\$ where
\\x\_{i}\\ and \\y\_{i}\\ are the values of the field estimate and its
TLS counterpart, respectively, corresponding to plot \\i\\ for \\i = 1,
\ldots, n\\.

## Value

- fixed.area.plot:

  If no “`fixed.area.plot`” element exists in `simulations` argument,
  missing; otherwise, the matrix will include estimates of the the
  relative bias under circular fixed area plot design between each field
  estimation specified in `variables` argument and its TLS
  counterpart(s) existing in the “`fixed.area.plot`” element in the
  `simulations` argument. Each row will correspond to a radius value,
  and the following columns will be included:

  - `radius`: radius (m) of the simulated plots used for computing the
    estimated relative bias.

  - Column(s) ‘`<x>.<y>`’: numeric column(s) containing estimated
    relative bias between ‘`<x>`’, a field estimation, and ‘`<y>`’, a
    TLS counterpart.

- k.tree.plot:

  If no “`k.tree.plot`” element exists in the `simulations` argument,
  missing; otherwise, the matrix will include the relative bias
  estimated in the k-tree plot design between each field estimation
  specified in `variables` argument and its TLS counterpart(s) existing
  in “`k.tree.plot`” element in `simulations` argument. Each row will
  correspond to a k value, and the following columns will be included:

  - `k`: number of trees (trees) of the simulated plots used for
    computing the estimated relative bias.

  - Column(s) ‘`<x>.<y>`’: numeric column(s) containing estimated
    relative bias between ‘`<x>`’, a field estimation, and ‘`<y>`’, a
    TLS counterpart.

- angle.count.plot:

  If no “`angle.count`” element exists in `simulations` argument,
  missing; otherwise, the matrix will contain estimated relative bias
  under angle-count plot design between each field estimation specified
  in the `variables` argument and its TLS counterpart(s) existing in the
  “`angle.count.plot`” element in the `simulations` argument. Each row
  will correspond to a BAF value, and the following columns will be
  included:

  - `BAF`: BAF (\\{m}^{2}/ha\\) of the simulated plots used for
    computing the estimated relative bias.

  - Column(s) ‘`<x>.<y>`’: numeric column(s) containing estimated
    relative bias between ‘`<x>`’, a field estimation, and ‘`<y>`’, a
    TLS counterpart.

## Output Files

During the execution, if the `save.result` argument is `TRUE`, the
function will print the matrices described in the ‘Value’ section to
files. These are written without row names in `dir.result` directory
using [`write.csv`](https://rdrr.io/r/utils/write.table.html) function
from the utils package. The pattern used for naming these files is
`RB.<plot design>.csv`, where `<plot design>` is equal to
“`fixed.area.plot`”, “`k.tree.plot`” or “`angle.count.plot`” is
according to the plot design.

Furthermore, if the `save.result` argument is `TRUE`, interactive line
charts graphically representing relative bias values will also be
created and saved in the `dir.result` directory by means of the
[saveWidget](https://rdrr.io/pkg/htmlwidgets/man/saveWidget.html)
function in the
[htmlwidgets](https://CRAN.R-project.org/package=htmlwidgets) package.
Generated widgets allow users to consult relative bias data directly on
the plots, select/deselect different sets of traces, to zoom and scroll,
and so on. The pattern used for naming these files is
`RB.<x>.<plot design>.html`, where `<plot design>` is indicated for the
previously described files, and `<x>` equals `N`, `G`, `V`, `d` and/or
`h` according to the `variables` argument. All relative biases related
to diameters are plotted in the same chart (files named as
`RB.d.<plot design>.html`), and the same applies to those related to
heights (files named as `RB.h.<plot design>.html`).

## Note

The results obtained using this function are merely descriptive, and
they do not guarantee any type of statistical accuracy in using TLS
metrics instead of field estimations in order to estimate forest
attributes of interest.

## Author

Juan Alberto Molina-Valero and Adela Martínez-Calvo.

## See also

[`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md),
[`correlations`](https://molina-valero.github.io/FORTLS/reference/correlations.md).

## Examples

``` r
# \donttest{

# Load variables estimated from field data, and TLS metrics
# corresponding to Rioja data set

data("Rioja.simulations")


# Establish directory where relative bias results corresponding to Rioja example
# will be saved. For instance, current working directory

dir.result <- tempdir()


# Compute relative bias between field-based estimates of TLS metrics
# corresponding to Rioja example
# Relative bias for variables by default

rb <- relative.bias(simulations = Rioja.simulations, dir.result = dir.result)
#> Computing relative bias for fixed area plots
#>  (0.85 secs)
#> Computing relative bias for k-tree plots
#>  (0.58 secs)
#> Computing relative bias for angle-count plots
#>  (0.47 secs)


# Relative bias for variable 'N'

rb <- relative.bias(simulations = Rioja.simulations, variables = "N",
                    dir.result = dir.result)
#> Computing relative bias for fixed area plots
#>  (0.24 secs)
#> Computing relative bias for k-tree plots
#>  (0.16 secs)
#> Computing relative bias for angle-count plots
#>  (0.11 secs)


# Relative bias corresponding to angle-count design for all available variables

rb <- relative.bias(simulations = Rioja.simulations["angle.count"],
                    variables <- c("N", "G", "V", "d", "dg", "dgeom", "dharm",
                                   "d.0", "dg.0", "dgeom.0", "dharm.0", "h",
                                   "hg", "hgeom", "hharm", "h.0", "hg.0",
                                   "hgeom.0", "hharm.0"),
                    dir.result = dir.result)
#> Computing relative bias for angle-count plots
#>  (0.82 secs)
  # }
```

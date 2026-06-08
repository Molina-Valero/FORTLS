# Correlation Between Field Estimations and TLS Metrics

Computes correlations between variables estimates from field data and
metrics derived from TLS data. Field estimates and TLS metrics for a
common set of plots are required in order to compute correlations. These
data must be obtained from any of the three different plot designs
currently available (fixed area, k-tree and angle-count), and correspond
to plots with incremental values for the plot design parameter (radius,
k and BAF, respectively). Two correlation measures are implemented:
Pearson's correlation coefficient and Spearman's *rho*. In addition to
estimating these measures, tests for association are also executed, and
interactive line charts graphically representing correlations are
generated.

## Usage

``` r
correlations(simulations,
             variables = c("N", "G", "d", "dg", "d.0", "h", "h.0"),
             method = c("pearson", "spearman"), save.result = TRUE,
             dir.result = NULL)
```

## Arguments

- simulations:

  List including estimated variables based on field data and metrics
  derived from TLS data. The structure and format must be analogous to
  output returned by the
  [`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md)
  function. Specifically, it must have at least one of the following
  named elements:

  - `fixed.area.plot`: data frame with field estimates and TLS metrics
    under a circular fixed area plot design. Each row corresponds to a
    (plot, radius) pair, and all or any of the following columns are
    included:

    Plot identification and radius:

    - `id`, `radius`: same description and format as indicated for same
      named columns of `fixed.area.plot` in
      [`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md)
      ‘Value’.

    Variables estimated on the basis of simulated field plots:

    - `N`, `G`, `V`, `d`, `dg`, `dgeom`, `dharm`, `h`, `hg`, `hgeom`,
      `hharm`, `d.0`, `dg.0`, `dgeom.0`, `dharm.0`, `h.0`, `hg.0`,
      `hgeom.0`, `hharm.0`: same description and format as indicated for
      same named columns of `fixed.area.plot` in
      [`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md)
      ‘Value’.

    TLS metrics derived from simulated TLS plots:

    - `N.tls`, `N.hn`, `N.hr`, `N.hn.cov`, `N.hr.cov`, `N.sh`,
      `num.points`,  
      `num.points.est`, `num.points.hom`, `num.points.hom.est`,
      `G.tls`,  
      `G.hn`, `G.hr`, `G.hn.cov`, `G.hr.cov`, `G.sh`, `V.tls`, `V.hn`,
      `V.hr`, `V.hn.cov`,  
      `V.hr.cov`, `V.sh`, `d.tls`, `dg.tls`, `dgeom.tls`, `dharm.tls`,
      `h.tls`,  
      `hg.tls`, `hgeom.tls`, `hharm.tls`, `d.0.tls`, `dg.0.tls`,
      `dgeom.0.tls`, `dharm.0.tls`, `h.0.tls`, `hg.0.tls`,
      `hgeom.0.tls`, `hharm.0.tls`, `P01`, `P05`, `P10`, `P20`, `P25`,
      `P30`, `P40`, `P50`, `P60`, `P70`, `P75`, `P80`, `P90`, `P95`,
      `P99`: same description and format as indicated for same named
      columns of `fixed.area.plot` in
      [`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md)
      ‘Value’.

    If the `fixed.area.plot` element is included in `simulations`
    argument, it must contain at least `id` and `radius` columns, one of
    the field estimates columns and one of the TLS metrics columns.

  - `k.tree.plot`: data frame with field estimates and TLS metrics under
    the k-tree plot design. Each row corresponds to a (plot, k) pair,
    and all or any of the following columns are included:

    Plot identification and k:

    - `id`, `k`: same description and format as indicated for same named
      columns of `k.tree.plot` in
      [`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md)
      ‘Value’.

    Variables estimates on the basis of simulated field plots:

    - `N`, `G`, `V`, `d`, `dg`, `dgeom`, `dharm`, `h`, `hg`, `hgeom`,
      `hharm`, `d.0`, `dg.0`, `dgeom.0`, `dharm.0`, `h.0`, `hg.0`,
      `hgeom.0`, `hharm.0`: same description and format as indicated for
      same named columns of `k.tree.plot` in
      [`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md)
      ‘Value’.

    TLS metrics derived from simulated TLS plots:

    - `N.tls`, `N.hn`, `N.hr`, `N.hn.cov`, `N.hr.cov`, `N.sh`,
      `num.points`,  
      `num.points.est`, `num.points.hom`, `num.points.hom.est`,
      `G.tls`,  
      `G.hn`, `G.hr`, `G.hn.cov`, `G.hr.cov`, `G.sh`, `V.tls`, `V.hn`,
      `V.hr`,  
      `V.hn.cov`, `V.hr.cov`, `V.sh`, `d.tls`, `dg.tls`, `dgeom.tls`,
      `dharm.tls`, `h.tls`, `hg.tls`, `hgeom.tls`, `hharm.tls`,
      `d.0.tls`, `dg.0.tls`,  
      `dgeom.0.tls`, `dharm.0.tls`, `h.0.tls`, `hg.0.tls`,
      `hgeom.0.tls`,  
      `hharm.0.tls`, `P01`, `P05`, `P10`, `P20`, `P25`, `P30`, `P40`,
      `P50`, `P60`, `P70`, `P75`, `P80`, `P90`, `P95`, `P99`: same
      description and format as indicated for same named columns of
      `k.tree.plot` in
      [`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md)
      ‘Value’.

    If a `k.tree.plot` element is included in the `simulations`
    argument, it must include at least `id` and `k` columns, one of the
    field estimate columns, and one of the TLS metrics columns.

  - `angle.count.plot`: data frame with field estimates and TLS metrics
    under the angle-count plot design. Each row corresponds to a (plot,
    BAF) pair, and all or one of the following columns are included:

    Plot identification and BAF:

    - `id`, `BAF`: same description and format as indicated for same
      named columns of `angle.count.plot` in
      [`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md)
      ‘Value’.

    Variables estimated on the basis of simulated field plots:

    - `N`, `G`, `V`, `d`, `dg`, `dgeom`, `dharm`, `h`, `hg`, `hgeom`,
      `hharm`, `d.0`, `dg.0`, `dgeom.0`, `dharm.0`, `h.0`, `hg.0`,
      `hgeom.0`, `hharm.0`: same description and format as indicated for
      same named columns of `angle.count.plot` in
      [`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md)
      ‘Value’.

    TLS metrics derived from simulated TLS plots:

    - `N.tls`, `N.pam`, `num.points`, `num.points.est`,
      `num.points.hom`,  
      `num.points.hom.est`, `G.tls`, `G.pam`, `V.tls`, `V.pam`, `d.tls`,
      `dg.tls`, `dgeom.tls`, `dharm.tls`, `h.tls`, `hg.tls`,
      `hgeom.tls`, `hharm.tls`,  
      `d.0.tls`, `dg.0.tls`, `dgeom.0.tls`, `dharm.0.tls`, `h.0.tls`,  
      `hg.0.tls`, `hgeom.0.tls`, `hharm.0.tls`, `P01`, `P05`, `P10`,
      `P20`, `P25`, `P30`, `P40`, `P50`, `P60`, `P70`, `P75`, `P80`,
      `P90`, `P95`, `P99`: same description and format as indicated for
      same named columns of `angle.count.plot` in
      [`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md)
      ‘Value’.

    If the `angle.count.plot` element is included in the `simulations`
    argument, it must contain at least `id` and `BAF` columns, one of
    the field estimates columns and one of the TLS metrics columns.

- variables:

  Optional character vector naming field estimates for which
  correlations between these and all the available TLS metrics will be
  computed. If this argument is specified by the user, it must include
  at least one of the following character strings: “`N`”, “`G`”, “`V`”,
  “`d`”, “`dg`”, “`dgeom`”, “`dharm`”, “`d.0`”, “`dg.0`”, “`dgeom.0`”,
  “`dharm.0`”, “`h`”, “`hg`”, “`hgeom`”, “`hharm`”, “`h.0`”, “`hg.0`”,
  “`hgeom.0`”, or “`hharm.0`”. If this argument is not specified by the
  user, it will be set to

  `c("N", "G", "V", "d", "dg", "d.0", "h", "h.0")` by default. In both
  cases all the elements in the `simulations` argument must include at
  least the columns corresponding to the field estimates specified in
  the `variables` argument.

- method:

  Optional character vector naming which correlation measurements will
  be used. If this argument is specified by the user, it must include at
  least one of the following character strings: “`pearson`” or
  “`spearman`”. If this argument is not specified by the user, it will
  be set to `c("pearson", "spearman")` by default.

- save.result:

  Optional logical indicating wheter or not the output files described
  in ‘Output Files’ section must be saved in `dir.result` or not. If
  this argument is not specified by the user, it will be set to `TRUE`
  by default and, as a consequence, the output files will be saved.

- dir.result:

  Optional character string naming the absolute path of an existing
  directory where files described in the ‘Output Files’ section will be
  saved. `.Platform$file.sep` must be used as the path separator in
  `dir.result`. If this argument is not specified by the user, and
  `save.result` is `TRUE`, it will be set to `NULL` by default and, as
  consequence, the current working directory of the R process will be
  assigned to `dir.result` during execution.

## Details

For each radius, k or BAF value (according to the currently available
plot designs: circular fixed area, k-tree and angle-count), this
function computes correlations between each variable estimated from
field data specified in the `variables` argument and all the metrics
derived from TLS data existing in the data frames included in the
`simulations` argument.

Two correlation measures are implemented at present: Pearson's
correlation coefficient and Spearman's *rho*. For each method, in
addition to the estimated measure, the p-value of a test for association
is also returned. The
[`cor.test`](https://rdrr.io/r/stats/cor.test.html) function from the
utils package is used to compute both the estimated correlations and the
p-values of the associated tests; more details about these measures and
their tests for association can be found in the corresponding
documentation. There cannot be missing data for three or more plots, and
there cannot be zero standard deviation, in order to prevent missing
correlation values for each (field estimation, TLS metric) pair and plot
design parameter (radius, k or BAF).

Apart from estimated correlations and their corresponding p-values, for
ecah method, the function also returns the plot design parameter and
field estimates, the value of the optimal correlation (i.e. the maximum
of the absolute value of available correlations) and the TLS metric to
which it corresponds.

## Value

- correlations:

  A list including the estimated correlations for each measure specified
  in the `method` argument will be generated. This will include all or
  any of the following named elements:

  - `pearson`: if “`pearson`” is not included in `method` parameter,
    missing; otherwise, the list will include the estimated Pearson's
    correlations for each plot design specified in the `simulations`
    argument. In the latter case, the list will include all or any of
    the following named elements:

    - `fixed.area.plot`: if the `simulations` argument does not have an
      element named “`fixed.area.plot`”, missing; otherwise, the matrix
      will include the estimated Pearson's correlations for a circular
      fixed area plot design. Each row will correspond to a radius
      value, and the following columns will be included:

      - `radius`: radius (m) of the simulated plots used for computing
        the estimated correlations.

      - Column(s) ‘`<x>.<y>`’: numeric column(s) containing estimated
        Pearson's correlations between ‘`<x>`’, a field estimate, and
        ‘`<y>`’, a TLS metric.

    - `k.tree.plot`: if the `simulations` argument does not include an
      element named “`k.tree.plot`”, missing; otherwise, the matrix will
      include the estimated Pearson's correlations for a k-tree plot
      design. Each row will correspond to a k value, and the following
      columns will be included:

      - `k`: number of trees (trees) of the simulated plots used for
        computing the estimated correlations.

      - Column(s) ‘`<x>.<y>`’: same description and format as indicated
        in `correlations$pearson$fixed.area.plot` element.

    - `angle.count.plot`: if the `simulations` argument does not have
      any element named “`angle.count`”, missing; otherwise, the matrix
      will include the estimated Pearson's correlations for the
      angle-count plot design. Each row will correspond to a BAF value,
      and the following columns will be included:

      - `BAF`: BAF (\\{m}^{2}/ha\\) of the simulated plots used for
        computing the estimated correlations.

      - Column(s) ‘`<x>.<y>`’: same description and format as indicated
        in `correlations$pearson$fixed.area.plot` element.

  - `spearman`: if “`spearman`” is not included in `method` parameter,
    missing; otherwise, the list will include the estimated Spearman's
    correlations for each plot design specified in `simulations`
    argument. In the latter case, the structure and format will be
    analogous to that indicated for the previous element but estimated
    Pearson's correlations will be replaced by Spearman's correlations.

- correlations.pval:

  List containing the p-value of the test for association corresponding
  to each measure specified in `method` argument. The structure and
  format will be the same as indicated for the previous element but
  estimated correlations will be replaced by p-values for their
  corresponding tests for association.

- opt.correlations:

  List containing the optimal correlations, and the names of the TLS
  metrics to which they correspond, for each measure specified in
  `method` argument. The list will will include all or any of the
  following named elements:

  - `pearson`: if “`pearson`” is not included in `method` parameter,
    missing; otherwise, the list will include the optimal Pearson's
    correlations, and the names of the TLS metrics to which they
    correspond, for each plot design specified in `simulations`
    argument. In the latter case, it will include all or any of the
    following named elements:

    - `fixed.area.plot`: if `simulations` argument does not have any
      element named “`fixed.area.plot`”, missing; otherwise, the data
      frame will include the optimal Pearson's correlations, and the
      names of the TLS metrics to which they correspond, for a circular
      fixed area plot design. Each row will correspond to a radius
      value, and the following columns will be included:

      - `radius`: radius (m) of the simulated plots used for computing
        the estimated correlations.

      - Columns ‘`<x>.cor`’ and ‘`<x>.metric`’: the former, numeric
        column(s) including optimal Pearson's correlations between
        ‘`<x>`’, a field estimate, and all the available TLS metrics;
        and the latter, character column(s) will include names of the
        TLS metrics to which they correspond.

    - `k.tree.plot`: if the `simulations` argument does not have any
      element named “`k.tree.plot`”, missing; otherwise, the data frame
      will include the optimal Pearson's correlations and the names of
      the TLS metrics to which they correspond for the k-tree plot
      design. Each row will correspond to a k value, and the following
      columns will be included:

      - `k`: number of trees (trees) of the simulated plots used for
        computing the estimated correlations.

      - Columns ‘`<x>.cor`’ and ‘`<x>.metric`’: same description and
        format as indicated in
        `opt.correlations$pearson$fixed.area.plot` element.

    - `angle.count.plot`: if the `simulations` argument does not have
      any element named “`angle.count`”, missing; otherwise, the data
      frame will include the optimal Pearson's correlations, and the
      names of the TLS metrics to which they correspond for the
      angle-count plot design. Each row will correspond to a BAF value,
      and the following columns will be included:

      - `BAF`: BAF (\\{m}^{2}/ha\\) of the simulated plots used for
        computing the estimated correlations.

      - Columns ‘`<x>.cor`’ and ‘`<x>.metric`’: same description and
        format as indicated in
        `opt.correlations$pearson$fixed.area.plot` element.

  - `spearman`: if “`spearman`” is not included in `method` parameter,
    missing; otherwise, the list will include the optimal Spearman's
    correlations, and the names of the TLS metrics to which they
    correspond, for each plot design specified in `simulations`
    argument. In the latter case, the structure and format will be
    analogous to that indicated for the previous element, but optimal
    Pearson's correlations will be replaced by Spearman's correlations.

## Output Files

During the execution, if the `save.result` argument is `TRUE`, the
function prints to files the matrices and data frames included in
`correlations` and `opt.correlations` elements described in ‘Value’.
Both are written without row names in `dir.result` directory by using
the [`write.csv`](https://rdrr.io/r/utils/write.table.html) function in
the utils package. The patterns used for naming these files are
`correlations.<plot design>.<method>.csv` and
`opt.correlations.<plot design>.plot.<method>.csv` for correlation
matrices and optimal correlation data frames, respectively, where
`<plot design>` is equal to “`fixed.area.plot`”, “`k.tree.plot`” or
“`angle.count.plot`” according to plot design, and `<method>` equals
“`pearson`” or “`spearman`” according to the correlation measure.

Furthermore, if the `save.result` argument is `TRUE`, interactive line
charts graphically representing correlations will also be created and
saved in the `dir.result` directory by means of
[saveWidget](https://rdrr.io/pkg/htmlwidgets/man/saveWidget.html)
function in the
[htmlwidgets](https://CRAN.R-project.org/package=htmlwidgets) package.
Generated widgets enable users to consult correlation data directly on
the plots, select/deselect different sets of traces, to zoom and scroll,
etc. The pattern used for naming these files is
`correlations.<x>.<plot design>.<method>.html`, where both
`<plot design>` and `<method>` are as indicated for the previous
described files, and `<x>` is equal to any of elements specified in the
`variables` argument.

## Note

This function is particularly useful for further steps related to
model-based and model-assisted approaches, as correlations measure the
strength of a relationship between two variables (linear for Pearson's
correlation, monotonic for Spearman's correlation).

## Author

Juan Alberto Molina-Valero and Adela Martínez-Calvo.

## See also

[`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md),
[`optimize.plot.design`](https://molina-valero.github.io/FORTLS/reference/optimize.plot.design.md).

[`cor.test`](https://rdrr.io/r/stats/cor.test.html) in utils package.

## Examples

``` r
# \donttest{

# Load field estimates and TLS metrics corresponding to Rioja data set

data("Rioja.simulations")


# Establish directory where correlation results corresponding to the Rioja example
# will be saved. For instance, current working directory

# dir.result <- tempdir()


# Compute correlations between field estimates and TLS metrics corresponding
# to Rioja example
# Pearson's and Spearman's correlations for variables by default

# corr <- correlations(simulations = Rioja.simulations, dir.result = dir.result)


# Pearson's and Spearman's correlations for variable 'N'

# corr <- correlations(simulations = Rioja.simulations, variables = "N",
#                      dir.result = dir.result)


# Only Pearson's correlations for variables by default

# corr <- correlations(simulations = Rioja.simulations, method = "pearson",
#                      dir.result = dir.result)


# Pearson's and Spearman's correlations corresponding to angle-count design for
# all available variables

# corr <- correlations(simulations = Rioja.simulations["angle.count"],
#                      variables <- c("N", "G", "V", "d", "dg", "dgeom", "dharm",
#                                     "d.0", "dg.0", "dgeom.0", "dharm.0", "h",
#                                     "hg", "hgeom", "hharm", "h.0", "hg.0",
#                                     "hgeom.0", "hharm.0"),
#                      dir.result = dir.result)

  # }
```

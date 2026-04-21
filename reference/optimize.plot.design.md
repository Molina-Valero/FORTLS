# Optimize Plot Design Based on Optimal Correlations

Generation of interactive heatmaps graphically represent the optimal
correlations between variables estimated from field data, and metrics
derived from TLS data. These data must be derived from any of the three
different plot designs currently available (circular fixed area, k-tree
and angle-count) and correspond to plots with incremental values for the
plot design parameter (radius, k and BAF, respectively). In addition,
correlation measures that are currently admissible are Pearson's
correlation coefficient and/or Spearman's *rho*.

## Usage

``` r
optimize.plot.design(correlations,
                     variables = c("N", "G", "d", "dg", "d.0", "h", "h.0"),
                     dir.result = NULL)
```

## Arguments

- correlations:

  List including the optimal correlations between field estimations and
  TLS metrics. The structure and format must be analogous to the
  `opt.correlations` element in the output returned by the
  [`correlations`](https://molina-valero.github.io/FORTLS/reference/correlations.md)
  function. In particular, the list must include at least one of the
  following named elements:

  - `pearson`: list containing the optimal Pearson's correlations, and
    the names of the TLS metrics to which they correspond. It must
    include at least one of the following named elements:

    - `fixed.area.plot`: data frame containing the optimal Pearson's
      correlations, and the names of the TLS metrics to which they
      correspond, under circular fixed area plot design. Each row will
      correspond to a radius value, and the data frame will include the
      following columns:

      - `radius`: radius (m) of the simulated plots used for computing
        the estimated correlations.

      - Columns ‘`<x>.cor`’ and ‘`<x>.metric`’: the former, numeric
        column(s) containing optimal Pearson's correlations between
        ‘`<x>`’, a field estimate, and the available TLS metrics when
        [`correlations`](https://molina-valero.github.io/FORTLS/reference/correlations.md)
        function was executed; and the latter, character column(s)
        containing names of the TLS metrics to which they correspond.

      If `fixed.area.plot` is included in `pearson` element, it must
      contain at least the `radius` column and a (‘`<x>.cor`’,
      ‘`<x>.metric`’) pair of columns corresponding to the same field
      estimation.

    - `k.tree.plot`: data frame including the optimal Pearson's
      correlations and the names of the TLS metrics to which they
      correspond, under the k-tree plot design. Each row will correspond
      to a k value, and the following columns will be included:

      - `k`: number of trees (trees) of the simulated plots used for
        computing the estimated correlations.

      - Columns ‘`<x>.cor`’ and ‘`<x>.metric`’: same description and
        format as indicated in `correlations$pearson$fixed.area.plot`
        element.

      If `k.tree.plot` is included in `pearson` element, it must contain
      at least `k` column, and a (‘`<x>.cor`’, ‘`<x>.metric`’) pair of
      columns corresponding to the same field estimation.

    - `angle.count.plot`: data frame including the optimal Pearson's
      correlations and the names of the TLS metrics to which they
      correspond, for the angle-count plot design. Each row will
      correspond to a BAF value, and the data frame will include the
      following columns:

      - `BAF`: BAF (\\{m}^{2}/ha\\) of the simulated plots used to
        compute the estimated correlations.

      - Columns ‘`<x>.cor`’ and ‘`<x>.metric`’: same description and
        format as indicated in `correlations$pearson$fixed.area.plot`
        element.

      If the `angle.count.plot` is included in the `pearson` element, it
      must include at least the `BAF` column and a (‘`<x>.cor`’,
      ‘`<x>.metric`’) pair of columns corresponding to the same field
      estimation.

  - `spearman`: list containing the optimal Spearman's correlations, and
    the names of the TLS metrics to which they correspond. The structure
    and format will be analogous to that indicated for the previous
    element but with optimal Pearson's correlations replaced by
    Spearman's correlations.

- variables:

  Optional character vector naming field estimations whose optimal
  correlations will be represented graphically in the heatmaps generated
  during the execution. If this argument is specified by the user, it
  must include at least one of the following character strings: “`N`”,
  “`G`”, “`V`”, “`d`”, “`dg`”, “`dgeom`”, “`dharm`”, “`d.0`”, “`dg.0`”,
  “`dgeom.0`”, “`dharm.0`”, “`h`”, “`hg`”, “`hgeom`”, “`hharm`”,
  “`h.0`”, “`hg.0`”, “`hgeom.0`”, or “`hharm.0`”. If this argument is
  not specified by the user, it will be set to
  `c("N", "G", "V", "d", "dg", "d.0", "h", "h.0")` by default. In both
  cases, all data frames in the `correlations` argument must have at
  least the (‘`<x>.cor`’, ‘`<x>.metric`’) pairs corresponding to the
  field estimations specified in the `variables` argument.

- dir.result:

  Optional character string naming the absolute path of an existing
  directory where files described in ‘Output Files’ section will be
  saved. `.Platform$file.sep` must be used as the path separator in
  `dir.result`. If this argument is not specified by the user, it will
  be set to `NULL` by default and, as consequence, the current working
  directory of the R process will be assigned to `dir.result` during the
  execution.

## Details

This function represents graphically, by means of interactive heatmaps,
the strongest correlations (positive or negative) for each plot design
and size simulated, between the estimated variables based on field data
specified in the `variables` argument, and metrics derived from TLS
data, under circular fixed area, k-tree and/or angle-count plot designs.

Two correlation measures are implemented at present: Pearson’s
correlation coefficient and Spearman’s rho. Hence, only optimal
correlations based on `correlations` arguments will be taken into
account during the execution.

For each correlation measure and plot design, at least one no missing
value for optimal correlations must be represented; otherwise, execution
will be stopped, and an error message will appear. In addition, at least
two different no missing values for optimal correlations are required to
ensure that the colour palette is correctly applied when the heatmap is
generated.

## Value

Invisible `NULL`.

## Output Files

During the execution, interactive heatmaps graphically representing
optimal correlations values between field estimations and TLS metrics
are created and saved in `dir.result` directory by means of the
[saveWidget](https://rdrr.io/pkg/htmlwidgets/man/saveWidget.html)
function in the
[htmlwidgets](https://CRAN.R-project.org/package=htmlwidgets) package.
The widgets generated allow users to consult optimal correlations values
and TLS metrics to which they correspond directly on the plots, to zoom
and scroll, and so on. The pattern used for naming these files is
`opt.correlations.<plot design>.<method>.html`, where `<plot design>`
equals “`fixed.area.plot`”, “`k.tree.plot`” or “`angle.count.plot`”
according to plot design, and `<method>` equals “`pearson`” or
“`spearman`” according to correlation measure.

## Note

This function is key to choosing the best possible plot design (in terms
of correlation measures) considering all variables of interest before
establishing definitive sampling design.

## Author

Juan Alberto Molina-Valero and Adela Martínez-Calvo.

## See also

[`correlations`](https://molina-valero.github.io/FORTLS/reference/correlations.md).

## Examples

``` r
# \donttest{

# Load field estimations and TLS metrics corresponding to Rioja data set

data("Rioja.simulations")


# Compute correlations between field estimations and TLS metrics corresponding
# to Rioja example, and select optimal correlations results

corr <- correlations(simulations = Rioja.simulations,
                     variables = c("N", "G", "d", "dg", "dgeom","dharm",
                                   "d.0", "dg.0", "dgeom.0", "dharm.0", "h",
                                   "hg", "hgeom", "hharm", "h.0", "hg.0",
                                   "hgeom.0", "hharm.0"),
                     save.result = FALSE)
#> Computing correlations for fixed area plots
#>  (214.39 secs)
#> Computing correlations for k-tree plots
#>  (54.7 secs)
#> Computing correlations for angle-count plots
#>  (42.24 secs)

opt.corr <- corr$opt.correlations


# Establish directory where optimal correlations heatmaps corresponding to Rioja
# example will be saved. For instance, current working directory

dir.result <- tempdir()


# Generate heatmaps for optimal correlations between field estimations and TLS
# metrics corresponding to Rioja example
# Optimal Pearson's and Spearman's correlations for variables by default

# optimize.plot.design(correlations = opt.corr, dir.result = dir.result)


# Optimal Pearson's and Spearman's correlations for variables 'N' and 'G'

optimize.plot.design(correlations = opt.corr, variables = c("N", "G"),
                     dir.result = dir.result)
#> Plotting heatmap(s) for optimal Pearson's correlations 
#>  (0.46 secs)
#> Plotting heatmap(s) for optimal Spearman's correlations 
#>  (0.26 secs)


  # }
```

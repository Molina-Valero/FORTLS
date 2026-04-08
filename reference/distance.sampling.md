# Distance Sampling Methods for Correcting Occlusions Effects

Calculation of the probability of detection of every tree by using
distance sampling methodologies (more specifically point transects
methods), by means of fitting detection functions to the histogram of
tree distribution according to their distance to TLS. Use both half
normal and hazard rate functions, without and with *dbh* as covariate.
These probabilities are used for correcting estimation bias caused by
lack of detection of trees due to occlusion.

## Usage

``` r
distance.sampling(tree.tls,
                  id.plots = NULL,
                  strata.attributes = NULL)
```

## Arguments

- tree.tls:

  Data frame with a list of trees detected and their *dbh* and
  horizontal distances from TLS with the same structure and format as
  [`tree.detection.single.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.single.scan.md)
  and
  [`tree.detection.several.plots`](https://molina-valero.github.io/FORTLS/reference/tree.detection.several.plots.md)
  ‘Value’.

- id.plots:

  Optional vector with plot identification encoded as character string
  or numeric for the plots considered. In this case, `tree.tls` argument
  must include a common column named ‘id’. If this argument is not
  specified by the user, it will be set to NULL by default, and as a
  consequence, all plots will be considered.

- strata.attributes:

  Optional data frame inluding plot radius considered at strata level.
  It must contain a column named ‘stratum’ (numeric) with encoding
  coinciding with that used in previous functions
  ([`normalize`](https://molina-valero.github.io/FORTLS/reference/normalize.md),
  [`tree.detection.single.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.single.scan.md)
  and
  [`tree.detection.several.plots`](https://molina-valero.github.io/FORTLS/reference/tree.detection.several.plots.md))
  for identifying strata. Therefore, strata must heve been included
  previously in ‘tree.tls’. Another column named ‘plot.radius’ (numeric)
  will be required to set maximum horizontal distance (m) considered for
  fitting detection probability functions. If this argument is not
  specified by the user, it will be set to NULL by default, and as a
  consequence, all trees will be included.

## Details

All internal functions related to distance sampling methodologies are
fitted with the [ds](https://rdrr.io/pkg/Distance/man/ds.html) function
included in the [Distance](https://CRAN.R-project.org/package=Distance)
package.

Detection functions are left-truncated at 1 m, according to Astrup et
al., (2014).

Same warning messages as [ds](https://rdrr.io/pkg/Distance/man/ds.html)
function are provided when fits do not converge or another warnings
occur.

For further details on these point transects methods and similar
sampling methodologies, as well as their application with R, see
Buckland et al., (2001); Marques & Buckland, (2003); Miller & Thomas,
(2015) and Clark (2016). Examples of distance sampling analyses, as well
as lectures, are available at
<https://distancesampling.org/resources/vignettes.html>.

## Value

List containing the following elements:

- tree:

  Data frame with detection probabilities for every tree and method.

&nbsp;

- `stratum`: stratum identification (coincident with strata of
  `tree.tls`). If there are not strata, it will be set as a single
  stratum encoded as 1 (numeric).

- `id`: plot identification (coincident with id of `tree.tls`).

- `tree`: tree numbering (coincident with `tree` of `tree.tls`).

- `P.hn`: tree detection probability according to half normal function.

- `P.hn.cov`: tree detection probability according to half normal
  function with *dbh* as covariate.

- `P.hr`: tree detection probability according to half rate function.

- `P.hr.cov`: tree detection probability according to half rate function
  with *dbh* as covariate.

&nbsp;

- parameters:

  Data frame with parameters estimated for detection functions (see
  references for understanding their meaning).

&nbsp;

- `P.hn.scale`: scale parameter for half normal function (sigma).

- `P.hn.cov.scale.intercept`: alpha.0 parameter of scale parameter for
  half normal function with *dbh* as covariate.

- `P.hn.cov.dbh`: alpha.1 parameter of scale parameter for half normal
  function with *dbh* as covariate.

- `P.hr.scale`: scale parameter for half rate function (sigma).

- `P.hr.shape`: shape parameter for half rate function (b).

- `P.hr.cov.scale.intercept`: alpha.0 parameter of scale parameter for
  half normal function with *dbh* as covariate.

- `P.hr.cov.dbh`: alpha.1 parameter of scale parameter for half normal
  function with *dbh* as covariate.

- `P.hr.cov.shape`: shape parameter for half rate function with *dbh* as
  covariate (b).

&nbsp;

- AIC:

  Data frame with Akaike information criterions (AIC) of every detection
  function fit.

&nbsp;

- `P.hn`: AIC of half normal function fit.

- `P.hn.cov`: AIC of half normal function with *dbh* as covariate fit.

- `P.hr`: AIC of half rate function fit.

- `P.hr.cov`: AIC of half rate function with *dbh* as covariate fit.

## Note

Although this step is optional for other functionalities of FORTLS, such
as obtaining metrics and assessing the best plot designs (implemented in
[`metrics.variables`](https://molina-valero.github.io/FORTLS/reference/metrics.variables.md),
[`correlations`](https://molina-valero.github.io/FORTLS/reference/correlations.md),
[`relative.bias`](https://molina-valero.github.io/FORTLS/reference/relative.bias.md)
and
[`optimize.plot.design`](https://molina-valero.github.io/FORTLS/reference/optimize.plot.design.md)),
its inclusion is highly recommended, especially with high rates of
occlusions.

Note that this function could be more useful after assessing the best
possible plot design with
[`estimation.plot.size`](https://molina-valero.github.io/FORTLS/reference/estimation.plot.size.md),
[`correlations`](https://molina-valero.github.io/FORTLS/reference/correlations.md),
[`relative.bias`](https://molina-valero.github.io/FORTLS/reference/relative.bias.md)
or
[`optimize.plot.design`](https://molina-valero.github.io/FORTLS/reference/optimize.plot.design.md)
functions.

## References

Astrup, R., Ducey, M. J., Granhus, A., Ritter, T., & von Lüpke, N.
(2014). Approaches for estimating stand-level volume using terrestrial
laser scanning in a single-scan mode. *Canadian Journal of Forest
Research*, **44(6)**, 666-676.
[doi:10.1139/cjfr-2013-0535](https://doi.org/10.1139/cjfr-2013-0535) .

Buckland, S. T., Anderson, D. R., Burnham, K. P., Laake, J. L.,
Borchers, D. L., & Thomas, L. (2001). *Introduction to distance
sampling: estimating abundance of biological populations*, Oxford,
United Kindown, Oxford University Press.

Clark, R. G. (2016). Statistical efficiency in distance sampling. *PloS
one*, **11(3)**, e0149298.
[doi:10.1371/journal.pone.0149298](https://doi.org/10.1371/journal.pone.0149298)
.

Marques, F. F., & Buckland, S. T. (2003). Incorporating covariates into
standard line transect analyses. *Biometrics*, **59(4)**, 924-935.
[doi:10.1111/j.0006-341X.2003.00107.x](https://doi.org/10.1111/j.0006-341X.2003.00107.x)
.

Miller, D. L., & Thomas, L. (2015). Mixture models for distance sampling
detection functions. *PloS one*, **10(3)**, e0118726.
[doi:10.1371/journal.pone.0118726](https://doi.org/10.1371/journal.pone.0118726)
.

## Author

Juan Alberto Molina-Valero and Adela Martínez-Calvo.

## See also

[`tree.detection.single.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.single.scan.md),
[`tree.detection.several.plots`](https://molina-valero.github.io/FORTLS/reference/tree.detection.several.plots.md),
[`metrics.variables`](https://molina-valero.github.io/FORTLS/reference/metrics.variables.md),
[`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md).

[ds](https://rdrr.io/pkg/Distance/man/ds.html) in
[Distance](https://CRAN.R-project.org/package=Distance) package.

## Examples

``` r
# \donttest{

# Loading example data

data(Rioja.data)

tree.tls <- Rioja.data$tree.tls

# Whithout considering maximum distance

ds <- distance.sampling(tree.tls)


  # }
```

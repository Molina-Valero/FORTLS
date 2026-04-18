# Compute Metrics and Variables for Terrestrial-Based Technologies Point Clouds

This function computes a set of metrics and variables from the
Terrestrial-Based Technologies point cloud data, which have a high
potential to be related or used as direct estimates (in the case of
variables) of forest attributes at plot level. These can be obtained for
different plot designs (circular fixed area, k-tree and angle-count
plots). This function also includes methodologies for correcting
occlusions generated in TLS single-scan point clouds.

## Usage

``` r
metrics.variables(tree.tls, tree.ds = NULL, tree.field = NULL,
                  plot.design = c("fixed.area", "k.tree", "angle.count"),
                  plot.parameters = data.frame(radius = 25, k = 10, BAF = 2),
                  scan.approach = "multi",
                  var.metr = list(tls = NULL, field = NULL),
                  v.calc = "parab", dbh.min = 4, h.min = 1.3,
                  max.dist = Inf, dir.data = NULL,
                  save.result = TRUE, dir.result = NULL)
```

## Arguments

- tree.tls:

  Data frame with information about trees detected from
  terrestrial-based technologies point clouds data in the same format as
  [`tree.detection.single.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.single.scan.md)
  or
  [`tree.detection.multi.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.multi.scan.md)
  ‘Values’.

- tree.ds:

  Optional list containing detection probabilities of trees with
  distance sampling methods. The format must be the same as the ‘Value’
  obtained with
  [`distance.sampling`](https://molina-valero.github.io/FORTLS/reference/distance.sampling.md).
  If this argument is not specified by the user, it will be set to NULL
  by default and, as a consequence, TLS metrics using distance sampling
  based correction will not be calculated for circular fixed area and
  k-tree plot designs.

- tree.field:

  Data frame with information about trees measured in the field plots.
  Each row must correspond to a (plot, tree) pair, and it must include
  at least the following columns:

  - `id`: plot identification encoded as character string or numeric.
    Plot identifications must coincide with those included in `id`
    column of the `tree.tls` argument.

  - `tree`: trees numbering.

  - `h.dist`: horizontal distance (m) from plot centre to tree centre.
    Centres of the field plots must coincide with centres of their
    corresponding Terrestrial-based technologies plots.

  - `dbh`: tree diameter (cm) at breast height (1.3 m).

  - `h`: tree total height (m).

  - `h.blc`: height based live crown (m) (optional).

  - `h.bdc`: height based dead crown (m) (optional).

  - `v.user`: tree volume (m^3) (optional).

  - `w.user`: tree biomass (Mg) (optional).

  - `dead`: integer value indicating for each tree if it is dead (1) or
    not (NA).

- plot.design:

  Vector containing the plot designs considered. By default, all plot
  designs will be considered (circular fixed area, k-tree and
  angle-count plots).

- plot.parameters:

  Data frame containing parameters for circular fixed area, k-tree and
  angle-count plot designs. If there is a `stratum` column in the
  `tree.list.tls` argument, it must have the same number of rows as
  strata values and they must be named using strata encoding. If plot
  parameters are not specified, the corresponding plot designs will not
  be considered in the function. If no parameter is specified, the
  function will stop giving an error message! The parameters are as
  follows:

&nbsp;

- `radius`: plot radius (m) considered for circular fixed area plots.
  Absence of this argument rules out this plot design.

- `k.tree`: number of trees (trees) considered for k-tree plots. Absence
  of this argument rules out this plot design.

- `BAF`: basal area factor (\\{m}^{2}/ha\\) considered for angle-count
  plots. Absence of this argument rules out this plot design.

- `num.trees`: number of dominant trees per ha (tree/ha), i.e. those
  with largest *dbh*, considered for calculating dominant diameters and
  heights. In the absence of this argument, the number will be assumed
  to be 100 trees/ha.

&nbsp;

- scan.approach:

  Character parameter indicating TLS single-scan (‘single’) or TLS
  multi-scan approach or SLAM point clouds (‘multi’) approaches. If this
  argument is not specified by the user, it will be set to ‘multi’
  approach.

- var.metr:

  Optional vector containing all the metrics and variables of interest.
  By default it will be set as NULL and thus, all the metrics and
  variables available will be generated.

- v.calc:

  Optional parameter to calculate volume when is not included in
  tree.tls input data.

- dbh.min:

  Optional minimum *dbh* (cm) considered for detecting trees. By default
  it will be set at 4 cm.

- h.min:

  Optional minimum *h* (m) considered for detecting trees. By default it
  will be set at 1.3 m.

- max.dist:

  Optional argument to specify the maximum horizontal distance
  considered in which trees will be included.

- dir.data:

  Optional character string naming the absolute path of the directory
  where TXT files containing TLS point clouds are located.
  `.Platform$file.sep` must be used as the path separator in `dir.data`,
  and TXT files in the directory must have the same description and
  format as indicated for TXT files in
  [`normalize`](https://molina-valero.github.io/FORTLS/reference/normalize.md)
  ‘Output Files’. If this argument is not specified by the user, it will
  be set to `NULL` by default and, as consequence, the current working
  directory of the R process will be assigned to `dir.dat` during the
  execution.

- save.result:

  Optional logical which indicates whether or not the output files
  described in ‘Output Files’ section must be saved in `dir.result`. If
  this argument is not specified by the user, it will be set to `TRUE`
  by default and, as consequence, the output files will be saved.

- dir.result:

  Optional character string naming the absolute path of an existing
  directory where files described in ‘Output Files’ section will be
  saved. `.Platform$file.sep` must be used as the path separator in
  `dir.result`. If this argument is not specified by the user, and
  `save.result` is `TRUE`, it will be set to `NULL` by default and, as a
  consequence, the current working directory of the R process will be
  assigned to `dir.result` during the execution.

## Details

This function also works for several plots. In this case, a column named
"id" to identify plots (string character or numeric) must be included in
the `tree.list.tls` database argument. This must coincide with the
corresponding "id" assigned in
[`normalize`](https://molina-valero.github.io/FORTLS/reference/normalize.md)
to TXT files saved in `dir.data` (for more details see
[`normalize`](https://molina-valero.github.io/FORTLS/reference/normalize.md)).
In addition, if there are several strata, they can be processed
separately according to `plot.parameters` values (where each row will
represent one stratum). If `tree.list.tls` does not include a specific
"stratum" column, it will be assumed to have only one stratum, which
will be encoded according to `rownames(plot.parameters)[1]`.

Using TLS data, this function computes metrics and estimations of
variables (see ‘Value’) for plot design specified in the
`plot.parameters` argument. The notation used for variables is based on
IUFRO (1959).

At this stage, three plot designs are available:

- Circular fixed area plots, simulated only if a `radius` value is
  specified in the `plot.parameters` argument.

- k-tree plots, simulated only if a `k.tree` value is specified in the
  `plot.parameters` argument.

- Angle-count plots, simulated only if a `BAF` value is specified in the
  `plot.parameters` argument.

Volume is estimated modelling stem profile as a paraboloid and
calculating the volumes of revolution; where trees *dbh* are estimated
in
[`tree.detection.single.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.single.scan.md)
and
[`tree.detection.multi.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.multi.scan.md),
and total heights are estimated as percentile 99 of z coordinate of
points delimited by Voronoi polygons.

Regarding occlusion corrections for estimating variables, apart from
distance sampling methods considered in
[`distance.sampling`](https://molina-valero.github.io/FORTLS/reference/distance.sampling.md),
another occlusion corection based on correcting the shadowing effect
(Seidel & Ammer, 2014) has been included to estimate some variables in
circular fixed area and k-tree plots. In the case of angle-count plots,
occlusion corrections are based on gap probability attenuation with
distance to TLS depending on Poisson model (Strahler et al., 2008;
Lovell et al., 2011).

## Value

List including TLS-based metrics and variables computed for plot designs
considered. The list will contain one element per plot design considered
(fixed.area.plot, k.tree.plot and angle.count.plot):

- fixed.area.plot:

  If no value for `plot.radius` is specified in the `plot.parameters`
  argument, `NULL`; otherwise, the data frame will include TLS metrics
  and variables estimated in a circular fixed area plot design (rows
  represent simulations). The data frame will have the following
  columns:

  Stratum, plot identification and radius:

  - `stratum`: stratum identification encoded as a character string or
    numeric. It must coincide with those included in the `stratum`
    column of `tree.list.tls`.

  - `id`: plot identification encoded as character string or numeric. It
    will coincide with those included in the `id` column of
    `tree.list.tls` and, if applicable, the `distance.sampling`
    argument.

  - `radius`: radius (m) of the simulated plot.

  Terrestrial-Based Technologies variables:

  - `N.tls`: stand density (trees/ha) without occlusion corrections.

  - `N.hn`, `N.hr`, `N.hn.cov`, `N.hr.cov`: stand density (trees/ha)
    with occlusion corrections based on distance sampling methodologies
    (see
    [`distance.sampling`](https://molina-valero.github.io/FORTLS/reference/distance.sampling.md)).
    These columns will not appear if the `distance.sampling` argument is
    `NULL`.

  - `N.sh`: stand density (trees/ha) with correction of the shadowing
    effect. respectively).

  - `G.tls`: stand basal area (\\{m}^{2}/ha\\) without occlusion
    corrections.

  - `G.hn`, `G.hr`, `G.hn.cov`, `G.hr.cov`: stand basal area
    (\\{m}^{2}/ha\\) with occlusion corrections based on distance
    sampling methodologies (see
    [`distance.sampling`](https://molina-valero.github.io/FORTLS/reference/distance.sampling.md)).
    These columns are missing if `distance.sampling` argument is `NULL`.

  - `G.sh`: stand basal area (\\{m}^{2}/ha\\) with correction of the
    shadowing effect.

  - `V.tls`: stand volume (\\{m}^{3}/ha\\) without occlusion
    corrections.

  - `V.hn`, `V.hr`, `V.hn.cov`, `V.hr.cov`: stand volume
    (\\{m}^{3}/ha\\) with occlusion corrections based on distance
    sampling methodologies (see
    [`distance.sampling`](https://molina-valero.github.io/FORTLS/reference/distance.sampling.md)).
    These columns will be missing if the `distance.sampling` argument is
    `NULL`.

  - `V.sh`: stand volume (\\{m}^{3}/ha\\) with correction of the
    shadowing effect.

  - `d.tls`, `dg.tls`, `dgeom.tls`, `dharm.tls`: mean tree
    diameters (cm) at breast height (1.3 m), calculated from the
    arithmetic mean, quadratic mean, geometric mean and harmonic mean,
    respectively.

  - `h.tls`, `hg.tls`, `hgeom.tls`, `hharm.tls`: mean tree heights (m),
    calculated from the arithmetic mean, quadratic mean, geometric mean
    and harmonic mean, respectively.

  - `d.0.tls`, `dg.0.tls`, `dgeom.0.tls`, `dharm.0.tls`: dominant mean
    tree diameters (cm) at breast height (1.3 m), calculated from the
    arithmetic mean, quadratic mean, geometric mean and harmonic mean,
    respectively.

  - `h.0.tls`, `hg.0.tls`, `hgeom.0.tls`, `hharm.0.tls`: dominant mean
    tree heights (m), calculated from the arithmetic mean, quadratic
    mean, geometric mean and harmonic mean, respectively.

  Terrestrial-Based Technologies metrics:

  - `n.pts`, `n.pts.est`, `n.pts.red`, `n.pts.red.est`: number of points
    and estimated number of points (points) belonging to trees normal
    sections (+/- 5 cm) in the original point cloud (`n.pts` and
    `n.pts.est`, respectively); and number of points and estimated
    number of points (points) belonging to trees normal sections (+/-
    5 cm) in the reduced point cloud (`n.pts.red` and `n.pts.red.est`,
    respectively).

  - `P01`, `P05`, `P10`, `P20`, `P25`, `P30`, `P40`, `P50`, `P60`,
    `P70`, `P75`, `P80`, `P90`, `P95`, `P99`: height percentiles derived
    from z coordinates of TLS point clouds relative to ground level (m).

  - `mean.arit.z/rho/phi`, `mean.qua.z/rho/phi`, `mean.geom.z/rho/phi`,
    `mean.harm.z/rho/phi`, `median.z/rho/phi`, `mode.z/rho/phi`: central
    tendency statistics for z, rho (horizontal distance) and r
    (euclidean distance) coordinates (arithmetic, quadratic, geometrical
    and harmonic means, median and mode, respectively).

  - `var.z/rho/phi`, `sd.z/rho/phi`, `CV.z/rho/phi`, `D.z/rho/phi`,
    `ID.z/rho/phi`, `max.z/rho/phi`, `min.z/rho/phi`: dispersion
    tendency statistics for z, rho and r coordinates (variance, standar
    deviation, coefficient of variation, range, interquartile range,
    maximum and minimum, respectively).

  - `kurt.z/rho/phi`, `skw.z/rho/phi`: curtosis and skewness,
    respectively.

  - `L1.z/rho/phi`, `L2.z/rho/phi`, `L3.z/rho/phi`, `L4.z/rho/phi`,
    `L-CV.z/rho/phi`, `MAD-median.z/rho/phi`, `MAD-mode.z/rho/phi`:
    L-moments of order 2, 3 and 4, median of the absolute deviations
    from the overall median and mode of the absolute deviations from the
    overall mode, respectively.

  - `L3.mu.z/rho/phi`, `L4.mu.z/rho/phi`: third and fourth central
    moments, respectively.

  - `PA.2m`, `PA.mean.z`, `PA.mode.z`: percentage of laser returns above
    2 m, mean and mode, respectively.

  - `PB.2m`, `PB.mean.z`, `PB.median.z`, `PB.mode.z`: percentage of
    laser returns below 2 m, mean and mode, respectively.

  - `weibull.b.z/rho/phi`, `weibull.c.z/rho/phi`: scale and shape
    parameters, respectively, for Weibull distribution fitted for z
    coordinates of TLS point clouds relative to ground level.

- k.tree.plot:

  If no value for `k.tree` is specified in the `plot.parameters`
  argument, `NULL`; otherwise, the data frame will include TLS metrics
  and variables estimated in the k-tree plot design (rows represent
  simulations). The data frame will include the following columns (same
  description and format as indicated in `fixed.area.plot` element):

  Stratum, plot identification and k:

  - `stratum`: stratum identification encoded as character string or
    numeric. It must coincide with those included in the `stratum`
    column of `tree.list.tls`.

  - `id`: plot identification encoded as character string or numeric
    value. It will coincide with those included in `id` column of
    `tree.list.tls` or, if applicable, the `distance.sampling` argument.

  - `k`: number of trees in the simulated plot.

  Terrestrial-Based Technologies variables:

  - `N.tls`, `N.hn`, `N.hr`, `N.hn.cov`, `N.hr.cov`, `N.sh`,  
    `G.tls`, `G.hn`, `G.hr`, `G.hn.cov`, `G.hr.cov`, `G.sh`,  
    `V.tls`, `V.hn`, `V.hr`, `V.hn.cov`, `V.hr.cov`, `V.sh`,  
    `d.tls`, `dg.tls`, `dgeom.tls`, `dharm.tls`,  
    `h.tls`, `hg.tls`, `hgeom.tls`, `hharm.tls`,  
    `d.0.tls`, `dg.0.tls`, `dgeom.0.tls`, `dharm.0.tls`,  
    `h.0.tls`, `hg.0.tls`, `hgeom.0.tls`, `hharm.0.tls`

  Terrestrial-Based Technologies metrics:

  - `num.points`, `num.points.est`, `num.points.hom`,
    `num.points.hom.est`,  
    `P01`,`P05`, `P10`, `P20`, `P25`, `P30`, `P40`, `P50`, `P60`, `P70`,
    `P75`, `P80`, `P90`, `P95`, `P99`: same description and format as
    indicated in fixed.area.plot element.

- angle.count.plot:

  If no value for `BAF` is specified in the `plot.parameters` argument,
  `NULL`; otherwise, the data frame will include TLS metrics and
  variables estimated in the angle-count plot design (rows represent
  simulations). The data frame will include the following columns:

  Stratum, plot identification and BAF:

  - `stratum`: stratum identification encoded as character string or
    numeric. It must coincide with those included in `stratum` column of
    `tree.list.tls`.

  - `id`: plot identification encoded as character string or numeric. It
    will coincide with those included in the `id` column of
    `tree.list.tls`.

  - `BAF`: BAF (\\{m}^{2}/ha\\) of the simulated plot.

  Terrestrial-Based Technologies variables:

  - `N.tls`: same description and format as indicated in the
    `fixed.area.plot` element.

  - `N.pam`: stand density (trees/ha) with occlusion correction based on
    gap probability attenuation with distance to TLS. `fixed.area.plot`
    element.

  - `G.tls`: same description and format as indicated for the
    `fixed.area.plot` element.

  - `G.pam`: stand basal area (\\{m}^{2}/ha\\) with occlusion correction
    based on gap probability attenuation with distance from TLS.

  - `V.tls`: same description and format as indicated for the
    `fixed.area.plot` element.

  - `V.pam`: stand volume (\\{m}^{3}/ha\\) with occlusion correction
    based on gap probability attenuation with distance from TLS.

  - `d.tls`, `dg.tls`, `dgeom.tls`, `dharm.tls`,  
    `h.tls`, `hg.tls`, `hgeom.tls`, `hharm.tls`,  
    `d.0.tls`, `dg.0.tls`, `dgeom.0.tls`, `dharm.0.tls`,  
    `h.0.tls`, `hg.0.tls`, `hgeom.0.tls`, `hharm.0.tls`

  Terrestrial-Based Technologies metrics:

  - `num.points`, `num.points.est`, `num.points.hom`,
    `num.points.hom.est`,  
    `P01`, `P05`, `P10`, `P20`, `P25`, `P30`, `P40`, `P50`, `P60`,
    `P70`, `P75`, `P80`, `P90`, `P95`, `P99`: same description and
    format as indicated for the `fixed.area.plot` element.

## Output Files

After computing metrics and variables, if the `save.result` argument is
`TRUE`, the function will save the elements in the list described in
‘Value’ (`fixed.area.plot`, `k.tree.plot` and `angle.count.plot`), which
are different from `NULL` as CSV files. Data frames are written without
row names in the `dir.result` directory by using
[`write.csv`](https://rdrr.io/r/utils/write.table.html) function from
the utils package. The pattern used for naming these files is
‘metrics.variables.\<plot design\>.csv’, being ‘\<plot design\>’ equal
to “fixed.area.plot”, “k.tree.plot” or “angle.count.plot” according to
the plot design.

## Note

In order to optimize plot designs and, therefore, for better use of
`metrics.variables`, other functions such as
[`correlations`](https://molina-valero.github.io/FORTLS/reference/correlations.md),
[`relative.bias`](https://molina-valero.github.io/FORTLS/reference/relative.bias.md)
and
[`estimation.plot.size`](https://molina-valero.github.io/FORTLS/reference/estimation.plot.size.md)
should be used.

This function will be updated as new metrics are developed.

## References

IUFRO (1959). Standarization of symbols in forest mensuration. Vienna,
Austria, IUFRO.

Lovell, J. L., Jupp, D. L. B., Newnham, G. J., & Culvenor, D. S. (2011).
Measuring tree stem diameters using intensity profiles from ground-based
scanning lidar from a fixed viewpoint. *ISPRS Journal of Photogrammetry
and Remote Sensing*, **66(1)**, 46-55.
[doi:10.1016/j.isprsjprs.2010.08.006](https://doi.org/10.1016/j.isprsjprs.2010.08.006)

Seidel, D., & Ammer, C. (2014). Efficient measurements of basal area in
short rotation forests based on terrestrial laser scanning under special
consideration of shadowing. *iForest-Biogeosciences and Forestry*,
**7(4)**, 227.
[doi:10.3832/ifor1084-007](https://doi.org/10.3832/ifor1084-007)

Strahler, A. H., Jupp, D. L. B., Woodcock, C. E., Schaaf, C. B., Yao,
T., Zhao, F., Yang, X., Lovell, J., Culvenor, D., Newnham, G.,
Ni-Miester, W., & Boykin-Morris, W. (2008). Retrieval of forest
structural parameters using a ground-based lidar instrument (Echidna®).
*Canadian Journal of Remote Sensing*, **34(sup2)**, S426-S440.
[doi:10.5589/m08-046](https://doi.org/10.5589/m08-046)

## Author

Juan Alberto Molina-Valero, Adela Martínez-Calvo.

## See also

[`tree.detection.single.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.single.scan.md),
[`tree.detection.multi.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.multi.scan.md),
[`tree.detection.several.plots`](https://molina-valero.github.io/FORTLS/reference/tree.detection.several.plots.md),
[`distance.sampling`](https://molina-valero.github.io/FORTLS/reference/distance.sampling.md),
[`normalize`](https://molina-valero.github.io/FORTLS/reference/normalize.md).

## Examples

``` r
# \donttest{

# Establishment of working directories (optional)
# By default here we propose the current working directory of the R process

dir.data <- tempdir()
dir.result <- tempdir()


# Loading example data included in FORTLS

data("Rioja.data")
tree.tls <- Rioja.data$tree.tls
tree.tls <- tree.tls[tree.tls$id == "1", ]

# Download example of TXT file corresponding to plot 1 from Rioja data set

download.file(url = "https://www.dropbox.com/s/w4fgcyezr2olj9m/Rioja_1.txt?dl=1",
              destfile = file.path(dir.data, "1.txt"), mode = "wb")


# Considering distance sampling methods (only for single-scan approaches)

# ds <- distance.sampling(tree.tls)

met.var.TLS <- metrics.variables(tree.tls = tree.tls,
                                 # tree.ds = ds,
                                 plot.parameters = data.frame(radius = 10, k = 10, BAF = 2),
                                 dir.data = dir.data, dir.result = dir.result)
#> Computing metrics for plot: '1'
#>  (8.63 secs)


  # }
```

# Tree-Level Variables Estimation

Detects trees from point clouds corresponding to TLS multi-scan
approaches and SLAM devices. For each tree detected, the function
calculates the central coordinates and estimates the diameter at 1.3 m
above ground level (which is known as *dbh*, diameter at breast height)
and classifies it as fully visible or partially occluded. Finally, the
function obtains the number of points belonging to normal sections of
trees (those corresponding to *dbh* +/- 5 cm) and estimates them for
both original and reduced (with random selection process) point clouds.

## Usage

``` r
tree.detection.multi.scan(data, single.tree = NULL,
                          dbh.min = 4, dbh.max = 200, h.min = 1.3,
                          geo.dist = 0.1,
                          tls.precision = NULL,
                          stem.section = c(0.7, 3.5), stem.range = NULL, breaks = NULL,
                          slice = 0.1, understory = NULL, bark.roughness = 1,
                          den.type = 1, d.mer = NULL,
                          segmentation = NULL,
                          plot.attributes = NULL, plot = TRUE,
                          threads = 1,
                          dir.data = NULL, save.result = TRUE, dir.result = NULL)
```

## Arguments

- data:

  Data frame with same description and format as indicated for
  [`normalize`](https://molina-valero.github.io/FORTLS/reference/normalize.md)
  ‘Value’.

- single.tree:

  Optional argument to indicate if there is only one tree.

- dbh.min:

  Optional minimum *dbh* (cm) considered for detecting trees. By default
  it will be set at 4 cm.

- dbh.max:

  Optional maximum *dbh* (cm) considered for detecting trees. By default
  it will be set at 200 cm.

- h.min:

  Optional minimum *h* (m) considered for detecting trees. By default it
  will be set at 1.3 m.

- geo.dist:

  Local surface variation (also known as normal change rate, NCR). By
  default it will be set as 0.1. For better understanding of this
  argument see ‘Details’.

- tls.precision:

  Average point cloud precision in cm.

- stem.section:

  Section free of noise (shurb, branches, etc.) considered to detect
  trees. If not specified, an automatic internal algorithm will be
  applied (see ‘Details’).

- stem.range:

  Section considered to estimate straightness tree attributes.

- breaks:

  Height above ground level (m) of slices considered for detecting
  trees. By default it will be considered all possible sections from 0.4
  m to maximum height by 0.3 m intervals (+/- 5 cm).

- slice:

  Slice width considered for detecting trees. By default it will be
  considered as 0.1 m.

- understory:

  Optional argument to indicate if there is dense understory vegetation.

- bark.roughness:

  Bark roughness established in 3 degrees (1 \< 2 \< 3). By default it
  will be considered as 1.

- den.type:

  Numeric argument indicating the dendrometic type used to estimate
  volumen when there are not sections enough to fit a taper equation.
  Dendrometrics types available are the following: cylinder = 0,
  paraboloid = 1 (by default), cone = 2 and neiloid = 3.

- d.mer:

  Top stem diameter (cm) considered to estimate commercial timber
  volume.

- segmentation:

  Tree segmentation.

- plot.attributes:

  Data frame with attributes at plot level. It must contain a column
  named `id` (character string or numeric value) with encoding
  coinciding with that used in `id` argument of
  [`normalize`](https://molina-valero.github.io/FORTLS/reference/normalize.md)
  for identifying plots. If there are strata, another column named
  ‘stratum’ (numeric) will be required for other functionalities of
  FORTLS (see, for instance,
  [`estimation.plot.size`](https://molina-valero.github.io/FORTLS/reference/estimation.plot.size.md)
  or
  [`metrics.variables`](https://molina-valero.github.io/FORTLS/reference/metrics.variables.md)).
  If this argument is not specified by the user, it will be set to NULL
  by default and, as a consequence, the function will not add these
  possible plot attributes.

- plot:

  Optional logical which indicates whether or not the normalized point
  cloud will be plot. If this argument is not specified by the user, it
  will be set to `TRUE` by default and, as consequence, the normalized
  point cloud will be plot.

- threads:

  Number of threads.

- dir.data:

  Optional character string naming the absolute path of the directory
  where LAS files containing TLS point clouds are located.
  `.Platform$file.sep` must be used as the path separator in `dir.data`.
  If this argument is not specified by the user, it will be set to
  `NULL` by default and, as a consequence, the current working directory
  of the R process will be assigned to `dir.data` during the execution.

- save.result:

  Optional logical which indicates whether or not the output files
  described in ‘Output Files’ section should be saved in `dir.result`.
  If this argument is not specified by the user, it will be set to
  `TRUE` by default and, as a consequence, the output files will be
  saved.

- dir.result:

  Optional character string naming the absolute path of an existing
  directory where the files described in ‘Output Files’ section will be
  saved. `.Platform$file.sep` must be used as the path separator in
  `dir.result`. If this argument is not specified by the user, and the
  `save.result` is `TRUE`, it will be set to `NULL` by default and, as a
  consequence, the current working directory of the R process will be
  assigned to `dir.result` during the execution.

## Details

Slices determined by `breaks` argument are clustered using the DBSCAN
algorithm (Ester et al., 1996) on the horizontal plane according to
Cartesian coordinates (x, y). Before and after this process, several
algorithms are used to remove noisy points and apply classification
criteria to select the clusters of trees.

*dbh* is directly estimated for the section of 1.3 m above ground level,
and estimated from other sections using *dbh*~*breaks* linear
regression. Finally, the mean value of all estimates is provided in
‘Value’ as the *dbh* of the tree section.

Volume is estimated modelling stem profile as a paraboloid and
calculating the volumes of revolution; where trees *dbh* are estimated
in
[`tree.detection.single.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.single.scan.md),
and total heights are estimated as percentile 99 of z coordinate of
points delimited by Voronoi polygons.

The number of points corresponding to a normal section (+/- 5 cm) is
estimated in proportion to *dbh*, using the average number of points per
radius unit as reference. In this respect, only tree sections fully
visible at 1.3 m above ground level will be considered for estimating
the average number of points.

Local surface variation (also known as normal change rate ,NCR), is a
quantitative measure of curvature feature (Pauly et al., 2002). This is
useful for distinguishing points belonging to fine branches and foliage
(e.g. leaves, shrubs) and stem points (e.g. Jin et al., 2016; Zhang et
al., 2019). Just as we considered 5 cm as suitable for calculating local
surface variation for the stem separation in forests, according to other
authors (Ma et al., 2015; Xia et al., 2015), we also established the NCR
threshold as 0.1, according to Zhang et al. (2019). However, this
argument (`ncr.threshold`) may be modified in order to use more
appropriate values.

## Value

Data frame with the following columns for every tree detected (each row
corresponds to one tree detected):

- id:

  Optional plot identification encoded as a character string or numeric.
  If this argument is not specified by the user, it will be set to NULL
  by default and, as a consequence, the plot will be encoded as 1.

- file:

  Optional file name identification encoded as character string or
  numeric. If it is null, the file will be encoded as `id` by default.

- tree:

  tree numbering

- Coordinates:

  Cartesian (according to
  <https://en.wikipedia.org/wiki/Cartesian_coordinate_system> notation):

&nbsp;

- `x`: distance on x axis (m) of tree centre.

- `y`: distance on y axis (m) of tree centre.

Azimuthal angles:

- `phi`: angular coordinate (rad) of tree centre.

&nbsp;

- h.dist:

  horizontal distance (m) from plot centre to tree centre.

- *dbh*:

  estimated tree diameter (cm) at breast height (1.3 m).

- *h*:

  estimated tree total height (m).

- *h.com*:

  estimated commercial tree height (m) according to the top diameter
  defined in the argument `d.mer`.

- *v*:

  estimated tree stem volume (m^3).

- *v.com*:

  estimated commercial tree stem volume (m^3) according to the top
  diameter defined in the argument `d.mer`.

- *SS.max*:

  Maximum sagitta (see Prendes et al. (2022)).

- *sinuosity*:

  Sinuosity (see Prendes et al. (2022)).

- *lean*:

  Lean (see Prendes et al. (2022)).

- n.pts:

  number of points corresponding to a normal section (+/- 5 cm) in the
  original point cloud.

- n.pts.red:

  number of points corresponding to a normal section (+/- 5 cm) in the
  point cloud reduced by the point cropping process.

- n.pts.est:

  number of points estimated for a normal section (+/- 5 cm) in the
  original point cloud.

- n.pts.red.est:

  number of points estimated for a normal section (+/- 5 cm) in the
  point cloud reduced by the point cropping process.

- partial.occlusion:

  yes (1) or no (0)

## Output Files

At the end of the tree detection process, if the `save.result` argument
is `TRUE`, the function will save the data frame described in ‘Value’ as
a CSV file named ‘tree.tls.csv’. The data frame will be written without
row names in the `dir.result` directory by using
[`write.csv`](https://rdrr.io/r/utils/write.table.html) function from
the utils package.

## Note

Although `tree.detection.multi.scan` also works with reduced point
clouds, thus reducing the computing time, use of the original point
cloud is recommended in order to detect more trees. This will also
depend on forest conditions, especially those related to visibility. The
more distant the trees are, the lower the density of points will be, and
using reduced point clouds will therefore complicate detection of the
most distant trees.

Note that `dbh.min` and `dbh.max` are important for avoiding outlier
values when inventory data are used for reference purposes. Otherwise,
knowledge about the autoecology of species could be used for filtering
anomalous values of *dbh*.

The argument `breaks = 1.3` could be sufficient for detecting trees
visible at *dbh*, involving lower computational cost. However, those
trees not detected at *dbh*, may be estimated from lower and/or higher
sections. Considering the three default sections in the argument
`breaks = c(1.0, 1.3, 1.6)` maintains a good balance in the case study
of this package.

## References

Ester, M., Kriegel, H. P., Sander, J., & Xu, X. (1996). A density-based
algorithm for discovering clusters in large spatial databases with
noise. In Kdd (Vol. 96, No. 34, pp. 226-231).

Jin, S., Tamura, M., & Susaki, J. (2016). A new approach to retrieve
leaf normal distribution using terrestrial laser scanners. J. *Journal
of Forestry Research*, **27(3)**, 631-638.
[doi:10.1007/s11676-015-0204-z](https://doi.org/10.1007/s11676-015-0204-z)

Ma, L., Zheng, G., Eitel, J. U., Moskal, L. M., He, W., & Huang, H.
(2015). Improved salient feature-based approach for automatically
separating photosynthetic and nonphotosynthetic components within
terrestrial lidar point cloud data of forest canopies. *IEEE
Transactions Geoscience Remote Sensing*, **54(2)**, 679-696.
[doi:10.1109/TGRS.2015.2459716](https://doi.org/10.1109/TGRS.2015.2459716)

Pauly, M., Gross, M., & Kobbelt, L. P., (2002). Efficient simplification
of point-sampled surfaces. In IEEE Conference on Visualization. (pp.
163-170). Boston, USA.
[doi:10.1109/VISUAL.2002.1183771](https://doi.org/10.1109/VISUAL.2002.1183771)

Prendes, C., Canga, E., Ordoñez, C., Majada, J., Acuna, M., & Cabo, C.
(2022). Automatic assessment of individual stem shape parameters in
forest stands from tls point clouds: Application in pinus pinaster.
Forests, 13(3), 431.
[doi:10.3390/f13030431](https://doi.org/10.3390/f13030431)

Xia, S., Wang, C., Pan, F., Xi, X., Zeng, H., & Liu, H. (2015).
Detecting stems in dense and homogeneous forest using single-scan TLS.
*Forests*. **6(11)**, 3923-3945.
[doi:10.3390/f6113923](https://doi.org/10.3390/f6113923)

Zhang, W., Wan, P., Wang, T., Cai, S., Chen, Y., Jin, X., & Yan, G.
(2019). A novel approach for the detection of standing tree stems from
plot-level terrestrial laser scanning data. *Remote Sens*. **11(2)**,
211. [doi:10.3390/rs11020211](https://doi.org/10.3390/rs11020211)

## Author

Juan Alberto Molina-Valero and Adela Martínez-Calvo.

## See also

[`normalize`](https://molina-valero.github.io/FORTLS/reference/normalize.md),
[`tree.detection.single.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.single.scan.md),
[`tree.detection.several.plots`](https://molina-valero.github.io/FORTLS/reference/tree.detection.several.plots.md),
[`distance.sampling`](https://molina-valero.github.io/FORTLS/reference/distance.sampling.md),
[`estimation.plot.size`](https://molina-valero.github.io/FORTLS/reference/estimation.plot.size.md),
[`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md),
[`metrics.variables`](https://molina-valero.github.io/FORTLS/reference/metrics.variables.md)

## Examples

``` r
# \donttest{

# Establishment of working directories (optional)
# Establishment of working directories using tempdir() to avoid leaving files in
# user filespace

dir.data <- tempdir()
dir.result <- tempdir()


# Loading example data of TLS multi-scan approach point cloud (LAZ file) to dir.data

download.file(
"www.dropbox.com/scl/fi/es5pfj87wj0g6y8414dpo/PiceaAbies.laz?rlkey=ayt21mbndc6i6fyiz2e7z6oap&dl=1",
              destfile = file.path(dir.data, "PiceaAbies.laz"),
              mode = "wb")

# Normalizing the whole point cloud data without considering arguments

pcd <- normalize(las = "PiceaAbies.laz",

                 id = "PiceaAbies",

                 scan.approach = "multi",

                 max.dist = 7.5,

                 dir.data = dir.data, dir.result = dir.result)
#>                                                                                 


# Tree detection without considering arguments

tree.tls <- tree.detection.multi.scan(data = pcd,

                                      threads = 2,

                                      dir.data = dir.data, dir.result = dir.result)
#> Application of Statistical Outlier Removal (SOR) to the entire point cloud
#> Downloading uv...
#> Done!
#> Retention of points with high verticality & low surface variation
#> Computing geometric features...
#> [>                                                 ] 0% (0/41665)[=>                                                ] 2% (1000/41665)[==>                                               ] 4% (2000/41665)[===>                                              ] 7% (3000/41665)[====>                                             ] 9% (4000/41665)[======>                                           ] 12% (5000/41665)[=======>                                          ] 14% (6000/41665)[========>                                         ] 16% (7000/41665)[=========>                                        ] 19% (8000/41665)[==========>                                       ] 21% (9000/41665)[============>                                     ] 24% (10000/41665)[=============>                                    ] 26% (11000/41665)[==============>                                   ] 28% (12000/41665)[===============>                                  ] 31% (13000/41665)[================>                                 ] 33% (14000/41665)[==================>                               ] 36% (15000/41665)[===================>                              ] 38% (16000/41665)[====================>                             ] 40% (17000/41665)[=====================>                            ] 43% (18000/41665)[======================>                           ] 45% (19000/41665)[========================>                         ] 48% (20000/41665)[=========================>                        ] 50% (21000/41665)[==========================>                       ] 52% (22000/41665)[===========================>                      ] 55% (23000/41665)[============================>                     ] 57% (24000/41665)[==============================>                   ] 60% (25000/41665)[===============================>                  ] 62% (26000/41665)[================================>                 ] 64% (27000/41665)[=================================>                ] 67% (28000/41665)[==================================>               ] 69% (29000/41665)[====================================>             ] 72% (30000/41665)[=====================================>            ] 74% (31000/41665)[======================================>           ] 76% (32000/41665)[=======================================>          ] 79% (33000/41665)[========================================>         ] 81% (34000/41665)[==========================================>       ] 84% (35000/41665)[===========================================>      ] 86% (36000/41665)[============================================>     ] 88% (37000/41665)[=============================================>    ] 91% (38000/41665)[==============================================>   ] 93% (39000/41665)[================================================> ] 96% (40000/41665)[=================================================>] 98% (41000/41665)[==================================================] 100% (41665/41665)
#> Computing geometric features...
#> [>                                                 ] 0% (0/31916)[=>                                                ] 3% (1000/31916)[===>                                              ] 6% (2000/31916)[====>                                             ] 9% (3000/31916)[======>                                           ] 12% (4000/31916)[=======>                                          ] 15% (5000/31916)[=========>                                        ] 18% (6000/31916)[==========>                                       ] 21% (7000/31916)[============>                                     ] 25% (8000/31916)[==============>                                   ] 28% (9000/31916)[===============>                                  ] 31% (10000/31916)[=================>                                ] 34% (11000/31916)[==================>                               ] 37% (12000/31916)[====================>                             ] 40% (13000/31916)[=====================>                            ] 43% (14000/31916)[=======================>                          ] 46% (15000/31916)[=========================>                        ] 50% (16000/31916)[==========================>                       ] 53% (17000/31916)[============================>                     ] 56% (18000/31916)[=============================>                    ] 59% (19000/31916)[===============================>                  ] 62% (20000/31916)[================================>                 ] 65% (21000/31916)[==================================>               ] 68% (22000/31916)[====================================>             ] 72% (23000/31916)[=====================================>            ] 75% (24000/31916)[=======================================>          ] 78% (25000/31916)[========================================>         ] 81% (26000/31916)[==========================================>       ] 84% (27000/31916)[===========================================>      ] 87% (28000/31916)[=============================================>    ] 90% (29000/31916)[==============================================>   ] 93% (30000/31916)[================================================> ] 97% (31000/31916)[==================================================] 100% (31916/31916)
#> Detection of tree stem axes
#> Computing sections

  # }
```

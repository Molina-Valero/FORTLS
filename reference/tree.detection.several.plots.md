# Tree-Level Variables Estimation for Several Plots

This function integrates both, the
[`normalize`](https://molina-valero.github.io/FORTLS/reference/normalize.md)
and
[`tree.detection.single.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.single.scan.md)
or
[`tree.detection.multi.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.multi.scan.md)
functions, generating the same ‘Output Files’ as indicated for these,
and it returs the same ‘Value’ as described for
[`tree.detection.single.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.single.scan.md)
or
[`tree.detection.multi.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.multi.scan.md)
respectively. However, this function is designed for working with
several plots, producing a list of all scans considered automatically
from LAS files.

## Usage

``` r
tree.detection.several.plots(las.list, id.list = NULL, file = NULL,

                             scan.approach = "multi",
                             pcd.red = NULL, normalized = NULL,
                             center.coord = NULL,
                             x.side = NULL, y.side = NULL,
                             max.dist = NULL, min.height = NULL, max.height = 50,
                             algorithm.dtm = "knnidw", res.dtm = 0.2,
                             csf = list(cloth_resolution = 0.5),
                             intensity = NULL, RGB = NULL, voxel_size = NULL,

                             single.tree = NULL,
                             dbh.min = 4, dbh.max = 200, h.min = 1.3,
                             geo.dist = 0.1,
                             tls.resolution = NULL, tls.precision = NULL,
                             stem.section = c(0.7, 3.5), stem.range = NULL, breaks = NULL,
                             slice = 0.1, understory = NULL, bark.roughness = 1,
                             den.type = 1, d.mer = NULL,

                             plot.attributes = NULL, plot = NULL,

                             threads = 1,

                             dir.data = NULL, save.result = TRUE, dir.result = NULL,
                             save.las = NULL)
```

## Arguments

- las.list:

  Character vector containing the names of all LAS files for analysis
  and belonging to TLS point cloud, including .las extension (see
  ‘Examples’)

- id.list:

  Optional vector with plots identification encoded as character string
  or numeric. If this argument is not specified by the user, it will be
  set to NULL by default and, as a consequence, the plots will be
  encoded with correlative numbers from 1 to n plots.

- file:

  Optional vector containing files name identification encoded as
  character string or numeric value. If it is null, file will be encoded
  as `id` by default.

- scan.approach:

  Character parameter indicating TLS single-scan (‘single’) or TLS
  multi-scan approach or SLAM point clouds (‘multi’) approaches. If this
  argument is not specified by the user, it will be set to ‘multi’
  approach.

- pcd.red:

  Optional argument to indicate if point cloud density must be reduced
  to detect trees.

- normalized:

  Optional argument to establish as `TRUE` when point cloud is already
  normalized.

- center.coord:

  Planimetric x and y center coordinate of the plots. They has to be
  introduced as a data frame object with the following columns names:
  'id', 'x' and 'y'. They represent plot id, and center coordinates
  respectively.

- x.side:

  x-side (m) of the plot when the plot is square or rectangular.

- y.side:

  y-side (m) of the plot when the plot is square or rectangular.

- max.dist:

  Optional maximum horizontal distance (m) considered from the plot
  centre. All points farther than `max.dist` will be discarded after the
  normalization process. If this argument is not specified by the user,
  it will be set to NULL by default and, as a consequence, all points
  will be used in processing, with `max.dist` representing the farthest
  point.

- min.height:

  Optional minimum height (m) considered from ground level. All points
  below `min.height` will be discarded after the normalization process.
  If this argument is not specified by the user, it will be set to NULL
  by default and, as a consequence, all points will be used in
  processing, with `min.height` representing the lowest point.

- max.height:

  Optional maximum height (m) considered from ground level. All points
  above `max.height` will be discarded after the normalization process.
  If this argument is not specified by the user, it will be set to NULL
  by default and, as a consequence, all points will be used in
  processing, with `max.height` representing the highest point.

- algorithm.dtm:

  Algorithm used to generate the digital terrain model (DTM) from the
  TLS point cloud. There are two posible options based on spatial
  interpolation: ‘tin’ and ‘knnidw’ (see ‘Details’). If this argument is
  not specified by the user, it will be set to ‘knnidw’ algorithm.

- res.dtm:

  Numeric parameter. Resolution of the DTM generated to normalize point
  cloud (see ‘Details’). If this argument is not specified by the user,
  it will be set to 0.2 m.

- csf:

  List containing parameters of CSF algorithm:

&nbsp;

- `cloth_resolution`: by default 0.5.

&nbsp;

- intensity:

  Logical parameter useful when point clouds have intensity values. It
  may be useful in some internal process to filter data.

- RGB:

  Logical parameter useful when point clouds are colorized, thus
  including values of RGB colors. It is based on the Green Leaf
  Algorithm (GLA) (see ‘Details’).

- voxel_size:

  Defines the size of the 3D grid cells used for downsampling.

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

- tls.resolution:

  List containing parameters of TLS resolution. This can be defined by
  the angle aperture:

&nbsp;

- `horizontal.angle`: horizontal angle aperture (degrees).

- `vertical.angle`: vertical angle aperture (degrees).

- `point.dist`: distance (mm) between two consecutive points.

- `tls.dist`: distance (m) from TLS at which two consecutive points are
  separated by `point.dist`.

  If this argument is not specified by the user, it will be set to NULL
  by default and, as a consequence the function will stop giving an
  error message.

&nbsp;

- tls.precision:

  Optional argument indicating the average point cloud precision in cm.

- stem.section:

  Section free of noise (shurb, branches, etc.) considered to detect
  trees. If not specified, an automatic internal algorithm will be
  applied (see ‘Details’).

- breaks:

  Height above ground level (m) of slices considered for detecting
  trees. By default it will be considered all possible sections from 0.1
  m to maximum height by 0.3 m intervals (+/- 5 cm).

- stem.range:

  Section considered to estimate straightness tree attributes.

- slice:

  Slice width considered for detecting trees. By default it will be
  considered as 0.1 m.

- understory:

  Optional argument to indicate if there is dense understory vegetation.

- bark.roughness:

  Bark roughness established in 3 degrees (1 \< 2 \< 3). By default it
  will be considered as 2.

- den.type:

  Numeric argument indicating the dendrometic type used to estimate
  volumen when there are not sections enough to fit a taper equation.
  Dendrometrics types available are the following: cylinder = 0,
  paraboloid = 1 (by default), cone = 2 and neiloid = 3.

- d.mer:

  Top stem diameter (cm) considered to estimate commercial timber
  volume.

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
  `.Platform$file.sep` must be used as the path separator in `dir.dat`.
  If this argument is not specified by the user, it will be set to
  `NULL` by default and, as consequence, the current working directory
  of the R process will be assigned to `dir.dat` during the execution.

- save.result:

  Optional logical which indicates whether or not the output files
  described in ‘Output Files’ section should be saved in the
  `dir.result`. If this argument is not specified by the user, it will
  be set to `TRUE` by default and, as a consequence, the output files
  will be saved.

- dir.result:

  Optional character string naming the absolute path of an existing
  directory where files described in ‘Output Files’ section will be
  saved. `.Platform$file.sep` must be used as the path separator in
  `dir.result`. If this argument is not specified by the user, and
  `save.result` is `TRUE`, it will be set to `NULL` by default and, as a
  consequence, the current working directory of the R process will be
  assigned to `dir.result` during the execution.

- save.las:

  Optional logical which indicates whether or not the imput point cloud
  must be saved in `dir.result` as LAZ file.

## Details

See
[`normalize`](https://molina-valero.github.io/FORTLS/reference/normalize.md),
[`tree.detection.single.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.single.scan.md)
and
[`tree.detection.multi.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.multi.scan.md)
for further details.

## Value

Data frame with the same description and format as
[`tree.detection.single.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.single.scan.md)
and
[`tree.detection.multi.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.multi.scan.md)
‘Values’. In this case, the `id` of plots will be encoded with
correlative numbers from 1 to n, where n is the number of LAS files
included in `files` argument, and `file` column will be encoded as `id`,
but including .las extension.

## Output Files

At the end of the tree detection process, if the `save.result` argument
is `TRUE`, the function will save both, the reduced point clouds as TXT
files encoded according to `file` column of ‘Value’; and the data frame
with the tree list described in ‘Value’ as CSV file (see
[`normalize`](https://molina-valero.github.io/FORTLS/reference/normalize.md)
and
[`tree.detection.single.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.single.scan.md)
or
[`tree.detection.multi.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.multi.scan.md)
‘Output files’). All outputs are written without row names in the
`dir.result` directory using
[vroom_write](https://vroom.tidyverse.org/reference/vroom_write.html)
function from [vroom](https://CRAN.R-project.org/package=vroom) package.

## Note

This function has been developed for working with several plots, which
will be the most common situation in forest inventory approaches.
Nevertheless, several LAS files are not provided as examples due to
problems with memory capacity.

## Author

Juan Alberto Molina-Valero and Adela Martínez-Calvo.

## See also

[`normalize`](https://molina-valero.github.io/FORTLS/reference/normalize.md),
[`tree.detection.single.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.single.scan.md),[`tree.detection.multi.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.multi.scan.md),
[`distance.sampling`](https://molina-valero.github.io/FORTLS/reference/distance.sampling.md),
[`estimation.plot.size`](https://molina-valero.github.io/FORTLS/reference/estimation.plot.size.md),
[`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md),
[`metrics.variables`](https://molina-valero.github.io/FORTLS/reference/metrics.variables.md).

## Examples

``` r
# \donttest{

# Establishment of working directories (optional)
# Establishment of working directories using tempdir() to avoid leaving files in
# user filespace

dir.data <- tempdir()
dir.result <- tempdir()

# Loading example data (LAZ files) to dir.data

download.file(
"www.dropbox.com/scl/fi/hzzrt0a39crdy6uvcj9el/PinusSylve1.laz?rlkey=svpwvorkm8889fgbnj14ns1f2&dl=1",
              destfile = file.path(dir.data, "PinusSylvestris1.laz"),
              mode = "wb")

download.file(
"www.dropbox.com/scl/fi/zeszze31jh5m1g4o3ns1o/PinusSylve2.laz?rlkey=wx72bi6ggdc7wedwgzupekp9k&dl=1",
              destfile = file.path(dir.data, "PinusSylvestris2.laz"),
              mode = "wb")


# Tree detection (TLS single-scan aproach)

id <- c("PinusSylvestris1", "PinusSylvestris2")

center.coord <- data.frame(id = id,
                           x = rep(0, length(id)),
                           y = rep(0, length(id)))

tree.tls <- tree.detection.several.plots(las.list = c("PinusSylvestris1.laz",
                                                      "PinusSylvestris2.laz"),

                                         id.list = id,

                                         scan.approach = "single",

                                         center.coord = center.coord,

                                         tls.resolution = list(point.dist = 7.67, tls.dist = 10),

                                         max.dist = 7.5,

                                         threads = 2,

                                         dir.data = dir.data, dir.result = dir.result)
#> Computing plot: PinusSylvestris1
#> Normalizing
#> [===================================>              ] 71% ETA: 0s     [===================================>              ] 71% ETA: 0s     [===================================>              ] 71% ETA: 0s     [===================================>              ] 71% ETA: 0s     [===================================>              ] 71% ETA: 0s     [===================================>              ] 71% ETA: 0s     [===================================>              ] 71% ETA: 0s     [===================================>              ] 71% ETA: 0s     [====================================>             ] 72% ETA: 0s     [====================================>             ] 72% ETA: 0s     [====================================>             ] 72% ETA: 0s     [====================================>             ] 72% ETA: 0s     [====================================>             ] 72% ETA: 0s     [====================================>             ] 72% ETA: 0s     [====================================>             ] 72% ETA: 0s     [====================================>             ] 72% ETA: 0s     [====================================>             ] 72% ETA: 0s     [====================================>             ] 72% ETA: 0s     [====================================>             ] 72% ETA: 0s     [====================================>             ] 72% ETA: 0s     [====================================>             ] 72% ETA: 0s     [====================================>             ] 73% ETA: 0s     [====================================>             ] 73% ETA: 0s     [====================================>             ] 73% ETA: 0s     [====================================>             ] 73% ETA: 0s     [====================================>             ] 73% ETA: 0s     [====================================>             ] 73% ETA: 0s     [====================================>             ] 73% ETA: 0s     [====================================>             ] 73% ETA: 0s     [====================================>             ] 73% ETA: 0s     [====================================>             ] 73% ETA: 0s     [====================================>             ] 73% ETA: 0s     [====================================>             ] 73% ETA: 0s     [====================================>             ] 73% ETA: 0s     [====================================>             ] 73% ETA: 0s     [=====================================>            ] 74% ETA: 0s     [=====================================>            ] 74% ETA: 0s     [=====================================>            ] 74% ETA: 0s     [=====================================>            ] 74% ETA: 0s     [=====================================>            ] 74% ETA: 0s     [=====================================>            ] 74% ETA: 0s     [=====================================>            ] 74% ETA: 0s     [=====================================>            ] 74% ETA: 0s     [=====================================>            ] 74% ETA: 0s     [=====================================>            ] 74% ETA: 0s     [=====================================>            ] 74% ETA: 0s     [=====================================>            ] 74% ETA: 0s     [=====================================>            ] 74% ETA: 0s     [=====================================>            ] 74% ETA: 0s     [=====================================>            ] 75% ETA: 0s     [=====================================>            ] 75% ETA: 0s     [=====================================>            ] 75% ETA: 0s     [=====================================>            ] 75% ETA: 0s     [=====================================>            ] 75% ETA: 0s     [=====================================>            ] 75% ETA: 0s     [=====================================>            ] 75% ETA: 0s     [=====================================>            ] 75% ETA: 0s     [=====================================>            ] 75% ETA: 0s     [=====================================>            ] 75% ETA: 0s     [=====================================>            ] 75% ETA: 0s     [=====================================>            ] 75% ETA: 0s     [=====================================>            ] 75% ETA: 0s     [=====================================>            ] 75% ETA: 0s     [======================================>           ] 76% ETA: 0s     [======================================>           ] 76% ETA: 0s     [======================================>           ] 76% ETA: 0s     [======================================>           ] 76% ETA: 0s     [======================================>           ] 76% ETA: 0s     [======================================>           ] 76% ETA: 0s     [======================================>           ] 76% ETA: 0s     [======================================>           ] 76% ETA: 0s     [======================================>           ] 76% ETA: 0s     [======================================>           ] 76% ETA: 0s     [======================================>           ] 76% ETA: 0s     [======================================>           ] 76% ETA: 0s     [======================================>           ] 76% ETA: 0s     [======================================>           ] 76% ETA: 0s     [======================================>           ] 77% ETA: 0s     [======================================>           ] 77% ETA: 0s     [======================================>           ] 77% ETA: 0s     [======================================>           ] 77% ETA: 0s     [======================================>           ] 77% ETA: 0s     [======================================>           ] 77% ETA: 0s     [======================================>           ] 77% ETA: 0s     [======================================>           ] 77% ETA: 0s     [======================================>           ] 77% ETA: 0s     [======================================>           ] 77% ETA: 0s     [======================================>           ] 77% ETA: 0s     [======================================>           ] 77% ETA: 0s     [======================================>           ] 77% ETA: 0s     [======================================>           ] 77% ETA: 0s     [=======================================>          ] 78% ETA: 0s     [=======================================>          ] 78% ETA: 0s     [=======================================>          ] 78% ETA: 0s     [=======================================>          ] 78% ETA: 0s     [=======================================>          ] 78% ETA: 0s     [=======================================>          ] 78% ETA: 0s     [=======================================>          ] 78% ETA: 0s     [=======================================>          ] 78% ETA: 0s     [=======================================>          ] 78% ETA: 0s     [=======================================>          ] 78% ETA: 0s     [=======================================>          ] 78% ETA: 0s     [=======================================>          ] 78% ETA: 0s     [=======================================>          ] 78% ETA: 0s     [=======================================>          ] 78% ETA: 0s     [=======================================>          ] 79% ETA: 0s     [=======================================>          ] 79% ETA: 0s     [=======================================>          ] 79% ETA: 0s     [=======================================>          ] 79% ETA: 0s     [=======================================>          ] 79% ETA: 0s     [=======================================>          ] 79% ETA: 0s     [=======================================>          ] 79% ETA: 0s     [=======================================>          ] 79% ETA: 0s     [=======================================>          ] 79% ETA: 0s     [=======================================>          ] 79% ETA: 0s     [=======================================>          ] 79% ETA: 0s     [=======================================>          ] 79% ETA: 0s     [=======================================>          ] 79% ETA: 0s     [=======================================>          ] 79% ETA: 0s     [========================================>         ] 80% ETA: 0s     [========================================>         ] 80% ETA: 0s     [========================================>         ] 80% ETA: 0s     [========================================>         ] 80% ETA: 0s     [========================================>         ] 80% ETA: 0s     [========================================>         ] 80% ETA: 0s     [========================================>         ] 80% ETA: 0s     [========================================>         ] 80% ETA: 0s     [========================================>         ] 80% ETA: 0s     [========================================>         ] 80% ETA: 0s     [========================================>         ] 80% ETA: 0s     [========================================>         ] 80% ETA: 0s     [========================================>         ] 80% ETA: 0s     [========================================>         ] 81% ETA: 0s     [========================================>         ] 81% ETA: 0s     [========================================>         ] 81% ETA: 0s     [========================================>         ] 81% ETA: 0s     [========================================>         ] 81% ETA: 0s     [========================================>         ] 81% ETA: 0s     [========================================>         ] 81% ETA: 0s     [========================================>         ] 81% ETA: 0s     [========================================>         ] 81% ETA: 0s     [========================================>         ] 81% ETA: 0s     [========================================>         ] 81% ETA: 0s     [========================================>         ] 81% ETA: 0s     [========================================>         ] 81% ETA: 0s     [========================================>         ] 81% ETA: 0s     [=========================================>        ] 82% ETA: 0s     [=========================================>        ] 82% ETA: 0s     [=========================================>        ] 82% ETA: 0s     [=========================================>        ] 82% ETA: 0s     [=========================================>        ] 82% ETA: 0s     [=========================================>        ] 82% ETA: 0s     [=========================================>        ] 82% ETA: 0s     [=========================================>        ] 82% ETA: 0s     [=========================================>        ] 82% ETA: 0s     [=========================================>        ] 82% ETA: 0s     [=========================================>        ] 82% ETA: 0s     [=========================================>        ] 82% ETA: 0s     [=========================================>        ] 82% ETA: 0s     [=========================================>        ] 82% ETA: 0s     [=========================================>        ] 83% ETA: 0s     [=========================================>        ] 83% ETA: 0s     [=========================================>        ] 83% ETA: 0s     [=========================================>        ] 83% ETA: 0s     [=========================================>        ] 83% ETA: 0s     [=========================================>        ] 83% ETA: 0s     [=========================================>        ] 83% ETA: 0s     [=========================================>        ] 83% ETA: 0s     [=========================================>        ] 83% ETA: 0s     [=========================================>        ] 83% ETA: 0s     [=========================================>        ] 83% ETA: 0s     [=========================================>        ] 83% ETA: 0s     [=========================================>        ] 83% ETA: 0s     [=========================================>        ] 83% ETA: 0s     [==========================================>       ] 84% ETA: 0s     [==========================================>       ] 84% ETA: 0s     [==========================================>       ] 84% ETA: 0s     [==========================================>       ] 84% ETA: 0s     [==========================================>       ] 84% ETA: 0s     [==========================================>       ] 84% ETA: 0s     [==========================================>       ] 84% ETA: 0s     [==========================================>       ] 84% ETA: 0s     [==========================================>       ] 84% ETA: 0s     [==========================================>       ] 84% ETA: 0s     [==========================================>       ] 84% ETA: 0s     [==========================================>       ] 84% ETA: 0s     [==========================================>       ] 84% ETA: 0s     [==========================================>       ] 84% ETA: 0s     [==========================================>       ] 85% ETA: 0s     [==========================================>       ] 85% ETA: 0s     [==========================================>       ] 85% ETA: 0s     [==========================================>       ] 85% ETA: 0s     [==========================================>       ] 85% ETA: 0s     [==========================================>       ] 85% ETA: 0s     [==========================================>       ] 85% ETA: 0s     [==========================================>       ] 85% ETA: 0s     [==========================================>       ] 85% ETA: 0s     [==========================================>       ] 85% ETA: 0s     [==========================================>       ] 85% ETA: 0s     [==========================================>       ] 85% ETA: 0s     [==========================================>       ] 85% ETA: 0s     [==========================================>       ] 85% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 89% ETA: 0s     [============================================>     ] 89% ETA: 0s     [============================================>     ] 89% ETA: 0s     [============================================>     ] 89% ETA: 0s     [============================================>     ] 89% ETA: 0s     [============================================>     ] 89% ETA: 0s     [============================================>     ] 89% ETA: 0s     [============================================>     ] 89% ETA: 0s     [============================================>     ] 89% ETA: 0s     [============================================>     ] 89% ETA: 0s     [============================================>     ] 89% ETA: 0s     [============================================>     ] 89% ETA: 0s     [============================================>     ] 89% ETA: 0s     [============================================>     ] 89% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 97% ETA: 0s     [================================================> ] 97% ETA: 0s     [================================================> ] 97% ETA: 0s     [================================================> ] 97% ETA: 0s     [================================================> ] 97% ETA: 0s     [================================================> ] 97% ETA: 0s     [================================================> ] 97% ETA: 0s     [================================================> ] 97% ETA: 0s     [================================================> ] 97% ETA: 0s     [================================================> ] 97% ETA: 0s     [================================================> ] 97% ETA: 0s     [================================================> ] 97% ETA: 0s     [================================================> ] 97% ETA: 0s     [================================================> ] 97% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 99% ETA: 0s     [=================================================>] 99% ETA: 0s     [=================================================>] 99% ETA: 0s     [=================================================>] 99% ETA: 0s     [=================================================>] 99% ETA: 0s     [=================================================>] 99% ETA: 0s     [=================================================>] 99% ETA: 0s     [=================================================>] 99% ETA: 0s     [=================================================>] 99% ETA: 0s     [=================================================>] 99% ETA: 0s     [=================================================>] 99% ETA: 0s     [=================================================>] 99% ETA: 0s     [=================================================>] 99% ETA: 0s     [=================================================>] 99% ETA: 0s                                                                                     
#> Detecting trees
#> Application of Statistical Outlier Removal (SOR) to the entire point cloud
#> Retention of points with high verticality & low surface variation
#> Computing geometric features...
#> [>                                                 ] 0% (0/16797)[==>                                               ] 5% (1000/16797)[=====>                                            ] 11% (2000/16797)[========>                                         ] 17% (3000/16797)[===========>                                      ] 23% (4000/16797)[==============>                                   ] 29% (5000/16797)[=================>                                ] 35% (6000/16797)[====================>                             ] 41% (7000/16797)[=======================>                          ] 47% (8000/16797)[==========================>                       ] 53% (9000/16797)[=============================>                    ] 59% (10000/16797)[================================>                 ] 65% (11000/16797)[===================================>              ] 71% (12000/16797)[======================================>           ] 77% (13000/16797)[=========================================>        ] 83% (14000/16797)[============================================>     ] 89% (15000/16797)[===============================================>  ] 95% (16000/16797)[==================================================] 100% (16797/16797)
#> Computing sections
#> Computing plot: PinusSylvestris2
#> Normalizing
#> [==========================================>       ] 85% ETA: 0s     [==========================================>       ] 85% ETA: 0s     [==========================================>       ] 85% ETA: 0s     [==========================================>       ] 85% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 86% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [===========================================>      ] 87% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 88% ETA: 0s     [============================================>     ] 89% ETA: 0s     [============================================>     ] 89% ETA: 0s     [============================================>     ] 89% ETA: 0s     [============================================>     ] 89% ETA: 0s     [============================================>     ] 89% ETA: 0s     [============================================>     ] 89% ETA: 0s     [============================================>     ] 89% ETA: 0s     [============================================>     ] 89% ETA: 0s     [============================================>     ] 89% ETA: 0s     [============================================>     ] 89% ETA: 0s     [============================================>     ] 89% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 90% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [=============================================>    ] 91% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 92% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [==============================================>   ] 93% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 94% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [===============================================>  ] 95% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 96% ETA: 0s     [================================================> ] 97% ETA: 0s     [================================================> ] 97% ETA: 0s     [================================================> ] 97% ETA: 0s     [================================================> ] 97% ETA: 0s     [================================================> ] 97% ETA: 0s     [================================================> ] 97% ETA: 0s     [================================================> ] 97% ETA: 0s     [================================================> ] 97% ETA: 0s     [================================================> ] 97% ETA: 0s     [================================================> ] 97% ETA: 0s     [================================================> ] 97% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 98% ETA: 0s     [=================================================>] 99% ETA: 0s     [=================================================>] 99% ETA: 0s     [=================================================>] 99% ETA: 0s     [=================================================>] 99% ETA: 0s     [=================================================>] 99% ETA: 0s     [=================================================>] 99% ETA: 0s     [=================================================>] 99% ETA: 0s     [=================================================>] 99% ETA: 0s     [=================================================>] 99% ETA: 0s     [=================================================>] 99% ETA: 0s     [=================================================>] 99% ETA: 0s                                                                                     
#> Detecting trees
#> Application of Statistical Outlier Removal (SOR) to the entire point cloud
#> Retention of points with high verticality & low surface variation
#> Computing geometric features...
#> [>                                                 ] 0% (0/11583)[====>                                             ] 8% (1000/11583)[========>                                         ] 17% (2000/11583)[============>                                     ] 25% (3000/11583)[=================>                                ] 34% (4000/11583)[=====================>                            ] 43% (5000/11583)[=========================>                        ] 51% (6000/11583)[==============================>                   ] 60% (7000/11583)[==================================>               ] 69% (8000/11583)[======================================>           ] 77% (9000/11583)[===========================================>      ] 86% (10000/11583)[===============================================>  ] 94% (11000/11583)[==================================================] 100% (11583/11583)
#> Computing sections

  # }
```

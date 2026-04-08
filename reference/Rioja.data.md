# Inventoried Plots Data for a Stand Case Study in La Rioja

This list includes trees detected with TLS for 16 single scans
corresponding to plots located in La Rioja, a region of Spain, in the
north of the Iberian Peninsula (first element), as well as those
inventoried in the field for these 16 plots (second element). Plot
attributes related to stand stratum are also included (third element).

The elements of the list are as follows:

1.  `tree.tls`: data frame that includes the list of trees detected with
    [`tree.detection.single.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.single.scan.md)
    for 16 TLS single-scan sampling points. The following variables are
    provided for each tree (see
    [`tree.detection.single.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.single.scan.md)
    ‘Value’ for more details):

    |            |                                            |                   |
    |------------|--------------------------------------------|-------------------|
    | \[,1\]     | id                                         | character/numeric |
    | \[,2\]     | file                                       | character         |
    | \[,3\]     | tree                                       | numeric           |
    | \[,4\]     | x                                          | numeric           |
    | \[,5\]     | y                                          | numeric           |
    | \[,6:8\]   | phi, phi.left, phi.right                   | numeric           |
    | \[,9\]     | h.dist                                     | numeric           |
    | \[,10\]    | dbh                                        | numeric           |
    | \[,11\]    | h                                          | numeric           |
    | \[,12\]    | v                                          | numeric           |
    | \[,13:16\] | n.pts, n.pts.red, n.pts.est, n.pts.red.est | numeric           |
    | \[,15\]    | partial.occlusion                          | numeric           |

2.  `tree.list.field`: data frame that includes the list of trees
    measured in 16 circular fixed area plots of radius 20 m, whose
    centres coincide with TLS single-scans points. The following
    variables are provided for each tree:

    |         |        |         |                                                         |
    |---------|--------|---------|---------------------------------------------------------|
    | \[,1\]  | id     | numeric | plot identification (coincident to TLS scans)           |
    | \[,2\]  | tree   | numeric | trees numbering                                         |
    | \[,3\]  | Sp     | numeric | specie code according to NFI ()                         |
    | \[,4\]  | x      | numeric | x cartesian coordinate                                  |
    | \[,5\]  | y      | numeric | 7 cartesian coordinate                                  |
    | \[,6\]  | h.dist | numeric | horizontal distance (m) from plot center to tree center |
    | \[,7\]  | dbh    | numeric | tree diameter (cm) at breast height (1.3 m)             |
    | \[,8\]  | h      | numeric | tree total height (m)                                   |
    | \[,9\]  | dead   | numeric | dead (1) or not (NA)                                    |
    | \[,10\] | v.user | numeric | stem volume (m^3) estimated with allometric equations   |
    | \[,11\] | w.user | numeric | stem biomass (Mg) estimated with allometric equations   |

## Usage

``` r
data(Rioja.data)
```

## Format

List with 2 data frames containing 604 observations and 17 variables
(tree.tls) and 659 observations and 11 variables (tree.field).

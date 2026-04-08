# Simulated Metrics and Variables for a Stand Case Study in La Rioja

This list contains metrics and variables estimated from field data and
TLS data from
[`Rioja.data`](https://molina-valero.github.io/FORTLS/reference/Rioja.data.md).

The elements on this list correspond to
[`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md)
‘Value’, as follows:

1.  `fixed.area`: data frame with TLS metrics and variables estimated on
    the basis of simulated plots in a fixed area plot design with radius
    increment of 0.1 m (from smallest possible radius to 20 m). The
    following variables are provided for each pair (plot, radius) (see
    [`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md)
    ‘Value’ for more details):

    |            |                                                                           |                   |
    |------------|---------------------------------------------------------------------------|-------------------|
    | \[,1\]     | id                                                                        | character/numeric |
    | \[,2\]     | radius                                                                    | numeric           |
    | \[,3:5\]   | N, G, V                                                                   | numeric           |
    | \[,6:9\]   | d, dg, dgeom, dharm                                                       | numeric           |
    | \[,10:13\] | h, hg, hgeom, hharm                                                       | numeric           |
    | \[,14:17\] | d.0, dg.0, dgeom.0, dharm.0                                               | numeric           |
    | \[,18:21\] | h.0, hg.0, hgeom.0, hharm.0                                               | numeric           |
    | \[,22:27\] | N.tls, N.hn.tls, N.hr.tls, N.hn.cov.tls, N.hr.cov.tls, N.sh.tls           | numeric           |
    | \[,28:31\] | num.points, num.points.est, num.points.hom, num.points.hom.est            | numeric           |
    | \[,32:37\] | G.tls, G.hn.tls, G.hr.tls, G.hn.cov.tls, G.hr.cov.tls, G.sh.tls           | numeric           |
    | \[,38:43\] | V.tls, V.hn.tls, V.hr.tls, V.hn.cov.tls, V.hr.cov.tls, V.sh.tls           | numeric           |
    | \[,44:47\] | d.tls, dg.tls, dgeom.tls, dharm.tls                                       | numeric           |
    | \[,48:51\] | h.tls, hg.tls, hgeom.tls, hharm.tls                                       | numeric           |
    | \[,52:55\] | d.0, dg.0, dgeom.0, dharm.0                                               | numeric           |
    | \[,56:59\] | h.0, hg.0, hgeom.0, hharm.0                                               | numeric           |
    | \[,60:74\] | P01, P05, P10, P20, P25, P30, P40, P50, P60, P70, P75, P80, P90, P95, P99 | numeric           |

2.  `k.tree`: data frame with TLS metrics and variables estimated on the
    basis of simulated plots under k-tree plot design for incremental
    values of 1 tree (from 1 to largest number of trees in one plot).
    The following variables are provided for each pair (plot, k) (see
    [`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md)
    ‘Value’ for more details):

    |            |                                                                           |                   |
    |------------|---------------------------------------------------------------------------|-------------------|
    | \[,1\]     | id                                                                        | character/numeric |
    | \[,2\]     | k                                                                         | numeric           |
    | \[,3:5\]   | N, G, V                                                                   | numeric           |
    | \[,6:9\]   | d, dg, dgeom, dharm                                                       | numeric           |
    | \[,10:13\] | h, hg, hgeom, hharm                                                       | numeric           |
    | \[,14:17\] | d.0, dg.0, dgeom.0, dharm.0                                               | numeric           |
    | \[,18:21\] | h.0, hg.0, hgeom.0, hharm.0                                               | numeric           |
    | \[,22:27\] | N.tls, N.hn.tls, N.hr.tls, N.hn.cov.tls, N.hr.cov.tls, N.sh.tls           | numeric           |
    | \[,28:31\] | num.points, num.points.est, num.points.hom, num.points.hom.est            | numeric           |
    | \[,32:37\] | G.tls, G.hn.tls, G.hr.tls, G.hn.cov.tls, G.hr.cov.tls, G.sh.tls           | numeric           |
    | \[,38:43\] | V.tls, V.hn.tls, V.hr.tls, V.hn.cov.tls, V.hr.cov.tls, V.sh.tls           | numeric           |
    | \[,44:47\] | d.tls, dg.tls, dgeom.tls, dharm.tls                                       | numeric           |
    | \[,48:51\] | h.tls, hg.tls, hgeom.tls, hharm.tls                                       | numeric           |
    | \[,52:55\] | d.0, dg.0, dgeom.0, dharm.0                                               | numeric           |
    | \[,56:59\] | h.0, hg.0, hgeom.0, hharm.0                                               | numeric           |
    | \[,60:74\] | P01, P05, P10, P20, P25, P30, P40, P50, P60, P70, P75, P80, P90, P95, P99 | numeric           |

3.  `angle.count`: data frame with TLS metrics and variables estimated
    on the basis of simulated plots in an angle-count plot design. They
    plots are simulated for correlative angle-count plots and
    incremental values of 0.1 \\{m}^{2}/ha\\ for BAF. The following
    variables are provided for each pair (plot, BAF) (see
    [`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md)
    ‘Value’ for more details):

    |            |                                                                           |                   |
    |------------|---------------------------------------------------------------------------|-------------------|
    | \[,1\]     | id                                                                        | character/numeric |
    | \[,2\]     | BAF                                                                       | numeric           |
    | \[,3:5\]   | N, G, V                                                                   | numeric           |
    | \[,6:9\]   | d, dg, dgeom, dharm                                                       | numeric           |
    | \[,10:13\] | h, hg, hgeom, hharm                                                       | numeric           |
    | \[,14:17\] | d.0, dg.0, dgeom.0, dharm.0                                               | numeric           |
    | \[,18:21\] | h.0, hg.0, hgeom.0, hharm.0                                               | numeric           |
    | \[,22:23\] | N.tls, N.pam.tls                                                          | numeric           |
    | \[,24:27\] | num.points, num.points.est, num.points.hom, num.points.hom.est            | numeric           |
    | \[,28:29\] | G.tls, G.pam.tls                                                          | numeric           |
    | \[,30:31\] | V.tls, V.pam.tls                                                          | numeric           |
    | \[,32:35\] | d.tls, dg.tls, dgeom.tls, dharm.tls                                       | numeric           |
    | \[,48:51\] | h.tls, hg.tls, hgeom.tls, hharm.tls                                       | numeric           |
    | \[,36:39\] | d.0, dg.0, dgeom.0, dharm.0                                               | numeric           |
    | \[,40:43\] | h.0, hg.0, hgeom.0, hharm.0                                               | numeric           |
    | \[,44:62\] | P01, P05, P10, P20, P25, P30, P40, P50, P60, P70, P75, P80, P90, P95, P99 | numeric           |

## Usage

``` r
data(Rioja.simulations)
```

## Format

List with 3 data frames containing 2224 observations and 74 variables
(simulations.fixed.area.plot), 272 observations and 74 variables
(simulations.k.tree.plot), and 576 observations and 62 variables
(simulations.angle.count.plot).

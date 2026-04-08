# FORTLS: Automatic Processing of Terrestrial-Based Technologies Point Cloud Data for Forestry Purposes

Process automation of point cloud data derived from terrestrial-based
technologies such as Terrestrial Laser Scanner (TLS) or Mobile Laser
Scanner. 'FORTLS' enables (i) detection of trees and estimation of
tree-level attributes (e.g. diameters and heights), (ii) estimation of
stand-level variables (e.g. density, basal area, mean and dominant
height), (iii) computation of metrics related to important forest
attributes estimated in Forest Inventories at stand-level, and (iv)
optimization of plot design for combining TLS data and field measured
data. Documentation about 'FORTLS' is described in Molina-Valero et al.
(2022, \<doi:10.1016/j.envsoft.2022.105337\>).

## Details

Usage of FORTLS includes the following functionalities:

- Tree detection: this is the first and necessary step for the other
  functionalities of FORTLS. This can be achieved using the following
  functions:

  1.  [`normalize`](https://molina-valero.github.io/FORTLS/reference/normalize.md):
      mandatory first step for obtaining the relative coordinates of a
      TLS point cloud.

  2.  [`tree.detection.single.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.single.scan.md):
      detects as many trees as possible from a normalized TLS
      single-scan point clouds.

  3.  [`tree.detection.multi.scan`](https://molina-valero.github.io/FORTLS/reference/tree.detection.multi.scan.md):
      detects as many trees as possible from a normalized TLS
      multi-scan, SLAM, or similar terrestrial-based technologies point
      clouds.

  4.  [`tree.detection.several.plots`](https://molina-valero.github.io/FORTLS/reference/tree.detection.several.plots.md):
      includes the two previous functions for a better workflow when
      there are several plots to be sequentially analyzed.

- Estimation of variables when no field data are available: this is the
  main functionality of FORTLS and can be achieved using the following
  functions:

  1.  [`distance.sampling`](https://molina-valero.github.io/FORTLS/reference/distance.sampling.md):
      optional function which can be used for considering methodologies
      for correcting occlusion effects in estimating variables.

  2.  [`estimation.plot.size`](https://molina-valero.github.io/FORTLS/reference/estimation.plot.size.md):
      enables the best plot design to be determined on the basis of TLS
      data only.

  3.  [`metrics.variables`](https://molina-valero.github.io/FORTLS/reference/metrics.variables.md):
      is used for estimating metrics and variables potentially related
      to forest attributes at stand level.

- Estimation of variables when field data are available: this is the
  main and most desirable functionality of FORTLS and can be achieved
  using the following functions:

  1.  [`distance.sampling`](https://molina-valero.github.io/FORTLS/reference/distance.sampling.md):
      as before.

  2.  [`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md):
      computes simulations of TLS and field data for different plot
      designs. This is a prior step to the next functions.

  3.  [`relative.bias`](https://molina-valero.github.io/FORTLS/reference/relative.bias.md):
      uses
      [`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md)
      output to assess the accuracy of direct estimations of variables
      according to homologous TLS and field data.

  4.  [`correlations`](https://molina-valero.github.io/FORTLS/reference/correlations.md):
      uses
      [`simulations`](https://molina-valero.github.io/FORTLS/reference/simulations.md)
      output to assess correlations among metrics and variables obtained
      from TLS data, and variables of interest estimated from field
      data.

  5.  [`optimize.plot.design`](https://molina-valero.github.io/FORTLS/reference/optimize.plot.design.md):
      using
      [`correlations`](https://molina-valero.github.io/FORTLS/reference/correlations.md)
      output, represents the best correlations for variables of interest
      according to the plot design. It is thus possible to select the
      best plot design for estimating forest attributes from TLS data.

  6.  [`metrics.variables`](https://molina-valero.github.io/FORTLS/reference/metrics.variables.md):
      as before, but in this case plot parameters will be choosen on the
      basis of field data and better estimates will therefore be
      obtained.

## Author

**Maintainer**: Juan Alberto Molina-Valero
<juanalberto.molina.valero@usc.es> \[copyright holder\]

Authors:

- María José Ginzo Villamayor \[contributor\]

- Manuel Antonio Novo Pérez \[contributor\]

- Adela Martínez-Calvo \[contributor\]

- Juan Gabriel Álvarez-González \[contributor\]

- Fernando Montes \[contributor\]

- César Pérez-Cruzado \[contributor\]

## References

Molina-Valero J. A., Ginzo-Villamayor M. J., Novo Pérez M. A.,
Martínez-Calvo A., Álvarez-González J. G., Montes F., & Pérez-Cruzado C.
(2019). FORTLS: an R package for processing TLS data and estimating
stand variables in forest inventories. *The 1st International Electronic
Conference on Forests — Forests for a Better Future: Sustainability,
Innovation, Interdisciplinarity*.
[doi:10.3390/IECF2020-08066](https://doi.org/10.3390/IECF2020-08066)

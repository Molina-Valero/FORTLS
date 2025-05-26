<!-- README.md is generated from README.Rmd. Please edit that file -->

# FORTLS <img src="man/figures/logo.png" align="right" alt="" width="250" />

<!-- badges: start -->

![license](https://img.shields.io/badge/Licence-GPL--3-blue.svg)
[![CRAN Status](https://www.r-pkg.org/badges/version/FORTLS)](https://cran.r-project.org/package=FORTLS)
[![DOI](https://zenodo.org/badge/DOI/10.3390/IECF2020-08066.svg)](https://doi.org/10.1016/j.envsoft.2022.105337)
![](https://cranlogs.r-pkg.org/badges/grand-total/FORTLS)

<!-- badges: end -->

## Automatic Processing of Close-Range Technologies Point Cloud Data for Forestry Purposes
Process automation of point cloud data derived from terrestrial-based technologies such as Terrestrial Laser Scanner (TLS) or Mobile Laser Scanner (MLS). 'FORTLS' enables (i) detection of trees and estimation of tree-level attributes (e.g. diameters and heights), (ii) estimation of stand-level variables (e.g. density, basal area, mean and dominant height), (iii) computation of metrics related to important forest attributes estimated in Forest Inventories (FIs) at stand-level, and (iv) optimization of plot design for combining TLS data and field measured data. Documentation about 'FORTLS' is described in Molina-Valero et al. (2022, <https://doi.org/10.1016/j.envsoft.2022.105337>).

Get the lat stable version of FORTLS from GitHub (included in the master branch)

```r
remotes::install_github("Molina-Valero/FORTLS", ref = "devel", dependencies = TRUE)
```

# Taller de manejo de nubes de puntos forestales - 9CFE

## Instalación de FORTLS

```r
install.packages(“FORTLS”)
library(FORTLS)
```

## Establecimiento del directorio de trabajo

Por ejemplo: "C:\taller_FORTLS"

```r
setwd("C:/taller_FORTLS")

```

## Normalización de la nube de puntos

```r
pcd <- normalize(las = "HLS_LiGrip.laz",
                 id = "HLS_LiGrip",
                 max.dist = 12.5,
                 scan.approach = "multi")
```

## Variables de árbol individual (o dendrométricas)

```r
tree.tls <- tree.detection.multi.scan(data = pcd,
                                      understory = TRUE)
```

## Variables de masa (o dasométricas)

```r
met.var.TLS <- metrics.variables(tree.tls = tree.tls,
                                 scan.approach = "multi",
                                 plot.parameters = data.frame(radius = 10, k = 10, BAF = 2))
                                 
# Parcela circular de área fija (10 m de radio)
parcela.circular <- met.var.TLS$fixed.area

# Parceka k-tree (k = 10)
parcela.k.tree <- met.var.TLS$k.tree

# Parcela relascópica (BAF = 2)
parcela.relascopica <- met.var.TLS$angle.count
```

# Acknowledgements 

**FORTLS** it is being developed at [Czech University of Life Sciences Prague](https://www.czu.cz/en) and [University of Santiago de Compostela](https://www.usc.gal/en).

Development of the `FORTLS` package is being possible thanks to the following fellowships/projects:

* Climate Change Adaptation of Forests in the Brdy Highland [LIFE21-CCA-CZ-LIFE-Adapt-Brdy/101074426](https://adaptbrdy.czu.cz/en)
* Design of forest monitoring systems on a regional scale [ED431F 2020/02] supported by the Regional Government of Galicia
* Ramón Areces Foundation Grants for Postdoctoral Studies [XXXV Call for Expansion of Studies Abroad in Life and Matter Sciences](https://www.fundacionareces.es/fundacionareces/es/tratarAplicacionInvestigador.do?paginaActual=2&idConvocatoria=2770&tipo=2)

<p align = "center">
  <img src="https://www.dropbox.com/scl/fi/cy3cfikrgwl54eovz3ncn/CZU_logotype_V_ENG_green.png?rlkey=hbbggghvn93412oqa85m0fpm0&raw=1" height="100"> 
  <img src="https://www.dropbox.com/scl/fi/g7dyqq5yzzvg2vu2dk6jv/usc.png?rlkey=z3x7mwx1ebsioivrwg9fpgdyq&raw=1" height="50" hspace="10"> 
  <img src="https://www.dropbox.com/scl/fi/9ohh7hs6sd9imxzsfb768/ccefpu-positivo.png?rlkey=g71a5x4qejmyybpwwc2vavfzg&raw=1" height="50"> 
  <img src="https://www.dropbox.com/scl/fi/zk0ktudsu0caszlw2z3dm/logotipo-fra-color.jpg?rlkey=1fiee4ra7mm98pdrozirlyo1p&raw=1" height="100">
</p>

<img src="https://www.dropbox.com/scl/fi/ec1m3266bcoq8qrgjqykv/logolink-RGB_LAB-LIFE-NATURA-MZP_en-okraje-1000x156.jpg?rlkey=bi5018o95zq63rwhhwa6wfs2y&raw=1" align="center">

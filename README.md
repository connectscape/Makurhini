
<!-- README.md is generated from README.Rmd. Please edit that file -->
Makurhini
=========

<!-- badges: start -->
<!-- badges: end -->
The goal of Makurhini is to provide a set of functions to estimate landscape fragmentation and connectivity metrics

Installation
------------

You can install the released version of Makurhini from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("Makurhini")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
library(sf)
#> Warning: package 'sf' was built under R version 3.5.1
#> Linking to GEOS 3.6.1, GDAL 2.2.3, PROJ 4.9.3
## basic example code
```

``` r
library(Makurhini)
ruta <- system.file("extdata", "Habitat_Patches.shp", package = "Makurhini")
cores <- read_sf(ruta)
cores
#> Simple feature collection with 143 features and 2 fields
#> geometry type:  POLYGON
#> dimension:      XY
#> bbox:           xmin: 3246137 ymin: 322869.6 xmax: 3739484 ymax: 696540.5
#> epsg (SRID):    NA
#> proj4string:    +proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +datum=WGS84 +units=m +no_defs
#> # A tibble: 143 x 3
#>       id   Area                                                    geometry
#>    <int>  <dbl>                                               <POLYGON [m]>
#>  1  1347 4196.  ((3676911 589967.3, 3676931 589895.5, 3676948 589818.3, 36~
#>  2  1609   60.2 ((3558044 696202.5, 3557972 696280.9, 3557957 696291.6, 35~
#>  3  1632   48.9 ((3569169 687776.4, 3569146 687749.5, 3569096 687745.6, 35~
#>  4  1643   15.2 ((3547317 685713.2, 3547363 685573.9, 3547391 685401.4, 35~
#>  5  1650   33.3 ((3567471 684357.4, 3567380 684214.3, 3567303 684046, 3567~
#>  6  1659   53.1 ((3590569 672451.7, 3590090 672574.9, 3589912 672547.5, 35~
#>  7  1660   83.8 ((3570789 670959.4, 3570860 671015.4, 3570909 671019.3, 35~
#>  8  1661   17.4 ((3440118 666273.2, 3440372 666849.2, 3440584 667001.7, 34~
#>  9  1662   18.0 ((3451637 671232.4, 3451616 671287.1, 3451535 671315.1, 34~
#> 10  1663   36.3 ((3444396 671675.7, 3444715 671834.8, 3444873 672019, 3444~
#> # ... with 133 more rows
```

<img src="man/figures/README-cores-1.png" width="100%" />

``` r
nrow(cores)
#> [1] 143
IIC <- dConnectivity(nodes = cores, id = "id", attribute = "Area",
                     distance = list(type = "centroid"),
                     metric = "IIC", distance_thresholds = 30000,
                     LA = NULL, overall = FALSE)
plot(IIC["dIIC"], breaks = "jenks")
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

``` r
plot(IIC["dIICintra"], breaks = "jenks")
```

<img src="man/figures/README-unnamed-chunk-3-2.png" width="100%" />

``` r
plot(IIC["dIICflux"], breaks = "jenks")
```

<img src="man/figures/README-unnamed-chunk-3-3.png" width="100%" />

``` r
plot(IIC["dIICconnector"], breaks = "jenks")
```

<img src="man/figures/README-unnamed-chunk-3-4.png" width="100%" />

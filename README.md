
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/OscarGOGO/Makurhini?branch=master&svg=true)](https://ci.appveyor.com/project/OscarGOGO/Makurhini)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Makurhini

![LOGO\_MAKHURINI](https://user-images.githubusercontent.com/30152793/79033305-ce8c2280-7b72-11ea-8df7-b8b48409b818.png)

The goal of Makurhini is to provide a set of functions to estimate
landscape fragmentation and connectivity metrics

## Installation

You can install the released version of Makurhini from
[GitHub](https://github.com) with:

``` r
library(devtools)
install_github("OscarGOGO/Makurhini", dependencies = TRUE)
```

## Example

This is a basic example which shows you how to solve some common
problems:

### Protected Connected Land (ProtConn)

Protected areas:

``` r
data("Protected_areas", package = "Makurhini")
data("regions", package = "Makurhini")
region <- regions[1,]
```

<img src="man/figures/README-cores-1.png" width="100%" /><img src="man/figures/README-cores-2.png" width="100%" />

``` r
test <- MK_ProtConn(nodes = Protected_areas, region = region,
                    attribute = "Intersected area", area_unit = "ha",
                    distance = list(type= "centroid"),
                    distance_thresholds = 10000,
                    probability = 0.5, transboundary = 50000,
                    LA = NULL, plot = TRUE, dPC = FALSE,
                    write = NULL, SAGA = FALSE, intern = FALSE)
test$`Protected Connected (Viewer Panel)`
```

<table class="table table-condensed">

<thead>

<tr>

<th style="text-align:left;">

Index

</th>

<th style="text-align:center;">

Value

</th>

<th style="text-align:left;">

ProtConn indicator

</th>

<th style="text-align:center;">

Percentage

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold">EC(PC)</span>

</td>

<td style="text-align:center;">

130189.18

</td>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold">Unprotected </span>

</td>

<td style="text-align:center;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffa90d">92.540</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold">PC </span>

</td>

<td style="text-align:center;">

1.2324e-03

</td>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold">Prot </span>

</td>

<td style="text-align:center;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff8eb">7.460</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold"> </span>

</td>

<td style="text-align:center;">

</td>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold">ProtConn </span>

</td>

<td style="text-align:center;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffbf5">3.511</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold"> </span>

</td>

<td style="text-align:center;">

</td>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold">ProtUnconn </span>

</td>

<td style="text-align:center;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffbf4">3.950</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold"> </span>

</td>

<td style="text-align:center;">

</td>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold">RelConn </span>

</td>

<td style="text-align:center;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffd383">47.058</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold"> </span>

</td>

<td style="text-align:center;">

</td>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold">ProConn\_design </span>

</td>

<td style="text-align:center;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffbf4">3.950</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold"> </span>

</td>

<td style="text-align:center;">

</td>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold">ProConn\_Bound </span>

</td>

<td style="text-align:center;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffbf5">3.511</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold"> </span>

</td>

<td style="text-align:center;">

</td>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold">ProtConn\_Prot </span>

</td>

<td style="text-align:center;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffa500">97.512</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold"> </span>

</td>

<td style="text-align:center;">

</td>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold">ProtConn\_Trans </span>

</td>

<td style="text-align:center;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">0.000</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold"> </span>

</td>

<td style="text-align:center;">

</td>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold">ProtConn\_Unprot </span>

</td>

<td style="text-align:center;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf8">2.488</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold"> </span>

</td>

<td style="text-align:center;">

</td>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold">ProtConn\_Within </span>

</td>

<td style="text-align:center;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffa707">94.784</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold"> </span>

</td>

<td style="text-align:center;">

</td>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold">ProtConn\_Contig </span>

</td>

<td style="text-align:center;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf7">2.728</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold"> </span>

</td>

<td style="text-align:center;">

</td>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold">ProtConn\_Within\_land</span>

</td>

<td style="text-align:center;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffbf6">3.327</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold"> </span>

</td>

<td style="text-align:center;">

</td>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold">ProtConn\_Contig\_land</span>

</td>

<td style="text-align:center;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffefe">0.096</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold"> </span>

</td>

<td style="text-align:center;">

</td>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold">ProtConn\_Unprot\_land</span>

</td>

<td style="text-align:center;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffefe">0.087</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold"> </span>

</td>

<td style="text-align:center;">

</td>

<td style="text-align:left;">

<span style="color: #636363; font-weight: bold">ProtConn\_Trans\_land
</span>

</td>

<td style="text-align:center;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">0.000</span>

</td>

</tr>

</tbody>

</table>

``` r
test$`ProtConn Plot`
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

### Equivalent Connectivity (EC)

Example with old-growth vegetation fragments of four times
(?list\_forest\_patches).

``` r
data("list_forest_patches", package = "Makurhini")
data("study_area", package = "Makurhini")

Max_attribute <- unit_convert(gArea(study_area), "m2", "ha")
```

``` r
dECA_test <- MK_dECA(nodes= list_forest_patches, attribute = NULL, area_unit = "ha",
                  distance = list(type= "centroid"), metric = "PC",
                  probability = 0.05, distance_thresholds = 5000,
                  LA = Max_attribute, plot= c("1993", "2003", "2007", "2011"))
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

``` r
dECA_test
```

<table class="table table-condensed">

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:left;">

Scenary

</th>

<th style="text-align:right;">

Area (ha)

</th>

<th style="text-align:right;">

ECA (ha)

</th>

<th style="text-align:right;">

Distance

</th>

<th style="text-align:right;">

Normalized ECA

</th>

<th style="text-align:right;">

dA

</th>

<th style="text-align:right;">

dECA

</th>

<th style="text-align:right;">

dA/dECA comparisons

</th>

<th style="text-align:right;">

Type of change

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

T1

</td>

<td style="text-align:left;">

1

</td>

<td style="text-align:right;">

<span style="display: inline-block; direction: rtl; border-radius: 4px; padding-right: 2px; background-color: #94D8B1; width: 98.07%">91438.11</span>

</td>

<td style="text-align:right;">

50657.60

</td>

<td style="text-align:right;">

5000

</td>

<td style="text-align:right;">

<span>55.40%</span>

</td>

<td style="text-align:right;">

<span style="color: red">-67.246</span>

</td>

<td style="text-align:right;">

<span style="color: red">-81.854</span>

</td>

<td style="text-align:right;">

dECA \< dA \< 0

</td>

<td style="text-align:right;">

  - Connectivity loss
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    T2
    </td>
    <td style="text-align:left;">
    2
    </td>
    <td style="text-align:right;">
    <span style="display: inline-block; direction: rtl; border-radius: 4px; padding-right: 2px; background-color: #94D8B1; width: 100.00%">93238.91</span>
    </td>
    <td style="text-align:right;">
    53604.33
    </td>
    <td style="text-align:right;">
    5000
    </td>
    <td style="text-align:right;">
    <span>57.49%</span>
    </td>
    <td style="text-align:right;">
    <span style="color: green">1.969</span>
    </td>
    <td style="text-align:right;">
    <span style="color: green">5.817</span>
    </td>
    <td style="text-align:right;">
    dECA or dA gain
    </td>
    <td style="text-align:right;">
    Habitat or connectivity gain
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    T3
    </td>
    <td style="text-align:left;">
    3
    </td>
    <td style="text-align:right;">
    <span style="display: inline-block; direction: rtl; border-radius: 4px; padding-right: 2px; background-color: #94D8B1; width: 89.57%">83517.49</span>
    </td>
    <td style="text-align:right;">
    38756.64
    </td>
    <td style="text-align:right;">
    5000
    </td>
    <td style="text-align:right;">
    <span>46.41%</span>
    </td>
    <td style="text-align:right;">
    <span style="color: red">-10.426</span>
    </td>
    <td style="text-align:right;">
    <span style="color: red">-27.699</span>
    </td>
    <td style="text-align:right;">
    dECA \< dA \< 0
    </td>
    <td style="text-align:right;">
      - Connectivity loss
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">
        T4
        </td>
        <td style="text-align:left;">
        4
        </td>
        <td style="text-align:right;">
        <span style="display: inline-block; direction: rtl; border-radius: 4px; padding-right: 2px; background-color: #94D8B1; width: 89.94%">83859.71</span>
        </td>
        <td style="text-align:right;">
        40187.05
        </td>
        <td style="text-align:right;">
        5000
        </td>
        <td style="text-align:right;">
        <span>47.92%</span>
        </td>
        <td style="text-align:right;">
        <span style="color: green">0.410</span>
        </td>
        <td style="text-align:right;">
        <span style="color: green">3.691</span>
        </td>
        <td style="text-align:right;">
        dECA or dA gain
        </td>
        <td style="text-align:right;">
        Habitat or connectivity gain
        </td>
        </tr>
        </tbody>
        </table>

### Citing Makurhini package

A formal paper detailing this packe is forthcoming, but until it is
published, please use the something like the following to cite if you
use it in your work:

<code> <i> Godínez-Gómez, O. and Correa Ayram C.A. 2020. Makurhini: An R
package for analyzing landscape connectivity.
<https://github.com/OscarGOGO/Makurhini>, DOI: 10.5281/zenodo.3748095
</i> </code>

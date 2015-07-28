<!-- README.md is generated from README.Rmd. Please edit that file -->
The Tissot Indicatrix
=====================

The Tissot Indicatrix is used to characterize local distortions within map projections.

<https://en.wikipedia.org/wiki/Tissot%27s_indicatrix>

NOTE: this package is in early development and will change a lot, this is just a quick prototype and needs quite a bit more to be useful.

TODO
----

-   sort out the interface for inputting data, projections
-   semi-automate the scaling for the plot
-   provide control and better defaults for colours
-   check terminology (and spelling) is consistent

Installation
============

Can be installed with

``` r
devtools::install_github("mdsumner/tissot")
```

Minimal example
===============

``` r
library(tissot)
 library(rgdal)
 prj <- function(z, proj.in, proj.out) {
   z.pt <- SpatialPoints(coords=matrix(z, ncol=2), proj4string=proj.in)
   w.pt <- spTransform(z.pt, CRS=proj.out)
   w.pt@coords[1, ]
 }
 # Longitude, latitude, and reprojection function
 # NAD 27 in
 # World Robinson projection out
 r <- tissot(130, 54, prj,
             proj.in=CRS("+init=epsg:4267"),
             proj.out=CRS("+init=esri:54030"))

 i <- indicatrix(r, scale=10^4, n=71)
 plot(i$outline, type="n", asp=1, xlab="Easting", ylab="Northing")
 polygon(i$base, col=rgb(0, 0, 0, .025), border="Gray")
 lines(i$d.lambda, lwd=2, col="Gray", lty=2)
 lines(i$d.phi, lwd=2, col=rgb(.25, .7, .25), lty=2)
 lines(i$axis.major, lwd=2, col=rgb(.25, .25, .7))
 lines(i$axis.minor, lwd=2, col=rgb(.7, .25, .25))
 lines(i$outline, asp=1, lwd=2)
```

![](readmefigs/README-unnamed-chunk-3-1.png)

Derived from

<http://gis.stackexchange.com/questions/31651/an-example-tissot-ellipse-for-an-equirectangular-projection>

Polar example
=============

``` r
 library(tissot)
library(rgdal)
library(maptools)
#> Checking rgeos availability: TRUE
data(wrld_simpl)
prj <- function(z, proj.in, proj.out) {
  z.pt <- SpatialPoints(coords=matrix(z, ncol=2), proj4string=proj.in)
  w.pt <- spTransform(z.pt, CRS=proj.out)
  w.pt@coords[1, ]
}
library(raster)

buildandplot <- function(data, title = "", scale = 5e5, nr = 7, nc = 7) {
  ## grid of points
  gr <- rasterToPoints(raster(data, nrow = nr, ncol = nc), spatial = TRUE)
  grll <- spTransform(gr, CRS("+proj=longlat +ellps=WGS84"))
  tis <- vector("list", length(gr))
  for (i in seq_along(tis)) tis[[i]] <- tissot(coordinates(grll)[i, 1], coordinates(grll)[i, 2], prj,  
                                               proj.in = CRS(projection(wrld_simpl)), proj.out = projection(data))
plot(data, main = title)
for (j in seq_along(tis)) {
  indi <- indicatrix(tis[[j]], scale = scale, n = 71)
  plot(indi, add = TRUE)
}
invisible(NULL)
}
## choose a projection
ptarget1 <- "+proj=stere +lon_0=147 +lat_ts-71 +lat_0=-90 +ellps=WGS84"
w1 <- spTransform(subset(wrld_simpl, coordinates(wrld_simpl)[,2] < 10), CRS(ptarget1))

ptarget2 <- "+proj=laea +lon_0=147 +lat_0=-90 +ellps=WGS84"
w2 <- spTransform(subset(wrld_simpl, coordinates(wrld_simpl)[,2] < 10), CRS(ptarget2))

ptarget3 <- "+proj=omerc +lonc=147 +gamma=9 +alpha=9 +lat_0=-80 +ellps=WGS84"
w3 <- spTransform(subset(wrld_simpl, coordinates(wrld_simpl)[,2] < -12), CRS(ptarget3), scale = 3e5)

buildandplot(w1, "Polar Stereographic")
```

![](readmefigs/README-unnamed-chunk-4-1.png)

``` r
buildandplot(w2, "Lambert Azimuthal Equal Area")
```

![](readmefigs/README-unnamed-chunk-4-2.png)

``` r
buildandplot(w3, "Oblique Mercator")
```

![](readmefigs/README-unnamed-chunk-4-3.png)

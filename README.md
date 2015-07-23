<!-- README.md is generated from README.Rmd. Please edit that file -->
The Tissot Indicatrix
=====================

Can be installed with

``` r
devtools::install_github("mdsumner/tissot")
```

Minimal example

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
#> Note: no visible binding for global variable '.Data'

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



## choose a projection
ptarget <- "+proj=stere +lat_ts-71 +lat_0=-90 +ellps=WGS84"
w <- spTransform(subset(wrld_simpl, coordinates(wrld_simpl)[,2] < 0), CRS(ptarget))
#> Note: no visible binding for global variable 'plotOrder' 
#> Note: no visible binding for global variable 'plotOrder' 
#> Note: no visible binding for global variable 'plotOrder'

## grid of points
library(raster)
#> Note: no visible binding for global variable 'plotOrder' 
#> Note: no visible binding for global variable 'coords' 
#> Note: no visible binding for global variable 'coords.nrs' 
#> Note: no visible binding for global variable 'coords' 
#> Note: no visible binding for global variable 'coords'
gr <- rasterToPoints(raster(w, nrow = 7, ncol = 7), spatial = TRUE)
grll <- spTransform(gr, CRS(projection(wrld_simpl)))
tis <- vector("list", length(gr))
for (i in seq_along(tis)) tis[[i]] <- tissot(coordinates(grll)[i, 1], coordinates(grll)[i, 2], prj,  proj.in = CRS(projection(wrld_simpl)), proj.out = projection(w))
plot(w)
for (j in seq_along(tis)) {
  i <- indicatrix(tis[[j]], scale = 3e5, n = 71)
  polygon(i$base, col=rgb(0, 0, 0, .025), border="Gray")
  lines(i$d.lambda, lwd=2, col="Gray", lty=2)
  lines(i$d.phi, lwd=2, col=rgb(.25, .7, .25), lty=2)
  lines(i$axis.major, lwd=2, col=rgb(.25, .25, .7))
  lines(i$axis.minor, lwd=2, col=rgb(.7, .25, .25))
  lines(i$outline, asp=1, lwd=1)
  
}
```

![](readmefigs/README-unnamed-chunk-4-1.png)

<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->
[![R-CMD-check](https://github.com/hypertidy/tissot/workflows/R-CMD-check/badge.svg)](https://github.com/hypertidy/tissot/actions)
<!-- badges: end -->

# The Tissot Indicatrix

The [Tissot
Indicatrix](https://en.wikipedia.org/wiki/Tissot%27s_indicatrix) is used
to characterize local distortions within map projections.

I have derived the code in this package (with permission) from Bill
Huber’s wonderful online answer here:

<http://gis.stackexchange.com/questions/31651/an-example-tissot-ellipse-for-an-equirectangular-projection>

Also see

<https://gis.stackexchange.com/questions/5068/how-to-create-an-accurate-tissot-indicatrix>

# Installation

Can be installed with

``` r
remotes::install_github("hypertidy/tissot")
```

# Minimal example

``` r
library(tissot)
# NAD 27 in
# World Robinson projection out
r <- tissot(130, 54,
           proj.in= "EPSG:4267",  
           proj.out= "ESRI:54030")
i <- indicatrix(r, scale=10^4, n=71)
plot(i)
```

![](readmefigs/README-minimal-1.png)

Since an original port of whuber’s code we have now made it much easier
to create many indicatrixes and plot them in one step. Or we can still
just grab one and plot it on its own. Note that the scale is quite
different in these plots.

``` r
x <- seq(-175, 175, by = 20)
y <- seq(-82.5, 82.5, by = 15)
xy <- expand.grid(x, y)
r <- tissot(xy,
            proj.in= "OGC:CRS84",
            proj.out= "+proj=robin")

i <- indicatrix0(r[1, ], scale=10^4, n=71)
plot(i, add = FALSE)
```

![](readmefigs/README-bigger-example-1.png)

``` r
ii <- indicatrix(r, scale=3e5, n=71)
plot(ii, add = FALSE)
```

![](readmefigs/README-bigger-example-2.png)

Mollweide.

``` r
m <- tissot(xy,
            proj.in= "OGC:CRS84",
            proj.out= "+proj=moll")


plot(indicatrix(m, scale=3e5, n=71), add = FALSE)
```

![](readmefigs/README-mollweide-1.png)

Eckhert I

``` r
e <- tissot(xy,
            proj.in= "OGC:CRS84",
            proj.out= "+proj=eck1")


plot(indicatrix(e, scale=3e5, n=71), add = FALSE)
```

![](readmefigs/README-eckhert-1.png)

Sinusoidal

``` r
s <- tissot(xy,
            proj.in= "OGC:CRS84",
            proj.out= "+proj=sinu")


plot(indicatrix(s, scale=3e5, n=71), add = FALSE)
```

![](readmefigs/README-sinu-1.png)

# Polar example

``` r
p <- tissot(xy[xy[,2] < -40, ],
            proj.in= "OGC:CRS84",
            proj.out= "+proj=stere +lon_0=147 +lat_ts-71 +lat_0=-90 +datum=WGS84")

plot(indicatrix(p, scale = 3e5))
```

![](readmefigs/README-polar-stereo-1.png)

``` r
laea <- tissot(xy[xy[,2] < 20, ],
            proj.in= "OGC:CRS84",
            proj.out= "+proj=laea +lon_0=147 +lat_0=-90 +datum=WGS84")

plot(indicatrix(laea, scale = 3e5))
```

![](readmefigs/README-polar-laea-1.png)

Oblique Mercator

You get the idea … many projections need extra attention for real data.

``` r
mp0 <- do.call(cbind, maps::map(plot = FALSE)[1:2])
omerc <- "+proj=omerc +lonc=147 +gamma=9 +alpha=9 +lat_0=-80 +ellps=WGS84"
mp <- tissot:::.prj(mp0, proj.out = omerc, proj.in = "OGC:CRS84")
o <- tissot(xy,
            proj.in= "OGC:CRS84",
            proj.out= omerc)

plot(indicatrix(o, scale = 3e5))
lines(mp)
```

![](readmefigs/README-omerc-1.png)

VicGrid

``` r
vgrid <- "+proj=lcc +lat_1=-36 +lat_2=-38 +lat_0=-37 +lon_0=145 +x_0=2500000 +y_0=2500000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
mp <- tissot:::.prj(mp0, 
                    proj.out = vgrid, proj.in = "OGC:CRS84")
v <- tissot(as.matrix(expand.grid(seq(120, 165, by =5 ), 
                                          seq(-45, -35, by = 5))),
            proj.in= "OGC:CRS84",
            proj.out= vgrid)

plot(indicatrix(v, scale = 2e5))
lines(mp)
```

![](readmefigs/README-vicgrid-1.png)

UTM Zone 54 (Hobart)

``` r
utm <- "+proj=utm +zone=54 +south"
mp <- tissot:::.prj(mp0, 
                    proj.out = utm, proj.in = "OGC:CRS84")
u <- tissot(as.matrix(expand.grid(seq(108, 162, by =6 ), 
                                          seq(-65, 55, by = 15))),
            proj.in= "OGC:CRS84",
            proj.out= utm)

plot(indicatrix(u, scale = 2e5))
lines(mp)
```

![](readmefigs/README-utm54-1.png)

## Code of Conduct

Please note that the tissot project is released with a [Contributor Code
of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

.prj <- function(z, proj.out, ..., proj.in = NULL) {
  if (is.null(proj.in)) {
    message("assuming WGS84 for unprojected angular coordinates")
    proj.in <- "OGC:CRS84"
  }
  capture.output(out <- PROJ::proj_trans(matrix(z, ncol = 2), proj.out, source = proj.in))
   do.call(cbind, out)
}
.to.degrees <- function(x) x * 180 / pi
.to.radians <- function(x) x * pi / 180
.clamp <- function(x) min(max(x, -1), 1)                             # Avoids invalid args to asin
.norm <- function(x) sqrt(sum(x*x))


#' Tissot
#'
#' Create the Tissot Indicatrix.
#'
#' Compute properties of scale distortion and Tissot's indicatrix at location x = c(lambda, phi)
#' using prj as the projection.  A is the ellipsoidal semi-major axis (in meters) and f.inv is
#' the inverse flattening.  The projection must return a vector (x, y) when given a vector (lambda, phi).
#' (Not vectorized.)  Optional arguments ... are passed to prj.
#' Source: Snyder pp 20-26 (WGS 84 defaults for the ellipsoidal parameters).
#' All input and output angles are in degrees.
#' @param lambda longitude
#' @param phi latitude
#' @param degrees logical, work in degrees or radians
#' @param A ellipsoidal semi-major axis (meters)
#' @param f.inv the inverse flattening
#' @param ... passed to internal function
#' @param proj.in projection of input
#' @param proj.out projection of context
#' @return list with stuff as per below
#' @export
#' @examples
#' x <- seq(-175, 175, by = 15)
#' y <- seq(-82.5, 82.5, by = 15)
#' xy <- expand.grid(x, y)
#' r <- tissot(xy,
#'             proj.in= "OGC:CRS84",
#'             proj.out= "+proj=robin")
#'
#' j <- which.min(abs(135 - r$lon) + abs(54 - r$lat))
#' idx0 <- indicatrix0(r[j, ], scale=10^4, n=71)
#' op <- par(mfrow =   c(2, 1))
#' plot(idx0, add = FALSE)

#' idx <- indicatrix(r, scale=3e5, n=71)
#' plot(idx, add = FALSE)
#' tissot_abline(r$lon[j], r$lat[j])
#' par(op)
#'
#'
#' @importFrom grDevices grey rgb
#' @importFrom graphics lines plot polygon
#' @importFrom stats numericDeriv
#' @importFrom grDevices xy.coords
tissot <- function(lambda, phi = NULL,  degrees=TRUE, A = 6378137, f.inv=298.257223563, ..., proj.in, proj.out) {
 xy <- xy.coords(lambda, phi)
 lam <- xy[[1L]]
 phi <- xy[[2L]]
 out <- vector("list", length(lam))
 for (i in seq_along(lam)) {
   out[[i]] <- tissot0(lam[i], phi[i], degrees = degrees, A = A, f.inv = f.inv, proj.in = proj.in, proj.out = proj.out, ...)
 }
 do.call(rbind, out)
}
#' @importFrom tibble tibble
tissot0 <- function(lambda, phi,  degrees=TRUE, A = 6378137, f.inv=298.257223563, proj.out = NULL,
                    ..., proj.in = NULL) {
  if (is.null(proj.in)) {
    proj.in <- "OGC:CRS84"
  }
  if (is.null(proj.out)) {
    stop("'proj.out' must be specified")
  }
  #
  # Precomputation.
  #
  if (f.inv==0) {                                                     # Use f.inv==0 to indicate a sphere
    e2 <- 0
  } else {
    e2 <- (2 - 1/f.inv) / f.inv                                       # Squared eccentricity
  }
  if (degrees) phi.r <- .to.radians(phi) else phi.r <- phi
  cos.phi <- cos(phi.r)                                               # Convenience term
  e2.sinphi <- 1 - e2 * sin(phi.r)^2                                  # Convenience term
  e2.sinphi2 <- sqrt(e2.sinphi)                                       # Convenience term
  if (degrees) units <- 180 / pi else units <- 1                    # Angle measurement units per radian
  #
  # Lengths (the metric).
  #
  radius.meridian <- A * (1 - e2) / e2.sinphi2^3                      # (4-18)
  length.meridian <- radius.meridian                                  # (4-19)
  radius.normal <- A / e2.sinphi2                                     # (4-20)
  length.normal <- radius.normal * cos.phi                            # (4-21)
  #
  # The projection and its first derivatives, normalized to unit lengths.
  #
  x <- c(lambda, phi)
  d <- numericDeriv(quote(.prj(x, proj.out = proj.out, proj.in = proj.in)[1L, ]), theta="x")
  z <- d[1:2]                                                         # Projected coordinates
  names(z) <- c("x", "y")
  g <- attr(d, "gradient")                                            # First derivative (matrix)
  g <- g %*% diag(units / c(length.normal, length.meridian))          # Unit derivatives
  dimnames(g) <- list(c("x", "y"), c("lambda", "phi"))
  g.det <- det(g)                                                     # Equivalent to (4-15)
  #
  # Computation.
  #
  h <- .norm(g[, "phi"])                                               # (4-27)
  k <- .norm(g[, "lambda"])                                            # (4-28)
  a.p <- sqrt(max(0, h^2 + k^2 + 2 * g.det))                          # (4-12) (intermediate)
  b.p <- sqrt(max(0, h^2 + k^2 - 2 * g.det))                          # (4-13) (intermediate)
  a <- (a.p + b.p)/2                                                  # (4-12a)
  b <- (a.p - b.p)/2                                                  # (4-13a)
  omega <- 2 * asin(.clamp(b.p / a.p))                                 # (4-1a)
  theta.p <- asin(.clamp(g.det / (h * k)))                             # (4-14)
  conv <- (atan2(g["y", "phi"], g["x","phi"]) + pi / 2) %% (2 * pi) - pi # Middle of p. 21
  #
  # The indicatrix itself.
  # "svd" essentially redoes the preceding computation of "h", "k", and "theta.p".
  #
  m <- svd(g)
  axes <- zapsmall(diag(m$d) %*% apply(m$v, 1, function(x) x / .norm(x)))
  dimnames(axes) <- list(c("major", "minor"), NULL)

  # list(location=c(lambda, phi), projected=z,
  #             meridian_radius=radius.meridian, meridian_length=length.meridian,
  #             normal_radius=radius.normal, normal_length=length.normal,
  #             scale.meridian=h, scale.parallel=k, scale.area=g.det, max.scale=a, min.scale=b,
  #             .to.degrees(zapsmall(c(angle_deformation=omega, convergence=conv, intersection_angle=theta.p))),
  #             axes=axes, derivatives=g)

  dfconint <- .to.degrees(zapsmall(c(omega, conv, theta.p)))
  tibble::tibble(lon = lambda, lat = phi, x = z[1L], y = z[2L],
                 meridian_radius=radius.meridian, meridian_length=length.meridian,
                 normal_radius=radius.normal, normal_length=length.normal,
                 scale.meridian=h, scale.parallel=k, scale.area=g.det, max.scale=a, min.scale=b,
                 angle_deformation= dfconint[1L], convergence= dfconint[2L], intersection_angle= dfconint[3L],
                 axes_x_major = axes[1L,1L, drop = TRUE], axes_x_minor = axes[2L,1L, drop = TRUE],
                 axes_y_major = axes[1L,2L, drop = TRUE], axes_y_minor = axes[2L,2L, drop = TRUE],
                 lambda_dx = g[1L, 1L, drop = TRUE], lambda_dy = g[2L, 1L, drop = TRUE],
                 phi_dx = g[1L, 2L, drop = TRUE], phi_dy = g[2L, 2L, drop = TRUE],
                 proj.in = proj.in, proj.out = proj.out)

}

#' Indicatrix
#'
#' Indicatrix
#'
#' Reprocesses the output of \code{tissot} into convenient geometrical data.
#' @param x object from \code{tissot}
#'
#' @param scale scaling
#' @param ... arguments \code{n}, \code{from} and \code{to} passed to \code{ti_ellipse} function
#'
#' @export
#'
indicatrix <- function(x, scale = 1, ...) {
  structure(lapply(split(x, 1:nrow(x)), indicatrix0, scale = scale, ...),
            class = c("indicatrixes", "list"))
}

#' @name indicatrix
#' @export
plot.indicatrixes <- function(x, asp=1, xlab="x", ylab="y", add = FALSE, ...,
 col.base = rgb(0, 0, 0, .1),
 col.lambda = grey(0.75),
 #col.phi =  "#45A271", ## #5CBD92", #rgb(.25, .7, .25),  ## green
 #col.major = "#A782C3", ##"#7DB0DD", #rgb(.25, .25, .7), ## purple
 #col.minor = "#C87A8A", ##"#C7A76C", #rgb(.7, .25, .25), ## red
 col.phi = "#1b9e77",  ## colorbrewer, qualitative, color friendly
 col.major = "#7570b3",
 col.minor = "#d95f02",
 col.outline = "black") {
 if (!add) {
   center <- do.call(rbind, lapply(x, function(a) a$center))
   major <- do.call(rbind, lapply(x, function(a) a$axis.major))
   minor <- do.call(rbind, lapply(x, function(a) a$axis.minor))
   props <- list(x = center, type = "n", asp = asp, xlab = xlab, ylab = ylab, ...)
   if (is.null(props$xlim)) props$xlim <-range(c(major[,1L], minor[,1L]))
   if (is.null(props$ylim)) props$ylim <- range(c(major[,2L], minor[,2L]))

   do.call(plot, props)
   options(tissot.last.plot.proj = x[[1]]$proj.out)
 }
  lapply(x, plot, add = TRUE,
         col.base = col.base,
         col.lambda = col.lambda,
         col.phi = col.phi,
         col.major = col.major,
         col.minor = col.minor,
         col.outline = col.outline
         )
  invisible(NULL)
}
#' @name indicatrix
#' @export
indicatrix0 <- function(x, scale=1, ...) {
  o <- unlist(x[c("x", "y")])
  base <- ti_ellipse(o, matrix(c(1,0,0,1), 2), scale=scale, ...)             # A reference circle
  axes_major <- unlist(x[c("axes_x_major", "axes_y_major")])
  axis.major <- rbind(o + scale * axes_major, o - scale * axes_major)
  axes_minor <- unlist(x[c("axes_x_minor", "axes_y_minor")])
  axis.minor <- rbind(o + scale * axes_minor, o - scale * axes_minor)
  outline <- ti_ellipse(o, matrix(c(axes_major, axes_minor), 2L, byrow = T), scale=scale, ...)
  #dimnames(g) <- list(c("x", "y"), c("lambda", "phi"))
  lambda_d <- unlist(x[c("lambda_dx", "lambda_dy")])
  d.lambda <- rbind(o + scale * lambda_d, o - scale * lambda_d)
  phi_d <- unlist(x[c("phi_dx", "phi_dy")])
  d.phi <- rbind(o + scale * phi_d, o - scale * phi_d)

  i <- list(center=o, base = base, outline = outline,
              axis.major = axis.major, axis.minor = axis.minor,
              d.lambda = d.lambda, d.phi = d.phi,
            proj.in = x$proj.in, proj.out = x$proj.out)
  class(i) <- c("indicatrix0", "list")
  i
}


#' Methods for indicatrix
#'
#' plot indicatrix
#'
#' @param asp aspect ratio
#' @param xlab x-axis labels
#' @param ylab y-axis labels
#' @param add add to existing plot
#' @param col.base colour of base
#' @param col.lambda colour of lambda
#' @param col.phi colour of phi
#' @param col.major major axis colour
#' @param col.minor minor axis colour
#' @param col.outline outline colour
#'
#' @rdname indicatrix
#' @export
plot.indicatrix0 <- function(x, asp=1, xlab="Easting", ylab="Northing", add = FALSE, ...,
                            col.base = rgb(0, 0, 0, .1),
                            col.lambda = grey(0.75),
                            col.phi =  "#45A271", ## #5CBD92", #rgb(.25, .7, .25),
                            col.major = "#A782C3", ##"#7DB0DD", #rgb(.25, .25, .7),
                            col.minor = "#C87A8A", ##"#C7A76C", #rgb(.7, .25, .25),
                            col.outline = "black") {
  if (!add) plot(x$outline, type="n", asp = asp, xlab = xlab, ylab = ylab, ...)
  polygon(x$base, col= col.base, border="Gray")
  lines(x$d.lambda, lwd=2, col= col.lambda, lty=2)
  lines(x$d.phi, lwd=2, col= col.phi, lty=2)
  lines(x$axis.major, lwd=2, col= col.major)
  lines(x$axis.minor, lwd=2, col = col.minor)
  lines(x$outline, asp=1, lwd=1, col = col.outline)
 invisible(NULL)
}
#' Ellipse
#'
#' @param center center
#' @param axes axes
#' @param scale scale
#' @param n n
#' @param from from
#' @param to to
#'
#' @return matrix
#' @export
ti_ellipse <- function(center, axes, scale=1, n=36, from=0, to=2*pi) {
  # Vector representation of an ellipse at "center" with axes in the *rows* of "axes".
  # Returns an "n" by 2 array of points, one per row.
  theta <- seq(from=from, to=to, length.out=n)
  t((scale * t(axes))  %*% rbind(cos(theta), sin(theta)) + center)
}

# plot.indicatrix <- function(x, y, ...) {
#   add <- list(...)$add
#   if (is.null(add) || !add) {
#     plot(x$outline, type="n", asp=1, xlab="x", ylab="y")
#   }
#   polygon(x$base, col=rgb(0, 0, 0, .025), border="Gray")
#   lines(x$d.lambda, lwd=2, col="Gray", lty=2)
#   lines(x$d.phi, lwd=2, col=rgb(.25, .7, .25), lty=2)
#   lines(x$axis.major, lwd=2, col=rgb(.25, .25, .7))
#   lines(x$axis.minor, lwd=2, col=rgb(.7, .25, .25))
#   lines(x$outline, asp=1, lwd=2)
# }

# smooth.cloud.R
# Create a smooth circle around the pth portion of a data cloud (centered at mean(x), mean(y))
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created March 24, 2006

smooth.cloud <- function(x, y, q, sf.arc = 1, sf.theta = 1, n = 50, na.rm=FALSE)
{
  require(Hmisc)
  
  # remove NAs
  if(na.rm)
  {
    x.new <- x[!is.na(x) & !is.na(y)]
    y.new <- y[!is.na(x) & !is.na(y)]
    
    x <- x.new
    y <- y.new
  }

  # center the cloud at (0,0) and convert to polar coordinates
  xbar <- mean(x)
  ybar <- mean(y)

  x <- x - xbar
  y <- y - ybar
  
  # smooth right half of circle
  z <- to.polar(x, y)
  smooth.r <- .smc(z$theta, z$r, n, sf.arc, sf.theta, q)
  smooth.r <- subset(to.cartesian(smooth.r$eval.points, smooth.r$estimate), x >= 0)

  # smooth left half of circle
  z <- to.polar(-x, y)
  smooth.l <- .smc(z$theta, z$r, n, sf.arc, sf.theta, q)
  smooth.l <- subset(to.cartesian(smooth.l$eval.points, smooth.l$estimate), x >= 0)
  smooth.l$x <- smooth.l$x * -1

  # put two halves together ... add extra point on end to complete circle for use in lines()
  n <- length(smooth.l$x)

  smooth <- data.frame(x = c(smooth.r$x, smooth.l$x[n:1], smooth.r$x[1]),
                       y = c(smooth.r$y, smooth.l$y[n:1], smooth.r$y[1]))

  # recenter cloud and return
  smooth$x <- smooth$x + xbar
  smooth$y <- smooth$y + ybar
  
  return(smooth)
}

.smc <- function(theta, r, n, sf.arc, sf.theta, q)
{
  # caclulate each half separately to avoid boundry problems
  eval.points <- seq(-pi / 2, pi / 2, length=round(n/2))

  width <- eval.points[2] - eval.points[1]

  estimate <- sapply(eval.points, function(pts)
                   {
                     # weight by change in theta
                     theta.w <- pnorm(abs(theta - pts), sd = sf.theta, lower.tail=FALSE) * 2
                     # weight by change in arc length to line defined by theta = pts, and r = r
                     s.w <- pnorm(r*abs(theta - pts), sd = sf.arc, lower.tail=FALSE) * 2
                     return(wtd.quantile(r, theta.w * s.w, probs=q))
                   }
                  )

  return(data.frame(eval.points = eval.points,
                    estimate = estimate))
}

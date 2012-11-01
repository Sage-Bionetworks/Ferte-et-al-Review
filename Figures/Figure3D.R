require(lattice)

x <- seq(-2.2*pi, 1.7*pi, len = 20)
y <- seq(-.3*pi,pi,len=20)
g <- expand.grid(x = x, y = y)
g$z <- sin(sqrt(g$x^2) +sqrt(g$y^2)) +  c(1:400)/(60*pi)

wireframe(z ~ x *y, g, drape = TRUE, shade=TRUE,
                aspect = c(1,1),pretty=TRUE)

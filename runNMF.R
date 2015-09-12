setwd("./stochasticNMF/")
require(bmp)
source("sNMF.R")
source("bNMF.R")

# show imput image
im <- read.bmp("lena512.bmp")
image(t(im)[, 512:1], col=gray.colors(64), axes=FALSE)

res.s <- snmf(im, 30, n.itr=100, conv.l=1e-2, conv.g=1e-2)
res.b <- bnmf(im, 30, n.itr=100, conv.g=1e-2)

fx <- (res.b$params$g/res.b$params$h) %*% (res.b$params$a/res.b$params$b)
image(t(fx)[, 512:1], col=gray.colors(64), axes=FALSE)




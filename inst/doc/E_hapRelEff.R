## ----setup,include=FALSE------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(Haplin, quietly = TRUE)

## ----eval=TRUE----------------------------------------------------------------
hapRelEff(
  cases.comp = c(c=1), 
  controls.comp = c(c=1),
  cases.ref = c(mfc=1),
  haplo.freq = c(0.1,0.9),
  RR = c(1,1)
)

## -----------------------------------------------------------------------------
hapRelEff(
  cases.comp = c(mfc=1), 
  controls.comp = c(mfc=1),
  cases.ref = c(mfc=1),
  haplo.freq = c(0.2,0.8),
  RR = c(1,1)
)

## -----------------------------------------------------------------------------
hapRelEff(
  cases.comp = c(mfc=1), 
  controls.comp = c(mfc=1),
  cases.ref = c(mfc=1),
  haplo.freq = c(0.2,0.8),
  RRcm = c(1,1),
  RRcf = c(1,1)
)

## ----eval=TRUE----------------------------------------------------------------
hapRelEff(
  cases.comp = list(c(mc=1)), 
  cases.ref = list(c(mfc=1)),
  haplo.freq = c(0.1,0.9),
  RR = c(1,1),
  RR.mat=c(1,1)
)

## ----eval=TRUE----------------------------------------------------------------
hapRelEff(
  nall = c(2,2),
  cases.comp = c(c=1), 
  controls.comp = c(c=1),
  cases.ref = c(mfc=1),
  haplo.freq = c(0.1,0.2,0.3,0.4),
  RR = c(1,1,1,1)
)

## -----------------------------------------------------------------------------
hapRelEff(
  cases.comp = c(mfc=1),
  controls.comp = c(mfc=1), 
  cases.ref = c(mfc=1),
  haplo.freq = c(0.8,0.2), 
  RRcm = c(1,2),
  RRcf = c(1,1),
  xchrom = T, 
  sim.comb.sex = "double",
  BR.girls = 1
)



if (FALSE) {
  length(fr@gam$coefficients) # 10L, default number of knots `k`
  (fr_old = Qindex(PFS ~ marker, data = Ki67q_0, k = 40))
  length(fr_old@gam$coefficients) # 40L, user choice of `k`
  # misung: in old simulation it was `k = ceiling({# of pct} * .4)`
}

if (FALSE) {
  length(nlfr@gam$coefficients) # 20L (=(5-1)*5), default
  # misung: to reproduce your old result
  (nlfr_misung = Qindex(PFS ~ marker, data = Ki67q_0, nonlinear = TRUE, k = 10))
  length(nlfr_misung@gam$coefficients) # 90L (=(10-1) * 10)
  (nlfr_old2 = Qindex(PFS ~ marker, data = Ki67q_0, nonlinear = TRUE, k = 7))
  length(nlfr_old2@gam$coefficients) # 42 (=(7-1)*7)
  (nlfr_old3 = Qindex(PFS ~ marker, data = Ki67q_0, nonlinear = TRUE, k = c(5,10)))
  length(nlfr_old3@gam$coefficients) # 5*(10-1)
}


if (FALSE) {
  integrandSurface(nlfr, newid = NULL)
  vis.gam(nlfr@gam)
  #autoplot.vis_gam_(vis_gam_(nlfr@gam))
} # these are mathematically identical

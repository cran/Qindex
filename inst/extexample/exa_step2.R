
### Step 2 (after Step 1)

## Step 2a: Linear Sign-Adjusted Quantile Indices
(fr = Qindex(PFS ~ marker, data = Ki67q_0))
stopifnot(all.equal.numeric(c(fr), predict(fr)))
\donttest{integrandSurface(fr)}
fr_1 = predict(fr, newdata = Ki67q_1)
\donttest{integrandSurface(fr, newdata = Ki67q_1)}

## Step 2b: Non-Linear Sign-Adjusted Quantile Indices
(nlfr = Qindex(PFS ~ marker, data = Ki67q_0, nonlinear = TRUE))
stopifnot(all.equal.numeric(c(nlfr), predict(nlfr)))
\donttest{integrandSurface(nlfr)}
nlfr_1 = predict(nlfr, newdata = Ki67q_1)
\donttest{integrandSurface(nlfr, newdata = Ki67q_1)}

## view linear and non-linear sign-adjusted quantile indices together
\donttest{integrandSurface(fr, nlfr)}

\donttest{
### Step 2c: Optimal Dichotomizing
set.seed(14837); (m1 = optimSplit_dichotom(
  PFS ~ marker, data = Ki67q_0, nsplit = 20L, top = 2L)) 
predict(m1)
predict(m1, boolean = FALSE)
predict(m1, newdata = Ki67q_1)
}
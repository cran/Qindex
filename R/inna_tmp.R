#require(gdata)
#require(rpart)
#require(survival)
#require(rpart.plot)
#require(spatstat)

#library(mgcv)
#library(refund)

#devtools::load_all("/Users/ixc107/Box/Qindex")

# These are temperary functions from Inna
# They may make their way to Qindex package eventually

######## ######## ######## Compute FRi's based on NNDQ

#' @importFrom survival Surv
FRiNNDQ <- function(
    data,
    Nknots,
    Nboot, 
    ...
) {	
  FRind <- FRindex(Surv(RECFREESURV_MO, RECURRENCE) ~ Marker, data = data, knot.value = Nknots, ...)
  data$FRi <- FRind@FRi
  data$FRi_std <- std_IQR(data$FRi)
  
  ############### FRi_std as continuous predictor
  cfit <- coxph(Surv(RECFREESURV_MO, RECURRENCE) ~ FRi_std, data = data) # I think this is wrong
  # cfit <- coxph(Surv(RECFREESURV_MO, RECURRENCE) ~ FRi, data = data) # I think this is correct
  
  ############### FRi as dichotomized predictor
  Dfun <- rpartD(y = with(data, Surv(RECFREESURV_MO, RECURRENCE)), x = data$FRi)
  data$d_FRi <- Dfun(data$FRi)
  trfit <- coxph(Surv(RECFREESURV_MO, RECURRENCE) ~ d_FRi, data = data)
  
  ############### dichotomized FRi witt bias correction
  m1 = BBC_dichotom(Surv(RECFREESURV_MO, RECURRENCE) ~ 1, data = data, dichotom = ~ FRi, R = Nboot)
  
  output = rbind(summary(cfit)$coefficients, summary(trfit)$coefficients, summary(m1)$coefficients)
  output
}


 











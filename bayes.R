library(mlrMBO)
library(ggplot2)
library(DiceKriging)
library(rjson)
try({#when sourcing the file
  dr = getSrcDirectory()[1]
})
try({#for running code directly in rstudio
  dr = dirname(rstudioapi::getActiveDocumentContext()$path)
})

setwd(dr) #set for usage as a script/running lines
stronk=c(2,3,9,14,17,21)

source('functions.R')

spctr = fromJSON(file = "panel1_lymph_subtracted_fixed.json")
spectra = matrix(ncol = 64, nrow = length(spctr))
spec_names = c()
for (i in 1:length(spctr)) {
  spectra[i, ] = spctr[[i]]$spectrum$mS
  spec_names = c(spec_names, spctr[[i]]$antigen)
}
colnames(spectra) = spctr[[1]]$spectrum$channels
rm(spctr)
rownames(spectra) = spec_names
nspectra = nrow(spectra)
ndetects = ncol(spectra)

pheno = read.csv('phenotypes_noAFonly.csv',
                 sep = ',',
                 fileEncoding = "UTF-8-BOM")
todel = c()
pheno = recur_unpack(c('base'), pheno, todel)

set.seed(1337)
resa = generate_artificial_cytometry_data(spectra, pheno, 100000, stronk,uncertainty=1)
set.seed(42)
resb = generate_artificial_cytometry_data(spectra, pheno, 100000, stronk,uncertainty=1.5)
set.seed(69)
resc = generate_artificial_cytometry_data(spectra, pheno, 100000, c(),uncertainty=1)
set.seed(420)
resd = generate_artificial_cytometry_data(spectra, pheno, 100000, c(),uncertainty=1.5)

fn = makeSingleObjectiveFunction(
  name = "ngd_eval",
  fn = function(x) {
    mse = 0
    for (r in list(resa,resb,resc,resd)) { 
      ngd = trans(
        nougadmt(
          mixed = r[[1]],
          spectra = spectra,
          snw = x['snw'],
          spw = x['spw'],
          nw = x['nw'],
          start = x['start'],
          iters = 500,
          alpha = x['alpha'] * 0.001,
          accel = x['accel'] * 0.001,
          nthreads = 24
        )$unmixed
      )
      mse_curr = sum((ngd - trans(r[[2]])) ^ 2) / length(ngd)
      mse = mse + mse_curr
    }
     if (is.nan(mse)) {
       mse = 40
     }
     if (!is.finite(mse)) {
       mse = 40
     }
    mse
    
  },
  par.set = makeParamSet(
    makeIntegerParam("snw", lower = 1, upper = 300),
    makeIntegerParam("spw", lower = 1, upper = 300),
    makeIntegerParam("nw", lower = 1, upper = 2000),
    makeIntegerParam("start", lower = 1, upper = 500),
    makeIntegerParam("alpha", lower = 1, upper = 100),
    makeIntegerParam("accel", lower = 1, upper = 1000)
  )
)

print(fn)
surr.km = makeLearner(
  "regr.km",
  predict.type = "se",
  covtype = "gauss",
  control = list(trace = FALSE)
)

ctrl =  makeMBOControl()
ctrl = setMBOControlTermination(ctrl, iters = 2000, target.fun.value = 1e-20)
ctrl = setMBOControlInfill(ctrl, makeMBOInfillCritEI())
run2 = mbo(fn,
          learner = surr.km,
          control = ctrl,
          show.info = TRUE)



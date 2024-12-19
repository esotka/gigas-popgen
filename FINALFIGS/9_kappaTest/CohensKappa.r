### Cohens Kappa
library(psych)
dat = read.delim("forCohensKappa.txt")
hyp =  dat$Historical_numeral
inf =  dat$Genetic_numeral
print(cohen.kappa(x=cbind(hyp,inf)))

#Call: cohen.kappa1(x = x, w = w, n.obs = n.obs, alpha = alpha, levels = levels)

#Cohen Kappa and Weighted Kappa correlation coefficients and confidence boundaries 
#                  lower estimate upper
#unweighted kappa -0.098     0.26  0.62
#weighted kappa    0.139     0.51  0.88

#Number of subjects = 14 

# remove uncertain species (by genetic inference)

dat2 = dat[!dat$Genetic.Inference=="Unknown",]
hyp2 =  dat2$Historical_numeral
inf2 =  dat2$Genetic_numeral
print(cohen.kappa(x=cbind(hyp2,inf2)))

#Cohen Kappa and Weighted Kappa correlation coefficients and confidence boundaries 
#                 lower estimate upper
#unweighted kappa 0.021     0.34  0.65
#weighted kappa   0.147     0.56  0.97

 #Number of subjects = 11 


#Cohen suggested the Kappa result be interpreted as follows: values ≤ 0 as indicating no agreement and 0.01–0.20 as none to slight, 0.21–0.40 as fair, 0.41– 0.60 as moderate, 0.61–0.80 as substantial, and 0.81–1.00 as almost perfect agreement.

# Unweighted kappa, therefore, is inappropriate for ordinal scales Sim J, Wright CC (2005) The kappa statistic in reliability studies: use, interpretation, and sample size requirements. Phys Ther 85:257–268


# load in packages
library(tidyr)
library(dplyr)
library(ggplot2)
library(nlme)

# make our own function
baranyi_dose_response <- function(LOG10N0, mumax, conc, lag){
  return(LOG10N0 + mumax * conc/log(10) + log10(exp(-mumax * conc) * 
                                               (1 - exp(-mumax * lag)) + exp(-mumax * lag)))
}

# weibull model 1 from drc
weibull1 <- function(dose, b, c, d, e){
  return(c+(d - c)*exp(-exp(b*(log(dose) - log(e)))))
}
weibull1 <- function(dose, b, c, d){
  return(c*exp(b*(log(dose))) + d)
}

dose = 1:100

plot(weibull1(dose, 2, 4))

# Analysis of gentamycin ####
# load in data
# gentamycin
d_gent <- readxl::read_excel('data/Competition Fitness data sheet.xlsx', range = 'A1:I13') %>%
  gather(., 'conc', 'fitness', c(2:7)) %>%
  mutate(., conc = as.numeric(conc),
         # change conc 0 to 0.001
         conc = ifelse(conc == 0, 0.001, conc)) %>%
  data.frame()

# quick plot                             
ggplot(d_gent) +
  geom_point(aes(log10(conc), fitness, group = rep, col = community)) +
  facet_wrap(~ community)

# non-linear regression - using Baranyi without lag ####

# prep data - log10 the concentrations
# normalised the concentration so that 0 is the smallest concentration (essentially added three to everything)
d_gent_fit <- mutate(d_gent, log_conc = log10(conc),
                     log_conc_norm = log_conc + abs(min(log_conc)))

# fit global model  
fit_gent <- gnls(fitness ~ weibull1(conc, b, c, d),
                  data = d_gent_fit,
                  params = b + c + d ~ community,
                  start = c(0.7, 0, 1, 0, 0, 0),
                  control = nls.control(maxiter = 1000))

# smallest t value / largest p value is lag, take that out
fit_gent2 <- gnls(fitness ~ weibull1(LOG10N0, mumax, conc = log_conc_norm, lag),
                 data = d_gent_fit,
                 params = c(mumax + LOG10N0 ~ community,
                            lag ~ 1),
                 start = c(0.7, 0, 1, 0, 1),
                 control = nls.control(maxiter = 1000))

# take out mumax next
fit_gent3 <- gnls(fitness ~ baranyi_dose_response(LOG10N0, mumax, conc = log_conc_norm, lag),
                  data = d_gent_fit,
                  params = c(LOG10N0 ~ community,
                             lag + mumax ~ 1),
                  start = c(0.7, 0, 1, 1),
                  control = nls.control(maxiter = 1000))

# take out LOG10NO
fit_gent4 <- gnls(fitness ~ baranyi_dose_response(LOG10N0, mumax, conc = log_conc_norm, lag),
                  data = d_gent_fit,
                  params = c(LOG10N0 + lag + mumax ~ 1),
                  start = c(0.7, 1, 1),
                  control = nls.control(maxiter = 1000))

# do model comparison on these
anova(fit_gent, fit_gent2, fit_gent3, fit_gent4)
AIC(fit_gent, fit_gent2, fit_gent3, fit_gent4)

# so fit_gent3 is the best fit !!!

# need to think what the best model is and what the parameters in this model mean...
# or whether there is a different dose-response model you would like to fit

# create predictions ####
preds_gent <- data.frame(expand.grid(conc = seq(min(d_gent_fit$conc), max(d_gent_fit$conc), length.out = 1000), community = c('Y', 'N'), stringsAsFactors = FALSE)) %>%
  mutate(., pred = predict(fit_gent, .))

# plot ####
ggplot(d_gent) +
  geom_point(aes(log10(conc), fitness, group = rep, col = community)) + 
  geom_line(aes(log10(conc), pred, col = community), preds_gent) +
  theme_bw()

# do the same with kanamycin
# kanamycin
d_kan <- readxl::read_excel('data/Competition Fitness data sheet.xlsx', range = 'A17:I29') %>%
  gather(., 'conc', 'fitness', c(2:7)) %>%
  mutate(., conc = as.numeric(conc),
         conc = ifelse(conc == 0, 0.002, conc)) %>%
  data.frame()

ggplot(d_kan) +
  geom_point(aes(log10(conc), fitness, group = rep, col = community))

# non-linear regression - using Baranyi without lag ####

# prep data - log10 the concentrations
# normalised the concentration so that 0 is the smallest concentration (essentially added three to everything)
d_kan_fit <- mutate(d_kan, log_conc = log10(conc),
                     log_conc_norm = log_conc + abs(min(log_conc)))

# fit global model  
fit_kan <- gnls(fitness ~ baranyi_dose_response(LOG10N0, mumax, conc = log_conc_norm, lag),
                 data = d_kan_fit,
                 params = mumax + lag + LOG10N0 ~ community,
                 start = c(0.7, 0, 1, 0, 1, 0),
                 control = nls.control(maxiter = 1000))

# smallest t value / largest p value is LOG10N0, take that out
fit_kan2 <- gnls(fitness ~ baranyi_dose_response(LOG10N0, mumax, conc = log_conc_norm, lag),
                  data = d_kan_fit,
                  params = c(mumax + lag ~ community,
                             LOG10N0 ~ 1),
                  start = c(0.7, 0, 1, 0, 1),
                  control = nls.control(maxiter = 1000))

# take out mumax next
fit_kan3 <- gnls(fitness ~ baranyi_dose_response(LOG10N0, mumax, conc = log_conc_norm, lag),
                  data = d_kan_fit,
                  params = c(lag ~ community,
                             LOG10N0 + mumax ~ 1),
                  start = c(0.7, 0, 1, 1),
                  control = nls.control(maxiter = 1000))

# take out LOG10NO
fit_kan4 <- gnls(fitness ~ baranyi_dose_response(LOG10N0, mumax, conc = log_conc_norm, lag),
                  data = d_kan_fit,
                  params = c(LOG10N0 + lag + mumax ~ 1),
                  start = c(0.7, 1, 1),
                  control = nls.control(maxiter = 1000))

# do model comparison on these
anova(fit_kan, fit_kan2, fit_kan3)
AIC(fit_kan, fit_kan2, fit_kan3)

# so fit_kan2 is the best fit !!!

# need to think what the best model is and what the parameters in this model mean...
# or whether there is a different dose-response model you would like to fit

# create predictions ####
preds_kan <- data.frame(expand.grid(log_conc_norm = seq(min(d_kan_fit$log_conc_norm), max(d_kan_fit$log_conc_norm), length.out = 50), community = c('Y', 'N'), stringsAsFactors = FALSE)) %>%
  mutate(., pred = predict(fit_kan2, .))

# plot ####
ggplot(d_kan) +
  geom_point(aes(log10(conc), fitness, group = rep, col = community)) + 
  geom_line(aes(log_conc_norm + min(d_kan_fit$log_conc), pred, col = community), preds_kan) +
  theme_bw()

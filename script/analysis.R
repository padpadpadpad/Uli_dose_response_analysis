# load in packages
library(tidyr)
library(dplyr)
library(ggplot2)
library(nlme)
library(rstan)

# weibull model 1 from drc
weibull1 <- function(dose, b, c, d, e){
  return(c+(d - c)*exp(-exp(b*(log(dose) - log(e)))))
}

# edited version with no asymptote and an intercept
weibull1 <- function(dose, b, c, d){
  return(c*exp(b*(log(dose))) + d)
}

# function for working out x at a specific y
x_at_spec_y <- function(xmin, xmax, model, y, colname_x, fac, colname_fac){
  if(missing(colname_fac)){
    temp = data.frame(x = seq(xmin, xmax, length.out = 100000))
    colnames(temp) <- colname_x
    temp$z <- predict(model, temp)
    return(temp[which.min(abs(temp$z-y)),1])
  }
  if(!missing(colname_fac)){
    temp = data.frame(expand.grid(x = seq(xmin, xmax, length.out = 100000), fac = fac, stringsAsFactors = FALSE))
    colnames(temp) <- c(colname_x, colname_fac)
    temp$z <- predict(model, temp)
    temp$w <- abs(temp$z - y)
    temp <- group_by_(temp, colname_fac) %>%
      arrange(., w) %>%
      top_n(., -1, w) %>%
      data.frame()
    return(select_(temp, colname_x, colname_fac))
  }
}

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

# fit global model  
fit_gent <- gnls(fitness ~ weibull1(conc, b, c, d),
                  data = d_gent,
                  params = b + c + d ~ community,
                  start = c(0.7, 0, 1, 0, 0, 0),
                  control = nls.control(maxiter = 1000))

# fit d_gent where the intercept is set to 0
fit_gent2 <- gnls(fitness ~ weibull1(conc, b, c, d = 0),
                 data = d_gent_fit,
                 params = b + c~ community,
                 start = c(0.7, 0, 1, 0),
                 control = nls.control(maxiter = 1000))

# compare models
AIC(fit_gent,
    fit_gent2)
anova(fit_gent,
      fit_gent2)

# smallest t value / largest p value is lag, take that out
fit_gent2 <- gnls(fitness ~ weibull1(conc, b, c, d),
                  data = d_gent,
                  params = c(b + d ~ community,
                                c ~ 1),
                  start = c(0.3, 0.17, 0.11, 0.9, -0.15),
                  control = nls.control(maxiter = 1000))

# take out mumax next
fit_gent3 <- gnls(fitness ~ weibull1(conc, b, c, d),
                  data = d_gent,
                  params = c(d ~ community,
                             b + c ~ 1),
                  start = c(0.3, 0.11, 0.9, -0.15),
                  control = nls.control(maxiter = 1000))

# do model comparison on these
anova(fit_gent, fit_gent2, fit_gent3)
MuMIn::AICc(fit_gent, fit_gent2, fit_gent3)

# so fit_gent3 is the best fit !!!

# create predictions ####
preds_gent <- data.frame(expand.grid(conc = seq(min(d_gent$conc), max(d_gent$conc), length.out = 1000), community = c('Y', 'N'), stringsAsFactors = FALSE)) %>%
  mutate(., pred = predict(fit_gent3, .))

# Use a bayesian method to get prediction intervals around the models
# Will fit in rstan

# add dummy variable for community, 0 or 1
preds_gent <- mutate(preds_gent, comm_dummy = ifelse(community == 'Y', 1, 0))
d_gent <- mutate(d_gent, comm_dummy = ifelse(community == 'Y', 1, 0))

# create a data list for stan ####
stan_data_list <- list(N = nrow(d_gent),
                       conc = d_gent$conc,
                       comm = d_gent$comm_dummy,
                       fitness = d_gent$fitness,
                       Nnew = nrow(preds_gent),
                       conc_new = preds_gent$conc,
                       comm_new = preds_gent$comm_dummy
                       )

# initial values function ####
start <- function() list(b = 1, c = 0, d = 0, d_comm = 0, sigma = 1)

# 1. run stan model without extensive predictions ####
model_stan <- stan(file = 'script/bayesian_model_no_preds.stan',
                   data = stan_data_list, 
                   init = start, 
                   iter = 1e4, 
                   chains = 3)

# view parameter estimates
model_stan
# these are exactly the same as fit_gent3!!! This is good.

# plot traceplot - this also looks great
traceplot(model_stan)

# 2. run model while generating predictions
# fitting loads of predictions so that we get the value of x closest to a y of 1
model_stan_pred <- stan(file = 'script/bayesian_model.stan',
                   data = stan_data_list, 
                   init = start, 
                   iter = 1e4, 
                   chains = 3)

# extract samples, keep means and then compute mean, 95% credible intervals for new predictions
samples <- rstan::extract(model_stan_pred)
mu_samps <- data.frame(samples[['mu_new']]) %>%
  gather(., 'Nnew', 'samps') %>%
  mutate(., Nnew = as.numeric(gsub('X', '', Nnew))) %>%
  group_by(., Nnew) %>%
  summarise(., mu = quantile(samps, probs = 0.5),
            lwr_CI = quantile(samps, probs = 0.025),
            upr_CI = quantile(samps, probs = 0.975)) %>%
  data.frame()
sim_samps <- data.frame(samples[['sim_new']]) %>%
  gather(., 'Nnew', 'samps') %>%
  mutate(., Nnew = as.numeric(gsub('X', '', Nnew))) %>%
  group_by(., Nnew) %>%
  summarise(., lwr_PI = quantile(samps, probs = 0.025),
            upr_PI = quantile(samps, probs = 0.975)) %>%
  data.frame()

# bind with predictions dataframe
p_gent <- cbind(select(preds_gent, conc, community), select(mu_samps, -Nnew), select(sim_samps, -Nnew))

# find the x where fitness = 1
mu_at_fitness_1 <- mutate(p_gent, w = abs(mu - 1)) %>%
  group_by(., community) %>%
  top_n(., -1, w) %>%
  select(., conc, community) %>%
  data.frame() %>%
  rename(., mu = conc)
lwrCI_at_fitness_1 <- mutate(p_gent, w = abs(lwr_CI - 1)) %>%
  group_by(., community) %>%
  top_n(., -1, w) %>%
  select(., conc, community) %>%
  data.frame() %>%
  rename(., lwr_CI = conc)
uprCI_at_fitness_1 <- mutate(p_gent, w = abs(upr_CI - 1)) %>%
  group_by(., community) %>%
  top_n(., -1, w) %>%
  select(., conc, community) %>%
  data.frame() %>%
  rename(., upr_CI = conc)
x_at_fitness_1 <- merge(mu_at_fitness_1,
                        lwrCI_at_fitness_1,
                        by = 'community') %>%
  merge(., uprCI_at_fitness_1, by = 'community')

# plot model and predictions ####
# on a log10 scale
ggplot(d_gent) +
  geom_point(aes(log10(conc), fitness, group = rep, col = community)) + 
  geom_line(aes(log10(conc), mu, col = community), p_gent) +
  geom_ribbon(aes(log10(conc), ymin = lwr_CI, ymax = upr_CI, fill = community), alpha = 0.2, p_gent) +
  theme_bw()

# plot the intervals of where the lwr_ci and upr_ci and mu cross 1
ggplot(x_at_fitness_1) +
  geom_point(aes(community, mu, col = community), size = 3) +
  geom_errorbar(aes(community, ymin = lwr_CI, ymax = upr_CI, col = community), width = 0.1) +
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

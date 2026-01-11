rm(list = ls())


library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(odin2)
library(dust2)
library(monty)
library(openxlsx)



# library(ggplot2)


# ---- Odin model --------------------------------------------------------
gen_sir <- odin2::odin('sir.R')


# ---- Load & shape data ----------------------------------------------------
# Total cases (weekly) and typing (cv16, ev71, other, neg)
data_observed <- readRDS("combine_data_2009_2022_final.rds") %>%
  filter(Region == "Peninsular Malaysia") %>%
  filter(OnsetDate >= as.Date("2017-01-01") & OnsetDate <= as.Date("2022-12-31")) %>%
  mutate(week = lubridate::week(OnsetDate), year = lubridate::year(OnsetDate)) %>%
  mutate(week = ifelse(week == 53, 52, week)) %>%
  group_by(year, week) %>%
  summarise(cases_total = n(), .groups = 'drop') %>%
  group_by(year) %>%
  mutate(week = row_number()) %>%
  ungroup() %>%
  mutate(time = 1:n())

data_typing <- read.csv("weekly_serotype_number.csv") %>%
  select(-c(1)) 

# merge total + typing (fill missing typing with 0)
data <- data_observed %>%
  left_join(data_typing) %>%
  mutate(across(c(cv16, ev71, other_ev,total), ~ replace_na(., 0))) %>%
  arrange(time) %>%
  rowwise() %>% 
  mutate(cases = I(list(c_across(c(cv16, ev71,other_ev))))) %>%
  rename("hfmd_cases" = "cases_total") %>%
  select(time, hfmd_cases, cases)
  

# ---- Demography -----------------------------------------------------------
pop_data <-  read.csv('pop_clean.csv') %>% mutate(week = c(0, 52, 104, 156, 209, 261, 312))

# ---- Holiday & NPI series -------------------------------------------------
npi <- read.csv('npi_clean.csv') %>% mutate(week = row_number() - 1L)
holiday <- read.csv('holiday_clean.csv') %>% mutate(week = row_number() - 1L)


# ---- Fixed parameters -----------------------------------------------------
fixed_pars <- list(
  gamma = 1,
  # w = 0,
  # p_4 = 0, 
  N_0 = pop_data$total[1],
  birth_time = pop_data$week, birth_value = pop_data$weekly_birth_rate,
  death_time = pop_data$week, death_value = pop_data$weekly_death_rate,
  chi_time = npi$week, chi_value = npi$chi,
  holiday_time = holiday$week, holiday_value = holiday$value,
  covid_start = 168,   # school close 1:  as.Date("2020-03-18") - as.Date("2017-01-01")
  covid_end = 276  # school close 2:  as.Date("2022-04-10") - as.Date("2017-01-01")
)


# ---- Monty packer ---------------------------------------------------------
estimated_pars <- c(
  'beta_0_1', 's_0_1', 'i_0_1', 'p_1', 
  'c1', 's_0_2', 'i_0_2',  
  'c2', 's_0_3', 'i_0_3', 
  'r0','r1','h','w'
)

packer <- monty_packer(estimated_pars, fixed = fixed_pars)

# ---- Priors (match your mcstate ranges) -----------------------------------
prior <- monty_dsl({
  beta_0_1 ~ Uniform(0, 25)  
  s_0_1 ~ Uniform(0, 0.5)
  i_0_1 ~ Uniform(0, 0.02)
  p_1 ~ Uniform(0, 1)
  
  c1 ~ Uniform(0.1, 10)  
  # beta_0_2 ~ Uniform(0, 25)  
  s_0_2 ~ Uniform(0, 0.5)
  i_0_2 ~ Uniform(0, 0.02)
  # p_2 ~ Uniform(0, 1)
  
  c2 ~ Uniform(0.1, 10)   
  # beta_0_3 ~ Uniform(0, 25)   
  s_0_3 ~ Uniform(0, 0.5)
  i_0_3 ~ Uniform(0, 0.02)
  # p_3 ~ Uniform(0, 1)

  r0  ~ Uniform(0, 20)
  r1  ~ Uniform(0, 0.5)
  h   ~ Uniform(0.5, 1)
  # se  ~ Uniform(0, 1)
  w  ~ Uniform(0, 1000)
})

# ---- MCMC ---------------
# deterministic mcmc
initial_values <- c(10, 0.1, 1e-3, 0.01,  # 'beta_0_1', 's_0_1', 'i_0_1', 'p_1', 
                    1, 0.1, 1e-3,  # 'c1', 's_0_2', 'i_0_2', 'p_2',
                    1, 0.1, 1e-3,   # 'c2', 's_0_3', 'i_0_3', 'p_3',
                    4, 0.05, 0.8,10) # 'r0','r1','h', 'w'

vcv <- diag(c(0.1, 0.01, 1e-6, 1e-4,
              0.01, 0.01, 1e-6,
              0.01, 0.01, 1e-6, 
              0.1, 1e-4, 0.01, 1), 14)

filter <- dust_filter_create(gen_sir(), data = data, time_start = 0, n_particles = 2000, n_threads = 4)
likelihood <- dust_likelihood_monty(filter, packer, save_trajectories = TRUE, save_state = TRUE)
posterior <- likelihood + prior
runner <- monty_runner_callr(4)

# 1
sampler <- monty_sampler_random_walk(vcv, rerun_every = 50, rerun_random = TRUE)
samples <- monty_sample(posterior, sampler, initial = initial_values, n_steps = 5000, n_chains = 4, runner = runner)
write_rds(samples, "samples-1.rds")


# ########################
# # 2
# ########################
samples <- readRDS("samples-1.rds")
samples_df <- posterior::as_draws_df(samples)

vcv <- cov(samples_df[, 1:(ncol(samples_df)-3)])
mcmc_list <- coda::as.mcmc.list(samples)
if((1 - coda::rejectionRate(mcmc_list))[1] < 0.01) { vcv <- vcv * 0.01 }

initial_values <- posterior::summarise_draws(samples_df) %>% dplyr::pull(mean)
sampler <- monty_sampler_random_walk(vcv, rerun_every = 50, rerun_random = TRUE)
samples <- monty_sample(posterior, sampler, initial = initial_values,
                        n_steps = 5000, n_chains = 4, runner = runner)
write_rds(samples, "samples-2.rds")


# ########################
# # 3
# ########################
samples <- readRDS("samples-2.rds")
samples_df <- posterior::as_draws_df(samples)

vcv <- cov(samples_df[, 1:(ncol(samples_df)-3)])
mcmc_list <- coda::as.mcmc.list(samples)
if((1 - coda::rejectionRate(mcmc_list))[1] < 0.01) { vcv <- vcv * 0.01 }

initial_values <- posterior::summarise_draws(samples_df) %>% dplyr::pull(mean)
sampler <- monty_sampler_random_walk(vcv, rerun_every = 50, rerun_random = TRUE)
samples <- monty_sample(posterior, sampler, initial = initial_values,
                        n_steps = 5000, n_chains = 4, runner = runner)
write_rds(samples, "samples-3.rds")
# 
# 
# ########################
# # 4
# ########################
samples <- readRDS("samples-3.rds")
samples_df <- posterior::as_draws_df(samples)

vcv <- cov(samples_df[, 1:(ncol(samples_df)-3)])
mcmc_list <- coda::as.mcmc.list(samples)
if((1 - coda::rejectionRate(mcmc_list))[1] < 0.01) { vcv <- vcv * 0.01 }

initial_values <- posterior::summarise_draws(samples_df) %>% dplyr::pull(mean)
sampler <- monty_sampler_random_walk(vcv, rerun_every = 50, rerun_random = TRUE)
samples <- monty_sample(posterior, sampler, initial = initial_values,
                        n_steps = 5000, n_chains = 4, runner = runner)
write_rds(samples, "samples-4.rds")


########################
# 5
########################
samples <- readRDS("samples-4.rds")
samples_df <- posterior::as_draws_df(samples)

vcv <- cov(samples_df[, 1:(ncol(samples_df)-3)])
mcmc_list <- coda::as.mcmc.list(samples)
if((1 - coda::rejectionRate(mcmc_list))[1] < 0.01) { vcv <- vcv * 0.01 }

initial_values <- posterior::summarise_draws(samples_df) %>% dplyr::pull(mean)
sampler <- monty_sampler_random_walk(vcv, rerun_every = 50, rerun_random = TRUE)
samples <- monty_sample(posterior, sampler, initial = initial_values,
                        n_steps = 5000, n_chains = 4, runner = runner)
write_rds(samples, "samples-5.rds")


########################
# 6
########################
samples <- readRDS("samples-5.rds")
samples_df <- posterior::as_draws_df(samples)

vcv <- cov(samples_df[, 1:(ncol(samples_df)-3)])
mcmc_list <- coda::as.mcmc.list(samples)
if((1 - coda::rejectionRate(mcmc_list))[1] < 0.01) { vcv <- vcv * 0.01 }

initial_values <- posterior::summarise_draws(samples_df) %>% dplyr::pull(mean)
sampler <- monty_sampler_random_walk(vcv, rerun_every = 50, rerun_random = TRUE)
samples <- monty_sample(posterior, sampler, initial = initial_values,
                        n_steps = 5000, n_chains = 4, runner = runner)
write_rds(samples, "samples-6.rds")


########################
# 7
########################
samples <- readRDS("samples-6.rds")
samples_df <- posterior::as_draws_df(samples)

vcv <- cov(samples_df[, 1:(ncol(samples_df)-3)])
mcmc_list <- coda::as.mcmc.list(samples)
if((1 - coda::rejectionRate(mcmc_list))[1] < 0.01) { vcv <- vcv * 0.01 }

initial_values <- posterior::summarise_draws(samples_df) %>% dplyr::pull(mean)
sampler <- monty_sampler_random_walk(vcv, rerun_every = 50, rerun_random = TRUE)
samples <- monty_sample(posterior, sampler, initial = initial_values,
                        n_steps = 5000, n_chains = 4, runner = runner)
write_rds(samples, "samples-7.rds")


########################
# 8
########################
samples <- readRDS("samples-7.rds")
samples_df <- posterior::as_draws_df(samples)

vcv <- cov(samples_df[, 1:(ncol(samples_df)-3)])
mcmc_list <- coda::as.mcmc.list(samples)
if((1 - coda::rejectionRate(mcmc_list))[1] < 0.01) { vcv <- vcv * 0.01 }

initial_values <- posterior::summarise_draws(samples_df) %>% dplyr::pull(mean)
sampler <- monty_sampler_random_walk(vcv, rerun_every = 50, rerun_random = TRUE)
samples <- monty_sample(posterior, sampler, initial = initial_values,
                        n_steps = 5000, n_chains = 4, runner = runner)
write_rds(samples, "samples-8.rds")


########################
# 9
########################
samples <- readRDS("samples-8.rds")
samples_df <- posterior::as_draws_df(samples)

vcv <- cov(samples_df[, 1:(ncol(samples_df)-3)])
mcmc_list <- coda::as.mcmc.list(samples)
if((1 - coda::rejectionRate(mcmc_list))[1] < 0.01) { vcv <- vcv * 0.01 }

initial_values <- posterior::summarise_draws(samples_df) %>% dplyr::pull(mean)
sampler <- monty_sampler_random_walk(vcv, rerun_every = 50, rerun_random = TRUE)
samples <- monty_sample(posterior, sampler, initial = initial_values,
                        n_steps = 5000, n_chains = 4, runner = runner)
write_rds(samples, "samples-9.rds")


########################
# 10
########################
samples <- readRDS("samples-9.rds")
samples_df <- posterior::as_draws_df(samples)

vcv <- cov(samples_df[, 1:(ncol(samples_df)-3)])
mcmc_list <- coda::as.mcmc.list(samples)
if((1 - coda::rejectionRate(mcmc_list))[1] < 0.01) { vcv <- vcv * 0.01 }

initial_values <- posterior::summarise_draws(samples_df) %>% dplyr::pull(mean)
sampler <- monty_sampler_random_walk(vcv, rerun_every = 50, rerun_random = TRUE)
samples <- monty_sample(posterior, sampler, initial = initial_values,
                        n_steps = 5000, n_chains = 4, runner = runner)
write_rds(samples, "samples-10.rds")


########################
# 11
########################
samples <- readRDS("samples-10.rds")
samples_df <- posterior::as_draws_df(samples)

vcv <- cov(samples_df[, 1:(ncol(samples_df)-3)])
mcmc_list <- coda::as.mcmc.list(samples)
if((1 - coda::rejectionRate(mcmc_list))[1] < 0.01) { vcv <- vcv * 0.01 }

initial_values <- posterior::summarise_draws(samples_df) %>% dplyr::pull(mean)
sampler <- monty_sampler_random_walk(vcv, rerun_every = 50, rerun_random = TRUE)
samples <- monty_sample(posterior, sampler, initial = initial_values,
                        n_steps = 5000, n_chains = 4, runner = runner)
write_rds(samples, "samples-11.rds")


########################
# 12
########################
samples <- readRDS("samples-11.rds")
samples_df <- posterior::as_draws_df(samples)

vcv <- cov(samples_df[, 1:(ncol(samples_df)-3)])
mcmc_list <- coda::as.mcmc.list(samples)
if((1 - coda::rejectionRate(mcmc_list))[1] < 0.01) { vcv <- vcv * 0.01 }

initial_values <- posterior::summarise_draws(samples_df) %>% dplyr::pull(mean)
sampler <- monty_sampler_random_walk(vcv, rerun_every = 50, rerun_random = TRUE)
samples <- monty_sample(posterior, sampler, initial = initial_values,
                        n_steps = 5000, n_chains = 4, runner = runner)
write_rds(samples, "samples-12.rds")


########################
# 13
########################
samples <- readRDS("samples-12.rds")
samples_df <- posterior::as_draws_df(samples)

vcv <- cov(samples_df[, 1:(ncol(samples_df)-3)])
mcmc_list <- coda::as.mcmc.list(samples)
if((1 - coda::rejectionRate(mcmc_list))[1] < 0.01) { vcv <- vcv * 0.01 }

initial_values <- posterior::summarise_draws(samples_df) %>% dplyr::pull(mean)
sampler <- monty_sampler_random_walk(vcv, rerun_every = 50, rerun_random = TRUE)
samples <- monty_sample(posterior, sampler, initial = initial_values,
                        n_steps = 5000, n_chains = 4, runner = runner)
write_rds(samples, "samples-13.rds")


########################
# 14
########################
samples <- readRDS("samples-13.rds")
samples_df <- posterior::as_draws_df(samples)

vcv <- cov(samples_df[, 1:(ncol(samples_df)-3)])
mcmc_list <- coda::as.mcmc.list(samples)
if((1 - coda::rejectionRate(mcmc_list))[1] < 0.01) { vcv <- vcv * 0.01 }

initial_values <- posterior::summarise_draws(samples_df) %>% dplyr::pull(mean)
sampler <- monty_sampler_random_walk(vcv, rerun_every = 50, rerun_random = TRUE)
samples <- monty_sample(posterior, sampler, initial = initial_values,
                        n_steps = 5000, n_chains = 4, runner = runner)
write_rds(samples, "samples-14.rds")


########################
# 15
########################
samples <- readRDS("samples-14.rds")
samples_df <- posterior::as_draws_df(samples)

vcv <- cov(samples_df[, 1:(ncol(samples_df)-3)])
mcmc_list <- coda::as.mcmc.list(samples)
if((1 - coda::rejectionRate(mcmc_list))[1] < 0.01) { vcv <- vcv * 0.01 }

initial_values <- posterior::summarise_draws(samples_df) %>% dplyr::pull(mean)
sampler <- monty_sampler_random_walk(vcv, rerun_every = 50, rerun_random = TRUE)
samples <- monty_sample(posterior, sampler, initial = initial_values,
                        n_steps = 5000, n_chains = 4, runner = runner)
write_rds(samples, "samples-15.rds")


########################
# 16
########################
samples <- readRDS("samples-15.rds")
samples_df <- posterior::as_draws_df(samples)

vcv <- cov(samples_df[, 1:(ncol(samples_df)-3)])
mcmc_list <- coda::as.mcmc.list(samples)
if((1 - coda::rejectionRate(mcmc_list))[1] < 0.01) { vcv <- vcv * 0.01 }

initial_values <- posterior::summarise_draws(samples_df) %>% dplyr::pull(mean)
sampler <- monty_sampler_random_walk(vcv, rerun_every = 50, rerun_random = TRUE)
samples <- monty_sample(posterior, sampler, initial = initial_values,
                        n_steps = 5000, n_chains = 4, runner = runner)
write_rds(samples, "samples-16.rds")


########################
# 17
########################
samples <- readRDS("samples-16.rds")
samples_df <- posterior::as_draws_df(samples)

vcv <- cov(samples_df[, 1:(ncol(samples_df)-3)])
mcmc_list <- coda::as.mcmc.list(samples)
if((1 - coda::rejectionRate(mcmc_list))[1] < 0.01) { vcv <- vcv * 0.01 }

initial_values <- posterior::summarise_draws(samples_df) %>% dplyr::pull(mean)
sampler <- monty_sampler_random_walk(vcv, rerun_every = 50, rerun_random = TRUE)
samples <- monty_sample(posterior, sampler, initial = initial_values,
                        n_steps = 5000, n_chains = 4, runner = runner)
write_rds(samples, "samples-17.rds")


########################
# 18
########################
samples <- readRDS("samples-17.rds")
samples_df <- posterior::as_draws_df(samples)

vcv <- cov(samples_df[, 1:(ncol(samples_df)-3)])
mcmc_list <- coda::as.mcmc.list(samples)
if((1 - coda::rejectionRate(mcmc_list))[1] < 0.01) { vcv <- vcv * 0.01 }

initial_values <- posterior::summarise_draws(samples_df) %>% dplyr::pull(mean)
sampler <- monty_sampler_random_walk(vcv, rerun_every = 50, rerun_random = TRUE)
samples <- monty_sample(posterior, sampler, initial = initial_values,
                        n_steps = 5000, n_chains = 4, runner = runner)
write_rds(samples, "samples-18.rds")


########################
# 19
########################
samples <- readRDS("samples-18.rds")
samples_df <- posterior::as_draws_df(samples)

vcv <- cov(samples_df[, 1:(ncol(samples_df)-3)])
mcmc_list <- coda::as.mcmc.list(samples)
if((1 - coda::rejectionRate(mcmc_list))[1] < 0.01) { vcv <- vcv * 0.01 }

initial_values <- posterior::summarise_draws(samples_df) %>% dplyr::pull(mean)
sampler <- monty_sampler_random_walk(vcv, rerun_every = 50, rerun_random = TRUE)
samples <- monty_sample(posterior, sampler, initial = initial_values,
                        n_steps = 5000, n_chains = 4, runner = runner)
write_rds(samples, "samples-19.rds")


########################
# 20
########################
samples <- readRDS("samples-19.rds")
samples_df <- posterior::as_draws_df(samples)

vcv <- cov(samples_df[, 1:(ncol(samples_df)-3)])
mcmc_list <- coda::as.mcmc.list(samples)
if((1 - coda::rejectionRate(mcmc_list))[1] < 0.01) { vcv <- vcv * 0.01 }

initial_values <- posterior::summarise_draws(samples_df) %>% dplyr::pull(mean)
sampler <- monty_sampler_random_walk(vcv, rerun_every = 50, rerun_random = TRUE)
samples <- monty_sample(posterior, sampler, initial = initial_values,
                        n_steps = 5000, n_chains = 4, runner = runner)
write_rds(samples, "samples-20.rds")



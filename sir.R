
# gen_sir <- odin2::odin({

update(N[]) <- S[i] + I[i] + R[i]
update(S[]) <- S[i] + n_birth[i] - n_SI_event[i]
update(I[]) <- I[i] + n_SI[i] - n_IR_event[i]
update(R[]) <- R[i] + n_IR[i] - n_R_death[i]
update(incidence[]) <- incidence[i] + n_SI[i]
# update(X) <- N[1] - (I[1] + I[2] + I[3])


# transition probabilities of (death and I) from S

beta_0[1] <- beta_0_1
beta_0[2] <- c1*beta_0_1
beta_0[3] <- c2*beta_0_1

beta[] <- if (holiday == 1) h * beta_0[i] else if (time >= covid_start && time <= covid_end) cr * beta_0[i] else beta_0[i]

foi[] <- if (beta[i] * (I[i] + w) / N[i] < 0) 0 else beta[i] * (I[i] + w) / N[i]

# from S
p_SI_event[] <- 1 - exp(-(foi[i] + death_rate) * dt)  
p_SI_death[] <- death_rate / (foi[i] + death_rate)

n_SI_event[] <- if (S[i] < 1) 0 else Binomial(S[i], p_SI_event[i])
n_SI_death[] <- if (n_SI_event[i] < 1) 0 else Binomial(n_SI_event[i], p_SI_death[i])
n_SI[] <- n_SI_event[i] - n_SI_death[i]

# from I
p_IR_event <- 1 - exp(-(gamma + death_rate) * dt)  
p_IR_death <- death_rate / (death_rate + gamma)

n_IR_event[] <- if (I[i] < 1) 0 else Binomial(I[i], p_IR_event)
n_IR_death[] <- if (n_IR_event[i] < 1) 0 else Binomial(n_IR_event[i], p_IR_death)
n_IR[] <- n_IR_event[i] - n_IR_death[i]


# from R 
n_R_death[] <- if (R[i] < 1) 0 else Binomial(R[i], 1 - exp(-death_rate * dt)) 


# from N
n_birth[] <- if (N[i] < 1) 0 else Binomial(N[i], 1 - exp(-birth_rate * dt))


# Initial values
s_0[1] <- s_0_1
s_0[2] <- s_0_2
s_0[3] <- s_0_3
i_0[1] <- i_0_1
i_0[2] <- i_0_2
i_0[3] <- i_0_3

initial(N[]) <- N_0
initial(S[]) <- floor(N_0 * s_0[i])
initial(I[]) <- floor(N_0 * i_0[i])
initial(R[]) <- N_0 - floor(N_0 * s_0[i]) - floor(N_0 * i_0[i])
initial(incidence[], zero_every = 1) <- 0   
# initial(X, zero_every = 1) <- 0


gamma <- parameter()      
N_0 <- parameter()        
w <- parameter() 

beta_0_1 <- parameter()     
s_0_1 <- parameter()        
i_0_1 <- parameter()        
         
c1 <- parameter()
# beta_0_2 <- parameter()     
s_0_2 <- parameter()        
i_0_2 <- parameter()        

c2 <- parameter()
# beta_0_3 <- parameter()     
s_0_3 <- parameter()        
i_0_3 <- parameter()        

birth_rate <- interpolate(birth_time, birth_value, "spline")
birth_time <- parameter(constant = TRUE)
birth_value <- parameter(constant = TRUE)
dim(birth_time, birth_value) <- parameter(rank = 1)

death_rate <- interpolate(death_time, death_value, "spline")
death_time <- parameter(constant = TRUE)
death_value <- parameter(constant = TRUE)
dim(death_time, death_value) <- parameter(rank = 1)

chi <- interpolate(chi_time, chi_value, "constant")
chi_time <- parameter(constant = TRUE)
chi_value <- parameter(constant = TRUE)
dim(chi_time, chi_value) <- parameter(rank = 1)

holiday_series <- interpolate(holiday_time, holiday_value, "constant")
holiday_time <- parameter(constant = TRUE)
holiday_value <- parameter(constant = TRUE)
dim(holiday_time, holiday_value) <- parameter(rank = 1)

holiday <- if (holiday_series > 0.5) 1 else 0
cr <- 1 / (1 + exp(-(r0 - r1 * chi)))


p_1 <- parameter()
# p_2 <- parameter()
# p_3 <- parameter()
# p_4 <- parameter()

covid_start <- parameter()
covid_end   <- parameter()
r0 <- parameter()
r1 <- parameter()

h <- parameter()
# se <- parameter()

noise <- Exponential(rate = 1e6)

#- observational process -#
#- hfmd
hfmd_cases <- data()
modelled_cases <- incidence[1] * p_1 + incidence[2] * p_1 + incidence[3] * p_1
hfmd_cases ~ Poisson(modelled_cases)
#- serotype
prop[1] <- (incidence[1] * p_1  + noise) / (modelled_cases + noise)
prop[2] <- (incidence[2] * p_1  + noise) / (modelled_cases + noise)
prop[3] <- 1 - (prop[1] + prop[2])
# prop[3] <- (incidence[3] * p_1  + noise) / (modelled_cases + noise)
# prop[4] <- 1 - (prop[1] + prop[2] + prop[3])

cases <- data()
dim(cases) <- 3
cases[1:(length(cases) - 1)] ~ Binomial(sum(cases[i:length(cases)]), prop[i] / sum(prop[i:length(cases)]))




n_sero <- parameter(3)
dim(N) <- n_sero
dim(S, n_birth, n_SI_event, n_SI_death, n_SI, p_SI_event, p_SI_death) <- n_sero
dim(I, n_IR_event, n_IR_death, n_IR) <- n_sero
dim(R, n_R_death) <- n_sero
dim(beta, foi) <- n_sero
dim(beta_0, s_0, i_0, p) <- n_sero
dim(incidence) <- n_sero
dim(prop) <- n_sero

# })


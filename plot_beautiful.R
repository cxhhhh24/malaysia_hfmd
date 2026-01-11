options(scipen = 10)
## ========= 依赖 =========
library(dplyr)
library(tidyr)
library(lubridate)
library(reshape2)
library(ggplot2)
library(patchwork)
## ========= 读入并抽稀 =========
samples <- readRDS("samples-14.rds")
samples_df <- posterior::as_draws_df(samples)
posterior::summarise_draws(samples_df)

posterior::summarise_draws(samples_df)->tb
tb[, c(1,2,6,7)] %>% mutate(posterior = paste0(round(mean,3),"(",round(q5,3),"-",round(q95,3),")")) %>%
select(variable, posterior)

bayesplot::mcmc_trace(samples_df)
mcmc_list <- coda::as.mcmc.list(samples)
1 - coda::rejectionRate(mcmc_list)

samples <- monty_samples_thin(samples, burnin = 2000, thinning_factor = 8)

trj  <- samples$observations$trajectories  # [state=16, time=T, step=S, chain=C]
pars <- samples$pars                       # [param=P,  step=S,  chain=C]
stopifnot(identical(dim(trj)[3:4], dim(pars)[2:3]))
Tn <- dim(trj)[2]; S <- dim(trj)[3]; C <- dim(trj)[4]

## ========= A) overall =========

data1 <- data %>%
  ungroup() %>%
  mutate(cases = lapply(cases, function(v) { v <- as.integer(v); length(v) <- 3; v[is.na(v)] <- 0L; v })) %>%
  unnest_wider(cases, names_sep = "_") %>%
  rename(cv16 = cases_1, ev71 = cases_2, other_ev = cases_3)

par_names <- dimnames(pars)[[1]]
# if (is.null(par_names)) {
#   par_names <- c('beta_0_1','s_0_1','i_0_1','p_1',
#                  'beta_0_2','s_0_2','i_0_2','p_2',
#                  'beta_0_3','s_0_3','i_0_3','p_3',
#                  'r0','r1','h','se')
# }
i_p1 <- match("p_1", par_names) 
# i_p2 <- match("p_2", par_names)
# i_p3 <- match("p_3", par_names)
# i_se <- match("se",  par_names)
# stopifnot(all(!is.na(c(i_p1,i_p2,i_p3,i_se))))

p1 <- c(pars[i_p1, , ])    # 长度 = S*C
# p2 <- c(pars[i_p2, , ])
# p3 <- c(pars[i_p3, , ])
# se <- c(pars[i_se, , ])

# 状态顺序：N(1:3), S(1:3), I(1:3), R(1:3), incidence(1:3)
# matrix(trj[, 1, 1, 4])
inc1_mat <- matrix(trj[13, , , ], nrow = Tn, ncol = S*C)
inc2_mat <- matrix(trj[14, , , ], nrow = Tn, ncol = S*C)
inc3_mat <- matrix(trj[15, , , ], nrow = Tn, ncol = S*C)




# modelled_cases = inc1*p1 + inc2*p2 + inc3*p3 （p4=0，忽略 X）
yhat_time_draws <- sweep(inc1_mat, 2, p1, "*") +
  sweep(inc2_mat, 2, p1, "*") +
  sweep(inc3_mat, 2, p1, "*")      # [time, draws]
yhat_draws <- t(yhat_time_draws)                    # [draws, time]

yhat_long <- reshape2::melt(yhat_draws)   # Var1 = draw, Var2 = time索引, value = cases
colnames(yhat_long) <- c("draw", "time_idx", "cases")
yhat_long$time <- data1$time[yhat_long$time_idx]
start_date <- as.Date("2017-01-01")
data1$date <- start_date + weeks(data1$time - 1)
yhat_long$date <- start_date + weeks(yhat_long$time - 1)





covid_start <- as.Date("2020-07-01")
covid_end   <- as.Date("2021-11-01")

yhat_covid  <- yhat_long %>%
  filter(date >= covid_start, date <= covid_end)

data1_covid <- data1 %>%
  filter(date >= covid_start, date <= covid_end)

p_main <- ggplot() +
  annotate(
    "rect",
    xmin  = as.Date("2020-03-18"),
    xmax  = as.Date("2022-04-10"),
    ymin  = -Inf,
    ymax  = Inf,
    alpha = 0.4,
    fill  = "#FDB462"
  )+
  geom_line(
    data = yhat_long,
    aes(x = date, y = cases, group = draw, color = "Particle trajectories"),
    alpha = 0.3, linewidth = 1
  ) +
  geom_point(
    data = data1,
    aes(x = date, y = hfmd_cases, color = "Observed data"),
    size = 0.6
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Particle trajectories" = "grey60",
               "Observed data"        = "red")
  ) +
  labs(x = "Year", y = "Weekly HFMD cases") +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = c(0.01, 0.01)
  ) +
  annotate("text",
           x = as.Date("2021-03-15"),
           y = max(data1$hfmd_cases, na.rm = TRUE) * 0.45,
           label = "COVID-19 Pandemic") +
  scale_y_continuous(
    limits = c(0, 20000),
    breaks = seq(0, 20000, 5000)
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position  = "top",
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(vjust = 1, hjust = -0.5)
  )

ggplot() +
  geom_line(
    data = yhat_covid,
    aes(x = date, y = cases, group = draw, color = "Particle trajectories"),
    alpha = 0.3, linewidth = 0.6
  ) +
  geom_point(
    data = data1_covid,
    aes(x = date, y = hfmd_cases, color = "Observed data"),
    size = 1
  ) +
  geom_line(
    data = data1_covid,
    aes(x = date, y = hfmd_cases, color = "Observed data"),
    linewidth = 0.75
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Particle trajectories" = "grey60",
               "Observed data"        = "red")
  ) +
  scale_x_date(
    limits = c(covid_start, covid_end),
    breaks = as.Date(c("2020-07-01", "2021-01-01")), 
    labels = c("2020", "2021"),                       
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 250),                 
    breaks = seq(0, 250, 50)
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(hjust = -1, vjust = 1, size = 8),
    axis.text.y     = element_text(size = 8),
    plot.margin     = margin(1, 1, 1, 1),
    panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    axis.line       = element_blank()
  ) -> p_covid


p_main +
  inset_element(
    p_covid,
    left   = 0.1,   
    bottom = 0.3,
    right  = 0.5,
    top    = 0.95,
    align_to = "panel"
  )

## ========= B) serotype =========

base1 <- sweep(inc1_mat, 2, p1, "*")   # [time, draws]
base2 <- sweep(inc2_mat, 2, p1, "*")
base3 <- sweep(inc3_mat, 2, p1, "*")
den   <- base1 + base2 + base3
# den   <- pmax(den, .Machine$double.eps)  # 防 0 除


# num1 <- sweep(base1, 2, se, "*")
# num2 <- sweep(base2, 2, se, "*")
# num3 <- sweep(base3, 2, se, "*")

prop1_mat <- base1 / den
prop2_mat <- base2 / den
prop3_mat <- base3 / den
# prop3_mat <- 1 - prop1_mat - prop2_mat
# sum123    <- prop1_mat + prop2_mat + prop3_mat
# prop4_mat <- pmax(1 - sum123, 0)  # 兜底非负

# 数量：比例 × 每周 typing 总数
typing_total <- with(data1, cv16 + ev71 + other_ev)     # [time]
tot_mat      <- matrix(typing_total, nrow = Tn, ncol = ncol(prop1_mat))

y1_draws <- t(prop1_mat * tot_mat)   # CV16 [draws, time]
y2_draws <- t(prop2_mat * tot_mat)   # EV71
y3_draws <- t(prop3_mat * tot_mat)   # Other
# y4_draws <- t(prop4_mat * tot_mat)   # Neg

# 1) 展开每组 draws
y1_long <- melt(y1_draws, varnames = c("draw", "time_idx"), value.name = "count")
y1_long$category <- "CV16"
y1_long$time <- data1$time[y1_long$time_idx]

y2_long <- melt(y2_draws, varnames = c("draw", "time_idx"), value.name = "count")
y2_long$category <- "EV71"
y2_long$time <- data1$time[y2_long$time_idx]

y3_long <- melt(y3_draws, varnames = c("draw", "time_idx"), value.name = "count")
y3_long$category <- "Other EVs"
y3_long$time <- data1$time[y3_long$time_idx]

# y4_long <- melt(y4_draws, varnames = c("draw", "time_idx"), value.name = "count")
# y4_long$category <- "Negative"
# y4_long$time <- data1$time[y4_long$time_idx]

traj_long <- rbind(y1_long, y2_long, y3_long)

# 2) 观测数据长表
obs_long <- data.frame(
  time     = rep(data1$time, times = 3),
  count    = c(data1$cv16, data1$ev71, data1$other_ev),
  category = rep(c("CV16", "EV71", "Other EVs"), each = nrow(data1))
)
start_date <- as.Date("2017-01-01")

# 给 traj_long 和 obs_long 增加日期列
traj_long$date <- start_date + weeks(traj_long$time - 1)
obs_long$date  <- start_date + weeks(obs_long$time  - 1)

unique(traj_long$category)

traj_long$category <- factor(traj_long$category, levels = c("CV16","EV71","Other EVs"))
obs_long$category <- factor(obs_long$category, levels = c("CV16","EV71","Other EVs"))


#- 易感人数
## 状态顺序：N(1:3), S(1:3), I(1:3), R(1:3), incidence(1:3)
##          => S1 = index 4, S2 = 5, S3 = 6

S1_mat <- matrix(trj[4, , , ], nrow = Tn, ncol = S*C)  # CV16 的 S
S2_mat <- matrix(trj[5, , , ], nrow = Tn, ncol = S*C)  # EV71
S3_mat <- matrix(trj[6, , , ], nrow = Tn, ncol = S*C)  # Other EVs

S1_draws <- t(S1_mat)  # [draw, time]
S2_draws <- t(S2_mat)
S3_draws <- t(S3_mat)

library(reshape2)

S1_long <- melt(S1_draws, varnames = c("draw", "time_idx"), value.name = "S")
S1_long$category <- "CV16"
S1_long$time     <- data1$time[S1_long$time_idx]

S2_long <- melt(S2_draws, varnames = c("draw", "time_idx"), value.name = "S")
S2_long$category <- "EV71"
S2_long$time     <- data1$time[S2_long$time_idx]

S3_long <- melt(S3_draws, varnames = c("draw", "time_idx"), value.name = "S")
S3_long$category <- "Other EVs"
S3_long$time     <- data1$time[S3_long$time_idx]

sus_traj_long <- rbind(S1_long, S2_long, S3_long)

sus_traj_long$category <- factor(sus_traj_long$category,
                                 levels = c("CV16","EV71","Other EVs"))

start_date <- as.Date("2017-01-01")
sus_traj_long$date <- start_date + weeks(sus_traj_long$time - 1)


#- 双y轴 计算比例
max_cases <- max(obs_long$count, na.rm = TRUE)
max_S     <- max(sus_traj_long$S, na.rm = TRUE)

ratio <- max_cases / max_S * 0.8   # 0.8 让 S 线略低一点，不顶到最上面



# ggplot 绘图（x 改成 date）
covid_start <- as.Date("2020-03-18")
covid_end   <- as.Date("2022-04-10")

ggplot() +
  # 1) 阴影：COVID 时间段
  annotate(
    "rect",
    xmin  = covid_start,
    xmax  = covid_end,
    ymin  = -Inf,
    ymax  = Inf,
    alpha = 0.4,
    fill  = "#FDB462"
  ) +

  geom_line(
    data = traj_long,
    aes(x = date, y = count,
        group = interaction(category, draw),
        colour = "Particle trajectories"),
    alpha = 0.05, linewidth = 0.3
  ) +

  geom_point(
    data = obs_long,
    aes(x = date, y = count, colour = "Observed data"),
    size = 0.8
  ) +
  geom_line(
    data = obs_long,
    aes(x = date, y = count, colour = "Observed data"),
    linewidth = 0.4
  ) +
  geom_line(
    data = sus_traj_long,
    aes(x = date, y = S * ratio,
        group = interaction(category, draw),
        colour = "Simulated Susceptible"),
    alpha = 0.05, linewidth = 0.3
  ) +
  geom_text(
    data = data.frame(
      category = factor(c("CV16","EV71","Other EVs"),
                        levels = c("CV16","EV71","Other EVs")),
      x = as.Date("2021-03-01"),
      y = Inf
    ),
    aes(x = x, y = y, label = "COVID-19 Pandemic"),
    vjust = 7, size = 3
  ) +
  facet_wrap(~ category, scales = "free_y", ncol = 1) +
  scale_colour_manual(
    name = NULL,
    values = c(
      "Particle trajectories"   = "grey60",
      "Observed data"            = "red",
      "Simulated Susceptible" = "#80B1D3"
    )
  ) +
  # 8) x 轴按年份
  scale_x_date(
    limits = c(as.Date("2017-01-01"), as.Date("2022-12-31")),
    breaks = seq(as.Date("2017-01-01"), as.Date("2022-01-01"), by = "1 year"),
    labels = function(x) format(x, "%Y"),
    expand = c(0.01, 0.01)
  )+
  scale_y_continuous(
    name   = "Weekly serotyped cases",
    expand = c(0, 5),
    sec.axis = sec_axis(~ . / ratio,
                        name = "Simulated Susceptible")
  ) +
  labs(x = "Year") +
  theme_classic(base_size = 12) +
  theme(
    legend.position    = "top",
    panel.grid.minor   = element_blank(),
    strip.background   = element_rect(fill = "white"),
    strip.text         = element_text(face = "bold"),
    axis.title.y.right = element_text(margin = margin(l = 5)),
    axis.text.x      = element_text(vjust = 1, hjust = -1)
    ) +
  guides(
    colour = guide_legend(
      override.aes = list(
        alpha    = 1,              # legend 里的线不透明
        linewidth = c(0.7, 0.7, 0.7),
        linetype  = c("solid","solid","solid"),
        shape     = c(NA, 16, NA)  # 给 Observed 一个点，其它只是线
      )
    )
  )


## ========= C) I compartment (Infected) =========


# 1) 从 trj 导出 incidence 13-15（按血清型）
# inc1_mat <- matrix(trj[13, , , ], nrow = Tn, ncol = S*C)  # CV16 incidence
# inc2_mat <- matrix(trj[14, , , ], nrow = Tn, ncol = S*C)  # EV71
# inc3_mat <- matrix(trj[15, , , ], nrow = Tn, ncol = S*C)  # Other EVs

inc1_draws <- t(inc1_mat)  # [draw, time]
inc2_draws <- t(inc2_mat)
inc3_draws <- t(inc3_mat)

inc1_long <- melt(inc1_draws, varnames = c("draw","time_idx"), value.name = "inc")
inc1_long$category <- "CV16"
inc1_long$time <- data1$time[inc1_long$time_idx]

inc2_long <- melt(inc2_draws, varnames = c("draw","time_idx"), value.name = "inc")
inc2_long$category <- "EV71"
inc2_long$time <- data1$time[inc2_long$time_idx]

inc3_long <- melt(inc3_draws, varnames = c("draw","time_idx"), value.name = "inc")
inc3_long$category <- "Other EVs"
inc3_long$time <- data1$time[inc3_long$time_idx]

inc_traj_long <- rbind(inc1_long, inc2_long, inc3_long)
inc_traj_long$category <- factor(inc_traj_long$category,
                                 levels = c("CV16","EV71","Other EVs"))

# 2) date 列（和你原来一致）
start_date <- as.Date("2017-01-01")
inc_traj_long$date <- start_date + weeks(inc_traj_long$time - 1)

# 3) 右轴缩放比例：把 incidence 缩放到左轴的 cases 量级上画
max_cases <- max(obs_long$count, na.rm = TRUE)
max_inc   <- max(inc_traj_long$inc, na.rm = TRUE)
max_inc   <- pmax(max_inc, .Machine$double.eps)

ratio_inc <- max_cases / max_inc * 0.8

# 4) 画图（基本照抄你那段，只把 sus_traj_long 换成 inc_traj_long）
ggplot() +
  annotate(
    "rect",
    xmin  = covid_start,
    xmax  = covid_end,
    ymin  = -Inf,
    ymax  = Inf,
    alpha = 0.4,
    fill  = "#FDB462"
  ) +
  geom_line(
    data = traj_long,
    aes(x = date, y = count,
        group = interaction(category, draw),
        colour = "Particle trajectories"),
    alpha = 0.05, linewidth = 0.3
  ) +
  geom_point(
    data = obs_long,
    aes(x = date, y = count, colour = "Observed data"),
    size = 0.8
  ) +
  geom_line(
    data = obs_long,
    aes(x = date, y = count, colour = "Observed data"),
    linewidth = 0.4
  ) +
  geom_line(
    data = inc_traj_long,
    aes(x = date, y = inc * ratio_inc,
        group = interaction(category, draw),
        colour = "Simulated incidence"),
    alpha = 0.05, linewidth = 0.3
  ) +
  geom_text(
    data = data.frame(
      category = factor(c("CV16","EV71","Other EVs"),
                        levels = c("CV16","EV71","Other EVs")),
      x = as.Date("2021-03-01"),
      y = Inf
    ),
    aes(x = x, y = y, label = "COVID-19 Pandemic"),
    vjust = 7, size = 3
  ) +
  facet_wrap(~ category, scales = "free_y", ncol = 1) +
  scale_colour_manual(
    name = NULL,
    values = c(
      "Particle trajectories" = "grey60",
      "Observed data"         = "red",
      "Simulated incidence"   = "#80B1D3"
    )
  ) +
  scale_x_date(
    limits = c(as.Date("2017-01-01"), as.Date("2022-12-31")),
    breaks = seq(as.Date("2017-01-01"), as.Date("2022-01-01"), by = "1 year"),
    labels = function(x) format(x, "%Y"),
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    name   = "Weekly serotyped cases",
    expand = c(0, 5),
    sec.axis = sec_axis(~ . / ratio_inc,
                        name = "Simulated incidence (weekly)")
  ) +
  labs(x = "Year") +
  theme_classic(base_size = 12) +
  theme(
    legend.position    = "top",
    panel.grid.minor   = element_blank(),
    strip.background   = element_rect(fill = "white"),
    strip.text         = element_text(face = "bold"),
    axis.title.y.right = element_text(margin = margin(l = 5)),
    axis.text.x        = element_text(vjust = 1, hjust = -1)
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(
        alpha     = 1,
        linewidth = c(0.7, 0.7, 0.7),
        linetype  = c("solid","solid","solid"),
        shape     = c(NA, 16, NA)
      )
    )
  )

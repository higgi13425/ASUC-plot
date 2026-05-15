#load libraries and functons
library(tidyverse)
library(flowchart)

source("analysis-functions.R")

# load data
data = read.csv("ASUC-cohort-2024-12-10.csv")

hdata <-read.csv("higgins_data.csv") #through 2025

  

# data preprocessing ----------------------------------------------------------------

## fix date variables

data = data %>%
  mutate_at("ADMIT_DATETIME", function(x) as.Date(x, format = "%Y-%m-%d %H:%M"))

hdata = hdata %>%
  mutate_at("ADMIT_DATETIME", function(x) as.Date(x, format = "%m/%d/%Y %H:%M"))
## control for multiple admissions -> only consider first admission

first_admit_data = data %>%
  group_by(PATIENT_ID) %>%
  filter(ADMIT_DATETIME == min(ADMIT_DATETIME)) %>%
  ungroup()
# 780 unique admissions
# 
first_admit_data_h = hdata %>%
  group_by(PATIENT_ID) %>%
  filter(ADMIT_DATETIME == min(ADMIT_DATETIME)) %>%
  ungroup()
# 842 unique admissions

#  Colectomy rates by year, 2014-2024
colectomy_rates = data %>%
mutate(ADMIT_YEAR = lubridate::year(ADMIT_DATETIME)) %>%
  group_by(ADMIT_YEAR) %>%
  summarize(
    n_colectomy = sum(YN_COLECTOMY_90_DAY == 1),
    count = n(),
    colectomy_rate = n_colectomy / count,
  ) %>%
  mutate_at("colectomy_rate", function(x) round(x, digits = 3)) %>%
  select(ADMIT_YEAR, colectomy_rate, n_colectomy, count)

colectomy_rates = hdata %>%
  mutate(ADMIT_YEAR = lubridate::year(ADMIT_DATETIME)) %>%
  group_by(ADMIT_YEAR) %>% 
  filter(!is.na(ADMIT_YEAR)) |> 
  summarize(
    n_colectomy = sum(YN_COLECTOMY_90_DAY == 1),
    count = n(),
    colectomy_rate = n_colectomy / count,
  ) %>%
  mutate_at("colectomy_rate", function(x) round(x, digits = 3)) %>%
  select(ADMIT_YEAR, colectomy_rate, n_colectomy, count)

print(colectomy_rates)
write.csv(colectomy_rates, "colectomy-90-rates.csv", row.names = FALSE)

## prep colectomy rates - now in colect_rate_pct
colectomy_rates <- colectomy_rates |> rowwise() |>
  mutate(colect_rate_pct = n_colectomy/count * 100) |> mutate(lcb = 100*Hmisc :: binconf(n_colectomy, count, alpha = 0.05, method = "wilson")[2], ucb = 100*Hmisc :: binconf(n_colectomy, count, alpha = 0.05, method = "wilson")[3]) |>
  mutate(colect_rate_pct = n_colectomy/count * 100)

## plot colectomy rates / year

ASUC_plot <- colectomy_rates %>%
  ggplot(aes(x = ADMIT_YEAR, y = colect_rate_pct)) +
  annotate("rect", xmin = 2019.75, xmax = 2022.25, ymin = 5, ymax = 9, alpha = .1,fill = "red") +
  annotate("text", x = 2021, y =7, label = "COVID Times") +
  geom_line(color = "red", linewidth=1) +
  geom_errorbar(aes(ymin=lcb, ymax=ucb), width=.3) +
  geom_label(aes(label = round(colect_rate_pct,1), size = 44)) +
  geom_text(aes(x = ADMIT_YEAR, y = -.3, label = paste0(n_colectomy, "/", count)), size = 4, vjust = 1.5) +
  scale_x_continuous(n.breaks=12,limits = c(2013.7,2025.3)) +
  xlab("Admission Year") + ylab("90-Day Colectomy Percentage") +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(-3,42)) +
  annotate("text", x = 2014.5, y =3, label = "N Colect /\nN ASUC") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none")

ASUC_plot

ggsave(ASUC_plot, filename = "ASUC_plot1.png", width = 10, height = 4.5, units = "in", dpi = 300)


## dataset with colectomy rates and prior therapy rates by ADMIT_YEAR
## calc percent adv Rx
data2 <- hdata |>
  mutate(ADMIT_YEAR = lubridate::year(ADMIT_DATETIME)) |>  replace_na(list(HNP_PRIOR_ADALIMUMAB=0, HNP_PRIOR_INFLIXIMAB=0, HNP_PRIOR_GOLIMUMAB=0, HNP_PRIOR_VEDOLIZUMAB=0, HNP_PRIOR_USTEKINUMAB=0, HNP_PRIOR_TOFACITINIB=0, HNP_PRIOR_UPADACITINIB=0)) |>
  rowwise() |>
  mutate(prior_biol = sum(HNP_PRIOR_ADALIMUMAB, HNP_PRIOR_GOLIMUMAB, HNP_PRIOR_INFLIXIMAB, HNP_PRIOR_USTEKINUMAB, HNP_PRIOR_VEDOLIZUMAB)) |> 
  mutate(prior_jak = sum(HNP_PRIOR_TOFACITINIB, HNP_PRIOR_UPADACITINIB)) |> 
  mutate(prior_biol = ifelse(prior_biol > 0, 1, 0)) |>
  mutate(prior_jak = ifelse(prior_jak > 0, 1, 0)) |>
  mutate(group = case_when(
    prior_biol > 0 & prior_jak > 0 ~ "both",
    prior_biol == 1 ~ "biol",
    prior_jak == 1 ~ "jak",
    TRUE ~ "none")) |>
  ungroup() |>
  select(ADMIT_YEAR, prior_biol, prior_jak, group) |> group_by(ADMIT_YEAR) |>
  summarise(pct_adv_rx = round(100*sum(group != "none")/n(),2)) |> 
  right_join(colectomy_rates) 

# plot wth PCT prior advanced Rx
colect_adv1 <- data2 |>
  ggplot(aes(x = ADMIT_YEAR, y = colect_rate_pct)) +
  geom_area(aes(x = ADMIT_YEAR, y = pct_adv_rx), fill = "cadetblue", alpha = 0.2) +
  # annotate("rect", xmin = 2014.5, xmax = 2018.5, ymin = 1, ymax = 5, alpha = .1,fill = "green") +
  # annotate("text", x = 2016.5, y =3, label = "Accelerated IFX Era") +
  # annotate("rect", xmin = 2018.5, xmax = 2022.5, ymin = 1, ymax = 5, alpha = .1,fill = "blue") +
  # annotate("text", x = 2020.8, y =3, label = "Tofa Era") +
  # annotate("rect", xmin = 2022.5, xmax = 2024, ymin = 1, ymax = 5, alpha = .3,fill = "orange") +
  # annotate("text", x = 2023.1, y =3, label = "Upa Era") +
  annotate("rect", xmin = 2019.75, xmax = 2022.25, ymin = 5, ymax = 9, alpha = .4,fill = "red") +
  annotate("text", x = 2021, y =7, label = "COVID Times") +
  geom_line(aes( y = colect_rate_pct), color = "red", linewidth=1) +
  geom_errorbar(aes(ymin=lcb, ymax=ucb), width=.2) +
  geom_label(aes( y = colect_rate_pct, label = round(colect_rate_pct,1), size = 44)) +
  geom_text(aes(x = ADMIT_YEAR, y = -.3, label = paste0(n_colectomy, "/", count)), size = 4, vjust = 1.5) +
  scale_x_continuous(n.breaks=10,limits = c(2013.9,2025.1)) +
  xlab("Admission Year") + ylab("90-Day Colectomy Percentage") +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(-3,65)) +
  annotate("text", x = 2014.3, y =3, label = "N Colect /\nN ASUC") +
  annotate("text", x = 2019.3, y =45, label = "Percent Prior\nAdvanced Therapy", size = 6) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none")

colect_adv1

ggsave(colect_adv1, filename = "ASUC_plot2.png", width = 10, height = 4.5, units = "in", dpi = 300)


# plot adding treatment eras
colect_adv2 <- data2 |>
  ggplot(aes(x = ADMIT_YEAR, y = colect_rate_pct)) +
  geom_area(aes(x = ADMIT_YEAR, y = pct_adv_rx), fill = "cadetblue", alpha = 0.2) +
  annotate("rect", xmin = 2014.5, xmax = 2018.5, ymin = 1, ymax = 5, alpha = .1,fill = "green") +
  annotate("text", x = 2016.5, y =3, label = "Accelerated IFX Era") +
  annotate("rect", xmin = 2018.5, xmax = 2022.5, ymin = 1, ymax = 5, alpha = .1,fill = "blue") +
  annotate("text", x = 2020.5, y =3, label = "Tofa Era") +
  annotate("rect", xmin = 2022.5, xmax = 2025.1, ymin = 1, ymax = 5, alpha = .3,fill = "orange") +
  annotate("text", x = 2023.4, y =3, label = "Upa") +
  annotate("text", x = 2024.6, y =3, label = "Era") +
  annotate("rect", xmin = 2019.75, xmax = 2022.25, ymin = 5, ymax = 9, alpha = .4,fill = "red") +
  annotate("text", x = 2021, y =7, label = "COVID Times") +
  geom_line(aes( y = colect_rate_pct), color = "red", linewidth=1) +
  geom_errorbar(aes(ymin=lcb, ymax=ucb), width=.2) +
  geom_label(aes( y = colect_rate_pct, label = round(colect_rate_pct,1), size = 44)) +
  geom_text(aes(x = ADMIT_YEAR, y = -.3, label = paste0(n_colectomy, "/", count)), size = 4, vjust = 1.5) +
  scale_x_continuous(n.breaks=10,limits = c(2013.9,2025.1)) +
  xlab("Admission Year") + ylab("90-Day Colectomy Percentage") +
  labs(caption = "IFX = Infliximab, Tofa = Tofacitinib, Upa = Upadacitinib\nHigher Percent Prior Advanced Therapy = Harder to Treat") +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(-3,65)) +
  annotate("text", x = 2014, y =3, label = "N Colect /\nN ASUC", size = 4.5) +
  annotate("text", x = 2019.3, y =45, label = "Percent Prior\nAdvanced Therapy", size = 6) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none")

colect_adv2

ggsave(colect_adv2, filename = "ASUC_plot3.png", width = 10, height = 4.5, units = "in", dpi = 300)



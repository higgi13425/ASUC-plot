# README -----------------------------------------------------------------------
#
#   file: higgins-talk.R
#   goal: complete analytic requests for Dr. Higgins talk (2024-09-01)
#
# ------------------------------------------------------------------------------

# packages and data ------------------------------------------------------------

library(tidyverse)
library(flowchart)

  source("analysis-functions.R")

  data = read.csv("ASUC-cohort-2024-12-10.csv")


# preprocessing ----------------------------------------------------------------

  ## fix date variables

    data = data %>%
      mutate_at("ADMIT_DATETIME", function(x) as.Date(x, format = "%Y-%m-%d %H:%M"))

  ## control for multiple admissions -> only consider first admission

    first_admit_data = data %>%
      group_by(PATIENT_ID) %>%
      filter(ADMIT_DATETIME == min(ADMIT_DATETIME)) %>%
      ungroup()
      # 819 unique admissions


# REQ1: IV steroids success ----------------------------------------------------

  ## questions:
  ## - What is the rate of success of IVCS in getting people to discharge without a rescue therapy?
  ## - a.	Has it changed (gone down) over time (year to year) in the biologic era?

  ## define cohort + success metric (success = YN_COLECTOMY_90_DAY)

    ivcs = data %>%
      filter(YN_INFLIXIMAB == 0 & YN_TOFACITINIB == 0 & YN_UPADACITINIB == 0 & YN_CYCLOSPORINE == 0)

  ## get success rate

    n_success = ivcs %>% filter(YN_COLECTOMY_90_DAY == 0) %>% nrow()
    n_total   = ivcs %>% nrow()

    overall_success_rate = n_success / n_total
    print(overall_success_rate)

  ## get success rates / year

    ivcs_success = ivcs %>%
      mutate(ADMIT_YEAR = lubridate::year(ADMIT_DATETIME)) %>%
      group_by(ADMIT_YEAR) %>%
      summarize(
        n_success = sum(YN_COLECTOMY_90_DAY == 0),
        count = n(),
        success_rate = sum(YN_COLECTOMY_90_DAY == 0) / count,
      ) %>%
      mutate_at("success_rate", function(x) round(x, digits = 3)) %>%
      select(ADMIT_YEAR, success_rate, n_success, count)

    print(ivcs_success)
    write.csv(ivcs_success, "all-ivcs-success-rates.csv", row.names = FALSE)

  ## plot success rates / year

    ivcs_success %>%
      ggplot(aes(x = ADMIT_YEAR, y = success_rate)) +
      geom_line(color = "blue", linewidth=1) +
      geom_label(aes(label = round(100*success_rate,2), size = 44)) +
      scale_x_continuous(n.breaks=10) +
      xlab("Admission Year") + ylab("Success Rate") +
      ggtitle("Michigan ASUC Cohort Success Rates with IVCS by Year") +
      labs(subtitle = "Success = No Colectomy in 90 Days") +
      theme_bw(base_size = 16)  +
      geom_text(aes(x = ADMIT_YEAR, y = 0.70, label = paste0(n_success, "/", count)))  +
      annotate("text", x = 2022.6, y =0.74, label = "Falling Use \nof IVCS Monotherapy", size =5) +
      annotate("rect", xmin = 2020.7, xmax = 2024.5, ymin = 0.71, ymax = 0.77, alpha = .1,fill = "red") +
      scale_y_continuous(labels = scales::percent_format(scale = 100), limits = c(0.7,1)) +
      theme(legend.position = "none")


ggsave("ivcs-success-rates.png", width = 8, height = 4.5, units = "in", dpi = 300)

# REQ2: important risk factors for colectomy -----------------------------------

  ## questions:
  ## - what are the most important risk factors for colectomy?
  ## - what are their associated odds ratios?

  ## define + scale predictors

    predictors = data %>%
      select(
        AGE, SEX, BMI,                                      # demographic variables
        YN_OUTSIDE_HOSPITAL_TRANSFER,                       # admission characteristics
        contains(c("ADMIT_", "DAY_")) & !contains("DATE")   # lab values
        & !contains(c("4", "5", "6", "7"))
      ) %>%
      colnames()

    scaled_data = data %>%
      mutate(across(contains("PLT"), function(x) x / 100))

  ## univariable analysis

    uni_results = get_LR_summary(scaled_data, predictors, "YN_COLECTOMY_90_DAY")

    print(head(uni_results, 30))
    write.csv(uni_results, "colectomy-90-univariable.csv", row.names = FALSE)


# REQ3: multivariable analysis -------------------------------------------------

  ## questions:
  ## - are admission CRP, admission albumin, and outside hospital transfer
  ##   good predictors of 90-day colectomy in a combined model?

  ## multivariable analysis

    multi_model   = glm(YN_COLECTOMY_90_DAY ~ ADMIT_CRP + ADMIT_ALB + YN_OUTSIDE_HOSPITAL_TRANSFER, data = data, family = "binomial")
    multi_results = broom::tidy(multi_model) %>%
      filter(term != "(Intercept)") %>%
      mutate(
        OR = exp(estimate) %>% round(digits = 2),
        lower_ci = estimate - 1.96 * std.error,
        upper_ci = estimate + 1.96 * std.error,
        lower_OR = exp(lower_ci) %>% round(digits = 2),
        upper_OR = exp(upper_ci) %>% round(digits = 2),
        OR_summary = glue("{OR} ({lower_OR}, {upper_OR})"),
        p = ifelse(
          p.value < 0.001, "< 0.001",
          round(p.value, digits = 3) %>% as.character()
        )
      ) %>%
      select(variable = term, OR_summary, p)

    print(multi_results)
    write.csv(multi_results, "colectomy-90-multivariable-CRP-ALB-OSH.csv", row.names = FALSE)

  ## assess VIF scores

    car::vif(multi_model)  # LGTM - all predictors well below 5


# REQ4: colectomy rates by year ------------------------------------------------

  ## questions:
  ## - what are the rates of 90-day colectomy by year across the dataset?
  ## - do they show any specific trend over time?

  ## get colectomy rates / year

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

    print(colectomy_rates)
    write.csv(colectomy_rates, "colectomy-90-rates.csv", row.names = FALSE)

## prep colectomy rates
colectomy_rates <- colectomy_rates |> rowwise() |>
  mutate(colect_rate_pct = n_colectomy/count * 100) |> mutate(lcb = 100*Hmisc :: binconf(n_colectomy, count, alpha = 0.05, method = "wilson")[2], ucb = 100*Hmisc :: binconf(n_colectomy, count, alpha = 0.05, method = "wilson")[3]) |>
  mutate(colect_rate_pct = n_colectomy/count * 100)

  ## plot colectomy rates / year

ASUC_plot <- colectomy_rates %>%
      ggplot(aes(x = ADMIT_YEAR, y = colect_rate_pct)) +
      annotate("rect", xmin = 2014.5, xmax = 2018.5, ymin = 1, ymax = 5, alpha = .1,fill = "green") +
      annotate("text", x = 2016.5, y =3, label = "Accelerated IFX Era") +
      annotate("rect", xmin = 2018.5, xmax = 2022.5, ymin = 1, ymax = 5, alpha = .1,fill = "blue") +
      annotate("text", x = 2020.8, y =3, label = "Tofa Era") +
      annotate("rect", xmin = 2022.5, xmax = 2024.5, ymin = 1, ymax = 5, alpha = .3,fill = "orange") +
      annotate("text", x = 2023.1, y =3, label = "Upa Era") +
      annotate("rect", xmin = 2019.75, xmax = 2022.25, ymin = 5, ymax = 9, alpha = .1,fill = "red") +
      annotate("text", x = 2021, y =7, label = "COVID Times") +
      geom_line(color = "red", linewidth=1) +
      geom_errorbar(aes(ymin=lcb, ymax=ucb), width=.3) +
      geom_label(aes(label = round(colect_rate_pct,1), size = 44)) +
      geom_text(aes(x = ADMIT_YEAR, y = -.3, label = paste0(n_colectomy, "/", count)), size = 4, vjust = 1.5) +
      scale_x_continuous(n.breaks=10,limits = c(2013.7,2024.3)) +
      xlab("Admission Year") + ylab("90-Day Colectomy Percentage") +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(-3,42)) +
      annotate("text", x = 2014, y =3, label = "N Colect /\nN ASUC") +
      theme_bw(base_size = 18) +
  theme(legend.position = "none")

ASUC_plot

ggsave(ASUC_plot, filename = "ASUC_plot.png", width = 10, height = 4.5, units = "in", dpi = 300)


## dataset with colectomy rates and prior therapy rates by ADMIT_YEAR
##
data2 <- data |>
  mutate(ADMIT_YEAR = lubridate::year(ADMIT_DATETIME)) |>  replace_na(list(PRIOR_ADALIMUMAB=0, PRIOR_INFLIXIMAB=0, PRIOR_GOLIMUMAB=0, PRIOR_VEDOLIZUMAB=0, PRIOR_USTEKINUMAB=0, PRIOR_TOFACITINIB=0, PRIOR_UPADACITINIB=0)) |>
  rowwise() |>
  mutate(prior_biol = sum(PRIOR_ADALIMUMAB, PRIOR_GOLIMUMAB, PRIOR_INFLIXIMAB, PRIOR_USTEKINUMAB, PRIOR_VEDOLIZUMAB)) |>
  mutate(prior_jak = sum(PRIOR_TOFACITINIB, PRIOR_UPADACITINIB)) |>
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


colect_adv <- data2 |>
  ggplot(aes(x = ADMIT_YEAR, y = colect_rate_pct)) +
  #geom_area(aes(x = ADMIT_YEAR, y = pct_adv_rx), fill = "blue", alpha = 0.2) +
  annotate("rect", xmin = 2014.5, xmax = 2018.5, ymin = 1, ymax = 5, alpha = .1,fill = "green") +
  annotate("text", x = 2016.5, y =3, label = "Accelerated IFX Era") +
  annotate("rect", xmin = 2018.5, xmax = 2022.5, ymin = 1, ymax = 5, alpha = .1,fill = "blue") +
  annotate("text", x = 2020.8, y =3, label = "Tofa Era") +
  annotate("rect", xmin = 2022.5, xmax = 2024, ymin = 1, ymax = 5, alpha = .3,fill = "orange") +
  annotate("text", x = 2023.1, y =3, label = "Upa Era") +
  annotate("rect", xmin = 2019.75, xmax = 2022.25, ymin = 5, ymax = 9, alpha = .4,fill = "red") +
  annotate("text", x = 2021, y =7, label = "COVID Times") +
  #annotate("text", x = 2019.65, y =36, label = "Percent Advanced\nRx Experienced") +
  geom_line(aes( y = colect_rate_pct), color = "red", linewidth=1) +
  geom_errorbar(aes(ymin=lcb, ymax=ucb), width=.2) +
  geom_label(aes( y = colect_rate_pct, label = round(colect_rate_pct,1), size = 44)) +
  geom_text(aes(x = ADMIT_YEAR, y = -.3, label = paste0(n_colectomy, "/", count)), size = 4, vjust = 1.5) +
  scale_x_continuous(n.breaks=10,limits = c(2013.9,2024.1)) +
  xlab("Admission Year") + ylab("90-Day Colectomy Percentage") +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(-3,50)) +
  annotate("text", x = 2014.3, y =3, label = "N Colect /\nN ASUC") +
  theme_bw(base_size = 18) +
  theme(legend.position = "none")

colect_adv

ggsave(colect_adv, filename = "ASUC_plot.png", width = 10, height = 4.5, units = "in", dpi = 300)

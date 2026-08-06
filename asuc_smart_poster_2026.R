# data from ASUC SMART Trial inpatient

library(tidyverse)
# in percentages
data <- tibble::tribble(
  ~pt_id, ~arm, ~d7resp, ~avoid_rescue, ~colect_90d,
  1, 'upa', 1, 1, 0,
  2, 'upa', 1, 1, 0,
  3, 'upa', 1, 1, 0,
  4, 'upa', 1, 1, 0,
  5, 'upa', 1, 1, 0,
  6, 'upa', 1, 1, 0,
  7, 'upa', 1, 1, 0,
  8, 'upa', 1, 0, 0,
  9, 'upa', 0, 0, 0,
  10, 'upa', 0, 0, 0,
  11, 'upa', 0, 0, 0,
  12, 'upa', 0, 0, 0,
  13, 'upa', 0, 0, 1,
  14, 'upa', 0, 0, 1,
  15, 'upa+ivcs', 1, 1, 0,
  16, 'upa+ivcs', 1, 1, 0,
  17, 'upa+ivcs', 1, 1, 0,
  18, 'upa+ivcs', 1, 1, 0,
  19, 'upa+ivcs', 1, 1, 0,
  20, 'upa+ivcs', 1, 1, 0,
  21, 'upa+ivcs', 1, 1, 0,
  22, 'upa+ivcs', 1, 1, 0,
  23, 'upa+ivcs', 1, 1, 0,
  24, 'upa+ivcs', 1, 1, 0,
  25, 'upa+ivcs', 1, 1, 0,
  26, 'upa+ivcs', 1, 1, 0,
  27, 'upa+ivcs', 1, 1, 0,
  28, 'upa+ivcs', 0, 1, 1,
  29, 'upa+ivcs', 0, 0, 1,
  30, 'ivcs', 1, 1, 0,
  31, 'ivcs', 1, 1, 0,
  32, 'ivcs', 1, 1, 0,
  33, 'ivcs', 1, 1, 0,
  34, 'ivcs', 1, 1, 0,
  35, 'ivcs', 0, 1, 0,
  36, 'ivcs', 0, 1, 0,
  37, 'ivcs', 0, 0, 0,
  38, 'ivcs', 0, 0, 0,
  39, 'ivcs', 0, 0, 0,
  40, 'ivcs', 1, 0, 0,
  41, 'ivcs', 1, 0, 1,
  42, 'ivcs', 1, 0, 1,
)

plot_data <- data %>%
  pivot_longer(
    cols = c(d7resp, avoid_rescue, colect_90d),
    names_to = "outcome",
    values_to = "value"
  ) %>%
  group_by(arm, outcome) %>%
  summarize(pct = mean(value) * 100, .groups = "drop") %>%
  mutate(
    outcome = factor(
      outcome,
      levels = c("d7resp", "avoid_rescue", "colect_90d"),
      labels = c("Day 7\nClinical Response", "Avoided Need for\nRescue Therapy", "90-Day Colectomy")
    ),
    arm = factor(arm, levels = c("ivcs", "upa", "upa+ivcs"))
  )


# 3. Create Horizontal Lollipop Plot
pd <- position_dodge(width = 0.6)

ggplot(plot_data, aes(y = outcome, color = arm)) +
  geom_linerange(
    aes(xmin = 0, xmax = pct),
    position = pd,
    linewidth = 0.8
  ) +
  geom_point(
    aes(x = pct),
    position = pd,
    size = 3.5
  ) +
  geom_text(
    aes(x = pct, label = paste0(round(pct, 1), "%")),
    position = pd,
    fontface  = "bold",
    hjust = -0.3,
    size = 3.2,
    show.legend = FALSE
  ) +
  scale_x_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, 110),
    expand = c(0, 0)
  ) +
  scale_color_manual(values = c("green4", "#2297E6", "#DF536B"),
                     name = "Treatment Arm (N)",
                     limits = c("upa+ivcs", "upa", "ivcs"),
                     labels = c(
                       "upa+ivcs" = "Upadacitinib + IV Steroids (15)",
                       "upa"      = "Upadacitinib Monotherapy (14)",
                       "ivcs"     = "IV Steroids Alone (13)")
                     ) +
  labs(
    title = "ASUC SMART Trial: Outcomes by Treatment Arm",
    subtitle = "Randomized Open-label Pilot Trial in 42 ASUC Patients at the University of Michigan",
    x = "Percentage of Patients (%)",
    y = "Outcome",
    caption = "Presented at DDW 2026, Abstract Su1605\nASUC Inpatients treated with 30 mg bid or 45 mg po qd Upa + 1 mg/kg/d IV methylprednisolone"
  ) +
  annotate('text', x = 40, y =3, label = "p > 0.9", size = 5) +
  annotate('text', x = 80, y =2, label = "p = 0.017", size = 5) +
  annotate('text', x = 80, y =1, label = "p = 0.029", size = 5) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    
    # Position legend inside top right
    legend.position = c(0.95, 0.96),
    legend.justification = c(1, 1),
    legend.direction = "vertical",
    
    # add a transparent background box behind the legend for clarity
    legend.background = element_rect(color = "black", fill = "white")
  )

ggsave('smart_trial.png', width = 8, height = 6)



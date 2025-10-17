hist_combo2 <- bind_rows(
  tibble(
    group_id          = "trial_A",
    drug_A            = c(3, 4.5, 6, 8),
    num_patients      = c(3, 3, 6, 3),
    num_toxicities    = c(0, 0, 0, 2),
    cohort_time       = 0
  ),
  tibble(
    group_id = "trial_B",
    drug_B = c(33.3, 50, 100, 200, 400, 800, 1120),
    num_patients = c(3, 3, 4, 9, 15, 20, 17),
    num_toxicities = c(0, 0, 0, 0, 0, 2, 4),
    cohort_time = 0
  )
) %>%
  replace_na(list(drug_A = 0, drug_B = 0)) %>%
  mutate(group_id = factor(group_id, c("trial_A", "trial_B", "IIT", "trial_AB"))) %>%
  select(group_id, drug_A, drug_B, num_patients, num_toxicities, cohort_time)

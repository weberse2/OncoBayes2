dose_info_combo2 <- bind_rows(
  crossing(
    group_id = "trial_A",
    drug_A = c(3, 4.5, 6, 8),
  ),
  crossing(
    group_id = c("trial_AB", "IIT"),
    drug_A = c(0, 3, 4.5, 6, 8),
    drug_B = c(0, 400, 600, 800)
  )
) %>%
  filter(!(drug_A == 0 & drug_B == 0)) %>%
  replace_na(list(drug_A = 0, drug_B = 0)) %>%
  mutate(
    group_id = factor(group_id, c("trial_A", "trial_B", "IIT", "trial_AB")),
    dose_id = 1:n()
  )

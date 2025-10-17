hist_combo3 <- bind_rows(
  tibble(
    stratum_id = "BID",
    group_id = "Combo",
    drug_A = 400,
    drug_B = 800,
    drug_C = c(80, 160, 240),
    num_toxicities = c(0, 1, 2),
    num_patients = c(4, 8, 5)
  ),
  tibble(
    stratum_id = "QD",
    group_id = "HistAgent1",
    drug_A = c(0, 0, 0, 240, 400, 400, 400),
    drug_B = c(240, 400, 800, 80, 400, 800, 1000),
    num_toxicities = c(0, 0, 1, 0, 1, 0, 1),
    num_patients = c(5, 10, 11, 4, 6, 11, 5)
  ),
  tibble(
    stratum_id = "QD",
    group_id = "HistAgent2",
    drug_A = c(0, 0, 0, 0, 0, 400, 400, 400),
    drug_C = c(80, 160, 320, 480, 640, 160, 320, 240),
    num_toxicities = c(0, 0, 0, 0, 1, 1, 2, 1),
    num_patients = c(3, 3, 6, 3, 5, 6, 5, 1)
  )
) %>%
  replace_na(list(drug_A = 0, drug_B = 0, drug_C = 0)) %>%
  mutate(
    group_id = factor(group_id),
    stratum_id = factor(stratum_id)
  )

drug_info_combo2 <- tibble(
  drug_name = c("drug_A", "drug_B"),
  ##  dose_ref  = c(    300,     960), ## originally in paper
  dose_ref = c(6, 1500), ## adapted to actual scale of doses
  dose_unit = "mg",
  ## reference_p_dlt = 0.1  ## originally in paper
  reference_p_dlt = 0.2 ## adapted to updated best practice
)

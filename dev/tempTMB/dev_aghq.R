#' Develop the aghq models here!

load("data/mw.rda")
mw <- sf::st_as_sf(mw)

system.time({constant_aghq(mw)})
system.time({iid_aghq(mw)})
system.time({besag_aghq(mw)})
system.time({bym2_aghq(mw)})
system.time({fck_aghq(mw)})
system.time({fik_aghq(mw)})
system.time({ck_aghq(mw)})
system.time({ik_aghq(mw)})

#' Test that the cross-validation enabling code works
ii <- 0
cv_test <- iid_aghq(mw, ii = ii)

cv_summary <- summary(cv_test)

fig_one_left_out <- cv_summary$randomeffectsummary %>%
  filter(variable == "u") %>%
  tibble::rownames_to_column("index") %>%
  mutate(
    index = as.numeric(index),
    left_out = case_when(
      (index - 1) %in% ii ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  ggplot(aes(x = index, y = mean, ymax = `2.5%`, ymin = `97.5%`, col = left_out)) +
    geom_pointrange() +
    theme_minimal() +
    labs(x = "Index", y = "Spatial effect", col = "Left out?")

ii <- 1:5
cv_test <- iid_aghq(mw, ii = ii)

cv_summary <- summary(cv_test)

fig_many_left_out <- cv_summary$randomeffectsummary %>%
  filter(variable == "u") %>%
  tibble::rownames_to_column("index") %>%
  mutate(
    index = as.numeric(index),
    left_out = case_when(
      (index - 1) %in% ii ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  ggplot(aes(x = index, y = mean, ymax = `2.5%`, ymin = `97.5%`, col = left_out)) +
  geom_pointrange() +
  theme_minimal() +
  labs(x = "Index", y = "Spatial effect", col = "Left out?")

fig_one_left_out / fig_many_left_out

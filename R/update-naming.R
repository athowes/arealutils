#' A dictionary to convert from internal names to human readable names, used in the manuscript.
#'
#' @param df A dataframe with columns `geometry`, `sim_model` and `inf_model`.
#' @return A dataframe where the entries in those columns have been renamed.
#' @export
update_naming <- function(df) {
  mutate(
    df,
    geometry = recode_factor(
      geometry,
      "1" = "1",
      "2" = "2",
      "3" = "3",
      "4" = "4",
      "grid" = "Grid",
      "civ" = "Cote d'Ivoire",
      "tex" = "Texas"
    ),
    sim_model = recode_factor(
      sim_model,
      "iid" = "IID",
      "icar" = "Besag",
      "ik" = "IK"
    ),
    inf_model = recode_factor(
      inf_model,
      "constant_aghq" = "Constant",
      "iid_aghq" = "IID",
      "besag_aghq" = "Besag",
      "bym2_aghq" = "BYM2",
      "fck_aghq" = "FCK",
      "ck_aghq" = "CK",
      "fik_aghq" = "FIK",
      "ik_aghq" = "IK"
    )
  )
}
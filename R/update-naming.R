#' A dictionary to convert from internal names to human readable names, used in the manuscript.
#'
#' @param df A dataframe with columns `geometry`, `sim_model` and `inf_model`.
#' @return A dataframe where the entries in those columns have been renamed.
#' @export
update_naming <- function(df) {
  mutate(
    df,
    geometry = recode_factor(geometry,
                             "grid" = "Grid",
                             "civ" = "Cote d'Ivoire",
                             "tex" = "Texas"),
    sim_model = recode_factor(sim_model,
                              "iid" = "IID",
                              "icar" = "Besag",
                              "ik" = "IK"),
    inf_model = recode_factor(inf_model,
                              "constant_inla" = "Constant",
                              "iid_inla" = "IID",
                              "besag_inla" = "Besag",
                              "bym2_inla" = "BYM2",
                              "fck_inla" = "FCK",
                              "ck_stan" = "CK",
                              "fik_inla" = "FIK",
                              "ik_stan" = "IK")
  )
}
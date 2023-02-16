#' Bind the rows of a list together, adding an `id` column.
#'
#' @param list A list.
#' @return A dataframe.
#' @export
list_to_df <- function(list){
  data.frame(dplyr::bind_rows(list, .id = "replicate"))
}

#' Convert the rows of a dataframe to a list.
#'
#' @param df A dataframe.
#' @return A list.
#' @export
rows_to_list <- function(df) {
  x <- as.list((data.frame(t(df))))
  names(x) <- NULL
  return(x)
}

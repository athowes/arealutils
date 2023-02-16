#' Create height by width grid `sf` object
#'
#' @param height Number of grid cells height.
#' @param width Number of grid cells width.
#' @return `height` by `width` `sf` object.
#' @export
create_sf_grid <- function(height, width){
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(width, 0), c(width, height), c(0, 0)))))
  grid <- sf::st_make_grid(sfc, cellsize = 1, square = TRUE)
  return(grid)
}

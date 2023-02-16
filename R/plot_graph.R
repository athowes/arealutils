#' Produce a `ggplot2` visualisation of a (neighbourhood) graph.
#'
#' @template sf
#' @param add_geography Should a map of the geography be plotted alongside the graph?
#' @return A `ggplot2` plot.
#' @examples
#' plot_graph(mw, add_geography = TRUE)
#' @export
plot_graph <- function(sf, add_geography = FALSE) {
  nb <- sf_to_nb(sf)
  
  nb_sf <- spdep::nb2lines(nb, coords = sp::coordinates(as(sf, "Spatial"))) %>%
    as("sf") %>%
    sf::st_set_crs(sf::st_crs(sf))

  b <- ggplot2::ggplot(sf) +
    ggplot2::geom_sf(data = nb_sf) +
    ggplot2::theme_minimal() +
    ggplot2::labs(subtitle = "Graph") + 
    ggplot2::theme_void()
    
  if(add_geography) {
    a <- ggplot2::ggplot(sf) +
      ggplot2::geom_sf() +
      ggplot2::theme_minimal() +
      ggplot2::labs(subtitle = "Geography") + 
      ggplot2::theme_void()
  
    return(cowplot::plot_grid(a, b, ncol = 2))  
  } else {
    return(b)
  }
}

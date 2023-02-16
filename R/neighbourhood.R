#' Wrapper function for `spdep::poly2nb`.
#'
#' @template sf
#' @return A neighbourhood list object.
#' @examples
#' nb <- sf_to_nb(mw)
#' @export
sf_to_nb <- function(sf){
  nb <- spdep::poly2nb(as(sf, "Spatial"))
  return(nb)
}

#' Wrapper function for `spdep::mat2listw`.
#'
#' @param mat An adjacency matrix.
#' @return A neighbourhood list object.
#' @examples
#' adj <- rbind(c(0, 1, 0), c(0, 0, 1), c(0, 0, 0))
#' adj <- adj + t(adj)
#' nb <- mat_to_nb(adj)
#' @export
mat_to_nb <- function(mat){
  nb <- spdep::mat2listw(mat)$neighbours
  return(nb)
}

#' Create a graph from `nb`.
#'
#' @param nb A neighbourhood list object.
#' @return A list containing
#' * `n` The number of nodes in the graph
#' * `n_edges` To-do
#' * `node1` To-do
#' * `node2` To-do
#' @source From \href{https://github.com/stan-dev/example-models/blob/master/knitr/car-iar-poisson/nb_data_funs.R}{code} by Mitzi Morris.
#' @examples
#' nb <- sf_to_nb(mw)
#' nb_to_graph(nb)
#' @export
nb_to_graph <- function(nb){
  n <- length(nb)
  n_links <- 0
  for (i in 1:n) {
    if (nb[[i]][1] > 0) {
      n_links <- n_links + length(nb[[i]])
    }
  }
  n_edges <- n_links / 2;
  node1 <- vector(mode = "numeric", length = n_edges)
  node2 <- vector(mode = "numeric", length = n_edges)
  idx <- 0
  for (i in 1:n) {
    if (nb[[i]][1] > 0) {
      for (j in 1:length(nb[[i]])) {
        n2 <- unlist(nb[[i]][j])
        if (i < n2) {
          idx <- idx + 1
          node1[idx] <- i
          node2[idx] <- n2
        }
      }
    }
  }
  return(list("n" = n, "n_edges" = n_edges, "node1" = node1, "node2" = node2))
}

#' Create ICAR precision matrix from `nb`.
#'
#' @param nb A neighbourhood list object.
#' @param W A `length(nb)` by `length(nb)` symmetric weights matrix with zeros for non-neighbouring areas.
#' @param return_sparse Should the matrix returned be of class `sparseMatrix`?
#' @param scale_precision Should the matrix returned be scaled using `scale_gmrf_precision`?
#' @examples
#' nb <- sf_to_nb(mw)
#' nb_to_precision(nb)
#' @export
nb_to_precision <- function(nb, W = NULL, return_sparse = FALSE, scale_precision = FALSE){
  n <- length(nb)
  
  if(is.null(W)) {
    # Without W make Besag precision
    Q <- matrix(data = 0, nrow = n, ncol = n) # Empty matrix
    for(i in 1:n){
      if(nb[[i]][1] == 0){
        Q[i, i] = 0 
      }
      else {
        Q[i, i] <- length(nb[[i]]) 
      }
      Q[i, nb[[i]]] <- -1
    }
  } 
  else {
    # With W
    if(!identical(nb, mat_to_nb(W))) {
      stop("W should have a neighbourhood list corresponding to nb")
    }
    Q <- diag(rowSums(W)) - W
  }

  if(scale_precision) {
    Q <- scale_gmrf_precision(Q)$Q
  }
  
  if(return_sparse) {
   Q <- as(Q, "sparseMatrix")
  }
  
  return(Q)
}
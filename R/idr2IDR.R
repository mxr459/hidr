#' Calculates IDR value based on idr value
#'
#' The function Calculates global irreproducibility rate (IDR) given local irreproducibility rate.
#'
#' @param expr1 A n*p matrix or data frame of normalized gene expression data. The rows correspond samples (n), the columns correspond genes (p). 
#' 
#' @param expr2 The second gene expression data, should be in the same format, size as expr1.
#' 
#' @param s The tuning parameter for matrixes differences estimation, leave it as NULL to automatically select. 
#' 
#' @param s.seq The candidates for s selection.
#' 
#' @param alpha Prespecified level of false discovery rate. A relatively loose criterion is suggested for determines screening matrix.
#' 
#' @param verbose Set verbose to TURE to show details of s selection.
#' 
#' @details Please refer \bold{Yin et.al (2016). Testing differential networks with applications to the detection of gene-gene interactions. Biometrika(2015),pp. 1-20}
#' @export


idr2IDR <- function(idrv){
  o <- order(idrv)
  idrv.o <- idrv[o]
  idrv.rank <- rank(idrv.o, ties.method="max")
  
  top.mean <- function(index, x){mean(x[1:index])}
  
  # IDR of selected peaks
  IDRv.o <- sapply(idrv.rank, top.mean, idrv.o)
  IDRv <- rep(NA, length(IDRv.o))
  IDRv[o] <- IDRv.o
  return(IDRv)
}





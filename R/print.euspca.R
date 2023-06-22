print.euspca<-function(x, print.loadings=FALSE, ...){
  k<-x$k
  cat(paste(k,"uncorrelated sparse PCs",sep=" "), "\n")
  cat("% of explained var. :", format(round(x$p.ev, 2)), "\n")
  cat("% of non-zero loadings :", format(round(x$p.nz, 2)), "\n")
  if(print.loadings){
    rownames(x$loadings) = paste0("PC",1:k)
    cat("\n")
    cat("Sparse loadings \n")
    print(t(round(x$loadings,3)))
  }
  rownames(x$pc.cor) = paste0("PC",1:k)
  colnames(x$pc.cor) = paste0("PC",1:k)
  cat("\n")
  cat("Correlation of PCs \n")
  print(round(x$pc.cor,3))
  max.cor = max(abs( x$pc.cor - diag(k) ))
  cat("Max. abs. cor. :", format(round(max.cor, 3)), "\n")
}

Table_list_glm <- function(x) {
  dev <- unlist(lapply(x, deviance)) 
  dev <- signif(dev, digits = 9)
  ux <- unique(dev)
  frq.dev <- tabulate(match(dev, ux))
  odr <- order(frq.dev, decreasing = T)
  ux <- ux[odr]
  frq.dev <- frq.dev[odr]
  results <- list()
  for(i in 1:length(ux)){
    results[[i]] <- list(frequency = frq.dev[i], model = x[[which(match(dev, ux)==i)[1]]])
  }
  results
}
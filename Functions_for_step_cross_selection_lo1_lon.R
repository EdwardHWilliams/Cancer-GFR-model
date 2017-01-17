RMSE_step <- function (object, scope, scale = 0, direction = c("both", "backward", 
                                                               "forward"), 
                       trace = 1, keep = NULL, steps = 1000, frac_improv = 0.995, 
                       seed = seed, CV="lo1", ...) 
{
  mydeviance <- function(x, ...) {
    if(CV=="lo1"){
      leave_out_1_cross_val(x)
    } else if(CV=="lon"){
      leave_out_n_cross_val(x, seed=seed)
    }
    
  }
  cut.string <- function(string) {
    if (length(string) > 1L) 
      string[-1L] <- paste0("\n", string[-1L])
    string
  }
  re.arrange <- function(keep) {
    namr <- names(k1 <- keep[[1L]])
    namc <- names(keep)
    nc <- length(keep)
    nr <- length(k1)
    array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr, 
                                                           namc))
  }
  step.results <- function(models, fit, object, usingCp = FALSE) {
    change <- sapply(models, "[[", "change")
    rd <- sapply(models, "[[", "deviance")
    dd <- c(NA, abs(diff(rd)))
    rdf <- sapply(models, "[[", "df.resid")
    ddf <- c(NA, diff(rdf))
    RMSE <- sapply(models, "[[", "RMSE")
    heading <- c("Stepwise Model Path \nAnalysis of Deviance Table", 
                 "\nInitial Model:", deparse(formula(object)), "\nFinal Model:", 
                 deparse(formula(fit)), "\n")
    aod <- data.frame(Step = I(change), Df = ddf, Deviance = dd, 
                      `Resid. Df` = rdf, `Resid. Dev` = rd, RMSE = RMSE, 
                      check.names = FALSE)
    if (usingCp) {
      cn <- colnames(aod)
      cn[cn == "RMSE"] <- "Cp"
      colnames(aod) <- cn
    }
    attr(aod, "heading") <- heading
    fit$anova <- aod
    fit
  }
  Terms <- terms(object)
  object$call$formula <- object$formula <- Terms
  md <- missing(direction)
  direction <- match.arg(direction)
  backward <- direction == "both" | direction == "backward"
  forward <- direction == "both" | direction == "forward"
  if (missing(scope)) {
    fdrop <- numeric()
    fadd <- attr(Terms, "factors")
    if (md) 
      forward <- FALSE
  }
  else {
    if (is.list(scope)) {
      fdrop <- if (!is.null(fdrop <- scope$lower)) 
        attr(terms(update.formula(object, fdrop)), "factors")
      else numeric()
      fadd <- if (!is.null(fadd <- scope$upper)) 
        attr(terms(update.formula(object, fadd)), "factors")
    }
    else {
      fadd <- if (!is.null(fadd <- scope)) 
        attr(terms(update.formula(object, scope)), "factors")
      fdrop <- numeric()
    }
  }
  models <- vector("list", steps)
  if (!is.null(keep)) 
    keep.list <- vector("list", steps)
  n <- nobs(object, use.fallback = TRUE)
  fit <- object
  if(CV=="lo1"){
    bRMSE <- leave_out_1_cross_val(fit)
  } else if(CV=="lon"){
    bRMSE <- leave_out_n_cross_val(fit, seed=seed)
  }
  edf <- fit$df.null- fit$df.residual
  if (is.na(bRMSE)) 
    stop("RMSE is not defined for this model, so 'step' cannot proceed")
  if (bRMSE == -Inf) 
    stop("RMSE is -infinity for this model, so 'step' cannot proceed")
  nm <- 1
  if (trace) {
    cat("Start:  RMSE=", format(round(bRMSE, 2)), "\n", cut.string(deparse(formula(fit))), 
        "\n\n", sep = "")
    utils::flush.console()
  }
  models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n - 
                         edf, change = "", RMSE = bRMSE)
  if (!is.null(keep)) 
    keep.list[[nm]] <- keep(fit, bRMSE)
  usingCp <- FALSE
  while (steps > 0) {
    steps <- steps - 1
    RMSE <- bRMSE
    ffac <- attr(Terms, "factors")
    scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
    aod <- NULL
    change <- NULL
    if (backward && length(scope$drop)) {
      aod <- drop1_RMSE(fit, scope$drop, scale = scale, trace = trace, CV = CV, seed=seed)
      rn <- row.names(aod)
      row.names(aod) <- c(rn[1L], paste("-", rn[-1L], sep = " "))
      if (any(aod$Df == 0, na.rm = TRUE)) {
        zdf <- aod$Df == 0 & !is.na(aod$Df)
        change <- rev(rownames(aod)[zdf])[1L]
      }
    }
    if (is.null(change)) {
      if (forward && length(scope$add)) {
        aodf <- add1_RMSE(fit, scope$add, scale = scale, trace = trace, CV=CV, seed=seed)
        rn <- row.names(aodf)
        row.names(aodf) <- c(rn[1L], paste("+", rn[-1L], 
                                           sep = " "))
        aod <- if (is.null(aod)) 
          aodf
        else rbind(aod, aodf[-1, , drop = FALSE])
      }
      attr(aod, "heading") <- NULL
      nzdf <- if (!is.null(aod$Df)) 
        aod$Df != 0 | is.na(aod$Df)
      aod <- aod[nzdf, ]
      if (is.null(aod) || ncol(aod) == 0) 
        break
      nc <- match(c("Cp", "RMSE"), names(aod))
      nc <- nc[!is.na(nc)][1L]
      o <- order(aod[, nc])
      if (trace) 
        print(aod[o, ])
      if (o[1L] == 1) 
        break
      change <- rownames(aod)[o[1L]]
    }
    usingCp <- match("Cp", names(aod), 0L) > 0L
    fit_new <- update(fit, paste("~ .", change), evaluate = FALSE)
    fit_new <- eval.parent(fit_new)
    nnew <- nobs(fit, use.fallback = TRUE)
    if (all(is.finite(c(n, nnew))) && nnew != n) 
      stop("number of rows in use has changed: remove missing values?")
    Terms <- terms(fit_new)
    if(CV=="lo1"){
      bRMSE <- leave_out_1_cross_val(fit_new)
    } else if(CV=="lon"){
      bRMSE <- leave_out_n_cross_val(fit_new, seed=seed)
    }
    edf <- fit_new$df.null- fit_new$df.residual
    if (bRMSE >= RMSE * frac_improv) 
      break
    fit <- fit_new
    if (trace) {
      cat("\nStep:  RMSE=", format(round(bRMSE, 2)), "\n", 
          cut.string(deparse(formula(fit))), "\n\n", sep = "")
      utils::flush.console()
    }
    nm <- nm + 1
    models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n - 
                           edf, change = change, RMSE = bRMSE)
    if (!is.null(keep)) 
      keep.list[[nm]] <- keep(fit, bRMSE)
  }
  if (!is.null(keep)) 
    fit$keep <- re.arrange(keep.list[seq(nm)])
  step.results(models = models[seq(nm)], fit, object, usingCp)
}






add1_RMSE <- function (object, scope, scale = 0, x = NULL, seed = seed, CV=CV,  ...) 
{
  if (!is.character(scope)) 
    scope <- add.scope(object, update.formula(object, scope))
  if (!length(scope)) 
    stop("no terms in scope for adding to object")
  oTerms <- attr(object$terms, "term.labels")
  int <- attr(object$terms, "intercept")
  ns <- length(scope)
  dfs <- dev <- score <- numeric(ns + 1)
  names(dfs) <- names(dev) <- names(score) <- c("<none>", scope)
  add.rhs <- paste(scope, collapse = "+")
  add.rhs <- eval(parse(text = paste("~ . +", add.rhs), keep.source = FALSE))
  new.form <- update.formula(object, add.rhs)
  Terms <- terms(new.form)
  y <- object$y
  if (is.null(x)) {
    fc <- object$call
    fc$formula <- Terms
    fob <- list(call = fc, terms = Terms)
    class(fob) <- oldClass(object)
    m <- model.frame(fob, xlev = object$xlevels)
    offset <- model.offset(m)
    wt <- model.weights(m)
    x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    oldn <- length(y)
    y <- model.response(m)
    if (!is.factor(y)) 
      storage.mode(y) <- "double"
    if (NCOL(y) == 2) {
      n <- y[, 1] + y[, 2]
      y <- ifelse(n == 0, 0, y[, 1]/n)
      if (is.null(wt)) 
        wt <- rep.int(1, length(y))
      wt <- wt * n
    }
    newn <- length(y)
    if (newn < oldn) 
      warning(sprintf(ngettext(newn, "using the %d/%d row from a combined fit", 
                               "using the %d/%d rows from a combined fit"), 
                      newn, oldn), domain = NA)
  }
  else {
    wt <- object$prior.weights
    offset <- object$offset
  }
  n <- nrow(x)
  if (is.null(wt)) 
    wt <- rep.int(1, n)
  Terms <- attr(Terms, "term.labels")
  asgn <- attr(x, "assign")
  ousex <- match(asgn, match(oTerms, Terms), 0L) > 0L
  if (int) 
    ousex[1L] <- TRUE
  X <- x[, ousex, drop = FALSE]
  z <- glm.fit(X, y, wt, offset = offset, family = object$family, 
               control = object$control)
  dfs[1L] <- z$rank
  if(CV=="lo1"){
    dev[1L] <- leave_out_1_cross_val_glm.fit(object, X, y)
  } else if(CV=="lon"){
    dev[1L] <- leave_out_n_cross_val_glm.fit(object, X, y, seed=seed)
  }
  r <- z$residuals
  w <- z$weights
  sTerms <- sapply(strsplit(Terms, ":", fixed = TRUE), function(x) paste(sort(x), 
                                                                         collapse = ":"))
  for (tt in scope) {
    stt <- paste(sort(strsplit(tt, ":")[[1L]]), collapse = ":")
    usex <- match(asgn, match(stt, sTerms), 0L) > 0L
    X <- x[, usex | ousex, drop = FALSE]
    z <- glm.fit(X, y, wt, offset = offset, family = object$family, 
                 control = object$control)
    if(CV=="lo1"){
      dev[tt] <- leave_out_1_cross_val_glm.fit(object, X, y)
    } else if(CV=="lon"){
      dev[tt] <- leave_out_n_cross_val_glm.fit(object, X, y, seed=seed)
    }
    dfs[tt] <- z$rank
  }
  RMSE <- dev
  dfs <- dfs - dfs[1L]
  dfs[1L] <- NA
  aod <- data.frame(Df = dfs, RMSE = RMSE, row.names = names(dfs), 
                    check.names = FALSE)
  if (all(is.na(RMSE))) 
    aod <- aod[, -3]
  head <- c("Single term additions", "\nModel:", deparse(formula(object)), 
            if (scale > 0) paste("\nscale: ", format(scale), "\n"))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}








# Still not working with backwards selection


drop1_RMSE <- function (object, scope, scale = 0, seed = seed, CV=CV, ...) 
{
  
  x <- model.matrix(object)
  n <- nrow(x)
  asgn <- attr(x, "assign")
  tl <- attr(object$terms, "term.labels")
  if (missing(scope)) 
    scope <- drop.scope(object)
  else {
    if (!is.character(scope)) 
      scope <- attr(terms(update.formula(object, scope)), 
                    "term.labels")
    if (!all(match(scope, tl, 0L) > 0L)) 
      stop("scope is not a subset of term labels")
  }
  ndrop <- match(scope, tl)
  ns <- length(scope)
  dfs <- numeric(ns)
  dev <- numeric(ns)
  score <- numeric(ns)
  y <- object$y
  if (is.null(y)) {
    y <- model.response(model.frame(object))
    if (!is.factor(y)) 
      storage.mode(y) <- "double"
  }
  wt <- object$prior.weights
  if (is.null(wt)) 
    wt <- rep.int(1, n)
  for (i in seq_len(ns)) {
    ii <- seq_along(asgn)[asgn == ndrop[i]]
    jj <- setdiff(seq(ncol(x)), ii)
    z <- glm.fit(x[, jj, drop = FALSE], y, wt, offset = object$offset, 
                 family = object$family, control = object$control)
    X = x[, jj, drop = FALSE]
    dfs[i] <- z$rank
    if(CV=="lo1"){
      dev[i] <- leave_out_1_cross_val_glm.fit(object, X, y)
    } else if(CV=="lon"){
      dev[i] <- leave_out_n_cross_val_glm.fit(object, X, y, seed=seed)
    }
  }
  scope <- c("<none>", scope)
  dfs <- c(object$rank, dfs)
  if(CV=="lo1"){
    devc <- leave_out_1_cross_val(object)
  } else if(CV=="lon"){
    devc <- leave_out_n_cross_val(object, seed=seed)
  }
  dev <- c(devc, dev)
  RMSE <- dev
  dfs <- dfs[1L] - dfs
  dfs[1L] <- NA
  aod <- data.frame(Df = dfs,  RMSE = RMSE, row.names = scope, 
                    check.names = FALSE)
  if (all(is.na(RMSE))) 
    aod <- aod[, -3]
  
  head <- c("Single term deletions", "\nModel:", deparse(formula(object)), 
            if (!is.null(scale) && scale > 0) paste("\nscale: ", 
                                                    format(scale), "\n"))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}






leave_out_1_cross_val <- function(object){
  data <- object$data
  n <- dim(data)[1]
  res <- rep(0,n)
  for(i in 1:dim(data)[1]){
    fit <- glm(formula(object), data[-i,], family = family(object))
    res[i] <- (predict(fit, newdata=data[i,], type="response")-
                 data[i,][,which(colnames(data)==colnames(full$model)[1])])^2
  }
  sqrt(mean(res))
}
leave_out_1_cross_val_glm.fit <- function(object, X, y){
  data <- object$data
  n <- dim(data)[1]
  res <- rep(0,n)
  for(i in 1:n){
    fit <- glm.fit(X[-i,], y[-i],  family = family(object))
    X_new <- X[i,]
    prct <- (X_new %*% cbind(fit$coefficients))
    res[i] <- (prct-data[i,][,which(colnames(data)==colnames(full$model)[1])])^2
  }
  sqrt(mean(res))
}



leave_out_n_cross_val <- function(object, seed, n=5 ){
  set.seed(seed)
  data <- object$data
  x <- sample(1:dim(data)[1], dim(data)[1])
  split <- chunk(x,n)
  res <- rep(0,n)
  for(i in 1:n){
    fit <- glm(formula(object), data[-split[[i]],], family = family(object))
    res[i] <- rmse(predict(fit, newdata=data[split[[i]],], type="response"), 
                   data[split[[i]],][
                     ,which(colnames(data)==colnames(full$model)[1])])
  }
  mean(res)
}
leave_out_n_cross_val_glm.fit <- function(object, X, y, seed, n=5){
  set.seed(seed)
  data <- object$data
  x <- sample(1:dim(data)[1], dim(data)[1])
  split <- chunk(x,n)
  res <- rep(0,n)
  for(i in 1:n){
    fit <- glm.fit(X[-split[[i]],], y[-split[[i]]],  family = family(object))
    X_new <- X[split[[i]],]
    prct <- (X_new %*% cbind(fit$coefficients))
    res[i] <- rmse(prct, data[split[[i]],][
      ,which(colnames(data)==colnames(full$model)[1])])
  }
  mean(res)
}







# x1 = rnorm(100); x2 = rnorm(100); x3 = rnorm(100); x4 = rnorm(100); x5 = rnorm(100); x6 = rnorm(100)
# y = 5*x1+3*x3+2*x4+2*x1*x4+rnorm(100)
# df <- data.frame(y, x1,x2,x3,x4,x5,x6)
# null = glm(y~1, data=df, x=T, y=T)
# full = glm(y~.*., data=df, x=T, y=T)
# RMSE_step(null, scope = list(upper=full, lower=null), seed=sample(1:10000,1), trace=T, 
#           frac_improv = .999, CV = "lon")


# Select Build, Build and reload to build and lode into the R-session.

mylm <- function(formula, data = list(), contrasts = NULL, ...){
  # Extract model matrix & responses
  mf <- model.frame(formula = formula, data = data)
  X  <- model.matrix(attr(mf, "terms"), data = mf, contrasts.arg = contrasts)
  y  <- model.response(mf)
  terms <- attr(mf, "terms")


  # Add code here to calculate coefficients, residuals, fitted values, etc...
  # and store the results in the list est
  df <- nrow(X) - ncol(X)
  xtx_inv <- solve(t(X)%*% X)
  beta_hat <- xtx_inv %*% t(X) %*%y
  fitted_values <- X %*% beta_hat
  res <- drop(y- fitted_values)

  sse <- sum(res^2)
  sst <- sum((y-mean(y))^2)
  sigma2 <- sse / df
  vcov_mat <- as.numeric(sigma2) * xtx_inv
  r_squared <- 1 - (sse / sst)

  est <- list(terms = terms, model = mf)
  est$coefficients <- drop(beta_hat)
  est$fitted.values <- fitted_values
  est$residuals <- res
  est$df.residual <- df
  est$sse <- as.numeric(sse)
  est$sst <- as.numeric(sst)
  est$vcov <- vcov_mat
  est$r.squared <- as.numeric(r_squared)

  # Store call and formula used
  est$call <- match.call()
  est$formula <- formula

  # Set class name. This is very important!
  class(est) <- 'mylm'

  # Return the object with all results
  return(est)
}

print.mylm <- function(object, ...){
  # Code here is used when print(object) is used on objects of class "mylm"
  # Useful functions include cat, print.default and format
  cat('\nCall:\n')
  print(object$call)
  cat('\nCoefficients:\n')
  print(object$coefficients)
}

summary.mylm <- function(object, ...){
  # Code here is used when summary(object) is used on objects of class "mylm"
  # Useful functions include cat, print.default and format
  cat('Summary of object\n')
  se <- sqrt(diag(object$vcov))
  z_val <- object$coefficients / se
  p_val <- 2 * pnorm(-abs(z_val))
  coef_table <- cbind(Estimate = object$coefficients, `Std. Error` = se, `z value` = z_val, `Pr(>|z|)` = p_val)
  cat('\nCall:\n')
  print(object$call)
  cat('\nCoefficients:\n')
  printCoefmat(coef_table, P.values = TRUE, has.Pvalue = TRUE)
  cat(sprintf('\nMultiple R-squared:  %f\n', object$r.squared))
}

plot.mylm <- function(object, ...){
  # Code here is used when plot(object) is used on objects of class "mylm"
  library(ggplot2)
  labs(title = "Residuals vs Fitted Values", x = "Fitted values", y = "Raw residuals")
  library(ggplot2)
  df_plot <- data.frame(fitted = object$fitted.values, residuals = object$residuals)

  p <- ggplot(df_plot, aes(x = fitted, y = residuals)) + geom_point(alpha = 0.5) + geom_hline(yintercept = 0, col = "red", linetype = "dashed") + theme_minimal() +
    labs(title = "Residuals vs Fitted Values", x = "Fitted values", y = "Raw residuals")
  print(p)
}



# This part is optional! You do not have to implement anova
anova.mylm <- function(object, ...){
  # Code here is used when anova(object) is used on objects of class "mylm"

  # Components to test
  comp <- attr(object$terms, "term.labels")

  # Name of response
  response <- deparse(object$terms[[2]])

  # Fit the sequence of models
  txtFormula <- paste(response, "~", sep = "")
  model <- list()
  for(numComp in 1:length(comp)){
    if(numComp == 1){
      txtFormula <- paste(txtFormula, comp[numComp])
    }
    else{
      txtFormula <- paste(txtFormula, comp[numComp], sep = "+")
    }
    formula <- formula(txtFormula)
    model[[numComp]] <- mylm(formula = formula, data = object$model)
  }

  # Print Analysis of Variance Table
  cat('Analysis of Variance Table\n')
  cat(c('Response: ', response, '\n'), sep = '')
  cat('          Df  Sum sq X2 value Pr(>X2)\n')
  base_formula <- formula(paste(response, "~ 1"))
  base_model <- mylm(formula = base_formula, data = object$model)
  prev_sse <- base_model$sse
  prev_df <- base_model$df.residual
  for(numComp in 1:length(comp)){
    # Add code to print the line for each model tested
    curr_sse <- model[[numComp]]$sse
    curr_df <- model[[numComp]]$df.residual

    df_diff <- prev_df - curr_df
    ss_diff <- prev_sse - curr_sse

    chisq_val <- ss_diff / (object$sse / object$df.residual)
    p_val <- pchisq(chisq_val, df = df_diff, lower.tail = FALSE)

    cat(sprintf("%-9s %-3d %-7.2f %-8.2f %-8.4e\n",
                comp[numComp], df_diff, ss_diff, chisq_val, p_val))

    prev_sse <- curr_sse
    prev_df <- curr_df
  }

  return(model)

}

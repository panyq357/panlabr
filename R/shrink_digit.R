#' Shrink decimal columns in data.frame.
#'
#' @export
shrink_digit <- function(df, digits) {
  for (name in names(df)) {
    if (is.numeric(df[[name]]) && !is.integer(df[[name]])) {
      df[[name]] <- round(df[[name]], digits)
    }
  }
  return(df)
}

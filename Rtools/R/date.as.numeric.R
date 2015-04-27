date.as.numeric <- function(string = NULL, date.format = NULL, format = "%m/%d/%y")
{
  if(is.null(date.format))
  {
    if(is.null(string))
      stop("date or string must be specified")

    date.format <- as.Date(string, format = format)
  }

  date.numeric <- as.numeric(date.format) / 365.25 + 1970
  
  return(date.numeric)
}

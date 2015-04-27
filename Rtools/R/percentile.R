#Randy Johnson
#R code for computing percentiles
#Created March 18, 2004
#Last Modified June 18, 2004

percentile = function(x, distn)
{
    distn = subset(distn, !is.na(distn))
    n = length(distn)
    
    if(n<2 | is.null(distn) | max(distn)==min(distn))
    { 
        warning("Distribution vector is too short, NULL, or has only one value in it.")
        return(NA)
    }
    if(is.null(x))
    { 
        warning("Vector of quantiles is NULL.")
        return(NA)
    }

    y = findInterval(x, distn[order(distn)]) / n
    
    return(y)
}

# update.data.frame.R
# Randall Johnson
# George's update funtionality from SAS
# Created December 8, 2009
# Last Modified December 8, 2009

# object: the data.frame to be updated
# formula: a two sided formula, left hand side contains variables identifying each unique
#          row of data, right hand side contains variables to update
# newdata: the data.frame containing updated/new data

update.data.frame <- function(object, formula, ..., newdata)
{
  ######### Data #########
  
  # create variables we'll need from formula
  id <- strsplit(gsub(' + ', '#', as.character(formula)[2], fixed = TRUE), '#')[[1]]

  newvars <- strsplit(gsub(' + ', '#', as.character(formula)[3], fixed = TRUE), '#')[[1]]

  ######### Error Checking #########
  
  # check that the appropriate variables are present
  if(!all(id %in% names(object) & id %in% names(newdata)))
  {
    missing.vars <- id[!id %in% names(object) | !id %in% names(newdata)]
    stop(paste('Missing the following id variables:',
               paste(missing.vars, collapse = ', ')))
  }

  if(!all(newvars %in% names(newdata)))
  {
    missing.vars <- newvars[!newvars %in% names(newdata)]
    stop(paste('Missing the following replacement variables in newdata:',
               paste(missing.vars, collapse = ', ')))
  }

  # check that all id variable combinations of newdata are unique
  newdata.ids <- subset(newdata, select = id)
  newdata$newdata.ids <- apply(newdata.ids, 1, paste, collapse = '/')

  newdata.ids <- table(newdata$newdata.ids)
  if(any(newdata.ids > 1))
  {
    dups <- paste(names(newdata.ids)[newdata.ids > 1], collapse = ', ')
    if(nchar(dups) > 40)
      dups <- paste(substr(dups, 1, 40), '...', sep = '')
    
    stop(paste('The following ids / id combinations are duplicated in newdata:', dups))
  }
  
  # check that all id variable combinations of object are unique
  object.ids <- subset(object, select = id)
  object$object.ids <- apply(object.ids, 1, paste, collapse = '/')
  
  object.ids <- table(object$object.ids[object$object.ids %in% newdata$newdata.ids])
  if(any(object.ids > 1))
  {
    dups <- paste(names(object.ids)[object.ids > 1], collapse = ', ')
    if(nchar(dups) > 40)
      dups <- paste(substr(dups, 1, 40), '...', sep = '')
    
    stop(paste('The following ids / id combinations are duplicated in object:', dups))
  }

  ######### Merge Data #########
  
  # get changing rows
  updated <- subset(object, object.ids %in% newdata$newdata.ids)

  # drop variables that will be changing...keep the rest
  for(i in newvars)
  {
    updated[[i]] <- NULL
  }

  # merge static data with chaning data
  updated <- merge(updated, newdata)

  # get unchanging rows
  object <- subset(object, !object.ids %in% newdata$newdata.ids)

  # merge changing rows with unchanging rows
  object <- merge(object, updated, all = TRUE)

  # drop extra variables and return
  object$object.ids <- NULL
  object$newdata.ids <- NULL

  return(object)
}

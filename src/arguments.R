parserArgs <- function(args){
  parser <- list()
  name <- NULL
  
  num.arg <- length(args)
  
  index <- 1
  index.arg <- 1
  while(index <= num.arg){
    if(substr(args[index], 1, 1) == "-"){
      if(substr(args[index], 2, 2) == "-"){
        name.tmp <- substring(args[index], 3)
        
        parser.tmp <- NULL
        while(TRUE){
          index <- index + 1
          if(index > num.arg){
            break
          }
          
          if(substr(args[index], 1, 1) == "-"){
            index <- index - 1
            break
          }
          parser.tmp <- c(parser.tmp, args[index])
        }
        
        if(is.null(parser.tmp)){
          warning(paste("Argument", name.tmp, "is not set!"))
        }
        else{
          name <- c(name, name.tmp)
          parser[[index.arg]] <- parser.tmp
          index.arg <- index.arg + 1
        }
      }
      else{
        name.tmp <- substring(args[index], 2)
        
        parser.tmp <- NULL
        while(TRUE){
          index <- index + 1
          if(index > num.arg){
            break
          }
          if(substr(args[index], 1, 1) == "-"){
            index <- index - 1
            break
          }
          parser.tmp <- c(parser.tmp, args[index])
        }
        
        if(is.null(parser.tmp)){
          warning(paste("Argument", name.tmp, "is not set!"))
        }
        else{
          name <- c(name, name.tmp)
          parser[[index.arg]] <- parser.tmp
          index.arg <- index.arg + 1
        }
      }
    }
    else{
      warning(paste("Cannot recogonize argument:", args[index]))
    }
    index <- index + 1
  }
  
  names(parser) <- name
  return(parser)
}



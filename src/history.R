make.history = function(command, args) {
  
  time = format(Sys.time(), "%FT%XZ%z")
  args = paste(args, collapse = " ")
  history = paste0(time, ": ", command, " ", args)
  return(history)
  
}

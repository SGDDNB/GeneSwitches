#' @title Download example input data
#'
#' @description This function checks if example input files are in the current directory and if not download them
#'
#' @return
#'
#' @export
#'
get_example_inputData <- function(){
  logexpdata =  "./logexpdata.RData"
  if(!file.exists(logexpdata)) {
    res <- tryCatch(download.file("http://files.ddnetbio.com/logexpdata.RData",
                                  destfile = logexpdata,
                                  method = "auto"),
                    error = function(e) 1)
    if(res == 1) { print("Error: cannot download logexpdata.RData.") }
  }

  cardiac_monocle2 =  "./cardiac_monocle2.RData"
  if(!file.exists(cardiac_monocle2)) {
    res <- tryCatch(download.file("http://files.ddnetbio.com/cardiac_monocle2.RData",
                                  destfile = cardiac_monocle2,
                                  method = "auto"),
                    error = function(e) 1)
    if(res == 1) { print("Error: cannot download cardiac_monocle2.RData") }
  }
}

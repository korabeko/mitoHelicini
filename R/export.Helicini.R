#' Export Helicini dataset into your working directory.
#'
#' Exports the complete mitochondrial Helicini dataset into your working directory.
#' @param readme Boolean. Default FALSE, TRUE prints the readme file accompanying the data.
#' @importFrom utils packageVersion
#' @export
#' @return Creates an directory './Helicini_data' and copies the full mitochondrial Helicini dataset with metadata into this directory.
#' @examples

# version 28.01.2022

export.Helicini<-function(readme=FALSE){
exportDIR <- paste("./Helicini_data_",packageVersion("mitoHelicini"),sep="")
if(!dir.exists(file.path(exportDIR))){

  dir.create(exportDIR)
  file.copy(list.files(system.file("extdata",package = "mitoHelicini",mustWork=TRUE)),exportDIR,overwrite=FALSE)
  if(length(list.files(exportDIR))>0){print("Data sucessfully exported.")}else{warning("Data were not copied, but I do not know why.")}
  }else{
  stop(paste("The directory '",exportDIR,"' already exists and it will not be overwritten. Data were not exported.",sep=""))
}

if(readme==TRUE){
  print(strwrap(readLines(system.file("extdata","README.txt", package = "mitoHelicini",mustWork=TRUE)),width = 150),quote=FALSE)
}

} #end of function

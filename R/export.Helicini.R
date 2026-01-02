#' Export Helicini dataset into your working directory.
#'
#' Exports the complete mitochondrial Helicini dataset into your working directory.
#' @importFrom utils packageVersion
#' @export
#' @return Creates an directory './Helicini_data' and copies the full mitochondrial Helicini dataset with metadata into this directory.
#' @examples

# version 28.01.2022

export.Helicini<-function(readme=FALSE){
exportDIR <- paste("./Helicini_data_",packageVersion("mitoHelicini"),sep="")
if(!dir.exists(file.path(exportDIR))){

  dir.create(exportDIR)
  files<-list.files(system.file("extdata",package = "mitoHelicini",mustWork=TRUE))
  files<-paste0(system.file("extdata",package = "mitoHelicini",mustWork=TRUE),"/",files)
  file.copy(files,exportDIR,overwrite=FALSE)
  if(length(list.files(exportDIR))>0){print("Data sucessfully exported.")}else{warning("Data were not copied, but I do not know why.")}
  }else{
  stop(paste("The directory '",exportDIR,"' already exists and it will not be overwritten. Data were not exported.",sep=""))
}

} #end of function

#' Extract metadata for sequences in alignment.
#'
#' Selects and exports table rows containing metadata of sequences contained in an alignment.
#' @param tab data.frame. Table of localities data and sequence identifiers.
#' @param column integer. The position of the column containing sequence identifiers within tab. Default 1.
#' @param alignment vector or data.frame. Source alignment in sequential fasta format, where the sequence headers begin with an 8-letter sequence identifier or a list of sequence identifiers.
#' @param save Boolean. If true, the selected table rows are exported to a file in .ods format.
#' @param filename character. File name for the exported table
#' @importFrom readODS write_ods
#' @export
#' @return
#' @examples

# version 28.06.2022

tab.selector<-function(tab,alignment,column=1,save=TRUE,filename){
if(is.data.frame(alignment)==T){alignment<-alignment[,1]}
alignment<-alignment[!which(alignment=="")]
if(substr(alignment[1],1,1)==">"){
  identif<-alignment[seq.int(1L, length(alignment), 2L)]
}else{
  identif<-alignment
}
#tab<-cbind(tab,tab[,dim(tab)[2]])
#tab[,dim(tab)[2]]<-as.character(tab[,dim(tab)[2]])
#for (i in 1:length(identif)){
#	tab[which(tab[,column]==substr(identif[i], 2, 9)),dim(tab)[2]]<-"keep"
#}
#tab<-tab[which(tab[,dim(tab)[2]]=="keep"),1:dim(tab)[2]-1]
identif<-substr(identif,2,9)
tab<-tab[tab[,column]%in%identif,]
if(save==TRUE){write_ods(tab,filename,append=F)}
else{return(tab)}
}

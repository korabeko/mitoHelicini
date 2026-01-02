#' Select sequences from an alignment
#'
#' Selects sequences from alignment using a vector of sequence identifiers.
#' @param individuals character. List of sequence identifiers to be selected
#' @param alignment character vector or data.frame. Source alignment in sequential fasta format.
#' @param save Boolean. Save the selected sequences in a file? Default TRUE.
#' @param filename character. Name of the file to be saved.
#' @export
#' @examples
#' \dontrun{
#'  tab<-read_ods("210917_tabulka_mol_kopie.ods",col_names = TRUE,col_types = NA,sheet="Helicini")
#'  individuals<-tab$ID
#' }

# version 14.4.2022

selector<-function (individuals,alignment,save=T,filename){
if(class(alignment)=="data.frame"){alignment<-alignment[,1]}
vyber<-vector(length=length(individuals)*2,mode="character")
for (i in 1:length(individuals)){
	#print(paste(i,": ",individuals[i]))
	if(individuals[i]%in%substr(alignment,2,9)){ #11,18 could be set to take lab codes
	vyber[2*i-1]<-alignment[which(substr(alignment,2,9)==individuals[i])]
	if(length(vyber[2*i-1])!=length(alignment[which(substr(alignment,2,9)==individuals[i])])){
	  warning(paste("There is a duplicity in the alignment:",alignment[which(substr(alignment,2,9)==individuals[i])]))
	}
	vyber[2*i]<-alignment[which(substr(alignment,2,9)==individuals[i])+1]
	}
	}
if(save==T){write.table(vyber,filename,col.names=F,row.names=F,quote=F)}
return(vyber)
}

#' Which sequences in the alignment still do not have GenBank accession numbers?
#'
#' For each indvidual in an alignment, the function checks if it has a GenBank accession number in the specified column of a table. It is assumed that the accession contain 8 characters, which is not valid for some old accessions but does not matter for Helicini as yet.
#' @param tab data.frame. Table with columns containing sequence (individual) identifiers and accession numbers for those sequences already submitted to GenBank.
#' @param alignment data.frame. Fasta sequence alignment loaded in as data.frame by read.table. When loading, apply stringsAsFactors=F. Sequence headers begin with an 8-letter identifier.
#' @param IDs integer. Column number in tab containing the sequence identifiers, default 1.
#' @param acc.nos integer. Column number in tab containing the GenBank accession numbers.
#' @param filestem character. File name stem for the renamed alignment, default "submit_these".
#' @export
#' @return Alignment in sequential fasta format as a character vector, and saves it in a file along with the respective rows of tab.
#' @examples
#' \dontrun{
#'  tab<-read.table("new_ids.csv",header=T,as.is=T,sep="\t")
#'  alignment<-read.table("new.fas",header=F,as.is=T,sep="\t")
#'  IDs<-c(2,3,4)#columns with identifiers in the table
#'  newID<-1 #column in the table containing the new IDs
#' }

# version 10.03.2023

which.submit<-function(alignment,tab,IDs=1,acc.nos,filestem="submit_these"){

  if(is.data.frame(alignment)==T){alignment<-alignment[,1,drop=TRUE]}
  identifiers<-alignment[seq(1,length(alignment)-1,by=2)]
  identifiers<-substr(identifiers,2,9)
  submit<-sapply(identifiers,FUN=function(x){
    nchar(tab[which(tab[,IDs]==x),acc.nos])
  })
  submit<-submit[is.na(submit)==TRUE|submit!=8]
  submit<-names(submit)

  #extract tab entries for the sequences without accession number in the tab
  selected.rows<-tab[tab[,IDs]%in%submit,]
  write_ods(selected.rows,paste0(filestem,".ods"),append=F)

  #save the corresponding alignment entries
  vyber<-selector(individuals=submit,alignment=alignment,save=FALSE)
  write.table(vyber,paste0(filestem,".fas"),col.names=F,row.names=F,quote=F)

}

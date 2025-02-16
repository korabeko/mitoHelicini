#' Give sequences new identifiers
#'
#' Assigns new sequence identifiers from a table using original 8-letter identifier from the fasta sequence header.
#' @param tab data.frame. Table with columns containing the original and new identifiers. The original identifiers are 8 characters long.
#' @param alignment data.frame. Fasta sequence alignment loaded in as data.frame by read.table. When loading, apply stringsAsFactors=F. Sequence headers begin with an 8-letter identifier.
#' @param IDs integer. Position(s) of columns in "tab" with original identifier(s), in the order they should appear in sequence header. Multiple collumns may be specified e.g. when some identifiers are isolate codes, and other GenBank accession numbers.
#' @param newID integer. Position of collumn in "tab" containing the new identifiers
#' @param filename character. File name for the renamed alignment, default "new_ids_output.fas".
#' @export
#' @return Alignment in sequential fasta format as a character vector, and optionaly saves it in a file.
#' @examples
#' \dontrun{
#'  tab<-read.table("new_ids.csv",header=T,as.is=T,sep="\t")
#'  alignment<-read.table("dnes_12S.fas",header=F,as.is=T,sep="\t")
#'  IDs<-c(2,3,4)#columns with identifiers in the table
#'  newID<-1 #column in the table containing the new IDs
#' }

# version 17.11.2022

assign.ids<-function(tab,alignment,IDs,newID,filename="new_ids_output.fas"){

IDs<-c(IDs) #makes sure that IDs is a vector even if it has just one value
#make sure that alignment is a vector, not a data.frame
if(is.data.frame(alignment)==T){alignment<-alignment[,1,drop=TRUE]}
#extract original identifiers from alignment
identifiers<-alignment[seq(1,length(alignment)-1,by=2)]
identifiers<-substr(identifiers,2,9)

#check that each original identifier is only once in the alignment
if(length(which(duplicated(identifiers,incomparables=NA)==TRUE))>0){
stop(paste("the following sequence identifiers in the alignment are duplicates:",identifiers[duplicated(identifiers)],". This is a fatal error, the execution will stop!",sep=" "))
}

for(i in IDs){
  if(length(which(duplicated(tab[,i],incomparables=c(NA,""))==TRUE))>0){
    print(paste("in table column number ",i," the element on the row ",which(duplicated(tab[,i],incomparables=NA)==TRUE)," ",tab[which(duplicated(tab[,i],incomparables=NA)==TRUE),i]," is a duplicate of an earlier row",sep=""))
    stop("this is a fatal error, the execution will stop!")
  }
}

c<-rep("",length(identifiers)*2)
e<-1

a<-1
repeat{ #master repeat over individuals
  individual<-tab[(a+1)/2,IDs]
  individual<-data.frame(individual) #creates a data.frame with all IDs of a give individual
  d<-individual%in%identifiers #is any of the IDs in the alignment?
  IDsfound<-which(d==TRUE)
  if(length(IDsfound)==0){ #if not, skip and print the individual
    print((a+1)/2)
    print(paste("individual *",paste0(individual,collapse=" "),"* not found in the alignment",sep="",collapse=""))
    a<-a+2
  }
  else{
    z<-1
    repeat{
      c[e]<-paste0(">",tab[(a+1)/2,newID],"_",substr(alignment[which(substr(alignment,2,9)==individual[,d])],2,nchar(alignment[which(substr(alignment,2,9)==individual[,d])])),collapse="")
      c[e+1]<-alignment[which(substr(alignment,2,9)==individual[,d])+1]
      e<-e+2
      a<-a+2
      z<-z+1
      if(z>length(IDsfound)){break}
    }
  }
  if(a/2>dim(tab)[1]){break}
} #end master repeat over individuals

'%notin%'<-Negate('%in%')
not_in_table<-identifiers%notin%substr(c[seq(1,length(c)-1,by=2)],11,18)
message("assign.ids: The following sequences were not processed, because they were not found in tab:")
message(paste(identifiers[not_in_table],collapse=", ",sep=""),collapse="",sep=", ")
write.table(c,filename,col.names=F,row.names=F,quote=F)

}


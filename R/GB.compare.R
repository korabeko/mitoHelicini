#' Compare your record for a sequence with that in GenBank
#'
#' Compares the stored version of a sequence with the version available from GenBank. If the two sequences are different, it stores the both versions in a fasta file and prints the corresponding GenBank accession number to another file. Useful when checking that your published sequences were submitted correctly...
#' @param acc.nos data.frame. Table including accession numbers and individual identifiers (both assumed to be exactly 8 characters); the identifiers are assumed to be in the first column
#' @param loci integer. Vector of column positions of acc.nos in which the accession numbers are stored.
#' @param alignments list, data.frame or character. (List of) fasta sequence alignments loaded in as data.frame by read.table. When loading, apply stringsAsFactors=F. The order of alignments must be the same as a the order of loci in acc.nos.
#' @param file.seqs character. File name for output containing the differing pairs of sequences.
#' @param file.numbers character. File name for output containing the accession numbers of sequences where your record does not match GenBank.
#' @export
#' @return Textual output on the screen and two plain text files.
#' @importFrom stringr str_remove_all
#' @importFrom utils write.table
#' @importFrom ape read.GenBank
#' @examples
#' \dontrun{
#'  loci<-which(names(tab)=="16S_acc_nos"|names(tab)=="COX1_acc_nos") #specifies columns listing acc. numbers for 16S and COX1
#'  seq1<-read.table("211227_Helicini_16S.fas",header=F,stringsAsFactors=F)  #example, always change file name to match desired version
#'  seq2<-read.table("211227_Helicini_COX1.fas",header=F,stringsAsFactors=F)  #example, always change file name to match desired version
#'  alignments<-list(seq1,seq2)
#' }

# version 14.4.2022

GB.compare<-function(acc.nos,loci,alignments,file.seqs,file.numbers){
#require(ape)

pairs<-character()
dif.accessions<-character()

if(is.list(alignments)==TRUE){
alignments<-lapply(alignments,FUN=function(X){
  if(is.data.frame(X)==TRUE){
  X<-X[,1,drop=TRUE]}else{X<-X}
})

}else{
  if(is.data.frame(alignments)==TRUE){
      alignments<-alignments[,1,drop=TRUE]
  }
  alignments<-list(alignments)
}

for (j in 1:length(loci)){ #iterates over loci, extracting data
IDs.in.align<-substr(alignments[[j]][seq(1,length(alignments[[j]]),by=2)],2,9)
accessions<-acc.nos[IDs.in.align%in%acc.nos[,1],loci[j]]
accessions<-accessions[!is.na(accessions)]
accessions<-accessions[accessions!=""]
tryCatch(gb<-ape::read.GenBank(access.nb=accessions,as.character = T),
  error=function(cond){
  message("If the error message from ape::read.GenBank (see below) reads 'cannot open the connection to 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?...' but the connection URL includes valid accession numbers, check your internet connection.")
  message(cond)
  }
)
gb<-lapply(gb,FUN=function(x){paste(x,collapse="")})
gb<-c(unlist(gb))
gb.order<-names(gb)
gb<-toupper(gb)

for (i in 1:length(gb.order)){ #iterates over individual, comparing data
my.id<-which(acc.nos[,loci[j]]==gb.order[i])[1]
my<-character()
my<-alignments[[j]][which(substr(alignments[[j]],2,9)==acc.nos[my.id,1])+1]
my<-toupper(my)
my<-stringr::str_remove_all(my, "[?-]")
if(length(my)>0){
  if((my==gb[i])==F){
  print(paste(c(acc.nos[my.id,1]," does not match ",acc.nos[my.id,loci[j]]," (locus in column ",loci[j],")"),collapse=""))
  if(nchar(my)>nchar(gb[i])){print(paste(c(acc.nos[my.id,1],": the present record is longer than the GenBank version "),collapse=""))}
  pairs<-c(pairs,paste(">",acc.nos[my.id,1],"_",gb.order[i],"_your_version",sep="",collapse=""),my,paste(">",acc.nos[my.id,1],"_",gb.order[i],"_GB",sep="",collapse=""),gb[i])
  dif.accessions<-c(dif.accessions,accessions[i])
  } #end inner if
}else{message(paste("While there is the accession number ",gb.order[i]," in the table, corresponding sequence is not present in the alignment. Did you input your alignments in the same order as the corresponding collumns in acc.nos? j=",j,", i=",i,collapse = "",sep=""))}
} #end iteration over accession numbers in a locus
} #end iteration over loci

utils::write.table(dif.accessions,file=file.numbers,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
utils::write.table(pairs,file=file.seqs,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
} #end of function

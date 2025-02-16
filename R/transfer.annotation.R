#' Transfer annotation from a reference sequence to sequence(s) aligned with it.
#'
#' Tranfers annotation from a reference sequence to a sequence aligned with the reference. Not yet implemented.
#' @param ref ref character or data.frame. Fasta sequence loaded in as data.frame by read.table or as a vector of length 2, containing the reference sequence to which annotation applies. The reference MUST be aligned with the alignment.
#' @param annotation data.frame. Table with annotation, must contain columns "start" and "end" (integer) for the features to be exported in individual alignments and "attributes" with names of those features conforming the GFF3 format (https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
#' @param your.seq data.frame or character vector. The sequence to be annotated in fasta format.
#' @importFrom utils read.table
#' @export
#' @return
#' @examples

# version 17.08.2022

transfer.annotation<-function(ref=NULL,annotation=NULL,your.seq){
message("This function is not finished, it currently accepts only one sequence to be annotated at a time.")

#prepare data
  if(is.null(annotation)==TRUE|is.null(ref)==TRUE){
    message("At least one of annotation and ref is not specified; Helix pomatia mitogenome data (individual HE000505) from the package 'mitoHelicini' will be used.")
    files<-list.files(system.file("extdata",package = "mitoHelicini"))
    annotation.file<-files[grep(pattern="_annotations.gff$",files)]
    annotation<-read.table(paste(system.file("extdata",package = "mitoHelicini"),"/",annotation.file,sep=""),header=F,as.is=TRUE,sep="\t")
    colnames(annotation)<-c("seqid","source","type","start","end","score","strand","phase","attributes")
    annotation<-annotation[which(substr(annotation$seqid,1,8)=="HE000505"&(annotation$type=="CDS"|annotation$type=="rRNA"|annotation$type=="tRNA")&annotation$source=="current"),]
    alignment.file<-alignment.file[grep(pattern="_Helicini_mitogenomes.fas$",files)]
    alignment.default<-read.table(paste(system.file("extdata",package = "mitoHelicini"),"/",alignment.file,sep=""),header=F,as.is=TRUE)
    ref<-alignment.default[which(substr(alignment.default[,1],2,9)=="HE000505"):(which(substr(alignment.default[,1],2,9)=="HE000505")+1),1]
  }else{
    colnames(annotation)<-c("seqid","source","type","start","end","score","strand","phase","attributes")
  }

if(is.data.frame(your.seq)==TRUE){your.seq<-your.seq[,1,drop=TRUE]}
if(is.data.frame(ref)==TRUE){ref<-ref[,1,drop=TRUE]}

##from here: calculate start and end of features of the reference sequence
table.genes<-sub(";.*", "", annotation$attributes)
table.genes<-sub(".*name=", "", table.genes)

#the following block takes annotation and recalculates feature starts and ends in the aligned reference sequence, to be used later to backwards calculate the start and end in the new sequence
ref<-ref[2]
nuc<-unlist(strsplit(ref,split=""))
gaps<-which(nuc=="-")
ref.pos<-c(1:(nchar(ref)-length(gaps)))
print("gap positions in reference:")
for (i in 1:length(gaps)){
  if(gaps[i]<=max(ref.pos)){
    add.gap<-which(ref.pos==gaps[i])
    ref.pos[add.gap:length(ref.pos)]<-sapply(ref.pos[add.gap:length(ref.pos)],FUN=function(x){x+1})
    print(add.gap)
  }else{
    add.gap<-which(ref.pos==gaps[i])
    print(gaps[i])
  }
}

annotation$al.start<-ref.pos[annotation$start]
annotation$al.end<-ref.pos[annotation$end]

annotation$al.start[1]<-1
annotation$al.end[length(annotation$al.end)]<-nchar(your.seq[2])

table.genes<-cbind(annotation,table.genes)

#now calculate the positions of all gaps in the aligned new sequence and use that to correct the start and ends.
new.annotation<-annotation
new.annotation$start<-new.annotation$al.start
new.annotation$end<-new.annotation$al.end
your.gaps<-which(unlist(strsplit(your.seq[2],split=""))=="-")

if(length(your.gaps)>0){

for(j in length(your.gaps):1){
  new.annotation$start[new.annotation$start>your.gaps[j]]<-new.annotation$start[new.annotation$start>your.gaps[j]]-1
  new.annotation$end[new.annotation$end>=your.gaps[j]]<-new.annotation$end[new.annotation$end>=your.gaps[j]]-1
}
} #end if

new.annotation<-new.annotation[,1:9]
new.annotation[,1]<-rep(substr(your.seq[1],2,nchar(your.seq[1])),times=length(new.annotation[,1]))

return(new.annotation)
} #end of function

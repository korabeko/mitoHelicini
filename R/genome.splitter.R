#' Splits a mitogenome alignment into blocks specified by annotation file.
#'
#' Splits alignment (for example whole genomes) according to annotation file in GFF3 format referring to one of the aligned sequences (reference sequence). Meant for example for export of a single gene alignment from an alignment of mitogenomes.
#' @param alignment data.frame. Fasta sequence alignment loaded in as data.frame by read.table. When loading, apply as.is=TRUE. If not specified, the function prints a warning and uses internal data instead, ignoring parameters ref and annotation.
#' @param ref data.frame/character/integer The reference sequence to which the annotation file refers. Can be supplied as fasta sequence loaded in as data.frame by read.table or provided as a character vector of length 2. Alternatively, it may be an integer indicating its position in the alignment file. The reference MUST be aligned with the alignment. Typically included in the alignment file. If not specified, the function prints a warning and uses internal data instead, ignoring parameters alignment and annotation.
#' @param annotation data.frame. Table with annotation of the reference sequence, must contain columns "start" and "end" (integer) for the features to be exported in individual alignments and "attributes" with names of those features conforming the GFF3 format (https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
#' @param save Boolean. Save results in files?
#' @param file.stem character. File name stem for the exported alignments. Default "testfile".
#' @param overhang Boolean. Tells whether positions at the beginning and end of the alignment, unaligned to reference, are to be included in the exported alignments. Default TRUE
#' @export
#' @return A list of alignments in fasta format. If either annotation or ref are not specified, the use reference sequence and annotation file from the package 'mitoHelicini'.
#' @examples
#' \dontrun{
#'  #example:
#'  table<-read.table("HE000505_mitogenome_annotations.gff",header=F,sep="\t",stringsAsFactors=F) #read annotations file, which contains the positions of the single gene alignments relative to a reference Helix pomatia mitogenome.
#'  colnames(table)<-c("seqid","source","type","start","end","score","strand","phase","attributes")
#'  names<-colnames(table)
#'  names<-c(names,"al.start","al.end")
#'  table<-cbind(table,numeric(length(table[,1])),numeric(length(table[,1])))
#'  colnames(table)<-names
#'  annotation<-table[which(table$type=="experimental_feature"&table$source=="current"),]
#'  alignment<-read.table("211227_Helicini_mitogenomes.fas",header=F,stringsAsFactors=F)
#'  ref<-alignment[which(substr(alignment[,1],2,9)=="HE000505"):(which(substr(alignment[,1],2,9)=="HE000505")+1),1]
#'  file.stem<-"test"
#'  #end of example
#' }

# version 09.10.2022

genome.splitter<-function(alignment=NULL,ref=NULL,annotation=NULL,save=F,file.stem="testfile",overhang=T){

  message("genome.splitter: Not yet adapted to handle vectors, alignment must be data.frame.")

#if annotation file or reference sequence are not specified, the function will use internal data.
if(is.null(annotation)==TRUE|is.null(ref)==TRUE){
  message("At least one of annotation and ref is not specified in genome.splitter; Helicini mitogenome data from the package 'mitoHelicini' will be used!")
  
  files<-list.files(system.file("extdata",package = "mitoHelicini"))
  annotation.file<-files[grep(pattern="_annotations.gff$",files)]
  
 #annotation<-read.table(system.file("extdata","220930_mitogenome_annotations.gff", package = "mitoHelicini",mustWork=TRUE),header=F,as.is=TRUE,sep="\t")
  annotation<-read.table(paste(system.file("extdata",package = "mitoHelicini"),"/",annotation.file,sep=""),header=F,as.is=TRUE,sep="\t")
  colnames(annotation)<-c("seqid","source","type","start","end","score","strand","phase","attributes")
  annotation<-annotation[which(substr(annotation$seqid,1,8)=="HE000505"&annotation$type=="experimental_feature"&annotation$source=="current"),]

  alignment.file<-files[grep(pattern="_Helicini_mitogenomes.fas$",files)]
  alignment.default<-read.table(paste(system.file("extdata",package = "mitoHelicini"),"/",alignment.file,sep=""),header=F,as.is=TRUE,sep="\t")
  alignment.default<-alignment.default[,1,drop=TRUE]
  ref<-alignment.default[which(substr(alignment.default,2,9)=="HE000505"):(which(substr(alignment.default,2,9)=="HE000505")+1)]

  if(is.null(alignment)==TRUE){
    alignment<-alignment.default
  }
}  #else{

#if(is.null(alignment)==TRUE){
#  warning("genome.splitter: Alignment not specified, Helicini mitogenome data from the package 'mitoHelicini' will be used but may not correspond to user-provided reference sequence and annotation table")
#  alignment.file<-list.files(system.file("extdata",package = "mitoHelicini"))
#  alignment.file<-alignment.file[grep(pattern="_Helicini_mitogenomes.fas$",alignment.file)]
#  alignment.default<-read.table(paste(system.file("extdata",package = "mitoHelicini"),"/",alignment.file,sep=""),header=F,as.is=TRUE)
#  alignment<-alignment.default
#  }
#}
  #end loading default data if none specified

#ensure alignment and ref are vectors
if(is.data.frame(alignment)==T){alignment<-alignment[,1,drop=TRUE]}
if(is.data.frame(ref)==T){ref<-ref[,1,drop=TRUE]}
  
#when ref is specified as an x-th sequence in the alignment
  #there is no check that the number for ref is interger...
if(is.numeric(ref)==TRUE){ref.seq<-alignment[(ref*2-1):(ref*2)]}else{ref.seq<-ref}


#prepare names for the single gene alignments

colnames(annotation)<-c("seqid","source","type","start","end","score","strand","phase","attributes")
table.genes<-sub(";.*", "", annotation$attributes)
table.genes<-sub(".*name=", "", table.genes)

#the following block takes annotation and recalculates feature starts and ends in the aligned reference sequence
ref.seq<-ref.seq[2]
nuc<-unlist(strsplit(ref.seq,split=""))
gaps<-which(nuc=="-")
ref.pos<-c(1:(nchar(ref.seq)-length(gaps)))
print("gaps were introduced in the reference sequence by alignment with others at the following positions:")
if(length(gaps)>0){
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
}

annotation$al.start<-ref.pos[annotation$start]
annotation$al.end<-ref.pos[annotation$end]

if(overhang==T){
annotation$al.start[1]<-1
annotation$al.end[length(annotation$al.end)]<-nchar(alignment[2])
}

table.genes<-cbind(annotation,table.genes)

#produce single gene alignments
single.alignments<-vector(mode="list",length(table.genes$table.genes))
for (k in 1:length(table.genes$table.genes)){
g.alignment<-alignment
g.alignment[seq(2,length(g.alignment),by=2)]<-substr(g.alignment[seq(2,length(g.alignment),by=2)],table.genes$al.start[k],table.genes$al.end[k])
if(save==TRUE){
write.table(g.alignment,file=paste(file.stem,"_",table.genes$table.genes[k],".fas",collapse="",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
}
single.alignments[[k]]<-g.alignment
names(single.alignments)[k]<-table.genes$table.genes[k]
}

#name the single alignments within the list:

return(single.alignments)
} #end of function


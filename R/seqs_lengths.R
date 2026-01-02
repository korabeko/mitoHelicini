#' Calculate and plot the distribution of sequence lengths in an alignment
#'
#' Extracts the lengths of sequences in an alignment, returns a named vector of the lengths and plots their distribution
#' @param alignment vector or data.frame. Alignment in sequential fasta format.
#' @param ambiguous Boolean. Specifies if ambiguous bases should be included in the sequence length, default FALSE.
#' @param verbose Boolean. If TRUE, length statistics are printed to the screen and the length distribution plotted; default TRUE.
#' @importFrom stats median
#' @export
#' @return Produces a named vector of sequence lengths.

# version 10.1.2024

seqs_lengths<-function(alignment,ambiguous=FALSE,verbose=TRUE){
  if(is.data.frame(alignment)==T){alignment<-alignment[,1,drop=TRUE]}
  seq.names<-substr(alignment[seq(1,length(alignment),by=2)],2,9)
  alignment<-alignment[seq(2,length(alignment),by=2)]
  if(ambiguous==TRUE){
    b<-sapply(alignment,function(x){if(x==""){0}else{ #gregexpr returns -1 for x==""
      length(unlist(gregexpr(pattern ='A|C|G|T|W|S|M|K|R|Y|B|D|H|V|N',x)))}
      },USE.NAMES = F)
  }else{
    b<-sapply(alignment,function(x){if(x==""){0}else{
      length(unlist(gregexpr(pattern ='A|C|G|T',x)))}
      },USE.NAMES = F)
  }
  if(length(b)!=length(seq.names)){stop("There is a problem with the alignment supplied to seqs.lengths, number of FASTA headers does not match the number of sequences (odd number of rows).")}
  names(b)<-seq.names
  if(verbose==TRUE){
    plot(b[order(b)],ylab="sequence length",xlab="sequence rank")
    print(paste("median lenght:",median(b)))
    print(paste("mean lenght:",mean(b)))
    print(paste("min. lenght:",min(b)))
    print(paste("max. lenght:",max(b)))
  }
  return(b)
}

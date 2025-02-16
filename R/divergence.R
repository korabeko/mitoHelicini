#' XXXXX
#'
#' Xxxxx
#' @param alignment vector or data.frame. Alignment in sequential fasta format as read in by read.table.
#' @param slindow integer. Size of a sliding window in bp, default 101.
#' @param circular Boolean. Is the sequence circular? Default TRUE.
#' @param type character. Nucleotide ("nt") or amino acid ("aa") sequence.
#' @details The current implementations uses step 1 and allows for up to 25% missing data in the window. "-" is interpreted as gap and considered a valid character state. Code missing data as "?"! Ambuguity codes interpreted as missing data.
#' @export
#' @return

# version 14.4.2022 #cvicny.fas se muze pouzit

divergence<-function(alignment,slindow=101,circular=TRUE,type=c("nt","aa")){
  
  #write a check for type!!!!!
  message("function not tested yet")
  
  nucleotides<-c("A","C","T","G","a","c","t","g","-")
  aminoacids<-c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","-")
  
  #transform alignment from data.frame to vector, if necessary
  if(is.data.frame(alignment)==T){alignment<-alignment[,1,drop=TRUE]}
  #extract sequence names
  seq.names<-substr(alignment[seq(1,length(alignment),by=2)],2,9)
  #extract sequences
  alignment<-alignment[seq(2,length(alignment),by=2)]
  lengths<-vapply(alignment,FUN=nchar,FUN.VALUE = numeric(length=1),USE.NAMES=FALSE)
  if(sum(lengths)!=length(lengths)*lengths[1]){stop("Sequences are of different lengths!")}
  
  #calculate the number of pairwise comparisons
  size<-(length(seq.names)*(length(seq.names)-1))/2
  if(size>1000)message("The number of pairwise comparisons is high, consider splitting your job.") #should better be based on the product of sequence length and the number of sequences.
  
  #create list for results
  mate1<-character(length = size)
  mate2<-character(length = size)
  divergences<-matrix(NA,nrow=size,ncol=nchar(alignment[1]))
  
  #main cycle, running all pairwise comparisons
  for(i in 1:(length(alignment)-1)){
    for(j in (i+1):length(alignment)){
      
      mate1[((i-1)*length(alignment))-((i-1)*((i-1)+1)/2)+(j-i)]<-seq.names[i]
      mate2[((i-1)*length(alignment))-((i-1)*((i-1)+1)/2)+(j-i)]<-seq.names[j]
      print(paste(i, seq.names[i], j, seq.names[j]))
      comparison<-rbind(unlist(strsplit(alignment[i],split="")),unlist(strsplit(alignment[j],split="")))
      
      if(type=="nt"){
        d<-apply(comparison,MARGIN=2,FUN=function(x){
          if(x[1]%in%nucleotides&x[2]%in%nucleotides){
          x[1]==x[2]
          }else{
            NA
          }
        })
      }
      
      if(type=="aa"){
        d<-apply(comparison,MARGIN=2,FUN=function(x){
          if(x[1]%in%aminoacids&x[2]%in%aminoacids){
            x[1]==x[2]
          }else{
            NA
          }
        })
      }
      
      # calculate mean in sliding window                   
      for(k in 1:(nchar(alignment[1])-slindow)){
        positions<-c(k:(k+slindow-1))
        dpercent<-1-mean(d[positions],na.rm=TRUE)
        if(sum(is.na(d[positions]))>slindow/4){dpercent<-NA}
        divergences[((i-1)*length(alignment))-((i-1)*((i-1)+1)/2)+(j-i),positions[slindow%/%2]]<-dpercent
      }
      
      if(circular==TRUE){
        for(k in c((nchar(alignment[1])-slindow+1):nchar(alignment[1]))){
          positions<-c(k:nchar(alignment[1]),1:((nchar(alignment[1])-(k+slindow-1))*(-1)))
          dpercent<-1-mean(d[positions],na.rm=TRUE)
          if(sum(is.na(d[positions]))>slindow/4){dpercent<-NA}
          divergences[((i-1)*length(alignment))-((i-1)*((i-1)+1)/2)+(j-i),positions[slindow%/%2]]<-dpercent
        }
      } #end if circular
      
    } #end inner for
  } #end outer for
  
  
  results<-cbind.data.frame(mate1,mate2,divergences)
  print(class(results[,3]))
  return(results)
  
} #function end


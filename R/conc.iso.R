#' Concatenate sequences based on identifiers in the fasta header.
#'
#' Performs concatenation of sequences using 8-letter identifiers at the beginning of fasta headers.
#' @param alignments list. List containing fasta sequence alignments loaded in as data.frame or character vector by read.table.
#' @param save Boolean. If TRUE, the results will be saved into a fasta file.
#' @param file character. Name of the file containing the concatenated alignment, if save=TRUE)
#' @param uknown character. The character used to denote missing data when the individual is not represented by a given locus. Does affect coding of missing data in the concatenated alignments. Default "?", other meaningful options are "-" (same as gap) or "N" (it is assumed that there is a nucleotide).
#' @export
#' @return Returns concatenated alignment in a sequential FASTA format as an character vector.

# version 17.04.2022

conc.iso<-function (alignments,save=T,file="concatenate.fas",uknown="?"){

#starts by taking the first alignment file and recording which individuals from the other files are not present in the first one.
differences<-vector(mode="character")
uniqs<-vector(mode="character")

#standardize input
alignments<-lapply(alignments,FUN=function(x){
  if(is.data.frame(x)){
    x[,1,drop=TRUE]
  }else{
    x
  }
})

#scans all alignments for sample identifiers
names<-unlist(lapply(alignments,FUN=function(x){
  if(length(x)>=2){
    a<-x[seq(from=1,to=length(x),by=2)]
    a<-a[a!=""]
    a<-unlist(sapply(a,FUN=(function(z){
      substr(a,1,9)
      }))) #end inner lapply
  } #end if
}))
names<-unique(names)


#scans all alignments 2 to n for individuals not present in the first one.
#for(k in 2:length(alignments)){
# uniqs<-setdiff(substr(alignments[[k]][2*(1:(length(alignments[[k]])/2))-1],1,9),substr(alignments[[1]][2*(1:(length(alignments[[1]])/2))-1],1,9))
# differences<-c(differences,uniqs)
#}
#unique2<-unique(differences)

#appends the sequence IDs from the first alignment to those from the previous step to have a complete list of all individuals represented in alignments
#names<-substr(alignments[[1]][seq(1,length(alignments[[1]]),by=2)],1,9)
#names<-c(names,unique2)

#creates an empty concatenated aligment
n.seq<-length(names)
concatenated<-rep("",n.seq*2)


for(j in 1:length(alignments)){ #iterates over single gene alignments
  if(length(alignments[[j]])>=2){

 #print(paste("alignment:",j)) #for testing

 gaps1<-paste(rep(uknown,times=nchar(alignments[[j]][2])),collapse="") #creates gaps for sequences missing in alignment j

 for (i in 1:n.seq){

 #print(paste("seq:",i)) #for testing

 acc.seq1<-names[i]

 if (acc.seq1%in%substr(alignments[[j]],1,9)){
  	concatenated[2*i]<-paste(concatenated[2*i],alignments[[j]][which(substr(alignments[[j]],1,9)==acc.seq1)+1],sep="")
  	}else{
    	concatenated[2*i]<-paste(concatenated[2*i],gaps1,sep="",collapse="")
  	} #end if else

 if (acc.seq1%in%substr(alignments[[j]],1,9)){
  if(nchar(concatenated[(2*i)-1])==0){
    	concatenated[(2*i)-1]<-paste(concatenated[(2*i)-1],alignments[[j]][which(substr(alignments[[j]],1,9)==acc.seq1)],sep="")
	  }
  } #end if

 } #end inner loop
}else{
  print(paste("conc.iso: alignment ",names(alignments)[[j]] ," was empty and is not represented in the final concatenated alignment"))
} #end outer if
} #end outer loop

if(save==T){write.table(concatenated,file,col.names=F,row.names=F,quote=F)}
return(concatenated)
#end of function
}

#' Try to locate mitochondrial origins of replication following https://doi.org/10.1016/j.mito.2014.05.009
#'
#' Uses the distribution of GC skew along the heavy strand to locate the origins of replication.
#' @param alignment character. Alignment in sequential FASTA format. No empty rows allowed. May be also data.frame with one column.
#' @param alignment integer. Size of the sliding window, recommended 150-350
#' @details The function assumes that a H strand alignment is provided. The GC skew is calculated in non-overlapping 251 bp windows (250 used originally in the cited paper) for all the sequences in the window together. Dssh is calculated for different combinations of possible OH a OL locations (locations of the centers of the windows are used for OL, OH and p). For each combination of OH and OL, linear regression model of the observed skew vs. the calculated Dssh for each position p is fitted and the combination with the lowest sum of squared errors (SSEs) is chosen as the result. Significance could be assessed by randomization, but we did not need it.
#' @export
#' @return Numeric vector containing the approximate locations of OH and OL (in this order).

# version 14.1.2024
#alignment<-read.table("mitogenomes_cvicny.fas") # data used in development

replication_origin<-function(alignment,window){

# preparing the data
if(is.data.frame(alignment)==TRUE){alignment<-alignment[,1,drop=TRUE]}

# in the cited paper, they used alignments; here we tried it also sequence by sequence and examined the distributions of the inferred OH and OL positions
#alignment.source<-alignment
#alignment.nogaps<-alignment
#alignment.nogaps[seq(from=2,to=length(alignment),by=2)]<-gsub("-",replacement="",alignment[seq(from=2,to=length(alignment),by=2)])
#for (j in seq(from=1,to=length(alignment.source),by=2)){
#alignment<-alignment.nogaps[j:(j+1)]

# create windows (they used 250, but an odd number has a clearly defined middle nucleotide)
windows<-seq(from=(window%/%2)+1,to=nchar(alignment[2]),by=window)
nwindows<-length(windows)

# iterate over windows, calculate GC skew
# in Deuterostomes, most genes are coded on the L strand, so that is the strand in sequence databases and the one they calculated the GC skew for
# in land snails, most genes are coded on the H strand, so our alignments are also the H strand
# we thus calculate the complementary GC skew, to calculate it for the L strand like in the cited paper
skews<-numeric()
for (i in 1:nwindows){
  #extract sequences within window, count Cs and Gs
  strings<-substr(alignment[seq(from=2,to=length(alignment),by=2)],windows[i]-(window%/%2),windows[i]+(window%/%2))
  Gs<-lengths(regmatches(strings, gregexpr("G", strings)))
  Cs<-lengths(regmatches(strings, gregexpr("C", strings)))
  # GC skew = (G − C)/(G + C) for the L strand, we work with the H strand
  skew<-(sum(Cs)-sum(Gs))/(sum(Gs)+sum(Cs))
  skews[i]<-skew
} # end iterating over windows


#specify the formula to calculate Dssh
formula<-function(OL,OH, p, L, windows){ 
  # the mitogenome is circular... this situation needs a fake coordinate system with respect to OH to calculate the distances
  # OH, OL, and p can only take values specified in windows  
    offset<-OH-min(windows) # offset is by how much is OH shifted with respect to the starting OH
    OH<-min(windows) # we put OH at position 1 in windows
    if(OL<OH){OL<-OL-offset+max(windows)}else{OL<-OL-offset}
    if(p<=OH){p<-p-offset+max(windows)}else{p<-p-offset}
    
    if((p-OH)<=(OL-OH)){Dssh<-((2*(OL-p))/L)} #checked
    if((p-OH)>(OL-OH)&(OL-OH)<L/2){Dssh<-((2*(p-OH))/L)} #result must be positive, so OH-p would not work
    if((p-OH)>(OL-OH)&(OL-OH)>=L/2){Dssh<-((L-(2*(p-OL)))/L)} #checked
  
  return(Dssh)
}

# prepare variables
L<-nchar(alignment[2])
results<-matrix(nrow=(nwindows*nwindows*(nwindows-1)),ncol=5)
combination<-character(length=(nwindows*nwindows*(nwindows-1)))

#calculate the Dssh for each OH and OL combination and all genome positions
a<-1
for(p in windows){
  for(OH in windows){
    for(OL in windows[which(windows!=OH)]){
      Dssh<-formula(p=p,OH=OH,OL=OL,L=L,windows=windows)
      results[a,1]<-p
      results[a,2]<-skews[which(windows==(p))]
      results[a,3]<-OH
      results[a,4]<-OL
      results[a,5]<-Dssh
      combination[a]<-paste0("OH",OH,"OL",OL)
      a<-a+1
    }
  }
}

# put results in one data.frame
results<-cbind.data.frame(results,combination,stringsAsFactors = TRUE)
colnames(results)<-c("p","skew","OH","OL","Dssh","combination")

# for each OH and OL combination, get the SSE and slope
correlations<-numeric(length=length(levels(results$combination)))
slopes<-numeric(length=length(levels(results$combination)))
b<-1
for(i in levels(results$combination)){
  s<-results$skew[which(results$combination==i)]
  D<-results$Dssh[which(results$combination==i)]
  model<-lm(s~D)
  correlations[b]<-sum((fitted(model)-s)^2) # calculate SSE
  slopes[b]<-model$coefficients[2]
  b<-b+1
}

correlations<-cbind.data.frame(levels(results$combination),correlations,slopes) # put results into one datasets, using just cbind would make it all character as the first vector
correlations<-correlations[correlations[,3]<0,]# remove negative correlations


location<-correlations[which(correlations[,2]==min(correlations[,2])),1] # saves the locations of OH and OL
origins<-as.numeric(c(substr(location,start=regexpr("OH",location)+2,stop=regexpr("OL",location)-1),substr(location,start=regexpr("OL",location)+2,stop=nchar(location))))
# in the paper, they analyzed the H strand but needed to return coordinates on the L strand, so I assume that I must switch the values in the end to get correct coordinates

#summarize results
names(origins)<-c("OL","OH")
print(alignment[1])
print(paste0("OH: ",origins[2]," OL: ",origins[1]))
return(c(origins[2],origins[1]))
# } #end of iteration over sequences, see lines just below data loading

} # end of function

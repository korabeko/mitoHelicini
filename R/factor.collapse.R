#' Collapse identical sequences by factor, keeping the most complete one.
#'
#' Identifies identical sequences for each factor level (locality, species, ...) and keeps only the most complete one. The first is kept in case of a tie. Indels are considered missing data and not taken into account in comparisons.
#' @param alignment data.frame or character vector.
#' @param tab data.frame. Table including sequence identifiers (8-letter) and the factor used for alignment thinning. Default NULL, see parameter fact.
#' @param fact integer. Column in 'tab' containing the factor. Default NULL causes the function to remove all duplicated sequences.
#' @param id integer. Column in 'tab' containing sequence IDs. Default 1
#' @param prefix character. File name stem for the resulting files.
#' @details level has three levels, corresponding to Table S1 in Korábek et al. 2022: 'genus' for genus, 'clade' for a roughly species-level clades (these clades do not always correspond to species borders), and 'subclade' for intraspecific clades used in figs. S4-62.
#' @export
#' @return Writes the thinned alignment along with the corresponding thinned table into files.
#' @examples

# version 15.4.2022 UNFINISHED: it remains to be tested and debugged.

#setwd("~/helix-clanek/BEAST_phylogeo")
#alignment<-read.table("pomatia.fas",header=F, as.is=T)
#tab<-read.table("localities_pomatia_copy.csv",sep="\t",header=T, as.is=T)
#fact<-4
#id<-1
#prefix<-"pomatia"

factor.collapse<-function(alignment,tab=NULL,fact=NULL,id=1,prefix){

message("UNFINISHED: it remains to be tested and debugged.")

if(is.data.frame(alignment)==TRUE){alignment<-alignment[,1,drop=TRUE]}

#creates data frame from alignment where sequence headers are on the same row as sequence
alig2<-cbind(alignment[seq(from=1,to=length(alignment),by=2)],alignment[seq(from=2,to=length(alignment),by=2)])
alig2[,1]<-sapply(alig2[,1],FUN=substr,start=2,stop=9)

if(is.null(fact)==TRUE){
  tab<-cbind(substr(alignment[seq(from=1,to=length(alignment),by=2)],2,9),rep("a",times=length(alignment)/2))
  tab[,2]<-as.factor(tab[,2])
  fact<-2
}

#identifies sets of sequences for each factor level
duplicated.fact<-tab[duplicated(tab[,fact]),fact]
multi<-tab[tab[,fact] %in% duplicated.fact,]
"%not_in%"<-Negate("%in%")
single<-tab[tab[,fact] %not_in% duplicated.fact,]

facts<-levels(base::as.factor(multi[,fact]))

# for each factor level, the function takes the first sequence in table and compares it with the next. If they are identical, it discards the second or the shorter one, and updates the table. If they are different by at least one resolved site, it moves to the second and repeats the procedure.

for (i in 1:length(facts)){
site<-multi[which(multi[,fact]==facts[i]),]
alig3<-alig2[alig2[,1]%in%site[,id],]

a<-1
repeat{
j<-a+1
repeat{
len1<-length(unlist(gregexpr(pattern ='A|C|G|T',alig3[a,2])))
len2<-length(unlist(gregexpr(pattern ='A|C|G|T',alig3[j,2])))
overlap<-c(unlist(gregexpr(pattern ='A|C|G|T',alig3[a,2])),unlist(gregexpr(pattern ='A|C|G|T',alig3[j,2])))
overlap<-overlap[duplicated(overlap)]
print(paste("compared",alig3[a,1],alig3[j,1]))
print(paste("i =",i,"a =",a,"j =",j))
if (all(unlist(strsplit(alig3[a,2],""))[overlap]==unlist(strsplit(alig3[j,2],""))[overlap])){
	if(len1==len2){
	site<-site[which(site[,id]!=alig3[j,1]),]
	print(paste("removed",alig3[j,1]))
	}
	if(len1!=len2){
	shorter<-which(c(len1,len2)==min(len1,len2))
	site<-site[which(site[,id]!=alig3[c(a,j)[shorter],1]),]#i = 152 a = 3 j = 6"
	print(paste("removed",alig3[c(a,j)[shorter],1]))
	}
alig3<-alig2[alig2[,1]%in%site[,id],]
j<-j
}
else{j<-j+1}
if(j>=dim(site)[1]){break}
}
a<-a+1
if(a>=dim(site)[1]){break}
}
single<-rbind(single,site)
}



alig2<-alig2[alig2[,1]%in%single[,id],]
alig2[,1]<-sapply(alig2[,1],function(x){paste(">",x,sep="",collapse=NULL)})
fin.alig<-character(length(alig2))
b<-1
d<-2
repeat{
fin.alig[b]<-alig2[(b+1)/2,1]
fin.alig[d]<-alig2[d/2,2]
b<-b+2
d<-d+2
if(d>length(fin.alig)){break}
}

write.table(fin.alig, file=paste(prefix,"_reduced.fas",sep="",collapse=NULL), quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(single, file=paste("localities_",prefix,"_reduced.csv",sep="",collapse=NULL), quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
}




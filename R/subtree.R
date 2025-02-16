#' Extracting subtree from a phylogeny while keeping posterior supports from BEAST.
#'
#' Attaches trees to a backbone phylogeny by replacing subtrees with less complete sampling with larger trees containing the same clade, stores posterior probability read by read.annotated.nexus in a new slot and writes them to file as a Newick comment.
#' @param tree phylo. Rooted phylogeny (is.rooted test not yet implemented) stored as the phylo class of ape package.
#' @param clade character vector or single integer. The subtree to be extracted. May be defined by a vector containing names of taxa defining the subtree, or by a node number.
#' @param save Boolean. If true, file "supertree.tre" is saved in translated nexus format.
#' @importFrom ape getMRCA extract.clade node.depth.edgelength
#' @export
#' @return Returns a tree of the class phylo, and optionally saves it in a file.


#version 13.03.2023 ### resolve dependencies to ape - see check() output

#tree<-read.annotated.nexus("sub_albescens.tree")

subtree<-function(tree,clade,save=FALSE){

  #if it exists, extract posterior
  if(!is.null(tree$annotations)&length(tree$annotations)>0){
    if(!is.null(tree$annotations[[1]]$posterior)){
      treeposterior<-unlist(lapply(tree$annotations,FUN=function(x){if(is.null(x$posterior)==FALSE){x$posterior}else{1}}))
      tree$posterior<-treeposterior
    }
  }else{
    if(!is.null(tree$posterior)){
      treeposterior<-tree$posterior
    }
  }

  #finding the subclade to be extracted, first using node number
    #checks follow to see that it is correctly specified
  if(is.integer(clade)){
    if(length(clade)==1){
      if(clade>length(tree$tip.label)&clade<2*length(tree$tip.label)-1){
        detachpoint<-clade
      }else{
        stop("Please, specify an INTERNAL node.")
      }
    }else{
      stop("More than one subclade specified.")
    }
  }

    #checks that the taxa specified are included in the tree
  if(is.character(clade)){
    if(sum(clade%in%tree$tip.label)==length(clade)){
      detachpoint<-getMRCA(tree,clade)
    }else{
      stop("Some of the specified taxa were not found in the tree.")
    }
  }

  #extract the subtree with ape
  subtree<-extract.clade(tree,detachpoint,root.edge=0,collapse.singles=TRUE)
  #strip subtree of annotations that are now no longer valid
  subtree<-subtree[1:4] #beware, this also removes the S3 class attribute
  class(subtree)<-"phylo"

  if(exists("treeposterior")){
    subtree$posterior<-treeposterior[which(tree$edge[,2]==which(tree$tip.label==subtree$tip.label[1])):which(tree$edge[,2]==which(tree$tip.label==subtree$tip.label[length(subtree$tip.label)]))]
  }
  #write data
  if(save==TRUE){write_supertree(subtree)}
  return(subtree)
}

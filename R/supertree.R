#' Constructing supertree by replacing subtrees of a backbone phylogeny with new trees (presumably containing more samples).
#'
#' Attaches trees to a backbone phylogeny by replacing subtrees with less complete sampling with larger trees containing the same clade, stores posterior probability read by read.annotated.nexus in a new slot and writes them to file as a Newick comment.
#' @param backbone phylo. Rooted backbone phylogeny (is.rooted test not yet implemented) stored as the phylo class of ape package.
#' @param tree phylo. Rooted phylogeny (is.rooted test not yet implemented) to replace a subtree in backbone, stored as the phylo class of ape package.
#' @param save Boolean. If true, file "supertree.tre" is saved in translated nexus format.
#' @importFrom ape drop.tip getMRCA bind.tree extract.clade node.depth.edgelength
#' @export
#' @return Returns a tree of the class phylo, and optionally saves it in a file.


#version 16.02.2023 ### resolve dependencies to ape - see check() output


#backbone<-read.annotated.nexus("backbone.tre")
#tree<-read.annotated.nexus("tree.tre")

#backbone<-read.annotated.nexus("sub_backbone.tree")
#tree<-read.annotated.nexus("sub_albescens.tree")

supertree<-function(backbone,tree,save=FALSE){

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

  #if it exists, extract posterior
  if(!is.null(backbone$annotations)&length(backbone$annotations)>0){
    if(!is.null(backbone$annotations[[1]]$posterior)){
      backposterior<-unlist(lapply(backbone$annotations,FUN=function(x){if(is.null(x$posterior)==FALSE){x$posterior}else{1}}))
      backbone$posterior<-backposterior
    }
  }else{
    if(!is.null(backbone$posterior)){
      backposterior<-backbone$posterior
    }
  }

#write test that the overlap is at least two tips and that the mrca of those in the subtree is the root
  overlap<-tree$tip.label[tree$tip.label%in%backbone$tip.label]
  if(length(overlap)<2){stop("Tree cannot be attached, because it is represented by less than two branches in the backbone.")}
  #getMRCA(tree,tree$tip.label[tree$tip.label%in%backbone$tip.label])==length(tree$tip.label)+1

  #find node at which the subtree is attached in the original backbone
  detachpoint<-ape::getMRCA(backbone,tree$tip.label[tree$tip.label%in%backbone$tip.label])
  branchlength<-backbone$edge.length[which(backbone$edge[,2]==detachpoint)]
  removed<-ape::extract.clade(backbone,detachpoint)

  if(length(removed$tip.label)>length(tree$tip.label)){message("The replaced subtree has more tips than the replacement. Not a supported scenario, the function may fail or produce spurious results.")}
  #test that we are the replaced and replacement subtrees are concordant
  if(identical(removed$tip.label[order(removed$tip.label)],overlap[order(overlap)])==FALSE){stop("The replaced subtree contains some tips that are not present in the replacement tree.")}

  backbone2<-ape::drop.tip(backbone,removed$tip.label,trim.internal = FALSE)

  if(length(which(backbone2$tip.label=="NA"))>1){
    repeat{
      backbone2<-ape::drop.tip(backbone2,c("NA"),trim.internal = FALSE) #repeat until all are removed
      if(length(which(backbone2$tip.label=="NA"))<2){break}
    }
  }

  attachpoint<-which(backbone2$tip.label=="NA")

  #for some reason, ape screws the length of the branch from where subtree is removed, so we have saved it earlier from unaltered backbone tree and correct it here after dropping subtree
  backbone2$edge.length[which(backbone2$edge[,2]==attachpoint)]<-branchlength

  #adapt tree length to the backbone
  target.heigth<-abs(ape::node.depth.edgelength(backbone2)-max(ape::node.depth.edgelength(backbone2)))[attachpoint]
  subtree.root.height<-abs(ape::node.depth.edgelength(tree)-max(ape::node.depth.edgelength(tree)))[length(tree$tip.label)+1]
  attach.ratio<-target.heigth/subtree.root.height
  tree$edge.length<-tree$edge.length*attach.ratio

  #attach tree
  super<-ape::bind.tree(backbone2,tree,attachpoint)

#backposterior1<-backposterior[1:which(backbone$edge[,2]==detachpoint)]
#backposterior2<-backposterior[(which(backbone$edge[,2]==detachpoint)+1):length(backbone$edge.length)]

#c(backbone$edge.length[1:which(backbone$edge[,2]==detachpoint)],tree$edge.length,backbone$edge.length[(1+which(backbone$edge[,2]==detachpoint)+length(removed$edge.length)):length(backbone$edge.length)])
  if(exists("treeposterior")&exists("backposterior")){
    super$posterior<-c(backposterior[1:which(backbone$edge[,2]==detachpoint)],treeposterior,backposterior[(1+which(backbone$edge[,2]==detachpoint)+length(removed$edge.length)):length(backbone$edge.length)])
  }
  #write data
  if(save==TRUE){write_supertree(super)}
  return(super)
}

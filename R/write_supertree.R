#' Write tree with posterior probabilities (crudely adapted from ape)
#'
#' Writes tree in nexus format with posterior probabilities stored as comments.
#' @param phylo phylo. Phylogenetic tree
#' @param file character. Filename of the tree saved, default="supertree.tre".
#' @param translate Booloean. Save the tree in translated form or with taxon names? Default translate=TRUE

#' @return
#' @examples
#' \dontrun{
#'  seq1<-read.table("SSU-mafft-vyber-21-09-16.fas",header=F,stringsAsFactors=F)
#'  seq2<-read.table("COI-mafft-vyber-21-09-16.fas",header=F,stringsAsFactors=F)
#'  alignments<-list(seq1,seq2)
#' }

# version 14.02.2023



## write.tree.R (2019-03-01)

##   Write Tree File in Parenthetic Format

## Copyright 2002-2019 Emmanuel Paradis, Daniel Lawson, and Klaus Schliep

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

checkLabel <- function(x)
{
  ## delete all leading and trailing spaces and tabs, and
  ## the leading left and trailing right parentheses:
  ## (the syntax will work with any mix of these characters,
  ##  e.g., "    ( ( ((  " will correctly be deleted)
  x <- gsub("^[[:space:]\\(]+", "", x)
  x <- gsub("[[:space:]\\)]+$", "", x)
  ## replace all spaces and tabs by underscores:
  x <- gsub("[[:space:]]", "_", x)
  ## replace commas, colons, and semicolons with dashes:
  x <- gsub("[,:;]", "-", x)
  ## replace left and right parentheses with dashes:
  x <- gsub("[\\(\\)]", "-", x)
  x
}

write.tree<-function(phy,file="",append=FALSE,digits=10,tree.names=FALSE){

    if (!(inherits(phy, c("phylo", "multiPhylo"))) &&
        !all(vapply(phy, inherits, logical(1), 'phylo')))
      stop("object \"phy\" must contain trees")

    if (inherits(phy, "phylo")) phy <- c(phy)
    N <- length(phy)
    res <- character(N)

    if (is.logical(tree.names)) {
      if (tree.names) {
        tree.names <-
          if (is.null(names(phy))) character(N)
        else names(phy)
      } else tree.names <- character(N)
    }

    ## added by KS (2019-03-01):
    check_tips <- TRUE
    if (inherits(phy, "multiPhylo")) {
      if (!is.null(attr(phy, "TipLabel"))) {
        attr(phy, "TipLabel") <- checkLabel(attr(phy, "TipLabel"))
        check_tips <- FALSE
      }
    }

    ## added by EP (2019-01-23):
    phy <- .uncompressTipLabel(phy)
    class(phy) <- NULL

    for (i in 1:N)
      res[i] <- .write.tree2(phy[[i]], digits = digits,
                             tree.prefix = tree.names[i], check_tips)

    if (file == "") return(res)
    else cat(res, file = file, append = append, sep = "\n")
  }



.write.tree2 <- function(phy, digits = 10, tree.prefix = "", check_tips)
{
  brl <- !is.null(phy$edge.length)
  nodelab <- !is.null(phy$node.label)
  if (check_tips) phy$tip.label <- checkLabel(phy$tip.label)
  if (nodelab) phy$node.label <- checkLabel(phy$node.label)
  f.d <- paste("%.", digits, "g", sep = "")
  cp <- function(x){
    STRING[k] <<- x
    k <<- k + 1
  }
  add.internal <- function(i) {
    cp("(")
    desc <- kids[[i]]
    for (j in desc) {
      if (j > n) add.internal(j)
      else add.terminal(ind[j])
      if (j != desc[length(desc)]) cp(",")
    }
    cp(")")
    if (nodelab && i > n) cp(phy$node.label[i - n]) # fixed by Naim Matasci (2010-12-07)
    if (brl) {
      ############add node annotation if present
      cp(paste0("[&posterior=",phy$posterior[ind[i]],"]"))
      ############add
      cp(":")
      cp(sprintf(f.d, phy$edge.length[ind[i]]))
    }

  }
  add.terminal <- function(i) {
    cp(phy$tip.label[phy$edge[i, 2]])
    if (brl) {
      cp(":")
      cp(sprintf(f.d, phy$edge.length[i]))
    }
  }

  n <- length(phy$tip.label)

  ## borrowed from phangorn:
  parent <- phy$edge[, 1]
  children <- phy$edge[, 2]
  kids <- vector("list", n + phy$Nnode)
  for (i in 1:length(parent))
    kids[[parent[i]]] <- c(kids[[parent[i]]], children[i])

  ind <- match(1:max(phy$edge), phy$edge[, 2])

  LS <- 4*n + 5
  if (brl) LS <- LS + 4*n
  if (nodelab)  LS <- LS + n
  STRING <- character(LS)
  k <- 1
  cp(tree.prefix)
  cp("(")
  getRoot <- function(phy)
    phy$edge[, 1][!match(phy$edge[, 1], phy$edge[, 2], 0)][1]
  root <- getRoot(phy) # replaced n+1 with root - root has not be n+1
  desc <- kids[[root]]
  for (j in desc) {
    if (j > n) add.internal(j)
    else add.terminal(ind[j])
    if (j != desc[length(desc)]) cp(",")
  }

  if (is.null(phy$root.edge)) {
    cp(")")
    if (nodelab) cp(phy$node.label[1])
    cp(";")
  }
  else {
    cp(")")
    if (nodelab) cp(phy$node.label[1])
    cp(":")
    cp(sprintf(f.d, phy$root.edge))
    cp(";")
  }
  paste(STRING, collapse = "")
}













#function (..., file = "", translate = TRUE)
#{
#  obj <- .getTreesFromDotdotdot(...)
#  ntree <- length(obj)

write_supertree<-function(tree,file="supertree.tre",translate=TRUE){
    obj<-list(tree)
    ntree <- length(obj)

    cat("#NEXUS\n", file = file)
    cat(paste("[R-package APE, ", date(), "]\n\n", sep = ""),
        file = file, append = TRUE)

    N <- length(obj[[1]]$tip.label)

    cat("BEGIN TAXA;\n", file = file, append = TRUE)
    cat(paste("\tDIMENSIONS NTAX = ", N, ";\n", sep = ""),
        file = file, append = TRUE)
    cat("\tTAXLABELS\n", file = file, append = TRUE)
    cat(paste("\t\t", obj[[1]]$tip.label, sep = ""),
        sep = "\n", file = file, append = TRUE)
    cat("\t;\n", file = file, append = TRUE)
    cat("END;\n", file = file, append = TRUE)

    cat("BEGIN TREES;\n", file = file, append = TRUE)
    if (translate) {
      cat("\tTRANSLATE\n", file = file, append = TRUE)
      obj <- .compressTipLabel(obj)
      X <- paste("\t\t", 1:N, "\t", attr(obj, "TipLabel"), ",", sep = "")
      ## We remove the last comma:
      X[length(X)] <- gsub(",", "", X[length(X)])
      cat(X, file = file, append = TRUE, sep = "\n")
      cat("\t;\n", file = file, append = TRUE)
      class(obj) <- NULL
      for (i in 1:ntree)
        obj[[i]]$tip.label <- as.character(1:N)
    } else {
      if (is.null(attr(obj, "TipLabel"))) {
        for (i in 1:ntree)
          obj[[i]]$tip.label <- checkLabel(obj[[i]]$tip.label)
      } else {
        attr(obj, "TipLabel") <- checkLabel(attr(obj, "TipLabel"))
        obj <- .uncompressTipLabel(obj)
      }
    }

    title <- names(obj)
    if (is.null(title))
      title <- rep("UNTITLED", ntree)
    else {
      if (any(s <- title == "")) title[s] <- "UNTITLED"
    }

    for (i in 1:ntree) {
      if (class(obj[[i]]) != "phylo") next
      root.tag <- if (is.rooted(obj[[i]])) "= [&R] " else "= [&U] "
      cat("\tTREE *", title[i], root.tag, file = file, append = TRUE)
      cat(write.tree(obj[[i]], file = ""),
          "\n", sep = "", file = file, append = TRUE)
    }
    cat("END;\n", file = file, append = TRUE)
  }

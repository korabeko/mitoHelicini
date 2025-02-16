#' Export representative mitochondrial alignment for the Helicini tribe allowing for choosing divergence level, loci, and individuals.
#'
#' Finds the individual with most complete data for each clade at the specified divergence level, concatenates the specified loci and outputs the alignment. Useful when you want to show the phylogenetic position of your new sample relative to what is already known, but want to reduce the data. In output it provides some useful info about the data.
#' @param level character. Divergence level to be represented. See details for allowed values. Only a single level is accepted. Using ID would export all individuals (see also "individuals").
#' @param genes character. Loci to be included. See details for allowed values.
#' @param genus character. Genera to be considered. If not present in the data, they will be ignored without warning.
#' @param species character. Species to be considered. If not present in the data, they will be ignored without warning.
#' @param individuals character. Optional, vector listing individual IDs to be considered in the selection. Can also be a data.frame, but then it is assumed that the IDs are in the first column.
#' @param weight character. Either "max" or "average" (use "max" as default). Specifies the scoring mechanism used, see details.
#' @param maxloci Boolean. If TRUE, the individual with best coverage (see details) is selected as usual, but then additional individuals are added if these help to increase the number of loci available for each of the groups. Default FALSE.
#' @param save Boolean. Save the final alignment to file? Default TRUE.
#' @param file character. Name for the file containing the exported representative alignment.
#' @param save_single Boolean. Specifies whether to save single gene alignments ("best.alignments" in the function output. Default FALSE.).
#' @details You may choose from a list of alignment files to output only selected loci. By default, all are exported in the order as they appear in the mitogenomes ("COX1","TRNV","16S","TRNL1_TRNA","ND6","TRNP","ND5","ND1_ND4L","CYTB","TRND-TRNF","COX2","TRNY-TRNL2","ATP8","TRNN","ATP6","TRNR_TRNE","12S_TRNM","ND3","TRNS2_TRNT","COX3","TRNS1","ND4","TRNI","ND2","TRNK"). Running with "loci" at default value outputs complete mitogenomes where available (with gaps and question marks for missing data in part of or the whole locus, respectively).
#' @details "level" has three levels, corresponding to Table S1 in Korábek et al. 2022: 'genus' for genus, 'clade' for a roughly species-level clades (these clades do not always correspond to species borders), and 'subclade' for intraspecific clades used in figs. S4-62.
#' @details The individuals are selected as follows: Samples within each clade specified by "level" are ranked for the sequence length of each locus. Loci are then weighted by their alignment length or by average non-zero sequence length and the individual with the greatest weighted sequence length is chosen as representative of its clade.
#' @details The "maxloci" option, if on, outputs the individual with best coverage overall, but if some locus is not available for this individual, it looks for an individual that has it (the one with the longest sequence) and adds also this individual to the final alignment.
#' @importFrom readODS read_ods
#' @export
#' @return Returns a list of the following objects: "per.group.inds","lengths.ordered.groups","missing.loci","best.alignments","final.alignment". "per.group.inds" gives the number of individuals in each group specified by "level" parameter. 'missing.loci' has three values: 0=locus present in the best scoring individual, 1=locus missing in the best scoring individual ,2=locus missing for the group altogether. "lengths.ordered.groups" contains all the sequence lengths structured by group and locus.
#' @examples

# version 10.1.2024
# calls selector, seq_lengths, conc.iso

representative.alignment<-function(level=c("ID","subclade","clade","genus"),genes=c("COX1","TRNV","16S","TRNL1_TRNA","ND6","TRNP","ND5","ND1_ND4L","CYTB","TRND-TRNF","COX2","TRNY-TRNL2","ATP8","TRNN","ATP6","TRNR_TRNE","12S_TRNM","ND3","TRNS2_TRNT","COX3","TRNS1","ND4","TRNI","ND2","TRNK"),genus=NULL,species=NULL,individuals=NULL,weight=c("max","average"),maxloci=FALSE,save=TRUE,file="representative.fas",save_single=FALSE){
  #temporary warnings and errors - none
  #complete or nearly complete mitogenomes: a<-representative.alignment(level="ID",individuals=c("HE003273", "HE003274", "HE003275", "HE003276", "HE003277", "HE002644", "HE003132", "HE001286", "HE001246", "HE003254", "HE003253", "HE000505", "HE002946", "HE003257", "HE003148", "HE003251", "HE000808", "HE001725", "HE002405", "HE000991", "HE002487","HE003272","HE003268","HE002964","HE001273","HE002954","HE002865","HE002843","HE002600","HE001338","HE001230","HE002250","HE002944"),weight="max",maxloci=FALSE,file ="lucorum.fas")
  
  # check that arguments are OK
  level<-match.arg(level)
  if(length(weight)!=1){stop("Weight must be specified! See help for details")}
  if(maxloci==TRUE){weight<-match.arg(weight)}else{weight<-"max"}
  genes<-match.arg(genes,several.ok = TRUE)
  if(level=="ID"&maxloci==TRUE){stop("If exporting data for selected or all individuals (level=\"ID\"), maxloci must be set to FALSE")}

  #read the master table with sequence metadata
  tab.file<-list.files(system.file("extdata",package = "mitoHelicini"))
  tab.file<-tab.file[grep(pattern="_Helicini_table_mol.ods$",tab.file)]
  tab<-read_ods(path=paste(system.file("extdata",package = "mitoHelicini"),"/",tab.file[1],sep=""),sheet = 1)
  #remove any NA rows
  tab<-tab[!is.na(tab$ID),]

  #find alignments of the chosen loci
  alignment.file<-list.files(system.file("extdata",package = "mitoHelicini"))
  alignment.file<-unlist(sapply(genes, FUN=function(x){alignment.file[grep(pattern=paste("Helicini_",x,sep=""),alignment.file)]}))

  #read the alignments
  alignments<-lapply(alignment.file,FUN=function(x){
    read.table(paste(system.file("extdata",package = "mitoHelicini"),"/",x,sep=""),header=F,as.is=TRUE,sep="\t")[,1,drop=TRUE]
  })

  #extract the genera wanted, if specified
  if(!is.null(genus)){
    message("representative_alignment: Only representatives of specified genera will be considered.")
    if(is.data.frame(genus)==T){genus<-genus[,1,drop=TRUE]}
    tab<-tab[tab$genus%in%genus,]
    alignments<-lapply(alignments,FUN=function(x){
      selector(tab$ID,x,save=FALSE)
    })
  }
  
  #extract the genera wanted, if specified
  if(!is.null(species)){
    message("representative_alignment: Only representatives of specified species will be considered.")
    if(is.data.frame(species)==T){species<-species[,1,drop=TRUE]}
    tab<-tab[tab$species%in%species,]
    alignments<-lapply(alignments,FUN=function(x){
      selector(tab$ID,x,save=FALSE)
    })
  }
  
  #extract the individuals wanted, if specified
  if(!is.null(individuals)){
    message("representative_alignment: Only specified individuals will be considered.")
    if(is.data.frame(individuals)==T){individuals<-individuals[,1,drop=TRUE]}
    alignments<-lapply(alignments,FUN=function(x){
      selector(individuals,x,save=FALSE)
    })
    tab<-tab[tab$ID%in%individuals,]
  }

  #calculate the lengths of individuals sequences in the alignments (would be sequence completeness if divided by alignment length)
  lengths<-lapply(alignments,FUN=function(w){
    seqs_lengths(w,verbose=FALSE)
    })

  #order the sequence lengths from longest to shortest
  lengths.ordered<-lapply(lengths,FUN=function(v){v[order(v,decreasing=TRUE)]})

  #extract individuals for each clade specified by "level" into a named list
  column<-which(names(tab)==level)
  IDs.clades<-lapply(levels(as.factor(tab[,column,drop=TRUE])),FUN=function(u){
    tab$ID[which(tab[,column,drop=TRUE]==u)]
  })
  names(IDs.clades)<-lapply(levels(as.factor(tab[,column,drop=TRUE])),FUN=function(u){
    u
  })

  #calculates the number of individuals per clade
  per.group.inds<-as.data.frame(summary(IDs.clades)[,1],optional=TRUE)
  names(per.group.inds)<-c("n.ind")

  #for each group, extract the sorted lengths for each locus
  #creates two-level named list, with groups as the upper level and loci in the lower, nested level
  lengths.ordered.groups<-lapply(IDs.clades,FUN=function(s){
    lapply(lengths.ordered,FUN=function(t){
      t[names(t)%in%s]
    })
  })

  ####select the best individual for each group, which is the one with best average rank in each locus weighted by the overall number of individuals for that locus
  
  #here the selection procedure starts
  #calculate weights for each alignment file
  if(weight=="max"){
    gene.weights<-lapply(alignments,FUN=function(r){
      if(is.data.frame(r)==T){nchar(r[2,1])} #if and else added as a workaround, because rbind when merging alignments outputs vector instead of data.frame, do not know why...
      else{
        #nchar(r[2]) #failed if the first position in alignment was empty, which happened when alignments were filtered by individuals vector
        max(nchar(r[seq(2,length(r),by=2)]))
        }
      })
  }
    
  if(weight=="average"){
    gene.weights<-lapply(lengths,FUN=function(r){
      mean(r[r>0])
    })
  }

  gene.weights<-unlist(gene.weights)/max(unlist(gene.weights)) #gene weight in the per-individual score is determined by gene length in the full alignment, mean length of actually sequenced fragment per gene would possibly be better, or per-group maximum length

  if(maxloci==FALSE){
    #calculate weighted alignment length for each individual
    ind.lengths<-numeric(length(tab$ID))
    names(ind.lengths)<-tab$ID
    for (i in 1:length(lengths.ordered.groups)){ #outer loop iterates over groups
      for (j in 1:length(lengths.ordered.groups[[i]])){ #inner loop iterates over loci
        for(k in 1:length(lengths.ordered.groups[[i]][[j]])){ #iterates over individuals
          a<-which(names(ind.lengths)==names(lengths.ordered.groups[[i]][[j]][k]))
          ind.lengths[a]<-ind.lengths[a]+(gene.weights[j]*lengths.ordered.groups[[i]][[j]][k])
        }
      }
    }

  
  #create a list of weighted alignment length for each individual, grouped by groups
  ind.lengths.groups<-lapply(levels(as.factor(tab[,column,drop=TRUE])),FUN=function(u){
    ind.lengths[names(ind.lengths)%in%tab$ID[which(tab[,column,drop=TRUE]==u)]]
  })
  names(ind.lengths.groups)<-lapply(levels(as.factor(tab[,column,drop=TRUE])),FUN=function(u){
    u
  })
  
  } #end if maxloci=FALSE

  #################################
  if(maxloci==FALSE){  #find the best individual in each group, report which loci are missing from final alignment

    #find the best individual in each group
    best<-unlist(lapply(ind.lengths.groups,FUN=function(o){
      names(o)[which(o==max(o))[1]] #takes the first ==max to resolve ties
    }))

    #which gene is missing in which group? result saved in data.frame missing.loci
    missing.loci<-as.data.frame(matrix(data=0,nrow=length(IDs.clades),ncol=length(genes)))
    #to.add<-vector(mode = "list", (length(IDs.clades)*length(genes))) #stores group (order) and locus (order) that should be added
    colnames(missing.loci)<-genes
    rownames(missing.loci)<-names(IDs.clades)
    for (i in 1:length(lengths.ordered.groups)){ #outer loop iterates over groups
      for (j in 1:length(lengths.ordered.groups[[i]])){ #inner loop iterates over loci
        if(!any(names(lengths.ordered.groups[[i]][[j]])%in%best)){#condition asks if any of the individuals for which the locus exists are present in 'best' (i.e. among the individuals chosen in the first round)
         missing.loci[i,j]<-1
        }
       if(length(lengths.ordered.groups[[i]][[j]])==0){
         missing.loci[i,j]<-2
        }
      }
    }
  
  } #end if maxloci==FALSE

  #################################

  if(maxloci==TRUE){ #find the best individual in each group and locus, report which loci are missing from final alignment

    #which loci are represented by which individuals?
    loci.inds<-lapply(alignments,FUN=function(x){
      substr(x[seq(from=1,to=length(x),by=2)],2,9)
    })
    
    
    #which gene is missing in which group? result saved in data.frame missing.loci. differs by only considering loci missing altogether
    missing.loci<-as.data.frame(matrix(data=0,nrow=length(IDs.clades),ncol=length(genes)))
    colnames(missing.loci)<-genes
    rownames(missing.loci)<-names(IDs.clades)
    for (i in 1:length(lengths.ordered.groups)){ #outer loop iterates over groups
      for (j in 1:length(lengths.ordered.groups[[i]])){ #inner loop iterates over loci
        if(length(lengths.ordered.groups[[i]][[j]])==0){
          missing.loci[i,j]<-2
        }
      }
    }

    #take group then locus and extract the longest sequence(s): loci.best is a list per group per locus of the IDs of the longest sequence or sequences (if there is a tie)
    loci.best<-vector(mode="list",length=length(IDs.clades)) #produces list where best individuals for each group and locus are grouped first by locus, then by group
    names(loci.best)<-names(IDs.clades)
    for(i in 1:length(IDs.clades)){
      for(j in 1:length(alignments)){
        if(length(lengths.ordered.groups[[i]][[j]])>0){longest<-max(lengths.ordered.groups[[i]][[j]])}
        loci.best[[i]][[j]]<-names(lengths.ordered.groups[[i]][[j]])[which((longest-lengths.ordered.groups[[i]][[j]])<longest/100)] #find longest sequences for each locus, introduce tolerance to small length variations
      }
    names(loci.best[[i]])<-names(alignments)
    }
    #so now i have all the best individuals for each group and locus


    best<-vector(mode="list",length=length(IDs.clades))
    for(i in 1:length(loci.best)){ #iterates over groups
      xx<-table(unlist(loci.best[[i]])) #names in xx are the best IDs combined from across loci in a group, xx is the number of loci in these individuals
      if(length(xx)>0){best.in.group<-names(xx[which(xx==max(xx))])[1]}else{best.in.group<-character()} #best in group is, among the longest sequences as identified for each locus, the individual with the greatest number of loci available
      loci.remaining<-loci.best[[i]] #creates a list of loci to be treated when lookingfor loci missing in the best individual but possibly present in another
      repeat{ #iterates over loci
        remove.ID<-lapply(loci.inds,FUN=function(x){!any(x%in%best.in.group)}) #in which loci occurs the best ID from the first locus?
        loci.deleted<-loci.remaining[which(names(loci.remaining)%in%names(remove.ID[remove.ID==FALSE])==TRUE)]
        loci.remaining<-loci.remaining[which(names(loci.remaining)%in%names(remove.ID[remove.ID==FALSE])==FALSE)] #remove these loci ######### assumes that there between any pair of loci an ID is shared / test for that is needed

        if(length(unlist(loci.deleted))==0){
          warning(paste("representative.alignment: no individual for any of the loci (",paste(names(loci.remaining),collapse = " ",sep=" "),") in group", names(loci.best)[i],", group drops out."))
        }# should warn if in no individual has any of the loci in the group
        
        if(
                   length(names(loci.remaining))>0 #check that some loci remain to be processed
                  &length(intersect(unlist(lapply(lengths.ordered.groups[[i]][unlist(lapply(lengths.ordered.groups[[i]],FUN=function(x){any(names(x)%in%best.in.group)}))],FUN=function(z){names(z)})), unlist(loci.best[[i]][unlist(lapply(loci.best[[i]],FUN=function(x){!any(x%in%best.in.group)}))]) ) ) ==0 #checks that no sample from the best set of the remaining loci (second part) is present in the already processed loci (first part)
                  &length(unlist(loci.best[[i]][unlist(lapply(loci.best[[i]],FUN=function(x){!any(x%in%best.in.group)}))]))>0 #checks that the remaining loci are not empty
          ){
          warning(paste("representative.alignment: no shared individual between loci (",paste(names(loci.remaining),collapse = " ",sep=" "),") and (",paste(names(loci.deleted),collapse = " ",sep=" "),") in group", names(loci.best)[i],"."))
          }# should warn if in any group two sets of individuals not overlapping in any locus are created
       
        if(length(loci.remaining)==0|length(unlist(loci.remaining))==0){break} #if nothing's left, stop

        xx<-table(unlist(loci.remaining))
        best.in.group<-c(best.in.group,names(xx[which(xx==max(xx))])[1]) #add the best of the remaining loci to the selection
      }
      best[[i]]<-best.in.group
    }

    best<-unlist(best)


  } #end if maxloci==TRUE

  #################################

  best<-unique(best)

  #extract best individuals from the alignments
  best.alignments<-lapply(alignments,FUN=function(m){selector(best,m,save=FALSE)})
  best.alignments<-lapply(best.alignments,FUN=function(n){n[n!=""]})

  #concatenate the best individuals
  final.alignment<-conc.iso(best.alignments,save=save,file = file)

  #prepare final output
  results<-list(per.group.inds,lengths.ordered.groups,missing.loci,best.alignments,final.alignment)
  names(results)<-c("per.group.inds","lengths.ordered.groups","missing.loci","best.alignments","final.alignment")

  #save single gene alignments
  if(save_single==TRUE){
    for(i in 1:length(best.alignments)){
      write.table(best.alignments[i],file=paste0(names(best.alignments)[i],".fas"),quote=FALSE,col.names=FALSE,row.names=FALSE)
    }
  }

  #return output
  return(results)

} #end of function


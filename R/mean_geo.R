#' Reconstructing rough geographic locations for internal nodes in a tree based on georeferenced tips.
#'
#' The function was developed for a very crude but fast approximation of the continuous phylogeographic analysis implemented in BEAST 1.10
#' @param tree phylo. Rooted phylogenetic tree, loaded as the ape's 'phylo' class. For the 'mean' method the tree should be ultrametric, but it would not fail.
#' @param tab data.frame. Table containing the tip labels and the corresponding GPS coordinates, with column names as required for BEAST 1.10: "taxa", "lat", "lon".
#' @param method character. Name of the method used: "mean" calculates node location as located in between the locations of its daughter nodes, but proportionally closer to that with the shorter branch; "proximity" assumes common origin where sister clades live near each other. For details see the relevant publication (not yet available), but the "mean" method is biologically less realistic. See Details.
#' @param threshold numeric. Number in the range <0,1>, specifying the support value above which there is a as good evidence for a given branch in the tree (default 0.95).
#' @param weighted.proximity Boolean. Whether to average 'proximity' witm 'mean' method or not when used with method="support_weighted".
#' @importFrom geosphere distVincentyEllipsoid bearing geodesic
#' @details The function is meant to be a rough approximation useful to get a first idea about the data where running a continuous phylogeographic analysis in BEAST 1.x.x would be unnecessarily time consuming. It is not meant to be a replacement of that method.
#' @details The method 'mean' calculates node location as located in between the locations of its two daughter nodes, but proportionally closer to that with the shorter branch ("weighted mean" thereafter).
#' @details The following methods always consider a subtree with 3-4 tips: a parent branch, its two daughter branches, and the branches descending from one (if the other is a terminal of the whole tree) or both these daughter nodes.
#' @details The method 'proximity' takes the daughter nodes and calculates the geographic distances between nodes descending from one of them on the one hand and the other on the other hand. These are then compared. The location of the daughter node is set to be the location of that of its descendants, which is is geographically the closest to the closest descendant of the other daughter node.
#' @details The method 'support_weighted' combines the two above (and a simple mean) in three variants depending on the support for the subtree (parental node, its two daughter nodes, and their daughter nodes) evaluated.
#' @details In variant 1, the tree is resolved enough for the proximity method to make sense: we can be sure that we are considering nodes that are related. The 'proximity' method works best if after split between the daughter lineages these expand in opposite directions. It is preferred, because it handles cases on unequal range sizes or non-adjacent ranges of the daughter clades better. If weighted.proximity==TRUE, the locations for the daughter nodes are calculated also using the 'simple mean' or 'weighted mean' method and the result is a mean of this and the 'proximity' method results, weighted by branch lengths. The branch length between the parent and the daughter node is compared to the mean branch length of the two branches descending from the daughter node. The shorter the former compared to the latter, the greater weight is given to the 'proximity' estimate, because less time for range change elapsed before the daughter branches split further, so their close relationship is more informative. If the ratio is the opposite, the locations of the nodes descending from them are more informative than their common descent from the parent node of the subtree.
#' @details In variant 2, we do not have enough certainty regarding the subtree topology (relationships, branch lengths), so a simple mean is calculated as the best estimate when we do not even known if the daughter node is supported.
#' @details In variant 3, we cannot be sure that the two daughter nodes are related, but they themselves are well supported. So for calculating the locations of the daughter nodes we can reasonably take into account also the lengths of there descendant branches and the 'weighted mean' method is applied. If only one daughter node is supported, 'simple mean' is used for the other.
#' @export
#' @return Returns a data frame containing tree branches and the associated geographic locations of their terminal nodes.

#version 07.07.2022


mean_geo<-function(tree,tab,method=c("mean","proximity","support_weighted"),threshold=0.90,weighted.proximity=FALSE){

method<-match.arg(method)
#tab<-read.table("localities_atrolabiata.csv",header=TRUE) #tab must have header like for BEAST, ie. columns "taxa", "lon", "lat"

#extract info from the tree and tab
branches<-tree$edge
branches<-rbind(branches, c(0,length(tree$tip.label)+1))	#adds the root, assumes that the root is the first node numbered after the tips
branches<-cbind(branches,rep(NA,dim(branches)[1]),rep(NA,dim(branches)[1]),rep(NA,dim(branches)[1]),rep(NA,dim(branches)[1]),rep(NA,dim(branches)[1]),rep(NA,dim(branches)[1]))	#add columns for branch terminal node height, location 1&2
colnames(branches)<-c("parent","terminal","terminal.height","lat","lon","tip","support","length")
branches<-data.frame(branches)

#a<-1
#repeat{
#  branches$terminal.height[a]<-tree$annotations[[a]]$height
#  a<-a+1
#  if(a==dim(branches)[1]){break}
#}
#branches$terminal.height[a]<-tree$root.annotation$height

#get node depths
x<-ape::node.depth.edgelength(tree)
x<-cbind(x,c(1:length(x)))
for(i in 1:length(x)){
  branches$terminal.height[which(branches$terminal==i)]<-x[i]
}

#calculate branch lengths ############################ run tests!!!!!!!
a<-1
repeat{
  branches$length[a]<-abs(branches$terminal.height[a]-branches$terminal.height[which(branches$terminal==branches$parent[a])])
  a<-a+1
  if(a==dim(branches)[1]){break}
}

#get tips
a<-1
repeat{
  branches$tip[which(branches$terminal==a)]<-tree$tip.label[a]
  # branches$tip[which(branches$terminal==a)]<-substr(tree$tip.label[a],1,8) ######substring only for my data, otherwise take it out
  a<-a+1
  if(a==length(tree$tip.label)+1){break}
}

#get tip locations
tips<-tree$tip.label
a<-1
repeat{
  branches$lon[which(branches$tip==tips[a])]<-tab$lon[which(tab$taxa==tips[a])]
  branches$lat[which(branches$tip==tips[a])]<-tab$lat[which(tab$taxa==tips[a])]
  a<-a+1
  if(a==length(tips)+1){break}
}

#get branch supports
if(method=="support_weighted"){
  print("This method requires support values as posterior probabilities (named 'posterior'), read in with read.annotated.nexus. Succesful reading of the values will be confirmed before running the calculations.")
a<-1
repeat{
  if(length(tree$annotations[[a]]$posterior)>0){
    branches$support[a]<-tree$annotations[[a]]$posterior
  }
  a<-a+1
  if(a==dim(branches)[1]){break}
}
branches$support[a]<-0
branches$support[branches$terminal<=length(tree$tip.label)]<-1
print("Reading supports finished")
}


############## method 1
if(method=="mean"){
repeat{
  if(!any(is.na(branches[,4]))){break} #ends loop if all values are filled
   #checks that the branch has the necessary data
  a<-length(tree$tip.label)+1
  repeat{
    parent<-which(branches[,2]==a) #takes node
    offspring<-which(branches[,1]==branches[parent,2]) #finds its two daugthers
    test1<-!any(is.na(branches[offspring,4]))
    test2<-is.na(branches[parent,4])
    if(test1&test2){
      #branches[parent,4]<-1000 #used for testing
      loc1<-c(branches[offspring[1],5],branches[offspring[1],4])
      loc2<-c(branches[offspring[2],5],branches[offspring[2],4])
      branch1<-abs(branches[parent,3]-branches[offspring[1],3])
      branch2<-abs(branches[parent,3]-branches[offspring[2],3])
      #calculates the location at baseline time as a point on the ellipsoid path from origin to terminal at distance given by the position of baseline relative to the branch length
      distance.meters<-geosphere::distVincentyEllipsoid(loc1,loc2)
      distance.azi<-geosphere::bearing(loc1,loc2) #calculates initial bearing from loc1 towards loc2
      ratio<-branch1/(branch2+branch1) #the ratio of the branch lengths to the offspring nodes
      distance.share<-distance.meters*ratio #the distance from offspring 1 to the parent
      parent.location<-geosphere::geodesic(p=loc1,azi=distance.azi,d=distance.share)
      branches[parent,4]<-parent.location[2]
      branches[parent,5]<-parent.location[1]
      print(a)
      }
    a<-a+1
    if(a>length(branches[,1])){break}
  }
} #end main repeat of method 1
} #end method 1

############## method 2
if(method=="proximity"){
  repeat{
    if(!any(is.na(branches[branches$parent!=0,4]))){break} #ends loop if all values are filled except for the root
    #checks that the branch has the necessary data
    a<-length(tree$tip.label)+1
    repeat{
      subtree<-ape::extract.clade(tree,a)
      parent<-which(branches[,2]==a) #takes node
      offspring<-which(branches[,1]==branches[parent,2]) #finds its two daughters
      offspring2a<-which(branches[,1]==branches[offspring[1],2]) #row in branches
        if(length(offspring2a)==0){offspring2a<-offspring[1]}
      offspring2b<-which(branches[,1]==branches[offspring[2],2])
        if(length(offspring2b)==0){offspring2b<-offspring[2]}


      # make tests to see if this node is suitable for the processing
      test1<-!any(is.na(branches[offspring2a,4]))&!any(is.na(branches[offspring2b,4])) #controls that I have all the necessary data for the node
      test2<-any(is.na(branches[offspring,4])) #controls that the node is not already done
      test3<-length(subtree$tip.label)>=3 #controls that it has at least three daughter nodes
      #dopsat podminku, aby se to zastavilo i kdyz neni spocitano pro root
      if(test1&test2&test3){
        #ted musin najit pro oba offspring jejich offspring, spocitat vydalenosti a vybrat nejbliydi par
        distances<-cbind(expand.grid(offspring2a,offspring2b,stringsAsFactors=FALSE),rep(NA,length(offspring2a)*length(offspring2b)))
        colnames(distances)<-c("offspring2a","offspring2b","distance")

        for(i in 1:length(distances[,1])){
         loc1<-c(branches[distances$offspring2a[i],5],branches[distances$offspring2a[i],4])
         loc2<-c(branches[distances$offspring2b[i],5],branches[distances$offspring2b[i],4])
         distances$distance[i]<-geosphere::distVincentyEllipsoid(loc1,loc2)
        }
        b<-which(distances$distance==min(distances$distance))[1]  #which point pair from bith offsping clades is the closest? takes the first if there is a tie

        #now put data back to branches
        branches$lat[offspring[1]]<-branches$lat[distances$offspring2a[b]]
        branches$lat[offspring[2]]<-branches$lat[distances$offspring2b[b]]
        branches$lon[offspring[1]]<-branches$lon[distances$offspring2a[b]]
        branches$lon[offspring[2]]<-branches$lon[distances$offspring2b[b]]

        print(a)
        }
      a<-a+1
      if(a>length(branches[,1])){break}
    }
  }
  } #method 2 end

#################method 3
  if(method=="support_weighted"){

    #functions for var2 and var3 (simple and weighted "mean")
    mean_location<-function(loc1,loc2,branches,m.method=c("simple.mean","weighted.mean"),which_offspring=NULL) {
      distance.meters<-geosphere::distVincentyEllipsoid(loc1,loc2)
      distance.azi<-geosphere::bearing(loc1,loc2) #calculates initial bearing from loc1 towards loc2
      if(m.method=="simple.mean"){ratio<-0.5}
      if(m.method=="weighted.mean"){
        if(which_offspring==1){
          branch1<-abs(branches$terminal.height[offspring[1]]-branches$terminal.height[offspring2a[1]])
          branch2<-abs(branches$terminal.height[offspring[1]]-branches$terminal.height[offspring2a[2]])
        }
        if(which_offspring==2){
          branch1<-abs(branches$terminal.height[offspring[2]]-branches$terminal.height[offspring2b[1]])
          branch2<-abs(branches$terminal.height[offspring[2]]-branches$terminal.height[offspring2b[2]])
        }
        #calculates the location as a point on the ellipsoid path from origin to terminal at distance given by the position of baseline relative to the branch length
        ratio<-branch1/(branch2+branch1)
      }
      distance.share<-distance.meters*ratio #the distance from offspring 1 to the parent
      parent.location<-geosphere::geodesic(p=loc1,azi=distance.azi,d=distance.share)
      return(parent.location)
    }

    repeat{
      if(!any(is.na(branches[branches$parent!=0,4]))){break} #ends loop if all values are filled except for the root
      #checks that the branch has the necessary data

      a<-length(tree$tip.label)+1 #a gives you a node number, variables parent, offspring and offspring2a+offspring2b give row(s) in the branches dataframe
      #offspring2a is daughter of offspring[1] and offspring2b is daughter of offspring[2]

      repeat{
        subtree<-ape::extract.clade(tree,a)
        parent<-which(branches$terminal==a) #takes node, the variable stores the row in branches where is its terminal
        offspring<-which(branches$parent==branches$terminal[parent]) #finds its two daughters
        offspring2a<-which(branches$parent==branches$terminal[offspring[1]]) #row in branches with the offspring of first offspring
        if(length(offspring2a)==0){offspring2a<-offspring[1]}
        offspring2b<-which(branches$parent==branches$terminal[offspring[2]])
        if(length(offspring2b)==0){offspring2b<-offspring[2]}

        # make tests to see if this node is suitable for the processing
        test1<-!any(is.na(branches[offspring2a,4]))&!any(is.na(branches[offspring2b,4])) #controls that I have all the necessary data for the node
        test2<-any(is.na(branches[offspring,4])) #controls that the node is not already done
        test3<-length(subtree$tip.label)>=3 #controls that it has at least three daughter nodes

        if(test1&test2&test3){
          #notes supports for parent and offspring to decide which variant to use on that subtree
          supp.parent<-branches$support[parent]
          supp.offspring.a<-branches$support[offspring[1]]
          supp.offspring.b<-branches$support[offspring[2]]

          #tests
          var1test1<-(supp.parent>=threshold&(supp.offspring.a>=threshold|supp.offspring.b>=threshold)) #checks that the support of parent and at least one offspring are over threshold
          var2test1<-supp.offspring.a<threshold & supp.offspring.b<threshold #testing that both offspring branches are unsupported (below threshold)
          var2test2<-length(offspring2a)<2 & supp.offspring.b<threshold & supp.parent<threshold #testing that parent is unsupported and if offspring a is supported it is a terminal
          var2test3<-length(offspring2b)<2 & supp.offspring.a<threshold & supp.parent<threshold #testing that, when offspring b (i.e. offspring[2]) is terminal, the other one is unsupported
          var3test1<-(supp.parent<threshold)

          #variant 1 always proximity_______________________________________________________________________________

          if(var1test1){ #if supports variant 1

            #creates table of distances between all combinations of offspring2a and offspring2b
            distances<-cbind(expand.grid(offspring2a,offspring2b,stringsAsFactors=FALSE),rep(NA,length(offspring2a)*length(offspring2b)))
            colnames(distances)<-c("offspring2a","offspring2b","distance")

            for(i in 1:length(distances[,1])){
              loc1<-c(branches[distances$offspring2a[i],5],branches[distances$offspring2a[i],4])
              loc2<-c(branches[distances$offspring2b[i],5],branches[distances$offspring2b[i],4])
              distances$distance[i]<-geosphere::distVincentyEllipsoid(loc1,loc2)
            }
            b<-which(distances$distance==min(distances$distance))[1]  #which point pair from both offspring clades is the closest? takes the first if there is a tie

            #records the offspring locations calculated with the proximity method
            offspring1lat.p<-branches$lat[distances$offspring2a[b]]
            offspring2lat.p<-branches$lat[distances$offspring2b[b]]
            offspring1lon.p<-branches$lon[distances$offspring2a[b]]
            offspring2lon.p<-branches$lon[distances$offspring2b[b]]

            if(weighted.proximity==TRUE){ ##########default off, trying what it does

            #calculates the offspring locations using the weighted or simple mean method
            if(length(offspring2a)>1&branches$support[offspring[1]]>=threshold){ #check that offspring 1 branches and is above threshold
              loc1<-c(branches$lon[offspring2a[1]],branches$lat[offspring2a[1]])
              loc2<-c(branches$lon[offspring2a[2]],branches$lat[offspring2a[2]])
              parent.location<-mean_location(loc1,loc2,branches,m.method = "weighted.mean",which_offspring=1)
              offspring1lat.m<-parent.location[2]
              offspring1lon.m<-parent.location[1]
            }
            if(length(offspring2b)>1&branches$support[offspring[2]]>=threshold){ #check that the other offspring 1 branches  and is above threshold
              loc1<-c(branches$lon[offspring2b[1]],branches$lat[offspring2b[1]])
              loc2<-c(branches$lon[offspring2b[2]],branches$lat[offspring2b[2]])
              parent.location<-mean_location(loc1,loc2,branches,m.method = "weighted.mean",which_offspring=2)
              offspring2lat.m<-parent.location[2]
              offspring2lon.m<-parent.location[1]
            }

           #if the offspring branch is not supported, calculate simple mean
            if(length(offspring2a)>1&branches$support[offspring[1]]<threshold){ #check that offspring 1 branches and is not supported
              loc1<-c(branches$lon[offspring2a[1]],branches$lat[offspring2a[1]])
              loc2<-c(branches$lon[offspring2a[2]],branches$lat[offspring2a[2]])
              parent.location<-mean_location(loc1,loc2,branches,m.method = "simple.mean",which_offspring=1)
              offspring1lat.m<-parent.location[2]
              offspring1lon.m<-parent.location[1]
            }
            if(length(offspring2b)>1&branches$support[offspring[2]]<threshold){ #check that offspring 1 branches and is not supported
              loc1<-c(branches$lon[offspring2b[1]],branches$lat[offspring2b[1]])
              loc2<-c(branches$lon[offspring2b[2]],branches$lat[offspring2b[2]])
              parent.location<-mean_location(loc1,loc2,branches,m.method = "simple.mean",which_offspring=2)
              offspring2lat.m<-parent.location[2]
              offspring2lon.m<-parent.location[1]
            }


            #calculating the position between the proximity-estimated and mean-estimated offspring position
            #now for the first offspring branch
            if(length(offspring2a)>1){ #check that offspring a branches
              loc1<-c(offspring1lon.m,offspring1lat.m) #loc1 is what is estimated by the mean method, preferred when branch.offspring is long and branch.offspring2 is short
              loc2<-c(offspring1lon.p,offspring1lat.p)
              distance.meters<-geosphere::distVincentyEllipsoid(loc1,loc2)
              distance.azi<-geosphere::bearing(loc1,loc2) #calculates initial bearing from loc1 towards loc2
              branch.offspring2<-(branches$length[offspring2a[1]]+branches$length[offspring2a[2]])/2
              branch.offspring<-branches$length[offspring[1]]
              ratio<-branch.offspring2/(branch.offspring+branch.offspring2) #ratio is large when offspring branch is short relative to offspring2, meaning that we get far from loc1 (mean method estimate) towards loc2 (proximity method estimate)
              distance.share<-distance.meters*ratio #the distance from offspring 1 to the parent
              parent.location<-geosphere::geodesic(p=loc1,azi=distance.azi,d=distance.share)
              offspring1lat.p<-parent.location[2]
              offspring1lon.p<-parent.location[1]
            }

            if(length(offspring2b)>1){ #check that offspring b branches
              loc1<-c(offspring2lon.m,offspring2lat.m)
              loc2<-c(offspring2lon.p,offspring2lat.p)
              distance.meters<-geosphere::distVincentyEllipsoid(loc1,loc2)
              distance.azi<-geosphere::bearing(loc1,loc2) #calculates initial bearing from loc1 towards loc2
              branch.offspring2<-(branches$length[offspring2b[1]]+branches$length[offspring2b[2]])/2
              branch.offspring<-branches$length[offspring[2]]
              ratio<-branch.offspring2/(branch.offspring+branch.offspring2) #ratio is large when offspring branch is short relative to offspring2, meaning that we get far from loc1 (mean method estimate) towards loc2 (proximity method estimate)
              distance.share<-distance.meters*ratio #the distance from offspring 1 to the parent
              parent.location<-geosphere::geodesic(p=loc1,azi=distance.azi,d=distance.share)
              offspring2lat.p<-parent.location[2]
              offspring2lon.p<-parent.location[1]
            }

            } #end if for weighted.proximity ####in testing

            #now put data back to branches data.frame
            branches$lat[offspring[1]]<-offspring1lat.p
            branches$lat[offspring[2]]<-offspring2lat.p
            branches$lon[offspring[1]]<-offspring1lon.p
            branches$lon[offspring[2]]<-offspring2lon.p

          } #end if supports variant 1

          #variant 2 always simple mean_______________________________________________________________________________
          if(var2test1|var2test2|var2test3){ #if supports variant 2
            if(length(offspring2a)>1){ #check that offspring a branches
              #simplemean<-function(loc1,loc2){} #make it a function, output will be parent.location
              #simplemean(loc1,loc2)
              loc1<-c(branches$lon[offspring2a[1]],branches$lat[offspring2a[1]])
              loc2<-c(branches$lon[offspring2a[2]],branches$lat[offspring2a[2]])
              parent.location<-mean_location(loc1,loc2,branches,m.method = "simple.mean",which_offspring=1)
              branches$lat[offspring[1]]<-parent.location[2]
              branches$lon[offspring[1]]<-parent.location[1]
            }
            if(length(offspring2b)>1){ #check that offspring b branches
              loc1<-c(branches$lon[offspring2b[1]],branches$lat[offspring2b[1]])
              loc2<-c(branches$lon[offspring2b[2]],branches$lat[offspring2b[2]])
              parent.location<-mean_location(loc1,loc2,branches,m.method = "simple.mean",which_offspring=2)
              branches$lat[offspring[2]]<-parent.location[2]
              branches$lon[offspring[2]]<-parent.location[1]
            }
          } #end if supports variant 2

          #variant 3 - chooses between weighted mean and simple mean_______________________________________________________________________________

          if(var3test1){ #if supports variant 3 / var3 was defined by no support for parental branch, but when there is no support for one of the branhces, move it to var2
            #if the offspring branch is not supported, calculate weighted mean

            if(length(offspring2a)>1&branches$support[offspring[1]]>=threshold){ #check that offspring 1 branches and is above threshold
              loc1<-c(branches$lon[offspring2a[1]],branches$lat[offspring2a[1]])
              loc2<-c(branches$lon[offspring2a[2]],branches$lat[offspring2a[2]])
              parent.location<-mean_location(loc1,loc2,branches,m.method = "weighted.mean",which_offspring=1)
              branches$lat[offspring[1]]<-parent.location[2]
              branches$lon[offspring[1]]<-parent.location[1]
            }
            if(length(offspring2b)>1&branches$support[offspring[2]]>=threshold){ #check that the other offspring 1 branches  and is above threshold
              loc1<-c(branches$lon[offspring2b[1]],branches$lat[offspring2b[1]])
              loc2<-c(branches$lon[offspring2b[2]],branches$lat[offspring2b[2]])
              parent.location<-mean_location(loc1,loc2,branches,m.method = "weighted.mean",which_offspring=2)
              branches$lat[offspring[2]]<-parent.location[2]
              branches$lon[offspring[2]]<-parent.location[1]
            }
            #if the offspring branch is not supported, calculate simple mean
            if(length(offspring2a)>1&branches$support[offspring[1]]<threshold){ #check that offspring 1 branches
              loc1<-c(branches$lon[offspring2a[1]],branches$lat[offspring2a[1]])
              loc2<-c(branches$lon[offspring2a[2]],branches$lat[offspring2a[2]])
              parent.location<-mean_location(loc1,loc2,branches,m.method = "simple.mean",which_offspring=1)
              branches$lat[offspring[1]]<-parent.location[2]
              branches$lon[offspring[1]]<-parent.location[1]
            }
            if(length(offspring2b)>1&branches$support[offspring[2]]<threshold){ #check that offspring 1 branches
              loc1<-c(branches$lon[offspring2b[1]],branches$lat[offspring2b[1]])
              loc2<-c(branches$lon[offspring2b[2]],branches$lat[offspring2b[2]])
              parent.location<-mean_location(loc1,loc2,branches,m.method = "simple.mean",which_offspring=2)
              branches$lat[offspring[2]]<-parent.location[2]
              branches$lon[offspring[2]]<-parent.location[1]
            }

          } #end if supports variant 3

          #for development - check that all possible combinations of supports are covered
          if(!var1test1&!var2test1&!var2test2&!var2test3&!var3test1) {stop(paste(parent, "tadytahle kombinace ti ve funkci chybi"))}

        print(a)

        }

        a<-a+1

        ################put here a clean-up step, removing all values generated so that they cannot be recycled if something fails
        #if(a==262){break} #for development only!!!!!!!!!!!!!!!!!!!!
        if(a>length(branches[,1])){break} # could be temporarily disabled for development
        } #end inner repeat

    }  #end outer repeat
  } #method 3 end
#################end methods

# calculate mean for root
if(method=="support_weighted"){
  print("Calculating location for root:")
  parent<-which(branches[,1]==0) #finds row with the root node in branches
  print(parent)
  offspring<-which(branches[,1]==branches[parent,2]) #finds its two daughters
  print(offspring)
  test4<-!any(is.na(branches[offspring,4]))
  if(test4){
    #branches[parent,4]<-1000 #used for testing
    loc1<-c(branches[offspring[1],5],branches[offspring[1],4])
    loc2<-c(branches[offspring[2],5],branches[offspring[2],4])
    supp.offspring.ra<-branches$support[offspring[1]]
    print(branches$support[offspring[1]])
    supp.offspring.rb<-branches$support[offspring[2]]
    print(branches$support[offspring[2]])
    if(supp.offspring.ra>=threshold & supp.offspring.rb>=threshold){
      distance.meters<-geosphere::distVincentyEllipsoid(loc1,loc2)
      distance.azi<-geosphere::bearing(loc1,loc2) #calculates initial bearing from loc1 towards loc2
      branch1<-abs(branches$terminal.height[parent]-branches$terminal.height[offspring[1]])
      branch2<-abs(branches$terminal.height[parent]-branches$terminal.height[offspring[2]])
      ratio<-branch1/(branch2+branch1)
      distance.share<-distance.meters*ratio #the distance from offspring 1 to the parent
      parent.location<-geosphere::geodesic(p=loc1,azi=distance.azi,d=distance.share)
    }else{
      parent.location<-mean_location(loc1,loc2,branches,m.method = "simple.mean")
    }
    branches[parent,4]<-parent.location[2]
    branches[parent,5]<-parent.location[1]
  }
}else{
    loc1<-c(branches[offspring[1],5],branches[offspring[1],4])
    loc2<-c(branches[offspring[2],5],branches[offspring[2],4])
    distance.meters<-geosphere::distVincentyEllipsoid(loc1,loc2)
    distance.azi<-geosphere::bearing(loc1,loc2) #calculates initial bearing from loc1 towards loc2
    branch1<-abs(branches$terminal.height[parent]-branches$terminal.height[offspring[1]])
    branch2<-abs(branches$terminal.height[parent]-branches$terminal.height[offspring[2]])
    ratio<-0.5
    distance.share<-distance.meters*ratio #the distance from offspring 1 to the parent
    parent.location<-geosphere::geodesic(p=loc1,azi=distance.azi,d=distance.share)
    branches[parent,4]<-parent.location[2]
    branches[parent,5]<-parent.location[1]
} #end if


return(branches)
}#function end

#tree<-read.annotated.nexus("vind1f.tree")
#tab<-read.table("vind-tabulka1.csv",header=TRUE)

setwd("D:/TXT/Work/Recherche/PUBLICATIONS/2023-1 - Hand et al. - Vielase bat/Vielasia-sigei_Hand-et-al/Statistical treatment")
invisible(lapply(c("ape","caper","doParallel","foreach","geiger","geomorph","lattice","LOST","MASS","mda","nnet","parallel","phylolm","phytools","scales","sensiPhy","sf","treeio","ULT"), library, character.only = TRUE))
setwd("../test")
# load("D:/TXT/Work/Recherche/PUBLICATIONS/2023-1 - Hand et al. - Vielase bat/Vielasia-sigei_Hand-et-al/Statistical treatment/Vielasia stat analyses.RData")

# Custom functions (with additional details in the Supp Data 1 for the align_rep_ppca and conv.pgls in their script part as #comments)
{
  # Simple function to replace a phylogenetic subtree by another phylogenetic tree
  subs.tree.to.ref<-function(back.tree,front.tree,front.tree.subset=NULL){
    if(length(front.tree)==4&class(front.tree)!="multiPhylo"){
      front.tree<-list(as.phylo(front.tree))
    }
    
    if(all(!is.na(front.tree.subset)&!is.null(front.tree.subset))){
      back.tree<-keep.tip(back.tree,front.tree.subset[which(front.tree.subset%in%back.tree$tip.label)])
      for (i in 1:length(front.tree)){
        front.tree[[i]]<-keep.tip(front.tree[[i]],front.tree.subset[which(front.tree.subset%in%front.tree[[i]]$tip.label)])
      }
    }
    
    common.clade<-getClade(back.tree,which(back.tree$tip.label%in%front.tree[[1]]$tip.label))
    
    back.tree.ages<-max(nodeHeights(back.tree))-nodeHeights(back.tree)
    back.stem.common.age<-back.tree.ages[which(back.tree$edge[,2]==common.clade),2]
    
    conc.topo<-list()
    if(length(front.tree)>1){
      class(conc.topo)<-"multiPhylo"
    }
    else{
      class(conc.topo)<-"phylo"
    }
    j<-1
    for (i in 1:length(front.tree)){
      curr.front.stem.age<-max(nodeHeights(front.tree[[i]]))
      if(curr.front.stem.age>back.stem.common.age){next}
      else{
        conc.topo[[j]]<-subs.tree(back.tree,front.tree[[i]],common.clade,node.age=curr.front.stem.age)
        j<-j+1
      }
    }
    return(conc.topo)
  }
  
  Arbour_et_al_treatment<-function(lmk_data,specimen_species,pairs,slides,omit.lmk=NA){
    n_specimens<-dim(lmk_data)[1]
    n_lmk<-dim(lmk_data)[2]/3
    col_lmk<-setNames(c(1:dim(lmk_data)[2]),as.character(rep(c(1:n_lmk),each=3)))
    geomorph_lmk<-arrayspecs(lmk_data,n_lmk,3,n_specimens)
    
    for (i in 1:n_specimens){
      if(!any(is.na(geomorph_lmk[,,i]))){next}
      else{
        geomorph_lmk[,,i]<-flipped(specimen=geomorph_lmk[,,i],land.pairs=pairs)
      }
    }
    
    if(!any(is.na(omit.lmk))){
      col_lmk<-col_lmk[-which(as.numeric(names(col_lmk))%in%omit.lmk)]
      old_n_lmk<-n_lmk
      n_lmk<-n_lmk-length(omit.lmk)
      
      pairs<-pairs[!apply(apply(pairs,c(1,2),function(x){x%in%omit.lmk}),1,any),]
      for(i in 1:nrow(pairs)){
        for(j in 1:ncol(pairs)){
          pairs[i,j]<-pairs[i,j]-length(which(pairs[i,j]>omit.lmk))
        }
      }
      
      geomorph_lmk<-arrayspecs(matrix(two.d.array(geomorph_lmk),nrow=n_specimens,ncol=old_n_lmk*3,dimnames = list(rownames(lmk_data),colnames(lmk_data)))[,col_lmk],n_lmk,3,n_specimens)
      slmk<-setNames(unlist(slides),c(1:length(unlist(slides))))
      slmk<-sort(slmk)
      if(any(!omit.lmk%in%slmk)){
        for(i in 1:length(omit.lmk)){
          slmk[slmk>omit.lmk[i]]<-slmk[slmk>omit.lmk[i]]-1
        }
      }
      bounds<-c()
      for(i in 1:length(slides)){
        bounds[2*i-1]<-ifelse(i==1,0,bounds[2*(i-1)])+1
        bounds[2*i]<-bounds[2*i-1]+length(slides[[i]])-1
      }
      omit.slmk<-slmk[!slmk%in%omit.lmk]
      omit.slmk[1]<-setNames(slmk[1],omit.slmk[1])
      for(i in 2:length(omit.slmk)){
        if((omit.slmk[i]-omit.slmk[i-1])>1){
          omit.slmk[i:length(omit.slmk)]<-omit.slmk[i:length(omit.slmk)]-1
        }
      }
      for(i in 1:length(slides)){
        slides[[i]]<-unname(omit.slmk[as.numeric(names(omit.slmk))>=bounds[2*i-1]&as.numeric(names(omit.slmk))<=bounds[2*i]])
      }
    }
    
    nslides<-length(unlist(slides))-2*length(slides)
    curves<-matrix(ncol=3,nrow=0,NA)
    for (i in 1:length(slides)){
      for (j in 2:(length(slides[[i]])-1)){
        curves<-rbind(curves,t(c(slides[[i]][j-1],slides[[i]][j],slides[[i]][j+1])))
      }
    }
    colnames(curves)<-c("before","slide","after")
    
    specimens_gpa<-gpagen(MissingGeoMorph(geomorph_lmk,method="BPCA"),ProcD=TRUE,curves=curves)
    
    species<-levels(factor(specimen_species))
    gpa<-matrix(nrow=length(species),ncol=n_lmk*3,NA)
    for (i in 1:length(species)){
      gpa[i,]<-c(t(mshape(specimens_gpa$coords[,,c(which(specimen_species==species[i]))])))
    }
    rownames(gpa)<-species
    
    def_gpa<-matrix(nrow=length(species),ncol=(n_lmk-length(pairs[,1]))*3,NA)
    for (i in 1:length(species)){
      curr_spec<-matrix(ncol=3,nrow=n_lmk,gpa[i,],byrow=T)
      curr_spec[,2]<-abs(curr_spec[,2])
      temp_def<-matrix(ncol=3,nrow=0,NA)
      for(j in 1:n_lmk){
        if(!j%in%pairs){
          temp_def<-rbind(temp_def,curr_spec[j,])
        }
        else{
          if(j%in%pairs[,2]){next}
          else{
            which_pair<-which(pairs[,1]==j)
            temp_def<-rbind(temp_def,apply(curr_spec[pairs[which_pair,],],2,mean))
          }
        }
      }
      temp_def<-c(t(temp_def))
      def_gpa[i,]<-temp_def
    }
    rownames(def_gpa)<-species
    
    return(def_gpa)
  }
  
  # General remark for the align_rep_ppca function.
  # After the pPCAs on Procrustes residuals, individuals are distributed along each component no matter the 
  # “direction” of that component (i.e., species are distributed from negative to positive values or the reverse 
  # depending on the analysis). For instance, the mostly retrieved orientation of the two first pPC axes of our pPCAs 
  # are inverted compared with those of Arbour et al. (2019). 
  # Thus, while the sign inversion itself is not an issue, each species is likely to have both positive and negative 
  # values on each principal component while gathering pPCA results over all phylogenies; averaging them would mostly 
  # yield close-to-zero values for all species. 
  # Therefore, before averaging values for each species and consider these pPCA results, a step is needed to check 
  # which axes are inverted for which iteration (= phylogenetic tree), and to de-invert them. 
  # This custom function (align_rep_ppca) deals with this issue and tentatively solves it.
  # There are comments for each step of this function within the code.
  
  align_rep_ppca<-function(data,threshold=seq(0,1,0.1),dec.prio=TRUE,print.threshold=FALSE,print.not.inverted=FALSE){
    threshold<-rev(threshold)
    not_inverted<-c(1:dim(data[[1]])[2])
    pPC_threshold<-matrix(ncol=length(not_inverted),nrow=3,c(rep(0,length(not_inverted)),rep(NA,length(not_inverted)*2)),byrow=TRUE)
    colnames(pPC_threshold)<-colnames(data[[1]])
    rownames(pPC_threshold)<-c("threshold","% non-fully agreeing taxa","mean value non-fully agreeing taxa")
    for(t in 1:(length(threshold)-1)){
      cl<-makeCluster(detectCores()-1)
      registerDoParallel(cl)
      
      # This function first considers the variation of the values for each species and for each pPC axis over all 
      # ‘iterations’ (i.e., all phylogenies considered). For a given species and a given pPC axis, one could expect 
      # the values to be floating around a given value and therefore around its negative counterpart, with an 
      # “intermediate blank range” between the two clouds of points. 
      # The align_rep_ppca function computes such intermediate range at the line with a ###1 (as being the gap between the closest 
      # positive and negative values to zero) and expresses it relatively to the whole range of values (i.e., the gap 
      # between the most positive and negative values) for each pPC axis and for each species. This relative range is 
      # necessarily comprised between zero (no intermediate range) and one (all points are exactly equal to one value 
      # and its negative counterpart). 
      # Such method is rather conservative because one extremely negative or positive or null point will impact, 
      # alone, the value of the proportion of intermediate blank range; high proportions of intermediate blank range 
      # will therefore reveal two clearly separated clouds of points.
      
      diff<-foreach(i=1:dim(data[[1]])[1],.combine="rbind",.packages="foreach")%dopar%{
        foreach(j=not_inverted,.combine="c")%dopar%{
          local<-mapply(function(x){x[i,j]},data)
          if(all(local>0)|all(local<0)){
            NA
          }
          else{
            (min(local[local>0])-max(local[local<0]))/(range(local)[2]-range(local)[1]) ###1
          }
        }
      }
      stopCluster(cl)
      
      colnames(diff)<-not_inverted
      rownames(diff)<-rownames(data[[1]])
      
      # The function then takes a degressive “definitive” threshold value between one and zero: if a pPC is 
      # considered at a time, it will not be considered after, avoiding errors (and allowing the function to process 
      # faster and faster). 
      
      n_diff_over_threshold<-apply(diff,2,function(x){length(which(x>threshold[t]))})
      pPC_to_mod<-as.numeric(names(which(n_diff_over_threshold>(dim(diff)[1]/2))))
      prop_not_full<-numeric(length=length(pPC_to_mod))
      value_not_full<-numeric(length=length(pPC_to_mod))
      if(length(pPC_to_mod)>0){
        for(p in pPC_to_mod){
          
          # For each threshold, the function identifies (lines below), for each pPC, which species display an 
          # “intermediate blank range” over the considered threshold; such species are likely to show sign inversions 
          # on the given pPC, and are therefore selected. 
          
          tax_ref<-diff[which(diff[,which(colnames(diff)==p)]>threshold[t]),which(colnames(diff)==p)]
          red_dat<-matrix(ncol=length(data),nrow=0)
          for (r in names(tax_ref)){
            red_dat<-rbind(red_dat,mapply(function(x){x[r,p]},data))
          }
          rownames(red_dat)<-names(tax_ref)
          
          # To be able to go further, there must be more than the half of all species showing an inversion (which is 
          # sort of a majority rule). 
          # Then, the function (line just below) identifies which sign is mostly shared by each species on this pPC 
          # axis (i.e., considering all iterations): each species will generally have a negative or positive value 
          # for a given pPC axis, and the positive or negative values indicate iterations to invert. 
          
          major_sign<-apply(red_dat,1,function(x){ifelse(length(which(x>0))>length(which(x<0)),1,-1)})
          
          # The function (line just below) then multiplies all species values by its major sign in order to get 
          # (temporarily) mostly positive values for all species: values are multiplied by one for species that have 
          # mostly positive values and by minus one in the reverse case. 
          
          major_sign_red_dat<-red_dat*major_sign
          
          # Therefore, negative values indicate inversions (line below): negative values for species with mostly 
          # positive values have been multiplied by one and remain negative, and positive values for species with 
          # mostly negative values have been multiplied by minus one and become negative. 
          
          inversions_in_reps<-apply(major_sign_red_dat,2,function(x){length(which(x<0))})
          
          # These negative values, for all species, thus serve to identify which iterations should be 
          # inverted. All species however do not necessarily agree; the function only interprets (line below) as an 
          # inversion the iterations where more than the half of the considered species (i.e., not all species, only 
          # those that show an intermediate blank range above the currently considered threshold) agree, to avoid 
          # too “local” cases (i.e., iterations with only some species indicating an inversion). 
          
          rep_to_invert<-which(inversions_in_reps>(length(tax_ref)/2))
          
          prop_not_full[which(pPC_to_mod==p)]<-round(100*length(which(inversions_in_reps[rep_to_invert]!=length(tax_ref)))/length(tax_ref),2)
          if(prop_not_full[which(pPC_to_mod==p)]==0){
            value_not_full[which(pPC_to_mod==p)]<-NA
          }
          else{
            value_not_full[which(pPC_to_mod==p)]<-round(100*mean(inversions_in_reps[rep_to_invert][inversions_in_reps[rep_to_invert]!=length(tax_ref)])/length(tax_ref),2)
          }
          
          # Finally (lines below), all values are inverted for the given pPC on the iterations that apparently 
          # yielded inverted values. 
          
          for(r in rep_to_invert){
            data[[r]][,p]<--data[[r]][,p]
          }
        }
        not_inverted<-not_inverted[-which(not_inverted%in%pPC_to_mod)]
        pPC_threshold[1,pPC_to_mod]<-threshold[t]
        pPC_threshold[2,pPC_to_mod]<-prop_not_full
        pPC_threshold[3,pPC_to_mod]<-value_not_full
      }
      if(length(not_inverted)==0){break}
    }
    
    # The function returns the inputted datasets whose pPC axes have been modified for each iteration showing an 
    # apparent inversion, but it can also return the number of the pPC axes that have not been inverted (if there 
    # was a too narrow “intermediate range” between positive and negative values; this may imply that there was no 
    # inversion at all for all iterations) and the threshold used for each inverted pPC axis by turning to TRUE the 
    # parameter print.not.inverted and print.threshold respectively.
    
    if(!print.threshold&!print.not.inverted){
      return(data)
    }
    else{
      results<-list("data"=data)
      if(print.threshold){
        results<-c(results,list("threshold"=pPC_threshold))
      }
      if(print.not.inverted){
        results<-c(results,list("not.inverted"=not_inverted))
      }
      return(results)
    }
  }
  
  # Replicate operations like mean, sd (etc) over data from align_rep_ppca
  rep_ppca_operations<-function(data,ops){
    res<-list()
    for (k in 1:length(ops)){
      cl<-makeCluster(detectCores()-1)
      registerDoParallel(cl)
      res[[k]]<-foreach(i=1:dim(data[[1]])[1],.combine="rbind",.packages="foreach")%dopar%{
        foreach(j=1:dim(data[[1]])[2],.combine="c")%dopar%{
          temp_res<-do.call(ops[k],list(mapply(function(x){x[i,j]},data)))
          if(grepl("test",ops[k],fixed=TRUE)){
            temp_res<-temp_res[["p.value"]]
            final_res<-"n"
            if(temp_res<0.01){final_res<-"***"}
            else if(temp_res<0.05){final_res<-"**"}
            else if (temp_res<0.1){final_res<-"*"}
          }
          else{
            final_res<-temp_res
          }
          final_res
        }
      }
      stopCluster(cl)
      colnames(res[[k]])<-colnames(data[[1]])
      rownames(res[[k]])<-rownames(data[[1]])
    }
    names(res)<-ops
    return(res)
  }
  
  tree_phylm_full_results<-function(x,y,xynames,multiphy,models=c("BM", "OUrandomRoot","OUfixedRoot", "lambda", "kappa", "delta", "EB"),output.all.models=TRUE,output="all"){
    if(length(models)==1){
      output.all.models<-TRUE
    }
    
    if(output=="all"){
      output<-c("all.results","sensi.estimates","all.stats","residuals")
    }
    
    data<-data.frame(y,x)
    rownames(data)<-xynames
    
    cl<-makeCluster(detectCores()-1)
    registerDoParallel(cl)
    
    all.results<-foreach(i=1:length(models),.packages="foreach")%dopar%{
      foreach(j=1:length(multiphy),.packages="phylolm")%dopar%{
        phylolm(data[,1]~data[,2],phy=multiphy[[j]],data=data,model=models[i])
      }
    }
    
    sensi.estimates<-foreach(i=1:length(models),.packages="foreach",.combine="rbind")%dopar%{
      foreach(j=1:length(multiphy),.combine = "rbind")%dopar%{
        coeffs<-coefficients(summary(all.results[[i]][[j]]))
        data.frame(j,
                   coeffs[1,1],coeffs[1,2],coeffs[1,4],
                   coeffs[2,1],coeffs[2,2],coeffs[2,4],
                   models[i],all.results[[i]][[j]]$aic,
                   ifelse(models[i]=="BM",NA,all.results[[i]][[j]]$optpar))
      }
    }
    colnames(sensi.estimates)<-c("n.tree","intercept","se.intercept","pval.intercept","estimate","se.estimate","pval.estimate","model","aic","optpar")
    
    all.stats<-foreach(i=1:length(models))%dopar%{
      all.stats<-data.frame(min = apply(sensi.estimates[sensi.estimates$model==models[i],-c(1,8)], 2, min), 
                            max = apply(sensi.estimates[sensi.estimates$model==models[i],-c(1,8)], 2, max), 
                            mean = apply(sensi.estimates[sensi.estimates$model==models[i],-c(1,8)], 2, mean), 
                            sd_tree = apply(sensi.estimates[sensi.estimates$model==models[i],-c(1,8)], 2, sd))
      all.stats$CI_low <- all.stats$mean - qt(0.975, df = length(multiphy) - 1) * all.stats$sd_tree/sqrt(length(multiphy))
      all.stats$CI_high <- all.stats$mean + qt(0.975, df = length(multiphy) - 1) * all.stats$sd_tree/sqrt(length(multiphy))
      all.stats
    }
    
    residuals<-foreach(i=1:length(models))%dopar%{
      t(mapply(function(x){x$residuals[order(match(names(x$residuals),multiphy[[1]]$tip.label))]},all.results[[i]]))
    }
    names(all.results)<-names(all.stats)<-names(residuals)<-models
    
    stopCluster(cl)
    
    models.summary<-NULL
    
    if(!output.all.models){
      best.model<-which.min(mapply(function(x){x[7,3]},all.stats))
      
      models.summary<-mapply(function(x){x[,3]},all.stats)
      rownames(models.summary)<-rownames(all.stats[[1]])
      
      all.results<-all.results[[best.model]]
      sensi.estimates<-sensi.estimates[sensi.estimates$model==models[best.model],-8]
      all.stats<-all.stats[[best.model]]
      residuals<-residuals[[best.model]]
    }
    
    to_return<-list("all.results"=all.results,"sensi.estimates"=sensi.estimates,"all.stats"=all.stats,"residuals"=residuals)
    to_return<-to_return[which(names(to_return)%in%output)]
    if(!output.all.models){
      to_return<-c(to_return,list("models.summary"=models.summary))
    }
    return(to_return)
  }
  
  compare.tree_phylm<-function(x,y,xynames,multiphy,subset=NULL,subset.split.name=NULL,tree_phylm_full_results.opt=NULL,output=NULL){
    if(is.null(subset)){
      subset<-xynames
    }
    if(is.null(output)){
      output<-c("all_regs","red_regs","TT","ablines")
    }
    
    comp_multiphy<-suppressWarnings(list(lapply(multiphy,function(x){keep.tip(x,subset)}),lapply(multiphy,function(x){drop.tip(x,subset)})))
    comp_multiphy<-comp_multiphy[!foreach(i=1:2,.combine="c")%do%{any(mapply(is.null,comp_multiphy[[i]]))}]
    n_comp<-c(1:length(comp_multiphy))
    class(comp_multiphy[n_comp])<-"multiPhylo"
    comp_xynames<-list(xynames[xynames%in%subset],xynames[!xynames%in%subset])
    comp_xynames<-comp_xynames[n_comp]
    comp_x<-list(x[xynames%in%subset],x[!xynames%in%subset])
    comp_x<-comp_x[n_comp]
    comp_y<-list(y[xynames%in%subset],y[!xynames%in%subset])
    comp_y<-comp_y[n_comp]
    
    if(is.null(subset.split.name)){
      subset.split.name<-c("x","y")
    }
    subset.split.name<-subset.split.name[n_comp]
    
    out_res<-list()
    
    all_regs<-list()
    for(i in n_comp){
      all_regs[[i]]<-do.call("tree_phylm_full_results",c(list(comp_x[[i]],comp_y[[i]],comp_xynames[[i]],comp_multiphy[[i]]),tree_phylm_full_results.opt))
    }
    all_regs<-setNames(all_regs,subset.split.name)
    if("all_regs"%in%output){
      out_res<-c(out_res,"all_regs"=ifelse(length(all_regs)==1,list(all_regs[[1]]),list(all_regs)))
    }
    
    red_regs<-list()
    for(i in n_comp){
      red_regs[[i]]<-mapply(function(x){x[c(1:7),3]},all_regs[[i]]$all.stats)
      rownames(red_regs[[i]])<-c("intercept","intercept SE","intercept p-value","slope","slope SE","slope p-value","AIC")
    }
    red_regs<-setNames(red_regs,subset.split.name)
    if("red_regs"%in%output){
      out_res<-c(out_res,"red_regs"=ifelse(length(red_regs)==1,list(red_regs[[1]]),list(red_regs)))
    }
    
    if(any(c("TT_data","TT")%in%output)){
      TT_data<-setNames(foreach(i=c(1,4))%do%{
        num.data<-numeric()
        names<-character()
        for(j in n_comp){
          mean<-red_regs[[j]][i,which.min(red_regs[[j]][7,])]
          n<-dim(all_regs[[j]]$residuals[[1]])[2]
          var<-red_regs[[j]][i+1,which.min(red_regs[[j]][7,])]*sqrt(n)
          names<-c(names,subset.split.name[j])
          num.data<-c(num.data,mean,var,n)
        }
        list("num.data"=num.data,"names"=names)
      },c("intercept","slope"))
      if("TT_data"%in%output){
        out_res<-c(out_res,"TT_data"=ifelse(length(TT_data)==1,list(TT_data[[1]]),list(TT_data)))
      }
      if(all("TT"%in%output&max(n_comp)>1)){
        TT<-c(foreach(i=c(1,2))%do%{
          do.call("param.t.test",c(setNames(as.list(TT_data[[i]]$num.data),c("mean_x","var_x","n_x","mean_y","var_y","n_y")),
                                   setNames(as.list(TT_data[[i]]$names),c("name_x","name_y"))))
        },
        list(setNames(foreach(j=n_comp,.combine="c")%do%{red_regs[[j]][3,which.min(red_regs[[j]][7,])]},subset.split.name)),
        list(setNames(foreach(j=n_comp,.combine="c")%do%{red_regs[[j]][6,which.min(red_regs[[j]][7,])]},subset.split.name)))[c(1,3,2,4)]
        names(TT)<-c("intercepts t-test","intercepts p-values","slopes t-test","slopes p-values")
        out_res<-c(out_res,"TT"=list(TT))
      }
    }
    if("ablines"%in%output){
      ablines<-foreach(i=n_comp)%do%{
        setNames(red_regs[[i]][c(1,4),which.min(red_regs[[i]][7,])],c("a","b"))
      }
      ablines<-setNames(ablines,subset.split.name)
      out_res<-c(out_res,"ablines"=ifelse(length(ablines)==1,list(ablines[[1]]),list(ablines)))
    }
    
    out_res
  }
  
  conv.pgls<-function(x,y,phy,names=NULL,tol=1e-6,thinning=0.1,output="all",trace.conv=TRUE,silent=TRUE){
    custom.pgls.confint<-function (ML,fix,lower.b,upper.b,param.CI,x,y,V, which = c("kappa", "lambda", "delta")) {
      which <- match.arg(which)
      whichNum <- which(names(fix) == which)
      opt <- fix[whichNum]
      fix <- fix[-whichNum]
      belowML <- c(lower.b[[which]], opt)
      aboveML <- c(opt, upper.b[[which]])
      MLdelta <- (qchisq(param.CI, 1)/2)
      offset <- (-ML) + MLdelta
      lowerBound.ll <- pgls.likelihood(structure(lower.b[[which]], names = which), 
                                       fix, y, x, V, optim.output = TRUE)
      upperBound.ll <- pgls.likelihood(structure(upper.b[[which]], names = which), 
                                       fix, y, x, V, optim.output = TRUE)
      lrt0 <- 2 * (ML - lowerBound.ll)
      lrt1 <- 2 * (ML - upperBound.ll)
      lowerBound.p <- 1 - pchisq(lrt0, 1)
      upperBound.p <- 1 - pchisq(lrt1, 1)
      ll.fun <- function(opt) {
        pg <- pgls.likelihood(opt, fix, y, x, V, optim.output = TRUE, 
                              names.optim = which)
        ll <- pg + offset
        return(ll)
      }
      lowerCI <- if (lowerBound.ll < (ML - MLdelta)) 
        uniroot(ll.fun, interval = belowML)$root
      else lower.b[[which]]
      upperCI <- if (upperBound.ll < (ML - MLdelta)) 
        uniroot(ll.fun, interval = aboveML)$root
      else upper.b[[which]]
      return(c(lowerCI, upperCI))
    }
    
    if(length(output)==1){
      if(output=="all"){
        output<-c("fits","residuals","vals","AIC_table","estimates.comparison")
      }
    }
    print_fits<-any(output=="fits")
    print_residuals<-any(output=="residuals")
    print_vals<-any(output=="vals")
    print_AIC_table<-any(output=="AIC_table")
    print_estimates.comparison<-any(output=="estimates.comparison")
    
    ny<-length(y[1,])
    if(length(colnames(y))==0|length(colnames(y))<ny){
      colnames(y)<-paste(rep("y",ny),1:ny,sep="")
    }
    if(is.null(names)){ # Without names specified, data is assumed to be sorted according to the phylogeny tip labels
      names<-phy$tip.label
    }
    old_x<-x
    old_y<-y
    
    y<-list()
    for(i in 1:ny){
      data<-comparative.data(phy,data.frame(x=old_x,old_y,names=names),names.col="names",vcv.dim=3)
      formula<-paste0(colnames(old_y)[i],"~x")
      m<-model.frame(formula,data$data)
      y[[i]]<-m[,1]
    }
    x<-model.matrix(as.formula(formula), m)
    V<-VCV.array(data$phy,dim=3)
    
    bounds<-list(kappa = c(1e-06, 3), lambda = c(1e-06, 1), delta = c(1e-06, 3))
    lower.b <- sapply(bounds, "[", 1)
    upper.b <- sapply(bounds, "[", 2)
    optimPar<-rep(list(unlist(lapply(bounds,mean))),ny)
    fixedPar<-lapply(optimPar,function(x){x[-c(1:ny)]})
    bounds.gap<-mapply(function(x){x[2]-x[1]},bounds)
    trace.LB<-matrix(ncol=3,nrow=0)
    trace.UB<-matrix(ncol=3,nrow=0)
    trace.par<-matrix(ncol=3,nrow=0)
    trace.allpar<-matrix(ncol=3*ny,nrow=0)
    ML_par<-NULL
    ML_CI<-NULL
    
    cl<-makeCluster(detectCores()-1)
    registerDoParallel(cl)
    while(any(bounds.gap>tol)){
      local<-foreach(i=1:ny,.packages="caper")%dopar%{optim(par = optimPar[[i]], fn = pgls.likelihood, 
                                                            method = "L-BFGS-B", control = list(fnscale = -1), upper = upper.b, 
                                                            lower = lower.b, V=V, y=y[[i]], x=x, fixedPar = fixedPar[[i]], 
                                                            optim.output = TRUE)}
      local.optimPar<-lapply(local,function(x){x$par})
      if(is.null(ML_par)){ML_par<-local.optimPar}
      if(length(bounds.gap)==3&&all(bounds.gap==mapply(function(x){x[2]-x[1]},bounds))){
        local.ll<-lapply(local,function(x){x$value})
        new.bounds<-foreach(i=1:ny,.packages="foreach")%dopar%{
          setNames(foreach(k=names(optimPar[[i]]))%dopar%{
            custom.pgls.confint(ML=local.ll[[i]],fix=local.optimPar[[i]],lower.b,upper.b,param.CI=0.95,x=x,y=y[[i]],V=V,which=k)
          },names(optimPar[[i]]))
        }
        if(is.null(ML_CI)){ML_CI<-new.bounds}
        new.bounds<-setNames(foreach(i=1:length(optimPar[[1]]))%do%{c(max(mapply(function(x){x[[i]]},new.bounds)[1,]),min(mapply(function(x){x[[i]]},new.bounds)[2,]))},names(optimPar[[i]]))
        new.bounds<-lapply(new.bounds,function(x){c(min(x),max(x))})
        lower.b<-sapply(new.bounds,"[",1)
        upper.b<-sapply(new.bounds,"[",2)
      }
      else{
        old.lower.b<-lower.b
        old.upper.b<-upper.b
        
        # Additional details to what is provided in the Supplementary file 4.
        # The thinning is not equal on both sides of the interval, it is weighted accordingly with the average 
        # parameter value for all regressions: for a 0-1 interval and a 0.9 average parameter value, the interval 
        # will be narrowed of (0.9-0)/(1-0) = 0.9 = 90% on its inferior bound and of (1-0.9)/(1-0) = 0.1 = 10% on 
        # its superior one; with a total 10% of interval narrowing, the new interval will go 
        # from 0+0.9*(1-0)/10 = 0.09 to 1-0.1*(1-0)/10 = 0.99. 
        # A convergence is achieved when the average width of the interval of each parameter is below a user-given 
        # criterion (set by default to 10-4), and the retained value for each parameter is the mean of this interval. 
        
        lower.b<-foreach(i=1:length(old.lower.b),.combine="c")%do%{old.lower.b[i]+thinning*bounds.gap[names(old.lower.b)[i]]*((mean(mapply(function(x){x[names(old.lower.b)[i]]},optimPar))-old.lower.b[i])/(old.upper.b[i]-old.lower.b[i]))}
        upper.b<-foreach(i=1:length(old.upper.b),.combine="c")%do%{old.upper.b[i]-thinning*bounds.gap[names(old.upper.b)[i]]*((old.upper.b[i]-mean(mapply(function(x){x[names(old.upper.b)[i]]},optimPar)))/(old.upper.b[i]-old.lower.b[i]))}
      }
      
      temp_bounds.gap<-c(upper.b-lower.b)
      lower.b<-lower.b[temp_bounds.gap>tol]
      upper.b<-upper.b[temp_bounds.gap>tol]
      
      optimPar <- lapply(local.optimPar,function(x){x[which(temp_bounds.gap>tol)]})
      newfixedPar <- lapply(local.optimPar,function(x){x[which(temp_bounds.gap<tol)]})
      if(length(newfixedPar[[1]])>0){
        temp<-matrix(ncol=ny,nrow=length(newfixedPar[[1]]),mapply(function(x){x},as.list(newfixedPar)))
        rownames(temp)<-names(newfixedPar[[1]])
        newfixedPar<- rep(list(apply(temp,1,mean)),ny)
      }
      fixedPar<-foreach(i=1:ny)%do%{c(newfixedPar[[i]],fixedPar[[i]])[c("kappa","lambda","delta")[which(!c("kappa","lambda","delta")%in%unique(c(na.omit(mapply(names,optimPar)))))]]}
      allPars<-foreach(i=1:ny)%do%{c(optimPar[[i]],fixedPar[[i]])[c("kappa","lambda","delta")]}
      bounds.gap<-temp_bounds.gap
      final_vals<-setNames(foreach(i=1:3,.combine="c")%do%{mean(mapply(function(x){x[i]},allPars))},c("kappa","lambda","delta"))
      sd_vals<-setNames(foreach(i=1:3,.combine="c")%do%{sd(mapply(function(x){x[i]},allPars))},c("kappa","lambda","delta"))
      if(!silent){print(bounds.gap)}
      if(trace.conv){
        all.LB<-if(length(lower.b)==3){lower.b}else{c(setNames(rep(0,3-length(lower.b)),c("kappa","lambda","delta")[which(!c("kappa","lambda","delta")%in%names(lower.b))]),lower.b)[c("kappa","lambda","delta")]}
        all.UB<-if(length(upper.b)==3){upper.b}else{c(setNames(rep(0,3-length(upper.b)),c("kappa","lambda","delta")[which(!c("kappa","lambda","delta")%in%names(upper.b))]),upper.b)[c("kappa","lambda","delta")]}
        trace.LB<-rbind(trace.LB,all.LB)
        trace.UB<-rbind(trace.UB,all.UB)
        trace.par<-rbind(trace.par,final_vals)
        trace.allpar<-rbind(trace.allpar,unlist(allPars))
      }
    }
    if(trace.conv){
      traces<-foreach(i=1:3)%do%{cbind(trace.par[trace.LB[,i]>0,i],trace.LB[trace.LB[,i]>0,i],trace.UB[trace.LB[,i]>0,i])}
      par(pty="s",mfrow=c(1,3))
      for (i in 1:3){
        plot(1,1,type="n",xlab="iterations",ylab=paste0(c("kappa","lambda","delta")[i]," value"),ylim=c((min(trace.LB[trace.LB>0])),(max(trace.UB[trace.UB>0]))),xlim=c(1,max(mapply(function(x){dim(x)[1]},traces))))
        for (j in 2:3){
          points(traces[[i]][,j],type="l",lwd=1)
        }
        text(55,traces[[i]][length(traces[[i]][,1]),1]+0.025,c("kappa","lambda","delta")[i])
        polygon(c(c(1:length(traces[[i]][,2])),c(length(traces[[i]][,2]):1)),c(traces[[i]][,2],rev(traces[[i]][,3])),col=alpha("black",0.25),border=NA)
        for(j in 1:ny){
          points((trace.allpar[-1,(i+3*(j-1))]),type="l",col=contrasting.palette(ny)[j],lwd=2)
        }
        points(traces[[i]][,1],type="l",lwd=2)
      }
      par(pty="m",mfrow=c(1,1))
    }
    stopCluster(cl)
    
    data<-comparative.data(phy,data.frame(x=old_x,old_y,names=names),names.col="names",vcv.dim=3)
    
    if(print_AIC_table|print_estimates.comparison){
      indep_fits<-foreach(i=1:ny)%do%{
        do.call("pgls",c(list(formula(paste0(colnames(old_y)[i],"~x")),data=data),ML_par[[i]]))
      }
    }
    
    if(print_fits|print_residuals|print_AIC_table|print_estimates.comparison){
      final_fits<-foreach(i=1:ny)%do%{
        do.call("pgls",c(list(formula(paste0(colnames(old_y)[i],"~x")),data=data),final_vals))
      }
      if(print_residuals){
        final_residuals<-foreach(i=1:ny)%do%{residuals(final_fits[[i]])}
        final_residuals<-as.data.frame(final_residuals)
        colnames(final_residuals)<-colnames(old_y)
      }
    }
    
    if(print_AIC_table){
      AIC_table<-as.table(foreach(i=1:ny,.combine="cbind")%do%{t(AIC(final_fits[[i]],indep_fits[[i]]))})
      colnames(AIC_table)<-paste(rep(colnames(old_y),each=2),rep(c("opt","ML"),ny))
    }
    
    if(print_estimates.comparison){
      estimates<-as.table(cbind(foreach(i=1:ny,.combine="cbind")%do%{
        foreach(j=1:3,.combine="rbind")%do%{
          i_val<-ML_par[[i]][[j]]
          i_range<-ML_CI[[i]][[j]]
          paste0(round(i_val,4),paste(rep(" ",7-nchar(round(i_val,4))),collapse=""),"(",round(i_range[1],4)," : ",round(i_range[2],4),")")
        }},
        round(final_vals,4)))
      dimnames(estimates)<-list("parameters"=c("kappa","lambda","delta"),"models"=c(paste(colnames(old_y),c("value, 95%IC"),sep=": "),"opt"))
    }
    
    to_output<-list("fits"=if(print_fits){final_fits},
                    "residuals"=if(print_residuals){final_residuals},
                    "vals"=if(print_vals){final_vals},
                    "AIC_table"=if(print_AIC_table){AIC_table},
                    "estimates.comparison"=if(print_estimates.comparison){estimates})
    to_output<-to_output[which(!mapply(is.null,to_output))]
    
    return(to_output)
  }
  
  # Wrapper of the pFDA function to perform it over more than one phylogeny and to output organized results
  pfda_LSchmitz<-function(data,phy,groups,lambda="opt",pfda.opt=NULL,res.reorder=TRUE){
    all_results<-list("confusion"=list(),"test.results"=list(),"training.results"=list(),"DA.scores"=list(),"percent.explained"=list(),"opt.L"=c())
    for (i in 1:length(phy)){
      if(is.list(data)&&length(data)==length(phy)){
        i_data<-data[[i]]
      }
      else{
        i_data<-data
      }
      curr_data<-i_data[match(phy[[i]]$tip.label,rownames(i_data)),]
      curr_groups<-groups[names(groups)%in%rownames(curr_data)]
      curr_groups<-curr_groups[match(rownames(curr_data),names(curr_groups))]
      curr_taxa<-rownames(curr_data)
      
      testtaxa <- rownames(curr_data)[curr_groups=="unknown"] # specifying taxa with unknown group, e.g. fossils
      testtaxan <- row(curr_data)[curr_groups=="unknown",1]
      trainingtaxa <- rownames(curr_data[-testtaxan,]) # creating a dataframe that only contains taxa with known group affiliation
      X <- curr_data[-testtaxan,]
      dd <- curr_data[-testtaxan,]
      g <- curr_groups[-testtaxan]
      tre <- drop.tip(phy[[i]], testtaxa)
      
      if(is.list(lambda)&&length(lambda)==length(phy)){
        i_lambda<-lambda[[i]]
      }
      else{
        i_lambda<-lambda
      }
      
      if(i_lambda=="opt"|i_lambda<0|i_lambda>1){
        ol1 <- optLambda(X,g,tre,plot=FALSE,export=FALSE)$optlambda
        all_results[[6]]<-c(all_results[[6]],ol1)
      }
      else{
        ol1 <- i_lambda
      }
      
      if(!is.null(pfda.opt)){
        if(all(mapply(length,pfda.opt)==length(phy))){
          i_pfda.opt<-list(val_k=pfda.opt$val_k[[i]],val_d=pfda.opt$val_d[[i]],KLD=pfda.opt$KLD[[i]])
        }
        else{
          i_pfda.opt<-pfda.opt
        }
      }
      
      pfda <- do.call("phylo.fda.pred",c(list(curr_data,curr_groups,curr_taxa,phy[[i]],testtaxan,val=as.numeric(ol1)),i_pfda.opt)) # Warning message re: priors can be ignored.
      
      all_results[[1]][[i]]<-pfda$confusion # Misclassified poroportion listed as "error".
      all_results[[5]][[i]]<-c(pfda$percent.explained[1],diff(pfda$percent.explained))
      
      # pFDA performance on tested taxa (i.e., those with unknown ecology)
      test.class <- as.character(predict(pfda, pfda$DATAtest, type="class"))
      test.variates <- predict(pfda, pfda$DATAtest, type="variates")
      test.prob <- predict(pfda, pfda$DATAtest, type="posterior")
      all_results[[2]][[i]] <- cbind(test.class, test.prob, test.variates)
      colnames(all_results[[2]][[i]]) <- c("predicted class", paste0("P(",levels(factor(g)),")"), paste0("DA",col(test.variates)[1,]))
      rownames(all_results[[2]][[i]]) <- testtaxa
      
      #pFDA performance on trained taxa (i.e., those with known ecology)
      training.class <- as.character(predict(pfda, pfda$DATA, type="class"))
      training.variates <- predict(pfda, pfda$DATA, type="variates")
      training.prob <- predict(pfda, pfda$DATA, type="posterior")
      all_results[[3]][[i]] <- cbind(as.character(g), training.class, training.prob, training.variates)
      colnames(all_results[[3]][[i]]) <- c("true class", "predicted class", paste0("P(",levels(factor(g)),")"), paste0("DA",col(training.variates)[1,]))
      rownames(all_results[[3]][[i]]) <- trainingtaxa
      
      # Combining DA scores
      training <- as.data.frame(cbind(training.variates, as.character(curr_groups[-testtaxan])))
      colnames(training) <- c(paste0("DA",col(training.variates)[1,]), "groups")
      rownames(training) <- curr_taxa[-testtaxan]
      unknown <- as.character(rep("unknown", times=length(testtaxan)))
      test <- as.data.frame(cbind(test.variates, unknown))
      colnames(test) <- c(paste0("DA",col(test.variates)[1,]), "groups")
      rownames(test) <- testtaxa
      scatter <- rbind(training, test)
      scatter[, 1:2] <- lapply(scatter[,1:2], as.character)
      scatter[, 1:2] <- lapply(scatter[,1:2], as.numeric)
      scatter<-scatter[match(phy[[1]]$tip.label,rownames(scatter)),]
      all_results[[4]][[i]]<-scatter
    }
    
    all_results[[5]]<-matrix(nrow=length(all_results[[5]]),ncol=length(all_results[[5]][[1]]),unlist(all_results[[5]]),byrow=TRUE)
    
    if(is.null(all_results[[6]])){
      all_results<-all_results[-6]
    }
    if(res.reorder){
      for (i in 2:4){
        all_results[[i]]<-lapply(all_results[[i]],function(x){as.matrix(x[match(rownames(all_results[[i]][[1]]),rownames(x)),])})
      }
    }
    return(all_results)
  }
  
  # Wrapper of the lda function to perform it over more than one phylogeny and to output organized results
  rep_lda<-function(data,groups,test,res.reorder=TRUE){
    all_results<-list("confusion"=list(),"test.results"=list(),"training.results"=list(),"DA.scores"=list(),"percent.explained"=list())
    for (i in 1:length(data)){
      curr_data<-data[[i]][match(names(groups),rownames(data[[i]])),]
      subset<-which(!rownames(curr_data)%in%test)
      curr_groups<-groups
      curr_groups[-subset]<-NA
      curr_groups<-factor(curr_groups)
      curr_lda<-lda(x=curr_data,grouping=curr_groups,subset=subset)
      curr_lda_CV<-lda(x=curr_data,grouping=curr_groups,subset=subset,CV=TRUE)
      VIE_prediction<-predict(curr_lda,curr_data[-subset,])
      colnames(curr_lda_CV$posterior)<-colnames(VIE_prediction$posterior)<-paste0("P(",colnames(VIE_prediction$posterior),")")
      
      all_results[[1]][[i]]<-mda:::confusion.default(curr_lda_CV$class,curr_groups[subset])
      all_results[[2]][[i]]<-cbind("predicted class"=as.character(VIE_prediction$class),
                                   VIE_prediction$posterior,
                                   VIE_prediction$x)
      all_results[[3]][[i]]<-cbind("true class"=as.character(curr_groups[subset]),
                                   "predicted class"=as.character(curr_lda_CV$class),
                                   curr_lda_CV$posterior,
                                   as.matrix(curr_data[subset,])%*%curr_lda$scaling)
      full_LDA_x <- rbind(predict(curr_lda)$x, VIE_prediction$x)
      full_LDA_x<-full_LDA_x[match(names(groups),rownames(full_LDA_x)),]
      full_LDA_x<-as.data.frame(full_LDA_x)
      full_LDA_x$groups<-groups
      all_results[[4]][[i]]<-as.matrix(full_LDA_x)
      all_results[[5]][[i]]<-100*curr_lda$svd/sum(curr_lda$svd)
    }
    
    all_results[[5]]<-matrix(nrow=length(all_results[[5]]),ncol=length(all_results[[5]][[1]]),unlist(all_results[[5]]),byrow=TRUE)
    
    if(res.reorder){
      for (i in 2:4){
        all_results[[i]]<-lapply(all_results[[i]],function(x){as.matrix(x[match(rownames(all_results[[i]][[1]]),rownames(x)),])})
      }
    }
    return(all_results)
  }
  
  # Simple function to get average/majority results from lists of results
  averagize<-function(x){
    if(suppressWarnings(all(is.na(as.numeric(x))))){
      if(length(table(x))>1){
        paste(paste0(names(table(x)),": ",round(100*table(x)/sum(table(x)),2),"%"),collapse="; ")
      }
      else{
        names(table(x))
      }
    }
    else{
      mean(as.numeric(x))
    }
  }
  
  # Function to plot groups (morphospaces) of data and to add the Vielasia point, with several graphical options
  Vielasia_graph<-function(x,y,to_rm=NA,groups,Vielasia,cols,Vielasia_col="red",opt,xlab,ylab,cex,pch=21,new=TRUE,legend=TRUE,leg.pos){
    x_V<-x[Vielasia]
    y_V<-y[Vielasia]
    if(!all(is.na(to_rm))&!all(is.null(to_rm))&all(to_rm!=0)){
      if(!Vielasia%in%to_rm){
        to_rm<-c(to_rm,Vielasia)
      }
      if(length(pch)==length(x)){
        pch<-pch[-to_rm]
      }
      x<-x[-to_rm]
      y<-y[-to_rm]
      if(!is.factor(groups)){
        groups<-factor(groups[-to_rm])
      }
      else{
        groups<-droplevels(groups[-to_rm])
      }
      
    }
    else{
      if(length(pch)==length(x)){
        pch<-pch[-Vielasia]
      }
      x<-x[-Vielasia]
      y<-y[-Vielasia]
    }
    if(!is.factor(groups)){
      groups<-factor(groups[-Vielasia])
    }
    else{
      groups<-droplevels(groups[-Vielasia])
    }
    
    if(opt==0){
      xlim<-range(x)
      ylim<-range(y)
    }
    if(opt==1|opt==2){
      if(opt==1){
        morphospaces<-morphospace(x,y,groups,output=NA)
      }
      if(opt==2){
        morphospaces<-morphospace(x,y,groups,output=NA,smoothing.method="spline")
      }
      xlim<-range(na.omit(mapply(function(x){x[,1]},morphospaces)))
      ylim<-range(na.omit(mapply(function(x){x[,2]},morphospaces)))
    }
    par(pty="s")
    if(new){plot(x,y,xlab=xlab,ylab=ylab,type="n",xlim=xlim,ylim=ylim)}
    if(missing(cex)){
      cex<-rep(1,length(x))
      cex_V<-2
    }
    else{
      cex_V<-cex[Vielasia]
      cex<-cex[-Vielasia]
    }
    if(length(pch)>1){
      if(length(pch)==length(x)){
        pch<-lapply(c(1:nlevels(groups)),function(x){pch[groups==levels(groups)[x]]})
      }
      if(!is.list(pch)){
        pch<-as.list(pch)
      }
    }
    for (i in 1:nlevels(groups)){
      points(x[groups==levels(groups)[i]],y[groups==levels(groups)[i]],pch=if(length(pch)>1){pch[[i]]}else{pch},col=NA,bg=cols[i],cex=cex[groups==levels(groups)[i]])
      if((opt==1|opt==2)&table(groups)[i]>1){
        points(rbind(morphospaces[[which(which(table(groups)>1)==i)]],morphospaces[[which(which(table(groups)>1)==i)]][1,]),type="l",lwd=2,col=cols[i])
      }
    }
    points(x_V,y_V,pch=23,col="black",bg=Vielasia_col,cex=cex_V)
    # arrows(x0=x_V+ifelse(x_V<mean(xlim),1,-1)*0.125*diff(xlim),
    #        y0=y_V+ifelse(y_V<mean(ylim),1,-1)*0.125*diff(ylim),
    #        x1=x_V+ifelse(x_V<mean(xlim),1,-1)*0.025*diff(xlim),
    #        y1=y_V+ifelse(y_V<mean(ylim),1,-1)*0.025*diff(ylim),
    #        length=0.1,col=Vielasia_col,lwd = 5)
    if(legend){
      if(missing(leg.pos)){
        leg.pos<-c("topleft","topright","bottomleft","bottomright")[which.min(c(
          length(which(x<mean(xlim)&y>mean(ylim))),
          length(which(x>mean(xlim)&y>mean(ylim))),
          length(which(x<mean(xlim)&y<mean(ylim))),
          length(which(x>mean(xlim)&y<mean(ylim)))))]
      }
      legend(leg.pos[1],lwd=2,col=cols,legend=levels(groups),bty="n")
    }
    par(pty="m")
  }
  
  # Compute the average position and the 'full' and '95%' morphospaces (see Data S1 for explanations) for given data.
  plots.variation<-function(list_data=NA,chull_data=NA,avg_data,data.name="analysis",axes=c(1,2),axes.contrib=c(50,50),axes.sign=c(1,1),morpho.groups,morpho.cols=NA,morpho.pch=21,VIE.col="red",out=FALSE,out.names=NA,single.out=1,return.chull.data=FALSE){
    get_chulls<-function(list_data,axes){
      cl<-makeCluster(detectCores()-1)
      registerDoParallel(cl)
      chull_full<-foreach(i=1:length(list_data[[1]][,1]),.packages="ULT")%dopar%{
        x_data<-as.numeric(mapply(function(x){x[i,axes[1]]},list_data))
        y_data<-as.numeric(mapply(function(x){x[i,axes[2]]},list_data))
        morphospace(x_data,y_data,output=NA)[[1]]
      }
      chull_95<-foreach(i=1:length(list_data[[1]][,1]),.packages=c("ULT","sf"))%dopar%{
        x_data<-as.numeric(mapply(function(x){x[i,axes[1]]},list_data))
        y_data<-as.numeric(mapply(function(x){x[i,axes[2]]},list_data))
        centroid<-st_centroid(st_polygon(list(as.matrix(rbind(chull_full[[i]],chull_full[[i]][1,])))))
        x_center<-centroid[1]
        y_center<-centroid[2]
        each_square_dist<-sqrt((x_data-x_center)^2+(y_data-y_center)^2)
        kept_ones<-which(each_square_dist<=quantile(each_square_dist,0.95))
        x_data_95<-x_data[kept_ones]
        y_data_95<-y_data[kept_ones]
        morphospace(x_data_95,y_data_95,output=NA)[[1]]
      }
      stopCluster(cl)
      chull_data<-list("full"=chull_full,"95"=chull_95)
      return(chull_data)
    }
    
    plot_chulls<-function(chull_data,avg_data,data.name,axes,axes.contrib,axes.sign,morpho.groups,morpho.cols,morpho.pch,VIE.col,out,out.names,single.out){
      xlab<-paste0(data.name," ",axes[1]," (",axes.contrib[1],"%)")
      ylab<-paste0(data.name," ",axes[2]," (",axes.contrib[2],"%)")
      
      chull_data<-lapply(chull_data,function(x){lapply(x,function(y){t(t(as.matrix(y))*axes.sign)})})
      avg_data<-t(t(as.matrix(avg_data))*axes.sign)
      
      par(mar=rep(3,4),mgp=c(2,0.5,0),pty="s")
      
      # Plotting "full" and "95%" morphospaces for each individual
      for(c in 1:2){
        plot(x=0,y=0,type="n",xlab=xlab,ylab=ylab,xlim=range(unlist(mapply(function(x){x[,1]},chull_data[[c]]))),ylim=range(unlist(mapply(function(x){x[,2]},chull_data[[c]]))))
        for (i in 1:length(chull_data[[c]])){
          points(c(chull_data[[c]][[i]][,1],chull_data[[c]][[i]][1,1]),c(chull_data[[c]][[i]][,2],chull_data[[c]][[i]][1,2]),type="l",lwd=1,col=c(morpho.cols,VIE.col)[morpho.groups][i])
        }
        for (i in 1:(nlevels(morpho.groups)-1)){
          morphospace(do.call("rbind",chull_data[[c]][morpho.groups==levels(morpho.groups)[i]]),col=c(morpho.cols)[i],lwd=2,plot.points = FALSE)
        }
        morphospace(do.call("rbind",chull_data[[c]][which(names(morpho.groups)=="Vielasia_sigei")]),col=c(VIE.col),lwd=2,plot.points = FALSE)
      }
      
      # Plotting average position of the species
      
      Vg.args<-list(x=as.numeric(avg_data[,1]),
                    y=as.numeric(avg_data[,2]),
                    Vielasia=which(names(morpho.groups)=="Vielasia_sigei"),
                    groups = morpho.groups,
                    cols=c(morpho.cols),
                    opt=1,
                    xlab=xlab,
                    ylab=ylab,
                    pch=morpho.pch,
                    legend=FALSE,
                    Vielasia_col = VIE.col)
      
      if(single.out){
        postscript(out.names[1])
        do.call("Vielasia_graph",Vg.args)
        dev.off()
      }
      do.call("Vielasia_graph",Vg.args)
      if(out){
        dev.off()
      }
    }
    
    type<-ifelse(length(data.name)>1|length(axes)>2,2,1)
    
    if(type==1){
      if(any(is.na(chull_data))){
        chull_data<-get_chulls(list_data,axes)
      }
      
      single.out<-ifelse(out,TRUE,FALSE)
      if(out&&length(out.names)==2){
        postscript(out.names[2])
      }
      par(mfrow=c(1,3))
      plot_chulls(chull_data,avg_data,data.name,axes,axes.contrib,morpho.groups,morpho.cols,morpho.pch,VIE.col,out,out.names,single.out)
    }
    
    if(type==2){
      if(length(axes)>2){
        if(length(axes)==3){
          if(any(is.na(chull_data))){
            chull_data_12<-get_chulls(list_data,axes[c(1,2)])
            chull_data_13<-get_chulls(list_data,axes[c(1,3)])
            chull_data<-setNames(list(chull_data_12,chull_data_13),c("axes 1 and 2","axes 1 and 3"))
          }
          
          avg_data<-c(list(avg_data[,axes[c(1,2)]]),list(avg_data[,axes[c(1,3)]]))
          axes<-c(list(axes[c(1,2)]),list(axes[c(1,3)]))
          axes.contrib<-c(list(axes.contrib[c(1,2)]),list(axes.contrib[c(1,3)]))
          if(length(axes.sign)==2){
            axes.sign<-rep(list(axes.sign),2)
          }
          else{
            axes.sign<-c(list(axes.sign[c(1,2)]),list(axes.sign[c(1,3)]))
          }
        }
        if(length(axes)==4){
          if(any(is.na(chull_data))){
            chull_data_12<-get_chulls(list_data,axes[c(1,2)])
            chull_data_34<-get_chulls(list_data,axes[c(3,4)])
            chull_data<-setNames(list(chull_data_12,chull_data_34),c("axes 1 and 2","axes 3 and 4"))
          }
          
          avg_data<-c(list(avg_data[,axes[c(1,2)]]),list(avg_data[,axes[c(3,4)]]))
          axes<-c(list(axes[c(1,2)]),list(axes[c(3,4)]))
          axes.contrib<-c(list(axes.contrib[c(1,2)]),list(axes.contrib[c(3,4)]))
          if(length(axes.sign)==2){
            axes.sign<-rep(list(axes.sign),2)
          }
          else{
            axes.sign<-c(list(axes.sign[c(1,2)]),list(axes.sign[c(3,4)]))
          }
        }
        data.name<-rep(data.name,2)
        morpho.groups<-c(rep(list(morpho.groups),2))
        morpho.pch<-c(rep(list(morpho.pch),2))
      }
      else if(length(data.name)>1){
        if(any(is.na(chull_data))){
          chull_data_1<-get_chulls(list_data[[1]],axes)
          chull_data_2<-get_chulls(list_data[[2]],axes)
          chull_data<-setNames(list(chull_data_1,chull_data_2),data.name)
        }
        else{
          names(chull_data)<-data.name
        }
        
        avg_data<-lapply(avg_data,function(x){x[,1:2]})
        axes<-c(list(axes),list(axes))
        axes.contrib<-lapply(axes.contrib,function(x){x[1:2]})
        if(length(axes.sign)==2){
          axes.sign<-rep(list(axes.sign),2)
        }
        else{
          axes.sign<-c(list(axes.sign[c(1,2)]),list(axes.sign[c(3,4)]))
        }
      }
      
      if(out&&length(out.names)==2){
        postscript(out.names[2])
      }
      par(mfrow=c(2,3))
      for(d in 1:2){
        plot_chulls(chull_data[[d]],avg_data[[d]],data.name[d],axes[[d]],axes.contrib[[d]],axes.sign[[d]],morpho.groups[[d]],morpho.cols,morpho.pch[[d]],VIE.col,out=ifelse(d==2,TRUE,FALSE),out.names,single.out=ifelse(d==single.out,TRUE,FALSE))
      }
    }
    
    if(return.chull.data){
      return(chull_data)
    }
  }
}

# Functions for phylogenetic flexible discriminant analysis (pFDA) after Schmitz & Motani (2011), retrieved from https://github.com/lschmitz/phylo.fda
{
  ## phylo.FDA.v0.2
  
  require(nnet)
  require(mda)
  require(ape)
  require(geiger)
  require(lattice)
  
  ###----------------------------------------------------------------------
  ### Internal function from the package mda
  ###----------------------------------------------------------------------
  "contr.fda" <-
    function (p = rep(1, d[1]), contrast.default = contr.helmert(length(p)))
    {
      d <- dim(contrast.default)
      sqp <- sqrt(p/sum(p))
      x <- cbind(1, contrast.default) * outer(sqp, rep(1, d[2] +
                                                         1))
      qx <- qr(x)
      J <- qx$rank
      qr.qy(qx, diag(d[1])[, seq(2, J)])/outer(sqp, rep(1, J -
                                                          1))
    }
  
  ###----------------------------------------------------------------------
  ### Associated functions modified from the package mda
  ###----------------------------------------------------------------------
  "predict.phylo.fda" <-
    function (object, newdata, type = c("class", "variates", "posterior",
                                        "hierarchical", "distances"), prior, dimension = J - 1, ...)
    {
      dist <- function(x, mean, m = ncol(mean)) (scale(x, mean,
                                                       FALSE)^2) %*% rep(1, m)
      type <- match.arg(type)
      means <- object$means
      Jk <- dim(means)
      J <- Jk[1]
      k <- Jk[2]
      if (type == "hierarchical") {
        if (missing(dimension))
          dimension.set <- seq(k)
        else {
          dimension.set <- dimension[dimension <= k]
          if (!length(dimension.set))
            dimension.set <- k
          dimension <- max(dimension.set)
        }
      }
      
      else dimension <- min(max(dimension), k)
      if (missing(newdata))
        y <- predict(object$fit)
      else {
        if (inherits(newdata, "data.frame") || is.list(newdata)) {
          Terms <- delete.response(terms(object))
          attr(Terms, "intercept") <- 0
          newdata <- model.matrix(Terms, newdata)
        }
        y <- predict(object$fit, newdata)
      }
      y <- y %*% object$theta[, seq(dimension), drop = FALSE]
      lambda <- object$values
      alpha <- sqrt(lambda[seq(dimension)])
      sqima <- sqrt(1 - lambda[seq(dimension)])
      newdata <- scale(y, FALSE, sqima * alpha)
      if (missing(prior))
        prior <- object$prior
      else {
        if (any(prior < 0) | round(sum(prior), 5) != 1)
          stop("innappropriate prior")
      }
      means <- means[, seq(dimension), drop = FALSE]
      switch(type, variates = return(newdata), class = {
        n <- nrow(newdata)
        prior <- 2 * log(prior)
        mindist <- dist(newdata, means[1, ], dimension) - prior[1]
        pclass <- rep(1, n)
        for (i in seq(2, J)) {
          ndist <- dist(newdata, means[i, ], dimension) - prior[i]
          l <- ndist < mindist
          pclass[l] <- i
          mindist[l] <- ndist[l]
        }
        ## 2001-10-27: Need to provide levels or else if we get an error
        ## if the predicted classes do no contain all possible classes.
        ## Reported by Greg Jefferis <jefferis@stanford.edu>, fix by
        ## Bj/orn-Helge Mevik <bjorn-helge.mevik@matforsk.no>.
        return(factor(pclass, levels = seq(J),
                      labels = dimnames(means)[[1]]))
      }, posterior = {
        pclass <- matrix(0, nrow(newdata), J)
        for (i in seq(J)) pclass[, i] <- exp(-0.5 * dist(newdata, means[i,
        ], dimension)) * prior[i]
        dimnames(pclass) <- list(dimnames(newdata)[[1]], dimnames(means)[[1]])
        return(pclass/drop(pclass %*% rep(1, J)))
      }, hierarchical = {
        prior <- 2 * log(prior)
        Pclass <- vector("list", length(dimension.set))
        names(Pclass) <- paste("D", dimension.set, sep = "")
        for (ad in seq(along = dimension.set)) {
          d <- dimension.set[ad]
          dd <- seq(d)
          
          mindist <- dist(newdata[, dd, drop = FALSE], means[1, dd, drop = FALSE],
                          d) - prior[1]
          pclass <- rep(1, nrow(newdata))
          for (i in seq(2, J)) {
            ndist <- dist(newdata[, dd, drop = FALSE], means[i, dd,
                                                             drop = FALSE], d) - prior[i]
            l <- ndist < mindist
            pclass[l] <- i
            mindist[l] <- ndist[l]
          }
          levels(pclass) <- dimnames(means)[[1]]
          Pclass[[ad]] <- pclass
        }
        rownames <- dimnames(newdata)[[1]]
        if (is.null(rownames))
          rownames <- paste(seq(nrow(newdata)))
        return(structure(Pclass, class = "data.frame", row.names = rownames,
                         dimensions = dimension.set))
      }, distances = {
        dclass <- matrix(0, nrow(newdata), J)
        for (i in seq(J)) dclass[, i] <- dist(newdata, means[i, ],
                                              dimension)
        dimnames(dclass) <- list(dimnames(newdata)[[1]], dimnames(means)[[1]])
        return(dclass)
      })
    }
  "predict.polyreg.modified" <-
    function (object, newdata, ...)
    {
      if (missing(newdata)) {
        z <- fitted(object)
        if (is.null(z))
          stop("need to supply newdata")
        else return(z)
      }
      degree <- object$degree
      monomial <- object$monomial
      newdata %*% object$coef
    }
  "polyreg.modified" <-
    function (x, y, w, degree = 1, monomial = FALSE, ...)
    {
      #x <- polybasis(x, degree, monomial)
      y <- as.matrix(y) # just making sure ...
      if (iswt <- !missing(w)) {
        if (any(w <= 0))
          stop("only positive weights")
        w <- sqrt(w)
        y <- y * w
        x <- x * w
      }
      qrx <- qr(x)
      coef <- as.matrix(qr.coef(qrx, y))
      
      fitted <- qr.fitted(qrx, y)
      if ((df <- qrx$rank) < ncol(x))
        coef[qrx$pivot, ] <- coef
      if (iswt)
        fitted <- fitted/w
      structure(list(fitted.values = fitted, coefficients = coef,
                     degree = degree, monomial = monomial, df = df), class = "polyreg.modified")
    }
  "print.phylo.fda" <-
    function (x, ...)
    {
      if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
      }
      cat("\nDimension:", format(x$dimension), "\n")
      cat("\nPercent Between-Group Variance Explained:\n")
      print(round(x$percent, 2))
      error <- x$confusion
      df <- x$fit
      if (!is.null(df))
        df <- df$df
      if (!is.null(df)) {
        cat("\nDegrees of Freedom (per dimension):", format(sum(df)),
            "\n")
      }
      if (!is.null(error)) {
        n <- as.integer(sum(error))
        error <- format(round(attr(error, "error"), 5))
        cat("\nTraining Misclassification Error:", error, "( N =",
            n, ")\n")
      }
      invisible(x)
    }
  "plot.phylo.fda" <- function(pfdamodel,gfactor=pfdamodel$g,prdfactor=pfdamodel$prd)
  {
    pfdavar <- predict(pfdamodel, type="variate")
    lim1x <- c(min(pfdavar[,1]),max(pfdavar[,1]))
    lim1y <- c(min(pfdavar[,2]),max(pfdavar[,2]))
    m1 <- 4;m2 <- 1
    oldpar<-par(no.readonly=FALSE);on.exit(par(oldpar));x11(height=8,width=14);par(mfrow=c(1,2),mar=c(m1,m1,m1,m2),oma=c(m2,m2,m2,m2));
    matplot(pfdavar[gfactor==levels(gfactor)[1],1], pfdavar[gfactor==levels(gfactor)[1],2], xlab="pFDA1",ylab="pFDA2", xlim=lim1x, ylim=lim1y, pch=1, col=1, main="True Classes",sub=paste("lambda = ",pfdamodel$val," intrcpt=",pfdamodel$intercept," eqprior=",pfdamodel$eqprior,sep=""))
    for (i in 2:nlevels(gfactor)) matplot(pfdavar[gfactor==levels(gfactor)[i],1], pfdavar[gfactor==levels(gfactor)[i],2], add=TRUE, pch=i, col=i)
    legend(min(lim1x),max(lim1y),levels(gfactor), pch=1:nlevels(gfactor), col=1:nlevels(gfactor))
    
    legend(min(lim1x),min(lim1y)+(max(lim1y)-min(lim1y))*0.1,paste("lambda = ",pfdamodel$val," intrcpt=",pfdamodel$intercept," eqprior=",pfdamodel$eqprior," ",sep=""))
    addEllipseGrp(pfdavar[,1],pfdavar[,2],gfactor, pval=0.95, num=30)
    matplot(pfdavar[prdfactor==levels(prdfactor)[1],1], pfdavar[prdfactor==levels(prdfactor)[1],2], xlab="pFDA1",ylab="pFDA2", xlim=lim1x, ylim=lim1y, pch=1, col=1, main="Predicted Classes",sub=paste("lambda = ",pfdamodel$val," intercept=",pfdamodel$intercept," eqprior=",pfdamodel$eqprior,sep=""))
    for (i in 2:nlevels(prdfactor)) matplot(pfdavar[prdfactor==levels(prdfactor)[i],1], pfdavar[prdfactor==levels(prdfactor)[i],2], add=TRUE, pch=i, col=i)
    legend(min(lim1x),max(lim1y),levels(prdfactor), pch=1:nlevels(prdfactor), col=1:nlevels(prdfactor))
    legend(min(lim1x),min(lim1y)+(max(lim1y)-min(lim1y))*0.1,paste(levels(prdfactor),"=",pfdamodel$prior," ",sep=""))
    legend(max(lim1x)-(max(lim1x)-min(lim1x))*0.2,max(lim1y),signif(attr(pfdamodel$confusion,"error"),4))
    invisible()
  }
  
  ###----------------------------------------------------------------------
  ### Main pFDA function with training data only
  ###----------------------------------------------------------------------
  "phylo.fda" <-function (data,grp,tretre,val=1,
                          dimension = J - 1, eps = .Machine$double.eps,
                          keep.fitted = (n * dimension < 1000), method=polyreg.modified,intercept=TRUE,eqprior=FALSE,priin=1)
  {
    this.call <- match.call()
    if(intercept) data <- cbind(Intercept=rep(1,nrow(data)),data)
    data <- as.matrix(data)
    tretre <- geiger:::rescale.phylo(tretre,"lambda", val)
    g <- as.factor(grp)
    ng <- nlevels(g)
    W <- vcv.phylo(tretre)
    invW<-solve(W)
    invW.eig <- eigen(invW)
    N <- invW.eig$vectors %*% diag(sqrt(invW.eig$values)) %*% solve(invW.eig$vectors)
    divnum <-det(N)^(1/nrow(N))
    N <- N/divnum
    DATA <- N%*%data #Rao (4,57); transforming the data to linear
    n <- nrow(DATA)
    y <- matrix(0,nrow(data),ng)
    for (i in 1:nrow(data)){y[i,g[i]] <- 1}
    Y <- N%*%y #Dummy matrix with phylo bias removed
    x <- DATA
    fg <- factor(g)
    prior <- colSums(Y)/sum(colSums(Y))
    if(eqprior) prior <- c(rep(1/ng,ng))
    if(priin != 1) prior<-priin
    cnames <- levels(fg)
    g <- as.numeric(fg)
    J <- length(cnames)
    weights <- rep(1, n)
    dp <- tapply(weights, g, sum)/n
    theta <- contr.helmert(J)
    theta <- contr.fda(dp, theta)
    
    Theta <- Y%*%theta #fda p.7, above eq2
    fit <- method(x, Theta, weights)
    rss <- t(Theta-fit$fitted) %*% (Theta-fit$fitted)
    ssm <- t(Theta) %*% fitted(fit)/n
    ed <- svd(ssm, nu = 0)
    thetan <- ed$v
    lambda <- ed$d
    lambda[lambda > 1 - eps] <- 1 - eps
    discr.eigen <- lambda/(1 - lambda)
    pe <- (100 * cumsum(discr.eigen))/sum(discr.eigen)
    dimension <- min(dimension, sum(lambda > eps))
    if (dimension == 0) {
      warning("degenerate problem; no discrimination")
      return(structure(list(dimension = 0, fit = fit, call = this.call),
                       class = "phylo.fda"))
    }
    thetan <- thetan[, seq(dimension), drop = FALSE]
    pe <- pe[seq(dimension)]
    alpha <- sqrt(lambda[seq(dimension)])
    sqima <- sqrt(1 - lambda[seq(dimension)])
    vnames <- paste("v", seq(dimension), sep = "")
    means <- scale(theta %*% thetan, FALSE, sqima/alpha)
    dimnames(means) <- list(cnames, vnames)
    names(lambda) <- c(vnames, rep("", length(lambda) - dimension))
    names(pe) <- vnames
    frml <- "grp~"
    nc <- ncol(data)
    varnam <- colnames(data)
    for(i in 1:(nc-1)) frml <- paste(frml,varnam[i],"+", sep="")
    frml <- paste(frml,varnam[nc], sep="")
    frml <- as.formula(frml)
    dset <- as.data.frame(cbind(grp,DATA))
    Terms <- as.call(fda(formula = frml, data = dset, weights = weights))
    obj <- structure(list(percent.explained = pe, values = lambda,
                          means = means, theta.mod = thetan, dimension = dimension,
                          prior = prior, fit = fit, call = this.call, terms = Terms),
                     class = "phylo.fda")
    obj$confusion <- confusion(predict(obj), fg)
    obj$prd <- predict(obj)
    obj$g <- as.factor(grp)
    obj$val <- val
    obj$rss <- sum(diag(rss))
    obj$intercept <- intercept
    obj$eqprior <- eqprior
    if (!keep.fitted)
      obj$fit$fitted.values <- NULL
    obj
  }
  
  ###----------------------------------------------------------------------
  ### Main pFDA function with training and test data
  ###----------------------------------------------------------------------
  
  "phylo.fda.pred" <-function (dataA,grpA,taxtaxA,tretreA,testlistn,val=1,
                               method=polyreg.modified,
                               eps = .Machine$double.eps, intercept=TRUE,eqprior=FALSE,priin=1,KLD=FALSE,val_k=NULL,val_d=NULL)
  {
    ## Preparing data
    this.call <- match.call()
    if(intercept) dataA <- cbind(Intercept=rep(1,nrow(dataA)),dataA)
    dataA <- as.data.frame(dataA)
    nA <- nrow(dataA)
    testlist <- taxtaxA[testlistn]
    traininglist <- taxtaxA[-testlistn]
    rownames(dataA) <- taxtaxA
    tretre <- drop.tip(tretreA,testlistn)
    grp <- grpA[-testlistn]
    grp <- grp[grp %in% names(table(grp))[table(grp) > 0], drop=TRUE]
    g <- as.factor(grp)
    ng <- nlevels(g)
    grpA <- as.factor(grpA)
    ntest <- length(testlist)
    dataA <- as.matrix(dataA)
    if(KLD){ # JM modification : to be able to account for an estimation of kappa, delta, and lambda at the same time, or only lambda (as originally)
      tretreA<-geiger:::rescale.phylo(tretreA,"kappa",val_k) # Kappa transformation
      
      tretreA<-geiger:::rescale.phylo(tretreA,"lambda",val) # Lambda transformation
      
      # Same than geiger:::heights.phylo, but using phytools::nodeHeights, since heights.phylo leads to errors in solving the vcv phylo matrix
      ht<-((max(nodeHeights(tretreA))-nodeHeights(tretreA))[order(tretreA$edge[,2]),])
      ht<-rbind(ht[1:Ntip(tretreA),],rep(max(nodeHeights(tretreA)),2),ht[(Ntip(tretreA)+1):(Ntip(tretreA)+Nnode(tretreA)-1),])
      ht<-as.data.frame(ht)
      colnames(ht)<-c("start","end")
      
      N = Ntip(tretreA)
      Tmax = ht$start[N + 1]
      ht$t = Tmax - ht$end
      ht$e = ht$start - ht$end
      ht$a = ht$t - ht$e
      bl = (ht$a + ht$e)^val_d - ht$a^val_d
      tretreA$edge.length = bl[tretreA$edge[, 2]]
    }
    else{tretreA <- geiger:::rescale.phylo(tretreA,"lambda", val)}
    W <- vcv.phylo(tretreA)
    invW<-solve(W)
    invW.eig <- eigen(invW)
    N <- invW.eig$vectors %*% diag(sqrt(invW.eig$values)) %*% solve(invW.eig$vectors)
    divnum <-det(N)^(1/nrow(N))
    N <- N/divnum
    invN <- solve(N)
    y <- matrix(0,nA,nlevels(grpA))
    for (i in 1:nA){y[i,grpA[i]] <- 1}
    Y <- N%*%y #Dummy matrix with phylo bias removed
    Y <- Y[-testlistn,1:ng]
    DATAA <- N%*%as.matrix(dataA) #Rao (4,57); transforming the data to linear
    DATA <- DATAA[-testlistn,]
    DATAtest <- DATAA[testlistn,]
    n<-nrow(DATA)
    m<-nrow(DATAtest)
    x <- DATA
    fg <- factor(g)
    prior <- colSums(Y)/sum(colSums(Y))
    if(eqprior) prior <- c(rep(1/ng,ng))
    if(priin != 1) prior<-priin
    cnames <- levels(fg)
    g <- as.numeric(fg)
    J <- length(cnames)
    dimension = J - 1
    keep.fitted = (n * dimension < 1000)
    weights <- rep(1, n)
    dp <- tapply(weights, g, sum)/n
    theta <- contr.helmert(J)
    theta <- contr.fda(dp, theta)
    Theta <- Y%*%theta #fda p.7, above eq2
    
    fit <- method(x, Theta, weights)
    rss <- t(Theta-fit$fitted) %*% (Theta-fit$fitted)
    ssm <- t(Theta) %*% fitted(fit)/n
    ed <- svd(ssm, nu = 0)
    thetan <- ed$v
    lambda <- ed$d
    lambda[lambda > 1 - eps] <- 1 - eps
    discr.eigen <- lambda/(1 - lambda)
    pe <- (100 * cumsum(discr.eigen))/sum(discr.eigen)
    dimension <- min(dimension, sum(lambda > eps))
    if (dimension == 0) {
      warning("degenerate problem; no discrimination")
      return(structure(list(dimension = 0, fit = fit, call = this.call),
                       class = "fda"))
    }
    thetan <- thetan[, seq(dimension), drop = FALSE]
    pe <- pe[seq(dimension)]
    alpha <- sqrt(lambda[seq(dimension)])
    sqima <- sqrt(1 - lambda[seq(dimension)])
    vnames <- paste("v", seq(dimension), sep = "")
    means <- scale(theta %*% thetan, FALSE, sqima/alpha)
    dimnames(means) <- list(cnames, vnames)
    names(lambda) <- c(vnames, rep("", length(lambda) - dimension))
    names(pe) <- vnames
    frml <- "grp~"
    nc <- ncol(dataA)
    varnam <- colnames(dataA)
    for(i in 1:(nc-1)) frml <- paste(frml,varnam[i],"+", sep="")
    frml <- paste(frml,varnam[nc], sep="")
    frml <- as.formula(frml)
    dset <- as.data.frame(cbind(grp,DATA))
    Terms <- as.call(fda(formula = frml, data = dset, weights = weights))
    obj <- structure(list(percent.explained = pe, values = lambda,
                          means = means, theta.mod = thetan, dimension = dimension,
                          prior = prior, fit = fit, call = this.call, terms = Terms),
                     class = "phylo.fda")
    obj$confusion <- confusion(predict(obj), fg)
    obj$prd <- predict(obj)
    obj$x<-x
    obj$g <- as.factor(grp)
    obj$val <- val
    obj$rss <- sum(diag(rss))
    obj$intercept <- intercept
    obj$eqprior <- eqprior
    obj$DATAtest <- DATAtest
    obj$DATA <- DATA
    tpred <- predict(obj,DATAtest)
    tpredn <- as.numeric(tpred)
    tpred <- as.matrix(tpred)
    rownames(tpred) <- testlist
    obj$testprediction <- tpred
    obj$testprediction_numeral <- tpredn
    if (!keep.fitted)
      
      obj$fit$fitted.values <- NULL
    obj
  }
  
  
  ###----------------------------------------------------------------------
  ### Function for optimal lambda value search
  ###----------------------------------------------------------------------
  
  "phylo.RSS"<-function (datain,grp,tretre,val=1)
  {
    datainO <- as.matrix(datain)
    datainI <- cbind(Intercept=rep(1,nrow(datainO)),datainO)
    tretre <- geiger:::rescale.phylo(tretre,"lambda", val)
    n <- nrow(datain)
    g <- as.factor(grp)
    ng <- nlevels(g)
    W <- vcv.phylo(tretre)
    invW<-solve(W)
    y <- matrix(0,n,ng) #Dummy matrix without phylo bias
    for (i in 1:n){y[i,g[i]] <- 1}
    invW.eig <- eigen(invW)
    N <- invW.eig$vectors %*% diag(sqrt(invW.eig$values)) %*% solve(invW.eig$vectors)
    Y <- N %*% y # Pretending that there is no phylogenetic bias in y; otherwise Y <- N%*%y
    DATAI <- N%*%datainI
    bhatI <- solve(t(datainI)%*%invW%*%datainI)%*%t(datainI)%*%invW%*%y #Rohlf (9) -- data biased still Rao (4,64)
    yhatI <- datainI%*%bhatI #Rohlf (11)
    RSSyI <- t(y-yhatI) %*% invW %*% (y-yhatI) #Martins and Hansen 1997 (9)
    l0I<- lm(Y~DATAI-1)
    list(RSS=sum(diag(RSSyI))) 
  }
  
  ##dataA=XA;grpA=gA;taxtaxA=taxaA;tretreA=treA;testlistn=testtaxan;val=0;treetrans=lambdaTree
  
  "phylo.RSS.pred" <-function (dataA,grpA,taxtaxA,tretreA,testlistn,val=1)
  {
    dataA <- as.data.frame(dataA)
    nA <- nrow(dataA)
    testlist <- taxtaxA[testlistn]
    traininglist <- taxtaxA[-testlistn]
    rownames(dataA) <- taxtaxA
    tretre <- drop.tip(tretreA,testlistn)
    grp <- grpA[-testlistn]
    grp <- grp[grp %in% names(table(grp))[table(grp) > 0], drop=TRUE]
    g <- as.factor(grp)
    ng <- nlevels(g)
    grpA <- as.factor(grpA)
    icptA <- rep(1,nA)
    dataA <- cbind(icptA,dataA)
    ntest <- length(testlist)
    tretreA <- geiger:::rescale.phylo(tretreA,"lambda", val)
    
    W <- vcv.phylo(tretreA)
    invW<-solve(W)
    invW.eig <- eigen(invW)
    N <- invW.eig$vectors %*% diag(sqrt(invW.eig$values)) %*% solve(invW.eig$vectors)
    invN <- solve(N)
    y <- matrix(0,nA,nlevels(grpA))
    for (i in 1:nA){y[i,grpA[i]] <- 1}
    Y <- N%*%y #Dummy matrix with phylo bias removed
    Y <- Y[-testlistn,1:ng]
    DATAA <- N%*%as.matrix(dataA) #Rao (4,57); transforming the data to linear
    DATA <- DATAA[-testlistn,]
    BHAT <- solve(t(DATA)%*%DATA)%*%t(DATA)%*%Y
    YHAT <- DATA%*%BHAT
    l0<- lm(Y~DATA-1)
    RSSY <- t(Y-YHAT) %*% (Y-YHAT)
    list(RSS=sum(diag(RSSY))) 
  }
  
  ##measurements=X;grps=g;mytree=tre;idc=filename_stem
  
  "optLambda" <- function(measurements,grps,mytree,idc="default",sstep=0.01,srange=c(0,1),fldr="./",plot=TRUE,export=TRUE)
  {
    lambdalist <- seq(min(srange),max(srange),sstep)
    segnum <- length(lambdalist)
    rslt<-matrix(,segnum,2)
    colnames(rslt) <- c("Lambda","RSS") 
    
    for(i in 1:segnum){
      lambdaval <- lambdalist[i]
      rss <- phylo.RSS(measurements,grps,mytree,val=lambdaval)
      rslt[i,] <- c(lambdaval,rss$RSS) 
    }
    
    optlambda <- matrix(,1,1);colnames(optlambda)<- "RSS"
    optlambda[1,1]<-max(rslt[which(rslt[,2]==min(rslt[,2])),1])
    
    if(plot){
      matplot(rslt[,1],rslt[,2],type="l",xlab=expression(lambda),ylab="RSS",main="RSS",lty=1,col=1)
      abline(v=optlambda[1,1],col=2,lty=2);mtext(paste("Optimal Lambda = ",optlambda[1,1],sep=""))
    }
    
    if(export){
      pdf(height=11,width=6,file=paste(fldr,idc,".optLambda.pdf",sep=''));layout(matrix(c(1,2),2,1))
      matplot(rslt[,1],rslt[,2],type="l",xlab=expression(lambda),ylab="RSS",main="RSS",lty=1,col=1)
      abline(v=optlambda[1,1],col=2,lty=2);mtext(paste("Optimal Lambda = ",optlambda[1,1],sep=""))
      dev.off()
    }
    
    list(optlambda=optlambda,rslt=rslt)
  }
  
  "optLambda.pred" <- function(measurementsA,grpsA,taxaA,mytreeA,testn,idc="default",sstep=0.01,srange=c(0,1),fldr="./",plot=TRUE,export=TRUE)
    
  {
    lambdalist <- seq(min(srange),max(srange),sstep)
    segnum <- length(lambdalist)
    rslt<-matrix(,segnum,2)
    colnames(rslt) <- c("Lambda","RSS")
    
    for(i in 1:segnum){
      lambdaval <- lambdalist[i]
      rss <- phylo.RSS.pred(measurementsA,grpsA,taxaA,mytreeA,testn,val=lambdaval)
      rslt[i,] <- c(lambdaval,rss$RSS,rss$AICY) #replaced lLy with AICY
    }
    
    optlambda <- matrix(,1,1);colnames(optlambda)<- "RSS"
    optlambda[1,1]<-max(rslt[which(rslt[,2]==min(rslt[,2])),1])
    if(plot){
      matplot(rslt[,1],rslt[,2],type="l",xlab=expression(lambda),ylab="RSS",main="RSS",lty=1,col=1)
      abline(v=optlambda[1,1],col=2,lty=2);mtext(paste("Optimal Lambda = ",optlambda[1,1],sep=""))
    }
    
    if(export){
      pdf(height=11,width=6,file=paste(fldr,idc,".optLambda.pred.pdf",sep=''));layout(matrix(c(1,2),2,1))
      matplot(rslt[,1],rslt[,2],type="l",xlab=expression(lambda),ylab="RSS",main="RSS",lty=1,col=1)
      abline(v=optlambda[1,1],col=2,lty=2);mtext(paste("Optimal Lambda = ",optlambda[1,1],sep=""))
      dev.off()
    }
    list(optlambda=optlambda,rslt=rslt)
  }
}

# Defining colors and symbols

col_HF<-"darkgoldenrod"
col_LF<-contrasting.palette()[3]
col_LE<-"blue"
col_nasal<-"#F67BF6"
col_oral<-"#8B1EF0"
col_nonLE<-col_ptero<-"cyan"
col_nonbat<-"black"
col_VIE_fig3<-"#F37ECF"
col_VIE<-"firebrick1"
pch_nasal<-24
pch_oral<-25
pch_nonLE<-pch_yang<-22
pch_ptero<-25
pch_rhino<-24
pch_nonbat<-21

############### I - PHYLOGENETIC DATA
############### I - 1 - This study trees (hereafter called 'our trees')

# Important the trees and converting them to 'phylo' object to deal with lighter objects
mrb_trees<-read.mrbayes("1 Vielasia_post-burnin_trees.tre")
trees<-list()
for (i in 1:length(mrb_trees)){
  trees[[i]]<-as.phylo(mrb_trees[[i]])
}

# Pruning the trees to keep those we have further (i.e., morphometric) data on
to_rm_tips<-trees[[1]]$tip.label[c(1:3,7:8,15:16,18:30)]
pruned_trees<-list()
for (i in 1:length(trees)){
  pruned_trees[[i]]<-read.tree(text=write.tree(ladderize(drop.tip(trees[[i]],to_rm_tips))))
}

# Getting the ages for the Vielase bat and the stem bat node for each tree
stem_bat_tax_age<-numeric(length=length(pruned_trees))
stem_bat_node_age<-numeric(length=length(pruned_trees))
crown_bat_age<-numeric(length=length(pruned_trees))
for (i in 1:length(pruned_trees)){
  ages<-max(nodeHeights(pruned_trees[[i]]))-nodeHeights(pruned_trees[[i]])
  stem_bat_tax<-which(pruned_trees[[i]]$tip.label=="Vielase_bat")
  stem_bat_tax_age[i]<-ages[which(pruned_trees[[i]]$edge[,2]==stem_bat_tax),2]
  stem_bat_node<-pruned_trees[[i]]$edge[which(pruned_trees[[i]]$edge[,2]==which(pruned_trees[[i]]$tip.label=="Vielase_bat")),1]
  stem_bat_node_age[i]<-ages[which(pruned_trees[[i]]$edge[,2]==stem_bat_node),2]
  bats<-pruned_trees[[i]]$tip.label[getDescendants(pruned_trees[[i]],stem_bat_node)[getDescendants(pruned_trees[[i]],stem_bat_node)<=Ntip(pruned_trees[[i]])]]
  crown_bat_node<-getClade(pruned_trees[[i]],bats[-which(bats=="Vielase_bat")])
  crown_bat_age[i]<-ages[which(pruned_trees[[i]]$edge[,2]==crown_bat_node),2]                
}

############### I - 2 - Time-calibrated bat phylogenies of Shi & Rabosky 2015 and of Amador et al 2018

# Importing the dated phylogenies
SR_phylo<-read.tree("2 SR_dated_tree.tre")
AM_phylo<-read.tree("3 Amador_dated_tree.tre")

# In the phylogeny of Amador et al. 2018, correct one of the node which is younger than one of its two descendants 
# by giving to it the same age than its descendant + the half of the minimal internodal time gap of the tree
AM_ages<-max(nodeHeights(AM_phylo))-nodeHeights(AM_phylo)
node_neg_edge<-AM_phylo$edge[which((round(AM_ages[,1]-AM_ages[,2],10)-round(AM_phylo$edge.length[AM_phylo$edge.length<0],10))==0),1]
min_internode_dist<-min(c(AM_ages[,1]-AM_ages[,2])[which(c(AM_ages[,1]-AM_ages[,2])>0&AM_phylo$edge[,2]>Ntip(AM_phylo))])
value_to_set<-AM_ages[which(AM_phylo$edge[,1]==node_neg_edge),2][AM_ages[which(AM_phylo$edge[,1]==node_neg_edge),1]<AM_ages[which(AM_phylo$edge[,1]==node_neg_edge),2]]+min_internode_dist/2
AM_ages[which(AM_phylo$edge[,1]==node_neg_edge),1]<-AM_ages[which(AM_phylo$edge[,2]==node_neg_edge),2]<-value_to_set
AM_phylo$edge.length<-AM_ages[,1]-AM_ages[,2]

############### I - 3 - Mixing the two previous molecular scaffolds and our information about stem bat taxon & node ages

# 1st way, adding Vielasia to the molecular trees when their crown bat node age is younger than our stem bat node age
crown_bats_SR_phylo<-getClade(SR_phylo,SR_phylo$tip.label[-c(813:815)])
crown_bats_AM_phylo<-getClade(AM_phylo,AM_phylo$tip.label[-c(805:812)])
ages_SR_phylo<-max(nodeHeights(SR_phylo))-nodeHeights(SR_phylo)
ages_AM_phylo<-max(nodeHeights(AM_phylo))-nodeHeights(AM_phylo)
crown_bat_ages_SR<-ages_SR_phylo[which(SR_phylo$edge[,2]==crown_bats_SR_phylo),2]
crown_bat_ages_AM<-ages_AM_phylo[which(AM_phylo$edge[,2]==crown_bats_AM_phylo),2]

conc_topo_SR<-list()
conc_topo_AM<-list()
j<-0
k<-0
for (i in 1:length(pruned_trees)){
  if(stem_bat_node_age[i]>=crown_bat_ages_SR){
    j<-j+1
    conc_topo_SR[[j]]<-bind.tip(tree=SR_phylo,tip.label="Vielasia_sigei",edge.length=stem_bat_node_age[i]-stem_bat_tax_age[i],position=stem_bat_node_age[i]-crown_bat_ages_SR,where=crown_bats_SR_phylo)
    conc_topo_SR[[j]]<-drop.tip(conc_topo_SR[[j]],SR_phylo$tip.label[c(813:815)])
  }
  if(stem_bat_node_age[i]>=crown_bat_ages_AM){
    k<-k+1
    conc_topo_AM[[k]]<-bind.tip(tree=AM_phylo,tip.label="Vielasia_sigei",edge.length=stem_bat_node_age[i]-stem_bat_tax_age[i],position=stem_bat_node_age[i]-crown_bat_ages_AM,where=crown_bats_AM_phylo)
    conc_topo_AM[[k]]<-drop.tip(conc_topo_AM[[k]],AM_phylo$tip.label[c(805:812)])
  }
}

# 2nd way, consider our crown node age as true and rescale all node ages of both molecular scaffolds for each of our trees

rescaled_conc_topo_SR<-list()
rescaled_conc_topo_AM<-list()
for (i in 1:length(pruned_trees)){
  temp_SR_phylo<-drop.tip(SR_phylo,SR_phylo$tip.label[c(813:815)])
  temp_ages_SR_phylo<-max(nodeHeights(temp_SR_phylo))-nodeHeights(temp_SR_phylo)
  rescaled_ages_SR_phylo<-scales::rescale(temp_ages_SR_phylo,to=c(min(temp_ages_SR_phylo),crown_bat_age[i]),from=c(min(temp_ages_SR_phylo),max(temp_ages_SR_phylo)))
  temp_SR_phylo$edge.length<-rescaled_ages_SR_phylo[,1]-rescaled_ages_SR_phylo[,2]
  rescaled_conc_topo_SR[[i]]<-read.tree(text=paste0("(",strsplit(write.tree(temp_SR_phylo),split=";")[[1]],":",stem_bat_node_age[i]-crown_bat_age[i],",Vielasia_sigei:",stem_bat_node_age[i]-stem_bat_tax_age[i],");"))
  
  temp_AM_phylo<-drop.tip(AM_phylo,AM_phylo$tip.label[c(805:812)])
  temp_ages_AM_phylo<-max(nodeHeights(temp_AM_phylo))-nodeHeights(temp_AM_phylo)
  rescaled_ages_AM_phylo<-scales::rescale(temp_ages_AM_phylo,to=c(min(temp_ages_AM_phylo),crown_bat_age[i]),from=c(min(temp_ages_AM_phylo),max(temp_ages_AM_phylo)))
  temp_AM_phylo$edge.length<-rescaled_ages_AM_phylo[,1]-rescaled_ages_AM_phylo[,2]
  rescaled_conc_topo_AM[[i]]<-read.tree(text=paste0("(",strsplit(write.tree(temp_AM_phylo),split=":0;")[[1]],":",stem_bat_node_age[i]-crown_bat_age[i],",Vielasia_sigei:",stem_bat_node_age[i]-stem_bat_tax_age[i],");"))
}
all_conc_SR<-c(conc_topo_SR,rescaled_conc_topo_SR[-which(mapply(function(x){any(x$edge.length==0)},rescaled_conc_topo_SR))])
all_conc_AM<-c(conc_topo_AM,rescaled_conc_topo_AM[-which(mapply(function(x){any(x$edge.length==0)},rescaled_conc_topo_AM))])

############### II - BONY LABYRINTH PARAMETERS INVESTIGATIONS USING DATA FROM DAVIES ET AL. 2013A,B
############### II - 1 - Get phylogenetic data at the Mammalia scale

# Input mammalian-scale phylogeny of Alvarez-Carretero et al. 2021
mamm_phylo<-read.tree("4 Alvarez-Carretero et al 2021 Nature_mammalian phylogeny.tre")
for (i in 1:Ntip(mamm_phylo)){
  mamm_phylo$tip.label[i]<-paste(toupper(substr(mamm_phylo$tip.label[i], 1, 1)), substr(mamm_phylo$tip.label[i], 2, nchar(mamm_phylo$tip.label[i])), sep="")
}

# Drop subspecies for species of the datasets
mamm_phylo<-drop.tip(mamm_phylo,which(substr(mamm_phylo$tip.label,1,16)=="Saimiri_sciureus")[1:4])
mamm_phylo$tip.label[which(substr(mamm_phylo$tip.label,1,16)=="Saimiri_sciureus")]<-"Saimiri_sciureus"
mamm_phylo<-drop.tip(mamm_phylo,which(substr(mamm_phylo$tip.label,1,21)=="Cryptomys_hottentotus")[1:6])
mamm_phylo$tip.label[which(substr(mamm_phylo$tip.label,1,21)=="Cryptomys_hottentotus")]<-"Cryptomys_hottentotus"

# Get phylogeny ages in My
mamm_phylo$edge.length<-mamm_phylo$edge.length*100

############### II - 2 - Investigations on basilar membrane, number of cochlear turns, and body mass (after Davies et al. 2013a)

# Get the data
ear_data<-read.table("5 ear_data.txt",header=TRUE,dec=",",sep="\t")

############### II - 2 - a - Mixing phylogenetic data at the Mammalia and Chiroptera scales for the given dataset

# -> The goal here is to substitute the Chiroptera subtree of the mammalian-scale phylogeny of 
#    Alvarez-Carretero et al. 2021 by the phylogenies built in the section I

# In the 'original' topologies for both Chiroptera scale molecular scaffolds
conc_topo_SR_ear<-subs.tree.to.ref(mamm_phylo,conc_topo_SR,ear_data$Species)
conc_topo_AM_ear<-subs.tree.to.ref(mamm_phylo,conc_topo_AM,ear_data$Species)

# In the 'rescaled' topologies for both Chiroptera scale molecular scaffolds
rescaled_conc_topo_SR_ear<-subs.tree.to.ref(mamm_phylo,rescaled_conc_topo_SR,ear_data$Species)
rescaled_conc_topo_AM_ear<-subs.tree.to.ref(mamm_phylo,rescaled_conc_topo_AM,ear_data$Species)

# Assemble 'original' and 'rescaled' topologies for both Chiroptera scale molecular scaffolds
all_conc_SR_ear<-c(if(length(conc_topo_SR_ear)>0){conc_topo_SR_ear},
                   rescaled_conc_topo_SR_ear[-which(mapply(function(x){any(x$edge.length==0)},rescaled_conc_topo_SR_ear))])
all_conc_AM_ear<-c(if(length(conc_topo_AM_ear)>0){conc_topo_AM_ear},
                   rescaled_conc_topo_AM_ear[-which(mapply(function(x){any(x$edge.length==0)},rescaled_conc_topo_AM_ear))])

# Keep only taxa presents in the dataset, in the mammalian-scale scaffold for non-bats mammals, and in the two
# chiropteran-scale scaffolds for bats
both_SR_AM_taxa_ear<-intersect(unique(unlist(mapply(function(x){x$tip.label},all_conc_SR_ear))),unique(unlist(mapply(function(x){x$tip.label},all_conc_AM_ear))))
ear_data<-ear_data[which(ear_data$Species%in%both_SR_AM_taxa_ear),]
all_conc_SR_ear<-lapply(all_conc_SR_ear,function(x){keep.tip(x,both_SR_AM_taxa_ear)})
all_conc_AM_ear<-lapply(all_conc_AM_ear,function(x){keep.tip(x,both_SR_AM_taxa_ear)})
ear_phylo_full<-c(all_conc_SR_ear,all_conc_AM_ear)
class(ear_phylo_full)<-"multiPhylo"

############### II - 2 - b - Compute PGLS between log basilar membrane length and log body masses for all phylogenies

# -> The goal here is to compute these PGLS for all phylogenies under different evolutionary models, to retain the
#    set of PGLS under the best-fitting model, and then to contrast average regressions between groups we want to compare

# First, the PGLS considering the whole dataset
reg_full<-compare.tree_phylm(x=log(ear_data[,5]),y=log(ear_data[,3]),xynames=ear_data[,2],multiphy=ear_phylo_full,subset.split.name=c("all"),output=c("all_regs","red_regs","TT_data","ablines"))

# Second, contrasting chiropteran and non-chiropteran mammals
only_bats<-extract.from.factor(ear_data$Species,intersect(all_conc_SR[[1]]$tip.label,all_conc_AM[[1]]$tip.label),ear_data$Species)
reg_bats<-compare.tree_phylm(x=log(ear_data[,5]),y=log(ear_data[,3]),xynames=ear_data[,2],multiphy=ear_phylo_full,subset=only_bats,subset.split.name=c("bats","non-bats"),output=c("red_regs","TT_data","TT","ablines"))

# Third, contrasting echolocating bats and other mammals (i.e., non-bat mammals + non-echolocating bats)
only_echobats<-extract.from.factor(ear_data$Echolocation.call.type,c("Nasal","Oral"),ear_data$Species)
reg_echobats<-compare.tree_phylm(x=log(ear_data[,5]),y=log(ear_data[,3]),xynames=ear_data[,2],multiphy=ear_phylo_full,subset=only_echobats,subset.split.name=c("echolocating bats","other mammals"),output=c("red_regs","TT_data","TT","ablines"))

# Fourth, contrasting high-frequency hearing mammals (i.e., echolocating bats + odontocetes) and other mammals
only_HF<-extract.from.factor(ear_data$Echolocation.call.type,c("Nasal","Oral","HF-cet"),ear_data$Species)
reg_HF<-compare.tree_phylm(x=log(ear_data[,5]),y=log(ear_data[,3]),xynames=ear_data[,2],multiphy=ear_phylo_full,subset=only_HF,subset.split.name=c("HF-hearing mammals","other mammals"),output=c("red_regs","TT_data","TT","ablines"))

# Fifth, contrasting high- and low-frequency hearing mammals (i.e., mammals with particular acoustic adaptation,
# namely echolocating bats + cetaceans) and other mammals
only_HFLF<-extract.from.factor(ear_data$Echolocation.call.type,c("Nasal","Oral","HF-cet","LF-cet"),ear_data$Species)
reg_HFLF<-compare.tree_phylm(x=log(ear_data[,5]),y=log(ear_data[,3]),xynames=ear_data[,2],multiphy=ear_phylo_full,subset=only_HFLF,subset.split.name=c("HF-.and.LF-hearing mammals","other mammals"),output=c("red_regs","TT_data","TT","ablines"))

# Summary of the t-tests in each opposition regarding PGLS regression slopes parameters
groups_PGLS_TT<-rbind(
  foreach(i=1:4,.combine="rbind")%do%{
    test_pval<-foreach(j=c(1,3),.combine="c")%do%{
      c(list(reg_bats$TT),list(reg_echobats$TT),list(reg_HF$TT),list(reg_HFLF$TT))[[i]][[j]]$p.value
    }
    test_pval<-t(test_pval)
    rownames(test_pval)<-paste(strsplit(c(list(reg_bats$TT),list(reg_echobats$TT),list(reg_HF$TT),list(reg_HFLF$TT))[[i]][[j]]$data.name," and ")[[1]],collapse=" vs ")
    test_pval
  },
  foreach(i=1:4,.combine="rbind")%do%{
    test_pval<-foreach(j=1:2,.combine="rbind")%do%{
      foreach(k=1:2,.combine="cbind")%do%{
        do.call("param.t.test",c(setNames(c(as.list(reg_full$TT_data[[k]]$num.data),as.list(c(list(reg_bats$TT_data),list(reg_echobats$TT_data),list(reg_HF$TT_data),list(reg_HFLF$TT_data))[[i]][[k]]$num.data[c(1:3)+3*(j-1)])),
                                          c("mean_x","var_x","n_x","mean_y","var_y","n_y")),
                                 setNames(c(as.list(reg_full$TT_data[[k]]$names),as.list(c(list(reg_bats$TT_data),list(reg_echobats$TT_data),list(reg_HF$TT_data),list(reg_HFLF$TT_data))[[i]][[k]]$names[j])),
                                          c("name_x","name_y"))))$p.value
      }
    }
    curr_name<-c(list(reg_bats$TT_data),list(reg_echobats$TT_data),list(reg_HF$TT_data),list(reg_HFLF$TT_data))[[i]][[k]]$names[1]
    rownames(test_pval)<-paste("whole dataset vs ",c(curr_name,paste0("non-",curr_name)),sep="")
    test_pval
  })
colnames(groups_PGLS_TT)<-c("intercepts","slopes")

# write.table(groups_PGLS_TT,file="Table S1H - Comparisons between PGLS regression line parameters.txt",append=FALSE,quote=FALSE,sep="\t")

############### II - 2 - c - Regression graphics with the different groups regression lines

echo_pch_reg<-setNames(c(pch_nasal,pch_oral,pch_nonLE),c("Nasal","Oral","Non"))
ear_bat_groups<-ear_bat_groups<-ear_data$Echolocation.call.type[ear_data$Species%in%only_bats]
for (i in 1:5){
  if(i==1){
    # postscript("Fig. 3B - Log basilar membrane length (mm) versus log body mass (g) in extant mammals.eps")
    # postscript("Fig. 3B bis - Log basilar membrane length (mm) versus log body mass (g) in extant mammals.eps")
  }
  if(i==2){
    # postscript("Fig. S3A - Allometric relationship of log basilar membrane length and log body mass across placental mammals.eps")
    par(mfrow=c(2,2),mar=rep(3,4),mgp=c(2,0.5,0))
  }
  
  reg_groups_chu<-reg_groups_pch<-as.character(ear_bat_groups)
  reg_groups_chu[reg_groups_chu%in%c("Nasal","Oral")]<-"LE"
  reg_groups_chu<-factor(reg_groups_chu,c("LE","Non","unknown"))
  reg_groups_pch<-factor(reg_groups_pch,c("Nasal","Oral","Non","unknown"))
  levels(reg_groups_pch)[levels(reg_groups_pch)%in%names(echo_pch_reg)]<-echo_pch_reg[match(levels(reg_groups_pch)[levels(reg_groups_pch)%in%names(echo_pch_reg)],names(echo_pch_reg))]
  levels(reg_groups_pch)[levels(reg_groups_pch)=="unknown"]<-NA
  reg_groups_pch<-as.numeric(as.character(reg_groups_pch))
  
  par(pty="s")
  plot(log(ear_data[,5]),log(ear_data[,3]),xlab="log body mass (g)",ylab="log basilar membrane length (mm)",type="n")
  points(log(ear_data[!ear_data$Species%in%only_bats,5]),log(ear_data[!ear_data$Species%in%only_bats,3]),pch=21,col=NA,bg=col_nonbat,cex=ear_data[!ear_data$Species%in%only_bats,4]/2)
  points(log(ear_data[which(ear_data[,6]=="LF-cet"),5]),log(ear_data[which(ear_data[,6]=="LF-cet"),3]),pch=21,col=NA,bg=col_LF,cex=ear_data[which(ear_data[,6]=="LF-cet"),4]/2)
  points(log(ear_data[which(ear_data[,6]=="HF-cet"),5]),log(ear_data[which(ear_data[,6]=="HF-cet"),3]),pch=21,col=NA,bg=col_HF,cex=ear_data[which(ear_data[,6]=="HF-cet"),4]/2)
  if(i==1){
    do.call("abline",c(reg_full$ablines,list(lwd=2,col="black")))
    
    reg_PI<-list()
    for (j in 1:length(ear_phylo_full)){ # Prediction interval using the function gls.pi from evomap, and providing a variance-covariance matrix of the tree taking into account the value of Pagel's delta (the best-fitting model) used for each PGLS
      curr_phy<-ear_phylo_full[[j]]
      mod<-reg_full$all_regs$all.results[[which.min(reg_full$red_regs["AIC",])]][[j]]
      X<-mod$X[,2]
      Y<-mod$y
      val_d<-mod$optpar
      
      ht<-((max(nodeHeights(curr_phy))-nodeHeights(curr_phy))[order(curr_phy$edge[,2]),])
      ht<-rbind(ht[1:Ntip(curr_phy),],rep(max(nodeHeights(curr_phy)),2),ht[(Ntip(curr_phy)+1):(Ntip(curr_phy)+Nnode(curr_phy)-1),])
      ht<-as.data.frame(ht)
      colnames(ht)<-c("start","end")
      N = Ntip(curr_phy)
      Tmax = ht$start[N + 1]
      ht$t = Tmax - ht$end
      ht$e = ht$start - ht$end
      ht$a = ht$t - ht$e
      bl = (ht$a + ht$e)^val_d - ht$a^val_d
      curr_phy$edge.length = bl[curr_phy$edge[, 2]]
      Sigma<-vcv(curr_phy)
      
      k<-1e6
      n <- length(X)
      tr <- sum(diag(Sigma))
      Sigma <- n * Sigma/tr
      invSigma <- solve(Sigma)
      Yhat<-mod$fitted.values
      Yresid = Y - Yhat
      
      XX <- cbind(rep(1, n), X)
      C <- solve(t(XX) %*% invSigma %*% XX)
      new_n<-1000
      new_X<-cbind(rep(1, new_n), seq(min(X)-max(X),max(X)*2,length.out=new_n))
      SEYhat <- sqrt((diag(new_X %*% C %*% t(new_X)) + c(1/k)) %*% ((t(Yresid) %*% invSigma %*% Yresid)/(n - 2)))
      new_Yhat<-new_X[,2]*mod$coefficients[2]+mod$coefficients[1]
      PI_low <- (new_Yhat - qt(0.975, n) * SEYhat)[order(new_X[,2])]
      PI_up <- (new_Yhat + qt(0.975, n) * SEYhat)[order(new_X[,2])]
      reg_PI[[j]]<-cbind(sort(new_X[,2]),PI_low,PI_up)
    }
    avg_PI<-reg_PI[[1]]
    for(j in 1:nrow(avg_PI)){
      for(k in 1:ncol(avg_PI)){
        avg_PI[j,k]<-mean(mapply(function(x){x[j,k]},reg_PI))
      }
    }
    polygon(c(avg_PI[,1],rev(avg_PI[,1])),c(avg_PI[,2],rev(avg_PI[,3])),col=NA,density=-10)
    
  }
  else if(i==2){
    for(j in 1:2){do.call("abline",c(reg_bats$ablines[[j]],list(lwd=2,col="black",lty=j)))}
    title("bats vs non-bats")
  }
  else if(i==3){
    for(j in 1:2){do.call("abline",c(reg_echobats$ablines[[j]],list(lwd=2,col="black",lty=j)))}
    title("echolocating bats vs others")
  }
  else if(i==4){
    for(j in 1:2){do.call("abline",c(reg_HF$ablines[[j]],list(lwd=2,col="black",lty=j)))}
    title("HF-hearing mammals vs others")
  }
  else{
    for(j in 1:2){do.call("abline",c(reg_HFLF$ablines[[j]],list(lwd=2,col="black",lty=j)))}
    title("HF- and LF-hearing mammals vs others")
  }
  Vielasia_graph(x=log(ear_data[ear_data$Species%in%only_bats,5]),
                 y=log(ear_data[ear_data$Species%in%only_bats,3]),
                 to_rm=0,
                 groups=reg_groups_chu,
                 Vielasia=which(ear_data[ear_data$Species%in%only_bats,1]=="stem bat"),
                 opt=1,
                 cols=c(col_LE,col_nonLE),
                 xlab="log body mass",
                 ylab="log basilar membrane length",
                 cex=ear_data[ear_data$Species%in%only_bats,4]/2,
                 pch=reg_groups_pch,
                 new=FALSE,
                 legend=FALSE,
                 Vielasia_col = col_VIE_fig3)
  # points(log(c(16.09,23.71)),rep(log(ear_data[ear_data$Species%in%only_bats,3][which(ear_data[ear_data$Species%in%only_bats,1]=="stem bat")]),2),type="l",pch=21,col="red",bg="red",lwd=3)
  if(i==1|i==5){
    # dev.off()
    par(pty="m",mfrow=c(1,1),mgp=c(3,1,0),mar=c(5.1,4.1,4.1,2.1))
  }
}

############### II - 2 - d - Calculating low and high frequency limits for Vielasia sigei

freq_ear_data<-read.table("6 frequency limits.txt",header=TRUE,dec=",",sep="\t")
freq_groupcols<-c("#F8A400",rep("#444CA9",3),rep("#13A629",10),rep("#FFCDDC",2),rep("#8F3E14",18),rep("#C6E5F0",2),rep("#FFFFFF",2))

freq_ear_with<-freq_ear_data[c(1:14,17:34),]
freq_ear_nobats<-freq_ear_data[c(17:34),]

Vielasia_logrbm<-log10(ear_data[ear_data[,2]=="Vielasia_sigei",3]/c(ear_data[ear_data[,2]=="Vielasia_sigei",5]^0.33))
Vielasia_logrbmturns<-log10(ear_data[ear_data[,2]=="Vielasia_sigei",3]*ear_data[ear_data[,2]=="Vielasia_sigei",4])

freqs_Vielasia<-freqs_davies<-matrix(ncol=4,nrow=2)
colnames(freqs_Vielasia)<-colnames(freqs_davies)<-c("Low limit (30dB)","Low limit (60dB)","High limit (30dB)","High limit (60dB)")
rownames(freqs_Vielasia)<-rownames(freqs_davies)<-c("mammals with bats","mammals without bats")

davies_estimates<-c(11.75,2.45,19.22)

# postscript("Fig. S4 - Estimated hearing limits of Vielasia sigei.eps")
par(pty="s",mfrow=c(2,2),mar=c(3,2,1,1),mgp=c(2,1,0))

for(i in 1:4){
  curr_lm_full<-lm(log10(freq_ear_with[,c(4,2,5,3)[i]])~log10(freq_ear_with[,6]*freq_ear_with[,ifelse(i==1|i==3,8,7)]^(ifelse(i==1|i==3,-1/3,1))))
  curr_lm_nobats<-lm(log10(freq_ear_nobats[,c(4,2,5,3)[i]])~log10(freq_ear_nobats[,6]*freq_ear_nobats[,ifelse(i==1|i==3,8,7)]^(ifelse(i==1|i==3,-1/3,1))))
  plot(log10(freq_ear_data[,c(4,2,5,3)[i]])~log10(freq_ear_data[,6]*freq_ear_data[,ifelse(i==1|i==3,8,7)]^(ifelse(i==1|i==3,-1/3,1))),xlim=as.numeric(strsplit(ifelse(i==1|i==3,c("-0.5%1"),c("1%2.5")),"%")[[1]]),ylim=as.numeric(strsplit(ifelse(i==1|i==3,c("-0.5%2.5"),c("-2%2")),"%")[[1]]),pch=21,col="black",bg=freq_groupcols,cex=2,xlab=ifelse(i==1|i==3,"Log relative basilar membrane","Log basilar membrane x turns"),ylab=paste0("Log ",ifelse(i==1|i==3,"high","low")," frequency limit (",ifelse(i<3,"3","6"),"0dB)"))
  abline(curr_lm_full,lwd=2)
  abline(curr_lm_nobats,lwd=2,lty="43")
  curr_x_Vielasia<-ifelse(i==1|i==3,Vielasia_logrbm,Vielasia_logrbmturns)
  curr_y_Vielasia<-c(curr_x_Vielasia*coefficients(curr_lm_full)[2]+coefficients(curr_lm_full)[1],
                     curr_x_Vielasia*coefficients(curr_lm_nobats)[2]+coefficients(curr_lm_nobats)[1])
  points(rep(curr_x_Vielasia,2),curr_y_Vielasia,col=col_VIE,bg=col_VIE,pch=23,cex=1)
  
  curr_x_davies<-log10(davies_estimates[1]*davies_estimates[ifelse(i==1|i==3,3,2)]^(ifelse(i==1|i==3,-1/3,1)))
  curr_y_davies<-c(curr_x_davies*coefficients(curr_lm_full)[2]+coefficients(curr_lm_full)[1],
                   curr_x_davies*coefficients(curr_lm_nobats)[2]+coefficients(curr_lm_nobats)[1])
  segments(curr_x_davies,-10,curr_x_davies,max(curr_y_davies),lty="32")
  segments(curr_x_Vielasia,-10,curr_x_Vielasia,max(curr_y_Vielasia),lty="32",col=col_VIE)
  for(j in 1:2){
    segments(-10,curr_y_davies[j],curr_x_davies,curr_y_davies[j],lty="32")
    segments(-10,curr_y_Vielasia[j],curr_x_Vielasia,curr_y_Vielasia[j],lty="32",col=col_VIE)
  }
  
  freqs_Vielasia[,c(3,1,4,2)[i]]<-exp(log(10)*curr_y_Vielasia)
  freqs_davies[,c(3,1,4,2)[i]]<-exp(log(10)*curr_y_davies)
}
# dev.off()

h_estimates<-rbind(freqs_davies,freqs_Vielasia)
rownames(h_estimates)<-paste0(c(rep("Davies et al. 'bat ancestor'",2),rep("Vielasia sigei",2))," - ",rownames(h_estimates))
# write.table(h_estimates,file="Table S1L - Limits of hearing frequencies at 30 and 60 dB.txt",append=FALSE,quote=FALSE,sep="\t")

############### II - 3 - Investigations on the semicircular canals sizes vs cochlea size relationships

# Get the data
bony_lab<-read.table("7 bony_lab.txt",header=TRUE,sep="\t",dec=",")
bony_lab<-na.omit(bony_lab)

############### II - 3 - a - Mixing phylogenetic data at the Mammalia and Chiroptera scales for the given dataset

# -> Again, the goal is to substitute the Chiroptera subtree of the mammalian-scale phylogeny of 
#    Alvarez-Carretero et al. 2021 by the phylogenies built in the section I

# In the 'original' topologies for both Chiroptera scale molecular scaffolds
conc_topo_SR_BL<-subs.tree.to.ref(mamm_phylo,conc_topo_SR,bony_lab$Species)
conc_topo_SR_AM<-subs.tree.to.ref(mamm_phylo,conc_topo_AM,bony_lab$Species)

# In the 'rescaled' topologies for both Chiroptera scale molecular scaffolds
rescaled_conc_topo_SR_BL<-subs.tree.to.ref(mamm_phylo,rescaled_conc_topo_SR,bony_lab$Species)
rescaled_conc_topo_AM_BL<-subs.tree.to.ref(mamm_phylo,rescaled_conc_topo_AM,bony_lab$Species)

# Assemble 'original' and 'rescaled' topologies for both Chiroptera scale molecular scaffolds
all_conc_SR_BL<-c(if(length(conc_topo_SR_BL)>0){conc_topo_SR_BL},
                  rescaled_conc_topo_SR_BL[-which(mapply(function(x){any(x$edge.length==0)},rescaled_conc_topo_SR_BL))])
all_conc_AM_BL<-c(if(length(conc_topo_SR_AM)>0){conc_topo_SR_AM},
                  rescaled_conc_topo_AM_BL[-which(mapply(function(x){any(x$edge.length==0)},rescaled_conc_topo_AM_BL))])

# Keep only taxa presents in the dataset, in the mammalian-scale scaffold for non-bats mammals, and in the two
# chiropteran-scale scaffolds for bats
both_SR_AM_taxa_BL<-intersect(unique(unlist(mapply(function(x){x$tip.label},all_conc_SR_BL))),unique(unlist(mapply(function(x){x$tip.label},all_conc_AM_BL))))
bony_lab<-bony_lab[which(bony_lab$Species%in%both_SR_AM_taxa_BL),]
all_conc_SR_BL<-lapply(all_conc_SR_BL,function(x){keep.tip(x,both_SR_AM_taxa_BL)})
all_conc_AM_BL<-lapply(all_conc_AM_BL,function(x){keep.tip(x,both_SR_AM_taxa_BL)})
bony_lab_phylo<-c(all_conc_SR_BL,all_conc_AM_BL)
bony_lab_phylo<-lapply(bony_lab_phylo,ladderize)
class(bony_lab_phylo)<-"multiPhylo"

############### II - 3 - b - Simply reporting boxplots of semicircular canal size and cochlea size for extant bats
###############               vs the value of the stem bat Vielasia sigei

only_bats_BL<-extract.from.factor(bony_lab$Species,intersect(all_conc_SR[[1]]$tip.label,all_conc_AM[[1]]$tip.label),bony_lab$Species)
bony_lab_bats<-bony_lab[bony_lab$Species%in%only_bats_BL,]
bony_lab_bats$SuborderEcho<-bony_lab_bats$Echolocation.call.type
bony_lab_bats$Echolocation.call.type<-mapply(function(x){x[[1]]},strsplit(bony_lab_bats$Echolocation.call.type,";"))
bony_lab$Echolocation.call.type<-mapply(function(x){x[[1]]},strsplit(bony_lab$Echolocation.call.type,";"))
bony_lab_bats$Echolocation.call.type<-factor(bony_lab_bats$Echolocation.call.type,levels=levels(factor(bony_lab_bats$Echolocation.call.type)))
bony_lab$Echolocation.call.type<-factor(bony_lab$Echolocation.call.type,levels=levels(factor(bony_lab$Echolocation.call.type)))
bony_lab_bats$SuborderEcho<-mapply(function(x){x[[length(x)]]},strsplit(bony_lab_bats$SuborderEcho,";"))
bony_lab_bats$SuborderEcho<-factor(bony_lab_bats$SuborderEcho,levels=levels(factor(bony_lab_bats$SuborderEcho))[c(3,4,1,2)])

# postscript("Fig. 3C - Semicircular canals and cochlea size in extant bats and in Vielasia.eps")
layout(matrix(ncol=6,nrow=1,c(3,1,1,1,2,3),byrow=TRUE))
par(pty="m",mar=c(2,5.1,1,5.1),las=1,cex.axis=2)

bp_groups<-as.character(bony_lab_bats$Echolocation.call.type)
bp_groups[bp_groups%in%c("Nasal","Oral")]<-"LE"
bp_groups<-factor(bp_groups,c("LE","Non","unknown"))

boxgroups(bony_lab_bats[,2:4],bp_groups,box.width=0.3,points.opt=list("col"=c(col_LE,col_nonLE,col_VIE_fig3),"bg"=c(col_LE,col_nonLE,col_VIE_fig3),"pch"=21,"lwd"=3,cex=2),names=c("ASC","LSC","PSC"),box.opt=list("frame"=FALSE,"yaxt"="n",cex=2))
axis(2,cex.axis=1)
axis(2,mean(range(unlist(bony_lab_bats[,2:4]))),labels="Radius",tick=FALSE,cex.axis=2)
axis(4,mean(range(unlist(bony_lab_bats[,2:4]))),labels="Size",tick=FALSE,hadj=0.25,cex.axis=2)
par(pty="m",mar=c(2,0.1,1,2.1))
boxgroups(as.matrix(bony_lab_bats[,5]),bp_groups,box.width=0.3,x.gap=c(0.025,0.1),points.opt=list("col"=c(col_LE,col_nonLE,col_VIE_fig3),"bg"=c(col_LE,col_nonLE,col_VIE_fig3),"pch"=21,"lwd"=3,cex=2),names=c("Cochlea"),box.opt=list("frame"=FALSE,"yaxt"="n",cex=2))
axis(2,at=seq(1,5,1),cex.axis=1)
par(pty="m",mfrow=c(1,1),mgp=c(3,1,0),mar=c(5.1,4.1,4.1,2.1),cex.axis=1)
# dev.off()

############### II - 3 - c - Compute PGLS between each log semicircular canal size and log cochlea size for each 
###############               phylogeny, and make the PGLS parameters of the three regressions to converge

# -> The goal here is to obtain a single value for the PGLS parameters (Pagel's kappa, lambda, and delta) for all
#    three PGLS regressions of semicircular canal size against cochlea size. This allows then to perform a pFDA
#    (phylogenetic flexible discriminant analysis) encompassing the residuals of the three PGLS regressions satisfying
#    the required unique kappa/lambda/delta value for the whole dataset.

# Performing the PGLS regressions for each phylogeny and make the PGLS parameters values to converge
to_use<-list()
for (i in 1:length(bony_lab_phylo)){
  start<-Sys.time()
  to_use[[i]]<-conv.pgls(x=log(bony_lab[,5]),y=log(bony_lab[,2:4]),phy=bony_lab_phylo[[i]],names=bony_lab[,1],output=c("AIC_table","estimates.comparison","residuals","vals"),tol=1e-4,trace.conv=FALSE,silent=TRUE)
  end<-Sys.time()
  writeLines(paste(i,as.numeric(end-start)%/%1,round((as.numeric(end-start)-as.numeric(end-start)%/%1)*60,0),sep="_____"))
}

############### II - 3 - d - Perform a discriminant analysis taking into account phylogeny (pFDA) and taking not it (LDA)

# Defining ecological groups, based on echolocation call type
echo_groups<-factor(setNames(bony_lab$Echolocation.call.type,bony_lab$Species))

# Performing the pFDA using the unique PGLS parameters calculated before for each phylogeny
pfda_optL<-pfda_LSchmitz(lapply(to_use,function(x){x$residuals}),bony_lab_phylo,groups=echo_groups,lambda=lapply(to_use,function(x){round(x$vals["lambda"],5)}),pfda.opt=list(val_k=lapply(to_use,function(x){round(x$vals["kappa"],5)}),val_d=lapply(to_use,function(x){round(x$vals["delta"],5)}),KLD=rep(list(TRUE),length(bony_lab_phylo))))

# Performing a simple LDA for each phylogeny without accounting for phylogenetic signal (PGLS residuals for each
# phylogeny are used, but these are not decorrelated from phylogenetic relationships)
pfda_null<-rep_lda(lapply(to_use,function(x){x$residuals}),echo_groups,"Vielasia_sigei")

############### II - 3 - e - Get average results for each discriminant analysis

# Computing the average confusion matrix
avg_confusion_optL<-as.table(matrix(nrow=(nlevels(echo_groups)-1),ncol=(nlevels(echo_groups)-1),byrow=F,apply(mapply(function(x){x},pfda_optL[[1]]),1,mean)))
avg_confusion_null<-as.table(matrix(nrow=(nlevels(echo_groups)-1),ncol=(nlevels(echo_groups)-1),byrow=F,apply(mapply(function(x){x},pfda_null[[1]]),1,mean)))
avg_confusion_all<-as.table(foreach(i=1:dim(avg_confusion_optL)[1],.combine="rbind")%do%{
  foreach(j=1:dim(avg_confusion_optL)[2],.combine="cbind")%do%{
    avg_LDA<-round(avg_confusion_null[i,j],2)
    avg_pFDA<-round(avg_confusion_optL[i,j],2)
    cell<-ifelse(avg_LDA==0,ifelse(avg_pFDA==0,"",paste0(avg_pFDA," (pFDA)")),ifelse(avg_pFDA==0,paste0(avg_LDA," (LDA)"),paste0(avg_LDA," - ",avg_pFDA)))
  }
})
dimnames(avg_confusion_all)<-dimnames(avg_confusion_optL)<-dimnames(avg_confusion_null)<-list("predicted"=levels(factor(bony_lab$Echolocation.call.type))[1:(nlevels(echo_groups)-1)],
                                                                                              "true"=levels(factor(bony_lab$Echolocation.call.type))[1:(nlevels(echo_groups)-1)])
export_avg_confusion_all<-apply(avg_confusion_all,c(1,2),function(x){if(any(unlist(strsplit(x," - "))==0)){if(all(unlist(strsplit(x," - "))==0)){paste0("")}else{paste0(unlist(strsplit(x," - "))[unlist(strsplit(x," - "))!=0],ifelse(unlist(strsplit(x," - "))[1]==0," (pFDA)"," (LDA)"))}}else{x}})
# write.table(export_avg_confusion_all,file="Table S1I - Average confusion matrix for the performed discriminant analyses for the six placental groups shown.txt",append=FALSE,quote=FALSE,sep="\t")

# Contrasting these confusions to a 'no confusion' matrix
no_confusion<-matrix(nrow=(nlevels(echo_groups)-1),ncol=(nlevels(echo_groups)-1),0)
diag(no_confusion)<-table(pfda_optL$training.results[[1]][,1])

# Printing the average number of taxa erroneously assigned to a category (positive values) or that should have been 
# assigned (negative values)
avg_confusion_optL-no_confusion
avg_confusion_null-no_confusion

# Compute the average results (category assignation, probability for each category, and discriminant axes scores) for
# the 'test' sample, that is Vielasia sigei
avg_test.results_optL<-as.data.frame(t(apply(simplify2array(pfda_optL$test.results),c(1,2),averagize)))
avg_test.results_null<-as.data.frame(t(apply(simplify2array(pfda_null$test.results),c(1,2),averagize)))
avg_test.results_all<-rbind(avg_test.results_null[,1:7],avg_test.results_optL[,1:7])
avg_test.results_all[,2:7]<-apply(avg_test.results_all[,2:7],c(1,2),function(x){signif(as.numeric(x),4)})
rownames(avg_test.results_all)<-c("LDA","pFDA")
# write.table(avg_test.results_all,file="Table S1J - Average group pertaining prediction and probabilities for Vielasia for both types of discriminant analysis.txt",append=FALSE,quote=FALSE,sep="\t")

# Compute the average results (idem than previously) for the 'training' sample, that are all other mammals of the dataset
avg_training.results_optL<-as.data.frame(apply(simplify2array(pfda_optL$training.results),c(1,2),averagize))
avg_training.results_null<-as.data.frame(apply(simplify2array(pfda_null$training.results),c(1,2),averagize))

# Compute the average discriminant axes scores for the whole dataset and the average percents of explained variance
avg_DA_scores_optL<-as.data.frame(apply(simplify2array(pfda_optL$DA.scores),c(1,2),averagize))
avg_DA_scores_optL[,1:4]<-apply(avg_DA_scores_optL[,1:4],2,as.numeric)
avg_DA_scores_null<-as.data.frame(apply(simplify2array(pfda_null$DA.scores),c(1,2),averagize))
avg_DA_scores_null[,1:3]<-apply(avg_DA_scores_null[,1:3],2,as.numeric)
avg_pe_optL<-round(apply(pfda_optL$percent.explained,2,mean),1)
avg_pe_null<-round(apply(pfda_null$percent.explained,2,mean),1)

############### II - 3 - f - Get exhaustive information of each species position on the two first discriminant axes
###############               of each discriminant analysis by computing their 'full' and '95%' morphospaces and
###############               plot them

# Preparing groups and symbols (with LE bats grouped for groups, and Nasal/Oral-emitting bats distinguished by symbol) for both DA

echo_pch_da<-setNames(c(pch_nonbat,pch_nonbat,pch_nasal,pch_oral,pch_nonLE,pch_nonbat),c("HF-cet","LF-cet","Nasal","Oral","Non","non-bat"))

lda_groups_chu<-lda_groups_pch<-setNames(as.character(avg_DA_scores_null[,4]),rownames(avg_DA_scores_null))
lda_groups_chu[lda_groups_chu%in%c("Nasal","Oral")]<-"LE"
lda_groups_chu<-factor(lda_groups_chu,c("HF-cet","LF-cet","LE","Non","non-bat","unknown"))

lda_groups_pch<-factor(lda_groups_pch,c("HF-cet","LF-cet","Nasal","Oral","Non","non-bat","unknown"))
levels(lda_groups_pch)[levels(lda_groups_pch)%in%names(echo_pch_da)]<-echo_pch_da[match(levels(lda_groups_pch)[levels(lda_groups_pch)%in%names(echo_pch_da)],names(echo_pch_da))]
levels(lda_groups_pch)[levels(lda_groups_pch)=="unknown"]<-NA
lda_groups_pch<-as.numeric(as.character(lda_groups_pch))

pfda_groups_chu<-pfda_groups_pch<-setNames(as.character(avg_DA_scores_optL[,5]),rownames(avg_DA_scores_optL))
pfda_groups_chu[pfda_groups_chu%in%c("Nasal","Oral")]<-"LE"
pfda_groups_chu<-factor(pfda_groups_chu,c("HF-cet","LF-cet","LE","Non","non-bat","unknown"))

pfda_groups_pch<-factor(pfda_groups_pch,c("HF-cet","LF-cet","Nasal","Oral","Non","non-bat","unknown"))
levels(pfda_groups_pch)[levels(pfda_groups_pch)%in%names(echo_pch_da)]<-echo_pch_da[match(levels(pfda_groups_pch)[levels(pfda_groups_pch)%in%names(echo_pch_da)],names(echo_pch_da))]
levels(pfda_groups_pch)[levels(pfda_groups_pch)=="unknown"]<-NA
pfda_groups_pch<-as.numeric(as.character(pfda_groups_pch))

# Running the figures

LDA_pFDA_morphospaces<-plots.variation(list_data=list(pfda_null$DA.scores,pfda_optL$DA.scores),
                                       chull_data = if(exists("LDA_pFDA_morphospaces")){LDA_pFDA_morphospaces}else{NA},
                                       avg_data=c(list(avg_DA_scores_null),list(avg_DA_scores_optL)),
                                       data.name=c("DA","pFDA"),
                                       axes=c(1,2),
                                       axes.contrib=c(list(avg_pe_null),list(avg_pe_optL)),
                                       morpho.groups=c(list(lda_groups_chu),list(pfda_groups_chu)),
                                       morpho.cols=c(col_HF,col_LF,col_LE,col_nonLE,col_nonbat),
                                       morpho.pch=c(list(lda_groups_pch),list(pfda_groups_pch)),
                                       VIE.col = col_VIE_fig3,
                                       return.chull.data=TRUE,out=TRUE,single.out=1,
                                       out.names=c("Fig. 3D - Plot of mammal species including bats on first two axes of linear discriminant analysis (LDA) of PGLS regression residuals.eps",
                                                   "Fig. S3B - Individual variation in the discriminant analyses performed on the PGLS regression residuals.eps"))

############### III - GEOMETRIC MORPHOMETRICS OF CRANIA AND MANDIBLE AFTER ARBOUR ET AL. 2019
############### III - 1 - Input, prepare, and treat our perform the Procrustes superimposition following Arbour et al. 2019

# Crania landmarks and sliding semi-landmarks
lmk_crania<-read.table("8 lmk_crania.txt",header=T,dec=",",sep="\t")
crania_pairs<-rbind(cbind(c(seq(1,21,2)),c(seq(1,21,2))+1),matrix(ncol=2,nrow=length(c(38:53)),c(38:69),byrow=FALSE))
crania_slides<-list("1L"=c(27:37),"2L"=c(38:45),"2R"=c(54:61),"3L"=c(46:53),"3R"=c(62:69))
crania_gpa<-Arbour_et_al_treatment(lmk_crania[,3:dim(lmk_crania)[2]],lmk_crania$Shi.Tree.label,crania_pairs,crania_slides)

# Mandible landmarks and sliding semi-landmarks
lmk_mandible<-read.table("9 lmk_mandible.txt",header=T,dec=",",sep="\t")
mandible_pairs<-rbind(cbind(c(seq(2,12,2)),c(seq(2,12,2))+1),
                      matrix(ncol=2,nrow=length(c(14:21)),c(14:29),byrow=FALSE),
                      matrix(ncol=2,nrow=length(c(30:39)),c(c(30:39),c(50:41)),byrow=FALSE))
mandible_slides<-list("1L"=c(14:21),"1R"=c(22:29),"2L"=c(30:39),"2R"=c(50:41))
mandible_gpa<-Arbour_et_al_treatment(lmk_mandible[,3:dim(lmk_mandible)[2]],lmk_mandible$Shi.Tree.label,mandible_pairs,mandible_slides)

############### III - 2 - Perform pPCAs for both datasets on all phylogenies

cl<-makeCluster(detectCores()-1)
registerDoParallel(cl)
crania_ppca<-foreach(i=1:length(all_conc_SR),.packages = c("phytools","ape"))%dopar%{
  crania_keep<-rownames(crania_gpa)[rownames(crania_gpa)%in%all_conc_SR[[i]]$tip.label]
  phyl.pca(keep.tip(all_conc_SR[[i]],crania_keep),crania_gpa[crania_keep,],method="BM")$S
}
mandible_ppca<-foreach(i=1:length(all_conc_SR),.packages = c("phytools","ape"))%dopar%{
  mandible_keep<-rownames(mandible_gpa)[rownames(mandible_gpa)%in%all_conc_SR[[i]]$tip.label]
  phyl.pca(keep.tip(all_conc_SR[[i]],mandible_keep),mandible_gpa[mandible_keep,],method="BM")$S
}
stopCluster(cl)

############### III - 3 - Custom operations on the points resulting from the pPCAs

# 'Realign' points, so that there are no axes inversions throughout the replications (i.e., phylogenies)
aligned_crania_ppca<-align_rep_ppca(crania_ppca)
aligned_mandible_ppca<-align_rep_ppca(mandible_ppca)

# Get the % of variance explained on each pPCA axis following Arbour et al. 2019 (i.e., % of variance of an axis = R2 explained by the given axis while regressing morphometric data against that axis)
cl<-makeCluster(detectCores()-1)
registerDoParallel(cl)
crania_ppca_R2<-foreach(j=1:length(aligned_crania_ppca),.combine = "rbind",.packages="geomorph")%dopar%{
  procD.lm(crania_gpa~aligned_crania_ppca[[j]][,1]+
             aligned_crania_ppca[[j]][,2]+
             aligned_crania_ppca[[j]][,3]+
             aligned_crania_ppca[[j]][,4],
           data=list(crania_gpa,aligned_crania_ppca[[j]][,1:4]))$aov.table[1:4,4]
}
mandible_ppca_R2<-foreach(j=1:length(aligned_mandible_ppca),.combine = "rbind",.packages="geomorph")%dopar%{
  procD.lm(mandible_gpa~aligned_mandible_ppca[[j]][,1]+
             aligned_mandible_ppca[[j]][,2]+
             aligned_mandible_ppca[[j]][,3]+
             aligned_mandible_ppca[[j]][,4],
           data=list(mandible_gpa,aligned_mandible_ppca[[j]][,1:4]))$aov.table[1:4,4]
}
stopCluster(cl)

############### III - 4 - Get averages of each species position on each axis and of the % of variance explained by each axis

# Get the average value for each species on each pPCA axis. One can also ask for the median, the standard deviation, or the interquartile range by setting c("mean","sd","median","IQR").
avg_crania_ppca<-rep_ppca_operations(aligned_crania_ppca,"mean")$mean
avg_mandible_ppca<-rep_ppca_operations(aligned_mandible_ppca,"mean")$mean

# Get average % of variance explained
crania_ppca_R2<-round(100*apply(crania_ppca_R2,2,mean),1)
mandible_ppca_R2<-round(100*apply(mandible_ppca_R2,2,mean),1)

############### III - 5 - Get exhaustive information of each species position by computing their 
###############          'morphospaces' with 'full' and '95%' morphospaces and plot them

# Input ecological information (i.e., diet and echolocation call type) used by Arbour et al. 2019 and define colors
lmk_eco<-read.table("10 lmk_eco.txt",header=T,dec=",",sep="\t")

crania_species<-rownames(crania_gpa)
crania_eco<-lmk_eco[which(lmk_eco[,1]%in%crania_species),]
crania_eco<-crania_eco[match(rownames(avg_crania_ppca),crania_eco[,1]),]

mandible_species<-rownames(mandible_gpa)
mandible_eco<-lmk_eco[which(lmk_eco[,1]%in%mandible_species),]
mandible_eco<-mandible_eco[match(rownames(avg_mandible_ppca),mandible_eco[,1]),]

# For crania points
# NB: Graphics are displayed to resemble to those of Arbour et al. 2019.
#     Since on our crania pPCAs, pPC1 and pPC2 are inverted compared to the pPCA of Arbour et al. 2019, we plot opposite 
#     values for both axes (i.e., -pPC1 and -pPC2)
gm_groups_chu<-factor(setNames(crania_eco$Echolocation.Emission.Type,crania_eco$Species),c("Oral","Nasal","Non","unknown"))

echo_pch_gm<-setNames(c(pch_yang,pch_rhino,pch_ptero),c("Yang","Rhino","Ptero"))

gm_groups_pch<-crania_eco$Species
gm_groups_pch[gm_groups_pch%in%SR_phylo$tip.label[getDescendants(SR_phylo,819)[getDescendants(SR_phylo,819)<=Ntip(SR_phylo)]]]<-"Yang"
gm_groups_pch[gm_groups_pch%in%SR_phylo$tip.label[getDescendants(SR_phylo,1402)[getDescendants(SR_phylo,1402)<=Ntip(SR_phylo)]]]<-"Rhino"
gm_groups_pch[gm_groups_pch%in%SR_phylo$tip.label[getDescendants(SR_phylo,1522)[getDescendants(SR_phylo,1522)<=Ntip(SR_phylo)]]]<-"Ptero"
gm_groups_pch[!gm_groups_pch%in%c("Yang","Rhino","Ptero")]<-"unknown"
gm_groups_pch<-factor(gm_groups_pch,c("Yang","Rhino","Ptero","unknown"))
levels(gm_groups_pch)[levels(gm_groups_pch)%in%names(echo_pch_gm)]<-echo_pch_gm[match(levels(gm_groups_pch)[levels(gm_groups_pch)%in%names(echo_pch_gm)],names(echo_pch_gm))]
levels(gm_groups_pch)[levels(gm_groups_pch)=="unknown"]<-NA
gm_groups_pch<-as.numeric(as.character(gm_groups_pch))

crania_morphospaces<-plots.variation(list_data=aligned_crania_ppca,
                                       chull_data = if(exists("crania_morphospaces")){crania_morphospaces}else{NA},
                                     avg_data=avg_crania_ppca,
                                     data.name="pPC",
                                     axes=c(1,2,3),
                                     axes.contrib=crania_ppca_R2,
                                     axes.sign=c(-1,-1,1),
                                     morpho.groups=gm_groups_chu,
                                     morpho.cols=c(col_oral,col_nasal,col_nonLE),
                                     morpho.pch=gm_groups_pch,
                                     VIE.col = col_VIE,
                                     return.chull.data=TRUE,out=TRUE,single.out=1,
                                     out.names=c("Fig. 4A - Cranial shape in Vielasia sigei compared with extant bats.eps",
                                                 "Fig. S5A - Individual variation in pPCAs performed on cranium dataset.eps"))

# For mandible points
# NB: Graphics are displayed to resemble to those of Arbour et al. 2019.
#     Since on our mandible pPCAs, pPC3 and pPC4 are inverted compared to the pPCA of Arbour et al. 2019, we plot opposite 
#     values for both axes (i.e., -pPC3 and -pPC4)

mandible_morphospaces<-plots.variation(list_data=aligned_mandible_ppca,
                                       chull_data = if(exists("mandible_morphospaces")){mandible_morphospaces}else{NA},
                                       avg_data=avg_mandible_ppca,
                                       data.name="pPC",
                                       axes=c(1,2,3,4),
                                       axes.contrib=mandible_ppca_R2,
                                       axes.sign=c(1,1,-1,-1),
                                       morpho.groups=setNames(factor(mandible_eco$Diet),mandible_eco$Species),
                                       morpho.cols=c("#BA9234","#A11174","#007272","#77D377","#76A9CE","#8E2121"),
                                       VIE.col = col_VIE,
                                       return.chull.data=TRUE,out=TRUE,single.out=1,
                                       out.names=c("Fig. 4B - Mandible shape in Vielasia sigei compared with extant bats.eps",
                                                   "Fig. S5B - Individual variation in pPCAs performed on mandible dataset.eps"))

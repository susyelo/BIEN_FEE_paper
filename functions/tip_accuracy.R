#Taken from http://www.ecography.org/appendix/ecog-03480

tip_accuracy  <- function(Tree, Trait, runs = 500, threshold = 0.75, range.lambda = 0.05, method=c("PEM","Rphylopars"), verbose=TRUE) {
  
  require(spacodiR)
  require(Rphylopars)
  require(phytools)
  require(geiger)
  require(MPSEM)
  
  cophenetic(Tree)->Tree_dist
  Tree->Tree_Ori
  
  Trait_O <- Trait[is.na(Trait)==FALSE]
  Tree_signal <-drop.tip(Tree_Ori,Tree_Ori$tip.label[Tree_Ori$tip.label%in%names(Trait_O)==FALSE])
  Trait_O <-Trait_O[Tree_signal$tip.label]
  Trait_NA <- Trait[is.na(Trait)]
  
  phylosig(Tree_signal, Trait_O,"lambda")[[1]]->lambda
  
  if(length(Trait[is.na(Trait)])> length(Trait)/2) {print("caution, you have more missing data than observations!")} else
  {Sys.sleep(0.00000000001)}
  
  if(lambda>1){lambda<-1
  print("starting workflow, please wait")}else
  {if(lambda<0.9){print("caution, lambda statistic is less than 0.9, and may be too low for phylogenetic imputation")} else
  {print("starting workflow, please wait")}}
  
  Sys.sleep(2)
  
  LL <- 1:length(Tree$tip.label)
  names(LL)<-Tree$tip.label
  to_prun <-NULL
  checked <-NULL
  
  for(i in (length(Tree$tip.label)+1):((length(Tree$tip.label)+1)+(length(Tree$tip.label)-2))){
    
    LL[LL%in%get.descendants.of.node(i,Tree,tips=TRUE)]->LL_subset
    Trait[names(Trait)%in%names(LL_subset)]->Trait_Obs_subset
    
    if(sum(is.na(Trait_Obs_subset))<length(Trait_Obs_subset)){1} else
    {if(length(setdiff(names(Trait_Obs_subset),checked))==0){1}else
    {c(checked,names(Trait_Obs_subset))->checked
      c(to_prun,names(sample(Trait_Obs_subset,length(Trait_Obs_subset)-1)))->to_prun}
    }}
  
  if(length(to_prun)>0) {Tree <- drop.tip(Tree,to_prun)} else
  {1}
  
  Res_PEM<-vector(mode="numeric", length= length(Tree$tip.label))  
  Store_tip_PEM <- rep(1,length(Tree$tip.label)) 
  
  if(method == "PEM")
  {
    for(z in 1:runs){
      
      if(verbose==TRUE) print(paste("searching trait ",z,sep=""))
      
      for(r in 1:1000000){
        x<-fastBM(rescale(Tree,"lambda", lambda),a=0, sig2=1) 
        if(phylosig(Tree,x,"lambda")[[1]]>(lambda-range.lambda) & phylosig(Tree,x,"lambda")[[1]]<(lambda+range.lambda)) {break} else
        {if(r > 999999) {stop("failure to obtain trait of desired phylogenetic signal")} else
        {Sys.sleep(0.00000000001)}}
      }
      
      if(verbose==TRUE) print(paste("imputing trait ",z,sep=""))
      
      trait<-data.frame(as.factor(Tree$tip.label),x[names(x)%in%Tree$tip.label])
      colnames(trait)<-c("species","V1")
      for(i in 1:length(trait[,2])){
        if(trait[i,1]%in%names(Trait_NA)) {trait[i,2]<-NA} else
        {1}
      }
      x_prun <- trait[,2]
      names(x_prun) <- trait[,1]
      x_prun <- x_prun[is.na(x_prun)==FALSE]
      tree_prun <- drop.tip(Tree,Tree$tip.label[Tree$tip.label%in%names(x_prun)==FALSE])
      Ymodel <- matrix(0,length(tree_prun$tip.label),1L,dimnames=list(tree_prun$tip.label,"trait")) 
      grModel <- Phylo2DirectedGraph(tree_prun) 
      Ymodel[,1] <- x_prun 
      PEMAG_select<-PEM.fitSimple(Ymodel[,1],x=NULL,w=grModel,d = "distance", sp = "species",lower = 0, upper = 1) 
      lm_PEM<-lmforwardsequentialAICc(y=Ymodel[,1],object=PEMAG_select) 
      targ_ite <- trait[,2]
      names(targ_ite) <- trait[,1]
      targ<-names(targ_ite[is.na(targ_ite)==TRUE])
      Tree_targ <- Tree
      targ.loc <- getGraphLocations(Tree_targ,targets=targ)
      newdata<-subset(trait, is.na(trait[,2]))  
      newdata[,2]->newdata2					  
      names(newdata2)<-row.names(newdata)		  
      newdata<-as.matrix(newdata2)			  
      rm(newdata2)								
      colnames(newdata)[1]<-"trait"				
      pred <- predict(object=PEMAG_select, targets=targ.loc, lmobject=lm_PEM, newdata=newdata, "prediction",0.95)
      for(i in 1:length(Tree$tip.label)){
        if(Tree$tip.label[i]%in%row.names(pred$values)) {Res_PEM[i]<-pred$values[row.names(pred$values)%in%Tree$tip.label[i]]} else
        {Res_PEM[i]<-999999999999999999}
      }
      names(Res_PEM)<-row.names(trait)
      Res_PEM_NA <- Res_PEM[Res_PEM==999999999999999999]
      if(verbose==TRUE) {print(paste("Trait ",z," imputed",sep=""))}
      
      R_v_tip_PEM <- vector(mode="numeric", length=length(Tree$tip.label))  
      
      for (j in 1:length(Tree$tip.label)){
        if(Tree$tip.label[j]%in%names(Res_PEM_NA)){R_v_tip_PEM[j]<-999999999999999999} else
        {R_v_tip_PEM[j] <- 1-((sum((x[j]-Res_PEM[j])^2))/(var(x)))}
      }
      Store_tip_PEM <- data.frame(Store_tip_PEM,R_v_tip_PEM)
    }
    Store_tip_PEM <- Store_tip_PEM[,-1]
    row.names(Store_tip_PEM)<-names(Trait)
  }
  
  if(method == "Rphylopars") 
  {
    for(z in 1:runs){
      
      if(verbose==TRUE) print(paste("searching trait ",z,sep=""))
      
      for(r in 1:1000000){
        x<-fastBM(rescale(Tree,"lambda", lambda),a=0, sig2=1) 
        if(phylosig(Tree,x,"lambda")[[1]]>(lambda-range.lambda) & phylosig(Tree,x,"lambda")[[1]]<(lambda+range.lambda)) {break} else
        {if(r > 999999) {stop("failure to obtain trait of desired phylogenetic signal")} else
        {Sys.sleep(0.00000000001)}}
      }
      
      if(verbose==TRUE) print(paste("imputing trait ",z,sep=""))
      
      trait<-data.frame(as.factor(Tree$tip.label),x[names(x)%in%Tree$tip.label],x[names(x)%in%Tree$tip.label])
      colnames(trait)<-c("species","V1","V2")
      for(i in 1:length(trait[,2])){
        if(trait[i,1]%in%names(Trait_NA)) {trait[i,2]<-NA} else
        {1}
      }
      PE<-phylopars(trait, Tree, model = "lambda", pheno_error=FALSE,pheno_correlated = FALSE,phylo_correlated=FALSE) 
      PE_predict <- PE$anc_recon
      Predichos_pars <- PE_predict[c(1:length(Tree$tip.label)),1]
      
      Trait_NA_pars <- Trait_NA[names(Trait_NA)%in%names(Predichos_pars)]
      
      for(i in 1:length(Tree$tip.label)){
        if(Tree$tip.label[i]%in%names(Trait_NA_pars)) {Res_PEM[i]<-Predichos_pars[names(Predichos_pars)%in%Tree$tip.label[i]]} else
        {Res_PEM[i]<-999999999999999999}
      }
      names(Res_PEM)<-row.names(trait)
      
      Res_PEM_NA <- Res_PEM[Res_PEM==999999999999999999]
      if(verbose==TRUE) {print(paste("Trait ",z," imputed",sep=""))}
      
      R_v_tip_PEM <- vector(mode="numeric", length=length(Tree$tip.label))  
      
      for (j in 1:length(Tree$tip.label)){
        if(Tree$tip.label[j]%in%names(Res_PEM_NA)){R_v_tip_PEM[j]<-999999999999999999} else
        {R_v_tip_PEM[j] <- 1-((sum((x[j]-Res_PEM[j])^2))/(var(x)))}
      }
      Store_tip_PEM <- data.frame(Store_tip_PEM,R_v_tip_PEM)
    }
    Store_tip_PEM <- Store_tip_PEM[,-1]
    row.names(Store_tip_PEM)<-names(Trait)
  }
  
  P_Umbral<-vector(mode="numeric", length=length(Tree$tip.label))
  names(P_Umbral)<-Tree$tip.label
  
  for(f in 1:length(Tree$tip.label)){
    
    if(min(Store_tip_PEM[f,])>threshold | min(Store_tip_PEM[f,])==threshold) {P_Umbral[f]<-1}
    else {
      if(max(Store_tip_PEM[f,])<threshold | max(Store_tip_PEM[f,])==threshold) {P_Umbral[f]<-0}
      else {
        table(Store_tip_PEM[f,]>threshold | Store_tip_PEM[f,]==threshold)[[2]]/runs->P_Umbral[f]
      }
    }
  }
  
  Tree_dist <- Tree_dist[row.names(Tree_dist)%in%names(Trait_NA),colnames(Tree_dist)%in%names(Trait_NA)]
  P_Umbral_F <- NULL
  for(u in 1:length(Trait_NA)){
    if(names(Trait_NA)[u]%in%names(P_Umbral)){P_Umbral_F<-c(P_Umbral_F,P_Umbral[names(P_Umbral)%in%names(Trait_NA)[u]])} else
    {sort(Tree_dist[u,-u])->s_sort
      for(e in 1:length(s_sort)){
        if(names(s_sort)[e]%in%names(P_Umbral)){P_Umbral[names(P_Umbral)%in%names(s_sort)[e]]->Value
          names(Value)<-names(Tree_dist[u,])[u]
          P_Umbral_F<-c(P_Umbral_F,Value)
          break()} else
          {Sys.sleep(0.00000000001)}
      }}}
  
  P_Umbral_F<-P_Umbral_F[names(Trait_NA)]
  Store_tip_PEM_F <- Store_tip_PEM[names(Trait_NA),] 
  colnames(Store_tip_PEM_F) <- paste0("run", 1:runs)
  Store_tip_PEM_F <- t(Store_tip_PEM_F)
  
  if(method == "PEM") {
    
    trait<-data.frame(as.factor(Tree_Ori$tip.label),Trait)   
    colnames(trait)<-c("species","V1")
    x_prun <- Trait[names(Trait)%in%as.character(trait[,1][which(trait[,2]%in%NA)])==FALSE] 
    tree_prun <- drop.tip(Tree_Ori,as.character(trait[,1][which(trait[,2]%in%NA)]))
    Ymodel <- matrix(0,length(tree_prun$tip.label),1L,dimnames=list(tree_prun$tip.label,"trait")) 
    grModel <- Phylo2DirectedGraph(tree_prun) 
    Ymodel[,1] <- x_prun 
    PEMAG_select<-PEM.fitSimple(Ymodel[,1],x=NULL,w=grModel,d = "distance", sp = "species",lower = 0, upper = 1) 
    lm_PEM<-lmforwardsequentialAICc(y=Ymodel[,1],object=PEMAG_select) 
    targ<-names(Trait)[names(Trait)%in%names(x_prun)==FALSE]
    targ.loc <- getGraphLocations(Tree_Ori,targets=targ)
    newdata<-subset(trait, is.na(trait[,2]))  
    newdata[,2]->newdata2					  
    names(newdata2)<-row.names(newdata)		  
    newdata<-as.matrix(newdata2)			  
    rm(newdata2)								
    colnames(newdata)[1]<-"trait"				
    pred <- predict(object=PEMAG_select, targets=targ.loc, lmobject=lm_PEM, newdata=newdata, "prediction",0.95)
    Predictions2<-pred[[1]] 
    Predictions<-as.numeric(Predictions2[,1])
    names(Predictions)<-row.names(Predictions2)
  }
  
  if(method == "Rphylopars") {
    
    trait<-data.frame(as.factor(Tree_Ori$tip.label),Trait,Trait)
    colnames(trait)<-c("species","V1","V2")
    PE<-phylopars(trait, Tree_Ori, model = "lambda", pheno_error=FALSE,pheno_correlated = FALSE,phylo_correlated=FALSE) 
    PE_predict <- PE$anc_recon
    Predichos_pars <- PE_predict[c(1:length(Tree_Ori$tip.label)),1] 
    Predictions<-Predichos_pars[names(Predichos_pars)%in%names(Trait[is.na(Trait)==TRUE])] 
    
  } 
  
  return(list(lambda=lambda,tip_accuracy=P_Umbral_F,predictions=Predictions,accuracy_val=Store_tip_PEM_F)) 
  
}

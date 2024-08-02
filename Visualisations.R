# RefinedDAG ----

library(sharp)
library(igraph)
library(scales)
# source("plot_functions.R")

biomarkers = readRDS("../Data/Biomarkers.rds")
mytable = read.table("annot_biomarkers.txt")

mylabels = as.character(mytable[,2])
names(mylabels) = mytable[,1]

annot = as.character(mytable[,3])
names(annot) = mytable[,1]

panel.colour = c("#006FE6", "#ffa600", "#18a95d", "#9b0058", "#dc241f","#61686b", "#da4290")
names(panel.colour) = c("Growth","Metabolism","Inflammation", "Liver", "Kidney","Bone", "Reproductive system")

expo = readRDS("../Data/Exposures.rds")
colnames(expo) = c("Smoking", "Diet","PA","Deprivation","BP med","LL med","DM med","CS med")
expo_colours = c("#dc251f", "#0019a8", "#0019a8", "#dc251f", rep("#00722a",4))
names(expo_colours) = colnames(expo) 

asp = 9/16
width = 6
height = 4

# asp = 9/16
# width = 6
# height = 3.5

## Full DAG

for (instance in 0:3){
  print(paste0("Instance ", instance))
  for (cat in c("Female", "Male")){
    print(paste0("Sex: ",cat))
    
    ### Select biomarkers
    stab = readRDS(paste0("../Results/1-stab_",cat,"_Instance",instance,".rds"))
    selected = SelectedVariables(stab)
    selected = selected[grepl(paste0(colnames(biomarkers), collapse = "|"), names(selected))]
    X2 = biomarkers[,selected==1]
    annotsub = annot[colnames(X2)]
    
    # Create DAG
    pk = c(ncol(expo),ncol(X2),1)
    dag =  LayeredDAG(layers = pk)
    
    dag_refined = betas = selprop = posprop = dag
    for(k in 1:ncol(X2)){
      tmp = readRDS(paste0("../Results/2-stab_",colnames(X2)[k],"_",cat,"_Instance",instance,".rds")) # Layer 2 ~ Layer 1
      dag_refined[1:ncol(expo),ncol(expo)+k] = ifelse(SelectedVariables(tmp)==1, 1, 0)
      selprop[1:ncol(expo),ncol(expo)+k] = ifelse(SelectedVariables(tmp)==1, SelectionProportions(tmp), 0)
      
      # Computing average non-zero beta coefficient from models with calibrated lambda
      nconf = ncol(tmp$Beta)-length(SelectedVariables(tmp))
      beta_mu = apply(tmp$Beta[ArgmaxId(tmp)[1],-c(1:nconf),], 1, function(x) mean(x[x!=0]))
      betas[1:ncol(expo),ncol(expo)+k] = ifelse(SelectedVariables(tmp)==1, beta_mu, 0)
      
      a = apply(tmp$Beta[ArgmaxId(tmp)[1],-c(1:nconf),], 1, function(x) sum(x[x!=0]>0))
      b = apply(tmp$Beta[ArgmaxId(tmp)[1],-c(1:nconf),], 1, function(x) sum(x[x!=0]<0))
      posprop[1:ncol(expo),ncol(expo)+k] = ifelse(SelectedVariables(tmp)==1, a/(a+b), 0)
    }
    keep = c(1:ncol(expo),which(grepl(paste0(colnames(X2),collapse = "|"),names(SelectedVariables(stab)))))
    selected_all = SelectedVariables(stab)[keep]
    dag_refined[-nrow(dag),ncol(dag)] = selected_all
    selprop[-nrow(dag),ncol(dag)] = ifelse(selected_all==1, SelectionProportions(stab)[keep], 0)
    
    # Computing average non-zero beta coefficient from models with calibrated lambda
    nconf = ncol(stab$Beta)-length(SelectedVariables(stab))
    beta_mu = apply(stab$Beta[ArgmaxId(stab)[1],-c(1:nconf),], 1, function(x) mean(x[x!=0]))[keep]
    betas[-nrow(dag),ncol(dag)] = ifelse(selected_all==1, beta_mu, 0)
    
    a = apply(stab$Beta[ArgmaxId(stab)[1],-c(1:nconf),], 1, function(x) sum(x[x!=0]>0))[keep]
    b = apply(stab$Beta[ArgmaxId(stab)[1],-c(1:nconf),], 1, function(x) sum(x[x!=0]<0))[keep]
    posprop[-nrow(dag),ncol(dag)] = ifelse(selected_all==1,  a/(a+b), 0)
    
    ### Legend info
    panel.count = table(annotsub)
    
    pdf(paste0("../Figures/DAG_",cat,"_Instance",instance,".pdf"), width = width, height = height)
    par(mar = c(0,0,1,0))
    DrawDAG(
      adjacency = dag_refined, layers = c(rep(1,ncol(expo)), rep(2,sum(selected==1)),2*1.5),
      prop = posprop,
      edge.col = c("#dc251f","#0019a8","#000000"),
      node.col = c(alpha(expo_colours[colnames(expo)],0.5),alpha(panel.colour[annotsub],0.5),NA),
      node.size = rep(c(15,20,1),pk))
    legend("bottomright", legend = paste0(names(panel.count), " (N=", panel.count, ")"),
           bty = "n",
           pch = rep(19,length(panel.count)),
           col = panel.colour[names(panel.count)], cex = 0.6)
    mtext(paste0(cat, "s at Instance ", instance), side = 3, line = 0, adj = 0, font = 2)
    dev.off() 
    
    pdf(paste0("../Figures/DAG_",cat,"_Instance",instance,"_Direct.pdf"), width = width, height = height)
    par(mar = c(0,0,1,0))
    DrawDAG(
      adjacency = dag_refined, layers = c(rep(1,ncol(expo)), rep(2,sum(selected==1)),2*1.5),
      prop = posprop,
      edge.col = c("#dc251f","#0019a8","#000000"),
      node.col = c(alpha(expo_colours[colnames(expo)],0.5),alpha(panel.colour[annotsub],0.5),NA),
      node.size = rep(c(15,20,1),pk),
      direct = TRUE)
    legend("bottomright", legend = paste0(names(panel.count), " (N=", panel.count, ")"),
           bty = "n",
           pch = rep(19,length(panel.count)),
           col = panel.colour[names(panel.count)], cex = 0.6)
    mtext(paste0(cat, "s at Instance ", instance, " (Direct)"), side = 3, line = 0, adj = 0, font = 2)
    dev.off() 
  }
}

## SEM Parameters ----

library(sharp)

### Load data
biomarkers = readRDS("../Data/Biomarkers.rds")
for(instance in 0:3){
  print(paste0("Instance ",instance))
  for(cat in c("Female", "Male")){
    print(paste0("Sex: ",cat))
    stab = readRDS(paste0("../Results/1-stab_",cat,"_Instance",instance,".rds"))
    
    selected = SelectedVariables(stab)
    selected = selected[grepl(paste0(colnames(biomarkers), collapse = "|"),names(selected))]

    print(sum(selected==1))
  }
}

## Number of participants included ----
### Load data
covars = readRDS("../Data/Baseline.rds")
genetics = readRDS("../Data/Genetic_principal_components.rds")
outcome = readRDS("../Data/BMI_repeated.rds")
biomarkers = readRDS("../Data/Biomarkers.rds")
expo = readRDS("../Data/Exposures.rds")

### Subset data
common_ids <- Reduce(intersect, list(rownames(covars), rownames(biomarkers), rownames(expo)))

covars = covars[rownames(covars) %in% common_ids,]
genetics = genetics[rownames(genetics) %in% common_ids,]
outcome = outcome[rownames(outcome) %in% common_ids,]
biomarkers = biomarkers[rownames(biomarkers) %in% common_ids,]
expo = expo[rownames(expo) %in% common_ids,]

for(instance in 0:3){
  cat(paste0("Instance ",instance, "\n"))
  Y = outcome[, paste0("BMI",instance), drop=FALSE]
  
  tmp = cbind(covars, genetics, biomarkers, expo, Y)
  
  covars_sub = covars[complete.cases(tmp),]
  
  Y = Y[complete.cases(tmp),, drop=FALSE]
  for(cat in c("Female", "Male")){
    cat(paste0(cat, " "))
    Y = Y[covars_sub$Sex==cat,, drop=FALSE]
    cat(paste0("N=", nrow(Y), "\n"))
  }
}






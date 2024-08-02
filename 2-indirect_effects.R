# Structural Equation Model: Selecting indirect effects (exposure-biomarker) in DAG (Step 2)

library(sharp)

ifelse(dir.exists("../Results"), "", dir.create("../Results"))
ifelse(dir.exists("../Figures"), "", dir.create("../Figures"))

### Reading arguments
n_cores = 8

args = commandArgs(trailingOnly = TRUE)
k = as.numeric(args[1])
sex_id = as.numeric(args[2])
cat = c("Female", "Male")[sex_id + 1]
print(paste0("Sex: ",cat))

instance = as.numeric(args[3])

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

tmp = cbind(covars, genetics, biomarkers, expo, outcome)

covars = covars[complete.cases(tmp),]
genetics = genetics[complete.cases(tmp),]
biomarkers = biomarkers[complete.cases(tmp),]
expo = expo[complete.cases(tmp),]
outcome = outcome[complete.cases(tmp), ]

X = cbind(expo, biomarkers)[covars$Sex==cat,]

Age = covars$`Age-when-attended-assessment-centre`
G = genetics[,paste0("X22009.0.", 1:40)]
colnames(G) = paste0("Genetic PC", 1:40)
Time = outcome[, paste0("Time",instance), drop=FALSE]
Z = cbind(Age, G, Time)[covars$Sex==cat,]

# If adjusting for Ethnicity:
# keep = colnames(covars) %in% c("Age-when-attended-assessment-centre", "Ethnicity")
# Z = model.matrix(~ ., Z)[,-1]

## Select biomarkers
stab = readRDS(paste0("../Results/1-stab_",cat,"_Instance",instance,".rds"))

tmp = SelectedVariables(stab)
selected = tmp[grepl(paste0(colnames(biomarkers), collapse = "|"),names(tmp))]

Y = X[,colnames(biomarkers)]
Y = Y[,selected==1]
X = X[,colnames(expo)]

print(colnames(Y)[k])

# If Error in glmnet::glmnet(x = xdata[s, ], y = ydata[s, ], family = family, : the length of penalty.factor does not match the number of variables (run below)
# X = model.matrix(~ ., X)[,-1]

# If Error in CheckDataRegression(xdata = xdata, ydata = ydata, family = family,  : Arguments 'xdata' and 'ydata' are not compatible. They have different numbers of observations.
# apply(cbind(Z,X), 2, function(x) sum(is.na(x)))

t0=Sys.time()

# LASSO linear regression
stab = VariableSelection(cbind(Z,X), Y[,k], family = "gaussian",
                         n_cat=3, pi_list = seq(0.5, 0.99, by = 0.01),
                         K = 100, tau = 0.5,
                         penalty.factor = c(rep(0,ncol(Z)),rep(1,ncol(X))),
                         n_cores = n_cores)

pdf(paste0("../Figures/2-CalibrationPlot_",colnames(Y)[k],"_",cat,"_Instance",instance,".pdf"))
par(mar = c(7,5,7,6))
CalibrationPlot(stab)
dev.off()

saveRDS(stab, paste0("../Results/2-stab_",colnames(Y)[k],"_",cat,"_Instance",instance,".rds"))

t1=Sys.time()
print(t1-t0)  


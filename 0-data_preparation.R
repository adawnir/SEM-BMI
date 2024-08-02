# Data cleaning: Exposures ~ Biomarkers ~ BMI (Instance 0, 1, 2, and 3)

### Reading arguments
args <- commandArgs(trailingOnly = TRUE)
baseline_path <- as.character(args[1])
withdrawn <- as.character(args[2])
withdrawn = as.character(read.csv(withdrawn)[,1])
ukb_path <- as.character(args[3])
repeat_path <- as.character(args[4])
townsend_path <- as.character(args[5])
score_path <- as.character(args[6])

ifelse(dir.exists("../Data"), "", dir.create("../Data"))

library(dplyr)
library(data.table)

# Baseline ----
covars = readRDS(paste0(baseline_path, "Data/Baseline.rds"))

### Exclude participants who have withdrawn from study
print(paste0(sum(rownames(covars) %in% withdrawn), " participants withdrawn"))
covars = subset(covars[!rownames(covars) %in% withdrawn,])

## Ethnic background ----
input = readRDS(paste0(baseline_path, "Processed/Sociodemographics.rds"))
input = input[rownames(covars),]

ethnic_mapping = c("White" = 1,
                   "British" = 1,
                   "Irish" = 1,
                   "Any other white background" = 1,
                   "Prefer not to answer" = 1,
                   "Do not know" = 1,
                   "Indian" = 2,
                   "Pakistani" = 2,
                   "Bangladeshi" = 2,
                   "Asian or Asian British" = 2,
                   "Chinese" = 2,
                   "Any other Asian background" = 2,
                   "Black or Black British" = 3,
                   "Caribbean" = 3,
                   "African" = 3,
                   "Any other Black background" = 3,
                   "Mixed" = 4,
                   "White and Black Caribbean" = 4,
                   "White and Black African" = 4,
                   "White and Asian" = 4,
                   "Any other mixed background" = 4,
                   "Other ethnic group" = 5)

Ethnicity = factor(ethnic_mapping[as.character(input$`Ethnic-background`)],
                   levels = 1:5, labels = c("White", "Asian", "Black","Mixed", "Other"))

covars = cbind(covars, Ethnicity)
covars = covars[!is.na(covars$Ethnicity),]

saveRDS(covars, "../Data/Baseline.rds")

# Genetic PCs (Population stratification) ----

## Select fields based on category
myfields = paste0("220",c("09","04","20")) # Genetic principal components

# Extract dataset
mydata = data.frame(fread(ukb_path, nrows=1))

# Extracting the column ids 
column_id = grep("eid", colnames(mydata))
for (k in 1:length(myfields)){
  mygrep=grep(paste0("X",myfields[k],"."), fixed=TRUE, colnames(mydata))
  column_id=c(column_id, mygrep)
}

# Extracting required columns from dataset
extracted = data.frame(fread(ukb_path, select=column_id))
rownames(extracted) = extracted$eid
genotype = extracted[rownames(extracted) %in% rownames(covars),-1]

saveRDS(genotype, "../Data/Genetic_principal_components.rds")

# Outcomes (Repeated measurements) ----

outcome = NULL
for(j in 0:3){
  print(paste0("Instance ", j))
  baseline = readRDS("../Data/Baseline.rds")
  
  reception = readRDS(paste0(repeat_path, "Data/Reception_Instance", j,".rds"))
  reception = reception[rownames(baseline),]
  all(rownames(baseline) == rownames(reception))
  colnames(reception) = gsub("\\.[0-3]\\.0", "", colnames(reception))
  
  body = readRDS(paste0(repeat_path, "Processed/Body_measurements_Instance", j,".rds"))
  body = body[rownames(baseline),]
  
  ## Proportion of missing samples ----
  X = body$`Body-mass-index-BMI`
  naprop = sum(is.na(X), na.rm = T)/length(X)
  print(paste0("Proportion of missing samples: ", round(naprop,2)))
  print(paste0(naprop*length(X), " participants missing BMI"))

  
  ### Time and ID
  Time = reception$`Age-when-attended-assessment-centre` - baseline$`Age-when-attended-assessment-centre`
  
  outcome = cbind(outcome, X, Time)
}
outcome = as.data.frame(outcome)
rownames(outcome) = rownames(baseline)
colnames(outcome) = c(t(cbind(paste0("BMI", 0:3), paste0("Time", 0:3))))

saveRDS(outcome, "../Data/BMI_repeated.rds")

# Biomarker data ----
biomarkers = readRDS(paste0(baseline_path, "/Data/Blood_biochemistry.rds"))

## Proportion of missing samples per biomarker ----
X = biomarkers
naprop = apply(X, 2, function(x) sum(is.na(x), na.rm = T)/nrow(X))
prop = sort(naprop, decreasing = TRUE)
names(prop) = gsub("-"," ", names(prop))

ifelse(dir.exists("../Figures"), "", dir.create("../Figures"))

pdf("../Figures/Proportion_missing_Blood_biochemistry.pdf", width = 7, height = 6)
par(mar=c(16,5,2,1), pty = "m")
plot(prop, type = "h", ylab="Proportion missing", xlab="", xaxt = "n", lwd=1, ylim = c(0,1),
     col=ifelse(prop>=0.5, yes="tomato", no="black"), cex.lab=1.5)
abline(h=0.5, lty=2, col="darkred", lwd = 1)
for (i in 1:length(prop)){
  axis(side=1, at=i, labels=names(prop)[i], las=2,
       col=ifelse(prop[i]>=0.5, yes="tomato", no="black"),
       col.axis=ifelse(prop[i]>=0.5, yes="tomato", no="black"))
}
dev.off()

X = X[,naprop<0.5] # Exclude Rheumatoid factor and Oestradiol

## Proportion of missing biomarkers per participant----
naprop = apply(X, 1, function(x) sum(is.na(x), na.rm = T)/ncol(X))
summary(naprop)

print(paste0(sum(naprop>=0.5), " participants with 50% or more missing"))
X = X[naprop<0.5,]

## Imputation of missing values ---
X_imp = X
imp = impute::impute.knn(as.matrix(X[,sapply(X, class)=="numeric"]), k = 10, rowmax = 0.5, colmax = 0.5, maxp = 1500, rng.seed=362436069)
X_imp[,sapply(X, class)=="numeric"] = imp$data

rownames(X_imp) = rownames(X)
colnames(X_imp) = colnames(X)

saveRDS(X_imp, "../Data/Blood_biochemistry_imputed.rds")


print(paste0(sum(rownames(X_imp) %in% rownames(covars)), " participants with data"))
biomarkers = X_imp[rownames(X_imp) %in% rownames(covars),]

saveRDS(biomarkers, "../Data/Biomarkers.rds")

# Exposures ----

### Load data
lifestyle = readRDS(paste0(baseline_path, "Processed/Lifestyle.rds"))
medication = readRDS(paste0(baseline_path, "Processed/Medications_supplements.rds"))

all(rownames(lifestyle)==rownames(medication))

mydata = cbind(lifestyle, medication)[rownames(covars),]

### Townsend Deprivation Index ----
tmp = readRDS(townsend_path)
townsend = tmp$`189-0.0`
names(townsend) = tmp$eid

mydata$`Townsend-deprivation-index` = townsend[rownames(mydata)]

### DASH ----
dash = readRDS(paste0(score_path, "DASH_ukb673609.rds"))
dash = dash[rownames(mydata),]
mydata$DASH = dash$dash

### MET ----
met = readRDS(paste0(score_path, "MET_ukb673609.rds"))
mydata$MET = met[rownames(mydata)]

## Select exposures ----
tmp = mydata %>% select(`Pack-years-adult-smoking-as-proportion-of-life-span-exposed-to-smoking`,
                          DASH, MET, `Townsend-deprivation-index`,
                          `Blood pressure medication`, `Lipid lowering medication`,
                          `Diabetes mellitus medication`, Corticosteroids)

# Count rows with at least one "Prefer not to answer"
count <- apply(tmp, 1, function(row) any(row == "Prefer not to answer", na.rm = T))
# Summarize the count
print(paste0(sum(count), " participants responded 'prefer not to answer' to variable(s) of interest"))

# Filter out rows containing "Prefer not to answer"
mydata = mydata[!count,]  %>%  mutate_if(is.factor, droplevels)

### Smoking ----
current = mydata$`Current-tobacco-smoking`
pkyrs = mydata$`Pack-years-adult-smoking-as-proportion-of-life-span-exposed-to-smoking`
mydata$`Pack-years-adult-smoking-as-proportion-of-life-span-exposed-to-smoking` = ifelse(!is.na(current) & is.na(pkyrs), 0, pkyrs)

## Proportion of missing samples per feature ----
X = mydata %>% select(`Pack-years-adult-smoking-as-proportion-of-life-span-exposed-to-smoking`,
                      DASH, MET, `Townsend-deprivation-index`,
                      `Blood pressure medication`, `Lipid lowering medication`,
                      `Diabetes mellitus medication`, Corticosteroids)
colnames(X) = c("Pack years","Diet score", "Physical activity score",
                "Townsend deprivation index", "Blood pressure medication",
                "Lipid-lowering medication",
                "Diabetes mellitus medication", "Corticosteroids")
  
naprop = apply(X, 2, function(x) sum(is.na(x), na.rm = T)/nrow(X))
prop = sort(naprop, decreasing = TRUE)
print(prop)

saveRDS(X, "../Data/Exposures.rds")

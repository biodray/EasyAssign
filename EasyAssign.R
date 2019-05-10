# Info --------------------------------------------------------------------

# Perform assignment with assignPOP
# 
# Audrey Bourret
# 2019-05-10
#

# See https://alexkychen.github.io/assignPOP/index.html

# Library -----------------------------------------------------------------

# POPassign - ** Attention R >= 3.6.0 **
if(!require(assignPOP)){ install.packages("assignPOP") }
library(assignPOP)

if(!require(klaR)){ install.packages("klaR") }
library(klaR)

# readxl convert read gen data
if(!require(readxl)){ install.packages("readxl") }
library(readxl)

# readxl to play around xl filesconvert read gen data
if(!require(readxl)){ install.packages("readxl") }
library(readxl)

#tidyverse to play around xl filesconvert read gen data
if(!require(tidyverse)){ install.packages("tidyverse") }
library(tidyverse)

# Internal functions
for(i in 1:length( list.files("./04_Functions") )){
  source(file.path("./04_Functions",  list.files("./04_Functions")[i]))  
}


# Parameters --------------------------------------------------------------

# Set these parameters

locus <- c("SEB25", "SEB31", "SEB33", "SEB9")

ref.excel <- "Microsat_mentella_fasciatus_ Mux1_03042019.xls"

assign.excel <- "BIN-SEB_Maria_Senay_Mux1_tot_2019.xlsx"

# Nothing to change here
ref.file  <- ref.excel %>% str_replace(".xlsx", ".gen") %>% str_replace(".xls", ".gen")
ref.dir   <- file.path("01_Ref_Genotypes" ,ref.excel %>% str_remove(".xlsx") %>% str_remove(".xls"))

assign.file  <- assign.excel %>% str_replace(".xlsx", ".gen") %>% str_replace(".xls", ".gen")
assign.dir   <- file.path("03_Results" ,assign.excel %>% str_remove(".xlsx") %>% str_remove(".xls"))

# Prepare Genetic Data ----------------------------------------------------

# Step 1 : Excel to Dataframe format
ref.geno <- as.data.frame(read_excel(file.path("01_Ref_Genotypes",ref.excel),  col_types = rep("text", 2 + length(locus) * 2)))
ref.df <- merge.MSAT.alleles(data=ref.geno, locus=locus, na="NA")
ref.df


# Step 2 : Merge alleles
assign.geno <- as.data.frame(read_excel(file.path("02_Data_to_Assign",assign.excel),  col_types = rep("text", 1 + length(locus) * 2)))
assign.df <- merge.MSAT.alleles(data=assign.geno, locus=locus, na="NA")
assign.df

# Step 3: Save as genpop
write.genpop(fn = file.path("01_Ref_Genotypes",ref.file), 
             data = ref.df, 
             pop = "POP", 
             ind = "ID",
             locus = locus)

write.genpop(fn = file.path("02_Data_to_Assign",assign.file), 
             data = assign.df, 
             ind = "ID",
             locus = locus)

# Step 4: Upload in the rigth format

ref.gen    <- read.Genepop( file.path("01_Ref_Genotypes",ref.file), pop.names=unique(ref.df$SP), haploid = FALSE)
assign.gen <- read.Genepop( file.path("02_Data_to_Assign",assign.file), pop.names="POP1", haploid = FALSE)

# Evaluate baseline -------------------------------------------------------

# remove low variance loci

ref.gen.rd <- reduce.allele(ref.gen, p = 0.95)

# NOTE: it is not necessary to run this part each time ...

# Compute cross-validation statistics

# Make a directory
dir.create(ref.dir)

assign.MC(ref.gen.rd, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.5, 1),
           loci.sample="fst", iterations=30, model="svm", dir=paste0(file.path(ref.dir,"MC_cross-validation"),"/"))

assign.kfold(ref.gen.rd, k.fold=c(3, 4, 5), train.loci=c(0.5, 1), 
             loci.sample="random", model="lda", dir=paste0(file.path(ref.dir,"kfold_cross-validation"),"/"))

accuMC <- accuracy.MC(dir = paste0(file.path(ref.dir,"MC_cross-validation"),"/"))

accuracy.plot(accuMC, pop = "all")


accuKF <- accuracy.kfold(dir = paste0(file.path(ref.dir,"kfold_cross-validation"),"/"))

accuracy.plot(accuKF, pop = "all")


# Predict sources of unknown individuals ----------------------------------

# Make a directory
dir.create(assign.dir)

# 1.Perform assignment test using genetic data and naive Bayes
assign.X( x1=ref.gen.rd, x2=assign.gen, dir=paste0(file.path(assign.dir,"naiveBayes"),"/"), model="naiveBayes")

# 2.Perform assignment test using integrated data and decision tree
assign.X( x1=ref.gen.rd, x2=assign.gen, dir=paste0(file.path(assign.dir,"tree"),"/"),  model="tree")

# 3.Perform assignment test uisng non-genetic data and random forest
assign.X( x1=ref.gen.rd, x2=assign.gen, dir=paste0(file.path(assign.dir,"randomForest"),"/"),  model="randomForest")



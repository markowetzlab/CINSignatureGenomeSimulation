## Calculate signature activities on simulated samples

rm(list=ls(all=TRUE))

## Libraries
library(this.path)
library(data.table)
library(reshape2)
library(YAPSA)
library(ggplot2)
library(patchwork)
library(ggthemes)
library(lemon)
library(Polychrome)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))


#### Paths and functions
thisPath=dirname(this.path())
IN=file.path(dirname(thisPath), "input")
OUT=file.path(dirname(thisPath), "output")

## Input file
# Signatures were derived with the SignatureDiscovery repo
SEGTABLES = file.path(OUT, "Out_2_Segmentation_tables_20each_N240.rds")

## Calculate SxC matrix
INPUTMODELS=file.path(IN, "2_combined_mixmodels_merged_components.rds")
UNINFPRIOR="TRUE"

# Signatures
SIGS=file.path(IN, "Signatures_XK5TRD_normalised.rds")
SIGNAMES = c("W3"="HRDLST", "W5"="ecDNA", "W1"="CHR", "W4"="WGDearly", "W2"="WGDlate")


## Output Files
CXS=file.path(OUT, "Out_3_Step_3_CxS_Matrix_20each_N240.txt")
NAMECXS=file.path(OUT, "Out_5_Activities_20each_N240")


## Load data
lDFStep3 = readRDS(SEGTABLES)
source(file.path(thisPath, "CINGenomeSimulation_Functions.R"))
source(file.path(IN, "main_functions.R"))



#### Step 1: Extract features
ecnfStep3 = extractCopynumberFeatures(lDFStep3, cores = 6, prePath=file.path(IN, ""), rmNorm = TRUE)
ecnfStep3$copynumber = NULL
saveRDS(ecnfStep3, file.path(OUT, "Out_3_ECNF_20each_N240.rds"))


#### Step 2: Calculate input matrix
allMatStep3 = calcInputMatrix(ecnfStep3, INPUTMODELS, UNINFPRIOR, NAME = "Out_4_SxC_Matrix_20each_N240", 
                              THISPATH = file.path(OUT, ""), SAVE = TRUE)
pInputStep3 = plotInputMatrix(allMatStep3, lDFStep3, NAME = "Out_4_CxS_Matrix_20each_N240", THISPATH = OUT)


#### Step 3: Quantify signature activities
HCxS = calcSigActivities(CXS = allMatStep3, SIGS = SIGS, BASENAME = NAMECXS, SIGNAMES = SIGNAMES)


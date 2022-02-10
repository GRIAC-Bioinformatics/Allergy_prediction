# -----------------------------------------------------------------------------------------------------
## Function:    loadPackages()
## Input:       None
## Output:      None
## Description: Loads all the necessary packages for the remainder of this script, 
##              installs them if needed
# -----------------------------------------------------------------------------------------------------
loadPackages <- function(){
  
  if(!require(tictoc))        {install.packages("tictoc")};        library(tictoc)
  if(!require(flux))          {install.packages("flux")};          library(flux)
  if(!require(klaR))          {install.packages("klaR")};          library(klaR)
  if(!require(caret))         {install.packages("caret")};         library(caret)
  if(!require(caretEnsemble)) {install.packages("caretEnsemble")}; library(caretEnsemble)
  if(!require(glmnet))        {install.packages("glmnet")};        library(glmnet)
  if(!require(gbm))           {install.packages("gbm")};           library(gbm)
  if(!require(randomForest))  {install.packages("randomForest")};  library(randomForest)
  if(!require(foreign))       {install.packages("foreign")};       library(foreign)
  if(!require(R.utils))       {install.packages("R.utils")};       library(R.utils)
  if(!require(DMwR))          {install.packages("DMwR")};          library(DMwR)
  if(!require(doParallel))    {install.packages("doParallel")};    library(doParallel)
  if(!require(data.table))    {install.packages("data.table")};    library(data.table)
  if(!require(neuralnet))     {install.packages("neuralnet")};     library(neuralnet)
  if(!require(xgboost))       {install.packages("xgboost")};       library(xgboost)
  if(!require(e1071))         {install.packages("e1071")};         library(e1071)
  if(!require(nnet))          {install.packages("nnet")};          library(nnet)
  if(!require(RColorBrewer))  {install.packages("RColorBrewer")};  library(RColorBrewer)
  if(!require(fastAdaboost))  {install.packages("fastAdaboost")};  library(fastAdaboost)
  if(!require(ada))           {install.packages("ada")};           library(ada)
  if(!require(plyr))          {install.packages("plyr")};          library(plyr)
  if(!require(dplyr))         {install.packages("dplyr")};         library(dplyr)
  if(!require(rlist))         {install.packages("rlist")};         library(rlist)
  if(!require(cutpointr))     {install.packages("cutpointr")};     library(cutpointr)
  if(!require(tidyverse))     {install.packages("tidyverse")};     library(tidyverse)
  if(!require(naivebayes))    {install.packages("naivebayes")};    library(naivebayes)
  if(!require(klaR))          {install.packages("klaR")};          library(klaR)
  if(!require(bnlearn))       {install.packages("bnlearn")};       library(bnlearn)
  if(!require(ggplot2))       {install.packages("ggplot2")};       library(ggplot2)
  if(!require(ROSE))          {install.packages("ROSE")};          library(ROSE)
  if(!require(PRROC))         {install.packages("PRROC")};         library(PRROC)
  if(!require(precrec))       {install.packages("precrec")};       library(precrec)
  if(!require(limma))         {install.packages("limma")};         library(limma)
  if(!require(MLmetrics))     {install.packages("MLmetrics")};     library(MLmetrics)
  if(!require(PerfMeas))      {install.packages("PerfMeas")};      library(PerfMeas)
  if(!require(pROC))          {install.packages("pROC")};          library(pROC)
  if(!require(ranger))        {install.packages("ranger")};        library(ranger)
  if(!require(DescTools))     {install.packages("DescTools")};     library(DescTools)
  if(!require(R.utils))       {install.packages("R.utils")};       library(R.utils)
  if(!require(tidyr))         {install.packages("tidyr")};         library(tidyr)  
  if(!require(reshape2))      {install.packages("reshape2")};      library(reshape2)  
  
}
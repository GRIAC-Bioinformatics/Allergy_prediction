# Allergy_prediction

Rationale: 
Childhood allergic diseases, including asthma, rhinitis and eczema, are prevalent conditions that share strong genetic and environmental components. Diagnosis relies on clinical history and measurements of allergen-specific IgE. We hypothesize that a multi-omics model could accurately diagnose childhood allergic disease.

Link to the paper:

Contact: <br>
Merlijn van Breugel: MvanBreugel@micompany.nl <br>
Cancan Qi: tracyqican@gmail.com <br>
Chengjian Xu: Xu.Chengjian@mh-hannover.de <br>

## Content of the scripts used in this project
### 1. The scripts of model build (folder "Main")
* model_validation_main.R (generating the main model or apply 3-CpG site model to another dataset)
* support functions
* validation_model.rds (final 3-CpG site model, can be used to predict allergy/non-allergy in another dataset using code from model_validation_main.R)

### 2. The scripts of main figures (folder "Figures")
* boxplot of sub-phenotype analysis
* single-cell RNAseq figures

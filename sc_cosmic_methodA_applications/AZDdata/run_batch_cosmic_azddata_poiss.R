
library(reticulate)
library(anndata)
library(SparseM)
library(rhdf5)
library(testit)
library(data.table)
library(Seurat)
library(RcppHNSW)
library(rsvd)
library(pracma)
library(progress)
library(xgboost)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
context1 <- toString(args[1])
context2 <- toString(args[2])
celltype_to_run <- toString(args[3])
counterfac_method <- toString(args[4])
output_cell <- toString(args[5])
output_name <- toString(args[6])

# context1="Healthy"
# context2="Inflamed"
# celltype_to_run="TA"
# counterfac_method = "svd_knn"
# output_cell="/n/groups/price/kushal/Causal/output/ucdata"
# output_name = "uc_poiss"

h5ad_file = "/n/groups/price/karthik/scdata/alzheimers.h5ad"

scdata = h5read(h5ad_file, name = "layers/counts")$data
scdata[scdata > 500] = 500
indices_scdata = h5read(h5ad_file, name = "layers/counts")$indices
indptr_scdata = h5read(h5ad_file, name = "layers/counts")$indptr
metadata = h5read(h5ad_file, name = "obs")
genedata = h5read(h5ad_file, name = "var")
context_label="DiseaseStatus"
subject_label = "Subject"
celltypes_label = "annot_level_2"


###########################  Creating context metadata annotation for each cell #########################################

context_types = metadata[['__categories']][[paste0(context_label)]]
cat("The context labels are as follows:\n", context_types, "\n")

if(min(metadata[[paste0(context_label)]]) == 0){
        context_ids = metadata[[paste0(context_label)]]+1
}else{
        context_ids = metadata[[paste0(context_label)]]
}
context_metadata = context_types[context_ids]
context_metadata[context_metadata == "no-pathology"] = "Healthy"
context_metadata[context_metadata == "early-pathology"] = "Disease"
context_metadata[context_metadata == "late-pathology"] = "Disease"

###########################  Creating individual metadata annotation for each cell #########################################

indiv_types = metadata[['__categories']][[paste0(subject_label)]]
if(min(metadata[[paste0(subject_label)]]) == 0){
                indiv_ids = metadata[[paste0(subject_label)]]+1
}else{
                indiv_ids = metadata[[paste0(subject_label)]]
}
indiv_metadata = indiv_types[indiv_ids]

###########################  Creating cell type metadata annotation for each cell #########################################

cell_types = metadata[['__categories']][[paste0(celltypes_label)]]
if(min(metadata[[paste0(celltypes_label)]]) == 0){
        cell_ids = metadata[[paste0(celltypes_label)]]+1
}else{
        cell_ids = metadata[[paste0(celltypes_label)]]
}
celltype_metadata = cell_types[cell_ids]

gene_metadata = data.frame("Gene" = genedata[['_index']],
						   "score" = 1)

sample_metadata = data.frame("id" = metadata$TAG,
							 "indiv" = indiv_metadata,
							 "context" = context_metadata,
							 "celltype" = celltype_metadata)

S4data = new("matrix.csr",
             "ra" = as.numeric(scdata), ## entry in each column split by rows, must be numeric
             "ja" = as.integer(indices_scdata+1), ## column (gene) indices for each element
             "ia" = as.integer(indptr_scdata+1),  ## (length of this array will be number of samples + 1), +1 added because of move from python to R
             "dimension" = c(nrow(gene_metadata), nrow(sample_metadata)))

S4data_sparse = as(S4data, "sparseMatrix")

pack = list()
pack$input_data = S4data_sparse
pack$sample_metadata = sample_metadata
pack$gene_metadata = gene_metadata

source("/n/groups/price/kushal/Causal/codes/sc_cosmic/svd_knn_counterfactual.R")
source("/n/groups/price/kushal/Causal/codes/sc_cosmic/propensity_knn_counterfactual.R")
source("/n/groups/price/kushal/Causal/codes/sc_cosmic/propensity_svd_knn_counterfactual.R")
source("/n/groups/price/kushal/Causal/codes/sc_cosmic/svd_smooth_knn_counterfactual.R")
source("/n/groups/price/kushal/Causal/codes/sc_cosmic/propensity_svd_smooth_knn_counterfactual.R")
source("/n/groups/price/kushal/Causal/codes/sc_cosmic/utils.R")
source("/n/groups/price/kushal/Causal/codes/sc_cosmic/cosmic_poisson_run.R")


outt = cosmic_poisson_run(input_data = pack$input_data,
		  sample_metadata = pack$sample_metadata,
		  gene_metadata = pack$gene_metadata,
		  context1= context1,
		  context2 = context2,
		  NUMITER_VI = 500,
		  knn_k=5,
		  rank_svd=50,
		  prior_gamma0=1,
		  prior_gamma1=1,
		  CROSSVAL_BATCH_NUM = 5,
		  NUMITER_XGBoost = 200,
		  MIN_NUM_CELLTYPES=1000,
		  counterfac_method = counterfac_method,
		  poiss_lite=F,
		  geneselect = T,
		  estimate_mu = F,
		  celltypes_to_use = celltype_to_run,
		  output_cell = output_cell,
		  output_name = output_name)
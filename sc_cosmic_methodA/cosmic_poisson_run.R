#' sc-COSMIC: Poisson counterfactual inference model for single-cell RNA-seq data
#'
#' The \code{cosmic_poisson_run} function takes in a matrix/data frame (G x N), a metadata matrix N x 3 matrix with
#' columns representing context, individual and cell-type annotations for each sample. The data must be counts. The 
#' code performs Poisson DE model for counterfactual causal context-enriched genes. 
#' @param input_data matrix/sparse-matrix or data frame of RNA counts, with samples along columns and genes along rows.
#' @param sample_metadata NxM matrix with N samples along rows and at least 3 columns context, indiv and celltype
#'       representing context, individual and cell-type annotations for each sample
#' @param gene_metadata GxM matrix with G genes along rows and any other columns
#' @param select_genes List of genes from gene metadata to be used as baseline for causal inference. Default is NULL in which case, Variational DE genes 
#'                     are used. If select_genes are not NULL, these genes are added on top of Variational DE genes. 
#' @param context1 Value of 1st context in context column of sample_metadata to use in a 2-way comparison
#' @param context2 Value of 2nd context in context column of sample_metadata to use in a 2-way comparison
#' @param NUMITER_VI Number of iterations in the Mean field variational inference model. Default is 500. 
#' @param knn_k_bound The number of nearest neighbors selected for performing counterfactual imputation
#' @param rank_svd Rank of the SVD used for counterfactual imputation. Default is 30. 
#' @param prior_gamma0 The a0 paameter of the prior Gamma distribution on library size. Default is equal to 1.
#' @param prior_gamma1 The a1 paameter of the prior Gamma distribution on library size. Default is equal to 1.
#' @param caliper_bound_vec The set of caliper bounds to use for balance assessment in matching. The default is the vector c(0.10, 0.15, 0.20, 0.25).
#' @param CROSSVAL_BATCH_NUM Integer representing number of cross-validation batches in propensity score calculation. Default is 5.
#' @param NUMITER_XGBoost Number of iterations in the gradient boosting propensity score model. Default is 200. 
#' @param MIN_NUM_CELLTYPES Minimum number of cells in a cell-type required for consideration in sc-cosmic.
#' @param NUM_ITER_COUNTERFAC The number of iterated counterfactuals to generate to compare with observed data. 
#' @param counterfac_method Method for detecting counterfactual neighbors for each observation. Can either be 
#'								propensity_knn: kNN based on propensity score with a caliper from `caliper_bound_vec`,
#'							  svd_knn: kNN based on SVD on the data matrix
#'								propensity_svd_knn: kNN on SVD with propensity score caliper with a caliper bound from `caliper_bound_vec`
#'								svd_knn_smooth: `svd_knn` approach to counterfactual est. but using response surface smoothing. 
#'								propensity_svd_knn_smooth: `propensity_svd_knn` approach to counterfactual est. but using response surface smoothing. 
#' @param geneselect boolean, if TRUE, we select genes that are non-causal DE to perform causal analysis. Else, no gene filtering is done.
#' @param estimate_mu boolean, if TRUE, estimates confounder driven expression effects before computing causal effects. Default is FALSE. 
#' @param celltypes_to_use A vector with cell type names for the ones we want to report using COSMIC.
#' @param output_cell Diectory where output from the COSMIC model is stored
#' @param output_name Name of the output file for saving the output
#' @return Outputs two files - one for non-causal and other for causal DE analysis in the output_cell provided by the user
#' @export

cosmic_poisson_run = function(input_data,
							 sample_metadata,
							 gene_metadata,
							 select_genes,
							 context1,
							 context2,
							 NUMITER_VI=500,
							 knn_k_bound = 5,
							 rank_svd = 50,
							 prior_gamma0 = 1, 
							 prior_gamma1 = 1,
							 caliper_bound_vec = c(0.10, 0.15, 0.20, 0.25),
							 CROSSVAL_BATCH_NUM = 5,
							 NUMITER_XGBoost = 200,
							 MIN_NUM_CELLTYPES = 5000,
							 NUM_ITER_COUNTERFAC=100,
							 counterfac_method = c("propensity_knn", "svd_knn", "propensity_svd_knn"), 
							 geneselect = TRUE,
							 estimate_mu=FALSE,
							 celltypes_to_use=NA,
							 output_cell,
							 output_name){

	library(reticulate)
	library(anndata)
	library(SparseM)
	library(rhdf5)
	library(testit)
	library(data.table)
	library(Seurat)
	library(RcppHNSW)
	library(irlba)
	library(rsvd)
	library(pracma)
	library(progress)
	library(xgboost)
	library(MatchIt)

	if(!is.matrix(input_data) & !is.data.frame(input_data) & class(input_data)[1] != "dgCMatrix"){
		stop("The input data matrix must be a matrix, sparse dgCMatrix or a data frame")
	}else{
		input_data = as(input_data, "sparseMatrix")
	}

	
	if(!identical(ceiling(input_data), input_data) | length(which(input_data < 0)) > 0){
		stop("The input data matrix must be non-negative counts")
	}


	if(!is.matrix(sample_metadata) & !is.data.frame(sample_metadata)){
		stop("sample_metadata must be a matrix or data frame")
	}

	if(!is.matrix(gene_metadata) & !is.data.frame(gene_metadata)){
		stop("gene_metadata must be a matrix or data frame")
	}
	
	if(ncol(input_data) != nrow(sample_metadata)){
		stop("The number of columns of input data must match the number of rows in sample metadata")
	}

	if(nrow(input_data) != nrow(gene_metadata)){
		stop("The number of rows of input data must match the number of rows in gene metadata")
	}

	context_metadata = sample_metadata$context
	celltype_metadata = sample_metadata$celltype
	indiv_metadata = sample_metadata$indiv
	contingency_celltypes = table(celltype_metadata)
	select_cell_types = names(contingency_celltypes)[which(contingency_celltypes > MIN_NUM_CELLTYPES)]

	if(is.na(celltypes_to_use)){
		select_cell_types = select_cell_types
	}else{
		select_cell_types = intersect(select_cell_types, celltypes_to_use)
		if(length(select_cell_types) == 0){
			stop("The number of cell types selected after filtering for adequate number of cells and/or 
				user-specified list of cell types is 0")
		}
	}

	cat("Processing COSMIC for cell types\n")
	print(select_cell_types)
	cat("\n")
	gamma_a0 = prior_gamma0;
	gamma_b0 = prior_gamma1;
	gene_names = gene_metadata$Gene

	select_genes_list = list()
	
	for(counter1 in 1:length(select_cell_types)){
		cat("Performing Non-causal inference for cell type:", select_cell_types[counter1], "\n")

		## select cells from the given cell-type ##
		## filtered data and metadata restricting to the cells for a given cell type 

		S4data_sparse_celltype = input_data[, which(celltype_metadata == select_cell_types[counter1])]
		indiv_celltype_metadata = indiv_metadata[which(celltype_metadata == select_cell_types[counter1])]
		context_celltype_metadata = context_metadata[which(celltype_metadata == select_cell_types[counter1])]

		## select and reorder cells belonging to the two contexts of interest. ##

		idx2 = c(which(context_celltype_metadata  == context1),
				 which(context_celltype_metadata  == context2))
		context_celltype_metadata2 = context_celltype_metadata[idx2]
		indiv_celltype_metadata2 =  indiv_celltype_metadata[idx2]
		S4data_sparse_celltype2 = S4data_sparse_celltype[,idx2]

		merged_indiv_health_celltype_metadata2 = paste0(indiv_celltype_metadata2, "_", context_celltype_metadata2)
		ff = factor(merged_indiv_health_celltype_metadata2, levels = sort(unique(merged_indiv_health_celltype_metadata2)))
		mm = as(model.matrix(~ff-1), "sparseMatrix")
		colnames(mm) = sort(unique(merged_indiv_health_celltype_metadata2))

		size_fac = colSums(S4data_sparse_celltype2)     #### \sum_{g} Y_{gj} = M_{j} 
		YYindiv = S4data_sparse_celltype2 %*% mm  
	
		set.seed(10)
		rhos_start = rgamma(ncol(S4data_sparse_celltype2), shape=gamma_a0, rate=gamma_b0) 
		rhos = rhos_start
		pb <- progress_bar$new(total = 4*NUMITER_VI)
		for(niter in 1:4*NUMITER_VI){
			rhos_indiv = tapply(rhos, merged_indiv_health_celltype_metadata2, sum)
			lambdas_indiv = sweep(YYindiv+1, 2, rhos_indiv+1, "/")  ### \lambda_{g, i} 
			libsize_indiv = colSums(lambdas_indiv)
			libsize_indiv_cells = as.numeric(libsize_indiv[merged_indiv_health_celltype_metadata2])
	        rhos = (size_fac+1)/(libsize_indiv_cells+1)
			pb$tick()
			Sys.sleep(1 / NUMITER_VI)
		}
		
		xx_tabb = cbind(sapply(colnames(lambdas_indiv), function(x) return(strsplit(x, "_")[[1]][1])),
	      sapply(colnames(lambdas_indiv), function(x) return(strsplit(x, "_")[[1]][2])))

		idx11 = which(xx_tabb[,2] == context1)
		idx12 = which(xx_tabb[,2] == context2)
		p_diff1 = apply(lambdas_indiv, 1,
		function(yy){
			tt = wilcox.test(yy[idx12], yy[idx11], alternative = "two.sided")
			return(abs(tt$p.value))
		})

		if(geneselect){
			if(!is.null(select_genes)){
				select_genes = intersect(select_genes, gene_names)
				select_genes_list[[counter1]] = union(gene_names[which(p_diff1 < 0.1)], select_genes)
			}else{
				select_genes_list[[counter1]] =  gene_names[which(p_diff1 < 0.1)]
			}	
		}else{
			select_genes_list[[counter1]] = gene_names
		}

		is_select = rep(0, nrow(gene_metadata))
		is_select[match(intersect(select_genes, gene_metadata$Gene),gene_metadata$Gene)] = 1

		noncausal_VI_genes = cbind.data.frame(select_cell_types[counter1], gene_metadata$Gene, p_diff1, is_select)
		colnames(noncausal_VI_genes) = c( "Cell.type", "Gene", "p.value", "is.Select")
		write.table(noncausal_VI_genes, file = paste0(output_cell, "/", "Noncausal_VIDE_", output_name, "_",
					select_cell_types[counter1], "_",
					context1, "_", context2,  ".txt"),
					row.names=F, col.names=T, sep = "\t", quote=F)
		cat("Finished performing COSMIC non-causal DE analysis for cell type:", select_cell_types[counter1], "\n")
	}

	
	for(counter1 in 1:length(select_cell_types)){
		cat("Performing Causal inference for cell type:", select_cell_types[counter1], "\n")

		## select cells from the given cell-type ##
		S4data_sparse_celltype = input_data[, which(celltype_metadata == select_cell_types[counter1])]
		indiv_celltype_metadata = indiv_metadata[which(celltype_metadata == select_cell_types[counter1])]
		context_celltype_metadata = context_metadata[which(celltype_metadata == select_cell_types[counter1])]

		## select cells belonging to the two contexts of interest. ##
		idx2 = c(which(context_celltype_metadata  == context1),
				 which(context_celltype_metadata  == context2))
		context_celltype_metadata2 = context_celltype_metadata[idx2]
		indiv_celltype_metadata2 =  indiv_celltype_metadata[idx2]
		S4data_sparse_celltype2 = S4data_sparse_celltype[,idx2]

		merged_indiv_health_celltype_metadata2 = paste0(indiv_celltype_metadata2, "_", context_celltype_metadata2)
		S4data_sparse_celltype2 = S4data_sparse_celltype2[match(select_genes_list[[counter1]], gene_names), ]

		if(counterfac_method == "propensity_knn"){
			counterfac_neighbors_ll = propensity_knn_counterfactual(S4data_sparse_celltype2,
					context_celltype_metadata2,
					knn_k_bound = knn_k_bound,
					rank_svd = rank_svd,
					caliper_bound_vec = caliper_bound_vec,
					CROSSVAL_BATCH_NUM=CROSSVAL_BATCH_NUM,
					NUMITER_XGBoost=NUMITER_XGBoost)
		}else if (counterfac_method == "svd_knn"){
			counterfac_neighbors_ll = svd_knn_counterfactual(S4data_sparse_celltype2,
					context_celltype_metadata2,
					knn_k_bound = knn_k_bound,
					rank_svd = rank_svd)
		}else if (counterfac_method == "propensity_svd_knn"){
			counterfac_neighbors_ll = propensity_svd_knn_counterfactual(S4data_sparse_celltype2,
					context_celltype_metadata2,
					knn_k_bound = knn_k_bound,
					rank_svd = rank_svd,
					caliper_bound_vec = caliper_bound_vec,
					CROSSVAL_BATCH_NUM=CROSSVAL_BATCH_NUM,
					NUMITER_XGBoost=NUMITER_XGBoost)
		}else{
				stop("The counterfac method chosen does not match with the possible options.")
		}

		merged_zz_context = counterfac_neighbors_ll$counterfac_neighbors
		propen_weights_merged_zz_context = counterfac_neighbors_ll$counterfac_weights
		S4data_sparse_celltype3 = as(counterfac_neighbors_ll$X, "sparseMatrix")
		knn_k = ncol(merged_zz_context)

    ################################## Confounding and causal effects calculation. ##################################

		
		for(iter_counterfac in 1:NUM_ITER_COUNTERFAC){

			assert(nrow(merged_zz_context) == nrow(propen_weights_merged_zz_context))
			col_ids = apply(propen_weights_merged_zz_context, 1, function(x) return(sample(1:knn_k, size = 1, prob = x)))
			row_ids = 1:nrow(merged_zz_context)
			sampled_neighbors = merged_zz_context[cbind(row_ids, col_ids)]
			

			counterfac = S4data_sparse_celltype3[, sampled_neighbors]
			if(length(counterfac_neighbors_ll$cells_removed) > 0){
				observed = S4data_sparse_celltype3[, -counterfac_neighbors_ll$cells_removed]
				context_celltype_metadata3 = context_celltype_metadata2[-counterfac_neighbors_ll$cells_removed]
				context_celltype_metadata3_binary = rep(0, length(context_celltype_metadata3))
				context_celltype_metadata3_binary[which(context_celltype_metadata3 == context1)] = 0
				context_celltype_metadata3_binary[which(context_celltype_metadata3 == context2)] = 1
				
				context1_fulldata = cbind(observed[, which(context_celltype_metadata3 == context1)], counterfac[, which(context_celltype_metadata3 == context2)])
				context2_fulldata = cbind(counterfac[, which(context_celltype_metadata3 == context1)], observed[, which(context_celltype_metadata3 == context2)])
				context1_fulldata = as(context1_fulldata, "sparseMatrix")
				context2_fulldata = as(context2_fulldata, "sparseMatrix")
				merged_indiv_health_celltype_metadata3 = merged_indiv_health_celltype_metadata2[-counterfac_neighbors_ll$cells_removed]

			}else{
				observed = S4data_sparse_celltype3
				context_celltype_metadata3 = context_celltype_metadata2
				context_celltype_metadata3_binary = rep(0, length(context_celltype_metadata3))
				context_celltype_metadata3_binary[which(context_celltype_metadata3 == context1)] = 0
				context_celltype_metadata3_binary[which(context_celltype_metadata3 == context2)] = 1
				
				context1_fulldata = cbind(observed[, which(context_celltype_metadata3 == context1)], counterfac[, which(context_celltype_metadata3 == context2)])
				context2_fulldata = cbind(counterfac[, which(context_celltype_metadata3 == context1)], observed[, which(context_celltype_metadata3 == context2)])
				context1_fulldata = as(context1_fulldata, "sparseMatrix")
				context2_fulldata = as(context2_fulldata, "sparseMatrix")
				merged_indiv_health_celltype_metadata3 = merged_indiv_health_celltype_metadata2
			}
			

			
			ff = factor(merged_indiv_health_celltype_metadata3, levels = sort(unique(merged_indiv_health_celltype_metadata3)))
			model_indiv = as(model.matrix(~ff-1), "sparseMatrix")
			colnames(model_indiv) = sort(unique(merged_indiv_health_celltype_metadata3))
			
		
			if(estimate_mu){

				YY_counterfac_indiv = counterfac %*% model_indiv
				YY_observed_indiv = observed %*% model_indiv
				size_fac_observed = colSums(observed)
				size_fac_counterfac = colSums(counterfac)

				cat("Estimating shared causal effect shared across observed and counterfactual\n")
				rhos_start0 = rgamma(ncol(observed), gamma_a0, gamma_b0) 
				rhos_start1 = rgamma(ncol(counterfac), gamma_a0, gamma_b0) 
				rhos0 = rhos_start0
				rhos1 = rhos_start1
				pb <- progress_bar$new(total = NUMITER_VI)
				for(niter in 1:NUMITER_VI){
					rhos = rhos0 + rhos1
					rhos_indiv = tapply(rhos, merged_indiv_health_celltype_metadata2, sum)
					mu_indiv = sweep(YY_counterfac_indiv + YY_observed_indiv+1, 2, 
					                 rhos_indiv+1, "/")  ### \mu_{g, i} 
					libsize_indiv = colSums(mu_indiv)
					names(libsize_indiv) = sort(unique(merged_indiv_health_celltype_metadata2))
					libsize_indiv_cells = as.numeric(libsize_indiv[as.character(merged_indiv_health_celltype_metadata2)])
			        rhos0 = (size_fac_observed+1)/(libsize_indiv_cells+1)
			        rhos1 = (size_fac_counterfac+1)/(libsize_indiv_cells+1)
					pb$tick()
					Sys.sleep(1 / NUMITER_VI)
				  }

				cat("Estimating context-specific causal effect\n")
			  	rhos_start = rgamma(ncol(observed), gamma_a0, gamma_b0) 
				rhos = rhos_start
				pb <- progress_bar$new(total = NUMITER_VI)
				for(niter in 1:NUMITER_VI){
					rhos_indiv = tapply(rhos, merged_indiv_health_celltype_metadata2, sum)
					deltas_indiv = (1+YY_observed_indiv)/(1+ sweep(mu_indiv, 2, rhos_indiv, "*"))  ### \delta_{g, i} 
					denn = colSums(deltas_indiv*mu_indiv)+1
					names(denn) = sort(unique(merged_indiv_health_celltype_metadata2))
					numm = colSums(observed) + 1
					rhos = numm/as.numeric(denn[as.character(merged_indiv_health_celltype_metadata2)])
					pb$tick()
					Sys.sleep(1 / NUMITER_VI)
					}
				deltas_indiv = log(deltas_indiv+1e-08)

				}else{
					context1_indiv = context1_fulldata %*% model_indiv
					context2_indiv = context2_fulldata %*% model_indiv
					size_fac_context1 = colSums(context1_indiv) 
					size_fac_context2 = colSums(context2_indiv)

					cat("Estimating context-specific causal effect\n")
					  rhos_indiv = rgamma(ncol(context1_indiv), gamma_a0, gamma_b0) 
						pb <- progress_bar$new(total = NUMITER_VI)
						for(niter in 1:NUMITER_VI){
							context1_times_rho_indiv = sweep(context1_indiv, 2, rhos_indiv, "*") 
							deltas_indiv = (1+context2_indiv)/(1+ context1_times_rho_indiv) 
							rhos_indiv = (colSums(context2_indiv) + 1)/(colSums(context1_indiv*deltas_indiv) + 1)
							pb$tick()
							Sys.sleep(1 / NUMITER_VI)
					}
					deltas_indiv = log(deltas_indiv+1e-08)	
			}

			if(iter_counterfac == 1){
				deltas_indiv_full = deltas_indiv
			}else{
				deltas_indiv_full = deltas_indiv_full + deltas_indiv
			}
			cat("We are at counterfactual estimation step:", iter_counterfac, "\n")
		}

		deltas_indiv_full = deltas_indiv_full/NUM_ITER_COUNTERFAC

		
		############################  Identify causal genes based on causal effects across individuals  #################################

		p_diff2 = apply(as.matrix(deltas_indiv_full), 1,
				function(yy){
					tt = pnorm(abs(sqrt(length(yy))*(mean(yy)/sd(yy))), 0,  1, lower.tail=F)
					return(tt)
					})
		p_diff2[is.na(p_diff2)] = 1

		causal_VI_genes = cbind.data.frame(select_cell_types[counter1], select_genes_list[[counter1]], p_diff2)
		colnames(causal_VI_genes) = c("Cell.type", "Gene", "p.value.combo")
		write.table(causal_VI_genes, file = paste0(output_cell, "/", "Causal_VIDE_", output_name, "_",
					select_cell_types[counter1], "_",
					context1, "_", context2, "_", counterfac_method, ".twosided.txt"),
					row.names=F, col.names=T, sep = "\t", quote=F)
		cat("Finished performing COSMIC causal DE analysis for cell type:", select_cell_types[counter1], "\n")


		############################  Save the output from sc-cosmic   #################################

		outll = list("balance_before_match" = counterfac_neighbors_ll$balance_before_match,
			  "balance_after_match" = counterfac_neighbors_ll$balance_after_match,
			  "balance_sheet" = counterfac_neighbors_ll$balance_sheet,
			  "is_balance" = counterfac_neighbors_ll$is_balance,
			  "std_diff_before_match" = counterfac_neighbors_ll$std_diff_before_match,
			  "std_diff_after_match" = counterfac_neighbors_ll$std_diff_after_match,
			  "stdDE_gene_results" = noncausal_VI_genes,
			  "causalDE_gene_results" = causal_VI_genes,
			  "causal_eff_matrix" = deltas_indiv_full
			  )

		save(outll, file = paste0(output_cell, "/", "scCOSMIC_", output_name, "_",
					select_cell_types[counter1], "_",
					context1, "_", context2, "_", counterfac_method, ".rda"))
	}
}
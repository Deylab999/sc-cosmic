#' Find counterfactual neighbors based on SVD latent factor decomposition of cells in sc-COSMIC
#'
#' The \code{svd_knn_counterfactual} function takes in a matrix/data frame (G x N), a context metadata vector of length N 
#' comprising of two distint labels (corresponding to context1 and context2) and uses a SVD based latent factor decomposition nmethod
#' to generate a set of nearest neighbors counterfactuals of opposite context to each cell, along their 
#' relative weights and data matrix for downstream analysis.
#' @param Xmat matrix/sparse-matrix or data frame of RNA counts, with samples along columns (N) and genes along rows (G).
#' @param context_vec A vector of length N containing the context information for each cell.
#' @param knn_k_bound Upper bound on the number of nearest neighbors to assess for balance assessment. Default is 5. 
#' @param rank_svd Rank of the SVD used for counterfactual imputation. Default is 50. 
#' @return Outputs a list with three elements 
#'         `counterfac_neighbors`: A N x knn_k matrix reporting nearest knn_k neighbors for each cell n.
#'         `counterfac_weights`: A N x knn_k weight matrix reporting weights for the nearest neighbors
#'         `X`: A GxN matrix that is same as the input data matrix.
#' @export 

svd_knn_counterfactual = function(Xmat,
								 context_vec,
								 knn_k_bound=5,
								 rank_svd=50){

	if(!is.matrix(Xmat) & !is.data.frame(Xmat) & class(Xmat)[1] != "dgCMatrix"){
		stop("The input data matrix must be a matrix, sparse dgCMatrix or a data frame")
	}else{
		Xmat = as(Xmat, "sparseMatrix")
	}

	if(ncol(Xmat) != length(context_vec)){
		stop("The number of columns of input data must match the number of rows in sample metadata")
	}
	if(length(unique(context_vec)) !=2){
		stop("The number of labels (distinct elements) in context_vec must be 2; representing context1 and context2")
	}

	context1 = unique(context_vec)[1]
	context2 = unique(context_vec)[2]

	idx2 = c(which(context_vec == context1), which(context_vec == context2))
	context_vec2 = context_vec[idx2]
	Xmat2 = Xmat[,idx2]

	cat("Performing counterfactual imputation using SVD \n")

	context_vec2_binary = rep(0, length(context_vec2))
	context_vec2_binary[which(context_vec2 == context1)] = 0
	context_vec2_binary[which(context_vec2 == context2)] = 1

	svd_out = rsvd(Xmat2, k=rank_svd)
	svd_out_proj_combined = svd_out$v
	svd_out_proj_combined_t = t(svd_out$v)

	############################## Matchit balance criterion  ################################################

	ratio_vec = 2:knn_k_bound
	caliper_bound_vec = c(0.10, 0.15, 0.20, 0.25)

    balance_sheet = c()
    for(mm in 1:length(ratio_vec)){
    	for(nn in 1:length(caliper_bound_vec)){

    		#########   Assess balance in the covariates prior to performing matching  ####################

    		m.out0 <- matchit( as.factor(context_vec2_binary) ~ t(svd_out_proj_combined_t),
                 method = NULL, distance = "glm", caliper=caliper_bound_vec[nn], std.caliper=T, ratio = ratio_vec[mm], replace=T)
			balance_before_match = summary(m.out0)[3]$sum.all
			num_matched_0 = summary(m.out0)[2]$nn["Matched", 1]/summary(m.out0)[2]$nn["All", 1] ### should be equal to 1 
		 	num_imbalanced_0 = length(which(abs(balance_before_match[,3]) > 0.2)) ## 0.2 is the cut-off proposed by Rubin
		 	num_imbalanced_propvar_0 = length(which(abs(balance_before_match[,4]) < 0.5 | abs(balance_before_match[,4]) > 2))

			#########   Assess balance in the covariates prior after performing matching  ####################

			m.out1 <- matchit( as.factor(context_vec2_binary) ~ t(svd_out_proj_combined_t),
		                 method = "nearest", distance = "glm", caliper=caliper_bound_vec[nn], std.caliper=T, ratio = ratio_vec[mm], replace=T)
			balance_after_match = summary(m.out1)[4]$sum.matched
			num_matched_1 = summary(m.out1)[2]$nn["Matched", 1]/summary(m.out1)[2]$nn["All", 1]
			num_imbalanced_1 = length(which(abs(balance_after_match[,3]) > 0.2))
			num_imbalanced_propvar_1 = length(which(abs(balance_after_match[,4]) < 0.5 | abs(balance_after_match[,4]) > 2))

			#############  Summarize the results across calipers and knn  ##################################

			balance_sheet = rbind(balance_sheet,
				c(ratio_vec[mm], caliper_bound_vec[nn], num_matched_0, num_matched_1, num_imbalanced_0, num_imbalanced_1,
				num_imbalanced_propvar_0, num_imbalanced_propvar_1))
    	}
    }

    colnames(balance_sheet) = c("Ratio", "Caliper", "Num.match.pre", "Num.match.post", "Num.imbalance.pre", "Num.imbalance.post",
    	                        "Num.imbalance.prevar", "Num.imbalance.postvar")
	balance_sheet = as.data.frame(balance_sheet)

	idx_select = Reduce(intersect, 
		list(which(balance_sheet$Num.imbalance.post == min(balance_sheet$Num.imbalance.post)),
		which(balance_sheet$Num.imbalance.postvar == min(balance_sheet$Num.imbalance.postvar))))
	if(length(idx_select) > 0){
		idx_select2 = idx_select[which.max(balance_sheet$Num.match.post[idx_select])]
		caliper_choose = balance_sheet$Caliper[idx_select2]
		ratio_choose = balance_sheet$Ratio[idx_select2]
	}else{
		idx_select = which(balance_sheet$Num.imbalance.post == min(balance_sheet$Num.imbalance.post))
		idx_select2 = idx_select[which.max(balance_sheet$Num.match.post[idx_select])]
		caliper_choose = balance_sheet$Caliper[idx_select2]
		ratio_choose = balance_sheet$Ratio[idx_select2]
	}

	num_imbalanced_choose = pmax(balance_sheet$Num.imbalance.post[idx_select2], balance_sheet$Num.imbalance.postvar[idx_select2])
	if(num_imbalanced_choose >= ceiling(.1*nrow(svd_out_proj_combined_t))){
		warning("The number of imbalanced covariates after matching - using different calipers and knn - is still greater than 10% of covariates")
		is_balance=F
	}else{
		is_balance=T
	}


	#########################   Final match IT with chosen ratio and caliper  #############################################

	m.out0 <- matchit( as.factor(context_vec2_binary) ~ t(svd_out_proj_combined_t),
                 method = NULL, distance = "glm", caliper=caliper_choose, std.caliper=T, ratio = ratio_choose, replace=T)
	m.out1 <- matchit( as.factor(context_vec2_binary) ~ t(svd_out_proj_combined_t),
		                 method = "nearest", distance = "glm", caliper=caliper_choose, std.caliper=T, ratio = ratio_choose, replace=T)

	Xtbar = apply(svd_out_proj_combined_t, 1, function(x) return(mean(x[context_vec2 == context2])))
	Xcbar = apply(svd_out_proj_combined_t, 1, function(x) return(mean(x[context_vec2 == context1])))
	stbar = apply(svd_out_proj_combined_t, 1, function(x) return(sd(x[context_vec2 == context2])))
	std_diff_before_match = (Xtbar - Xcbar)/stbar


	knn_out = try(get_knn(svd_out_proj_combined, k=1000, distance = "euclidean"), silent=T)
	if(class(knn_out) == "try-error"){
		knn_out = try(get_knn(svd_out_proj_combined, k=1000, distance = "l2"), silent=T)
		cat("WARNING: Euclidean distance metric failed; We are using the L2 metric as distance in kNN computation \n")
	}
	if(class(knn_out) == "try-error"){
		knn_out = try(get_knn(svd_out_proj_combined, k=1000, distance = "ip"), silent=T)
		cat("WARNING: Euclidean distance metric failed; We are using the IP metric as distance in kNN computation \n")
	}
	if(class(knn_out) == "try-error"){
		knn_out = try(get_knn(svd_out_proj_combined, k=1000, distance = "cosine"), silent=T)
		cat("WARNING: Euclidean distance metric failed; We are using the cosine metric as distance in kNN computation \n")
	}
	if(class(knn_out) != "list"){
		stop("WARNING: Euclidean distance metric failed; The k nearest neighbor algorithm did not run successfully")
	}


	merged_zz_context = matrix(0, nrow(knn_out$idx), ratio_choose)
	propen_merged_zz_context = matrix(0, nrow(knn_out$idx), ratio_choose)

	for(i in 1:nrow(knn_out$idx)){
		if(context_vec2[i] == context1){
			fidx1 = knn_out$idx[i, which(knn_out$dist[i, ] < caliper_choose)]
			#fidx1 = knn_out$idx[i, which(knn_out$dist[i, ] < caliper_choose*sd(propensity_scores))]
			fidx2 = fidx1[which(context_vec2[fidx1] == context2)]

			if(length(fidx2) >= ratio_choose){
				merged_zz_context[i, ] = fidx2[1:ratio_choose]
				propen_merged_zz_context[i, ] = 1
			}else if (length(fidx2) >= 1){
				merged_zz_context[i, ] = c(fidx2, sample(fidx2, ratio_choose - length(fidx2), replace=T))
				propen_merged_zz_context[i, ] = c(rep(1, length(fidx2)), rep(0.01, ratio_choose - length(fidx2)))
			}else{
				merged_zz_context[i, ] = rep(NA, ratio_choose)
				propen_merged_zz_context[i, ] = 1
			}
		}else if (context_vec2[i] == context2){
			fidx1 = knn_out$idx[i, which(knn_out$dist[i, ] < caliper_choose)]
			#fidx1 = knn_out$idx[i, which(knn_out$dist[i, ] < caliper_choose*sd(propensity_scores))]
			fidx2 = fidx1[which(context_vec2[fidx1] == context1)]
			if(length(fidx2) >= ratio_choose){
				merged_zz_context[i, ] = fidx2[1:ratio_choose]
				propen_merged_zz_context[i, ] = 1
			}else if (length(fidx2) >= 1){
				merged_zz_context[i, ] = c(fidx2, sample(fidx2, ratio_choose - length(fidx2), replace=T))
				propen_merged_zz_context[i, ] = c(rep(1, length(fidx2)), rep(0.01, ratio_choose - length(fidx2)))
			}else{
				merged_zz_context[i, ] = rep(NA, ratio_choose)
				propen_merged_zz_context[i, ] = 1
			}
		}
	}

	propen_weights_merged_zz_context = t(apply(propen_merged_zz_context, 1, function(x) return(x/sum(x))))
	remove_cells = which(is.na(rowMeans(merged_zz_context)))
	if(length(remove_cells) >= 1){
		merged_zz_context = merged_zz_context[-remove_cells, ]
		propen_weights_merged_zz_context = propen_weights_merged_zz_context[-remove_cells, ]
		Xmat3 = Xmat2[, -remove_cells]
		context_vec3 = context_vec2[-remove_cells]
	}else{
		context_vec3 = context_vec2
	}

	res_bal = c()
	for(mm in 1:ratio_choose){
		tt1 = c(which(context_vec3 == context2), merged_zz_context[which(context_vec3 == context1), mm])
		cc1 = c(which(context_vec3 == context1), merged_zz_context[which(context_vec3 == context2), mm])
		Xtbar = apply(svd_out_proj_combined_t, 1, function(x) return(mean(x[tt1])))
		Xcbar = apply(svd_out_proj_combined_t, 1, function(x) return(mean(x[cc1])))
		stbar = apply(svd_out_proj_combined_t, 1, function(x) return(sd(x[tt1])))
	    res_bal = rbind(res_bal, (Xtbar - Xcbar)/stbar)
	}
	std_diff_after_match = apply(res_bal, 2, mean)


	ll = list("counterfac_neighbors" = merged_zz_context,
			  "counterfac_weights" = propen_weights_merged_zz_context,
			  "X" = Xmat2,
			  "balance_before_match" = m.out0,
			  "balance_after_match" = m.out1,
			  "balance_sheet" = balance_sheet,
			  "std_diff_before_match" = std_diff_before_match,
			  "std_diff_after_match" = std_diff_after_match,
			  "cells_removed" = remove_cells,
			  "is_balance" = is_balance)
	return(ll)
}
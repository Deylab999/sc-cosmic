#' Find counterfactual neighbors based on SVD latent factor decomposition of cells along with a response surface smoothing in sc-COSMIC
#'
#' The \code{svd_smooth_knn_counterfactual} function takes in a matrix/data frame (G x N), a context metadata vector of length N 
#' comprising of two distint labels (corresponding to context1 and context2) and uses a SVD based latent factor decomposition on 
#' response surface smoothed data to generate a set of nearest neighbors counterfactuals of opposite context to each cell, along their 
#' relative weights and data matrix for downstream analysis.
#' @param Xmat matrix/sparse-matrix or data frame of RNA counts, with samples along columns (N) and genes along rows (G).
#' @param context_vec A vector of length N containing the context information for each cell.
#' @param knn_k The number of nearest neighbors. Default is 5. 
#' @param rank_svd Rank of the SVD used for counterfactual imputation. Default is 50. 
#' @param transform The type of transformation used before smoothing based on nature of data: `poisson` for counts data and `gauss` for normal data.
#'                  Default is `poisson`.
#' @return Outputs a list with three elements 
#'         `counterfac_neighbors`: A N x knn_k matrix reporting nearest knn_k neighbors for each cell n.
#'         `counterfac_weights`: A N x knn_k weight matrix reporting weights for the nearest neighbors
#'         `X`: A GxN matrix that is a response surface smoothed version of input data matrix.
#' @export 

svd_smooth_knn_counterfactual = function(Xmat,
								 context_vec,
								 knn_k=5,
								 rank_svd=50,
								 transform = "poisson"){
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
	context_vec2_binary = rep(0, length(context_vec2))
	context_vec2_binary[which(context_vec2 == context1)] = 0
	context_vec2_binary[which(context_vec2 == context2)] = 1

	svd_out = rsvd(Xmat2, k=rank_svd)
	svd_out_proj_combined = svd_out$v

	context1_data = Xmat2[, which(context_vec2_binary == 0)]
	context2_data = Xmat2[, which(context_vec2_binary == 1)]
	svd_out_proj_combined_context1 = svd_out_proj_combined[which(context_vec2_binary == 0),]
	svd_out_proj_combined_context2 = svd_out_proj_combined[which(context_vec2_binary == 1),]

	if(transform == "poisson"){
		betahat_mat_context1 = solve(t(svd_out_proj_combined_context1) %*% 
		    svd_out_proj_combined_context1 + 1*diag(rank_svd)) %*% t(svd_out_proj_combined_context1) %*% t(log(context1_data+1))
		betahat_mat_context2 = solve(t(svd_out_proj_combined_context2) %*% 
			svd_out_proj_combined_context2 + 1*diag(rank_svd)) %*% t(svd_out_proj_combined_context2) %*% t(log(context2_data+1))
		
		Y_smoothed_context1 = exp(t(betahat_mat_context1) %*% t(svd_out_proj_combined_context1))
		Y_smoothed_context2 = exp(t(betahat_mat_context2) %*% t(svd_out_proj_combined_context2))
		
		Y_smoothed = cbind(Y_smoothed_context1, Y_smoothed_context2)
		Y_smoothed = as(Y_smoothed, "sparseMatrix")
		Y_smoothed[Y_smoothed > 500] = 500
	}else if (transform == "gauss"){
		betahat_mat_context1 = solve(t(svd_out_proj_combined_context1) %*% 
		    svd_out_proj_combined_context1 + 1*diag(rank_svd)) %*% t(svd_out_proj_combined_context1) %*% t(context1_data)
		betahat_mat_context2 = solve(t(svd_out_proj_combined_context2) %*% 
			svd_out_proj_combined_context2 + 1*diag(rank_svd)) %*% t(svd_out_proj_combined_context2) %*% t(context2_data)
		Y_smoothed_context1 = t(betahat_mat_context1) %*% t(svd_out_proj_combined_context1)
		Y_smoothed_context2 = t(betahat_mat_context2) %*% t(svd_out_proj_combined_context2)
		
		Y_smoothed = cbind(Y_smoothed_context1, Y_smoothed_context2)
		Y_smoothed = as(Y_smoothed, "sparseMatrix")
	}else{
		stop("transform must either be poisson or gauss")
	}
	

	knn_out = get_knn(svd_out_proj_combined, k=1000, distance = "euclidean")

	merged_zz_context = matrix(0, nrow(knn_out$idx), knn_k)
	propen_merged_zz_context = matrix(1/knn_k, nrow(knn_out$idx), knn_k)
	
	for(i in 1:nrow(knn_out$idx)){
		if(context_vec2[i] == context1){
			fidx1 = knn_out$idx[i,]
			fidx2 = fidx1[which(context_vec2[fidx1] == context2)]
			if(length(fidx2) >= knn_k){
				merged_zz_context[i, ] = fidx2[1:knn_k]
			}else{
				merged_zz_context[i, ] = c(fidx2, sample(setdiff(which(context_vec2== context2), fidx2), knn_k - length(fidx2), replace=T))
				propen_merged_zz_context[i, (length(fidx2)+1):knn_k] = 1e-04 
			}
		}else if (context_vec2[i] == context2){
			fidx1 = knn_out$idx[i,]
			fidx2 = fidx1[which(context_vec2[fidx1] == context1)]
			if(length(fidx2) >= knn_k){
				merged_zz_context[i, ] = fidx2[1:knn_k]
			}else{
				merged_zz_context[i, ] = c(fidx2, sample(setdiff(which(context_vec2 == context1), fidx2), knn_k - length(fidx2), replace=T))
				propen_merged_zz_context[i, (length(fidx2)+1):knn_k] = 1e-04 
			}
		}
	}
	propen_weights_merged_zz_context = t(apply(propen_merged_zz_context, 1, function(x) return(x/sum(x))))
	Xmat3 = Y_smoothed

	ll = list("counterfac_neighbors" = merged_zz_context,
			  "counterfac_weights" = propen_weights_merged_zz_context,
			  "X" = Xmat3)
	return(ll)
}


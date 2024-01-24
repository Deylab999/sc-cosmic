

tabb1 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/POPS-plus/POPS_0kb_qmatched_MAGMA.txt")
tabb2 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/PolyGene_plus/PolyGene_PPI_PoPS_justMAGMA_Mean.txt")
tabb3 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/MAGMA-v108/MAGMA_v108_GENE_0_ZSTAT.txt")
tabb3[is.na(tabb3)] = 0
tabb3[tabb3 < 0] = 0
common_genes = Reduce(intersect, list(rownames(tabb1), rownames(tabb2), rownames(tabb3)))

tabb1 = tabb1[match(common_genes, rownames(tabb1)), ]
tabb2 = tabb2[match(common_genes, rownames(tabb2)), ]
tabb3 = tabb3[match(common_genes, rownames(tabb3)), ]

NUMGENES=500
knn_idx = 3:7
traitname = "PASS_IBD_deLange2017"

select_cell_types = c("T_Lymphocytes")

enrlist = list()

for(counter1 in 1:length(select_cell_types)){
	tabb_noncausal = read.delim(paste0("/n/groups/price/kushal/Causal/output/ucdata2/Noncausal_VIDE_uc_poiss_", select_cell_types[counter1], "_Healthy_Inflamed.txt"))
	
	tabb_causal1 = read.delim(paste0("/n/groups/price/kushal/Causal/output/ucdata2/Causal_VIDE_uc_poiss_", select_cell_types[counter1], "_Healthy_Inflamed_propensity_knn.txt"))
	tabb_causal1 = tabb_causal1[which(tabb_causal1$Gene %in% common_genes), ]
	pvals1 = apply(tabb_causal1[,knn_idx], 1, function(x) {
		y = x
		y[y > 0.99] = 0.99
		return(CCT(y))})

	tabb_causal2 = read.delim(paste0("/n/groups/price/kushal/Causal/output/ucdata2/Causal_VIDE_uc_poiss_", select_cell_types[counter1], "_Healthy_Inflamed_svd_knn.txt"))
	tabb_causal2 = tabb_causal2[which(tabb_causal2$Gene %in% common_genes), ]
	pvals2 = apply(tabb_causal2[,knn_idx], 1, function(x) {
		y = x
		y[y > 0.99] = 0.99
		return(CCT(y))})

	tabb_causal3 = read.delim(paste0("/n/groups/price/kushal/Causal/output/ucdata2/Causal_VIDE_uc_poiss_", select_cell_types[counter1], "_Healthy_Inflamed_propensity_svd_knn.txt"))
	tabb_causal3 = tabb_causal3[which(tabb_causal3$Gene %in% common_genes), ]
	pvals3 = apply(tabb_causal3[,knn_idx], 1, function(x) {
		y = x
		y[y > 0.99] = 0.99
		return(CCT(y))})

	tabb_causal4 = read.delim(paste0("/n/groups/price/kushal/Causal/output/ucdata2/Causal_VIDE_uc_poiss_", select_cell_types[counter1], "_Healthy_Inflamed_svd_knn_smooth.txt"))
	tabb_causal4 = tabb_causal4[which(tabb_causal4$Gene %in% common_genes), ]
	pvals4 = apply(tabb_causal4[,knn_idx], 1, function(x) {
		y = x
		y[y > 0.99] = 0.99
		return(CCT(y))})

	tabb_causal5 = read.delim(paste0("/n/groups/price/kushal/Causal/output/ucdata2/Causal_VIDE_uc_poiss_", select_cell_types[counter1], "_Healthy_Inflamed_propensity_svd_knn_smooth.txt"))
	tabb_causal5 = tabb_causal5[which(tabb_causal5$Gene %in% common_genes), ]
	pvals5 = apply(tabb_causal5[,knn_idx], 1, function(x) {
		y = x
		y[y > 0.99] = 0.99
		return(CCT(y))})

	# genes_causal1 = tabb_causal1$Gene[which(pvals1 < 0.05)]
	# genes_causal2 = tabb_causal2$Gene[which(pvals2 < 0.05)]
	# genes_causal3 = tabb_causal3$Gene[which(pvals3 < 0.05)]
	# genes_causal4 = tabb_causal4$Gene[which(pvals4 < 0.05)]
	# genes_causal5 = tabb_causal5$Gene[which(pvals5 < 0.05)]

	genes_causal1 = tabb_causal1$Gene[order(pvals1, decreasing=F)[1:NUMGENES]]
	genes_causal2 = tabb_causal2$Gene[order(pvals2, decreasing=F)[1:NUMGENES]]
	genes_causal3 = tabb_causal3$Gene[order(pvals3, decreasing=F)[1:NUMGENES]]
	genes_causal4 = tabb_causal4$Gene[order(pvals4, decreasing=F)[1:NUMGENES]]
	genes_causal5 = tabb_causal5$Gene[order(pvals5, decreasing=F)[1:NUMGENES]]

	common_genes2 = Reduce(intersect, list(tabb_causal1$Gene, tabb_causal2$Gene, tabb_causal3$Gene,
		tabb_causal4$Gene, tabb_causal5$Gene))
	pvals_combo = apply(cbind(pvals1[match(common_genes2, tabb_causal1$Gene)], 
		                      pvals2[match(common_genes2, tabb_causal2$Gene)], 
		                      pvals3[match(common_genes2, tabb_causal3$Gene)], 
		                      pvals4[match(common_genes2, tabb_causal4$Gene)], 
		                      pvals5[match(common_genes2, tabb_causal5$Gene)]), 1, function(x) return(CCT(x)))
	#genes_causal6 = common_genes2[which(pvals_combo < 0.05)]
	genes_causal6 = common_genes2[order(pvals_combo, decreasing=F)[1:NUMGENES]]



	tabb_noncausal = tabb_noncausal[which(tabb_noncausal$Gene %in% common_genes), ]
	#pval_merge = apply(tabb_causal[,3:7], 1, function(x) return(CCT(x)))
	#pval_merge = tabb_causal[,3]
	#genes_causal = tabb_causal$Gene[order(pval_merge, decreasing=F)[1:100]]
	#genes_noncausal = tabb_noncausal$Gene[which(tabb_noncausal$p.value < 0.05)]
	genes_noncausal = tabb_noncausal$Gene[order(tabb_noncausal$p.value, decreasing=F)[1:NUMGENES]]

	enr1 = mean(tabb1[ match(genes_causal1, common_genes), traitname])/mean(tabb1[, traitname])
	enr2 = mean(tabb2[ match(genes_causal1, common_genes), traitname])/mean(tabb2[, traitname])
	enr3 = mean(tabb3[ match(genes_causal1, common_genes), traitname])/mean(tabb3[, traitname])

	enr4 = mean(tabb1[ match(genes_causal2, common_genes), traitname])/mean(tabb1[, traitname])
	enr5 = mean(tabb2[ match(genes_causal2, common_genes), traitname])/mean(tabb2[, traitname])
	enr6 = mean(tabb3[ match(genes_causal2, common_genes), traitname])/mean(tabb3[, traitname])

	enr7 = mean(tabb1[ match(genes_causal3, common_genes), traitname])/mean(tabb1[, traitname])
	enr8 = mean(tabb2[ match(genes_causal3, common_genes), traitname])/mean(tabb2[, traitname])
	enr9 = mean(tabb3[ match(genes_causal3, common_genes), traitname])/mean(tabb3[, traitname])

	enr10 = mean(tabb1[ match(genes_causal4, common_genes), traitname])/mean(tabb1[, traitname])
	enr11 = mean(tabb2[ match(genes_causal4, common_genes), traitname])/mean(tabb2[, traitname])
	enr12 = mean(tabb3[ match(genes_causal4, common_genes), traitname])/mean(tabb3[, traitname])

	enr13 = mean(tabb1[ match(genes_causal5, common_genes), traitname])/mean(tabb1[, traitname])
	enr14 = mean(tabb2[ match(genes_causal5, common_genes), traitname])/mean(tabb2[, traitname])
	enr15 = mean(tabb3[ match(genes_causal5, common_genes), traitname])/mean(tabb3[, traitname])

	enr16 = mean(tabb1[ match(genes_causal6, common_genes), traitname])/mean(tabb1[, traitname])
	enr17 = mean(tabb2[ match(genes_causal6, common_genes), traitname])/mean(tabb2[, traitname])
	enr18 = mean(tabb3[ match(genes_causal6, common_genes), traitname])/mean(tabb3[, traitname])

	enr19 = mean(tabb1[ match(genes_noncausal, common_genes), traitname])/mean(tabb1[, traitname])
	enr20 = mean(tabb2[ match(genes_noncausal, common_genes), traitname])/mean(tabb2[, traitname])
	enr21 = mean(tabb3[ match(genes_noncausal, common_genes), traitname])/mean(tabb3[, traitname])

	enrlist[[counter1]] = rbind(c(enr1, enr2, enr3),
		c(enr4, enr5, enr6),
		c(enr7, enr8, enr9),
		c(enr10, enr11, enr12),
		c(enr13, enr14, enr15),
		c(enr16, enr17, enr18),
		c(enr19, enr20, enr21))
	cat("We are at cell type:", counter1, "\n")
}

rownames(enrlist[[counter1]]) = c("Propensity_knn", "SVD_knn", "Propensity_SVD_knn", "SVD_knn_rs", "Propensity_SVD_knn_rs", "Combo.causal", "Noncausal")
colnames(enrlist[[counter1]]) = c("PoPS", "PolyGene", "MAGMA")
enrlist[[counter1]]





counter1=1
tabb_causal1 = read.delim(paste0("/n/groups/price/kushal/Causal/output/Causal_VIDE_uc_poiss_", select_cell_types[counter1], "_Healthy_Inflamed_propensity_knn.txt"))
tabb_causal2 = read.delim(paste0("/n/groups/price/kushal/Causal/output/Causal_VIDE_uc_poiss_", select_cell_types[counter1], "_Healthy_Inflamed_svd_knn.txt"))
tabb_causal3 = read.delim(paste0("/n/groups/price/kushal/Causal/output/Causal_VIDE_uc_poiss_", select_cell_types[counter1], "_Healthy_Inflamed_propensity_svd_knn.txt"))



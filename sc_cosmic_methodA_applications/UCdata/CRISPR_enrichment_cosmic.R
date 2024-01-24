

genes1 = read.table("/n/groups/price/kushal/extras/essentiality/core_essential.txt")[,1]
genes2 = read.table("/n/groups/price/kushal/extras/essentiality/core_essential_shrna.txt")[,1]
genes3 = read.table("/n/groups/price/kushal/extras/essentiality/mgi_essential.txt")[,1]

knn_idx = 3:7
traitname = "PASS_IBD_deLange2017"
NUMGENES=1000
select_cell_types = c("TA")

enrlist = list()

for(counter1 in 1:length(select_cell_types)){
	tabb_noncausal = read.delim(paste0("/n/groups/price/kushal/Causal/output/ucdata2/Noncausal_VIDE_uc_poiss_", select_cell_types[counter1], "_Healthy_Inflamed.txt"))
	
	tabb_causal1 = read.delim(paste0("/n/groups/price/kushal/Causal/output/ucdata2/Causal_VIDE_uc_poiss_", select_cell_types[counter1], "_Healthy_Inflamed_propensity_knn.txt"))
	pvals1 = apply(tabb_causal1[,knn_idx], 1, function(x) {
		y = x
		y[y > 0.99] = 0.99
		return(CCT(y))})

	tabb_causal2 = read.delim(paste0("/n/groups/price/kushal/Causal/output/ucdata2/Causal_VIDE_uc_poiss_", select_cell_types[counter1], "_Healthy_Inflamed_svd_knn.txt"))
	pvals2 = apply(tabb_causal2[,knn_idx], 1, function(x) {
		y = x
		y[y > 0.99] = 0.99
		return(CCT(y))})

	tabb_causal3 = read.delim(paste0("/n/groups/price/kushal/Causal/output/ucdata2/Causal_VIDE_uc_poiss_", select_cell_types[counter1], "_Healthy_Inflamed_propensity_svd_knn.txt"))
	pvals3 = apply(tabb_causal3[,knn_idx], 1, function(x) {
		y = x
		y[y > 0.99] = 0.99
		return(CCT(y))})

	tabb_causal4 = read.delim(paste0("/n/groups/price/kushal/Causal/output/ucdata2/Causal_VIDE_uc_poiss_", select_cell_types[counter1], "_Healthy_Inflamed_svd_knn_smooth.txt"))
	pvals4 = apply(tabb_causal4[,knn_idx], 1, function(x) {
		y = x
		y[y > 0.99] = 0.99
		return(CCT(y))})

	tabb_causal5 = read.delim(paste0("/n/groups/price/kushal/Causal/output/ucdata2/Causal_VIDE_uc_poiss_", select_cell_types[counter1], "_Healthy_Inflamed_propensity_svd_knn_smooth.txt"))
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

	#pval_merge = apply(tabb_causal[,3:7], 1, function(x) return(CCT(x)))
	#pval_merge = tabb_causal[,3]
	#genes_causal = tabb_causal$Gene[order(pval_merge, decreasing=F)[1:100]]
	#genes_noncausal = tabb_noncausal$Gene[which(tabb_noncausal$p.value < 0.05)]
	genes_noncausal = tabb_noncausal$Gene[order(tabb_noncausal$p.value, decreasing=F)[1:NUMGENES]]

	enr1 = (length(intersect(genes_causal1, genes1))*21500)/(length(genes_causal1)*length(genes1))
	enr2 = (length(intersect(genes_causal1, genes2))*21500)/(length(genes_causal1)*length(genes2))
	enr3 = (length(intersect(genes_causal1, genes3))*21500)/(length(genes_causal1)*length(genes3))

	enr4 = (length(intersect(genes_causal2, genes1))*21500)/(length(genes_causal2)*length(genes1))
	enr5 = (length(intersect(genes_causal2, genes2))*21500)/(length(genes_causal2)*length(genes2))
	enr6 = (length(intersect(genes_causal2, genes3))*21500)/(length(genes_causal2)*length(genes3))

	enr7 = (length(intersect(genes_causal3, genes1))*21500)/(length(genes_causal3)*length(genes1))
	enr8 = (length(intersect(genes_causal3, genes2))*21500)/(length(genes_causal3)*length(genes2))
	enr9 = (length(intersect(genes_causal3, genes3))*21500)/(length(genes_causal3)*length(genes3))

	enr10 = (length(intersect(genes_causal4, genes1))*21500)/(length(genes_causal4)*length(genes1))
	enr11 = (length(intersect(genes_causal4, genes2))*21500)/(length(genes_causal4)*length(genes2))
	enr12 = (length(intersect(genes_causal4, genes3))*21500)/(length(genes_causal4)*length(genes3))

	enr13 = (length(intersect(genes_causal5, genes1))*21500)/(length(genes_causal5)*length(genes1))
	enr14 = (length(intersect(genes_causal5, genes2))*21500)/(length(genes_causal5)*length(genes2))
	enr15 = (length(intersect(genes_causal5, genes3))*21500)/(length(genes_causal5)*length(genes3))

	enr16 = (length(intersect(genes_causal6, genes1))*21500)/(length(genes_causal6)*length(genes1))
	enr17 = (length(intersect(genes_causal6, genes2))*21500)/(length(genes_causal6)*length(genes2))
	enr18 = (length(intersect(genes_causal6, genes3))*21500)/(length(genes_causal6)*length(genes3))

	enr19 = (length(intersect(genes_noncausal, genes1))*21500)/(length(genes_noncausal)*length(genes1))
	enr20 = (length(intersect(genes_noncausal, genes2))*21500)/(length(genes_noncausal)*length(genes2))
	enr21 = (length(intersect(genes_noncausal, genes3))*21500)/(length(genes_noncausal)*length(genes3))

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
colnames(enrlist[[counter1]]) = c("Human-crispr", "Human-crispr-shRNA", "Mouse-crispr")
enrlist[[counter1]]

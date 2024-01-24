

celltypes = c("B_Lymphocytes", "Endothelial",  "Enterocytes",   "Fibroblast",
"Goblet", "Macrophages", "Mast", "T_Lymphocytes",  "TA")

methods = c("propensity_svd_knn", "propensity_knn",
            "svd_knn", "svd_knn_smooth", 
 	        "propensity_svd_knn_smooth")


outt1 = cbind.data.frame("Healthy", "Inflamed", rep(celltypes, length(methods)), rep(methods, each=length(celltypes)),
	            "/n/groups/price/kushal/Causal/output/ucdata2", "uc_poiss")
colnames(outt1) = c("context1", "context2", "Celltype", "method", "output_cell", "output_name")

outt2 = cbind.data.frame("Healthy", "Non-inflamed", rep(celltypes, length(methods)), rep(methods, each=length(celltypes)),
	            "/n/groups/price/kushal/Causal/output/ucdata2", "uc_poiss")
colnames(outt2) = c("context1", "context2", "Celltype", "method", "output_cell", "output_name")


combo_outt = rbind(outt1, outt2)
write.table(combo_outt, file = "/n/groups/price/kushal/Causal/data/batch_joblist2.txt",
	col.names=F, row.names=F, sep = "\t", quote=F)


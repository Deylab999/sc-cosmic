
celltypes = c("B_Lymphocytes", "Endothelial",  "Enterocytes",   "Fibroblast",
"Goblet", "Macrophages", "Mast", "T_Lymphocytes",  "TA")

methods = c("propensity_svd_knn", "propensity_knn",
            "svd_knn")

## context_label="health"
outt1 = cbind.data.frame("/n/groups/price/kushal/Causal/codes/UCdata/run_batch_cosmic_ucdata_gauss.R", "Healthy", "UC", rep(celltypes, length(methods)), rep(methods, each=length(celltypes)),
	            "/n/groups/price/kushal/Causal/output/ucdata_gauss", "uc_gauss")
colnames(outt1) = c("code", "context1", "context2", "Celltype", "method", "output_cell", "output_name")

outt2 = cbind.data.frame("/n/groups/price/kushal/Causal/codes/UCdata/run_batch_cosmic_ucdata_poiss.R", "Healthy", "UC", rep(celltypes, length(methods)), rep(methods, each=length(celltypes)),
	            "/n/groups/price/kushal/Causal/output/ucdata_poiss", "uc_poiss")
colnames(outt2) = c("code", "context1", "context2", "Celltype", "method", "output_cell", "output_name")

combo_outt = rbind(outt1, outt2)
write.table(combo_outt, file = "/n/groups/price/kushal/Causal/data/batch_joblist1.txt",
	col.names=F, row.names=F, sep = "\t", quote=F)

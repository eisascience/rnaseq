#Generates most variable matrix. voom_mat = log2 CPM counts, CPM_min = filter for minimum mean, n = max number of variable genes)
get_most_variable_matrix <- function(voom_mat, CPM_min, n)
{
	vars <- apply(voom_mat, 1, var);
	means <- apply(voom_mat, 1, mean);
	top_vars <- order(vars, decreasing = TRUE)[which(order(vars, decreasing = TRUE) %in% which(means > CPM_min))][1:n];
	return(list(voom_mat[top_vars,], top_vars));
}

#Generates WGCNA object. mat = matrix of log2 CPM counts (usually from above function)
get_net <- function(mat){
	power <- 1; 
	sft <- pickSoftThreshold(t(mat), networkType = "signed", powerVector = seq(1, 30, 1));
	pv <- which(sft$fitIndices[,2] >= .9)
	if(length(pv) > 1)
	{
		power <- sft$fitIndices[pv[1],1];
	}
	else
	{
		power <- sft$fitIndices[which(sft$fitIndices[,2] == max(sft$fitIndices[,2]))[1],1];
	}
	print(power);
	net <- blockwiseModules(t(mat), networkType = "signed", power = power, numericLabels = TRUE, minModuleSize = 100);
	return(list(net, sft));
}

#Generates eigengene object and boxplots. net = WGCNA object, mat = log2 CPM counts matrix, fname = jpeg name, groups = numeric sample groups, labels = group labels
generate_eigengene_plot <- function(net, mat, fname, groups, labels){
	mat_eigens <- moduleEigengenes(t(mat), net$colors)$eigengene;
	i <- 2;
	jpeg(fname, 1600, 1600, res = 300)
	par(mfrow = c(floor(sqrt(dim(mat_eigens)[2] - 1)), ceiling(sqrt(dim(mat_eigens)[2] - 1))));
	while(i <= max(net$colors) + 1)
	{
		to_plot <- list();
		j <- 1; 
		while(j <= max(groups))
		{
			to_plot[[j]] <- mat_eigens[which(groups == j),i];
			j <- j + 1;
		}
		boxplot(to_plot, main = paste("ME ", i - 1, sep = ''), axes = FALSE);
		axis(2);
		axis(1, 1:max(groups), labels);
		i <- i + 1;
	}
	dev.off();
	return(mat_eigens);
}

#Automates network construction based on top 5000 most variables genes. Pass a DGEList object (e.g. from edgeR scripts) as raw_mat for best results
run_WGCNA_top_5000 <- function(raw_mat, minCPM = 1, fname, groups, labels){
	voom_mat <- voom(raw_mat);
	most_var <- get_most_variable_matrix(voom_mat$E, minCPM, 5000);
	net_list <- get_net(most_var[[1]]);
	mat_eigens <- generate_eigengene_plot(net_list[[1]], most_var[[1]], fname, groups, labels);
	return(list(voom_mat, most_var, net_list, mat_eigens));
}
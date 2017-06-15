clinical_data <- read.table(pipe("awk 'NR>8' clinical_data.txt | awk -F \"\t\" '{print $2 \"\t\" $5}'"), fill =  TRUE, row.names = 1)
gene_data <- read.table(pipe("awk 'NR>16' GPL96-15653.txt | awk -F \"\t\" '{print $1 \"\t\" $11}'|sed 's/ //g' "), header = TRUE, fill = TRUE, row.names=1)
probe_data <- read.table(pipe("awk 'NR>287' GSE2034-22071.txt"), fill = TRUE, header = TRUE, row.names = 1)

clinical_data <- t(clinical_data)
# probe_data <- rbind(probe_data, nrow(probe_data) + 1)
probe_data <- rbind(rep(NA, ncol(probe_data)),probe_data)
probe_data <- cbind(rep(NA, nrow(probe_data)), probe_data)

rownames(probe_data)[1] <- "Relapse"
colnames(probe_data)[1] <- "Gene"

for ( i in 1 : ncol(clinical_data) ) {
	probe_data[1,colnames(clinical_data)[i]] <- clinical_data[1,colnames(clinical_data)[i]]
}

for ( i in 1 : nrow(gene_data) ) {
	probe_data[rownames(gene_data)[i],1] <- toString(gene_data[rownames(gene_data)[i],1])
}

write.table(probe_data, "wang_data.txt", row.names = TRUE, col.names = NA, sep="\t")

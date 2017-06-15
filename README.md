####The paper that originally published this dataset: http://www.sciencedirect.com/science/article/pii/S0140673605179471
####This directory is to apply statistical methods to identify differentially expressed genes between 179 lymph-node negative metastasis-free patients and 107 lymph-node negative patients that developed a subsequent metastasis.
####The Affymetrix Human Genome U133A Array is explored to find genetic markers associated with breast cancer metastasis. The detailed description of the microarray data and the purpose of the analysis can be found in the original paper.
####The raw data was obtained from the Gene Expression Omnibus (GEO) database
####(http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2034) and saved here as tab delimited GSE2034-22071.txt
####Clinical information associated with the patients here from http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=GSE2034&id=40089&db=GeoDb_blob26 saved as clinical_data.txt
####This links each patient to the GEO accession associated with each array and includes several pieces of information about each patient, including the relapse status.
####GPL96-39578.txt is the mapping between probes on the array and the corresponding genes in the human genome. This information for the Affymetrix Human Genome U133A Array was downloaded from http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL96

####import_wang_data.r was provided by Chad L Meyers and used to merge all of these files together into a single matrix
####to form the final tab-delimited data file wang_data.txt with sample IDs across the first row, the relapse status for each patient across the second row, and the probe values in the following rows. Each probe row contains the Affymetrix probe ID in the first column, and the corresponding human gene name in the second column.


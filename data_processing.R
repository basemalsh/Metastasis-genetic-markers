data = read.table('wang_data.txt', header = TRUE)
matrix = as.matrix(data)

#number of probes included in the dataset
number_of_probes= nrow(matrix) - 1 
#22283

#number of patient samples in the dataset
number_of_patients= ncol(matrix) - 2 
#286

#divide patients into two groups based on relapse status
relapse<- matrix[1,3:288]
relapsenum<- as.numeric(relapse)

#patients who had relapses (relapse=1)
relapse_patients = length(which(relapsenum==1.0)) #107
#patients relapse free (relapse=0, i.e. no metastases)
norelapse_patients= length(which(relapsenum==0.0)) #179

#Number of unique genes represented by the probes in the dataset
#In the future it would be better to average multiple probes mapping to same gene
#but here they are analyzed independently
repeated= length(which(duplicated(matrix[,2])))
unique = 22283- repeated #13212

#histogram of the data
all_data <- matrix[,3:288]
all_data <- as.numeric(all_data)
hist(all_data, breaks = 100) # distribution seems positively skewed

#Any values <= 0 replaced with a value of 1
all_data[all_data <= 0] <- 1
#log-transform (log 10) the entire data matrix & plot the log-transformed data
hist(log10(all_data), breaks = 100) # distribution seems relatively symmetrical

#Plot histograms for the first four arrays
first<- matrix[2:22284,3]
first<- as.numeric(first)
firstlog<- log10(first)
hist(firstlog, breaks = 100) # distribution still symmetrical

second<- matrix[,4]
second<-as.numeric(second)
secondlog<- log10(second)
hist(secondlog, breaks=100) # distribution still symmetrical

third<- matrix[,5]
third<- as.numeric(third)
thirdlog<- log10(third)
hist(thirdlog, breaks = 100) # distribution still symmetrical

fourth<- matrix[,6]
fourth<- as.numeric(fourth)
fourthlog<- log10(fourth)
hist(fourthlog, breaks = 100) # distribution still symmetrical

first_four<- matrix[,3:6]
first_four <- as.numeric(first_four)
firstfourlog<- log(first_four)
hist(firstfourlog, breaks = 100) # distribution still symmetrical

#Perform quantile normalization across all arrays (samples) such that each has the same empirical distribution

source("http://bioconductor.org/biocLite.R")
biocLite("preprocessCore")
library("preprocessCore")

copyfornormalize <- matrix[2:22284,3:288]
class(copyfornormalize)<- "numeric"
#is.numeric(copyfornormalize)
#copyfornormalize[copyfornormalize <= 0] <- 1
logdata<- log10(copyfornormalize)
normalized<- normalize.quantiles.use.target(logdata,firstlog,copy=TRUE,subset = NULL)

hist(normalized[,1], breaks = 100)
hist(normalized[,2], breaks = 100)
hist(normalized[,3], breaks = 100)
hist(normalized[,4], breaks = 100)

#Using the t-test and Wilcoxon rank-sum statistics to identify differentially expressed probes with a per-probe significance level of p < 0.05

#You should divide the gene expression data into two groups (metastasis vs. non-metastasis), 
#using the relapse variable in the clinical data, and test each probe independently
met<-which(relapsenum==1.0)
nonmet<-which(relapsenum==0.0)
metastasized_data<- normalized[,met]
nonmetastasized_data<- normalized[,nonmet]

df1<-as.data.frame(metastasized_data)
df0<-as.data.frame(nonmetastasized_data)

#ttest method 1

getPvalueRow<-function(rownumber){
  x<-as.numeric(df0[rownumber,])
  y<-as.numeric(df1[rownumber,])
  t<-t.test(x,y)
  df<- data.frame(genenum=rownumber, pvalue=t$p.value)
  return(df)
}

dfnew<-sapply(1:22283, getPvalueRow)
#df<-t(dfnew)

which.min(dfnew[2,]) #1852
dfnew[2,1852] #4.044459e-07
matrix[1853,2] #ACBD3

order(dfnew)

#ttest method 2
Pvalues<- apply(normalized,1, function(x) t.test(x[met],y=x[nonmet],var.equal = TRUE)$p.value)

ordered<-order(Pvalues)
ordered[1:10] #1852  8873 18842 21437  1297   897   706 14228   744  4168   
Pvalues[ordered[1:10]] #4.662076e-07 2.247677e-06 2.421236e-06 3.491345e-06 4.841266e-06 5.984451e-06 6.839090e-06  8.723646e-06 1.007707e-05 1.095083e-05
matrix[ordered[1:10]+1,2] #"ACBD3"   "ABCC5"   "WFDC1"   "RACGAP1" "CLINT1"  "ZFP36L2" "FBXO7"   "SHC1"    "ERP29"   "NEK2" 
matrix[ordered[1:10]+1,1] #"202324_s_at" "209380_s_at" "219478_at"   "222077_s_at" "201769_at"   "201369_s_at"  "201178_at"   "214853_s_at" "201216_at"   "204641_at" 

#4a
#The protein encoded by acyl-CoA binding domain containing 3 (ACBD3) gene 
#is involved in the maintenance of Golgi structure and function through its 
#interaction with the integral membrane protein giantin. It may also be 
#involved in the hormonal regulation of steroid formation. It is 
#among the top 10% under-expressed in liver cancer. It reacts with 
#PI4KB which also reacts with  Metastasis-associated protein MTA2. 
#ACBD3 may participate in PKARIA-mediated tumorigenesis and hormone-independent 
#hypercortisolism. Microarray data from published cancer-related databases 
#has shown that ACBD3 is involved in cell cycle control and is regulated by 
#retinoblastoma protein (pRB), a critical tumor suppressor.

#wilcox 

WPvalues<- apply(normalized,1, function(x) wilcox.test(x[met],y=x[nonmet],var.equal = TRUE)$p.value)
Wordered<-order(WPvalues)
Wordered[1:10] #1852 18842 21761  1297   897 21437  4168 14228 18065 12285
WPvalues[Wordered[1:10]] #7.349612e-07 1.780967e-06 2.411102e-06 3.634223e-06 4.676169e-06 7.092196e-06 7.166352e-06  8.143665e-06 8.172226e-06 9.060792e-06
matrix[Wordered[1:10]+1,2] #"ACBD3"   "WFDC1"   "BLZF1"   "CLINT1"  "ZFP36L2" "RACGAP1" "NEK2"    "SHC1"    "LACTB2"  "SEC24A" 
matrix[Wordered[1:10]+1,1] #"202324_s_at" "219478_at"   "32088_at"    "201769_at"   "201369_s_at" "222077_s_at"  "204641_at"   "214853_s_at" "218701_at"   "212900_at"  


#The number of selected probes for each test + plotting the histogram of the log-transformed (negative log10) p-values of all probes

log_transformed_TTest_P <- -log10(Pvalues)
hist(log_transformed_TTest_P , breaks = 100)
log_transformed_WTest_P <- -log10(WPvalues)
hist(log_transformed_WTest_P , breaks = 100)

#overlap between the set of probes deemed significant by the two different approaches

isDifferent<-function(rownumber){
  x<-as.numeric(df0[rownumber,])
  y<-as.numeric(df1[rownumber,])
  t<-t.test(x,y,var.equal = TRUE)
  return (t$p.value < 0.05)
}

ttestsignificantprobes<-sapply(1:22283, isDifferent)
Tsignificants<- which(ttestsignificantprobes== TRUE)
length(Tsignificants) #3250

WisDifferent<-function(rownumber){
  x<-as.numeric(df0[rownumber,])
  y<-as.numeric(df1[rownumber,])
  t<-wilcox.test(x,y,var.equal = TRUE)
  return (t$p.value < 0.05)
}
Wilcoxtestsignificantprobes<-sapply(1:22283, WisDifferent)
Wsignificants<- which(Wilcoxtestsignificantprobes== TRUE)
length(Wsignificants) #3290

length(intersect(Tsignificants,Wsignificants)) #2700

#Multiple hypothesis correction
#Bonferroni Correction

#Take the p-value of each gene and multiply it by the number of genes in the
#gene list. If the corrected p-value is still below the cutoff, the gene will be
#significant:
#  BC Corrected P-value= p-value * n (number of genes in test) <0.05 
Corrected <- WPvalues * 22283
which(Corrected<0.05) #1852 18842

#or
#New threshold = 0.05/22283 = 2.24386303e-6
#Bonferroni corrected ?? for these data is 0.05/22283 = 2.24386303e-6.
#If we had used this family-wise error rate in our individual hypothesis tests,
#we would have concluded only 2 results were significant (shown below, and above)

BCWisDifferent<-function(rownumber){
  x<-as.numeric(df0[rownumber,])
  y<-as.numeric(df1[rownumber,])
  t<-wilcox.test(x,y,var.equal = TRUE)
  return (t$p.value < 2.24386303e-6)
}
BCWilcoxtestsignificantprobes<-sapply(1:22283, BCWisDifferent)
BCWsignificants<- which(BCWilcoxtestsignificantprobes== TRUE)
length(BCWsignificants) #2


#the number of selected probes by the BH procedure

#pvalues entered from smallest to largest
pvals<- c(WPvalues[Wordered])
#enter the target FDR
FDR<- .05
#calculate the threshold for each p-value
threshold<- FDR*(1:length(pvals))/length(pvals)
#compare the p-value against its threshold and display results
compared<- cbind(pvals,threshold,pvals<=threshold)
which(compared[1:length(pvals),3]==1) #Starting with the smallest p-value and moving up 
#we find that the largest k for which the p-value is less than its threshold is k = 176.
#assuming independence of the genes

#Rank the probes by their p-values in ascending order and plot the p-values and the threshold used for adjusted p-values as a function of the index of the ranked probes. More specifically, with the x-axis as the index of the ranked probes from 1 to the number of probes, plot the first 500 probes (i on the x-axis, the p-value of the ith probe on the y-axis for each probe). As a separate color, plot the adjusted p-value threshold at each point (i on the x-axis, threshold for adjusted p-value on the y-axis), where i is the index of the genes and alpha=0.05
p_values<-compared[1:500]
plot(p_values, col="red")
lines(threshold[1:500],col="blue")
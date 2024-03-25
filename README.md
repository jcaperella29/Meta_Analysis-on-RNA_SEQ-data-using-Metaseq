# Meta_Analysis-on-RNA_SEQ-data-using-Metaseq
A script that performs Meta_analysis on RNA_seq data from multiple studies
#First we import the needed libraries .We then provide information on the studies  and phenotype. 
#Then we read in the sample data "Breasts cancer and attach the above information.
#Note , do perform this analysis with their own data , users can provide a dataframe with the counts data from  studies that they are pooling from as well as related phenotype and studies data.
#then we aquire the results of performing Onesided-NOISeq on the  using "data(Result.Meta)" 
#We then asssign the data to a the variable results
#Next we use the Fisher and Stouffer tests to detrimine the likehood that genes are up or down regulated with respect to the phenotype.
#Then we extract the genes most likely to up and down regulated repectively and save them as text files for futher analysis.
#Finally, we pull both of the above sets of genes togather and output that as a text files as well.
 

library("metaSeq")
library("snow")

 data(BreastCancer)

 
PhenoType <- c(1,1,1,0,0, 1,0, 1,1,1,1,1,1,1,0, 1,1,0)

study <- c("A","A","A","A","A", "B","B", "C","C","C","C","C","C","C","C", "D","D","D")

cds <- meta.readData(data = BreastCancer, factor = PhenoType, studies = study)
data(Result.Meta)
result <- Result.Meta

F <- Fisher.test(result)
 S <- Stouffer.test(result)
# These outputs are summalized as list whose length is 3. First member is the probability which means a gene is upper-regulated genes, and Second member is lower-regulated genes.
#Weight in each study is also saved as its third member (weight is used only by Stouer's
#method).
 
 #organizing results
 #gathering the likely up-regulated genes
 Upp_prob_by_Fisher<-data.frame(F$Upper)
 Upp_prob_by_Fisher$gene_names<-row.names(Upp_prob_by_Fisher)
Upp_prob_by_Fisher$test<-rep(list("Fisher"),23368)

 Upp_prob_by_Stouffer<-data.frame(S$Upper) 
 Upp_prob_by_Stouffer$gene_names<-row.names(Upp_prob_by_Stouffer)
 Upp_prob_by_Stouffer$test<-rep(list("Stouffer"),23368)

Upp_probDF<-merge(x= Upp_prob_by_Fisher,y = Upp_prob_by_Stouffer,by = "gene_names") 


Upp_probDF<-Upp_probDF[order(-Upp_probDF$F.Upper,-Upp_probDF$S.Upper),]

#these are the twenty genes that are most likely to be up-regulated
Top_twenty_upreg<-Upp_probDF[1:20,]
Top_twenty_upreg_genes<-Upp_probDF[1:20,"gene_names"]
write.table(Top_twenty_upreg_genes,file = "C:/Users/ccape/Downloads/up_regulated.txt")

#gathering the likely down-regulated genes
Lower_prob_by_Fisher<-data.frame(F$Lower)
Lower_prob_by_Fisher$gene_names<-row.names(Lower_prob_by_Fisher)
Lower_prob_by_Fisher$test<-rep(list("Fisher"),23368)

Lower_prob_by_Stouffer<-data.frame(S$Lower) 
Lower_prob_by_Stouffer$gene_names<-row.names(Lower_prob_by_Stouffer)
Lower_prob_by_Stouffer$test<-rep(list("Stouffer"),23368)

Lower_probDF<-merge(x= Lower_prob_by_Fisher,y = Lower_prob_by_Stouffer,by = "gene_names") 


Lower_probDF<-Lower_probDF[order(-Lower_probDF$F.Lower,-Lower_probDF$S.Lower),]

Top_twenty_down_reg<-Lower_probDF[1:20,]
Top_twenty_down_reg_genes<-Lower_probDF[1:20,"gene_names"]
write.table(Top_twenty_down_reg_genes,file = "C:/Users/ccape/Downloads/down_regulated.txt")
#lastly we group all differentially_expressed Genes
differentially_expressed_genes<-c(Top_twenty_upreg_genes,Top_twenty_down_reg_genes)

write.table(differentially_expressed_genes,file = "C:/Users/ccape/Downloads/differentially_expressed.txt")



# this was built off of https://www.bioconductor.org/packages/devel/bioc/vignettes/metaSeq/inst/doc/metaSeq.pdf.
#Added lines  for reporting results.

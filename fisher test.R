# Load canonical transcipts list ----
library(readr)
setwd("~/Documents/Research/all chrom")
stringsAsFactors= FALSE
canonical_transcripts<- read_csv("canonical transcript list.csv", col_names = FALSE)
colnames(canonical_transcripts)<- c("Chrom", "Start", "End", "Number", "Transcript", "Gene" )

m<- unlist(strsplit(as.character(canonical_transcripts$Transcript),".",fixed = TRUE))
canonical_transcripts$Transcript<- m[seq(1, length(m), 2)]

m<- unlist(strsplit(as.character(canonical_transcripts$Gene),".",fixed = TRUE))
canonical_transcripts$Gene<- m[seq(1, length(m), 2)]

# Autosomes ====
setwd("~/Documents/Research/all chrom/all_var_analysis/singletones")

var_synonymous <- read_csv("synonymous/new/all_chrom_new.csv")
var_synonymous<- var_synonymous[,-c(14:15)]

var_misssense <- read_csv("missense/new/all_chrom_new.csv")
var_misssense<- var_misssense[,-c(14:15)]

var_lof <- read_csv("lof/new/all_chrom_new.csv")
var_lof<- var_lof[,-c(14:15)]

edit_df<- function(df) {
  colnames(df)[2]<- "transcript"
  df<- df[which(df$transcript %in% canonical_transcripts$Transcript),]
  df<- df[,c(-1, -5, -6,-7)]
}

var_synonymous<- edit_df(var_synonymous)
var_misssense<- edit_df(var_misssense)
var_lof<- edit_df(var_lof)

#convet transcripts to gene symbol
gene_converter <- read_delim("~/Documents/gene_converter",  "\t", escape_double = FALSE, trim_ws = TRUE)

convert_genes<- function(df){
  m<- match(df$transcript, gene_converter$`#name`)
  df["gene_symbol"]<- gene_converter$value[m]
  return(df)
}

var_synonymous<- convert_genes(var_synonymous)
var_misssense<- convert_genes(var_misssense)
var_lof<- convert_genes(var_lof)


write.csv(var_synonymous, "synonymous_singleton_canonical.csv")
write.csv(var_misssense, "misssense_singleton_canonical.csv")
write.csv(var_lof, "lof_singleton_canonical.csv")

# chromosome X ====
#load synonymous 

var_synonymous_x <- read_csv("synonymous/singleX_synonymous.csv")
var_misssense_x <- read_csv("missense/singleX_missense.csv")
var_lof_x <- read_csv("lof/singleX_lof.csv")

edit_df_x<- function(df){
  colnames(df)[1]<- "transcript"
  df<- df[which(df$transcript %in% canonical_transcripts$Transcript),]
  df<- df[,c(-4:-6)]
}

var_synonymous_x<- edit_df_x(var_synonymous_x)
var_misssense_x<- edit_df_x(var_misssense_x)
var_lof_x<- edit_df_x(var_lof_x)

var_synonymous_x<- convert_genes(var_synonymous_x)
var_misssense_x<- convert_genes(var_misssense_x)
var_lof_x<- convert_genes(var_lof_x)

write.csv(var_synonymous_x, "synonymous_singleton_canonical_X.csv")
write.csv(var_misssense_x, "misssense_singleton_canonical_X.csv")
write.csv(var_lof_x, "lof_singleton_canonical_X.csv")

### Fisher Test within each mutation type ----

fisher_fun<- function(df){
  p_val <- apply(df[,6:9], 1, function(x) fisher.test(matrix(c(x[3:4], x[1]-x[3], x[2]-x[4]),2,2))$p.value)
  odds_ratio<- apply(df[,6:9], 1, function(x) fisher.test(matrix(c(x[3:4], x[1]-x[3], x[2]-x[4]),2,2))$estimate)
  
  plot(log(odds_ratio),-log10(p_val),pch=20, cex= 0.4,  main = "synonymous autosomes- fisher test", xlab ="log(odds ratio)", ylab = "-log(p.value)")
  
  return(cbind.data.frame(df[,-1], "p-value" = p_val,"fdr" = p.adjust(p_val, method = "fdr"), 
                                    "odds ratio" = odds_ratio))
}

# Autosomes
syn <- read.csv(file = "synonymous_singleton_canonical.csv")
#filter data with low sum varinats per gene 
s1 <- apply(syn[,8:9],1,sum) 
syn.f <- syn[s1>30,]

syn_fisher_res<- fisher_fun(syn.f)


mis <- read.csv(file = "misssense_singleton_canonical.csv")
#filter data with low sum varinats per gene 
m1 <- apply(mis[,8:9],1,sum) 
mis.f <- mis[m1>30,]

mis_fisher_res<- fisher_fun(mis.f)

lof <- read.csv(file = "lof_singleton_canonical.csv")
#filter data with low sum varinats per gene 
l1 <- apply(lof[,8:9],1,sum) 
lof.f <- lof[l1>30,]

lof_fisher_res<- fisher_fun(lof.f)

# X chromosome
syn_X <- read.csv(file = "synonymous_singleton_canonical_X.csv")
#filter data with low sum varinats per gene 
s1 <- apply(syn_X[,8:9],1,sum) 
syn_X.f <- syn_X[s1>30,]

syn_x_fisher_res<- fisher_fun(syn_X.f)


mis_X <- read.csv(file = "misssense_singleton_canonical_X.csv")
#filter data with low sum varinats per gene 
m1 <- apply(mis_X[,8:9],1,sum) 
mis_X.f <- mis_X[m1>30,]

mis_X_fisher_res<- fisher_fun(mis_X.f)

lof_X <- read.csv(file = "lof_singleton_canonical_X.csv")
#filter data with low sum varinats per gene 
l1 <- apply(lof_X[,8:9],1,sum) 
lof_X.f <- lof_X[l1>10,]

lof_X_fisher_res<- fisher_fun(lof_X.f)

# Fisher Test of non-synonymous/synonymous








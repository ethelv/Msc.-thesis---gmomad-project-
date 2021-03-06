#load packages ====
library("rafalib")
library("ggpubr")

#load data for autosomes ====
setwd("~/Documents/Research/all chrom/all_var_analysis/singletones/synonymous/new")
syn<- read.csv("all_chrom_new.csv")
setwd("~/Documents/Research/all chrom/all_var_analysis/singletones/missense/new")
mis<- read.csv("all_chrom_new.csv")
setwd("~/Documents/Research/all chrom/all_var_analysis/singletones/lof/new")
lof<- read.csv("all_chrom_new.csv")
#load data for chromosome X ====
setwd("~/Documents/Research/all chrom/all_var_analysis/singletones/synonymous")
synX<- read.csv("singleX_synonymousnew.csv")
setwd("~/Documents/Research/all chrom/all_var_analysis/singletones/missense")
misX<- read.csv("singleX_missensenew.csv")
setwd("~/Documents/Research/all chrom/all_var_analysis/singletones/lof")
lofX<- read.csv("singleX_lofnew.csv")

#drop the first column from the autosomes
drop_col<- function(df, col_num){
  return(df[,-col_num])
}

syn<- drop_col(syn,1)
mis<- drop_col(mis,1)
lof<- drop_col(lof,1)

#function to chane name of a column
col_name_change<- function(df, col_num, new_name) {
  colnames(df)[col_num]<- new_name
  return(df)
}

df_list<- list(syn, synX, mis, misX, lof, lofX)
res<- lapply(df_list, function(x) col_name_change(x, 1, "gene"))


# syn<- col_name_change(syn, 2, "gene")
# synX<- col_name_change(synX, 2, "gene")
# mis<- col_name_change(mis, 2, "gene")
# misX<- col_name_change(misX, 2, "gene")
# lof<- col_name_change(lof, 2, "gene")
# lofX<- col_name_change(lofX, 2, "gene")


#load list of canonical transcripts
setwd("~/Documents/Research/all chrom")
canonical_transcript_list<- read.csv("canonical transcript list.csv", header = F)
colnames(canonical_transcript_list)<- c("Chrom", "Start", "End", "Number", "Transcript", "Gene" )

m<- unlist(strsplit(as.character(canonical_transcript_list$Transcript),".",fixed = TRUE))
canonical_transcript_list$Transcript<- m[seq(1, length(m), 2)]

#choose only canonical transcripts for the analysis
choose_canonical<- function(df_to_filter, canonical_list){
  m<- df_to_filter$gene %in% canonical_list
  return(df_to_filter[m,])
}

res_canonical<- lapply(res, function(x) choose_canonical(x, canonical_transcript_list$Transcript))

#linear plot autosomes ====
lin_plot_syn<- ggscatter(main ="Synonymous", res_canonical[[1]], x = "AC_male", y = "AC_female", size = 1,
                         add = "reg.line", conf.int = TRUE,  color = "springgreen3",
                         cor.coef = TRUE, cor.method = "pearson",
                         xlab = "Male Observed", ylab = "Female Observed",
                         add.params = list(color = "grey60", fill = "lightgray", size = 0.5))#+ 
#  annotate("text", x = 260, y = 1300, label = paste0("n = ", nrow(res_canonical[[1]])))

lin_plot_mis<- ggscatter(main ="Missense", res_canonical[[3]], x = "AC_male", y = "AC_female", size = 1,
                         add = "reg.line", conf.int = TRUE, color = "orange1",
                         cor.coef = TRUE, cor.method = "pearson",
                         xlab = "Male Observed", ylab = "Female Observed",
                         add.params = list(color = "grey60", fill = "lightgray", size = 0.5))#+ 
#annotate("text", x = 500, y = 3600, label = paste0("n = ", nrow(res_canonical[[3]])))

lin_plot_lof<- ggscatter(main ="LoF", res_canonical[[5]], x = "AC_male", y = "AC_female", size = 1,
                         add = "reg.line", conf.int = TRUE, color = "brown3",
                         cor.coef = TRUE, cor.method = "pearson",
                         xlab = "Male Observed", ylab = "Female Observed",
                         add.params = list(color = "grey60", fill = "lightgray", size = 0.5))#+ 
# annotate("text", x = 25, y = 120, label = paste0("n = ", nrow(res_canonical[[5]])))

lin_plot<- ggarrange(lin_plot_syn, lin_plot_mis, lin_plot_lof, 
                     labels = c("A", "B", "C"), 
                     ncol = 3, nrow = 1)
lin_plot
annotate_figure(lin_plot, top = text_grob("Male-Female Correlarion Autosomes",size = 16))

#linear plot Chromosome X ====
lin_plot_synX<- ggscatter(main ="Synonymous", res_canonical[[2]], x = "AC_male", y = "AC_female", size = 1,
                          add = "reg.line", conf.int = TRUE,  color = "springgreen3",
                          cor.coef = TRUE, cor.method = "pearson",
                          xlab = "Male Observed", ylab = "Female Observed",
                          add.params = list(color = "grey60", fill = "lightgray", size = 0.5))

lin_plot_misX<- ggscatter(main ="Missense", res_canonical[[4]], x = "AC_male", y = "AC_female", size = 1,
                          add = "reg.line", conf.int = TRUE, color = "orange1",
                          cor.coef = TRUE, cor.method = "pearson",
                          xlab = "Male Observed", ylab = "Female Observed",
                          add.params = list(color = "grey60", fill = "lightgray", size = 0.5))

lin_plot_lofX<- ggscatter(main ="LoF", res_canonical[[6]], x = "AC_male", y = "AC_female", size = 1,
                          add = "reg.line", conf.int = TRUE, color = "brown3",
                          cor.coef = TRUE, cor.method = "pearson",
                          xlab = "Male Observed", ylab = "Female Observed",
                          add.params = list(color = "grey60", fill = "lightgray", size = 0.5))

lin_plotX<- ggarrange(lin_plot_synX, lin_plot_misX, lin_plot_lofX, 
                      labels = c("A", "B", "C"), 
                      ncol = 3, nrow = 1)
annotate_figure(lin_plotX, top = text_grob("Male-Female Correlarion Chromosome X",size = 16))
lin_plotX

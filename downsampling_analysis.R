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


 syn<- col_name_change(syn, 2, "gene")
 synX<- col_name_change(synX, 2, "gene")
 mis<- col_name_change(mis, 2, "gene")
 misX<- col_name_change(misX, 2, "gene")
 lof<- col_name_change(lof, 2, "gene")
 lofX<- col_name_change(lofX, 2, "gene")



#downsampling to see if the difference in r_2 arises from the difference in number of samples between the groups
downsampling_test<- function(df, y_axis) {
  x <- y_axis/ sum(df$count)
  r_sq<- c()
  for (val in x) {
    print(val)
    new_AC<- trunc(df$count* val)
    new_AC_female<- c()
    for (i in c(1:length(new_AC))){
      new_AC_female[i]<- rhyper(1,df$AC_female[i], df$AC_male[i],new_AC[i])
    }
  new_AC_male  = new_AC - new_AC_female
  r_sq[match(val,x)]<- summary(lm(new_AC_male ~ new_AC_female))$r.squared
}
return(r_sq)

}

y <- c(10000,25000,50000,75000,100000,150000,200000,300000,400000,500000,600000,700000,800000 )
y_lof<- c(10000,25000,50000,75000,100000,120663)

r_sq_syn<- downsampling_test(res_canonical[[1]],y)
r_sq_mis<- downsampling_test(res_canonical[[3]],y)
r_sq_lof<- downsampling_test(res_canonical[[5]], y_lof)

y_x<- c(1000,2500,3500,5000,8000,9500,11000,14500,17000, 22000, 25000)
y_x_lof <-  c(100,250,500,750,1000,1250,1500)

r_sq_syn_x<- downsampling_test(res_canonical[[2]],y_x)
r_sq_mis_x<- downsampling_test(res_canonical[[4]],y_x)
r_sq_lof_x<- downsampling_test(res_canonical[[6]], y_x_lof)

#plot all types of variants----
plot( y, r_sq_syn, ylim = c(0, 1), xlab = "sample size", ylab = "r-squared", type="b", col="red",lwd=1, pch=19 , xlim = c(0,800000),
      main = "r^2 male-female (downsampling)")
lines(y, r_sq_mis, type="b", col="green",lwd=1, pch=19)
lines(y_lof, r_sq_lof, type="b", col="blue",lwd=1, pch=19)
legend("bottomright",  legend = c("syn", "miss", "lof"), lty=c(1,1),col = c("red", "green", "blue"))
